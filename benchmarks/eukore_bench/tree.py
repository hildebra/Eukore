"""The species tree of a Eukore package's BUSCO proteins against NCBI's taxonomy (since 2026-11-12; compare.py runs
it where BUSCO's HMMs, hmmalign and FastTree are there; benchmarks/README.md): whether the markers' sequences place
each genome where its name says.

- The genomes: up to PER_GENUS species' genomes per genus, those with the most BUSCOs Single first (a genus is
  tested by two of its species; every species of all eukaryotes is more than a tree of this kind needs).
- The loci: the LOCI BUSCOs Single (and not of confidence low) in the most of those genomes; a genome with under
  MIN_LOCI of them is left out.
- The alignment: each locus' proteins (busco/proteins/) placed on BUSCO's profile by hmmalign (eukore/align.py: a
  residue on a match state is that column's; the build's own rows, WORK/align/, are read where given), the loci's
  columns joined, a column kept where OCCUPANCY of the genomes hold a residue.
- The tree: FastTree (-lg -gamma) of the joined columns.
- tree_taxa.tsv: each NCBI taxon of genus .. phylum with two genomes or more in the tree (and not all of them):
  whether its genomes are one branch of the tree (monophyletic; the tree is unrooted: one side of an edge), else the
  branch that fits it best (the F1 of the branch's genomes and the taxon's), the taxon's genomes outside it
  (`strays`) and the other taxa's inside (`intruders`).
- tree_strays.tsv: the genomes that stray from a taxon of three genomes or more whose other genomes hold together
  (STRAY_SHARE of it at most), with the highest rank they stray at and the genomes of the branch they sit on:
  placed apart from their name's relatives, the genomes to look at first (a misnamed genome, a contaminated one, or
  a taxon NCBI's taxonomy holds together and the genes do not).
- tree.nwk: the tree, its tips the genomes' accessions; tree_tips.tsv names them.

    python3 -m eukore_bench.tree --package OUT/Eukore_2026.11 --odb SOURCES/eukaryota_odb12 --out RESULTS [-t 16]
"""

import argparse
import collections
import concurrent.futures
import json
import os
import shutil
import subprocess
import sys
import tempfile

from . import prep
from eukore import align, benchmarks, ncbi as ncbi_mod, odb as odb_mod, package as package_mod

PER_GENUS = 2
LOCI = 50
MIN_LOCI = 0.5      # of the loci chosen: a genome with fewer Single is not placed
OCCUPANCY = 0.5     # a column held by this share of the genomes
STRAY_SHARE = 0.2   # a taxon's genomes outside its best branch, this share of it at most (one at least): strays
BESIDE = 3          # the genomes named of the branch a stray sits on
RANKS = ("genus", "family", "order", "class", "phylum")
TIP_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "loci") + RANKS
TAXON_COLUMNS = ("rank", "taxid", "name", "supergroup", "genomes", "monophyletic", "f1", "strays", "intruders")
STRAY_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "rank", "taxon", "genomes", "ranks", "beside")


# --- the genomes and loci ------------------------------------------------------------------------------------------

def choose(package, species, ncbi, per_genus=PER_GENUS, n_loci=LOCI):
    """The tree's genomes and loci -> ({accession: (species taxid, {rank: taxid}, its loci Single)}, [locus])."""
    genomes = {r["accession"]: r for r in package_mod.read_table(package, "genomes.tsv")
               if r["representative"] == "t" and r.get("source") != "no_assembly"}
    single = collections.defaultdict(set)
    for r in benchmarks.rows_of(os.path.join(package, "busco", "copies.tsv")):
        if r["accession"] in genomes and r["status"] == "Single" and r.get("confidence") != "low":
            single[r["accession"]].add(r["locus"])
    by_genus = collections.defaultdict(list)
    ranks = {}
    for acc, g in genomes.items():
        sp = int(g["species_taxid"])
        if g["species_taxid"] not in species or not single[acc]:
            continue
        ranks[acc] = {r: t for r, t in ncbi.rank_ids(sp).items() if r in RANKS}
        by_genus[ranks[acc].get("genus") or f"species {sp}"].append(acc)
    chosen = []
    for accs in by_genus.values():
        chosen += sorted(accs, key=lambda a: (-len(single[a]), a))[:max(1, per_genus)]
    held = collections.Counter(locus for a in chosen for locus in single[a])
    loci = sorted(sorted(held, key=lambda x: (-held[x], x))[:n_loci])
    tips = {a: (int(genomes[a]["species_taxid"]), ranks[a], single[a] & set(loci)) for a in chosen
            if len(single[a] & set(loci)) >= max(1, MIN_LOCI * len(loci))}
    return tips, loci


# --- the alignment -----------------------------------------------------------------------------------------------------

def locus_rows(package, locus, hmm, tips, hmmalign, cache, tmp):
    """One locus' proteins of the tips that hold it Single, each as its row on the profile's match states."""
    proteins = {a: p for a, p in package_mod.read_copies(package, locus, "proteins").items()
                if a in tips and locus in tips[a][2] and p}
    states = align.align_locus(hmm, proteins, hmmalign, cache, tmp)
    width = max([odb_mod.hmm_header(hmm)[1]] + [s for row in states.values() for s in row if s])
    rows = {}
    for acc, placed in states.items():
        row = ["-"] * width
        for residue, state in zip(proteins[acc].replace("*", "X"), placed):
            if state:
                row[state - 1] = residue.upper()
        rows[acc] = "".join(row)
    return locus, width, rows


def alignment(package, tips, loci, odb_folder, hmmalign, threads=1, cache=None, tmp=None, log=print):
    """The tips' proteins of the loci, each on its profile's match states and the loci joined, the columns OCCUPANCY
    of the tips hold kept -> ({accession: row}, the columns kept, the columns of all loci)."""
    hmm = {locus: os.path.join(odb_folder, "hmms", f"{locus}.hmm") for locus in loci}
    with concurrent.futures.ThreadPoolExecutor(max(1, threads)) as pool:
        got = list(pool.map(lambda x: locus_rows(package, x, hmm[x], tips, hmmalign, cache, tmp), loci))
    failed = [locus for locus, _w, rows in got if not rows]
    if failed:
        log(f"WARNING: tree: hmmalign placed no protein of {len(failed)} loci (e.g. {failed[:3]})")
    parts = {a: [] for a in tips}
    for _locus, width, rows in got:
        for a in tips:
            parts[a].append(rows.get(a) or "-" * width)
    joined = {a: "".join(p) for a, p in parts.items()}
    total = len(next(iter(joined.values()), ""))
    keep = [i for i in range(total) if sum(row[i] != "-" for row in joined.values()) >= OCCUPANCY * len(joined)]
    return {a: "".join(row[i] for i in keep) for a, row in joined.items()}, len(keep), total


def fasttree(rows, binary, work, threads=1):
    """FastTree's tree (-lg -gamma) of the rows -> Newick text."""
    aln = os.path.join(work, "tree.faa")
    with open(aln, "w") as fh:
        fh.writelines(f">{a}\n{row}\n" for a, row in sorted(rows.items()))
    env = dict(os.environ, OMP_NUM_THREADS=str(max(1, threads)))
    run = subprocess.run([binary, "-lg", "-gamma", "-quiet", aln], capture_output=True, text=True, env=env)
    if run.returncode or not run.stdout.strip():
        raise RuntimeError(f"FastTree failed (exit {run.returncode}): {run.stderr.strip()[-300:]}")
    return run.stdout.strip()


# --- the tree ----------------------------------------------------------------------------------------------------------

class Tree:
    """A Newick tree as lists by node: parent, children, name (a tip's), length (its branch's); node 0 the root."""

    def __init__(self, text):
        self.parent, self.children, self.name, self.length = [-1], [[]], [None], [0.0]
        node, i, text = 0, 0, text.strip()
        while i < len(text):
            c = text[i]
            if c == "(":
                node = self._add(node)
            elif c == ",":
                node = self._add(self.parent[node])
            elif c == ")":
                node = self.parent[node]
            elif c == ";":
                break
            else:
                j = i
                while j < len(text) and text[j] not in "(),;":
                    j += 1
                label, _colon, length = text[i:j].partition(":")
                if label.strip() and not self.children[node]:
                    self.name[node] = label.strip().strip("'")
                try:
                    self.length[node] = float(length) if length else 0.0
                except ValueError:
                    self.length[node] = 0.0
                i = j - 1
            i += 1
        self.tips = [n for n in range(len(self.parent)) if not self.children[n]]

    def _add(self, parent):
        self.parent.append(parent)
        self.children.append([])
        self.name.append(None)
        self.length.append(0.0)
        self.children[parent].append(len(self.parent) - 1)
        return len(self.parent) - 1

    def postorder(self):
        """The nodes, each after its children."""
        order, stack = [], [0]
        while stack:
            n = stack.pop()
            order.append(n)
            stack.extend(self.children[n])
        return order[::-1]

    def masks(self):
        """([each node's tips as a bit set], {tip name: its bit})."""
        bit = {self.name[n]: 1 << i for i, n in enumerate(self.tips)}
        out = [0] * len(self.parent)
        for n in self.postorder():
            out[n] = bit[self.name[n]] if not self.children[n] else 0
            for c in self.children[n]:
                out[n] |= out[c]
        return out, bit


def best_branch(mask, size, masks, everything, total):
    """The side of an edge (a node's tips, or all the others) that fits a taxon's tips (mask, size of them) best ->
    (F1, that side's tips as a bit set)."""
    best = (0.0, 0)
    for c in masks:
        inside = (c & mask).bit_count()
        n = c.bit_count()
        for shared, side_size, side in ((inside, n, c), (size - inside, total - n, everything ^ c)):
            if shared and side_size:
                f1 = 2 * shared / (side_size + size)
                if f1 > best[0]:
                    best = (f1, side)
    return best


def evaluate(tree, tips, ncbi, group_of):
    """tree_taxa.tsv's and tree_strays.tsv's rows. tips: {accession: (species taxid, {rank: taxid}, loci)};
    group_of: {accession: its group}."""
    masks, bit = tree.masks()
    everything, total = masks[0], len(tree.tips)
    seen = set(masks) | {everything ^ m for m in masks}
    name_of_bit = {b: a for a, b in bit.items()}
    members = collections.defaultdict(int)
    for acc, (_sp, ranks, _loci) in tips.items():
        if acc in bit:
            for rank, taxid in ranks.items():
                members[(rank, taxid)] |= bit[acc]
    parent_of_tip = {tree.name[n]: tree.parent[n] for n in tree.tips}
    rows, strays = [], collections.defaultdict(list)  # accession -> [(rank, taxid, genomes)]
    for (rank, taxid), mask in sorted(members.items(), key=lambda kv: (RANKS.index(kv[0][0]), kv[0][1])):
        size = mask.bit_count()
        if size < 2 or size == total:
            continue
        first = name_of_bit[mask & -mask]
        row = dict(rank=rank, taxid=taxid, name=ncbi.name.get(taxid, "-"), supergroup=group_of.get(first, "-"),
                   genomes=size)
        if mask in seen:
            rows.append(dict(row, monophyletic="t", f1="1.000", strays=0, intruders=0))
            continue
        f1, side = best_branch(mask, size, masks, everything, total)
        out, inside = mask & ~side, side & ~mask
        rows.append(dict(row, monophyletic="f", f1=f"{f1:.3f}", strays=out.bit_count(), intruders=inside.bit_count()))
        if size >= 3 and out.bit_count() <= max(1, int(STRAY_SHARE * size)):
            while out:
                low = out & -out
                strays[name_of_bit[low]].append((rank, taxid, size))
                out ^= low
    stray_rows = []
    for acc, found in sorted(strays.items()):
        rank, taxid, size = max(found, key=lambda f: RANKS.index(f[0]))  # the highest rank it strays at
        beside, node = [], parent_of_tip[acc]  # the genomes of the nearest branch that holds others
        while node != -1 and not beside:
            beside = [name_of_bit[b] for b in _bits(masks[node]) if name_of_bit[b] != acc]
            node = tree.parent[node]
        near = [ncbi.name.get(tips[a][0], a) if a in tips else a for a in beside[:BESIDE]]
        sp = tips[acc][0]
        stray_rows.append(dict(accession=acc, species_taxid=sp, name=ncbi.name.get(sp, "-"),
                               supergroup=group_of.get(acc, "-"), rank=rank, taxon=ncbi.name.get(taxid, "-"),
                               genomes=size, ranks=",".join(r for r in RANKS[::-1] if r in {f[0] for f in found}),
                               beside=", ".join(near) + (f" and {len(beside) - BESIDE} more"
                                                         if len(beside) > BESIDE else "") or "-"))
    return rows, stray_rows


def _bits(mask):
    while mask:
        low = mask & -mask
        yield low
        mask ^= low


def rank_summary(rows):
    """{rank: its taxa tested, those monophyletic, their share}."""
    out = {}
    for rank in RANKS:
        mine = [r for r in rows if r["rank"] == rank]
        if mine:
            mono = sum(r["monophyletic"] == "t" for r in mine)
            out[rank] = {"taxa": len(mine), "monophyletic": mono, "share": round(mono / len(mine), 4)}
    return out


# --- all of it -------------------------------------------------------------------------------------------------------

def run(package, out, ncbi, species, summary, odb_folder, hmmalign=None, binary=None, threads=1, cache=None, tmp=None,
        per_genus=PER_GENUS, n_loci=LOCI, log=print):
    """The tree of package and its taxa's tables into out, its headline figures into summary (tree); where a tool or
    BUSCO's profiles are missing, summary's left_out says so."""
    left = summary.setdefault("left_out", {})
    left.pop("tree", None)
    hmmalign = hmmalign or shutil.which("hmmalign")
    binary = binary or shutil.which("FastTreeMP") or shutil.which("FastTree") or shutil.which("fasttree")
    if not odb_folder or not os.path.isdir(os.path.join(odb_folder, "hmms")):
        left["tree"] = "no BUSCO dataset (its hmms/: prep.py's odb step, or --odb)"
    elif not hmmalign or not binary:
        left["tree"] = f"no {'hmmalign' if not hmmalign else 'FastTree'} (the benchmarks' environment)"
    if "tree" in left:
        return summary
    tips, loci = choose(package, species, ncbi, per_genus, n_loci)
    if len(tips) < 4:
        left["tree"] = f"{len(tips)} genomes to place: under 4"
        return summary
    work = tempfile.mkdtemp(prefix="tree.", dir=tmp or out)
    try:
        rows, kept, total = alignment(package, tips, loci, odb_folder, hmmalign, threads, cache, work, log)
        log(f"tree: {len(tips):,} genomes, {len(loci)} loci, {kept:,} of {total:,} columns; FastTree ..")
        newick = fasttree(rows, binary, work, threads)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    with open(os.path.join(out, "tree.nwk"), "w", newline="\n") as fh:
        fh.write(newick + "\n")
    group_of = {a: species.get(str(sp), {}).get("supergroup", "-") for a, (sp, _r, _l) in tips.items()}
    benchmarks.write(out, "tree_tips.tsv", TIP_COLUMNS, [
        dict(accession=a, species_taxid=sp, name=ncbi.name.get(sp, "-"), supergroup=group_of[a], loci=len(held),
             **{r: ncbi.name.get(ranks.get(r), "-") for r in RANKS}) for a, (sp, ranks, held) in sorted(tips.items())])
    taxa, strays = evaluate(Tree(newick), tips, ncbi, group_of)
    benchmarks.write(out, "tree_taxa.tsv", TAXON_COLUMNS, taxa)
    benchmarks.write(out, "tree_strays.tsv", STRAY_COLUMNS, strays)
    summary["tree"] = {"genomes": len(tips), "loci": len(loci), "columns": kept, "per_genus": per_genus,
                       "ranks": rank_summary(taxa), "strays": len(strays)}
    log("tree: " + ", ".join(f"{rank} {v['monophyletic']} of {v['taxa']} monophyletic"
                             for rank, v in summary["tree"]["ranks"].items()) + f"; {len(strays)} genomes stray")
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--out", required=True, help="the results' folder (compare.py's: its summary.json is added to)")
    p.add_argument("--bench", help="prep.py's folder (BENCH): its taxdump and BUSCO dataset")
    p.add_argument("--odb", help="BUSCO's dataset (eukaryota_odb12/, its hmms/; default: BENCH's)")
    p.add_argument("--taxdump", help="NCBI's whole taxdump (default: BENCH's, else the package's taxonomy/)")
    p.add_argument("--align_cache", help="the build's WORK/align (its proteins' rows on the profiles, read and added to)")
    p.add_argument("--hmmalign", help="hmmalign (default: on PATH)")
    p.add_argument("--fasttree", help="FastTree or FastTreeMP (default: on PATH)")
    p.add_argument("--per_genus", type=int, default=PER_GENUS, help="species' genomes per genus")
    p.add_argument("--loci", type=int, default=LOCI, help="the BUSCOs joined: those Single in the most genomes")
    p.add_argument("--tmp", help="the alignment's temporary files (default: in --out)")
    p.add_argument("-t", "--threads", type=int, default=8)
    a = p.parse_args(argv)
    from . import compare  # noqa: PLC0415
    compare.open_results(a.package, a.out)
    whole = prep.taxdump_of(a.bench, a.taxdump) if a.bench or a.taxdump else None
    ncbi = ncbi_mod.Ncbi(whole or os.path.join(a.package, "taxonomy"))
    summary = compare.load_summary(a.out)
    run(a.package, a.out, ncbi, compare.species_of(a.package), summary, a.odb or prep.odb_of(a.bench), a.hmmalign,
        a.fasttree, a.threads, a.align_cache, a.tmp, a.per_genus, a.loci,
        log=lambda text: print(text, file=sys.stderr, flush=True))
    compare.write_summary(a.package, a.out, summary)
    json.dump({k: summary[k] for k in ("tree",) if k in summary} or {"left_out": summary["left_out"].get("tree")},
              sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
