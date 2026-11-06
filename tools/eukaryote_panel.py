#!/usr/bin/env python3
"""How often Eukore's caller finds each BUSCO across eukaryotes (since 2026-11-05): a panel of NCBI genomes, one per order
in each supergroup (groups.GROUPS) up to --per_group (orders spread over the supergroup's classes), each the order's
smallest representative genome (genomes.representatives) of --max_mb or fewer that carries no flag (NCBI's, a partial
genome, small for its clade) and has MAX_CONTIGS_PER_MB or fewer contigs per Mb, in which every BUSCO of the dataset is called as a build calls
it (loci.call with the 30 proteins of the nearest species and the ancestral variants). BUSCO's own dataset holds few protists
(22 SAR, 6 Discoba, 6 Amoebozoa species of 837), so its species list overstates how universal a BUSCO is outside
animals, fungi and plants; this panel measures it in real genomes. tools/rank_loci.py --eukaryote_panel OUT reads it.

    python3 tools/eukaryote_panel.py --assemblies WORK/sources/ncbi_assemblies --taxdump WORK/sources/taxdump.tar.gz \\
        --odb WORK/sources/eukaryota_odb12 --out PANEL [--per_group 20 --max_mb 1500 -j 2 -t 3]

PANEL/panel.tsv lists the genomes chosen; PANEL/done/<accession>.tsv each one's calls (locus, status, length,
hmm_score); a genome is deleted once called (--keep: not). A rerun calls only those not done.
"""

import argparse
import collections
import functools
import glob
import multiprocessing
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import genomes, loci, ncbi as ncbi_mod, odb as odb_mod, placement, sources, util  # noqa: E402
from eukore.groups import group_of  # noqa: E402 (the supergroups, groups.GROUPS: the package's)

MAX_CONTIGS_PER_MB = 100  # more: contigs of ~10 kb, whose genes are cut at their ends


def choose(ncbi, assemblies, per_group, max_mb):
    """The panel -> [(group, class, order, species taxid, assembly row)]."""
    reps = genomes.representatives(assemblies, ncbi)
    by_order = collections.defaultdict(list)
    for sp, (row, _others) in reps.items():
        at = ncbi.ranks(sp)
        # a flagged genome (NCBI's flags, partial, small for its clade) or a fragmented one would count its gaps as
        # missing BUSCOs; of the rest, the smallest under max_mb (the run's cost)
        if not at.get("order") or row.get("flags") or not row.get("size") or row["size"] > max_mb * 1e6 \
                or (row.get("contigs") or 1 << 20) / (row["size"] / 1e6) > MAX_CONTIGS_PER_MB:
            continue
        by_order[(group_of(ncbi, sp), at.get("class") or "-", at["order"])].append((row["size"], sp, row))
    by_group = collections.defaultdict(list)
    for (grp, cls, order), rows in sorted(by_order.items()):
        _rank, sp, row = min(rows, key=lambda x: x[0])
        by_group[grp].append((grp, cls, order, sp, row))
    out = []
    for grp, rows in sorted(by_group.items()):
        if len(rows) <= per_group:
            out += rows
            continue
        # spread over the classes: round-robin over them, each class's orders in name order
        by_class = collections.defaultdict(list)
        for r in rows:
            by_class[r[1]].append(r)
        queues = [collections.deque(v) for _c, v in sorted(by_class.items())]
        while len([r for r in out if r[0] == grp]) < per_group and any(queues):
            for q in queues:
                if q and len([r for r in out if r[0] == grp]) < per_group:
                    out.append(q.popleft())
    return out


def call_one(item, settings):
    grp, cls, order, sp, row = item
    acc = row["accession"]
    done = os.path.join(settings["out"], "done", acc + ".tsv")
    if os.path.isfile(done):
        return acc, None
    path = os.path.join(settings["out"], "genomes", acc + ".fna.gz")
    if not os.path.isfile(path):
        why = util.fetch_file(placement.genome_url(row), path)
        if why:
            return acc, f"fetch failed: {why}"
    odb = odb_mod.Odb(settings["odb"], settings["loci"])
    got = loci.call(path, odb.queries(settings["lineage"][sp], settings["odb_lineages"]), odb,
                    threads=settings["threads"], tmp=os.path.dirname(done), table=settings["gencode"][sp])
    if got["failed"]:
        return acc, got["failed"]
    with open(done + ".part", "w", newline="\n") as fh:
        fh.write(f"#{acc}\t{grp}\t{cls}\t{order}\t{sp}\t{row['organism']}\n")
        fh.write("locus\tstatus\tlength\thmm_score\n")
        for locus in odb.loci:
            r = got["loci"][locus]
            fh.write(f"{locus}\t{r['status']}\t{len(r.get('copy', ''))}\t{r.get('hmm_score', 0):.1f}\n")
    os.replace(done + ".part", done)
    if not settings["keep"]:
        os.remove(path)
    return acc, None


def read(folder):
    """A panel's calls -> {accession: (group, order, {locus: status})}."""
    out = {}
    for path in glob.glob(os.path.join(folder, "done", "*.tsv")):
        with open(path) as fh:
            head = fh.readline().rstrip("\n").lstrip("#").split("\t")
            fh.readline()
            calls = {}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                calls[f[0]] = f[1]
        out[head[0]] = (head[1], head[3], calls)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--assemblies", required=True, help="NCBI's assembly tables (genbank_*.txt)")
    p.add_argument("--taxdump", required=True)
    p.add_argument("--odb", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--per_group", type=int, default=20)
    p.add_argument("--max_mb", type=float, default=1500)
    p.add_argument("--keep", action="store_true", help="keep the genomes once called")
    p.add_argument("-j", "--jobs", type=int, default=2)
    p.add_argument("-t", "--threads", type=int, default=3)
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(opts.taxdump)
    assemblies = sources.read_assemblies(sorted(glob.glob(os.path.join(opts.assemblies, "genbank_*.txt"))))
    items = choose(ncbi, assemblies, opts.per_group, opts.max_mb)
    os.makedirs(os.path.join(opts.out, "done"), exist_ok=True)
    os.makedirs(os.path.join(opts.out, "genomes"), exist_ok=True)
    with open(os.path.join(opts.out, "panel.tsv"), "w", newline="\n") as fh:
        fh.write("accession\tgroup\tclass\torder\tspecies_taxid\torganism\tsize\n")
        fh.writelines(f"{r['accession']}\t{g}\t{c}\t{o}\t{sp}\t{r['organism']}\t{r['size']}\n" for g, c, o, sp, r in items)
    with open(os.path.join(opts.odb, "scores_cutoff")) as fh:
        all_loci = [line.split()[0] for line in fh if line.strip()]
    odb = odb_mod.Odb(opts.odb, all_loci)
    settings = dict(out=opts.out, odb=opts.odb, loci=all_loci, threads=opts.threads, keep=opts.keep,
                    odb_lineages=odb_mod.species_lineages(odb, ncbi),
                    lineage={sp: list(ncbi.ancestors(sp)) for _g, _c, _o, sp, _r in items},
                    gencode={sp: ncbi.code_of(int(r.get("taxid") or sp)) for _g, _c, _o, sp, r in items})
    counts = collections.Counter(g for g, *_ in items)
    print(f"panel: {len(items)} genomes, {sum(r['size'] for *_x, r in items) / 1e9:.1f} Gb: {dict(counts)}; "
          f"{len(all_loci)} BUSCOs", flush=True)
    failed = []
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for i, (acc, why) in enumerate(pool.imap_unordered(functools.partial(call_one, settings=settings), items), 1):
            if why:
                failed.append((acc, why))
            if i % 10 == 0 or i == len(items):
                print(f"{i}/{len(items)} called ({len(failed)} failed)", flush=True)
    for acc, why in failed:
        print(f"failed: {acc}: {why}", flush=True)
    print("DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
