"""What a Eukore package holds, as a resource (since 2026-11-12): the tables a paper describes a database by, made
from a finished package and NCBI's whole taxonomy (compare.py runs them; benchmarks/README.md):

- breadth_ranks.tsv: per group and rank (species .. phylum) NCBI's named taxa, those NCBI holds a genome of (the
  build's selection_taxa.tsv), and Eukore's (with a genome; with any record): the breadth above the species;
- content.tsv: per group and marker kind (BUSCO genes, 18S, ITS, 28S, a whole operon, mitochondrial and plastid genes,
  and their combinations) the species holding it, and those holding it from their genome's own assembly: the markers
  of one species that are linked through one assembly;
- panels.tsv: the organisms of named lists (benchmarks/panels.tsv: WHO's fungal priority pathogens, human parasites,
  model organisms, the mock communities' eukaryotes), each found in NCBI's taxonomy by its name, and what the package
  holds of it.

    python3 -m eukore_bench.resource --package OUT/Eukore_2026.11 --bench BENCH --out RESULTS
"""

import argparse
import collections
import json
import os
import sys

from . import ROOT
from . import prep
from eukore import benchmarks, groups, labels, ncbi as ncbi_mod, rrna
from eukore.genomes import HYBRID

RANKS = ("species", "genus", "family", "order", "class", "phylum")
ABOVE = RANKS[1:]
RANK_COLUMNS = ("supergroup", "rank", "ncbi_named", "ncbi_genome", "eukore_genome", "eukore")
# selection_taxa.tsv's fates that are no species NCBI holds a genome of
NO_SPECIES = ("not a species", "placeholder", "hybrid", "not in the taxdump")
EUKARYOTA = 2759
PANELS = os.path.join(ROOT, "benchmarks", "panels.tsv")
PANEL_COLUMNS = ("panel", "entry", "taxid", "ncbi_name", "rank", "status", "eukore_genome", "eukore", "ncbi_genome",
                 "loci", "rrna_records", "organelle_genes", "note")
CONTENT_COLUMNS = ("supergroup", "of", "marker", "species", "same_assembly")
ITS = ("ITS1", "5.8S", "ITS2")
MARKERS = ("BUSCO", "18S", "ITS", "28S", "operon", "mitochondrion", "plastid")
# the combinations a reader links markers by: a genome's genes with its barcodes
COMBINATIONS = (("BUSCO + 18S", ("BUSCO", "18S")), ("BUSCO + operon", ("BUSCO", "operon")),
                ("BUSCO + operon + organelle", ("BUSCO", "operon", "organelle")))


class Lineages:
    """A taxon's group (groups.py's, None outside the eukaryotes) and its taxids at genus .. phylum, the taxa above
    the species kept (NCBI's whole taxonomy: two million species, a few hundred thousand taxa above them)."""

    def __init__(self, ncbi):
        self.ncbi, self.table, self.cache = ncbi, groups.rows(), {}

    def of(self, taxid):
        """(group or None, (genus, family, order, class, phylum taxids, None where NCBI names none))."""
        ncbi = self.ncbi
        taxid = ncbi.merged.get(taxid, taxid)
        path, t = [], taxid
        while t not in self.cache:
            parent = ncbi.parent.get(t)
            if parent is None or parent == t or len(path) > 200:  # the root, or a taxon the taxonomy does not hold
                self.cache[t] = (None, (None,) * len(ABOVE))
                break
            path.append(t)
            t = parent
        got = self.cache[t]
        for node in reversed(path):
            group, ids = got
            if node == EUKARYOTA:
                group = groups.OTHER
            if node in self.table and group is not None:
                group = self.table[node][2]
            rank = ncbi.rank.get(node)
            if rank in ABOVE:
                ids = tuple(node if r == rank else x for r, x in zip(ABOVE, ids))
            got = (group, ids)
            if rank != "species" and node != taxid:  # the taxa above the species alone
                self.cache[node] = got
        return got


def named_species(ncbi, taxid):
    """Whether an NCBI species is a named one: a binomial, no placeholder (uncultured ..), sp., cf. or hybrid."""
    name = ncbi.name.get(taxid, "")
    return labels.is_species_name(name) and not ncbi.is_placeholder(taxid) and not HYBRID.search(name)


def rank_rows(species, taxa, ncbi, lineages, whole=True):
    """breadth_ranks.tsv's rows. species: species.tsv's rows by taxid; taxa: selection_taxa.tsv's rows (none: a
    package before 2026-11-12, `ncbi_genome` "-"); whole: ncbi is NCBI's whole taxonomy (else `ncbi_named` "-": a
    package's own taxonomy holds its own lineages alone). A taxon above the species counts in every group one of its
    species is of (Chordata: Vertebrata and Metazoa), once in `all`."""
    held = {k: collections.defaultdict(set) for k in RANK_COLUMNS[2:]}

    def add(kind, group, taxid, is_species):
        _g, ids = lineages.of(taxid)
        for g in (group, "all"):
            if is_species:
                held[kind][(g, "species")].add(taxid)
            for rank, t in zip(ABOVE, ids):
                if t:
                    held[kind][(g, rank)].add(t)

    if whole:
        for taxid, rank in ncbi.rank.items():
            if rank == "species":
                group, _ids = lineages.of(taxid)
                if group is not None and named_species(ncbi, taxid):
                    add("ncbi_named", group, taxid, True)
    for r in taxa or ():
        if r["fate"] not in NO_SPECIES:
            add("ncbi_genome", r["supergroup"], ncbi.merged.get(int(r["taxid"]), int(r["taxid"])), True)
    for taxid, r in species.items():
        t = ncbi.merged.get(int(taxid), int(taxid))
        is_species = ncbi.rank.get(t) == "species"  # a published read labelled by its genus: no species
        add("eukore", r["supergroup"], t, is_species)
        if r["representative"] != "-":
            add("eukore_genome", r["supergroup"], t, is_species)
    names = {g for kind in held.values() for g, _r in kind} - {"all"}
    rows = []
    for g in sorted(names) + ["all"]:
        for rank in RANKS:
            rows.append(dict(supergroup=g, rank=rank,
                             ncbi_named=len(held["ncbi_named"][(g, rank)]) if whole else "-",
                             ncbi_genome=len(held["ncbi_genome"][(g, rank)]) if taxa else "-",
                             eukore_genome=len(held["eukore_genome"][(g, rank)]),
                             eukore=len(held["eukore"][(g, rank)])))
    return rows


def markers_held(folder, species):
    """{marker kind: (the species holding it from any source, those holding it from their genome's own assembly)}
    (taxids as species.tsv writes them): content_rows' sets, MARKERS', `organelle` and the COMBINATIONS'."""
    rep = {t: r["representative"] for t, r in species.items()}
    held = collections.defaultdict(lambda: [set(), set()])  # marker -> (species: any source, their own assembly)

    def add(marker, taxid, same):
        held[marker][0].add(taxid)
        if same:
            held[marker][1].add(taxid)

    for t, r in species.items():
        if benchmarks.count_of(r.get("loci")) > 0:
            add("BUSCO", t, True)
    for r in benchmarks.rows_of(os.path.join(folder, "rrna", "records.tsv")):
        t = r["species_taxid"]
        if t not in species:
            continue
        same = rep[t] != "-" and r.get("accession") == rep[t]
        regions = {name for name, _a, _b in rrna.parse_regions(r["regions"]) or ()}
        if r["kind"] == "18S":
            regions.add("18S")
        for marker, needed in (("18S", ("18S",)), ("ITS", ITS), ("28S", ("28S",))):
            if all(n in regions for n in needed):
                add(marker, t, same)
        if r["kind"] == "rDNA":  # the package's own operon (composition.tsv's with_operon)
            add("operon", t, same)
    for r in benchmarks.rows_of(os.path.join(folder, "organelle", "genes.tsv")):
        t = r["species_taxid"]
        if t in species and r["organelle"] in ("mitochondrion", "plastid"):
            same = rep[t] != "-" and r.get("accession") == rep[t]
            add(r["organelle"], t, same)
            add("organelle", t, same)
    for name, parts in COMBINATIONS:
        for i in (0, 1):
            held[name][i] = set.intersection(*(held[p][i] for p in parts))
    return held


def content_rows(folder, species):
    """content.tsv's rows: per group (and `all`) and marker kind the species holding it (from any source: their
    genome, a RefSeq organelle genome, a published read) and those holding it from their genome's own assembly
    (rrna/records.tsv's and organelle/genes.tsv's `accession` the species' representative). BUSCO: a BUSCO gene
    called in the genome; ITS: ITS1, 5.8S and ITS2 on one record; operon: a record of kind rDNA (a unit from its
    18S to its 28S, as composition.tsv counts them)."""
    held = markers_held(folder, species)
    group = {t: r["supergroup"] for t, r in species.items()}
    rows = []
    for g in sorted(set(group.values())) + ["all"]:
        mine = {t for t, x in group.items() if g in (x, "all")}
        for marker in MARKERS + tuple(n for n, _p in COMBINATIONS):
            rows.append(dict(supergroup=g, of=len(mine), marker=marker, species=len(held[marker][0] & mine),
                             same_assembly=len(held[marker][1] & mine)))
    return rows


# --- the named lists -------------------------------------------------------------------------------------------------

def read_panels(path=PANELS):
    """benchmarks/panels.tsv -> [(panel, name, taxid or None, note)]: its lines of panel, name, an NCBI taxid where
    the name alone does not find the taxon, a note."""
    out = []
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if line.startswith("#") or len(f) < 2 or f[0] == "panel" or not f[1]:
                continue
            out.append((f[0], f[1], int(f[2]) if len(f) > 2 and f[2].isdigit() else None, f[3] if len(f) > 3 else ""))
    return out


def name_key(name):
    """A name as the lists and NCBI's names are compared: labels.key, NCBI's brackets ([Candida] auris) dropped."""
    return labels.key(name.replace("[", "").replace("]", ""))


def resolve(names, ncbi, lineages):
    """{name: NCBI taxid} of the eukaryote taxa named so (a scientific name or a synonym; NCBI's whole taxonomy read
    once for these names alone). A name of several taxa is the one's that holds the others (a genus and its nominate
    subgenus), else a scientific name's before a synonym's; unrelated taxa's still: not found."""
    wanted = {name_key(n): n for n in names}
    found = collections.defaultdict(list)  # key -> [(taxid, by its scientific name)]
    for taxid, name in ncbi.name.items():
        for i, n in enumerate([name] + ncbi.synonyms.get(taxid, [])):
            k = name_key(n)
            if k in wanted and lineages.of(taxid)[0] is not None:
                found[k].append((taxid, i == 0))
    out = {}
    for k, hits in found.items():
        taxa = {t for t, _s in hits}
        if len(taxa) > 1:
            holder = set.intersection(*(set(ncbi.ancestors(t)) for t in taxa)) & taxa
            scientific = {t for t, s in hits if s}
            taxa = holder or (scientific if len(scientific) == 1 else taxa)
        if len(taxa) == 1:
            out[wanted[k]] = taxa.pop()
    return out


def panel_rows(entries, species, taxa, ncbi, lineages):
    """panels.tsv's rows: each entry of the lists with its NCBI taxon and what the package holds of it: `status`
    (`genome`: a species of it with a genome; `records`: with rRNA records alone; `absent`; `not found`: no NCBI
    taxon of that name), the package's species of it (with a genome; any record), NCBI's species of it with a genome
    (selection_taxa.tsv), and those species' BUSCO genes called (the most of one species), rRNA records and
    organelle genes."""
    found = resolve([name for _p, name, t, _n in entries if t is None], ncbi, lineages)
    wanted = {ncbi.merged.get(t, t) for _p, name, t, _n in entries for t in [t or found.get(name)] if t}
    under = collections.defaultdict(list)  # a list's taxon -> the package's rows of it
    for taxid, r in species.items():
        for t in ncbi.ancestors(ncbi.merged.get(int(taxid), int(taxid))):
            if t in wanted:
                under[t].append(r)
    with_genome = collections.Counter()
    for r in taxa or ():
        if r["fate"] not in NO_SPECIES:
            t0 = int(r["taxid"])
            for t in ncbi.ancestors(ncbi.merged.get(t0, t0)):
                if t in wanted:
                    with_genome[t] += 1
    rows = []
    for panel, name, taxid, note in entries:
        taxid = taxid or found.get(name)
        taxid = ncbi.merged.get(taxid, taxid) if taxid else None
        if not taxid or taxid not in ncbi.parent:
            rows.append(dict(panel=panel, entry=name, taxid="-", ncbi_name="-", rank="-", status="not found",
                             eukore_genome=0, eukore=0, ncbi_genome="-", loci=0, rrna_records=0, organelle_genes=0,
                             note=note or "-"))
            continue
        mine = under.get(taxid, [])
        genomes = [r for r in mine if r["representative"] != "-"]
        rows.append(dict(
            panel=panel, entry=name, taxid=taxid, ncbi_name=ncbi.name.get(taxid, "-"), rank=ncbi.rank.get(taxid, "-"),
            status="genome" if genomes else "records" if mine else "absent", eukore_genome=len(genomes),
            eukore=len(mine), ncbi_genome=with_genome[taxid] if taxa else "-",
            loci=max((benchmarks.count_of(r.get("loci")) for r in mine), default=0),
            rrna_records=sum(benchmarks.count_of(r.get("rrna_records")) for r in mine),
            organelle_genes=sum(benchmarks.count_of(r.get("organelle_genes")) for r in mine), note=note or "-"))
    return rows


def panel_summary(rows):
    """{panel: its entries, those with a genome, with records alone, absent, not found}."""
    out = {}
    for r in rows:
        s = out.setdefault(r["panel"], {"entries": 0, "genome": 0, "records": 0, "absent": 0, "not found": 0})
        s["entries"] += 1
        s[r["status"]] += 1
    return out


# --- all of it -------------------------------------------------------------------------------------------------------

def run(package, out, ncbi, species, summary, whole=True, panels=PANELS, log=print):
    """The resource's tables of package into out, their headline figures into summary (breadth_ranks, content,
    panels). species: compare.species_of's rows; whole: ncbi is NCBI's whole taxonomy."""
    lineages = Lineages(ncbi)
    taxa = list(benchmarks.rows_of(os.path.join(out, "selection_taxa.tsv")))
    if not taxa:
        summary.setdefault("left_out", {})["ncbi_genome"] = \
            "no selection_taxa.tsv (a package built before 2026-11-12): NCBI's taxa with a genome not counted by rank"
    rows = rank_rows(species, taxa, ncbi, lineages, whole)
    benchmarks.write(out, "breadth_ranks.tsv", RANK_COLUMNS, rows)
    summary["breadth_ranks"] = {r["rank"]: {k: r[k] for k in RANK_COLUMNS[2:]} for r in rows
                                if r["supergroup"] == "all"}
    rows = content_rows(package, species)
    benchmarks.write(out, "content.tsv", CONTENT_COLUMNS, rows)
    summary["content"] = {r["marker"]: {"species": r["species"], "same_assembly": r["same_assembly"]}
                          for r in rows if r["supergroup"] == "all"}
    entries = read_panels(panels)
    if entries:
        rows = panel_rows(entries, species, taxa, ncbi, lineages)
        benchmarks.write(out, "panels.tsv", PANEL_COLUMNS, rows)
        summary["panels"] = panel_summary(rows)
        missing = [r["entry"] for r in rows if r["status"] == "not found"]
        if missing:
            log(f"WARNING: the lists' names NCBI's taxonomy does not hold (benchmarks/panels.tsv: give their taxid): "
                f"{', '.join(missing)}")
    log("resource: " + ", ".join(
        f"{rank} {v['eukore_genome']:,} of {v['ncbi_genome']} with a genome in NCBI"
        for rank, v in summary["breadth_ranks"].items() if rank in ("species", "genus", "family")))
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--bench", help="prep.py's folder (BENCH): its taxdump")
    p.add_argument("--out", required=True, help="the results' folder (compare.py's: its summary.json is added to)")
    p.add_argument("--taxdump", help="NCBI's whole taxdump (default: BENCH/downloads/taxdump.tar.gz)")
    p.add_argument("--panels", default=PANELS, help="the named lists (default: benchmarks/panels.tsv)")
    a = p.parse_args(argv)
    from . import compare  # noqa: PLC0415
    os.makedirs(a.out, exist_ok=True)
    compare.open_results(a.package, a.out)
    whole = prep.taxdump_of(a.bench, a.taxdump) if a.bench or a.taxdump else None
    ncbi = ncbi_mod.Ncbi(whole or os.path.join(a.package, "taxonomy"))
    summary = compare.load_summary(a.out)
    run(a.package, a.out, ncbi, compare.species_of(a.package), summary, bool(whole), a.panels,
        log=lambda text: print(text, file=sys.stderr, flush=True))
    compare.write_summary(a.package, a.out, summary)
    json.dump({k: summary[k] for k in ("breadth_ranks", "content", "panels") if k in summary}, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
