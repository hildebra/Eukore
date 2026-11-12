"""How good a Eukore package's records are, as headline numbers (since 2026-11-12; compare.py runs it;
benchmarks/README.md), from the package's own tables:

- quality.tsv: per group and level (genome, BUSCO copy, rRNA record, organelle gene) each check and confidence tier:
  the genomes or records it flags, of how many. A genome's checks: NCBI's and its FCS-GX and ANI reports' flags, copies
  another clade's, a duplicated genome, a genetic code that differs from NCBI's;
- gencode.tsv: per group the genomes whose genetic code, checked against their BUSCO genes, agrees with NCBI's,
  differs or is undecided, and the tables that fit those that differ;
- rdna_published.tsv: the genomes' rDNA copies by read depth (genomes.tsv's, where tools/rdna_depth.py measured them)
  beside the published numbers of benchmarks/rdna_published.tsv.

    python3 -m eukore_bench.checks --package OUT/Eukore_2026.11 --out RESULTS
"""

import argparse
import collections
import json
import os
import sys

from . import ROOT
from eukore import benchmarks

QUALITY_COLUMNS = ("supergroup", "level", "check", "flagged", "of")
GENCODE_COLUMNS = ("supergroup", "genomes", "agrees", "differs", "undecided", "unchecked", "best_tables")
RDNA_COLUMNS = ("species", "species_taxid", "accession", "low", "high", "reference", "copies_reads",
                "units_assembled", "within")
RDNA_PUBLISHED = os.path.join(ROOT, "benchmarks", "rdna_published.tsv")
LEVELS = ("genome", "BUSCO copy", "rRNA record", "organelle gene")
DUPLICATED = ("haplotypes", "haplotigs", "wgd_or_polyploid")
ANY = "any flag"


def flags_of(text):
    """A flags column's flags (`-` for none; comma-separated)."""
    return [f.strip() for f in (text or "-").split(",") if f.strip() not in ("", "-")]


def quality_rows(package, species):
    """quality.tsv's rows: per group (and `all`), level and check the genomes or records flagged, of the level's in
    the group. species: species.tsv's rows by taxid, their `supergroup` the group."""
    genomes = {r["accession"]: r for r in benchmarks.rows_of(os.path.join(package, "genomes.tsv"))
               if r["representative"] == "t"}
    group_of = {a: species.get(g["species_taxid"], {}).get("supergroup", "-") for a, g in genomes.items()}
    counts, totals = collections.Counter(), collections.Counter()

    def add(level, group, checks):
        for g in (group, "all"):
            totals[(g, level)] += 1
            for c in checks:
                counts[(g, level, c)] += 1

    for acc, g in genomes.items():
        if g.get("source") == "no_assembly":
            continue
        checks = flags_of(g.get("flags"))
        if checks:
            checks.append(ANY)
        if benchmarks.count_of(g.get("foreign_copies")) > 0:
            checks.append("copies another clade's")
        if g.get("duplication") in DUPLICATED:
            checks.append(f"duplication: {g['duplication']}")
        if (g.get("gencode_check") or "").startswith("differs"):
            checks.append("genetic code differs from NCBI's")
        add("genome", group_of[acc], checks)
    for level, table, group in (
            ("BUSCO copy", "busco/copies.tsv", lambda r: group_of.get(r["accession"])),
            ("rRNA record", "rrna/records.tsv", lambda r: species.get(r["species_taxid"], {}).get("supergroup")),
            ("organelle gene", "organelle/genes.tsv",
             lambda r: species.get(r["species_taxid"], {}).get("supergroup"))):
        for r in benchmarks.rows_of(os.path.join(package, *table.split("/"))):
            g = group(r)
            if g is None:  # another assembly's copy, a record of no species of the package
                continue
            add(level, g, flags_of(r.get("flags")) + [f"confidence: {r.get('confidence') or '-'}"])
    rows = []
    for (g, level), n in sorted(totals.items(), key=lambda kv: (kv[0][0] == "all", kv[0][0], LEVELS.index(kv[0][1]))):
        checks = sorted((c for (g2, l2, c) in counts if (g2, l2) == (g, level)),
                        key=lambda c: (c.startswith("confidence"), -counts[(g, level, c)], c))
        rows += [dict(supergroup=g, level=level, check=c, flagged=counts[(g, level, c)], of=n) for c in checks]
    return rows


def gencode_rows(package, species):
    """gencode.tsv's rows: per group (and `all`) the species' genomes by their code's check (genomes.tsv's
    gencode_check: agrees, differs, undecided; `-`: not checked) and, of those that differ, the tables that fit
    (`best_table` with its genomes: `6:3; 10:1`)."""
    counts, tables = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
    for r in benchmarks.rows_of(os.path.join(package, "genomes.tsv")):
        if r["representative"] != "t" or r.get("source") == "no_assembly":
            continue
        check = (r.get("gencode_check") or "-").split(":")[0]
        check = check if check in ("agrees", "differs", "undecided") else "unchecked"
        for g in (species.get(r["species_taxid"], {}).get("supergroup", "-"), "all"):
            counts[g][check] += 1
            if check == "differs" and r.get("best_table", "-") != "-":
                tables[g][r["best_table"]] += 1
    return [dict(supergroup=g, genomes=sum(c.values()), agrees=c["agrees"], differs=c["differs"],
                 undecided=c["undecided"], unchecked=c["unchecked"],
                 best_tables="; ".join(f"{t}:{n}" for t, n in sorted(tables[g].items(), key=lambda kv: -kv[1])) or "-")
            for g, c in sorted(counts.items(), key=lambda kv: (kv[0] == "all", kv[0]))]


def rdna_rows(package, species, published=RDNA_PUBLISHED):
    """rdna_published.tsv's rows: each species of the published list the package holds a genome of, with its rDNA
    copies by read depth (`-` where not measured), the units its assembly holds, and whether the reads' number lies
    within the published range (`t`, `f`, `-`)."""
    if not os.path.isfile(published):
        return []
    with open(published) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    by_name = {r.get("name"): (t, r) for t, r in species.items() if r.get("representative", "-") != "-"}
    genomes = {r["accession"]: r for r in benchmarks.rows_of(os.path.join(package, "genomes.tsv"))}
    rows = []
    for f in lines[1:]:
        name, low, high, reference = (f + ["", "", "", ""])[:4]
        if name not in by_name:
            continue
        taxid, sp = by_name[name]
        g = genomes.get(sp["representative"], {})
        reads = g.get("rdna_copies_reads", "-") or "-"
        try:
            within = "t" if float(low) <= float(reads) <= float(high) else "f"
        except ValueError:
            within = "-"
        rows.append(dict(species=name, species_taxid=taxid, accession=sp["representative"], low=low, high=high,
                         reference=reference or "-", copies_reads=reads, units_assembled=g.get("rdna_units", "-"),
                         within=within))
    return rows


def run(package, out, species, summary, log=print):
    """The checks' tables of package into out, their headline figures into summary (quality, gencode, rdna)."""
    rows = quality_rows(package, species)
    benchmarks.write(out, "quality.tsv", QUALITY_COLUMNS, rows)
    summary["quality"] = {}
    for r in rows:
        if r["supergroup"] == "all":
            level = summary["quality"].setdefault(r["level"], {"of": r["of"], "checks": {}})
            level["checks"][r["check"]] = r["flagged"]
    rows = gencode_rows(package, species)
    if any(r["agrees"] or r["differs"] or r["undecided"] for r in rows):
        benchmarks.write(out, "gencode.tsv", GENCODE_COLUMNS, rows)
        summary["gencode"] = next({k: r[k] for k in GENCODE_COLUMNS[1:]} for r in rows if r["supergroup"] == "all")
    rows = rdna_rows(package, species)
    if rows:
        benchmarks.write(out, "rdna_published.tsv", RDNA_COLUMNS, rows)
        measured = [r for r in rows if r["within"] != "-"]
        summary["rdna_published"] = {"species": len(rows), "measured": len(measured),
                                     "within": sum(r["within"] == "t" for r in measured)}
    g = summary["quality"].get("genome", {})
    log(f"checks: {g.get('checks', {}).get(ANY, 0)} of {g.get('of', 0)} genomes flagged by NCBI; "
        f"{(summary.get('gencode') or {}).get('differs', 0)} genetic codes differ")
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--out", required=True, help="the results' folder (compare.py's: its summary.json is added to)")
    a = p.parse_args(argv)
    from . import compare  # noqa: PLC0415
    compare.open_results(a.package, a.out)
    summary = compare.load_summary(a.out)
    run(a.package, a.out, compare.species_of(a.package), summary,
        log=lambda text: print(text, file=sys.stderr, flush=True))
    compare.write_summary(a.package, a.out, summary)
    json.dump({k: summary[k] for k in ("quality", "gencode", "rdna_published") if k in summary}, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
