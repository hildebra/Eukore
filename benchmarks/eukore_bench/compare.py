"""A finished Eukore package against other references, and the read benchmark (since 2026-11-12; until then the
build's package stage, with --read_compare and --read_genomes): into an output folder of its own, never the
package's (its manifest checks the package's files):

    python3 -m eukore_bench.compare --package OUT/Eukore_2026.11 --bench BENCH --out OUT_BENCH [-t 32]

OUT_BENCH gets the package's own benchmarks (its benchmarks/ tables, summary.json and summary.md, copied), then:
- breadth.tsv: the build's (NCBI's species with an assembly, Eukore's) with the first other reference's columns
  (BENCH/compare/references.tsv's first: EukDetect2 where there);
- overlap.tsv: the build's rows (PR2's labels, EUKARYOME, KSGP, ROD, BUSCO's dataset) and one per other reference;
- reads.tsv, detection.tsv: the read benchmark (reads.py) of Eukore and every other reference, on the test genomes
  prep.py placed in BENCH/read_genomes and its bacterial background (BENCH/read_background);
- summary.json and summary.md with their figures added (breadth, overlap, reads, detection).
figures.py then draws every figure from OUT_BENCH. benchmarks/eukore-benchmark.sh runs this in the benchmarks'
environment.
"""

import argparse
import collections
import glob
import json
import os
import shutil
import sys
import tempfile

from . import ROOT  # noqa: F401  (the checkout on the path)
from . import prep, reads
from eukore import benchmarks, ncbi as ncbi_mod, package as package_mod


def read_summary(rows):
    """reads.tsv's headline: per reference and mode the marker reads' precision and recall at species and genus (held
    out: of the test genomes with another species of their genus in the package, the only ones a held-out read can be
    placed right to its genus), and the background reads assigned per million, and to a phylum or deeper."""
    out = collections.defaultdict(dict)
    sums = collections.defaultdict(collections.Counter)
    for r in rows:
        if r["mode"] == "held_out" and not int(r["congeners"]):
            continue
        key = (r["reference"], r["mode"], r["kind"], r["reads"])
        for c in ("n", "assigned", "species_correct", "species_wrong", "genus_correct", "genus_wrong",
                  "phylum_wrong"):
            sums[key][c] += int(r[c])
        sums[key]["censored"] += int(r.get("censored") or 0)  # since 2026-11-10
    for (ref, mode, kind, which), s in sums.items():
        if kind == "eukaryote" and which == "genome" and mode != "in" and s["n"]:
            # since 2026-11-10: the share of the genome reads whose secondary hits were all to records left out
            out[ref].setdefault("censored", {})[mode] = round(s["censored"] / s["n"], 4)
        if kind == "eukaryote" and which == "marker":
            out[ref][mode] = {f"{r}_{m}": round(v, 4) if v is not None else None for r in ("species", "genus")
                              for m, v in (("recall", s[f"{r}_correct"] / s["n"] if s["n"] else None),
                                           ("precision", s[f"{r}_correct"] / (s[f"{r}_correct"] + s[f"{r}_wrong"])
                                            if s[f"{r}_correct"] + s[f"{r}_wrong"] else None))}
        elif kind == "background":
            out[ref]["background_per_million"] = round(1e6 * s["assigned"] / s["n"], 1) if s["n"] else None
            out[ref]["background_to_a_phylum_per_million"] = round(1e6 * s["phylum_wrong"] / s["n"], 1) \
                if s["n"] else None
    return dict(out)


def species_of(package):
    """The package's species.tsv: {taxid: row, its supergroup set (format 5's group)}."""
    return {r["taxid"]: dict(r, supergroup=r.get("group") or r.get("supergroup", "-"))
            for r in package_mod.read_table(package, "species.tsv")}


def breadth_with(rows, species, theirs, ncbi):
    """The build's breadth.tsv rows with the compared reference's columns (theirs: its NCBI species)."""
    cache = {}
    mine = {int(t): r["supergroup"] for t, r in species.items() if r["representative"] != "-"}
    their = {sp: benchmarks.supergroup_of(ncbi, sp, cache) for sp in theirs}
    out = []
    groups = [r["supergroup"] for r in rows]
    for g in sorted(set(their.values()) - set(groups) - {"all"}):  # a supergroup the build had none of
        rows = rows[:-1] + [dict(supergroup=g, ncbi_species=0, eukore=0, eukore_rrna_only=0)] + rows[-1:]
    for r in rows:
        g = r["supergroup"]
        a = {sp for sp, x in mine.items() if g in (x, "all")}
        b = {sp for sp, x in their.items() if g in (x, "all")}
        out.append(dict(r, compared=len(b), both=len(a & b), eukore_only=len(a - b), compared_only=len(b - a)))
    return out


def run(package, bench, out, threads=1, minimap2=None, log=print, tmp=None, taxdump=None):
    """The comparisons and the read benchmark of package (BENCH: prep.py's folder) into out -> summary.json's
    figures. taxdump: NCBI's whole taxonomy (default BENCH's download; else the package's taxonomy/, which holds
    its own lineages alone: a species another reference holds alone is then not counted)."""
    os.makedirs(out, exist_ok=True)
    own = os.path.join(package, benchmarks.FOLDER)
    for path in glob.glob(os.path.join(own, "*")):  # the package's own benchmarks, its figures left out
        if os.path.isfile(path):
            shutil.copy2(path, out)
    try:
        with open(os.path.join(out, "summary.json")) as fh:
            summary = json.load(fh)
    except (OSError, ValueError):
        summary = {"species": 0, "busco_loci": 0, "left_out": {}}
    summary.setdefault("left_out", {})
    whole = prep.taxdump_of(bench, taxdump)
    if not whole:
        log("WARNING: no NCBI taxdump in BENCH (prep.py's taxdump step) or given: the package's taxonomy, the "
            "other references' species outside its lineages not counted")
    ncbi = ncbi_mod.Ncbi(whole or os.path.join(package, "taxonomy"))
    species = species_of(package)
    compare = prep.read_references(bench)
    if not compare:
        log(f"WARNING: no other reference in {os.path.join(bench, 'compare', prep.REFERENCES)} (prep.py): the read "
            "benchmark of Eukore alone, no comparisons")
    theirs = {name: benchmarks.compared_species(table, ncbi) for name, _fasta, table in compare}
    if theirs:
        first = compare[0][0]
        rows = list(benchmarks.rows_of(os.path.join(out, "breadth.tsv")))
        if rows:
            rows = breadth_with(rows, species, theirs[first], ncbi)
            benchmarks.write(out, "breadth.tsv", benchmarks.BREADTH_COLUMNS, rows)
            summary["breadth"] = dict(next(({k: v for k, v in r.items() if k != "supergroup"} for r in rows
                                            if r["supergroup"] == "all"), {}), compared_name=first)
        added = benchmarks.overlap_rows(species, theirs, ncbi)
        kept = [r for r in benchmarks.rows_of(os.path.join(out, "overlap.tsv")) if r["database"] not in theirs]
        rows = kept + added
        benchmarks.write(out, "overlap.tsv", benchmarks.OVERLAP_COLUMNS, rows)
        summary["overlap"] = {r["database"]: {k: int(r[k]) for k in ("species", "both", "eukore_only",
                                                                      "database_only")}
                              for r in rows if r["supergroup"] == "all"}
        log(f"compared: {', '.join(f'{n} {len(s):,} species' for n, s in theirs.items())}")
    tests_file = os.path.join(bench, "read_genomes", "tests.tsv")
    tests = [r["accession"] for r in benchmarks.rows_of(tests_file)] if os.path.isfile(tests_file) else []
    minimap2 = minimap2 or shutil.which("minimap2")
    summary["left_out"].pop("reads", None)
    if tests and minimap2:
        work = tempfile.mkdtemp(prefix="read_benchmark.", dir=tmp or out)
        detections = []
        try:
            rows = reads.run(package, ncbi, os.path.join(bench, "read_genomes"), tests, minimap2, threads, work,
                             compare, reads.read_background(prep.READ_BACKGROUND,
                                                            os.path.join(bench, "read_background")),
                             None, log, detections=detections)
            benchmarks.write(out, "reads.tsv", reads.COLUMNS, rows)
            benchmarks.write(out, "detection.tsv", reads.DETECT_COLUMNS, detections)
            summary["reads"] = read_summary(rows)
            summary["detection"] = reads.detection_summary(detections)
        finally:
            shutil.rmtree(work, ignore_errors=True)
    else:
        summary["left_out"]["reads"] = ("no test genomes (prep.py --package)" if not tests
                                        else "no minimap2 (the benchmarks' environment)")
    with open(os.path.join(out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    composition = list(benchmarks.rows_of(os.path.join(out, "composition.tsv")))
    with open(os.path.join(out, "summary.md"), "w", newline="\n") as fh:
        fh.write(benchmarks.summary_md(summary, os.path.basename(os.path.normpath(package)), None,
                                       composition or None))
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--bench", required=True, help="prep.py's folder (BENCH)")
    p.add_argument("--out", required=True, help="the results' folder (not the package's)")
    p.add_argument("--minimap2", default=None, help="minimap2 (default: on PATH)")
    p.add_argument("--taxdump", help="NCBI's whole taxdump (default: BENCH/downloads/taxdump.tar.gz)")
    p.add_argument("--tmp", help="the read benchmark's temporary files (a node's local disk; default: in --out)")
    p.add_argument("-t", "--threads", type=int, default=8)
    a = p.parse_args(argv)
    summary = run(a.package, a.bench, a.out, a.threads, a.minimap2,
                  log=lambda text: print(text, file=sys.stderr, flush=True), tmp=a.tmp, taxdump=a.taxdump)
    json.dump({k: summary[k] for k in ("breadth", "overlap", "reads", "left_out") if k in summary}, sys.stdout,
              indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
