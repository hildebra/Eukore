"""A finished Eukore package against other references, and the read benchmark (since 2026-11-12; until then the
build's package stage, with --read_compare and --read_genomes): into an output folder of its own, never the
package's (its manifest checks the package's files):

    python3 -m eukore_bench.compare --package OUT/Eukore_2026.11 --bench BENCH --out OUT_BENCH [-t 32]

OUT_BENCH gets the package's own benchmarks (its benchmarks/ tables, summary.json and summary.md, copied), then:
- breadth.tsv: the build's (NCBI's species with an assembly, Eukore's) with the first other reference's columns
  (BENCH/compare/references.tsv's first: EukDetect2 where there);
- overlap.tsv: the build's rows (PR2's labels, EUKARYOME, KSGP, ROD, BUSCO's dataset) and one per other reference;
- breadth_ranks.tsv, content.tsv, panels.tsv: the package as a resource (resource.py): its breadth at every rank
  against NCBI's, the markers its species hold from one assembly, the named lists' organisms;
- databases.tsv, upset.tsv, novel.tsv, agreement.tsv: the package beside the other databases (databases.py): the
  comparison table, the species several databases share, the species a genome is the first sequence of a marker
  for, the same species' markers in two databases;
- quality.tsv, gencode.tsv, rdna_published.tsv: how good its records are (checks.py): each check's share per
  group, the genetic codes against NCBI's, the rDNA copies by reads against published numbers;
- thresholds.tsv: the identity that parts taxa per marker and rank (guidance.py);
- sizes.tsv, growth.tsv, changes.tsv: the release among releases (history.py): its parts' sizes, its genomes by the
  year NCBI released them, and with --previous the species and genomes changed since an earlier package;
- tree.nwk, tree_tips.tsv, tree_taxa.tsv, tree_strays.tsv: the tree of its BUSCO proteins against NCBI's taxonomy
  (tree.py; where BENCH holds BUSCO's dataset and hmmalign and FastTree are on the path);
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
from . import checks, databases, guidance, history, prep, reads, resource, tree as tree_mod
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


def copy_own(package, out):
    """The package's own benchmarks (its benchmarks/ tables, summary.json and summary.md; its figures left out)
    copied into out."""
    os.makedirs(out, exist_ok=True)
    for path in glob.glob(os.path.join(package, benchmarks.FOLDER, "*")):
        if os.path.isfile(path):
            shutil.copy2(path, out)


def open_results(package, out):
    """A results folder for one step run alone (resource.py, databases.py, ..): the package's own benchmarks copied
    into it where it is new; one that holds a summary.json already is kept as it is, so that a step adds to what the
    others wrote."""
    if not os.path.isfile(os.path.join(out, "summary.json")):
        copy_own(package, out)


def load_summary(out):
    """out's summary.json (the package's own, copied), or a bare one."""
    try:
        with open(os.path.join(out, "summary.json")) as fh:
            summary = json.load(fh)
    except (OSError, ValueError):
        summary = {"species": 0, "busco_loci": 0, "left_out": {}}
    summary.setdefault("left_out", {})
    return summary


def write_summary(package, out, summary):
    """summary.json and summary.md of out."""
    with open(os.path.join(out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    composition = list(benchmarks.rows_of(os.path.join(out, "composition.tsv")))
    with open(os.path.join(out, "summary.md"), "w", newline="\n") as fh:
        fh.write(benchmarks.summary_md(summary, os.path.basename(os.path.normpath(package)), None,
                                       composition or None))


def run(package, bench, out, threads=1, minimap2=None, log=print, tmp=None, taxdump=None, tree=True,
        previous=None,
        align_cache=None):
    """The comparisons and the read benchmark of package (BENCH: prep.py's folder) into out -> summary.json's
    figures. taxdump: NCBI's whole taxonomy (default BENCH's download; else the package's taxonomy/, which holds
    its own lineages alone: a species another reference holds alone is then not counted)."""
    copy_own(package, out)
    summary = load_summary(out)
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
    # the package as a resource (since 2026-11-12): its breadth at every rank, its linked markers, the named lists
    resource.run(package, out, ncbi, species, summary, bool(whole), log=log)
    minimap2 = minimap2 or shutil.which("minimap2")
    # beside the other databases (since 2026-11-12): the comparison table, the species several share, the species a
    # genome is the first sequence of a marker for, the same species' markers in two databases
    databases.run(package, bench, out, ncbi, species, summary, theirs, compare, minimap2, threads, tmp, log)
    # how good its records are (since 2026-11-12): each check's share per group, the genetic codes, the rDNA copies;
    # the tree of its BUSCO proteins against NCBI's taxonomy (where BUSCO's profiles, hmmalign and FastTree are)
    checks.run(package, out, species, summary, log)
    guidance.run(package, out, summary, log)  # the identity that parts taxa, per marker and rank
    history.run(package, out, species, summary, previous, log)  # its sizes, the genomes by year, the changes
    if tree:
        try:
            tree_mod.run(package, out, ncbi, species, summary, prep.odb_of(bench), threads=threads,
                         cache=align_cache, tmp=tmp, log=log)
        except (RuntimeError, OSError) as e:  # a tool that failed: the other benchmarks are not lost
            summary["left_out"]["tree"] = f"failed ({e})"
            log(f"WARNING: the tree failed ({e}): left out")
    else:
        summary["left_out"]["tree"] = "not asked (--no_tree)"
    tests_file = os.path.join(bench, "read_genomes", "tests.tsv")
    tests = [r["accession"] for r in benchmarks.rows_of(tests_file)] if os.path.isfile(tests_file) else []
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
    write_summary(package, out, summary)
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--bench", required=True, help="prep.py's folder (BENCH)")
    p.add_argument("--out", required=True, help="the results' folder (not the package's)")
    p.add_argument("--minimap2", default=None, help="minimap2 (default: on PATH)")
    p.add_argument("--taxdump", help="NCBI's whole taxdump (default: BENCH/downloads/taxdump.tar.gz)")
    p.add_argument("--tmp", help="the read benchmark's temporary files (a node's local disk; default: in --out)")
    p.add_argument("--previous", help="an earlier package: the changes from it (history.py's changes.tsv)")
    p.add_argument("--no_tree", action="store_true", help="no species tree (tree.py: hours for all eukaryotes)")
    p.add_argument("--align_cache", help="the build's WORK/align: its proteins' rows on BUSCO's profiles, for the tree")
    p.add_argument("-t", "--threads", type=int, default=8)
    a = p.parse_args(argv)
    summary = run(a.package, a.bench, a.out, a.threads, a.minimap2,
                  log=lambda text: print(text, file=sys.stderr, flush=True), tmp=a.tmp, taxdump=a.taxdump,
                  tree=not a.no_tree, align_cache=a.align_cache, previous=a.previous)
    json.dump({k: summary[k] for k in ("breadth", "overlap", "reads", "left_out") if k in summary}, sys.stdout,
              indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
