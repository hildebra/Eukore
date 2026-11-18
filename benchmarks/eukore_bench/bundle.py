"""The benchmark step's evaluation bundle (since 2026-11-18): what it takes to judge a benchmark run away from the cluster,
as one OUT/Eukore_<version>_benchmarks_evaluation.tar.gz (the build's final job writes the like for the build:
eukore_build.evaluation_bundle):

    python3 -m eukore_bench.bundle PACKAGE RESULTS [--bench BENCH] [--out FILE] [--logs LOG:LOG]

Inside `<name>/`:
- BUNDLE.md: what is in it, what was left out or cut short, the references' versions;
- results/: the results' tables (summary.json and summary.md, reads.tsv, detection.tsv, background_hits.tsv,
  missing.tsv, ...) except those the package's own benchmarks/ held unchanged (compare.py copies them: the build's
  evaluation bundle has them); a table over BUNDLE_FILE_BYTES is cut to its first BUNDLE_LINES lines;
- results/figures/: the figures (each under FIGURE_BYTES, together under FIGURES_BYTES) with index.md;
- bench/: the references' list (compare/references.tsv) and the test genomes (read_genomes/tests.tsv);
- attracting/: the package's rows (rrna/records.tsv, genomes.tsv) of the records background_hits.tsv names, so that
  what attracts the bacterial reads can be told without the package: the records' species, source, kind, confidence
  and flags, and their genomes' tier, assembly level and NCBI's flags;
- logs/: the last LOG_BYTES of each log named (tools/eukore.sh gives the job's own).
A failure is said, never the job's end.
"""

import argparse
import collections
import contextlib
import datetime
import io
import os
import sys
import tarfile
import tempfile
import time

from . import ROOT  # noqa: F401  (the checkout on the path)
from eukore import benchmarks, package as package_mod, util

BUNDLE_FILE_BYTES = 30_000_000   # a table over this is cut to its first BUNDLE_LINES lines
BUNDLE_LINES = 20_000
FIGURE_BYTES = 8_000_000
FIGURES_BYTES = 80_000_000
LOG_BYTES = 20_000_000           # of each log, its last bytes
ATTRACTING_RECORDS = 5_000       # records of background_hits.tsv looked up (the most reads first)
FIGURE_SUFFIXES = (".png", ".svg", ".pdf", ".md")


def unchanged_own(package, name, path):
    """Whether a results file is a copy of the package's own benchmarks/ file (compare.copy_own: shutil.copy2, so the
    same size and modification time) that no step rewrote."""
    own = os.path.join(package, benchmarks.FOLDER, name)
    if not os.path.isfile(own):
        return False
    a, b = os.stat(own), os.stat(path)
    return a.st_size == b.st_size and int(a.st_mtime) == int(b.st_mtime)


def cut_head(path, lines=None):
    """A text file's first `lines` (default BUNDLE_LINES) lines -> bytes."""
    lines = BUNDLE_LINES if lines is None else lines
    out = []
    with open(path, "rb") as fh:
        for i, line in enumerate(fh):
            if i >= lines:
                break
            out.append(line)
    return b"".join(out)


def attracting(package, results, log=print):
    """The package's rows of the records and genomes background_hits.tsv names -> ({table name: (columns, rows)}),
    {} where the table is not there. A record named accession|locus is a BUSCO copy (its genome's row); the others
    are rRNA records (rrna/records.tsv's row, and the genome's it came from)."""
    hits = list(benchmarks.rows_of(os.path.join(results, "background_hits.tsv")))
    if not hits:
        return {}
    reads = collections.Counter()
    for r in hits:
        reads[r["record"]] += int(r["reads"])
    records = {rec for rec, _n in reads.most_common(ATTRACTING_RECORDS)}
    rrna, accessions = [], set()
    for rec in records:
        if rec.partition("|")[0].startswith(("GCA_", "GCF_")):
            accessions.add(rec.partition("|")[0])
    path = os.path.join(package, "rrna", "records.tsv")
    for row in benchmarks.rows_of(path) if os.path.isfile(path) else ():
        if row.get("record") in records:
            rrna.append(dict(row, reads_in_background=reads[row["record"]]))
            if row.get("accession", "-") not in ("-", ""):
                accessions.add(row["accession"])
    out = {}
    if rrna:
        out["records.tsv"] = (list(rrna[0]), sorted(rrna, key=lambda r: -r["reads_in_background"]))
    try:
        genomes = [r for r in package_mod.read_table(package, "genomes.tsv") if r["accession"] in accessions]
    except OSError:
        genomes = []
    if genomes:
        out["genomes.tsv"] = (list(genomes[0]), genomes)
    log(f"  attracting records: {len(records):,} named, {len(rrna):,} rRNA records, {len(genomes):,} genomes")
    return out


def make(package, results, bench=None, out=None, logs=(), say=print):
    """The bundle of a benchmark run (package, results: compare.py's folder, bench: prep.py's) -> its path (default
    RESULTS_evaluation.tar.gz beside RESULTS), or None where it failed."""
    results = os.path.normpath(results)
    base = os.path.basename(results) + "_evaluation"
    path = out or os.path.join(os.path.dirname(results), base + ".tar.gz")
    started = time.time()
    notes = {"left out (a copy of the package's own benchmarks/)": [], "cut to its first lines": [],
             "too large": []}
    try:
        with tempfile.TemporaryDirectory(prefix="bundle.", dir=os.path.dirname(path) or None) as tmp, \
                tarfile.open(util.part_of(path), "w:gz") as tar:
            def add_bytes(name, data):
                info = tarfile.TarInfo(f"{base}/{name}")
                info.size, info.mtime = len(data), int(time.time())
                tar.addfile(info, io.BytesIO(data))
            names = sorted(n for n in os.listdir(results)
                           if os.path.isfile(os.path.join(results, n)) and not n.startswith("."))
            for name in names:
                full = os.path.join(results, name)
                if name not in ("summary.json", "summary.md") and unchanged_own(package, name, full):
                    notes["left out (a copy of the package's own benchmarks/)"].append(name)
                elif os.path.getsize(full) > BUNDLE_FILE_BYTES and name.endswith(".tsv"):
                    add_bytes(f"results/{name}", cut_head(full))
                    notes["cut to its first lines"].append(f"{name} ({os.path.getsize(full) / 1e6:.0f} MB, "
                                                           f"{BUNDLE_LINES:,} lines kept)")
                elif os.path.getsize(full) > BUNDLE_FILE_BYTES:
                    notes["too large"].append(f"{name} ({os.path.getsize(full) / 1e6:.0f} MB)")
                else:
                    tar.add(full, f"{base}/results/{name}")
            figures = os.path.join(results, "figures")
            total = 0
            for name in sorted(os.listdir(figures)) if os.path.isdir(figures) else ():
                full = os.path.join(figures, name)
                size = os.path.getsize(full) if os.path.isfile(full) else None
                if size is None or not name.endswith(FIGURE_SUFFIXES):
                    continue
                if size > FIGURE_BYTES or total + size > FIGURES_BYTES:
                    notes["too large"].append(f"figures/{name} ({size / 1e6:.0f} MB)")
                    continue
                total += size
                tar.add(full, f"{base}/results/figures/{name}")
            for sub in (os.path.join("compare", "references.tsv"), os.path.join("read_genomes", "tests.tsv")):
                if bench and os.path.isfile(os.path.join(bench, sub)):
                    tar.add(os.path.join(bench, sub), f"{base}/bench/{os.path.basename(sub)}")
            for name, (columns, rows) in attracting(package, results, say).items():
                benchmarks.write(tmp, name, columns, rows)
                tar.add(os.path.join(tmp, name), f"{base}/attracting/{name}")
            for log in filter(None, logs):
                if os.path.isfile(log):
                    with open(log, "rb") as fh:
                        fh.seek(max(0, os.path.getsize(log) - LOG_BYTES))
                        add_bytes(f"logs/{os.path.basename(log)}", fh.read())
            add_bytes("BUNDLE.md", description(package, results, bench, base, notes).encode())
        os.replace(util.part_of(path), path)
    except OSError as e:
        say(f"WARNING: no benchmark bundle ({e})")
        with contextlib.suppress(OSError):
            os.remove(util.part_of(path))
        return None
    say(f"benchmark bundle: {path} ({os.path.getsize(path) / 1e6:.1f} MB, {time.time() - started:.0f} s): the results' "
        "tables and figures, what attracts the background reads, the logs; copy it into local/ of the checkout to "
        "have it looked at")
    return path


def description(package, results, bench, base, notes):
    """BUNDLE.md's text."""
    lines = [f"# {base}", "", f"Made {datetime.datetime.now().isoformat(timespec='seconds')} from the package "
             f"`{package}` and the results `{results}`.", "",
             "- `results/`: the results' tables (`summary.json` and `summary.md` first), `results/figures/`: the figures "
             "and `index.md`;",
             "- `bench/`: the other references (`references.tsv`) and the read benchmark's test genomes (`tests.tsv`);",
             "- `attracting/`: the package's rows of the records `background_hits.tsv` names (`records.tsv`: rRNA "
             "records, `genomes.tsv`: their genomes), so that what attracts the bacterial reads can be told;",
             "- `logs/`: the jobs' logs, the last bytes of each.", ""]
    for what, items in notes.items():
        if items:
            lines += [f"Not in the bundle, {what}:", ""] + [f"- {x}" for x in items] + [""]
    if bench:
        refs = os.path.join(bench, "compare", "references.tsv")
        if os.path.isfile(refs):
            with open(refs) as fh:
                lines += ["The other references (`bench/references.tsv`):", "", "```", fh.read().rstrip(), "```", ""]
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("package", help="the Eukore package (OUT/Eukore_<version>)")
    p.add_argument("results", help="compare.py's results folder (OUT/Eukore_<version>_benchmarks)")
    p.add_argument("--bench", help="prep.py's folder (BENCH): its references' list and the test genomes")
    p.add_argument("--out", help="the bundle's path (default: RESULTS_evaluation.tar.gz beside RESULTS)")
    p.add_argument("--logs", default="", help="logs to add, colon-separated (the last bytes of each)")
    a = p.parse_args(argv)
    if not os.path.isdir(a.results):
        p.error(f"{a.results}: no results (eukore-benchmark.sh makes them)")
    make(a.package, a.results, a.bench, a.out, a.logs.split(":"), say=lambda t: print(t, file=sys.stderr, flush=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
