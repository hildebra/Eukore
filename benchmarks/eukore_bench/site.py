"""The benchmarks' results for Eukore's website (since 2026-11-12): the build's release (OUT/release, eukore/release.py)
holds what the build knows of its package; the benchmarks run after it, so their tables reach the site by this step:

    python3 -m eukore_bench.site RESULTS --release OUT/release [--figures RESULTS/figures]

It writes OUT/release/data/benchmarks.json (the website's Benchmarks page reads it: website/pages/benchmarks.php) and
copies the figures the page shows into OUT/release/data/figures/. What the page shows is what a user of the database
looks up: what each database beside Eukore is and holds, the breadth at every rank, the markers one species links,
the named lists' organisms, the checks' shares, the same species' markers in two databases, the identity that parts
taxa, the tree's agreement with NCBI's taxonomy, the changes since the release before. Nothing of the build's own
files there (stats.json, taxonomy.json, downloads.json) is touched.
"""

import argparse
import datetime
import json
import os
import shutil
import sys

from . import figures as figures_mod

FILE = "benchmarks.json"
FIGURES = ("breadth_ranks", "content", "upset", "novel", "agreement", "quality", "thresholds", "tree", "growth",
           "changes")
TABLES = ("databases", "breadth_ranks", "content", "panels", "novel", "agreement", "thresholds", "gencode",
          "changes", "sizes")
ALL_ROWS = ("databases", "panels", "agreement", "sizes")  # tables without a group column, or shown whole
SHOWN_MARKERS = ("BUSCO (all)", "18S")  # thresholds.tsv's markers shown (each BUSCO's own stay in the table)


def rows_of(results, name):
    """A results table's rows for the page: those of `all` (the groups' stay in the table itself), each whole where
    the table has no groups; thresholds.tsv's without each BUSCO's own."""
    rows = figures_mod.read(os.path.join(results, f"{name}.tsv"))
    if name == "thresholds":
        return [r for r in rows if r["marker"] in SHOWN_MARKERS or r["source"].startswith("published")]
    if name in ALL_ROWS:
        return rows
    return [r for r in rows if r.get("supergroup") == "all"]


def data(results, figures=None):
    """benchmarks.json's content of a results folder (compare.py's), and the figures to copy [(name, path)]."""
    summary = figures_mod.load_summary(results)
    out = {"package": None,  # write(): the one given, else the release's own
           "made": datetime.datetime.now().isoformat(timespec="seconds"),
           "tables": {name: rows for name in TABLES for rows in [rows_of(results, name)] if rows},
           "summary": {k: summary[k] for k in ("upset", "tree", "quality", "rdna_published", "growth", "cost",
                                               "changes", "breadth_ranks") if k in summary},
           "left_out": summary.get("left_out") or {}, "figures": []}
    captions = {name: caption for name, _fn, caption, _tables in figures_mod.FIGURES}
    found = []
    for name in FIGURES:
        path = os.path.join(figures or os.path.join(results, figures_mod.FOLDER), f"{name}.png")
        if os.path.isfile(path):
            out["figures"].append({"name": name, "file": f"figures/{name}.png", "caption": captions.get(name, "")})
            found.append((name, path))
    return out, found


def write(results, release, figures=None, package=None, log=print):
    """benchmarks.json and the figures into release/data -> the JSON's path."""
    folder = os.path.join(release, "data")
    if not os.path.isdir(folder):
        raise SystemExit(f"{release}: no data/ in it (a build's OUT/release: eukore_build.py's release step)")
    out, found = data(results, figures)
    if package:
        out["package"] = os.path.basename(os.path.normpath(package))
    elif not out["package"]:  # the release's own package (the build's stats.json)
        try:
            with open(os.path.join(folder, "stats.json")) as fh:
                out["package"] = json.load(fh).get("package")
        except (OSError, ValueError):
            pass
    if found:
        os.makedirs(os.path.join(folder, "figures"), exist_ok=True)
        for name, path in found:
            shutil.copyfile(path, os.path.join(folder, "figures", f"{name}.png"))
    path = os.path.join(folder, FILE)
    with open(path + ".part", "w", newline="\n") as fh:
        json.dump(out, fh, indent=1)
    os.replace(path + ".part", path)
    log(f"site: {path} ({len(out['tables'])} tables, {len(found)} figures)")
    return path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("results", help="the benchmarks' results (eukore-benchmark.sh's RESULTS)")
    p.add_argument("--release", required=True, help="the build's release folder (OUT/release)")
    p.add_argument("--figures", help="the figures' folder (default RESULTS/figures: eukore-figures.sh's, as PNG)")
    p.add_argument("--package", help="the package's folder or name (default: the release's own)")
    a = p.parse_args(argv)
    write(a.results, a.release, a.figures, a.package, log=lambda text: print(text, file=sys.stderr))
    return 0


if __name__ == "__main__":
    sys.exit(main())
