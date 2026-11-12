"""Which identity parts taxa, per marker and rank (since 2026-11-12; compare.py runs it; benchmarks/README.md): the
numbers a user of the package chooses a marker and an identity cutoff by.

thresholds.tsv, per marker (all BUSCOs together, each BUSCO, the 18S) and rank (species, genus, family, order):
- `within`: sequences that should stay together at that rank, and `between`: those that should be apart;
- `threshold`: the identity (%) that parts the two best (the highest of those with the best balanced accuracy);
- `within_kept` and `between_apart`: the shares on the right side of it; `accuracy`: their mean.

A species: `within` is the intraspecific panel's species (intraspecific/barcode_gap.tsv: each species' least alike
pair of its own assemblies), `between` every copy's nearest copy of another species (the package's
benchmarks/resolution_hist.tsv). A genus, family or order: `within` the copies whose nearest other species shares
that rank and no deeper one, `between` those whose nearest shares only a higher rank or none. A package built before
2026-11-12 has no resolution_hist.tsv: nothing is written.

The published thresholds of benchmarks/thresholds_published.tsv follow as rows of `source` published.

    python3 -m eukore_bench.guidance --package OUT/Eukore_2026.11 --out RESULTS
"""

import argparse
import collections
import json
import os
import sys

from . import ROOT
from eukore import benchmarks

COLUMNS = ("marker", "rank", "threshold", "within", "within_kept", "between", "between_apart", "accuracy", "source")
PUBLISHED = os.path.join(ROOT, "benchmarks", "thresholds_published.tsv")
RANKS = ("species", "genus", "family", "order")
ABOVE = {"genus": ("family", "order", "class", "phylum", "none"), "family": ("order", "class", "phylum", "none"),
         "order": ("class", "phylum", "none")}
STEP = 0.5
FLOOR = 70.0     # the lowest identity resolution.tsv's search reports: a copy without a hit is under it
MIN_VALUES = 5   # sequences on each side, at least, for a threshold to be said


def best_threshold(within, between):
    """within, between: {identity: count} (a between of None: under FLOOR) -> (threshold, within kept, between
    apart, balanced accuracy): the identity T at which the most of within are T or more and of between under T; of
    several as good the highest."""
    n_in, n_out = sum(within.values()), sum(between.values())
    best = None
    t = 100.0
    while t >= FLOOR:
        kept = sum(c for v, c in within.items() if v is not None and v >= t) / n_in
        apart = sum(c for v, c in between.items() if v is None or v < t) / n_out
        if best is None or (kept + apart) / 2 > best[3] + 1e-9:
            best = (t, kept, apart, (kept + apart) / 2)
        t -= STEP
    return best


def read_hist(path):
    """resolution_hist.tsv -> {marker: {shared rank: {identity or None: copies}}}."""
    out = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for r in benchmarks.rows_of(path):
        value = None if r["identity"] == "-" else float(r["identity"])
        out[r["marker"]][r["shared_rank"]][value] += int(r["copies"])
    return out


def read_gap(path):
    """intraspecific/barcode_gap.tsv -> {marker: {intra_min: species}}: each panel species' least alike pair of its
    own assemblies; the BUSCOs also together (benchmarks.ALL_BUSCOS)."""
    out = collections.defaultdict(collections.Counter)
    for r in benchmarks.rows_of(path):
        try:
            value = float(r["intra_min"])
        except (KeyError, ValueError):
            continue
        out[r["marker"]][value] += 1
        if r["marker"].endswith("at2759"):
            out[benchmarks.ALL_BUSCOS][value] += 1
    return out


def threshold_rows(hist, gap):
    """thresholds.tsv's rows of the package's own markers."""
    rows = []
    for marker in sorted(hist, key=lambda m: (m != benchmarks.ALL_BUSCOS, m != "18S", m)):
        by_rank = hist[marker]
        for rank in RANKS:
            if rank == "species":
                within = gap.get(marker, {})
                between = collections.Counter()
                for shared, values in by_rank.items():
                    between.update(values)
            else:
                within = by_rank.get(rank, {})
                between = collections.Counter()
                for shared in ABOVE[rank] + ("unmatched",):
                    between.update(by_rank.get(shared, {}))
            if sum(within.values()) < MIN_VALUES or sum(between.values()) < MIN_VALUES:
                continue
            t, kept, apart, accuracy = best_threshold(within, between)
            rows.append(dict(marker=marker, rank=rank, threshold=f"{t:.1f}", within=sum(within.values()),
                             within_kept=f"{kept:.3f}", between=sum(between.values()), between_apart=f"{apart:.3f}",
                             accuracy=f"{accuracy:.3f}", source="package"))
    return rows


def published_rows(path=PUBLISHED):
    """benchmarks/thresholds_published.tsv's thresholds as thresholds.tsv's rows (`source`: published, its taxa and
    reference)."""
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    rows = []
    for f in lines[1:]:
        marker, rank, threshold, taxa, reference = (f + [""] * 5)[:5]
        rows.append(dict(marker=marker, rank=rank, threshold=threshold, within="-", within_kept="-", between="-",
                         between_apart="-", accuracy="-", source=f"published: {taxa}; {reference}"))
    return rows


def run(package, out, summary, log=print):
    """thresholds.tsv of package into out, its headline figures (the BUSCOs together and the 18S) into summary."""
    hist = read_hist(os.path.join(out, "resolution_hist.tsv"))
    if not hist:
        summary.setdefault("left_out", {})["thresholds"] = \
            "no resolution_hist.tsv (a package built before 2026-11-12, or without vsearch)"
        return summary
    rows = threshold_rows(hist, read_gap(os.path.join(package, "intraspecific", "barcode_gap.tsv")))
    benchmarks.write(out, "thresholds.tsv", COLUMNS, rows + published_rows())
    summary["thresholds"] = {}
    for r in rows:
        if r["marker"] in (benchmarks.ALL_BUSCOS, "18S"):
            summary["thresholds"].setdefault(r["marker"], {})[r["rank"]] = {
                "threshold": float(r["threshold"]), "accuracy": float(r["accuracy"])}
    log("thresholds: " + "; ".join(f"{m} " + ", ".join(f"{rank} {v['threshold']:g}%" for rank, v in by.items())
                                   for m, by in summary["thresholds"].items()))
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--out", required=True, help="the results' folder (compare.py's: its summary.json is added to)")
    a = p.parse_args(argv)
    from . import compare  # noqa: PLC0415
    compare.open_results(a.package, a.out)
    summary = compare.load_summary(a.out)
    run(a.package, a.out, summary, log=lambda text: print(text, file=sys.stderr, flush=True))
    compare.write_summary(a.package, a.out, summary)
    json.dump(summary.get("thresholds") or {"left_out": summary["left_out"].get("thresholds")}, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
