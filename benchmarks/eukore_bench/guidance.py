"""Which identity parts taxa, per marker and rank (since 2026-11-12; compare.py runs it; benchmarks/README.md): the
numbers a user of the package chooses a marker and an identity cutoff by.

thresholds.tsv, per group of species (`supergroup`: `all`, and since 2026-11-15 each group, for the markers but each
BUSCO's own), marker (all BUSCOs together, each BUSCO, the 18S, ...) and rank (species, genus, family, order):
- `within`: sequences that should stay together at that rank, and `between`: those that should be apart;
- `threshold`: the identity (%) that parts the two best (the highest of those with the best balanced accuracy);
- `within_kept` and `between_apart`: the shares on the right side of it; `accuracy`: their mean. An accuracy near
  0.5 is a coin's: the marker does not part that rank's taxa at any identity (USABLE), whatever threshold is named;
- `nearest`, `nearest_median`, `nearest_under_99` (since 2026-11-15): how alike the nearest other species is, with
  no threshold in between. A species' row: every copy's nearest copy of another species; a genus', family's or
  order's: the copies whose nearest other species shares that rank and no deeper one (a genus: the nearest
  congener). The median is a step of half a point (`99.5`: 99.5 to under 100; `<70`: no copy 70% alike).

A species: `within` is the intraspecific panel's species (intraspecific/barcode_gap.tsv: each species' least alike
pair of its own assemblies), `between` every copy's nearest copy of another species (the package's
benchmarks/resolution_hist.tsv). A genus, family or order: `within` the copies whose nearest other species shares
that rank and no deeper one, `between` those whose nearest shares only a higher rank or none. A package built before
2026-11-12 has no resolution_hist.tsv: nothing is written. One built before 2026-11-15 has no groups in it: the rows
of `all` alone. A group's species row has a threshold only where the panel holds MIN_VALUES species of the group;
its `nearest` columns are given without.

Read on the first full build (2026-11-15): the 18S' species row named 97% with an accuracy of 0.56, and the page
showed it as a cutoff, though 75% of the species' 18S have another species' within 97% (and 60% within 99%).

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

COLUMNS = ("supergroup", "marker", "rank", "threshold", "within", "within_kept", "between", "between_apart",
           "accuracy", "nearest", "nearest_median", "nearest_under_99", "source")
ALL = "all"
USABLE = 0.65    # a balanced accuracy under it: the marker does not part that rank's taxa (the page and summary say so)
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
    """resolution_hist.tsv -> {group: {marker: {shared rank: {identity or None: copies}}}}: ALL (every row's), and
    each group of its `supergroup` column (a table without one, or `-`: no group's)."""
    out = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(collections.Counter)))
    for r in benchmarks.rows_of(path):
        value = None if r["identity"] == "-" else float(r["identity"])
        for group in {ALL, r.get("supergroup") or "-"} - {"-"}:
            out[group][r["marker"]][r["shared_rank"]][value] += int(r["copies"])
    return out


def read_gap(path):
    """intraspecific/barcode_gap.tsv -> {group: {marker: {intra_min: species}}}: each panel species' least alike pair
    of its own assemblies, of ALL and of its `group`; the BUSCOs also together (benchmarks.ALL_BUSCOS)."""
    out = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for r in benchmarks.rows_of(path):
        try:
            value = float(r["intra_min"])
        except (KeyError, ValueError):
            continue
        for group in {ALL, r.get("group") or "-"} - {"-"}:
            out[group][r["marker"]][value] += 1
            if r["marker"].endswith("at2759"):
                out[group][benchmarks.ALL_BUSCOS][value] += 1
    return out


def nearest_figures(values):
    """{identity or None: copies} -> (copies, their median's step as text, the share under 99%)."""
    n = sum(values.values())
    if not n:
        return 0, "-", "-"
    seen, median = 0, None
    for v, c in sorted(values.items(), key=lambda kv: (kv[0] is not None, kv[0] or 0.0)):  # None (under FLOOR) first
        seen += c
        if seen * 2 >= n:
            median = v
            break
    under = sum(c for v, c in values.items() if v is None or v < 99.0) / n
    return n, f"<{FLOOR:g}" if median is None else f"{median:.1f}", f"{under:.3f}"


def threshold_rows(hist, gap):
    """thresholds.tsv's rows of the package's own markers: of ALL species, then of each group (the pooled markers
    alone: each BUSCO's own per group would be 129 rows a group and rank)."""
    rows = []
    for group in sorted(hist, key=lambda g: (g != ALL, g)):
        for marker in sorted(hist[group], key=lambda m: (m != benchmarks.ALL_BUSCOS, m != "18S", m)):
            if group != ALL and marker.endswith("at2759"):
                continue
            by_rank = hist[group][marker]
            for rank in RANKS:
                if rank == "species":
                    within = gap.get(group, {}).get(marker, {})
                    between = collections.Counter()
                    for shared, values in by_rank.items():
                        between.update(values)
                    described = between
                else:
                    within = by_rank.get(rank, {})
                    between = collections.Counter()
                    for shared in ABOVE[rank] + ("unmatched",):
                        between.update(by_rank.get(shared, {}))
                    described = within
                n_in, n_out = sum(within.values()), sum(between.values())
                if n_out < MIN_VALUES or (n_in < MIN_VALUES and (group == ALL or rank != "species")):
                    continue  # too few; a group's species row stays, its panel species too few for a threshold
                row = dict(supergroup=group, marker=marker, rank=rank, threshold="-", within=n_in, within_kept="-",
                           between=n_out, between_apart="-", accuracy="-", source="package")
                if n_in >= MIN_VALUES:
                    t, kept, apart, accuracy = best_threshold(within, between)
                    row.update(threshold=f"{t:.1f}", within_kept=f"{kept:.3f}", between_apart=f"{apart:.3f}",
                               accuracy=f"{accuracy:.3f}")
                row["nearest"], row["nearest_median"], row["nearest_under_99"] = nearest_figures(described)
                rows.append(row)
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
        rows.append(dict(supergroup=ALL, marker=marker, rank=rank, threshold=threshold, within="-", within_kept="-",
                         between="-", between_apart="-", accuracy="-", nearest="-", nearest_median="-",
                         nearest_under_99="-", source=f"published: {taxa}; {reference}"))
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
        if r["supergroup"] == ALL and r["marker"] in (benchmarks.ALL_BUSCOS, "18S") and r["threshold"] != "-":
            summary["thresholds"].setdefault(r["marker"], {})[r["rank"]] = {
                "threshold": float(r["threshold"]), "accuracy": float(r["accuracy"])}
    groups = sorted({r["supergroup"] for r in rows} - {ALL})
    log("thresholds: " + "; ".join(
        f"{m} " + ", ".join(f"{rank} {v['threshold']:g}%" + ("" if v["accuracy"] >= USABLE else " (no usable cutoff)")
                            for rank, v in by.items()) for m, by in summary["thresholds"].items())
        + (f"; {len(groups)} groups' too" if groups else "; no groups (a resolution_hist.tsv before 2026-11-15)"))
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
