#!/usr/bin/env python3
"""Draw a Eukore package's benchmarks (eukore/benchmarks.py; docs/benchmarks.md) as figures, and their headline
figures as a Markdown table.

    python3 docs/benchmarks/plot_benchmarks.py OUT/Eukore_2026.11                       # into docs/benchmarks/figures
    python3 docs/benchmarks/plot_benchmarks.py OUT/Eukore_2026.11 --out figs --format svg
    python3 docs/benchmarks/plot_benchmarks.py PKG --benchmarks other/benchmarks          # tables measured elsewhere

Needs matplotlib alone. Reads PACKAGE/benchmarks/ (or --benchmarks) and PACKAGE/busco/loci.tsv; a figure whose table
is missing is left out (the build measures what its inputs allow). Writes:
- completeness.png: the BUSCOs Single, Duplicated, Fragmented and Missing per supergroup, and each genome's share;
- completeness_vs_assembly.png: a genome's completeness against its assembly's contiguity;
- loci_call_rates.png: each BUSCO's call rate per supergroup (busco/loci.tsv), the loci a user chooses among;
- orthodb.png: Eukore's calls against OrthoDB's own gene models of the same species;
- resolution.png: the identity of each copy's nearest other species by the rank they share, BUSCOs and 18S; per
  marker the copies under 99%;
- confidence.png: the confidence tiers of the BUSCO copies and of the rRNA records per source;
- rrna.png: per rRNA source the share of records holding each region;
- labels.png: how the species were named in PR2, EUKARYOME and KSGP;
- runtime.png: each genome's run time against its size, and the package stage's steps;
- summary.md: the headline figures (summary.json) as a table, for the README.
"""

import argparse
import collections
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# the palette (the dataviz reference palette: categorical slots in fixed order, one blue ramp for ordered tiers)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
BLUE = {100: "#cde2fb", 250: "#86b6ef", 400: "#3987e5", 550: "#1c5cab", 700: "#0d366b"}
NONE = "#d6d5d0"  # missing, unnamed: a gray, not a series
TIER_COLORS = {"high": BLUE[550], "medium": BLUE[400], "low": BLUE[250]}
SEQUENTIAL = LinearSegmentedColormap.from_list("blue", ["#f3f7fd", BLUE[100], BLUE[250], BLUE[400], BLUE[550],
                                                        BLUE[700]])
SUPERGROUPS = ("Vertebrata", "Metazoa", "Fungi", "Viridiplantae", "Rhodophyta", "SAR", "Discoba", "Metamonada",
               "Amoebozoa", "Haptista", "Cryptophyceae", "other", "other eukaryotes")
RANKS = ("genus", "family", "order", "class", "phylum", "none")
BUSCOS = "BUSCO (all)"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 11,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.labelsize": 9.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": False, "axes.axisbelow": True, "grid.color": GRID,
    "grid.linewidth": 0.7, "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5, "text.color": INK, "legend.frameon": False, "legend.fontsize": 8.5,
    "font.family": "DejaVu Sans", "savefig.dpi": 160, "savefig.bbox": "tight"})


def read(path):
    """A TSV's rows (dicts), [] where there is none."""
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def num(text, default=None):
    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def group_order(names):
    """Supergroups in the tree's order, unknown ones last."""
    rank = {g: i for i, g in enumerate(SUPERGROUPS)}
    return sorted(set(names), key=lambda g: (rank.get(g, len(SUPERGROUPS)), g))


def save(fig, out, name, fmt):
    path = os.path.join(out, f"{name}.{fmt}")
    fig.savefig(path)
    plt.close(fig)
    return path


def legend_below(ax, handles=None, labels=None, ncol=4):
    kw = {} if handles is None else dict(handles=handles, labels=labels)
    ax.legend(**kw, loc="upper left", bbox_to_anchor=(0, -0.12), ncol=ncol, handlelength=1.2, columnspacing=1.2)


def stacked(ax, rows, parts, colors, labels=None):
    """Horizontal 100% stacked bars: rows [(label, {part: count})], a 2-px surface gap between segments."""
    labels = labels or parts
    for y, (_name, counts) in enumerate(rows):
        total = sum(counts.get(p, 0) for p in parts) or 1
        left = 0.0
        for p, c in zip(parts, colors):
            share = counts.get(p, 0) / total
            ax.barh(y, share, left=left, height=0.62, color=c, edgecolor=SURFACE, linewidth=1.2)
            left += share
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.grid(axis="x")
    handles = [matplotlib.patches.Patch(color=c) for c in colors]
    return handles, list(labels)


# --- the figures -------------------------------------------------------------------------------------------------

def completeness(bench, out, fmt):
    rows = read(os.path.join(bench, "completeness.tsv"))
    if not rows:
        return None
    by = collections.defaultdict(list)
    for r in rows:
        by[r["supergroup"]].append(r)
    groups = group_order(by)
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 0.42 * len(groups) + 1.9), sharey=True,
                               gridspec_kw=dict(width_ratios=(1.25, 1), wspace=0.06))
    parts = ("single", "duplicated", "fragmented", "missing")
    bars = [(f"{g} ({len(by[g])})", {p: sum(int(r[p]) for r in by[g]) for p in parts}) for g in groups]
    handles, labels = stacked(a, bars, parts, SERIES + [NONE], ("Single", "Duplicated", "Fragmented", "Missing"))
    a.set_title("BUSCOs per supergroup (species)")
    legend_below(a, handles, labels)
    for y, g in enumerate(groups):
        values = sorted(float(r["complete"]) for r in by[g])
        jitter = [((i * 0.618) % 1 - 0.5) * 0.36 for i in range(len(values))]
        b.scatter(values, [y + j for j in jitter], s=14, color=SERIES[0], alpha=0.55, linewidths=0)
        med = values[len(values) // 2]
        b.plot([med, med], [y - 0.3, y + 0.3], color=INK, linewidth=2, solid_capstyle="round")
    b.set_xlim(-0.02, 1.02)
    b.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    b.grid(axis="x")
    b.set_title("Each genome's complete BUSCOs (bar: median)")
    b.tick_params(axis="y", length=0)
    return save(fig, out, "completeness", fmt)


def completeness_vs_assembly(bench, out, fmt):
    rows = [r for r in read(os.path.join(bench, "completeness.tsv"))
            if num(r["contigs"]) and num(r["genome_size"])]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    x = [num(r["contigs"]) / (num(r["genome_size"]) / 1e6) for r in rows]
    y = [num(r["complete"]) for r in rows]
    ax.scatter(x, y, s=22, color=SERIES[0], alpha=0.6, linewidths=0.8, edgecolors=SURFACE)
    ax.set_xscale("log")
    ax.set_xlabel("contigs per Mb of assembly")
    ax.set_ylabel("complete BUSCOs")
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_ylim(-0.02, 1.02)
    ax.grid(axis="both")
    ax.set_title(f"Completeness against contiguity ({len(rows)} genomes)")
    return save(fig, out, "completeness_vs_assembly", fmt)


def loci_call_rates(package, out, fmt):
    rows = read(os.path.join(package, "busco", "loci.tsv"))
    if not rows:
        return None
    cols = [c for c in rows[0] if c.startswith("call_") and c != "call_rate"]
    groups = group_order(c[len("call_"):] for c in cols)
    rows.sort(key=lambda r: -num(r["call_rate"], 0))
    matrix = [[num(r["call_rate"], 0) for r in rows]] + [[num(r[f"call_{g}"], 0) for r in rows] for g in groups]
    fig, ax = plt.subplots(figsize=(max(7.0, 0.085 * len(rows) + 2.4), 0.36 * len(matrix) + 1.6))
    im = ax.imshow(matrix, aspect="auto", cmap=SEQUENTIAL, vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks(range(len(matrix)), ["all species"] + groups)
    ax.axhline(0.5, color=SURFACE, linewidth=3)
    ax.set_xticks([])
    ax.set_xlabel(f"the {len(rows)} BUSCOs of eukaryota_odb12, by their call rate over all species")
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    bar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    bar.ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, colors=INK2)
    ax.set_title("Each BUSCO's call rate per supergroup")
    return save(fig, out, "loci_call_rates", fmt)


def orthodb(bench, out, fmt):
    rows = read(os.path.join(bench, "orthodb.tsv"))
    if not rows:
        return None
    verdicts = ("concordant", "similar", "discordant", "missed")
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        by[r["name"]][r["verdict"]] += 1
    names = sorted(by, key=lambda n: (-by[n]["concordant"] / sum(by[n].values()), n))
    if len(names) > 25:  # many species: by supergroup instead
        by = collections.defaultdict(collections.Counter)
        for r in rows:
            by[r["supergroup"]][r["verdict"]] += 1
        names = group_order(by)
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, max(3.4, 0.36 * len(names) + 1.9)),
                               gridspec_kw=dict(width_ratios=(1.3, 1), wspace=0.45))
    handles, labels = stacked(a, [(f"{n} ({sum(by[n].values())})", by[n]) for n in names], verdicts,
                              [BLUE[550], BLUE[400], BLUE[250], NONE],
                              ("concordant (both ≥ 90%)", "similar (≥ 50%)", "discordant", "missed"))
    a.set_title("Calls against OrthoDB's gene models (BUSCOs)")
    for t in a.get_yticklabels():
        t.set_fontstyle("italic" if len(names) <= 25 and " " in t.get_text() else "normal")
    legend_below(a, handles, labels, ncol=2)
    called = [r for r in rows if r["recall"] != "-"]
    b.scatter([num(r["recall"]) for r in called], [num(r["precision"]) for r in called], s=16, color=SERIES[0],
              alpha=0.5, linewidths=0)
    b.add_patch(matplotlib.patches.Rectangle((0.9, 0.9), 0.12, 0.12, fill=False, edgecolor=INK2, linewidth=0.8,
                                             linestyle=(0, (3, 2))))
    b.set_xlim(-0.02, 1.02)
    b.set_ylim(-0.02, 1.02)
    b.set_xlabel("OrthoDB protein's 5-mers in the call (recall)")
    b.set_ylabel("call's 5-mers in OrthoDB's protein (precision)")
    b.grid(axis="both")
    b.set_title(f"Each called BUSCO ({len(called)})")
    return save(fig, out, "orthodb", fmt)


def resolution(bench, out, fmt):
    rows = read(os.path.join(bench, "resolution.tsv"))
    if not rows:
        return None
    at = {(r["marker"], r["shared_rank"]): r for r in rows}
    markers = [m for m in (BUSCOS, "18S") if any(k[0] == m for k in at)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw=dict(width_ratios=(1, 1.25), wspace=0.22))
    width = 0.8 / max(1, len(markers))
    for i, m in enumerate(markers):
        stats, positions = [], []
        for j, rank in enumerate(RANKS):
            r = at.get((m, rank))
            if not r or r["median"] == "-":
                continue
            stats.append(dict(med=num(r["median"]), q1=num(r["q25"]), q3=num(r["q75"]), whislo=num(r["q05"]),
                              whishi=num(r["q95"]), label=rank))
            positions.append(j + (i - (len(markers) - 1) / 2) * width)
        if stats:
            c = SERIES[i]
            a.bxp(stats, positions=positions, widths=width * 0.8, showfliers=False, patch_artist=True,
                  boxprops=dict(facecolor=c, edgecolor=c, alpha=0.85), medianprops=dict(color=SURFACE, linewidth=2),
                  whiskerprops=dict(color=c, linewidth=1.5), capprops=dict(color=c, linewidth=1.5))
            a.scatter(positions, [s["med"] for s in stats], s=10, color=c, zorder=3)  # a one-copy box still shows
    def count(m, rank):  # a box of one copy is a line: its count says so
        r = at.get((m, rank))
        return f"{int(r['copies']):,}" if r and r["median"] != "-" else "0"
    a.set_xticks(range(len(RANKS)), [f"{r}\n" + " | ".join(count(m, r) for m in markers)
                                     for r in RANKS])
    a.set_xlabel("rank shared with the nearest other species (copies: " +
                 " | ".join("BUSCO" if m == BUSCOS else m for m in markers) + ")")
    a.set_ylabel("identity to the nearest other species, %")
    a.axhline(99, color=INK2, linewidth=0.8, linestyle=(0, (3, 2)))
    a.text(len(RANKS) - 0.5, 99.2, "99%", color=INK2, fontsize=8, ha="right", va="bottom")
    a.grid(axis="y")
    a.set_title("Nearest other species, by the rank they share")
    a.legend([matplotlib.patches.Patch(color=SERIES[i]) for i in range(len(markers))],
             ["BUSCO coding sequences" if m == BUSCOS else "18S" for m in markers], loc="lower left")
    loci = sorted(((r["marker"], int(r["under_99"]) / int(r["copies"])) for r in rows
                   if r["shared_rank"] == "all" and r["marker"] not in (BUSCOS, "18S") and int(r["copies"]) >= 3),
                  key=lambda kv: -kv[1])
    if loci:
        b.bar(range(len(loci)), [v for _m, v in loci], width=0.82, color=SERIES[0], edgecolor=SURFACE,
              linewidth=0.4)
        ssu = at.get(("18S", "all"))
        if ssu:
            v = int(ssu["under_99"]) / int(ssu["copies"])
            b.axhline(v, color=SERIES[1], linewidth=2)
            b.text(len(loci) - 0.5, v + 0.015, f"18S: {v:.0%}", color=INK, fontsize=8.5, ha="right", va="bottom",
                   bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5))
        b.set_xticks([])
        b.set_xlabel(f"the {len(loci)} BUSCOs, most resolving first")
        b.set_ylim(0, 1)
        b.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
        b.grid(axis="y")
        b.set_title("Copies under 99% identical to any other species")
    else:
        b.set_axis_off()
    return save(fig, out, "resolution", fmt)


def confidence(bench, out, fmt):
    rows = []
    copies = read(os.path.join(bench, "completeness.tsv"))
    if copies:
        rows.append(("BUSCO copies", {t: sum(int(r[t]) for r in copies) for t in TIER_COLORS}))
    by = collections.defaultdict(collections.Counter)
    for r in read(os.path.join(bench, "rrna.tsv")):
        for t in TIER_COLORS:
            by[r["source"]][t] += int(r[t])
    rows += [(f"rRNA: {s}", c) for s, c in sorted(by.items(), key=lambda kv: -sum(kv[1].values()))]
    if not rows:
        return None
    rows = [(f"{n} ({sum(c.values()):,})", c) for n, c in rows]
    fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(rows) + 1.5))
    handles, labels = stacked(ax, rows, tuple(TIER_COLORS), list(TIER_COLORS.values()))
    legend_below(ax, handles, labels, ncol=3)
    ax.set_title("Confidence tiers (docs/FORMAT.md: the checks that set them)")
    return save(fig, out, "confidence", fmt)


def rrna(bench, out, fmt):
    rows = read(os.path.join(bench, "rrna.tsv"))
    if not rows:
        return None
    parts = ("18S", "ITS1", "5.8S", "ITS2", "28S", "lsu_complete", "5S", "introns")
    names = ("18S", "ITS1", "5.8S", "ITS2", "28S", "whole 28S", "inline 5S", "18S intron")
    labels = [f"{r['source']}, {r['kind']} ({int(r['records']):,})" for r in rows]
    matrix = [[int(r[p]) / max(1, int(r["records"])) for p in parts] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(rows) + 1.6))
    ax.imshow(matrix, aspect="auto", cmap=SEQUENTIAL, vmin=0, vmax=1)
    for i, line in enumerate(matrix):
        for j, v in enumerate(line):
            ax.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=8, color=SURFACE if v > 0.55 else INK)
    ax.set_xticks(range(len(names)), names)
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(labels)), labels)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("rRNA records holding each region", pad=26)
    return save(fig, out, "rrna", fmt)


def labels(bench, out, fmt):
    rows = read(os.path.join(bench, "labels.tsv"))
    if not rows:
        return None
    parts = ("species", "synonym", "genus", "none")
    overall = [(r["reference"], {p: int(r[p]) for p in parts}) for r in rows if r["supergroup"] == "all"]
    groups = group_order(r["supergroup"] for r in rows if r["supergroup"] != "all")
    refs = [r for r, _c in overall]
    share = {(r["supergroup"], r["reference"]): (int(r["species"]) + int(r["synonym"])) /
             max(1, sum(int(r[p]) for p in parts)) for r in rows}
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, max(2.8, 0.36 * len(groups) + 1.6)),
                               gridspec_kw=dict(width_ratios=(1.4, 1), wspace=0.45))
    handles, names = stacked(a, overall, parts, [BLUE[550], BLUE[400], BLUE[250], NONE],
                             ("by its name", "by a synonym", "by its genus", "not named"))
    legend_below(a, handles, names)
    a.set_title("Species named in each reference")
    matrix = [[share.get((g, ref), 0) for ref in refs] for g in groups]
    b.imshow(matrix, aspect="auto", cmap=SEQUENTIAL, vmin=0, vmax=1)
    for i, line in enumerate(matrix):
        for j, v in enumerate(line):
            b.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=8, color=SURFACE if v > 0.55 else INK)
    b.set_xticks(range(len(refs)), refs)
    b.set_yticks(range(len(groups)), groups)
    b.tick_params(length=0)
    for side in ("left", "bottom"):
        b.spines[side].set_visible(False)
    b.set_title("Named to the species, per supergroup")
    return save(fig, out, "labels", fmt)


def runtime(bench, out, fmt):
    rows = [r for r in read(os.path.join(bench, "runtime.tsv")) if num(r["seconds"]) and num(r["genome_size"])]
    stages = [(r["step"], num(r["seconds"], 0)) for r in read(os.path.join(bench, "stages.tsv"))]
    if not rows and not stages:
        return None
    fig, axes = plt.subplots(1, 2 if rows and stages else 1, figsize=(11 if rows and stages else 6.4, 4.2),
                             gridspec_kw=dict(wspace=0.5) if rows and stages else None)
    axes = list(axes) if rows and stages else [axes]
    if rows:
        ax = axes.pop(0)
        x = [num(r["genome_size"]) / 1e6 for r in rows]
        y = [num(r["seconds"]) / 60 for r in rows]
        ax.scatter(x, y, s=24, color=SERIES[0], alpha=0.65, linewidths=0.8, edgecolors=SURFACE)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("genome size, Mb")
        ax.set_ylabel("minutes")
        for axis in (ax.xaxis, ax.yaxis):  # plain numbers, not 6x10^0
            axis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1, 2, 5)))
            axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v:g}"))
            axis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.grid(axis="both", which="major")
        threads = sorted({r["threads"] for r in rows})
        ax.set_title(f"Each genome's run ({len(rows)}; {'/'.join(threads)} threads)")
    if stages:
        ax = axes.pop(0)
        ax.barh(range(len(stages)), [s / 60 for _n, s in stages], height=0.6, color=SERIES[0])
        ax.set_yticks(range(len(stages)), [n for n, _s in stages])
        ax.invert_yaxis()
        ax.set_xlabel("minutes")
        ax.grid(axis="x")
        for y, (_n, s) in enumerate(stages):
            ax.text(s / 60, y, f" {s / 60:.1f}", va="center", fontsize=8, color=INK2)
        ax.set_title("The package stage's steps")
    return save(fig, out, "runtime", fmt)


# --- the summary -------------------------------------------------------------------------------------------------

def pct(v):
    return "-" if v is None else f"{100 * v:.1f}%"


def summary_md(bench, package, out):
    path = os.path.join(bench, "summary.json")
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        s = json.load(fh)
    manifest = {}
    if os.path.isfile(os.path.join(package, "manifest.json")):
        with open(os.path.join(package, "manifest.json")) as fh:
            manifest = json.load(fh)
    lines = ["| benchmark | figure |", "|---|---|"]
    name = os.path.basename(os.path.normpath(package))
    partial = (manifest.get("inputs") or {}).get("partial")
    lines.append(f"| package | {name}{' (a trial: ' + str(len(partial['species_clades'])) + ' clades)' if partial else ''}"
                 f", {s['species']:,} species' genomes, {s['busco_loci']} BUSCOs |")
    c = s.get("completeness") or {}
    lines.append(f"| median complete BUSCOs per genome | {pct(c.get('median'))} "
                 f"({pct(c.get('at_least_90'))} of genomes at 90% or more) |")
    o = s.get("orthodb")
    if o:
        lines.append(f"| BUSCOs OrthoDB has a gene model of, called | {pct(o['called'])} of {o['pairs']:,} "
                     f"({o['species']} species) |")
        lines.append(f"| of those called, OrthoDB's gene model (5-mers ≥ 90% both ways) | "
                     f"{pct(o['concordant_of_called'])} |")
    for marker, r in (s.get("resolution") or {}).items():
        lines.append(f"| {'BUSCO' if marker.startswith('BUSCO') else marker} copies under 99% identical to any other "
                     f"species | {pct(r['under_99'])} of {r['copies']:,} |")
    r = s.get("rrna")
    if r:
        lines.append(f"| rRNA records; operons with a whole 28S | {r['records']:,}; "
                     f"{pct(r['lsu_complete_of_operons'])} of {r['operons']:,} |")
    for ref, v in (s.get("labels") or {}).items():
        lines.append(f"| species named in {ref} by name or synonym | {pct(v)} |")
    t = s.get("runtime")
    if t:
        per_gb = f", {t['seconds_per_gb'] / 60:.0f} min per Gb of genome" if t.get("seconds_per_gb") else ""
        lines.append(f"| genome run (BUSCOs, rRNA, tiara), median | {t['median_seconds'] / 60:.1f} min on "
                     f"{t.get('threads', '?')} threads{per_gb} |")
    for what, why in (s.get("left_out") or {}).items():
        lines.append(f"| {what} | not measured: {why} |")
    path = os.path.join(out, "summary.md")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


FIGURES = (completeness, completeness_vs_assembly, orthodb, resolution, confidence, rrna, labels, runtime)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("package", help="a Eukore package folder (Eukore_<version>/)")
    p.add_argument("--benchmarks", help="the benchmarks' tables (default PACKAGE/benchmarks)")
    p.add_argument("--out", default=os.path.join(HERE, "figures"), help="where the figures go (default "
                                                                        "docs/benchmarks/figures)")
    p.add_argument("--format", default="png", choices=("png", "svg", "pdf"))
    opts = p.parse_args(argv)
    bench = opts.benchmarks or os.path.join(opts.package, "benchmarks")
    if not os.path.isdir(bench):
        p.error(f"{bench}: no benchmarks (python3 -m eukore.benchmarks {opts.package} measures them)")
    os.makedirs(opts.out, exist_ok=True)
    made = {f.__name__: f(bench, opts.out, opts.format) for f in FIGURES}
    made["loci_call_rates"] = loci_call_rates(opts.package, opts.out, opts.format)
    made["summary"] = summary_md(bench, opts.package, opts.out)
    for path in filter(None, made.values()):
        print(path)
    left = [name for name, path in made.items() if not path]
    if left:
        print(f"left out (no table): {', '.join(left)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
