"""A Eukore package's figures (since 2026-11-10; until 2026-11-12 eukore/figures.py, drawn by the package stage into
the package): drawn from a benchmarks folder, compare.py's (the package's own benchmarks and the comparisons with
other references), or a package's own benchmarks/ alone (its figures that need no other reference):

    python3 -m eukore_bench.figures PKG --benchmarks OUT_BENCH                  # into OUT_BENCH/figures
    python3 -m eukore_bench.figures PKG --benchmarks OUT_BENCH --out DIR --format pdf

benchmarks/eukore-figures.sh runs it in the benchmarks' environment (benchmarks/env.yaml: matplotlib). A figure
whose tables are missing is left out, one that fails is logged and left
out: neither stops the others. Writes the figures of FIGURES (Figure 1, then the supplementary ones in their order),
index.md (each one's caption and tables) and summary.md (benchmarks.summary_md: the headline figures as tables).

Figure 1 (at a journal's two-column width, 180 mm), four panels, each a reviewer's question:
  a  How large? Species per supergroup: NCBI's with a genome (the gray track), Eukore's (with a genome; any record),
     and each other reference's (overlap.tsv: compare.py's references, PR2, EUKARYOME, KSGP).
  b  How far from BUSCO's dataset does a call stay right? Left, the dataset species' genomes called again without
     their genus', family's, order's dataset species (holdout.tsv; orthodb.tsv's own-protein calls at `own`): the
     share of OrthoDB's genes called, on the gene of the own-protein call, and OrthoDB's gene model. Right, each
     genome laid on its own annotation (annotation.tsv) at its distance: the share of its introns whose splice sites
     are the annotation's, per annotation source.
  c  Why genes and not the 18S? Per marker (each BUSCO a point; the 18S; the ITS): left, the copies under 99%
     identical to their nearest congener (resolution.tsv); right, the intraspecific panel's species told from their
     congeners, their own variation counted (intraspecific/barcode_gap.tsv).
  d  Does it place a metagenome's reads? Each test genome a sample of its reads with as many bacterial reads: F1 of
     detecting the right species (its species in the reference) and genus (held out) by the sample's depth, per
     reference (detection.tsv); the held-out reads censored by minimap2's list of secondary hits (reads.tsv).
A panel whose tables are missing says so in its place.
"""

import argparse
import collections
import json
import os
import random
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from eukore import benchmarks  # noqa: E402

FOLDER = "figures"  # in the package's benchmarks/

# the palette (the dataviz reference palette: categorical slots in fixed order, one blue ramp for ordered tiers)
SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
CATEGORICAL = SERIES + ["#eda100", "#e87ba4", "#008300"]  # slots 1-6; lines and bars alone past the third
BLUE = {100: "#cde2fb", 250: "#86b6ef", 400: "#3987e5", 550: "#1c5cab", 700: "#0d366b"}
NONE = "#d6d5d0"  # missing, unnamed: a gray, not a series
TIER_COLORS = {"high": BLUE[700], "unverified": BLUE[550], "medium": BLUE[400], "low": BLUE[250]}  # format 6
SEQUENTIAL = LinearSegmentedColormap.from_list("blue", ["#f3f7fd", BLUE[100], BLUE[250], BLUE[400], BLUE[550],
                                                        BLUE[700]])
# the groups in order (data/supergroups.tsv's, format 5), then a format-4 package's SAR and Haptista
SUPERGROUPS = ("Vertebrata", "Metazoa", "Fungi", "other Opisthokonta", "Apusozoa and Breviatea", "Amoebozoa",
               "Viridiplantae", "Rhodophyta", "Glaucophyta", "other Archaeplastida", "Stramenopiles", "Alveolata",
               "Rhizaria", "Telonemia", "other TSAR", "SAR", "Haptophyta", "Centroplasthelida", "other Haptista",
               "Haptista", "Cryptophyceae", "other Cryptista", "Discoba", "Metamonada", "other", "other eukaryotes")
RANKS = ("genus", "family", "order", "class", "phylum", "none")
BUSCOS = benchmarks.ALL_BUSCOS
EUKORE_REFS = ("Eukore", "Eukore exons", "Eukore 18S")  # the read benchmark's references of Eukore's own (reads.py)
SSU_REFERENCES = ("PR2", "EUKARYOME", "KSGP", "SILVA")  # 18S references; the others' markers are genomes' genes
NOT_IN_BREADTH = ("ROD", benchmarks.ODB_NAME)  # ROD's operons are Eukore's records, BUSCO's dataset its queries
TICK_LABELS = {"Eukore exons": "Eukore\nexons", "Eukore 18S": "Eukore\n18S", "EukDetect2": "Euk-\nDetect2",
               "MetaPhlAn4": "MetaPhlAn\n4", "PR2": "PR2\n(18S)"}
FATES = (("chosen", BLUE[550]), ("chosen, flagged", BLUE[400]), ("model species", BLUE[700]),
         ("no assembly (records)", BLUE[250]), ("every assembly excluded", "#b9b8b2"), ("placeholder", NONE),
         ("hybrid", NONE), ("not a species", NONE), ("not in the taxdump", NONE), ("no download", NONE))
MIN_READS_LABEL = "3 or more, 1% or more of those placed,"  # eukore/reads.py's MIN_READS, MIN_SHARE
DETECT_MODES = (("in", "species, its species in the reference", SERIES[0]),
                ("held_out", "genus, its species held out", BLUE[250]),
                ("novel", "class, its order held out", SERIES[2]))
VERDICTS = ("concordant", "similar", "discordant", "missed")
VERDICT_COLORS = [BLUE[550], BLUE[400], BLUE[250], NONE]
SOURCE_NAMES = {"ncbi": "NCBI", "ensembl": "Ensembl", "toga2": "TOGA2", "veupathdb": "VEuPathDB", "-": "NCBI"}

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


def save(fig, out, name, fmt, dpi=None):
    path = os.path.join(out, f"{name}.{fmt}")
    fig.savefig(path, **({"dpi": dpi} if dpi else {}))
    plt.close(fig)
    return path


def legend_below(ax, handles=None, labels=None, ncol=4):
    kw = {} if handles is None else dict(handles=handles, labels=labels)
    ax.legend(**kw, loc="upper left", bbox_to_anchor=(0, -0.12), ncol=ncol, handlelength=1.2, columnspacing=1.2)


def pct_axis(ax, axis="x"):
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))


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
    pct_axis(ax)
    ax.grid(axis="x")
    handles = [matplotlib.patches.Patch(color=c) for c in colors]
    return handles, list(labels)


def tick_label(ref):
    return TICK_LABELS.get(ref, ref)


def species_names(bench):
    """{accession: its species' name} of the package's genomes (completeness.tsv)."""
    return {r["accession"]: r["name"] for r in read(os.path.join(bench, "completeness.tsv"))}


def short_name(name):
    words = name.split()
    return f"{words[0][0]}. {' '.join(words[1:])}" if len(words) > 1 else name


def compared_refs(bench):
    """The other references of reads.tsv, in compare.py's order (BENCH/compare/references.tsv)."""
    return list(dict.fromkeys(r["reference"] for r in read(os.path.join(bench, "reads.tsv"))
                              if r["reference"] not in EUKORE_REFS))


def load_summary(bench):
    try:
        with open(os.path.join(bench, "summary.json")) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def detection(bench):
    """summary.json's detection: {reference: {"<mode>_<depth>": {tp, fn_absent, fn_missed, fp, f1, ...}}}."""
    return load_summary(bench).get("detection") or {}


def detection_refs(bench, det):
    refs = [r for r in EUKORE_REFS if r in det] + [r for r in compared_refs(bench) if r in det]
    return refs + [r for r in det if r not in refs]


def split_key(key):
    """detection's "<mode>_<depth>" -> (mode, depth)."""
    mode, depth = key.rsplit("_", 1)
    return mode, int(depth)


def ref_style(ref, others):
    """A read benchmark reference's line: Eukore's own in its blue (whole, exons dashed, 18S dotted), the others in
    the next categorical slots, each its own marker too."""
    if ref in EUKORE_REFS:
        return dict(color=SERIES[0], linestyle=("-", (0, (4, 1.6)), (0, (1, 1.4)))[EUKORE_REFS.index(ref)],
                    marker="o", markerfacecolor=(SERIES[0], SURFACE, SURFACE)[EUKORE_REFS.index(ref)])
    i = others.index(ref) if ref in others else len(others)
    return dict(color=CATEGORICAL[1 + i % (len(CATEGORICAL) - 1)], linestyle="-",
                marker=("s", "D", "^", "v", "P")[i % 5], markerfacecolor=CATEGORICAL[1 + i % (len(CATEGORICAL) - 1)])


# --- Figure 1 --------------------------------------------------------------------------------------------------------

FIG1_SIZE = (7.2, 7.4)  # inches: a journal's two-column width (~180 mm)
FIG1_RC = {"font.size": 7, "axes.titlesize": 7, "axes.titleweight": "normal", "axes.labelsize": 7,
           "xtick.labelsize": 6.3, "ytick.labelsize": 6.3, "legend.fontsize": 6.3, "axes.linewidth": 0.6,
           "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
           "lines.linewidth": 1.3, "lines.markersize": 3.5, "grid.linewidth": 0.5}
NEAREST = ("own", "genus", "family", "order", "class", "phylum", "none")  # the rank shared with BUSCO's dataset
NEAREST_LABELS = ("own sp.", "genus", "family", "order", "class", "phylum", "none")


def panel_title(ax, letter, text):
    ax.annotate(f"{letter}   {text}", xy=(0, 1), xycoords="axes fraction", xytext=(0, 26), textcoords="offset points",
                fontsize=8, fontweight="bold", ha="left", va="bottom", color=INK)


def not_measured(axes, why):
    """A panel's place, its tables missing: the reason, in its axes."""
    axes = list(axes)
    for ax in axes:
        ax.set_axis_off()
    axes[0].text(0.02, 0.5, f"not measured: {why}", transform=axes[0].transAxes, fontsize=6.5, color=MUTED,
                 va="center")


def nearest_axis(ax):
    ax.set_xticks(range(len(NEAREST)), NEAREST_LABELS, rotation=55, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.5, len(NEAREST) - 0.5)
    ax.set_ylim(-0.02, 1.04)
    pct_axis(ax, "y")
    ax.grid(axis="y")


def breadth_refs(bench, overlap):
    """Panel a's references: the read benchmark's (compare.py's, in its order), then the 18S references; five at
    most."""
    names = list(dict.fromkeys(r["database"] for r in overlap))
    compared = [n for n in compared_refs(bench) if n in names]
    return [n for n in compared + [n for n in names if n not in compared] if n not in NOT_IN_BREADTH][:5]


def breadth_marks(refs):
    """{reference: (color, marker)}: a genome-marker reference orange, an 18S reference green (three slots at most:
    dots meet in every pairing), each its own marker within its kind."""
    out, kinds = {}, collections.Counter()
    for ref in refs:
        kind = "ssu" if ref in SSU_REFERENCES else "genes"
        out[ref] = (SERIES[2] if kind == "ssu" else SERIES[1], ("s", "D", "^")[kinds[kind] % 3])
        kinds[kind] += 1
    return out


def panel_breadth(ax, bench):
    breadth = {r["supergroup"]: r for r in read(os.path.join(bench, "breadth.tsv"))}
    overlap = read(os.path.join(bench, "overlap.tsv"))
    panel_title(ax, "a", "Species per supergroup")
    if not breadth:
        not_measured([ax], "breadth.tsv")
        return False
    count = {(r["database"], r["supergroup"]): int(r["species"]) for r in overlap}
    every = {r["supergroup"]: int(r["both"]) + int(r["eukore_only"]) for r in overlap}  # Eukore's, any record
    refs = breadth_refs(bench, overlap)
    marks = breadth_marks(refs)
    groups = group_order(g for g in breadth if g != "all")
    rows = (["all"] if "all" in breadth else []) + groups
    ys = [0.0] + [i + 1.5 for i in range(len(groups))] if rows and rows[0] == "all" else list(range(len(rows)))
    top = 1
    for g, y in zip(rows, ys):
        r = breadth[g]
        ncbi = int(num(r["ncbi_species"], 0))
        genome = int(num(r["eukore"], 0))
        if ncbi:
            ax.plot([0, ncbi], [y, y], color=NONE, linewidth=4, solid_capstyle="round", zorder=1)
        for i, ref in enumerate(refs):
            color, marker = marks[ref]
            ax.scatter(count.get((ref, g), 0), y + (i - (len(refs) - 1) / 2) * 0.09, s=12, marker=marker, color=color,
                       edgecolors=SURFACE, linewidths=0.5, zorder=3)
        if every.get(g, genome) > genome:
            ax.scatter(every[g], y, s=26, facecolors=SURFACE, edgecolors=SERIES[0], linewidths=1.1, zorder=4)
        ax.scatter(genome, y, s=26, color=SERIES[0], edgecolors=SURFACE, linewidths=0.6, zorder=5)
        ax.text(1.01, y, f"{genome:,}", transform=ax.get_yaxis_transform(), fontsize=6, color=INK2, va="center")
        top = max(top, ncbi, genome, every.get(g, 0), *(count.get((ref, g), 0) for ref in refs))
    ax.text(1.01, ys[0] - 1.0, "Eukore", transform=ax.get_yaxis_transform(), fontsize=6, color=INK2, va="center")
    ax.set_yticks(ys, ["all eukaryotes" if g == "all" else g for g in rows])
    ax.set_ylim(ys[-1] + 0.7, ys[0] - 0.7)
    ax.tick_params(axis="y", length=0)
    ax.set_xscale("symlog", linthresh=1, linscale=0.35)
    ax.set_xlim(-0.15, top * 1.6)
    ax.xaxis.set_major_locator(matplotlib.ticker.SymmetricalLogLocator(base=10, linthresh=1))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v:,.0f}"))
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.grid(axis="x")
    ax.set_xlabel("species (log scale)")
    handles = [Line2D([], [], color=NONE, linewidth=4, solid_capstyle="round"),
               Line2D([], [], linestyle="", marker="o", color=SERIES[0], markersize=4.5),
               Line2D([], [], linestyle="", marker="o", markerfacecolor=SURFACE, markeredgecolor=SERIES[0],
                      markersize=4.5)]
    labels = ["NCBI: species with a genome", "Eukore: with a genome", "Eukore: any record (rRNA alone too)"]
    for ref in refs:
        color, marker = marks[ref]
        handles.append(Line2D([], [], linestyle="", marker=marker, color=color, markersize=4))
        labels.append(f"{ref}{' (18S)' if ref in SSU_REFERENCES else ''}")
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(-0.02, -0.17), ncol=2, handlelength=1.3,
              columnspacing=1.0, handletextpad=0.4)
    return True


def reach_points(bench):
    """Panel b's left: {rank shared with the nearest dataset species kept: Counter(loci, called, same_gene, concordant,
    genomes)}, holdout.tsv's levels pooled (a genome's levels that keep the same nearest rank, once: the first, the
    fewest left out) and orthodb.tsv's own-protein calls at `own`; and per genome {accession: {rank: Counter}}."""
    pooled, per_genome, seen = collections.defaultdict(collections.Counter), collections.defaultdict(dict), set()
    for r in read(os.path.join(bench, "holdout.tsv")):
        rank, key = r["nearest"], (r["accession"], r["nearest"])
        if rank not in NEAREST or key in seen:
            continue
        seen.add(key)
        c = collections.Counter({k: int(r[k]) for k in ("loci", "called", "same_gene", "concordant")}, genomes=1)
        pooled[rank].update(c)
        per_genome[r["accession"]][rank] = c
    own = collections.defaultdict(collections.Counter)
    for r in read(os.path.join(bench, "orthodb.tsv")):
        c = own[r["accession"]]
        c["loci"] += 1
        c["called"] += r["verdict"] != "missed"
        c["concordant"] += r["verdict"] == "concordant"
    for acc, c in own.items():
        pooled["own"].update(c)
        pooled["own"]["genomes"] += 1
        per_genome[acc]["own"] = c
    return pooled, per_genome


REACH = (("called", "called", BLUE[250]), ("same_gene", "on the gene its own proteins find", BLUE[400]),
         ("concordant", "OrthoDB's gene model", BLUE[700]))


def annotation_points(bench, min_introns=5):
    """Panel b's right: per genome with min_introns called introns laid on its annotation, (rank shared with the
    nearest dataset species, or own; annotation source; share of its introns with both splice sites the
    annotation's; introns)."""
    by, meta = collections.defaultdict(lambda: [0, 0]), {}
    for r in read(os.path.join(bench, "annotation.tsv")):
        if r["transcript"] == "-":
            continue
        by[r["accession"]][0] += int(r["introns"])
        by[r["accession"]][1] += int(r["exact_introns"])
        meta[r["accession"]] = ("own" if r["odb_self"] == "t" else r["odb_nearest"],
                                SOURCE_NAMES.get(r.get("source") or "-", r.get("source")))
    return [(meta[a][0], meta[a][1], exact / introns, introns) for a, (introns, exact) in by.items()
            if introns >= min_introns and meta[a][0] in NEAREST]


def boxes(ax, groups, color, width, marker="o"):
    """Per x: its values' points (jittered) and, with 5 or more, a box (quartiles; whiskers 5% and 95%) and median."""
    rng = random.Random(1)
    for x, values in groups.items():
        values = sorted(values)
        ax.scatter([x + rng.uniform(-width / 2.6, width / 2.6) for _v in values], values, s=4, color=color,
                   alpha=0.55, linewidths=0, marker=marker, zorder=2)
        if len(values) >= 5:
            q = [benchmarks.quantile(values, p) for p in (0.05, 0.25, 0.5, 0.75, 0.95)]
            ax.bxp([dict(whislo=q[0], q1=q[1], med=q[2], q3=q[3], whishi=q[4])], positions=[x], widths=width,
                   showfliers=False, patch_artist=True, manage_ticks=False, zorder=3,
                   boxprops=dict(facecolor="none", edgecolor=color, linewidth=0.9),
                   medianprops=dict(color=INK, linewidth=1.1), whiskerprops=dict(color=color, linewidth=0.8),
                   capprops=dict(color=color, linewidth=0.8))
        else:
            ax.plot([x - width / 2, x + width / 2], [values[len(values) // 2]] * 2, color=INK, linewidth=1.1,
                    zorder=3)


def panel_reach(left, right, bench):
    panel_title(left, "b", "Calls by distance to BUSCO's dataset")
    pooled, _per = reach_points(bench)
    points = annotation_points(bench)
    if not pooled and not points:
        not_measured([left, right], "holdout.tsv, orthodb.tsv, annotation.tsv")
        return False
    # the hold-out's levels alone: `own` (the own-protein calls) pools every dataset genome, the far ones too, and so
    # sits below the nearer levels it is no baseline of (the holdout figure draws each genome's line from it)
    for key, label, color in REACH:
        xs = [i for i, r in enumerate(NEAREST) if r != "own" and pooled.get(r, {}).get("loci")]
        ys = [pooled[NEAREST[i]][key] / pooled[NEAREST[i]]["loci"] for i in xs]
        if xs:
            left.plot(xs, ys, color=color, marker="o", label=label, zorder=3)
    for i, r in enumerate(NEAREST):
        if r != "own" and pooled.get(r, {}).get("genomes"):
            left.text(i, 0.01, f"{pooled[r]['genomes']:,}", ha="center", va="bottom", fontsize=5.5, color=MUTED)
    nearest_axis(left)
    left.set_title("the dataset's genomes,\nrelatives left out", fontsize=6.5, color=INK2)
    left.set_ylabel("share of OrthoDB's genes (left);\nof a genome's introns (right)", fontsize=6.3)
    sources = [s for s, _n in collections.Counter(s for _d, s, _v, _i in points).most_common()]
    width = 0.7 / max(1, len(sources))
    for k, source in enumerate(sources):
        color = CATEGORICAL[1 + k % (len(CATEGORICAL) - 1)]
        groups = collections.defaultdict(list)
        for dist, s, share, _n in points:
            if s == source:
                groups[NEAREST.index(dist) + (k - (len(sources) - 1) / 2) * width].append(share)
        boxes(right, groups, color, width * 0.85, ("o", "s", "D", "^")[k % 4])
    for i, r in enumerate(NEAREST):
        n = sum(1 for d, *_x in points if d == r)
        if n:
            right.text(i, 0.01, f"{n:,}", ha="center", va="bottom", fontsize=5.5, color=MUTED)
    nearest_axis(right)
    right.tick_params(axis="y", labelleft=False)
    right.set_title("each genome against its\nannotation: introns exact", fontsize=6.5, color=INK2)
    if not points:
        right.text(0.5, 0.5, "not measured:\nannotation.tsv", transform=right.transAxes, ha="center", fontsize=6,
                   color=MUTED)
    handles, labels = left.get_legend_handles_labels()
    for k, source in enumerate(sources):
        handles.append(Line2D([], [], linestyle="", marker=("o", "s", "D", "^")[k % 4],
                              color=CATEGORICAL[1 + k % (len(CATEGORICAL) - 1)], markersize=3.5))
        labels.append(f"{source}'s annotation (box: median, quartiles)")
    left.legend(handles, labels, loc="upper left", bbox_to_anchor=(-0.02, -0.33), ncol=1, handlelength=1.6)
    left.annotate("rank shared with the nearest dataset species (genomes at the foot)", xy=(1.06, -0.26),
                  xycoords="axes fraction", ha="center", va="top", fontsize=6.3, color=INK2)
    return True


def resolution_points(bench):
    """Panel c's left: {marker: (copies, share under 99% identical to the nearest congener)} (resolution.tsv's
    shared rank genus), each BUSCO apart, BUSCOS pooled."""
    out = {}
    for r in read(os.path.join(bench, "resolution.tsv")):
        if r["shared_rank"] == "genus" and int(r["copies"]):
            out[r["marker"]] = (int(r["copies"]), int(r["under_99"]) / int(r["copies"]))
    return out


def gap_points(package):
    """Panel c's right: {marker: (species measured, share resolved from their congeners)} (barcode_gap.tsv), each
    BUSCO apart, BUSCOS pooled."""
    counts = collections.defaultdict(collections.Counter)
    for r in read(os.path.join(package, "intraspecific", "barcode_gap.tsv")):
        if r["resolved"] in ("t", "f"):
            for m in (r["marker"], BUSCOS) if r["marker"].endswith("at2759") else (r["marker"],):
                counts[m][r["resolved"]] += 1
    return {m: (c["t"] + c["f"], c["t"] / (c["t"] + c["f"])) for m, c in counts.items()}


def own_loci(package):
    """Eukore's own marker loci (markers/loci.tsv: RPB1 .. H3; format 5), an empty set without."""
    return {r["locus"] for r in read(os.path.join(package, "markers", "loci.tsv"))}


MARKER_ROWS = ("BUSCOs", "Eukore's loci", "18S", "ITS")


def marker_rows(values, own):
    """{row: (each locus' share, of those with 3 or more; (n, share) of the row pooled, or None)} of MARKER_ROWS from
    {marker: (n, share)}: the BUSCOs (each a locus), Eukore's own loci (each a locus), the 18S, the ITS."""
    def pooled(names):
        got = [values[m] for m in names]
        n = sum(x for x, _s in got)
        return (n, sum(x * s for x, s in got) / n) if n else None
    busco = [m for m in values if m.endswith("at2759")]
    mine = [m for m in values if m in own]
    return {"BUSCOs": ([values[m][1] for m in busco if values[m][0] >= 3], values.get(BUSCOS) or pooled(busco)),
            "Eukore's loci": ([values[m][1] for m in mine if values[m][0] >= 3], pooled(mine)),
            "18S": ([], values.get("18S")), "ITS": ([], values.get("ITS"))}


def marker_strip(ax, rows, names):
    """One row per marker kind: each locus' share a point, the row's pooled share a bar with its value and n."""
    rng = random.Random(1)
    for y, name in enumerate(names):
        each, pooled = rows[name]
        if each:
            ax.scatter(each, [y + rng.uniform(-0.2, 0.2) for _v in each], s=5, color=BLUE[250], linewidths=0,
                       zorder=2)
        if pooled:
            n, share = pooled
            ax.plot([share, share], [y - 0.3, y + 0.3], color=INK if each else SERIES[0], linewidth=1.8,
                    solid_capstyle="butt", zorder=3)
            if not each:
                ax.scatter(share, y, s=18, color=SERIES[0], edgecolors=SURFACE, linewidths=0.6, zorder=4)
            ax.text(share, y - 0.34, f"{share:.0%} of {n:,}", ha="center" if 0.15 < share < 0.85 else
                    "left" if share <= 0.15 else "right", va="bottom", fontsize=5.6, color=INK)
        else:
            ax.text(0.5, y, "not measured", ha="center", va="center", fontsize=5.8, color=MUTED)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(len(names) - 0.45, -0.8)
    pct_axis(ax)
    ax.grid(axis="x")
    ax.tick_params(axis="y", length=0)


def panel_resolution(left, right, bench, package):
    panel_title(left, "c", "Species told apart, by marker")
    own = own_loci(package)
    res, gap = marker_rows(resolution_points(bench), own), marker_rows(gap_points(package), own)
    names = [n for n in MARKER_ROWS if res[n][1] or gap[n][1]]
    if not names:
        not_measured([left, right], "resolution.tsv, intraspecific/barcode_gap.tsv")
        return False
    marker_strip(left, res, names)
    marker_strip(right, gap, names)
    left.set_yticks(range(len(names)), [f"{n}\n(each a dot)" if res[n][0] or gap[n][0] else n for n in names])
    right.tick_params(axis="y", labelleft=False)
    left.set_title("copies <99% identical\nto the nearest congener", fontsize=6.5, color=INK2)
    right.set_title("species told from congeners\n(barcode gap: own variation)", fontsize=6.5, color=INK2)
    left.annotate("share (bar: the row's pooled, of n copies or species)", xy=(1.06, -0.13),
                  xycoords="axes fraction", ha="center", va="top", fontsize=6, color=INK2)
    return True


def censored_share(bench, mode):
    """Of the Eukore's genome reads in a mode, the share censored (reads.py: minimap2's secondary list all left out),
    or None (no such column: before 2026-11-10)."""
    n = c = 0
    for r in read(os.path.join(bench, "reads.tsv")):
        if (r["reference"], r["mode"], r["kind"], r["reads"]) == ("Eukore", mode, "eukaryote", "genome"):
            if "censored" not in r:
                return None
            n += int(r["n"])
            c += int(r["censored"] or 0)
    return c / n if n else None


def panel_detection(left, right, bench):
    panel_title(left, "d", "Detecting a sample's eukaryote from its reads")
    det = detection(bench)
    refs = detection_refs(bench, det)
    if not refs:
        not_measured([left, right], "detection.tsv (eukore_bench.compare)")
        return False
    others = [r for r in refs if r not in EUKORE_REFS]
    depths = sorted({split_key(k)[1] for r in refs for k in det[r]})
    for ax, (mode, label, _c) in zip((left, right), DETECT_MODES):
        samples = 0
        for ref in refs:
            pts = sorted((split_key(k)[1], v.get("f1")) for k, v in det[ref].items()
                         if split_key(k)[0] == mode and v.get("f1") is not None)
            samples = max([samples] + [v.get("samples", 0) for k, v in det[ref].items() if split_key(k)[0] == mode])
            if pts:
                ax.plot([d for d, _v in pts], [v for _d, v in pts], label=ref, markeredgewidth=0.8,
                        **ref_style(ref, others))
        ax.set_xscale("log")
        ax.set_xticks(depths, [f"{d / 1000:g}k" for d in depths])
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xlim(depths[0] / 1.5, depths[-1] * 1.5)
        ax.set_ylim(-0.03, 1.04)
        ax.grid(axis="y")
        ax.set_title(f"{label.split(', ', 1)[1].replace('its ', '')}:\nright {label.split(',')[0]} ({samples} "
                     f"samples)", fontsize=6.5, color=INK2)
        cens = censored_share(bench, mode) if mode != "in" else None
        if cens:
            ax.text(0.98, 0.04, f"reads censored: {cens:.1%}", transform=ax.transAxes, ha="right", fontsize=5.5,
                    color=MUTED)
    left.set_ylabel("F1 of detection", fontsize=6.5)
    right.tick_params(axis="y", labelleft=False)
    left.annotate("reads of the genome per sample (+ as many bacterial reads)", xy=(1.06, -0.13),
                  xycoords="axes fraction", ha="center", va="top", fontsize=6, color=INK2)
    handles, labels = left.get_legend_handles_labels()
    left.legend(handles, labels, loc="upper left", bbox_to_anchor=(-0.02, -0.25), ncol=3, handlelength=2.4,
                columnspacing=1.0)
    return True


def figure1(bench, package, out, fmt):
    with plt.rc_context(FIG1_RC):
        fig = plt.figure(figsize=FIG1_SIZE)
        outer = fig.add_gridspec(2, 2, width_ratios=(1, 1.08), height_ratios=(1.12, 1), hspace=0.95, wspace=0.62,
                                 left=0.14, right=0.97, top=0.94, bottom=0.17)
        b = outer[0, 1].subgridspec(1, 2, wspace=0.12)
        c = outer[1, 0].subgridspec(1, 2, wspace=0.24)
        d = outer[1, 1].subgridspec(1, 2, wspace=0.12)
        b1 = fig.add_subplot(b[0])
        c1 = fig.add_subplot(c[0])
        d1 = fig.add_subplot(d[0])
        drawn = [panel_breadth(fig.add_subplot(outer[0, 0]), bench),
                 panel_reach(b1, fig.add_subplot(b[1], sharey=b1), bench),
                 panel_resolution(c1, fig.add_subplot(c[1], sharey=c1), bench, package),
                 panel_detection(d1, fig.add_subplot(d[1], sharey=d1), bench)]
        if not any(drawn):
            plt.close(fig)
            return None
        return save(fig, out, "Fig1", fmt, dpi=300)


# --- the supplementary figures ---------------------------------------------------------------------------------------

def selection(bench, package, out, fmt):
    rows = read(os.path.join(bench, "selection.tsv"))
    if not rows:
        return None
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        by[r["supergroup"]][r["fate"]] += int(r["taxa"])
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(by) + 1.8))
    handles, labels = stacked(ax, [(f"{g} ({sum(by[g].values()):,})", by[g]) for g in group_order(by)],
                              [f for f, _c in FATES], [c for _f, c in FATES])
    legend_below(ax, handles, labels, ncol=3)
    ax.set_title("What became of NCBI's taxa with a GenBank assembly, per supergroup")
    return save(fig, out, "selection", fmt)


def overlap(bench, package, out, fmt):
    rows = read(os.path.join(bench, "overlap.tsv"))
    if not rows:
        return None
    # the databases by their species, most first; the supergroups as the other figures order them, `all` last
    names = sorted({r["database"] for r in rows},
                   key=lambda n: -next(int(r["species"]) for r in rows if r["database"] == n and r["supergroup"] == "all"))
    groups = group_order({r["supergroup"] for r in rows if r["supergroup"] != "all"}) + ["all"]
    at = {(r["database"], r["supergroup"]): r for r in rows}
    cols = min(3, len(names))
    nrow = -(-len(names) // cols)
    fig, axes = plt.subplots(nrow, cols, figsize=(4.2 * cols, (0.32 * len(groups) + 1.0) * nrow), sharey=True,
                             squeeze=False)
    for ax, name in zip(axes.flat, names):
        left = [0] * len(groups)
        for part, color, label in (("eukore_only", BLUE[250], "Eukore only"), ("both", BLUE[700], "in both"),
                                   ("database_only", SERIES[1], "the database only")):
            vals = [int(at[(name, g)][part]) if (name, g) in at else 0 for g in groups]
            ax.barh(range(len(groups)), vals, left=left, color=color, height=0.65, edgecolor=SURFACE, label=label)
            left = [a + b for a, b in zip(left, vals)]
        ax.set_yticks(range(len(groups)), groups)
        ax.grid(axis="x")
        ax.set_title(name, fontsize=10)
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4, integer=True))
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    axes[0][0].invert_yaxis()
    for ax in axes[-1]:
        ax.set_xlabel("species")
    fig.legend(*axes[0][0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3,
               frameon=False, fontsize=8.5)
    fig.suptitle("Species in Eukore (any record) and in each database, per supergroup", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    return save(fig, out, "overlap", fmt)


def completeness(bench, package, out, fmt):
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
    pct_axis(b)
    b.grid(axis="x")
    b.set_title("Each genome's complete BUSCOs (bar: median)")
    b.tick_params(axis="y", length=0)
    return save(fig, out, "completeness", fmt)


def completeness_vs_assembly(bench, package, out, fmt):
    rows = [r for r in read(os.path.join(bench, "completeness.tsv")) if num(r["contigs"]) and num(r["genome_size"])]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    x = [num(r["contigs"]) / (num(r["genome_size"]) / 1e6) for r in rows]
    y = [num(r["complete"]) for r in rows]
    ax.scatter(x, y, s=22 if len(rows) < 500 else 6, color=SERIES[0], alpha=0.6 if len(rows) < 500 else 0.3,
               linewidths=0.8 if len(rows) < 500 else 0, edgecolors=SURFACE)
    ax.set_xscale("log")
    ax.set_xlabel("contigs per Mb of assembly")
    ax.set_ylabel("complete BUSCOs")
    pct_axis(ax, "y")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(axis="both")
    ax.set_title(f"Completeness against contiguity ({len(rows):,} genomes)")
    return save(fig, out, "completeness_vs_assembly", fmt)


def loci_call_rates(bench, package, out, fmt):
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


def orthodb(bench, package, out, fmt):
    rows = read(os.path.join(bench, "orthodb.tsv"))
    if not rows:
        return None
    held = [r for r in rows if r.get("held_out_verdict", "-") != "-"]  # called without the species' own proteins
    group = "name" if len({r["name"] for r in rows}) <= 25 else "supergroup"  # many species: by supergroup

    def count(rs, column):
        by = collections.defaultdict(collections.Counter)
        for r in rs:
            by[r[group]][r[column]] += 1
        return by
    by, by_held = count(rows, "verdict"), count(held, "held_out_verdict")
    names = sorted(by, key=lambda n: (-by[n]["concordant"] / sum(by[n].values()), n)) if group == "name" \
        else group_order(by)
    items = []
    for n in names:
        items.append((f"{n} ({sum(by[n].values())})", by[n]))
        if n in by_held:
            items.append(("own proteins left out", by_held[n]))
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, max(3.4, 0.36 * len(items) + 1.9)),
                               gridspec_kw=dict(width_ratios=(1.3, 1), wspace=0.45))
    handles, labels = stacked(a, items, VERDICTS, VERDICT_COLORS,
                              ("concordant (both ≥ 90%)", "similar (≥ 50%)", "discordant", "missed"))
    a.set_title("Calls against OrthoDB's gene models (BUSCOs)")
    for t in a.get_yticklabels():
        text = t.get_text()
        t.set_fontstyle("italic" if group == "name" and " " in text and not text.startswith("own ") else "normal")
    legend_below(a, handles, labels, ncol=2)
    called = [r for r in rows if r["recall"] != "-"]
    held_called = [r for r in held if r["held_out_recall"] != "-"]
    if held_called:  # the held-out calls are the measure; the build's own, queried with OrthoDB's protein, behind
        b.scatter([num(r["recall"]) for r in called], [num(r["precision"]) for r in called], s=12, color=NONE,
                  linewidths=0, label="queried with its own proteins")
        b.scatter([num(r["held_out_recall"]) for r in held_called], [num(r["held_out_precision"]) for r in held_called],
                  s=16, color=SERIES[0], alpha=0.5, linewidths=0, label="its own proteins left out")
        b.legend(loc="lower left", fontsize=8, frameon=False)
        called = held_called
    else:
        b.scatter([num(r["recall"]) for r in called], [num(r["precision"]) for r in called], s=16, color=SERIES[0],
                  alpha=0.5, linewidths=0)
    b.add_patch(matplotlib.patches.Rectangle((0.9, 0.9), 0.12, 0.12, fill=False, edgecolor=INK2, linewidth=0.8,
                                             linestyle=(0, (3, 2))))
    b.set_xlim(-0.02, 1.02)
    b.set_ylim(-0.02, 1.02)
    b.set_xlabel("OrthoDB protein's 5-mers in the call (recall)")
    b.set_ylabel("call's 5-mers in OrthoDB's protein (precision)")
    b.grid(axis="both")
    b.set_title(f"Each called BUSCO ({len(called):,})")
    return save(fig, out, "orthodb", fmt)


def holdout(bench, package, out, fmt):
    """The clade hold-out: per measure, each genome's line over its levels (thin) and all genomes pooled (bold)."""
    pooled, per_genome = reach_points(bench)
    if not any(r != "own" for r in pooled):
        return None
    names = species_names(bench)
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.0), sharey=True, gridspec_kw=dict(wspace=0.08))
    for ax, (key, label, color) in zip(axes, REACH):
        for acc, ranks in per_genome.items():
            xs = [i for i, r in enumerate(NEAREST) if r in ranks and ranks[r]["loci"] and
                  (key != "same_gene" or r != "own")]
            if len(xs) > 1:
                ax.plot(xs, [ranks[NEAREST[i]][key] / ranks[NEAREST[i]]["loci"] for i in xs], color=MUTED,
                        linewidth=0.7, alpha=0.5 if len(per_genome) < 40 else 0.15, zorder=1)
                if len(per_genome) <= 12:
                    ax.annotate(short_name(names.get(acc, acc)), (xs[-1], ranks[NEAREST[xs[-1]]][key] /
                                                                  ranks[NEAREST[xs[-1]]]["loci"]),
                                xytext=(3, 0), textcoords="offset points", fontsize=6.5, color=INK2,
                                fontstyle="italic", va="center")
        xs = [i for i, r in enumerate(NEAREST) if pooled.get(r, {}).get("loci") and (key != "same_gene" or r != "own")]
        ax.plot(xs, [pooled[NEAREST[i]][key] / pooled[NEAREST[i]]["loci"] for i in xs], color=color, linewidth=2.4,
                marker="o", markersize=6, zorder=3)
        for i in xs:
            ax.text(i, 0.01, f"{pooled[NEAREST[i]]['genomes']:,}", ha="center", va="bottom", fontsize=7, color=MUTED)
        ax.set_xticks(range(len(NEAREST)), NEAREST_LABELS, rotation=40, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.5, len(NEAREST) - 0.5)
        ax.set_ylim(-0.02, 1.04)
        pct_axis(ax, "y")
        ax.grid(axis="y")
        ax.set_title(label[0].upper() + label[1:], fontsize=10)
    axes[0].set_ylabel("share of OrthoDB's genes of the genome")
    axes[1].set_xlabel("rank shared with the nearest BUSCO dataset species kept (genomes at the foot; bold: all "
                       "pooled)")
    fig.suptitle("The clade hold-out: each dataset species' genome called again without its genus', family's and "
                 "order's dataset species", x=0.01, ha="left", fontweight="bold", fontsize=10.5)
    return save(fig, out, "holdout", fmt)


def annotation(bench, package, out, fmt):
    rows = [r for r in read(os.path.join(bench, "annotation.tsv")) if r["transcript"] != "-"]
    if not rows:
        return None
    by, source = collections.defaultdict(collections.Counter), {}
    for r in rows:
        c = by[r["accession"]]
        for k in ("bases", "ann_bases", "shared_bases", "introns", "exact_introns", "exons", "exact_exons"):
            c[k] += int(r[k])
        c["copies"] += 1
        source[r["accession"]] = SOURCE_NAMES.get(r.get("source") or "-", r.get("source"))
    sources = [s for s, _n in collections.Counter(source.values()).most_common()]
    many = len(by) > 40
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.3 if many else max(4.3, 0.22 * len(by) + 1.5)),
                               gridspec_kw=dict(wspace=0.3))
    for k, s in enumerate(sources):
        mine = [c for acc, c in by.items() if source[acc] == s]
        a.scatter([c["shared_bases"] / c["ann_bases"] for c in mine], [c["shared_bases"] / c["bases"] for c in mine],
                  s=[(6 if many else 20) + c["copies"] / (8 if many else 2) for c in mine],
                  color=CATEGORICAL[1 + k % 5], alpha=0.6, edgecolors=SURFACE, linewidths=0.4, label=s)
    rec = [c["shared_bases"] / c["ann_bases"] for c in by.values()]
    prec = [c["shared_bases"] / c["bases"] for c in by.values()]
    lo = min(rec + prec + [0.8]) - 0.02
    a.set_xlim(lo, 1.005)
    a.set_ylim(lo, 1.005)
    pct_axis(a)
    pct_axis(a, "y")
    a.set_xlabel("the annotation's coding bases in the call (recall)")
    a.set_ylabel("the call's coding bases in the annotation's model (precision)")
    a.grid(axis="both")
    a.legend(loc="lower left", fontsize=8, title="annotation", title_fontsize=8)
    a.set_title(f"Coding bases, each genome ({len(by):,})")
    if many:  # per source and distance to BUSCO's dataset
        points = annotation_points(bench)
        width = 0.7 / max(1, len(sources))
        for k, s in enumerate(sources):
            groups = collections.defaultdict(list)
            for dist, src, share, _n in points:
                if src == s:
                    groups[NEAREST.index(dist) + (k - (len(sources) - 1) / 2) * width].append(share)
            boxes(b, groups, CATEGORICAL[1 + k % 5], width * 0.85)
        nearest_axis(b)
        b.set_xlabel("rank shared with the nearest BUSCO dataset species")
        b.set_title("Introns with both splice sites the annotation's, each genome")
    else:
        genomes = sorted(by, key=lambda acc: -by[acc]["introns"])
        b.barh(range(len(genomes)), [by[g]["exact_introns"] / by[g]["introns"] if by[g]["introns"] else 0
                                     for g in genomes], color=[CATEGORICAL[1 + sources.index(source[g]) % 5]
                                                               for g in genomes], height=0.6)
        for i, g in enumerate(genomes):
            b.text(0.01, i, f"{by[g]['introns']:,} introns", va="center", fontsize=7, color=SURFACE
                   if by[g]["introns"] and by[g]["exact_introns"] / by[g]["introns"] > 0.2 else INK)
        names = species_names(bench)
        b.set_yticks(range(len(genomes)), [names.get(g, g) for g in genomes], fontsize=7.5, fontstyle="italic")
        b.invert_yaxis()
        b.set_xlim(0, 1)
        pct_axis(b)
        b.grid(axis="x")
        b.set_title("Introns with both splice sites the annotation's")
    return save(fig, out, "annotation", fmt)


def introns(bench, package, out, fmt):
    rows = read(os.path.join(bench, "genes.tsv"))
    if not rows:
        return None
    by = collections.defaultdict(list)
    for r in rows:
        e, i = int(r["exon_bases"]), int(r["intron_bases"])
        if e + i:
            by[r["supergroup"]].append(i / (e + i))
    groups = group_order(by)
    fig, ax = plt.subplots(figsize=(7.5, 0.45 * len(groups) + 1.6))
    for y, g in enumerate(groups):
        vals = sorted(by[g])
        ax.scatter(vals, [y + ((k * 0.618) % 1 - 0.5) * 0.36 for k in range(len(vals))], s=16 if len(vals) < 300
                   else 4, color=SERIES[0], alpha=0.6 if len(vals) < 300 else 0.3, linewidths=0)
        ax.plot([vals[len(vals) // 2]] * 2, [y - 0.3, y + 0.3], color=INK, linewidth=2)
    ax.set_yticks(range(len(groups)), [f"{g} ({len(by[g]):,})" for g in groups])
    ax.invert_yaxis()
    ax.set_xlim(-0.02, 1.02)
    pct_axis(ax)
    ax.grid(axis="x")
    ax.set_title("Share of the BUSCO genes' bases in introns, each genome (bar: median)")
    return save(fig, out, "introns", fmt)


def resolution(bench, package, out, fmt):
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
    a.set_xticks(range(len(RANKS)), [f"{r}\n" + " | ".join(count(m, r) for m in markers) for r in RANKS])
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
        pct_axis(b, "y")
        b.grid(axis="y")
        b.set_title("Copies under 99% identical to any other species")
    else:
        b.set_axis_off()
    return save(fig, out, "resolution", fmt)


def barcode_gap(bench, package, out, fmt):
    """The barcode gap: per marker, each panel species' least alike pair of its own assemblies against its nearest
    other species; above the diagonal, told from it."""
    rows = [r for r in read(os.path.join(package, "intraspecific", "barcode_gap.tsv"))
            if num(r["intra_min"]) is not None and num(r["inter_max"]) is not None]
    if not rows:
        return None
    kinds = collections.defaultdict(list)
    for r in rows:
        kinds["BUSCOs" if r["marker"].endswith("at2759") else r["marker"]].append(r)
    names = ["BUSCOs"] + sorted(k for k in kinds if k != "BUSCOs") if "BUSCOs" in kinds else sorted(kinds)
    names = names[:6]
    fig, axes = plt.subplots(1, len(names), figsize=(3.3 * len(names), 3.6), squeeze=False)
    states = (("t", "told from its congener", SERIES[0]), ("f", "not told apart", SERIES[1]),
              ("-", "nearest other species not a congener", NONE))
    for ax, name in zip(axes[0], names):
        lo = min(min(num(r["intra_min"]), num(r["inter_max"])) for r in kinds[name])
        lo = max(0, min(97.0, lo - 1))
        for state, label, color in states:
            pts = [(num(r["inter_max"]), num(r["intra_min"])) for r in kinds[name] if r["resolved"] == state]
            if pts:
                ax.scatter([x for x, _y in pts], [y for _x, y in pts], s=10 if name == "BUSCOs" else 24, color=color,
                           alpha=0.7, linewidths=0, label=f"{label} ({len(pts):,})", zorder=2 if state != "-" else 1)
        ax.plot([lo, 100], [lo, 100], color=INK2, linewidth=0.8, linestyle=(0, (3, 2)))
        ax.set_xlim(lo, 100.3)
        ax.set_ylim(lo, 100.3)
        ax.set_aspect("equal")
        ax.grid(axis="both")
        ax.set_xlabel("identity to the nearest other species, %")
        ax.set_title(name, fontsize=10)
        ax.legend(loc="lower right", fontsize=6.5)
    axes[0][0].set_ylabel("least identity within the species, %")
    fig.suptitle("The barcode gap of the intraspecific panel (above the diagonal: the species' own assemblies nearer "
                 "to each other than to any other species)", x=0.01, ha="left", fontweight="bold", fontsize=10)
    return save(fig, out, "barcode_gap", fmt)


def markers(bench, package, out, fmt):
    rows = read(os.path.join(bench, "markers.tsv"))
    if not rows:
        return None
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        by[r["marker"]][r["verdict"]] += 1
    items = [(f"{m} ({sum(c.values())})", c) for m, c in sorted(by.items(), key=lambda kv: (
        -kv[1]["concordant"] / max(1, sum(kv[1].values())), kv[0]))]
    fig, ax = plt.subplots(figsize=(7.5, 0.36 * len(items) + 1.7))
    handles, labels = stacked(ax, items, VERDICTS, VERDICT_COLORS,
                              ("concordant (both ≥ 90%)", "similar (≥ 50%)", "discordant", "missed"))
    legend_below(ax, handles, labels)
    ax.set_title("Eukore's own marker loci against UniProt's protein of the same species (species)")
    return save(fig, out, "markers", fmt)


def flags(bench, package, out, fmt):
    rows = read(os.path.join(bench, "flags.tsv"))
    if not rows:
        return None
    combos = collections.Counter()
    tier = {}
    for r in rows:
        combos[r["flags"]] += int(r["copies"])
        tier[r["flags"]] = r["confidence"]
    top = combos.most_common(14)
    total = sum(combos.values())
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(top) + 1.5))
    ax.barh(range(len(top)), [n for _f, n in top], color=[TIER_COLORS.get(tier[f], NONE) for f, _n in top],
            height=0.65)
    for i, (_f, n) in enumerate(top):
        ax.text(n, i, f" {n:,} ({n / total:.0%})", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(top)), [("no flag" if f == "-" else f.replace(",", " + ")) for f, _n in top], fontsize=8)
    ax.invert_yaxis()
    ax.grid(axis="x")
    ax.set_xlabel("BUSCO copies")
    handles = [matplotlib.patches.Patch(color=c) for c in TIER_COLORS.values()]
    ax.legend(handles, list(TIER_COLORS), loc="lower right", title="tier", fontsize=7.5, title_fontsize=7.5)
    ax.set_title(f"The checks' flags and their combinations ({total:,} copies)")
    return save(fig, out, "flags", fmt)


def confidence(bench, package, out, fmt):
    rows = []
    copies = read(os.path.join(bench, "completeness.tsv"))
    if copies:
        rows.append(("BUSCO copies", {t: sum(int(r.get(t) or 0) for r in copies) for t in TIER_COLORS}))
    by = collections.defaultdict(collections.Counter)
    for r in read(os.path.join(bench, "rrna.tsv")):
        for t in TIER_COLORS:
            by[r["source"]][t] += int(r.get(t) or 0)
    rows += [(f"rRNA: {s}", c) for s, c in sorted(by.items(), key=lambda kv: -sum(kv[1].values()))]
    if not rows:
        return None
    rows = [(f"{n} ({sum(c.values()):,})", c) for n, c in rows]
    fig, ax = plt.subplots(figsize=(7.5, 0.42 * len(rows) + 1.5))
    handles, labels = stacked(ax, rows, tuple(TIER_COLORS), list(TIER_COLORS.values()))
    legend_below(ax, handles, labels, ncol=4)
    ax.set_title("Confidence tiers (docs/FORMAT.md: the checks that set them)")
    return save(fig, out, "confidence", fmt)


def rrna(bench, package, out, fmt):
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


# rrna_lengths' panels: (source, region, title); the organelle genes' from organelle/genes.tsv (those not partial)
LENGTH_PANELS = (("genome", "18S", "18S"), ("genome", "ITS1", "ITS1"), ("genome", "5.8S", "5.8S"),
                 ("genome", "ITS2", "ITS2"), ("genome", "28S", "28S"),
                 ("mitochondrion", "SSU", "mito. SSU\n(12S in animals)"),
                 ("mitochondrion", "COX1", "mito. COX1\n(with its introns)"))


def rrna_lengths(bench, package, out, fmt):
    rows = read(os.path.join(bench, "rrna_regions.tsv"))
    at = {(r["supergroup"], r["source"], r["region"]): r for r in rows}
    panels = [p for p in LENGTH_PANELS if any((r["source"], r["region"]) == p[:2] for r in rows)]
    if not panels:
        return None
    groups = group_order({g for g, s, x in at if (s, x) in {p[:2] for p in panels}})
    fig, axes = plt.subplots(1, len(panels), figsize=(2.3 * len(panels), 0.4 * len(groups) + 1.8), sharey=True,
                             squeeze=False)
    axes = axes[0]
    for ax, (source, reg, title) in zip(axes, panels):
        color = SERIES[0] if source == "genome" else SERIES[2]
        stats, pos = [], []
        for y, g in enumerate(groups):
            r = at.get((g, source, reg))
            if r:
                stats.append(dict(med=num(r["median"]), q1=num(r["q25"]), q3=num(r["q75"]), whislo=num(r["q05"]),
                                  whishi=num(r["q95"]), label=g))
                pos.append(y)
                ax.text(1.0, y, f" {r['records']}", transform=ax.get_yaxis_transform(), fontsize=6.5, color=INK2,
                        va="center", ha="left")
        if stats:
            ax.bxp(stats, positions=pos, orientation="horizontal", widths=0.55, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=color, edgecolor=color, alpha=0.85),
                   medianprops=dict(color=SURFACE, linewidth=2), whiskerprops=dict(color=color),
                   capprops=dict(color=color))
            # the median as a dot too: a group of a record or two, or of one length, draws no box
            ax.scatter([s["med"] for s in stats], pos, s=18, color=color, edgecolor=INK, linewidth=0.6, zorder=3)
        ax.set_title(title, fontsize=9.5)
        ax.grid(axis="x")
        ax.set_xlabel("bases")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
    axes[0].set_yticks(range(len(groups)), groups)
    axes[0].invert_yaxis()
    fig.suptitle("Lengths per supergroup: the genomes' rRNA operons' regions (blue) and the mitochondrial genes "
                 "(green); records at the right", x=0.01, y=1.08, ha="left", fontweight="bold", fontsize=10.5)
    fig.subplots_adjust(wspace=0.35)
    return save(fig, out, "rrna_lengths", fmt)


def rrna_copies(bench, package, out, fmt):
    rows = [r for r in read(os.path.join(bench, "rrna_copies.tsv")) if str(r["rrna_18s"]).isdigit()]
    if not rows:
        return None
    groups = group_order({r["supergroup"] for r in rows})
    chromosome = ("Complete Genome", "Chromosome")
    rng = random.Random(1)
    small = len(rows) < 2000
    fig, axes = plt.subplots(1, 2, figsize=(10, 0.5 * len(groups) + 1.8), sharey=True)
    for ax, (col, title) in zip(axes, (("rrna_18s", "18S copies"), ("rdna_units", "rDNA units (18S to 28S)"))):
        for whole, color, label in ((True, SERIES[0], "chromosome-level assembly"),
                                    (False, SERIES[1], "scaffolds or contigs")):
            pts = [(int(r[col]), groups.index(r["supergroup"]) + rng.uniform(-0.18, 0.18)) for r in rows
                   if (r["level"] in chromosome) == whole and str(r[col]).isdigit()]
            if pts:
                ax.scatter([x for x, _y in pts], [y for _x, y in pts], s=22 if small else 5, color=color,
                           alpha=0.8 if small else 0.35, label=label, edgecolor=SURFACE, linewidth=0.4 if small else 0)
        for y, g in enumerate(groups):
            vals = sorted(int(r[col]) for r in rows if r["supergroup"] == g and str(r[col]).isdigit())
            if vals:
                ax.plot([vals[len(vals) // 2]] * 2, [y - 0.3, y + 0.3], color=INK, linewidth=1.5)
        ax.set_xscale("symlog", linthresh=1)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v:,.0f}"))
        ax.set_xlim(left=-0.3)
        ax.grid(axis="x")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("copies in the assembly (per genome; bar: median)")
    axes[0].set_yticks(range(len(groups)), [f"{g} ({sum(r['supergroup'] == g for r in rows):,})" for g in groups])
    axes[0].invert_yaxis()
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2,
               frameon=False, fontsize=8.5)
    fig.suptitle("rRNA copies per genome, per supergroup (as assembled: tandem arrays collapse)", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    return save(fig, out, "rrna_copies", fmt)


def primers(bench, package, out, fmt):
    rows = read(os.path.join(bench, "primers.tsv"))
    if not rows:
        return None
    pairs = sorted({r["pair"] for r in rows})
    groups = ["all"] + group_order({r["supergroup"] for r in rows if r["supergroup"] != "all"})
    at = {(r["pair"], r["supergroup"]): r for r in rows}
    matrix = [[num(at[(p, g)]["share"], 0) if (p, g) in at else float("nan") for g in groups] for p in pairs]
    fig, ax = plt.subplots(figsize=(1.0 * len(groups) + 3, 0.36 * len(pairs) + 1.8))
    ax.imshow(matrix, aspect="auto", cmap=SEQUENTIAL, vmin=0, vmax=1)
    for i, p in enumerate(pairs):
        for j, g in enumerate(groups):
            r = at.get((p, g))
            if r:
                v = num(r["share"], 0)
                ax.text(j, i, f"{v:.0%}\n{int(r['amplified']):,}/{int(r['species']):,}", ha="center", va="center",
                        fontsize=6.5, color=SURFACE if v > 0.55 else INK)
    ax.set_xticks(range(len(groups)), groups, rotation=30, ha="left", rotation_mode="anchor")
    ax.xaxis.tick_top()
    ax.set_yticks(range(len(pairs)), pairs)
    ax.axvline(0.5, color=SURFACE, linewidth=3)
    ax.tick_params(length=0)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title("In-silico PCR: the species a primer pair amplifies, of those with a record it runs on", pad=40)
    return save(fig, out, "primers", fmt)


def organelle(bench, package, out, fmt):
    rows = read(os.path.join(bench, "organelle.tsv"))
    if not rows:
        return None
    genes = sorted({(r["organelle"], r["gene"]) for r in rows})
    groups = group_order({r["supergroup"] for r in rows})
    at = {(r["supergroup"], r["organelle"], r["gene"]): r for r in rows}
    matrix = [[int(at[(g, o, x)]["with_gene"]) / max(1, int(at[(g, o, x)]["species"])) if (g, o, x) in at else 0
               for o, x in genes] for g in groups]
    fig, ax = plt.subplots(figsize=(0.7 * len(genes) + 3, 0.42 * len(groups) + 1.6))
    ax.imshow(matrix, aspect="auto", cmap=SEQUENTIAL, vmin=0, vmax=1)
    for i, line in enumerate(matrix):
        for j, v in enumerate(line):
            ax.text(j, i, f"{v:.0%}", ha="center", va="center", fontsize=7, color=SURFACE if v > 0.55 else INK)
    ax.set_xticks(range(len(genes)), [f"{o[:4]}.\n{x}" for o, x in genes], fontsize=7)
    ax.set_yticks(range(len(groups)), groups)
    ax.tick_params(length=0)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title("Species with each organelle gene")
    return save(fig, out, "organelle", fmt)


def labels(bench, package, out, fmt):
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


def read_metrics(bench):
    """{reference: {metric: value}} of the marker reads, and background assignments per million."""
    sums = collections.defaultdict(collections.Counter)
    for r in read(os.path.join(bench, "reads.tsv")):
        if r["mode"] == "held_out" and not int(r["congeners"]):
            continue  # no other species of its genus in the package: its held-out reads cannot be placed to it
        key = (r["reference"], r["mode"], r["kind"], r["reads"])
        for c in r:
            if c in ("n", "assigned") or c.endswith(("_correct", "_wrong")):
                sums[key][c] += int(r[c])
    out = collections.defaultdict(dict)
    for (ref, mode, kind, reads), s in sums.items():
        if kind == "eukaryote" and reads == "marker" and s["n"]:
            for rank in ("species", "genus"):
                right, wrong = s[f"{rank}_correct"], s[f"{rank}_wrong"]
                out[ref][(mode, rank, "recall")] = right / s["n"]
                out[ref][(mode, rank, "precision")] = right / (right + wrong) if right + wrong else None
        elif kind == "background" and s["n"]:
            out[ref]["background"] = 1e6 * s["assigned"] / s["n"]
            out[ref]["background_phylum"] = 1e6 * s["phylum_wrong"] / s["n"]
    return out


READ_METRICS = ((("in", "species", "recall"), "species recall, species in the reference"),
                (("held_out", "genus", "recall"), "genus recall, species held out"),
                (("held_out", "genus", "precision"), "genus precision, species held out"))


def marker_reads(bench, package, out, fmt):
    m = read_metrics(bench)
    if not m:
        return None
    fig, ax = plt.subplots(figsize=(7.5, 5))
    refs = [r for r in EUKORE_REFS if r in m] + [r for r in compared_refs(bench) if r in m]
    w = 0.26
    for i, (key, label) in enumerate(READ_METRICS):
        xs = [j + (i - 1) * w for j in range(len(refs))]
        ax.bar(xs, [m[r].get(key) or 0 for r in refs], width=w * 0.92, color=[SERIES[0], BLUE[250], SERIES[2]][i],
               label=label)
    for j, r in enumerate(refs):
        bg = m[r].get("background_phylum")
        if bg is not None:
            ax.text(j, 1.03, f"{bg:,.0f}", ha="center", va="bottom", fontsize=7, color=INK2)
    ax.text(-0.5, 1.10, "background reads per million placed in a eukaryote phylum:", fontsize=7, color=INK2)
    ax.set_xticks(range(len(refs)), [tick_label(r) for r in refs])
    ax.set_ylim(0, 1.16)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    pct_axis(ax, "y")
    ax.grid(axis="y")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=7.5, ncol=1)
    ax.set_title("Simulated reads from Eukore's markers' loci (Eukore's own genes: other references hold others)",
                 fontsize=9.5)
    return save(fig, out, "marker_reads", fmt)


def detection_depth(bench, package, out, fmt):
    """Detection's F1, precision and recall against the sample's depth, per reference and mode (each a row)."""
    det = detection(bench)
    refs = detection_refs(bench, det)
    if not refs:
        return None
    others = [r for r in refs if r not in EUKORE_REFS]
    modes = [m for m in DETECT_MODES if any(split_key(k)[0] == m[0] for r in refs for k in det[r])]
    fig, axes = plt.subplots(len(modes), 3, figsize=(12, 3.2 * len(modes)), sharex=True, sharey=True, squeeze=False)
    for row, (mode, label, _c) in zip(axes, modes):
        for ax, what in zip(row, ("f1", "precision", "recall")):
            for ref in refs:
                pts = sorted((split_key(k)[1], v.get(what)) for k, v in det[ref].items()
                             if split_key(k)[0] == mode and v.get(what) is not None)
                if pts:
                    ax.plot([d for d, _v in pts], [v for _d, v in pts], markersize=5, linewidth=2, label=ref,
                            **ref_style(ref, others))
            ax.set_xscale("log")
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _p: f"{v:,.0f}"))
            ax.set_ylim(-0.03, 1.05)
            ax.grid(axis="y")
            ax.set_title(f"{what if what != 'f1' else 'F1'}: {label.split(',')[0]}", fontsize=9.5)
        row[0].set_ylabel(label.split(", ", 1)[1], fontsize=8.5)
        cens = censored_share(bench, mode) if mode != "in" else None
        if cens is not None:
            row[2].text(0.98, 0.04, f"Eukore's reads censored: {cens:.2%}", transform=row[2].transAxes, ha="right",
                        fontsize=7.5, color=INK2)
    for ax in axes[-1]:
        ax.set_xlabel("reads of the genome in its sample")
    fig.legend(*axes[0][0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.0),
               ncol=len(refs), frameon=False, fontsize=8.5)
    fig.suptitle(f"Detecting each test genome's taxon by sequencing depth (a taxon detected with "
                 f"{MIN_READS_LABEL} reads; as many bacterial reads in each sample)", x=0.01, ha="left",
                 fontweight="bold", fontsize=10.5)
    fig.tight_layout()
    return save(fig, out, "detection_depth", fmt)


def detection_genomes(bench, package, out, fmt):
    """Each test genome's reads placed to its right taxon per reference, at the deepest samples; gray where the
    reference does not hold it."""
    rows = read(os.path.join(bench, "detection.tsv"))
    if not rows:
        return None
    depth = max(int(r["depth"]) for r in rows)
    rows = [r for r in rows if int(r["depth"]) == depth]
    det = detection(bench)
    refs = detection_refs(bench, det) or sorted({r["reference"] for r in rows})
    names = species_names(bench)
    fig, axes = plt.subplots(1, 2, figsize=(1.1 * len(refs) * 2 + 4, 0.38 * len({r["test"] for r in rows}) + 2.2),
                             sharey=True)
    for ax, (mode, label, _c) in zip(axes, DETECT_MODES):
        mine = [r for r in rows if r["mode"] == mode and (mode == "in" or int(r["congeners"]))]
        tests = sorted({(r["supergroup"], r["test"]) for r in rows if r["mode"] == mode})
        at = {(r["test"], r["reference"]): r for r in mine}
        for y, (_g, t) in enumerate(tests):
            for x, ref in enumerate(refs):
                r = at.get((t, ref))
                if r is None:
                    color, text = SURFACE, "-"
                elif r["in_reference"] != "t":
                    color, text = NONE, "not in it"
                else:
                    color = BLUE[550] if r["detected"] == "t" else SERIES[1]
                    text = f"{int(r['true_reads']):,}"
                ax.add_patch(matplotlib.patches.Rectangle((x - 0.46, y - 0.42), 0.92, 0.84, color=color))
                ax.text(x, y, text, ha="center", va="center", fontsize=7,
                        color=SURFACE if color in (BLUE[550], SERIES[1]) else INK)
        ax.set_xlim(-0.5, len(refs) - 0.5)
        ax.set_ylim(len(tests) - 0.5, -0.5)
        ax.set_xticks(range(len(refs)), [tick_label(r) for r in refs], fontsize=8)
        ax.set_yticks(range(len(tests)), [f"{names.get(t, t)} ({g})" for g, t in tests], fontsize=7.5)
        ax.tick_params(length=0)
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(False)
        ax.set_title(f"the right {label}", fontsize=9.5)
    handles = [matplotlib.patches.Patch(color=BLUE[550], label="detected"),
               matplotlib.patches.Patch(color=SERIES[1], label="in the reference, not detected"),
               matplotlib.patches.Patch(color=NONE, label="not in the reference")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle(f"Reads placed to each test genome's right taxon, samples of {depth:,} reads (held out: genomes "
                 f"with congeners)", x=0.01, ha="left", fontweight="bold", fontsize=10.5)
    fig.tight_layout()
    return save(fig, out, "detection_genomes", fmt)


def runtime(bench, package, out, fmt):
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
        ax.scatter(x, y, s=24 if len(rows) < 500 else 5, color=SERIES[0], alpha=0.65 if len(rows) < 500 else 0.3,
                   linewidths=0.8 if len(rows) < 500 else 0, edgecolors=SURFACE)
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
        ax.set_title(f"Each genome's run ({len(rows):,}; {'/'.join(threads)} threads)")
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


# --- all of them -----------------------------------------------------------------------------------------------------

# (name, function, caption, its tables) in their order: Figure 1, then the supplementary figures S1 ..
FIGURES = (
    ("Fig1", figure1,
     "(a) Species per supergroup: NCBI's with a genome (gray track), Eukore's with a genome (dot) and with any record "
     "(ring), and the other references' (genome markers orange, 18S references green). (b) Left: BUSCO's dataset "
     "species' genomes called again without their genus', family's and order's dataset species, by the rank they "
     "share with the nearest dataset species kept (own: queried with their own proteins); the share of OrthoDB's genes "
     "called, on the gene the own-protein call found, and OrthoDB's gene model (5-mers ≥ 90% both ways); genomes at "
     "the foot. Right: each genome's calls laid on its own annotation, the share of its introns whose splice sites "
     "are the annotation's, at its distance to BUSCO's dataset. (c) Per marker: left, the copies under 99% identical "
     "to their nearest congener; right, the intraspecific panel's species whose own assemblies are nearer to each "
     "other than to any congener (each BUSCO a dot; bar: all BUSCOs pooled). (d) Each test genome a sample of its "
     "reads from anywhere in it, with as many bacterial reads; a taxon detected with 3 reads or more and 1% or more of "
     "those placed: F1 of detecting the right species (its species in the reference) and genus (held out) by depth, "
     "per reference.",
     "breadth, overlap, holdout, orthodb, annotation, resolution, intraspecific/barcode_gap, detection, reads"),
    ("selection", selection, "The genome set: what became of NCBI's taxa with an assembly, per supergroup.",
     "selection"),
    ("overlap", overlap, "Species in Eukore (any record) and in each other database the build read (PR2, EUKARYOME's "
                         "long reads, KSGP, ROD, BUSCO's OrthoDB dataset, the compared references), per supergroup.",
     "overlap"),
    ("completeness", completeness, "BUSCO completeness per supergroup and per genome.", "completeness"),
    ("completeness_vs_assembly", completeness_vs_assembly, "Completeness against the assembly's contiguity.",
     "completeness"),
    ("loci_call_rates", loci_call_rates, "Each BUSCO's call rate per supergroup (the loci a user chooses among).",
     "busco/loci.tsv"),
    ("orthodb", orthodb, "Calls against OrthoDB's gene models, with and without the species' own proteins.",
     "orthodb"),
    ("holdout", holdout, "The clade hold-out: each dataset species' genome called again without its genus', family's "
                         "and order's dataset species; per genome (thin) and pooled (bold), by the rank shared with "
                         "the nearest dataset species kept.", "holdout, orthodb"),
    ("annotation", annotation, "Calls laid on each genome's own annotation (NCBI's, Ensembl's, TOGA2's): coding bases "
                               "and splice sites, each genome.", "annotation"),
    ("introns", introns, "Intron share of the BUSCO genes' bases: what an exon-only reference leaves out.", "genes"),
    ("resolution", resolution, "Identity of each copy's nearest other species by the rank they share; BUSCOs and 18S.",
     "resolution"),
    ("barcode_gap", barcode_gap, "The intraspecific panel's barcode gap per marker: each species' least alike pair of "
                                 "its own assemblies against its nearest other species.", "intraspecific/barcode_gap"),
    ("markers", markers, "Eukore's own marker loci (RPB1 .. H3) against UniProt's protein of the same species.",
     "markers"),
    ("flags", flags, "Which checks fire, and together.", "flags"),
    ("confidence", confidence, "Confidence tiers of the BUSCO copies and the rRNA records.", "completeness, rrna"),
    ("rrna", rrna, "rRNA records holding each region, per source.", "rrna"),
    ("rrna_lengths", rrna_lengths, "rRNA regions' lengths in the genomes' operons and the mitochondrial SSU (12S in "
                                   "animals) and COX1 (introns included) per supergroup.", "rrna_regions"),
    ("rrna_copies", rrna_copies, "rRNA copies per genome as assembled (18S spans and whole rDNA units), per supergroup "
                                 "and assembly level: tandem arrays collapse in an assembly, so these are not the "
                                 "genomes' copy numbers.", "rrna_copies"),
    ("primers", primers, "In-silico PCR: per primer pair and supergroup, the species it amplifies of those with a "
                         "record it runs on.", "primers"),
    ("organelle", organelle, "Species with each organelle gene.", "organelle"),
    ("labels", labels, "Species named in PR2, EUKARYOME, KSGP, SILVA, UNITE and MIDORI2.", "labels"),
    ("marker_reads", marker_reads, "Reads simulated from Eukore's own markers' loci (BUSCO genes and rRNA) per "
                                   "reference: species recall (species in), genus recall and precision (held out). "
                                   "Eukore's genes, so the other references, which hold other genes, are measured on "
                                   "ground that is not theirs.", "reads"),
    ("detection_depth", detection_depth, "Detection's F1, precision and recall by the sample's depth, per reference "
                                         "and mode: the right species with it in the reference, the right genus with "
                                         "it held out, the right class with its order held out (novel); the reads "
                                         "minimap2's secondary list censored.", "detection, reads"),
    ("detection_genomes", detection_genomes, "Each test genome's reads placed to its right taxon per reference in the "
                                             "deepest samples: detected, held but missed, or not held.", "detection"),
    ("runtime", runtime, "Run time per genome and of the package stage.", "runtime, stages"),
)
NOT_COLLECTED = (
    "read assignment by a read profiler (protal) on held-out genomes: the end-to-end measure, outside the build",
    "the caller's gates against compleasm's on divergent genomes (a panel run, not a build's benchmark)",
    "the BUSCOs' prokaryote homologues (tools/prokaryote_panel.py, tools/prokaryote_separation.py)",
    "a contamination spike-in (to be designed: a eukaryote contaminant's copy where the native one is removed)",
)


def index_md(package, made, out):
    lines = [f"# Figures: {os.path.basename(os.path.normpath(package))}", "",
             "Drawn by `eukore/figures.py` from the package's `benchmarks/` (docs/benchmarks.md).", ""]
    number = 0
    for name, _fn, caption, tables in FIGURES:
        title = "Figure 1" if name == "Fig1" else f"S{number}: {name}"
        number += 1
        path = made.get(name)
        lines += [f"## {title}", "", caption, "", f"Tables: {tables}.", ""]
        lines += [f"![{name}]({os.path.basename(path)})", ""] if path else ["*Not drawn: its tables are missing.*", ""]
    lines += ["## Not measured by the build", ""] + [f"- {x}" for x in NOT_COLLECTED] + [""]
    path = os.path.join(out, "index.md")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(lines))
    return path


def draw(package, bench=None, out=None, fmt="png", log=print):
    """Every figure of FIGURES of a package (its benchmarks in bench, default PACKAGE/benchmarks) into out (default
    bench/figures), with index.md and summary.md -> {name: path, or None where left out}."""
    bench = bench or os.path.join(package, benchmarks.FOLDER)
    out = out or os.path.join(bench, FOLDER)  # OUT_BENCH/figures, or PACKAGE/benchmarks/figures
    os.makedirs(out, exist_ok=True)
    made = {}
    for name, fn, _caption, _tables in FIGURES:
        try:
            made[name] = fn(bench, package, out, fmt)
        except Exception as e:  # noqa: BLE001 (a figure's failure costs the others nothing)
            plt.close("all")
            log(f"WARNING: figure {name} failed ({type(e).__name__}: {e}): left out")
            made[name] = None
    made["index"] = index_md(package, made, out)
    summary = load_summary(bench)
    if summary:
        manifest = {}
        if os.path.isfile(os.path.join(package, "manifest.json")):
            with open(os.path.join(package, "manifest.json")) as fh:
                manifest = json.load(fh)
        partial = (manifest.get("inputs") or {}).get("partial")
        made["summary"] = os.path.join(out, "summary.md")
        with open(made["summary"], "w", newline="\n") as fh:
            fh.write(benchmarks.summary_md(summary, os.path.basename(os.path.normpath(package)),
                                           len(partial["species_clades"]) if partial else None,
                                           read(os.path.join(bench, "composition.tsv"))))
    return made


def main(argv=None, default_out=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("package", help="a Eukore package folder (Eukore_<version>/)")
    p.add_argument("--benchmarks", help="the benchmarks' tables (default PACKAGE/benchmarks)")
    p.add_argument("--out", default=default_out, help="where the figures go (default "
                                                      + (default_out or "PACKAGE/benchmarks/figures") + ")")
    p.add_argument("--format", default="png", choices=("png", "svg", "pdf"))
    opts = p.parse_args(argv)
    bench = opts.benchmarks or os.path.join(opts.package, benchmarks.FOLDER)
    if not os.path.isdir(bench):
        p.error(f"{bench}: no benchmarks (python3 -m eukore_bench.compare makes them)")
    made = draw(opts.package, bench, opts.out, opts.format, log=lambda t: print(t, file=sys.stderr))
    for path in filter(None, made.values()):
        print(path)
    left = [name for name, path in made.items() if not path]
    if left:
        print(f"left out (no table): {', '.join(left)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
