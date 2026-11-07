#!/usr/bin/env python3
"""The Eukore overview figure: (a) the build, from its sources to the package, and (b) the package, what each of its
parts makes available. Drawn with the paper's palette (docs/benchmarks/plot_benchmarks.py).

    python3 docs/figures/overview.py                              # docs/figures/overview.png, the numbers below
    python3 docs/figures/overview.py --package OUT/Eukore_2026.11 # the package's counts (its manifest) head panel b
    python3 docs/figures/overview.py --format svg --no_marks      # without the "added in 0.2" marks

The numbers are NUMBERS below (NCBI's tables and the sources' releases of 2026-11); an item that starts with "*" is a
part added in Eukore 0.2 (format 5) and gets a mark. After the first full build, run it with --package so that panel b
carries the package's own counts, and update NUMBERS where the sources moved.
"""

import argparse
import importlib.util
import json
import os
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("plot_benchmarks", os.path.join(HERE, "..", "benchmarks",
                                                                              "plot_benchmarks.py"))
pb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pb)
plt = pb.plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch  # noqa: E402 (after plot_benchmarks' rcParams)

INK, INK2, MUTED, GRID, SURFACE, BLUE = pb.INK, pb.INK2, pb.MUTED, pb.GRID, pb.SURFACE, pb.BLUE
NEW = pb.SERIES[1]            # the mark of a part added in 0.2
SOURCE_FILL = "#f3f7fd"       # the sources' and the package's boxes
STAGE_FILL = SURFACE

NUMBERS = {
    "ncbi_species": "~29,000", "ncbi_tb": "~6.5 TB", "toga2_species": "72", "tsa_species": "5,025",
    "eukprot_species": "75", "veupathdb_species": "71", "ensembl_species": "2,037", "toga2_models": "907",
    "buscos": "129", "primer_pairs": "18", "refseq_organelle": "~2.2 GB", "fcs_gx_db": "~470 GB",
    "sources": "21", "panel_species": "40", "panel_assemblies": "3", "benchmark_tables": "20",
    "rrna_primer_pairs": "10", "organelle_primer_pairs": "8", "tetranucleotides": "136", "references": "6",
    "id_systems": "7", "ncbi_date": "2026-11",
}

# layout, in mm
WIDTH = 180
PAD, LINE, TITLE, GAP, INDENT = 1.6, 2.6, 3.6, 2.2, 2.0
FS_BODY, FS_TITLE, FS_HEAD, FS_PANEL = 6.0, 6.8, 8.0, 10
EM = 0.54 * FS_BODY * 0.3528    # DejaVu Sans: a character's width in mm, on average


def fmt(text, n):
    return text.format(**n)


def wrap(items, w):
    """Each item (a string; "*" first marks it) as its wrapped lines [(text, marked, first)]."""
    lines = []
    for item in items:
        marked = item.startswith("*")
        text = item[1:].strip() if marked else item
        width = int((w - 2 * PAD - (INDENT if marked else 0)) / EM)
        parts = textwrap.wrap(text, width) or [""]
        lines += [(part, marked, i == 0) for i, part in enumerate(parts)]
    return lines


class Box:
    def __init__(self, title, items, w, min_h=0):
        self.title, self.w = title, w
        self.marked = title.startswith("*")
        if self.marked:
            self.title = title[1:].strip()
        self.lines = wrap(items, w)
        self.h = max(PAD * 2 + TITLE + LINE * len(self.lines), min_h)

    def draw(self, ax, x, y, fill, marks):
        """The box with its top-left corner at (x, y)."""
        ax.add_patch(FancyBboxPatch((x, y - self.h), self.w, self.h, boxstyle="round,pad=0,rounding_size=1.2",
                                    facecolor=fill, edgecolor=GRID, linewidth=0.7))
        tx = x + PAD
        if self.marked and marks:
            dot(ax, x + PAD + 0.6, y - PAD - TITLE * 0.55)
            tx += INDENT
        ax.text(tx, y - PAD - TITLE * 0.55, self.title, fontsize=FS_TITLE, fontweight="bold", color=INK,
                va="center", ha="left")
        yy = y - PAD - TITLE
        for text, marked, first in self.lines:
            indent = INDENT if (marked and marks) else 0
            if marked and marks and first:
                dot(ax, x + PAD + 0.6, yy - LINE * 0.55)
            ax.text(x + PAD + indent, yy - LINE * 0.55, text if marks else text.replace("*", ""), fontsize=FS_BODY,
                    color=INK, va="center", ha="left")
            yy -= LINE
        return y - self.h


def dot(ax, x, y):
    ax.plot([x], [y], marker="o", markersize=2.2, color=NEW, markeredgewidth=0, clip_on=False)


def column_height(boxes, gap=GAP):
    return sum(b.h for b in boxes) + gap * (len(boxes) - 1)


def column(ax, x, y, boxes, fill, marks, gap=GAP):
    """Boxes stacked from the top y down; returns [(top, bottom)] of each."""
    spans = []
    for b in boxes:
        bottom = b.draw(ax, x, y, fill, marks)
        spans.append((y, bottom))
        y = bottom - gap
    return spans


def head(ax, x, y, text, sub=None):
    ax.text(x, y, text, fontsize=FS_HEAD, fontweight="bold", color=INK, va="bottom", ha="left")
    if sub:
        ax.text(x, y - 3.0, sub, fontsize=FS_BODY, color=INK2, va="bottom", ha="left")


def arrow(ax, p, q, label=None, lw=0.9, color=INK2):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=7, linewidth=lw, color=color,
                                 shrinkA=0, shrinkB=0, clip_on=False))
    if label:
        ax.text((p[0] + q[0]) / 2, (p[1] + q[1]) / 2 + 1.0, label, fontsize=FS_BODY - 0.4, color=INK2,
                ha="center", va="bottom", bbox=dict(facecolor=SURFACE, edgecolor="none", pad=0.6))


def panel_a_content(n):
    xs, ws = 1, 50          # sources
    xf, wf = 55, 24         # fetch
    xp, wp = 83, 45         # process
    xk, wk = 132, 47        # package
    sources = [Box(t, items, ws) for t, items in [
        ("genomes", [
            fmt("NCBI assemblies (GenBank, RefSeq): {ncbi_species} species, {ncbi_tb}", n),
            fmt("TOGA2 assemblies: {toga2_species} vertebrates with no NCBI genome", n),
            fmt("* INSDC transcriptomes: {tsa_species} species with no genome", n),
            fmt("EukProt records: {eukprot_species} protists with no assembly", n)]),
        ("the genomes' own gene models", [
            fmt("NCBI · VEuPathDB ({veupathdb_species} sp.) · Ensembl ({ensembl_species} sp.) · "
                "TOGA2 ({toga2_models} sp.)", n)]),
        ("marker models", [
            fmt("BUSCO eukaryota_odb12: {buscos} HMMs and their proteins", n),
            fmt("RefSeq organelle genomes ({refseq_organelle})", n),
            fmt("* Rfam RF01960, the eukaryote SSU  ·  {primer_pairs} primer pairs", n)]),
        ("published rDNA", [
            "EUKARYOME 2.1 long reads  ·  ROD 1.2 operons"]),
        ("taxonomy and names", [
            "* NCBI taxdump, its nuclear and organelle codes",
            "PR2 5.1.1 · EUKARYOME 2.1 · KSGP 4.0 · SILVA 138.2* · UNITE 10.0* · MIDORI2*",
            "* Open Tree 3.7.3 (with GBIF, WoRMS, IRMNG, Index Fungorum ids) · BOLD BINs"]),
        ("quality screens", [
            fmt("NCBI's FCS-GX and ANI reports · FCS-GX database ({fcs_gx_db}) · tiara", n)]),
    ]]
    fetch = [Box("genome set", [
        "one best assembly per species: NCBI's flags, FCS-GX and ANI exclude; then RefSeq reference, "
        "representative, level, contiguity",
        fmt("{ncbi_tb} of genomes, pruned as it goes", n),
        fmt("* {panel_species} species' {panel_assemblies} further assemblies each", n),
        "every stage resumable"], wf, min_h=column_height(sources))]
    process = [Box(t, items, wp) for t, items in [
        ("BUSCO calls", [
            "miniprot aligns the 30 nearest proteins with the genome's genetic code; hmmsearch confirms",
            "rescue pass: each Missing HMM in the genome's ORFs",
            "the genome's own gene model where its HMM confirms it",
            "the gene kept whole: exons, introns, frameshifts"]),
        ("rDNA", [
            "barrnap, ITSx: 18S spans and whole units, 18S–ITS1–5.8S–ITS2–28S (an inline 5S)"]),
        ("organelle genes", [
            "RefSeq's record, else the assembly's organelle molecules: nhmmer, barrnap",
            "* miniprot for exact exons (an option)"]),
        ("assembly metadata", [
            "one pass over the genome while it is on disk; sourmash sketch; organelle molecules kept"]),
        ("contamination windows", [
            "50 kb around each copy: tiara, FCS-GX"]),
    ]]
    package = [Box(t, items, wk) for t, items in [
        ("checks across species", [
            "own-species agreement, 99% collapses, length outliers, rDNA chimera screen",
            "* a tier high / medium / low, with its reasons, on every copy, record and organelle gene"]),
        ("* genetics", [
            "hmmalign: intron phase, splice sites, U12, sites shared across species",
            "coding sequences and proteins; genetic-code check; duplication: haplotypes or a WGD"]),
        ("* rDNA regions and PCR", [
            "cmalign on RF01960: the 18S's V1–V9",
            fmt("in-silico PCR of {primer_pairs} primer pairs", n)]),
        ("names and ids", [
            fmt("lineage in {references} references, ids in {id_systems} systems", n),
            "* supergroup and group; empty ranks filled"]),
        ("* barcode gap", [
            "each marker's variation within a species beside its nearest other species"]),
        ("benchmarks, manifest", [
            fmt("{benchmark_tables} tables measured on the package; the manifest written last", n)]),
    ]]
    cols = [(xs, ws, sources, SOURCE_FILL), (xf, wf, fetch, STAGE_FILL), (xp, wp, process, STAGE_FILL),
            (xk, wk, package, STAGE_FILL)]
    return cols


def panel_a(ax, n, marks):
    """The build: sources -> fetch -> process (per genome) -> package (across species)."""
    cols = panel_a_content(n)
    top_pad = 13
    H = top_pad + 3.6 + max(column_height(b) for _x, _w, b, _f in cols) + 4
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    top = H - top_pad
    ax.text(0, H - 1, "a", fontsize=FS_PANEL, fontweight="bold", va="top", ha="left")
    ax.text(5, H - 1, "The build: every eukaryote genome NCBI holds, to one package", fontsize=FS_HEAD,
            fontweight="bold", va="top", ha="left")
    heads = [(fmt("sources: {sources}", n), fmt("NCBI's tables, the references of {ncbi_date}", n)),
             ("fetch", "the network stage"), ("process", "each genome once, in parallel"),
             ("package", "all species together")]
    spans = []
    for (x, w, boxes, fill), (title, sub) in zip(cols, heads):
        head(ax, x, top + 1.0, title, sub)
        spans.append(column(ax, x, top - 3.6, boxes, fill, marks))
    s_spans, f_spans, p_spans, k_spans = spans
    (xs, ws, _, _), (xf, wf, _, _), (xp, wp, _, _), (xk, _, _, _) = cols
    f_top, f_bot = f_spans[0]
    for t, b in s_spans:
        arrow(ax, (xs + ws, (t + b) / 2), (xf - 0.3, (t + b) / 2))
    p_mid = (p_spans[0][0] + p_spans[-1][1]) / 2
    arrow(ax, (xf + wf, (f_top + f_bot) / 2), (xp - 0.3, (f_top + f_bot) / 2))
    arrow(ax, (xp + wp, p_mid), (xk - 0.3, p_mid))
    if marks:
        dot(ax, WIDTH - 46, H - 2.6)
        ax.text(WIDTH - 44, H - 2.6, "and *: added in Eukore 0.2 (package format 5)", fontsize=FS_BODY,
                color=INK2, va="center", ha="left")
    return H


def panel_b_content(n):
    w = 42.6
    row1 = [Box(t, items, w) for t, items in [
        ("taxonomy and names", [
            "species.tsv: NCBI lineage, its empty ranks filled; lineage in PR2, EUKARYOME, KSGP, SILVA*, "
            "UNITE*, MIDORI2*",
            "* supergroup and group of every species",
            "* crossrefs.tsv: Open Tree, GBIF, WoRMS, IRMNG, Index Fungorum, UNITE SH, BOLD BIN",
            "genomes.tsv: accession, flags, rRNA copies, genetic code*, duplication*, type material*",
            "* taxonomy/: NCBI's, cut to the package, 3 codes per node"]),
        ("busco/", [
            fmt("{buscos} BUSCOs × every species", n),
            "copies/: each gene whole, introns and all",
            "* cds/, proteins/: by the genome's code",
            "features.tsv: exons, introns, frameshifts, masked runs, assembly gaps",
            "* introns.tsv, intron_sites.tsv: each intron's phase, sites, class; the sites shared across species",
            "copies.tsv: each copy's call, checks, confidence",
            "loci.tsv: each BUSCO's call rate per group, to choose loci by"]),
        ("rrna/", [
            "records: 18S spans and whole units (18S, ITS1, 5.8S, ITS2, 28S, 5S) of genomes, EUKARYOME and "
            "ROD, labelled by NCBI taxid",
            "checks, chimera screen, confidence",
            "* regions.tsv: V1–V9 of every 18S",
            fmt("* amplicons.tsv: in-silico PCR of {rrna_primer_pairs} primer pairs", n)]),
        ("organelle/", [
            "mitochondrion: SSU, LSU, 5S; COX1, CYTB, COX3*, nad1*, nad5*, atp6*",
            "plastid: SSU, LSU, 5S; rbcL, matK, tufA, psbA*, ycf1*, the trnH–psbA spacer*",
            "* cds, proteins by the organelle's code; checks and confidence per gene",
            fmt("* amplicons.tsv: {organelle_primer_pairs} primer pairs", n),
            "molecules.fna.gz: the organelles, whole"]),
    ]]
    row2 = [Box(t, items, w) for t, items in [
        ("metadata/", [
            "assembly.tsv.gz: size, contiguity, gaps, GC and its spread, masking, CpG, telomeres",
            fmt("{tetranucleotides} tetranucleotides; sourmash sketches", n),
            "ncbi.tsv.gz: NCBI's own summary",
            "* codons.tsv.gz: the genetic-code evidence"]),
        ("* intraspecific/", [
            fmt("{panel_species} species × {panel_assemblies} further assemblies: their copies, rDNA and "
                "organelle genes", n),
            "barcode_gap.tsv: per marker, identity within the species beside its nearest other species"]),
        ("benchmarks/", [
            fmt("{benchmark_tables} tables: completeness, agreement with OrthoDB and the genomes' annotation, "
                "resolution, rDNA, names, primers*, reads and detection, overlap with other databases, "
                "run time", n),
            "summary.json: the headline figures"]),
        ("manifest.json", [
            "format 5; the version of every input; the package's counts; every file's size and MD5; "
            "written last",
            "read by protal (its 30 BUSCOs), phylogenetics, primer design, contamination screens: plain TSV, "
            "by column header"]),
    ]]
    return w, row1, row2


def panel_b(ax, n, marks, counts):
    """The package: one row per species, and what each part holds."""
    w, row1, row2 = panel_b_content(n)
    gap = 2.6
    top_pad = 13
    h1, h2 = max(b.h for b in row1), max(b.h for b in row2)
    H = top_pad + h1 + gap + h2 + 2
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(0, H - 1, "b", fontsize=FS_PANEL, fontweight="bold", va="top", ha="left")
    ax.text(5, H - 1, "The package: plain TSV and gzipped FASTA, one row per species", fontsize=FS_HEAD,
            fontweight="bold", va="top", ha="left")
    if counts:
        line = "  ·  ".join(f"{k}: {v}" for k, v in counts)
    else:
        line = fmt("the first full build: {ncbi_species} species with a genome, {tsa_species} with a transcriptome, "
                   "{eukprot_species} with records; each with its {buscos} BUSCOs, rDNA, organelle genes", n)
    ax.text(5, H - 5.6, line, fontsize=FS_BODY, color=INK2, va="top", ha="left")
    top = H - top_pad
    xs = [1 + i * (w + gap) for i in range(4)]
    for x, b in zip(xs, row1):
        b.draw(ax, x, top, SOURCE_FILL, marks)
    for x, b in zip(xs, row2):
        b.draw(ax, x, top - h1 - gap, SOURCE_FILL, marks)
    return H


def manifest_counts(package):
    with open(os.path.join(package, "manifest.json")) as fh:
        m = json.load(fh)
    c = m["counts"]
    organelle = sum(sum(v.values()) for v in c.get("organelle_genes", {}).values())
    return [("species with a genome", f"{c['species']:,}"), ("with rDNA alone", f"{c.get('species_rrna_only', 0):,}"),
            ("BUSCO copies", f"{c['copies']:,}"), ("rDNA records", f"{c['rrna_records']:,}"),
            ("organelle genes", f"{organelle:,}"), ("taxa", f"{c['taxa']:,}")]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", help="a package whose manifest's counts head panel b")
    p.add_argument("--out", default=HERE)
    p.add_argument("--name", default="overview")
    p.add_argument("--format", default="png", choices=("png", "svg", "pdf"))
    p.add_argument("--no_marks", action="store_true", help="no marks on the parts added in 0.2")
    a = p.parse_args(argv)
    counts = manifest_counts(a.package) if a.package else None
    os.makedirs(a.out, exist_ok=True)
    # the panels' heights come from their content: lay them out once on throwaway axes
    probe = plt.figure()
    ha = panel_a(probe.add_subplot(211), NUMBERS, not a.no_marks)
    hb = panel_b(probe.add_subplot(212), NUMBERS, not a.no_marks, counts)
    plt.close(probe)
    fig = plt.figure(figsize=(WIDTH / 25.4, (ha + hb) / 25.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[ha, hb], hspace=0, left=0, right=1, top=1, bottom=0)
    panel_a(fig.add_subplot(gs[0]), NUMBERS, not a.no_marks)
    panel_b(fig.add_subplot(gs[1]), NUMBERS, not a.no_marks, counts)
    print(pb.save(fig, a.out, a.name, a.format))


if __name__ == "__main__":
    main()
