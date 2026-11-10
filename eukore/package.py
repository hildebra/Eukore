"""Eukore's package (docs/FORMAT.md; format 2 since 2026-11-05: Eukore's own locus calls in NCBI's genomes; format 3 since
2026-11-05: the organelle genes; format 4 since 2026-11-05: every BUSCO of the dataset, in busco/, with
busco/loci.tsv; each species' PR2, EUKARYOME and KSGP lineages, the genomes' rRNA counts, the copies' and records'
contig ends, gaps and confidence), written from the build's genome set and work folder, read and checked by its users
(protal's converter).
"""

import collections
import datetime
import gzip
import json
import os
import shutil
import sys

from . import gencode
from . import qc
from . import store
from . import util
from . import __version__

FORMAT = 6
# the formats a reader of this one takes (5 adds files and columns, and regroups supergroup; 6, since 2026-11-09, adds
# columns and the confidence tier `unverified`, which format 5 counted as medium)
READABLE = (4, 5, 6)
TIERS = ("high", "unverified", "medium", "low")  # confidence, best first (format 6; until then high, medium, low)
LOCUS_SETS = ("busco", "markers")  # the loci's folders: BUSCO's, Eukore's own markers (markers.py; since 2026-11-07)
MANIFEST = "manifest.json"
LOCI_COLUMNS = ("locus", "name", "cutoff", "hmm_length", "species", "called", "duplicated", "fragmented", "call_rate",
                "median_length", "median_cds_length", "median_protein_length")
GENOME_COLUMNS = ("accession", "fetched", "species_taxid", "taxid", "organism", "strain", "level", "genome_size",
                  "contigs", "category", "flags", "ftp", "source", "representative", "loci", "rrna_18s", "rdna_units",
                  "rdna_complete", "rdna_contigs", "local",
                  # since 2026-11-07 (format 5): the table its copies translate by; its duplicated BUSCOs; its
                  # BUSCOs' introns (introns.py); its code checked (codecheck.py)
                  "gencode", "duplicated_share", "dup_identity_median", "duplication", "introns",
                  "introns_per_kb_cds", "median_intron_length", "gencode_check", "gencode_codons", "best_table",
                  "type_material",
                  # since 2026-11-09 (format 6): what made its calls (genomes.GENOME_RECORD); its copies another
                  # clade's and the LCA of those clades (qc.cross_summary)
                  "genome_record", "foreign_copies", "foreign_source")
COPY_COLUMNS = ("accession", "locus", "status", "contig", "contig_length", "start", "end", "strand", "length",
                "cds_length", "exons", "identity", "coverage", "hmm_score", "cutoff", "hmm_coverage", "copies",
                "query", "rescued", "model", "protein_length", "window_class", "fcs", "own_species", "same_as", "length_check",
                "contig_edge", "gaps", "confidence", "flags",
                "internal_stops", "dup_identity",  # since 2026-11-07 (format 5)
                # since 2026-11-09 (format 6): its coding sequence against the protein confirmed (loci.cds_state);
                # the other clade's species it is shared with (qc.verdict); how it was chosen among a Duplicated
                # locus' copies: hmm (the best HMM score), lineage (another, whose lineage fits: collect)
                "cds_check", "cross", "choice")
DUP_LOCI_MIN = 20   # duplication_of: a genome with fewer BUSCOs called says nothing of its duplication
DUP_SHARE = 0.15    # this share of its called BUSCOs Duplicated: a duplicated genome
HAPLOTYPE_ID = 99.5  # its duplicated BUSCOs' copies this alike (median, %): an assembly's two haplotypes
# duplication_of's `haplotigs` (since 2026-11-10; the audit's B3): the alleles of an outbred diploid are under
# HAPLOTYPE_ID alike (1-4% heterozygosity: coding sequences 96-99.5%), as a recent polyploid's copies are; what tells
# them is where the copies lie: an assembly's second haplotype on contigs far shorter than the ones its collapsed
# (Single) genes lie on, a polyploid's on chromosomes like any gene.
HAPLOTIG_ID = 95.0   # alleles are this alike at least (a contaminant's copies, on short contigs too, are not)
SHORT_CONTIG = 0.25  # a contig under this share of the Single copies' median contig is short
SHORT_EXCESS = 0.25  # the Duplicated loci with a copy on a short contig, over the share chance gives: haplotigs
SHORT_NOISE = 3.0    # and that excess this many times its sampling error at least
SHORT_SINGLE_MIN = 10  # Single copies with a contig length, to say what a short contig is
FEATURE_COLUMNS = ("accession", "locus", "kind", "start", "end")
SPECIES_COLUMNS = ("taxid", "name", "representative", "source", "loci", "rrna_records", "organelle_genes",
                   "supergroup", "group", "ncbi_lineage", "ncbi_lineage_filled", "pr2", "pr2_match", "eukaryome",
                   "eukaryome_match", "ksgp", "ksgp_match", "silva", "silva_match", "unite", "unite_match",
                   "midori2", "midori2_match", "type_material")  # group, ncbi_lineage_filled, silva .. midori2_match,
# type_material: since format 5
REFERENCES = ("PR2", "EUKARYOME", "KSGP", "SILVA", "UNITE", "MIDORI2")  # the last three since format 5
GAP_MIN = 10        # a run of Ns this long in a copy: an assembly gap (features.tsv's gap)
EDGE = 1000         # a copy this close to its contig's end: maybe cut, or joined to another sequence
FCS_FOREIGN = {"EXCLUDE", "TRIM", "FIX"}
FOREIGN_WINDOW = {"bacteria", "archaea", "prokarya"}


SHARED_VALUES = ("locus", "status", "strand", "cutoff", "copies", "exons", "query", "window_class", "rescued", "model",
                 "dup_identity", "contig", "contig_length")  # copy rows' values that repeat: kept once (collect)


def gap_runs(seq, min_run=GAP_MIN):
    """[(start, end)] of the runs of min_run Ns or more in a copy (0-based, end excluded)."""
    out, i, seq = [], 0, seq.upper()
    while True:
        i = seq.find("N" * min_run, i)
        if i < 0:
            return out
        j = i
        while j < len(seq) and seq[j] == "N":
            j += 1
        out.append((i, j))
        i = j


def tier(low, medium, unverified):
    """A tier from its reasons (format 6): low where one says another organism's or not the gene, medium where one
    says something may be wrong, unverified where nothing is wrong but a check could not be made, else high -> (tier,
    flags: the reasons in that order, or -)."""
    name = "low" if low else "medium" if medium else "unverified" if unverified else "high"
    return name, ",".join(low + medium + unverified) or "-"


def copy_confidence(r, family=False, duplicated_genome=False):
    """A copy's tier and why (busco/copies.tsv's confidence and flags; tier()): low where a check calls it another
    organism's or not the locus (own_species f, FCS-GX's EXCLUDE/TRIM/FIX, tiara's prokaryote window, a marker copy
    whose best profile is a decoy's: `paralog`, since 2026-11-09); medium where something may be wrong (Fragmented,
    Duplicated, called by the rescue pass, its protein's length off its relatives', an assembly gap in it, at a
    contig's end, an organelle window, FCS-GX's REVIEW, a transcriptome's; since 2026-11-09 a stop inside its protein
    (`internal_stops`), a coding sequence not reading as the protein checked (`cds_mismatch`), a sequence shared with
    another clade's genome, its owner undecided (`cross_shared`)); unverified where nothing is wrong but its
    own-species check could not decide (no relative near enough; until 2026-11-09 medium, as if something were); else
    high. family: a multigene family's marker (markers.py's kind family: actin, tubulins, ...), whose Duplicated is
    its nature, no flag. duplicated_genome (since 2026-11-10; the audit's B3): the copy's genome holds its BUSCOs
    twice genome-wide (duplication_of's haplotypes, haplotigs or wgd_or_polyploid): two copies are that genome's
    truth, a Duplicated copy no flag there either (its status and `copies` still say so)."""
    low, medium, unverified = [], [], []
    if r.get("own_species") == "f":
        low.append("own_species")
    fcs = str(r.get("fcs", "-")).split(":")[0]
    if fcs in FCS_FOREIGN:
        low.append("fcs")
    if r.get("window_class", "-") in FOREIGN_WINDOW:
        low.append("window_prokaryote")
    if r.get("paralog") == "t":  # markers.recheck: its best profile a decoy's (EFL as EF-1a)
        low.append("paralog")
    if r.get("transcript") == "t":  # a transcriptome's copy (TSA): no introns, its Duplicated its isoforms
        medium.append("transcript")
    if r.get("status") == "Fragmented" or (r.get("status") == "Duplicated" and r.get("transcript") != "t"
                                           and not family and not duplicated_genome):
        medium.append(r["status"].lower())
    if r.get("rescued") == "t":  # loci.py's rescue pass: its span may fall short of the gene's ends
        medium.append("rescued")
    if r.get("length_check", "ok") != "ok":
        medium.append(f"length_{r['length_check']}")
    if int(r.get("gaps") or 0):
        medium.append("gap")
    if r.get("contig_edge") not in (None, "-") and int(r["contig_edge"]) < EDGE:
        medium.append("contig_edge")
    if r.get("window_class") == "organelle":
        medium.append("window_organelle")
    if fcs == "REVIEW":
        medium.append("fcs_review")
    if int(r.get("internal_stops") or 0) or r.get("cds_check") == "stops":
        medium.append("internal_stops")
    if r.get("cds_check") == "mismatch":
        medium.append("cds_mismatch")
    if r.get("own_species") in ("-", None):
        if r.get("cross") not in (None, "-", ""):
            medium.append("cross_shared")
        else:
            unverified.append("own_species_undecided")
    return tier(low, medium, unverified)


def type_material_of(row):
    """An assembly's relation to its species' type material (NCBI's relation_to_type_material: `assembly from type
    material`, `assembly from synonym type material`, `assembly from pathotype material`, `assembly designated as
    neotype`, ...; since 2026-11-07), `-` for none or an assembly NCBI's tables do not hold."""
    value = (row or {}).get("type_material") or "na"
    return "-" if value == "na" else value


COPY_NUMBER_COLUMNS = ("accession", "species_taxid", "set", "element", "copies", "partial", "basis")
COPY_NUMBERS = "copy_numbers.tsv.gz"


def write_gz(path, columns, rows):
    """A gzipped TSV table (- for an empty value; BGZF since 2026-11-09)."""
    with util.bgzf_open(util.part_of(path)) as fh:
        fh.write("\t".join(columns) + "\n")
        fh.writelines("\t".join("-" if r.get(c) in (None, "") else str(r[c]) for c in columns) + "\n" for r in rows)
    os.replace(util.part_of(path), path)


def read_gz(path):
    """A gzipped TSV table -> [row dicts]."""
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def organelle_molecules_of(folder):
    """A done genome's organelle molecules in its assembly (metadata.describe's stats.json) -> {mitochondrion: n,
    plastid: n}, {} where not measured."""
    from .metadata import STATS  # noqa: PLC0415
    try:
        stats = json.loads(store.text(folder, STATS) or "")
    except ValueError:
        return {}
    return {k: int(stats[f"{k}_molecules"]) for k in ("mitochondrion", "plastid")
            if str(stats.get(f"{k}_molecules", "")).isdigit()}


def copy_numbers(genomes, copy_rows, loci, set_of, organelle=()):
    """copy_numbers.tsv.gz's rows (since 2026-11-08): per species' genome and element, its copies and partial copies
    in the assembly. set busco or markers: each locus searched (complete copies; a Fragmented one 0 and 1 partial; a
    Missing one 0 and 0); rrna: rrna.element_counts's; mitochondrion and plastid: the molecules its FASTA titles name an
    organelle, and the organelle genes found on them (organelle/genes.tsv's of source assembly: their copies on the
    molecules, a plastid's inverted repeat 2). basis: assembly, or transcripts for a transcriptome (TSA: its copies are
    isoforms, not genes). tools/rdna_depth.py adds rows of basis reads. genomes: collect's (their _searched,
    _rrna_copies, _organelle_molecules)."""
    by = collections.defaultdict(dict)
    for r in copy_rows:
        by[r["accession"]][r["locus"]] = r
    genes = collections.defaultdict(list)
    for _n, _s, row in organelle or ():
        if row.get("source") == "assembly" and row.get("accession") not in (None, "-"):
            genes[row["accession"]].append(row)
    out = []
    for acc, g in sorted(genomes.items()):
        if g.get("representative") != "t" or g.get("source") == "no_assembly":  # EukDetect2's records: no genome
            continue
        base = dict(accession=acc, species_taxid=g["species_taxid"],
                    basis="transcripts" if g.get("source") == "tsa" else "assembly")
        for item in loci:
            locus = item["locus"]
            if locus not in g.get("_searched", ()):  # a genome run before the locus was (a marker): not 0, not known
                continue
            r = by[acc].get(locus)
            status = r["status"] if r else "Missing"
            out.append(dict(base, set=set_of.get(locus, "busco"), element=locus,
                            copies=int(r["copies"]) if status in ("Single", "Duplicated") else 0,
                            partial=int(status == "Fragmented")))
        for element, (n, partial) in (g.get("_rrna_copies") or {}).items():
            out.append(dict(base, set="rrna", element=element, copies=n, partial=partial))
        for kind, n in (g.get("_organelle_molecules") or {}).items():
            out.append(dict(base, set=kind, element="molecules", copies=n, partial=0))
        for row in sorted(genes[acc], key=lambda r: (r["organelle"], r["gene"])):
            out.append(dict(base, set=row["organelle"], element=row["gene"], copies=int(row.get("copies") or 1),
                            partial=int(row.get("partial") in (True, "t"))))
    return out


def duplication_of(rows, alternates=None):
    """A genome's duplicated BUSCOs (its copies.tsv rows; since 2026-11-07) -> dict(duplicated_share (of those called
    Single or Duplicated), dup_identity_median (the Duplicated ones' dup_identity, "<50" read as 50), duplication,
    dup_short_contigs): `haplotypes` where many are duplicated and their copies are HAPLOTYPE_ID% alike or more (an
    assembly holding both haplotypes of a diploid, or a recent hybrid), `wgd_or_polyploid` where many are duplicated
    and less alike (a whole genome duplication, a polyploid, an old hybrid), `low` where few are, `-` where too few
    are called. A heuristic: a reader weighs it with the genome's size and the species' biology.

    Since 2026-11-10 `haplotigs`, a class of its own: many duplicated, HAPLOTIG_ID to HAPLOTYPE_ID alike, and their
    copies on short contigs (under SHORT_CONTIG of the Single copies' median contig) SHORT_EXCESS more often than
    chance gives (the share of the Single copies on such contigs; with two copies seen, either's): an assembly not
    purged of an outbred diploid's second haplotype. Until then `wgd_or_polyploid`, which still is every genome whose
    copies are less alike than HAPLOTYPE_ID and not told so: a polyploid, and two haplotypes assembled equally long.
    alternates: the genome's other complete copies' rows (its alternates.tsv: locus, contig_length), where the caller
    has them; without them the chosen copies' contigs alone are seen, half of them the first haplotype's.
    dup_short_contigs: `<share of the Duplicated loci on a short contig>/<the share chance gives>`, - where not
    told (no column of genomes.tsv yet)."""
    called = [r for r in rows if r.get("status") in ("Single", "Duplicated")]
    ids = sorted(float(str(r["dup_identity"]).lstrip("<")) for r in called
                 if r["status"] == "Duplicated" and r.get("dup_identity") not in (None, "", "-"))
    if len(called) < DUP_LOCI_MIN:
        return dict(duplicated_share="-", dup_identity_median="-", duplication="-", dup_short_contigs="-")
    share = sum(r["status"] == "Duplicated" for r in called) / len(called)
    median = ids[len(ids) // 2] if ids else None

    def length(r):
        return int(r["contig_length"]) if str(r.get("contig_length", "")).isdigit() else None
    single = sorted(n for n in (length(r) for r in called if r["status"] == "Single") if n)
    others = collections.defaultdict(list)
    for a in alternates or ():
        if length(a):
            others[a["locus"]].append(length(a))
    seen = [[length(r)] + others.get(r.get("locus"), []) for r in called if r["status"] == "Duplicated" and length(r)]
    short, chance, noise = None, None, 0.0
    if len(single) >= SHORT_SINGLE_MIN and seen:
        limit = SHORT_CONTIG * single[len(single) // 2]
        by_chance = sum(n < limit for n in single) / len(single)
        short = sum(min(ns) < limit for ns in seen) / len(seen)
        chance = sum(1 - (1 - by_chance) ** len(ns) for ns in seen) / len(seen)
        # both shares' sampling error (a genome nearly all Duplicated has few Single copies to say what chance is)
        slope = sum(len(ns) * (1 - by_chance) ** (len(ns) - 1) for ns in seen) / len(seen)
        noise = (short * (1 - short) / len(seen) + slope ** 2 * by_chance * (1 - by_chance) / len(single)) ** 0.5
    kind = ("low" if share < DUP_SHARE else "-" if median is None else
            "haplotypes" if median >= HAPLOTYPE_ID else
            "haplotigs" if median >= HAPLOTIG_ID and short is not None
            and short - chance >= max(SHORT_EXCESS, SHORT_NOISE * noise) else "wgd_or_polyploid")
    return dict(duplicated_share=f"{share:.3f}", dup_identity_median="-" if median is None else f"{median:.1f}",
                duplication=kind, dup_short_contigs="-" if short is None else f"{short:.2f}/{chance:.2f}")


def _tsv(path, columns, rows):
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "-")) if r.get(c, "-") not in (None, "") else "-" for c in columns) + "\n")
    os.replace(util.part_of(path), path)


def done_table(work, accession, name):
    """A done genome's table (WORK/done/<accession>/<name>, or its packed genome's, store.py; header first) -> [row
    dicts]."""
    from .genomes import done_dir  # noqa: PLC0415
    return store.table(done_dir(work, accession), name)


def done_copies(work, accession, swaps=None):
    """A genome's copies (WORK/done/<accession>/copies.fna) -> {locus: copy}; swaps: {locus: alt}, the loci whose
    copy is one of its alternates (collect's choice by lineage, since 2026-11-09)."""
    from .genomes import done_dir  # noqa: PLC0415
    lines = [line.strip() for line in (store.text(done_dir(work, accession), "copies.fna") or "").splitlines()]
    out = {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 2)}
    if swaps:
        alts = {(a["locus"], a["alt"]): a for a in done_alternates(work, accession)}
        out.update({locus: alts[(locus, k)]["copy"] for locus, k in swaps.items() if (locus, k) in alts})
    return out


def done_alternates(work, accession):
    """A genome's alternates (genomes.ALTERNATES: a Duplicated locus' other complete copies whole; since 2026-11-09)
    -> [loci.tsv-like row dicts with alt (int), copy and features [(kind, start, end)]], [] for none."""
    from .genomes import ALTERNATES, done_dir  # noqa: PLC0415
    folder = done_dir(work, accession)
    rows = store.table(folder, ALTERNATES[0])
    if not rows:
        return []
    lines = [line.strip() for line in (store.text(folder, ALTERNATES[1]) or "").splitlines()]
    copies = {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 2)}
    feats = collections.defaultdict(list)
    for f in store.table(folder, ALTERNATES[2]):
        feats[(f["locus"], int(f["alt"]))].append((f["kind"], int(f["start"]), int(f["end"])))
    for r in rows:
        r["alt"] = int(r["alt"])
        r["copy"], r["features"] = copies.get(f"{r['locus']}|{r['alt']}", ""), feats[(r["locus"], r["alt"])]
    return rows


LINEAGE_RANK = {"t": 3, "-": 2, "f": 0}  # own_species, best first; an undecided copy shared with another clade: 1


def lineage_rank(check):
    """A copy's own-species verdict as a rank (collect's choice among a Duplicated locus' copies): t, then undecided,
    then undecided but shared with another clade's genome, then another clade's."""
    return 1 if check["own_species"] == "-" and check.get("cross") else LINEAGE_RANK[check["own_species"]]


def genome_part(args):
    """One done genome's share of collect (since 2026-11-10 a function of its own, so that a pool of processes reads
    the genomes, `threads` at once; until then one after another, ~30,000 of them from the shared file system, and
    their coding sequences, introns and proteins made in the one process): args (work, accession, its genetic code)
    -> None where the genome is not done, else dict(loci (its copies' count), rrna ((counts, element_counts) or None),
    molecules, searched, record, copies [(copy row, coding sequence, (protein, internal stops), intron rows)] in its
    loci.tsv's order, alternates [(locus, alt, coding sequence, protein length)]). Its copies whole and its features
    are not returned: write reads them again (the memory bound of 2026-11-08)."""
    from . import introns as introns_mod, rrna as rrna_mod  # noqa: PLC0415
    from .genomes import done_dir, done_summary, is_done, record_of  # noqa: PLC0415
    work, acc, table = args
    folder = done_dir(work, acc)
    got = done_copies(work, acc)
    # a done genome without a BUSCO called (a divergent protist) is still the species' genome: its rRNA and
    # organelle records name it (until 2026-11-06 it was left out of genomes.tsv, their accession dangling)
    if not got and not is_done(work, acc):
        return None
    rrna = None
    if store.has(folder, "rrna.tsv"):
        rrna_rows = rrna_mod.load(folder)
        rrna = (rrna_mod.counts(rrna_rows), rrna_mod.element_counts(rrna_rows))
    own = collections.defaultdict(list)  # this genome's features by locus
    for f in done_table(work, acc, "features.tsv"):
        own[f["locus"]].append((f["kind"], int(f["start"]), int(f["end"])))
    fcs = [dict(f, start=int(f["start"]), end=int(f["end"])) for f in done_table(work, acc, "fcs.tsv")]
    screened = store.has(folder, "fcs.tsv")
    searched = done_table(work, acc, "loci.tsv")
    copies = []
    for r in searched:
        if r["status"] != "Missing":
            row = copy_row(r, acc, fcs if screened else None, len(gap_runs(got.get(r["locus"], ""))))
            cds = qc.spliced(got.get(r["locus"], ""), own[r["locus"]])
            # its introns (introns.py), their HMM places set by genetics
            introns = introns_mod.of_copy(got[r["locus"]], own[r["locus"]]) if got.get(r["locus"]) else []
            copies.append((row, cds, gencode.protein_of(cds, table), introns))
    # a Duplicated locus' other complete copies (since 2026-11-09, genome record 2): their coding sequences checked
    # beside the copies, for the choice by lineage; their copies read again where one is taken
    alternates = [(a["locus"], a["alt"], qc.spliced(a["copy"], a["features"]), int(a["protein_length"]))
                  for a in done_alternates(work, acc)]
    return dict(loci=len(got), rrna=rrna, molecules=organelle_molecules_of(folder),
                searched={r["locus"] for r in searched}, record=record_of(done_summary(work, acc)) or "-",
                copies=copies, alternates=alternates)


def collect(ncbi, job_map, assemblies, work, vsearch=None, threads=1, kinds=None, paralogs_of=None):
    """What the package holds of the BUSCO genes -> dict(genomes {accession: row}, copy_rows [COPY_COLUMNS dicts],
    taxids, cds and proteins {(accession, locus): sequence}, intron_rows (introns.of_copy's, accession and locus set),
    work and done (the genomes whose copies and features write reads from WORK/done), copy_count, swaps {accession:
    {locus: alt}}: the copies taken by lineage among a Duplicated locus' alternates, since 2026-11-09). Each species'
    genome (its job) with the copies its run called; the species' other assemblies listed (source "assembly", no
    copies); the copies' checks across species (qc.check). kinds: {marker locus: clean or family} (markers.kinds): a
    genome's duplication is read on its BUSCOs alone, a family marker's Duplicated is no flag. paralogs_of: a function
    of the marker copies' coding sequences {(accession, locus): cds} and the genomes -> the keys whose best profile is
    a decoy's (markers.recheck; since 2026-11-09: the decoys a marker set gained since a genome's run, EFL's).
    threads: the genomes read `threads` at once by a pool of processes (genome_part; since 2026-11-10), the result
    the same as one after another; vsearch's threads too."""
    from . import introns as introns_mod  # noqa: PLC0415
    kinds = kinds or {}
    by_accession = {r["accession"]: r for r in assemblies}
    genomes, copy_rows, cds, intron_rows, done, count = {}, [], {}, [], [], 0
    alt_cds, alt_length = {}, {}  # a Duplicated locus' other copies' (since 2026-11-09)
    made = {}  # {(accession, locus): (protein, internal stops)} of the copies as read (genome_part)
    jobs = sorted(job_map.items())
    parts = util.ordered(genome_part, [(work, acc, job.get("gencode", 1)) for acc, job in jobs],
                         max(1, min(threads, len(jobs))))
    for (acc, job), part in zip(jobs, parts):
        if part is None:
            continue
        got_n = part["loci"]
        row = by_accession.get(acc, {})
        genomes[acc] = dict(accession=acc, fetched=job.get("fetch", acc), species_taxid=job["species_taxid"],
                            taxid=job.get("taxid") or job["species_taxid"],
                            organism=row.get("organism") or job.get("organism", ""),
                            strain=row.get("strain") or row.get("isolate") or "-", level=row.get("level") or job.get(
                                "level", ""), genome_size=row.get("size") or job.get("size"),
                            contigs=row.get("contigs"), category=row.get("category") or job.get("category", ""),
                            flags=",".join(job.get("flags") or ()) or "-", ftp=row.get("ftp"), source=job["why"][0],
                            representative="t", loci=got_n, gencode=job.get("gencode", 1),
                            type_material=type_material_of(row))
        if part["rrna"] is not None:
            n, elements = part["rrna"]
            genomes[acc].update(rrna_18s=n["18S"], rdna_units=n["rDNA"], rdna_complete=n["rDNA_complete"],
                                rdna_contigs=n["rDNA_contigs"])
            genomes[acc]["_rrna_copies"] = elements  # copy_numbers.tsv's (no column)
        genomes[acc]["_organelle_molecules"] = part["molecules"]
        # its copies and features read here and again by write, never all held (since 2026-11-08: every copy, introns
        # and all, and every feature as a dict were, ~1 KB a feature: tens of GB for all eukaryotes)
        done.append(acc)
        count += got_n
        name = sys.intern(acc)
        genomes[acc]["_searched"] = part["searched"]  # Missing ones too: copy_numbers.tsv's zeros
        genomes[acc]["genome_record"] = part["record"]  # since 2026-11-09
        for row, seq, protein, introns in part["copies"]:
            row["accession"] = name
            for k in SHARED_VALUES:  # a few values, many rows: one string each (in this process: a worker's are not)
                if isinstance(row.get(k), str):
                    row[k] = sys.intern(row[k])
            copy_rows.append(row)
            cds[(name, row["locus"])], made[(name, row["locus"])] = seq, protein
            intron_rows.extend(dict(i, accession=name, locus=row["locus"]) for i in introns)
        for locus, k, seq, length in part["alternates"]:
            alt_cds[(name, locus, k)], alt_length[(name, locus, k)] = seq, length
        for other in job.get("others", ()):
            o = by_accession.get(other)
            if o and other not in genomes:
                genomes[other] = dict(accession=other, fetched=other, species_taxid=job["species_taxid"],
                                      taxid=o["taxid"], organism=o["organism"],
                                      strain=o.get("strain") or o.get("isolate") or "-", level=o["level"],
                                      genome_size=o["size"], contigs=o["contigs"], category=o["category"],
                                      flags=",".join(o.get("flags") or ()) or "-", ftp=o["ftp"], source="assembly",
                                      type_material=type_material_of(o),
                                      representative="f", loci=0)
    lengths = {(r["accession"], r["locus"]): int(r["protein_length"]) for r in copy_rows}
    checks = qc.check(cds, {a: g["species_taxid"] for a, g in genomes.items()}, ncbi, vsearch, threads,
                      {**lengths, **alt_length}, alt_cds)
    # the choice by lineage (since 2026-11-09; until then the best HMM score's copy, a contaminant's where it scored
    # best): a Duplicated locus whose copy's own-species verdict ranks below one of its other copies' (lineage_rank:
    # another clade's, or undecided where the other agrees with the lineage) takes that one, its row, coding sequence
    # and introns, and write its copy and features (`choice` lineage)
    swaps, best_alt = collections.defaultdict(dict), {}
    for (acc, locus, k) in sorted(alt_cds):
        c = checks[(acc, locus, k)]
        if (acc, locus) in checks and (lineage_rank(c), -k) > best_alt.get((acc, locus), (-1, 0))[:2]:
            best_alt[(acc, locus)] = (lineage_rank(c), -k, c)
    taken = {key: (-negk, c) for key, (rank, negk, c) in best_alt.items() if rank > lineage_rank(checks[key])}
    if taken:
        from .genomes import done_dir  # noqa: PLC0415
        by_key = {(r["accession"], r["locus"]): i for i, r in enumerate(copy_rows)}
        new_introns = []
        for acc in sorted({a for a, _l in taken}):
            alts = {(a["locus"], a["alt"]): a for a in done_alternates(work, acc)}
            screened = store.has(done_dir(work, acc), "fcs.tsv")
            fcs = [dict(f, start=int(f["start"]), end=int(f["end"]))
                   for f in done_table(work, acc, "fcs.tsv")] if screened else None
            for locus, (k, c) in sorted((lo, v) for (a_, lo), v in taken.items() if a_ == acc):
                if (locus, k) not in alts:
                    continue
                alt, old = alts[(locus, k)], copy_rows[by_key[(acc, locus)]]
                row = copy_row(dict(alt, dup_identity=old.get("dup_identity", "-")), old["accession"], fcs,
                               len(gap_runs(alt["copy"])))
                row["choice"] = "lineage"
                copy_rows[by_key[(acc, locus)]] = row
                cds[(acc, locus)], checks[(acc, locus)] = alt_cds[(acc, locus, k)], c
                swaps[acc][locus] = k
                new_introns.extend(dict(i, accession=old["accession"], locus=locus)
                                   for i in introns_mod.of_copy(alt["copy"], alt["features"]))
        swapped = {(a, lo) for a, s in swaps.items() for lo in s}
        intron_rows = [i for i in intron_rows if (i["accession"], i["locus"]) not in swapped] + new_introns
    paralogs = paralogs_of({k: v for k, v in cds.items() if k[1] in kinds}, genomes) if paralogs_of and kinds else set()
    proteins = {}
    for r in copy_rows:
        key = (r["accession"], r["locus"])
        c = checks[key]
        r.update(own_species=c["own_species"], same_as=c["same_as"], length_check=c["length"],
                 cross=c["cross"] if c.get("cross") else "-",
                 nearest_id=c["nearest_id"], nearest_depth=c["nearest_depth"])  # the last two: benchmarks', no columns
        if genomes[r["accession"]].get("source") == "tsa":  # since 2026-11-07: a transcriptome's copy
            r["transcript"] = "t"
        if key in paralogs:  # a marker copy whose best profile is a decoy's (markers.recheck; since 2026-11-09)
            r["paralog"] = "t"
        # its protein before its tier (since 2026-11-09: a stop inside it is a flag); made by genome_part, again for a
        # copy taken by lineage
        swapped = r["locus"] in swaps.get(r["accession"], ())
        proteins[key], r["internal_stops"] = made[key] if not swapped else gencode.protein_of(
            cds[key], genomes[r["accession"]].get("gencode", 1))
        r["confidence"], r["flags"] = copy_confidence(r, kinds.get(r["locus"]) == "family")
    by_genome, all_rows = collections.defaultdict(list), collections.defaultdict(list)
    for r in copy_rows:
        all_rows[r["accession"]].append(r)
        if r["locus"] not in kinds:  # the BUSCOs' (single-copy by their choice): the markers' are not
            by_genome[r["accession"]].append(r)
    for acc, g in genomes.items():
        if g["representative"] == "t" and g.get("source") != "tsa":  # a transcriptome's isoforms are no duplication
            g.update(duplication_of(by_genome.get(acc, [])))
        if g["representative"] == "t":  # its copies another clade's (since 2026-11-09: cross-contamination)
            g["foreign_copies"], g["foreign_source"] = qc.cross_summary(all_rows.get(acc, []), ncbi)
    taxids = {g["taxid"] for g in genomes.values()} | {int(r["same_as"]) for r in copy_rows if r["same_as"] != "-"} \
        | {int(g["foreign_source"]) for g in genomes.values() if str(g.get("foreign_source", "-")).isdigit()}
    return dict(genomes=genomes, copy_rows=copy_rows, taxids=taxids, cds=cds, proteins=proteins,
                intron_rows=intron_rows, work=work, done=done, copy_count=count, swaps=dict(swaps))


def copy_row(r, accession, fcs, gaps):
    """A copy's row of busco/copies.tsv from its done genome's loci.tsv row (or an alternate's): its FCS-GX call (fcs:
    the genome's ranges, None where it was not screened), its distance to its contig's end, its assembly gaps."""
    from .genomes import fcs_of_copy  # noqa: PLC0415
    edge = min(int(r["start"]) - 1, int(r["contig_length"]) - int(r["end"])) \
        if str(r.get("contig_length", "")).isdigit() and str(r.get("start", "")).isdigit() else "-"
    out = {k: v for k, v in r.items() if k not in ("alt", "copy", "features")}
    return dict(out, accession=accession, fcs=fcs_of_copy(fcs, r["locus"], r["contig"], r["start"], r["end"])
                if fcs is not None else "-", contig_edge=edge, gaps=gaps, cds_check=r.get("cds_check") or "-",
                choice="hmm")


def genome_features(work, accession, copies=None, swaps=None):
    """A done genome's features (features.tsv, and its copies' assembly gaps: gap_runs) as features.tsv rows
    (FEATURE_COLUMNS), in the package's order (locus, start, kind). copies: done_copies' (read where not given);
    swaps: {locus: alt} whose copy is an alternate's (done_copies), its features the alternate's."""
    copies = done_copies(work, accession, swaps) if copies is None else copies
    rows = [(f["locus"], int(f["start"]), f["kind"], int(f["end"])) for f in done_table(work, accession, "features.tsv")
            if not swaps or f["locus"] not in swaps]
    if swaps:
        rows += [(a["locus"], s, kind, e) for a in done_alternates(work, accession)
                 if swaps.get(a["locus"]) == a["alt"] for kind, s, e in a["features"]]
    rows += [(locus, a, "gap", b) for locus, seq in copies.items() for a, b in gap_runs(seq)]
    return [dict(accession=accession, locus=lo, kind=k, start=s, end=e) for lo, s, k, e in sorted(rows)]


def genome_files(args):
    """A done genome's lines of write's files (since 2026-11-10 a function of its own: a pool of processes reads the
    genomes): args (work, accession, swaps: done_copies') -> ([(locus, its copy's FASTA record)], [(locus, its
    features.tsv lines joined)]), in write's order: genome_features' rows are by locus, a locus' lines one string (a
    few strings a genome to pass between the processes, not one a feature)."""
    work, acc, swaps = args
    got = done_copies(work, acc, swaps)
    feats = {}
    for f in genome_features(work, acc, got, swaps):
        feats.setdefault(f["locus"], []).append("\t".join(str(f[c]) for c in FEATURE_COLUMNS) + "\n")
    return ([(locus, f">{acc}\n{seq.upper()}\n") for locus, seq in got.items() if seq],
            [(locus, "".join(lines)) for locus, lines in feats.items()])


def genetics(collected, ncbi, hmm_of, hmmalign="hmmalign", work=None, log=print, threads=1):
    """The copies' proteins placed on their HMMs (align.py, `threads` loci at once), their introns (introns.py) and
    each genome's genetic code checked (codecheck.py; since 2026-11-07) -> collected, with `introns`, `intron_sites`,
    `intron_site_columns` and `codons` and the genomes' rows completed; without hmmalign, introns without their places
    and no code check."""
    from . import align as align_mod, codecheck, groups, introns as introns_mod  # noqa: PLC0415
    states = align_mod.align(collected, hmm_of, hmmalign, work, log, threads) if hmmalign else {}
    group_of = {a: groups.group_of(ncbi, int(g["species_taxid"])) for a, g in collected["genomes"].items()
                if g["representative"] == "t"}
    rows, site_rows, site_columns, per_genome = introns_mod.annotate(collected, states, group_of)
    codons, checked = codecheck.run(collected, states, hmm_of) if states else ([], {})
    for acc, g in collected["genomes"].items():
        g.update(per_genome.get(acc, {}))
        g.update(checked.get(acc, {}))
    collected.update(introns=rows, intron_sites=site_rows, intron_site_columns=site_columns, codons=codons)
    return collected


def loci_rows(loci, collected, ncbi):
    """busco/loci.tsv's rows: each BUSCO (loci: dicts of locus, name, cutoff, hmm_length) with how the package's
    species hold it: species with a genome searched for it, called (Single or Duplicated), duplicated, fragmented, the
    call rate over all and per supergroup (groups.GROUPS), the median copy, coding and protein lengths. A genome run
    before a locus was (a marker added later: collect's _searched) counts for it nowhere (until 2026-11-08 as one not
    calling it: its call rates too low).

    Since 2026-11-10 (the audit's B5) also the call rates over the genomes alone: `species`, `called`, `call_rate`
    and `call_<group>` count every species' record, a transcriptome (source tsa) and EukDetect2's records
    (no_assembly) with the assemblies, and a BUSCO a transcriptome misses is a gene not expressed (or not assembled),
    not one the species lost. Further columns, the old ones as they were: `genomes` and `genome_call_rate` (the
    species whose record is an assembly, searched for the locus; of them called), `genome_call_<group>` (per group
    with such a species), `transcriptomes` and `transcriptome_call_rate`; `-` where there is none to rate. None
    starts with `call_`: a reader takes those for the groups."""
    from . import groups  # noqa: PLC0415
    genomes = collected["genomes"]
    reps = {a: g for a, g in genomes.items() if g.get("representative") == "t"}
    group = {a: groups.group_of(ncbi, g["species_taxid"]) for a, g in reps.items()}
    names = sorted(set(group.values()))
    kind = {a: {"tsa": "transcriptome", "no_assembly": "records"}.get(g.get("source"), "genome")
            for a, g in reps.items()}
    genome_names = sorted({group[a] for a in reps if kind[a] == "genome"})

    def rate(n, of):
        return f"{n / of:.3f}" if of else "-"
    by_locus = collections.defaultdict(list)
    for r in collected["copy_rows"]:
        if r["accession"] in reps:
            by_locus[r["locus"]].append(r)
    out = []
    for item in loci:
        searched = [a for a, g in reps.items() if g.get("_searched") is None or item["locus"] in g["_searched"]]
        sizes = collections.Counter(group[a] for a in searched)
        rs = by_locus.get(item["locus"], [])
        called = [r for r in rs if r["status"] in ("Single", "Duplicated")]
        per_group = collections.Counter(group[r["accession"]] for r in called)

        def median(key):
            values = sorted(int(r[key]) for r in called if str(r.get(key, "")).isdigit())
            return values[len(values) // 2] if values else 0
        row = dict(item, species=len(searched), called=len(called),
                   duplicated=sum(r["status"] == "Duplicated" for r in rs),
                   fragmented=sum(r["status"] == "Fragmented" for r in rs),
                   call_rate=f"{len(called) / max(len(searched), 1):.3f}", median_length=median("length"),
                   median_cds_length=median("cds_length"), median_protein_length=median("protein_length"))
        for name in names:
            row[f"call_{name}"] = f"{per_group[name] / max(sizes[name], 1):.3f}"
        of = collections.Counter(kind[a] for a in searched)
        got = collections.Counter(kind[r["accession"]] for r in called)
        row.update(genomes=of["genome"], genome_call_rate=rate(got["genome"], of["genome"]),
                   transcriptomes=of["transcriptome"],
                   transcriptome_call_rate=rate(got["transcriptome"], of["transcriptome"]))
        of = collections.Counter(group[a] for a in searched if kind[a] == "genome")
        got = collections.Counter(group[r["accession"]] for r in called if kind[r["accession"]] == "genome")
        for name in genome_names:
            row[f"genome_call_{name}"] = rate(got[name], of[name])
        out.append(row)
    return out, [f"call_{name}" for name in names] + ["genomes", "genome_call_rate", "transcriptomes",
                                                      "transcriptome_call_rate"] + [f"genome_call_{name}"
                                                                                    for name in genome_names]


def unverified_alone(copy_rows, genomes, ncbi, set_of=None):
    """The copies of tier `unverified` and, of them, those of a lone lineage (since 2026-11-10; the audit's B4) ->
    {locus set: dict(copies, alone_in_family)}. A copy's own-species check passes only where a species of its family
    (or genus) holds the locus as its nearest (qc.verdict); a species alone in its family among the package's
    holders of the locus (its genus, where NCBI names it no family; neither: alone) has no one to agree with, and
    its copies are `unverified` by the package's sampling, not by anything in them. The manifest's
    `counts.copy_unverified`; release.stats' `unverified_alone`. genomes: {accession: row with species_taxid};
    set_of: {locus: busco or markers}."""
    cache, holders, unverified = {}, collections.defaultdict(set), []

    def clade(sp):
        if sp not in cache:
            ranks = ncbi.rank_ids(sp)
            cache[sp] = ranks.get("family") or ranks.get("genus")
        return cache[sp]
    for r in copy_rows:
        sp = int(genomes[r["accession"]]["species_taxid"])
        holders[(r["locus"], clade(sp))].add(sp)
        if r.get("confidence") == "unverified":
            unverified.append((r["locus"], sp))
    out = {}
    for locus, sp in unverified:
        n = out.setdefault((set_of or {}).get(locus, "busco"), dict(copies=0, alone_in_family=0))
        n["copies"] += 1
        n["alone_in_family"] += clade(sp) is None or len(holders[(locus, clade(sp))]) < 2
    return out


def write(out, version, ncbi, collected, loci, rrna=None, genomes_folder=None, inputs=None, organelle=None,
          references=(), before_manifest=None, metadata=None, tables=(), writers=(), threads=1):
    """The package into `out` (written beside it and renamed into place; the manifest last) -> its folder.
    before_manifest: a function of the unfinished package's folder, run once its tables are written (the build's
    benchmarks: benchmarks.run), whose result (a dict) the manifest keeps as `benchmarks`.
    loci: the BUSCOs called ([{locus, name, cutoff, hmm_length}], in the dataset's order); rrna: (records [(name,
    sequence, row dict of records.tsv)]) or None; organelle: the organelle genes' records (organelle.collect) or
    None; references: labels.Reference of PR2, EUKARYOME and KSGP, each species' lineage in them (species.tsv);
    metadata: a function of (the unfinished folder, the genomes' rows) writing the assemblies' metadata
    (metadata.write: metadata/, organelle/molecules.fna.gz), its counts the manifest's `metadata`; tables: further
    tables [(path in the package, columns, rows)] (since 2026-11-07: rrna/regions.tsv, the amplicons); writers:
    functions of the unfinished folder writing further files (since 2026-11-07: intraspecific/); threads: the done
    genomes' copies and features read `threads` at once by a pool of processes (genome_files; since 2026-11-10), the
    files the same as one after another."""
    part = util.part_of(out)
    shutil.rmtree(part, ignore_errors=True)
    # each locus' set: busco/ (BUSCO's), markers/ (Eukore's own, markers.py; since 2026-11-07), the same files in each
    set_of = {item["locus"]: item.get("set", "busco") for item in loci}
    sets = [s for s in LOCUS_SETS if s in set_of.values()] or ["busco"]
    for sub in ["taxonomy", "rrna", "organelle"] + [f"{s}/{k}" for s in sets for k in ("copies", "cds", "proteins")]:
        os.makedirs(os.path.join(part, sub))
    genomes = collected["genomes"]
    local = set()
    if genomes_folder and os.path.isdir(genomes_folder):
        local = {f[:-len(".fna.gz")] for f in os.listdir(genomes_folder) if f.endswith(".fna.gz")}
    for acc, r in genomes.items():
        r["local"] = f"genomes/{acc}.fna.gz" if acc in local else "-"
    _tsv(os.path.join(part, "genomes.tsv"), GENOME_COLUMNS,
         sorted(genomes.values(), key=lambda r: (str(r["species_taxid"]), r["representative"] != "t", r["accession"])))
    for s in sets:
        rows_of_loci, call_columns = loci_rows([i for i in loci if set_of[i["locus"]] == s], collected, ncbi)
        _tsv(os.path.join(part, s, "loci.tsv"), LOCI_COLUMNS + tuple(call_columns), rows_of_loci)
    order = {item["locus"]: i for i, item in enumerate(loci)}
    by_locus = {"cds": collections.defaultdict(dict), "proteins": collections.defaultdict(dict)}
    for sub in by_locus:  # once, not a pass over every copy per locus (until 2026-11-08)
        for (a, lo), s in collected.get(sub, {}).items():
            by_locus[sub][lo][a] = s
    for item in loci:
        locus = item["locus"]
        # each copy's coding sequence (exons joined, frameshifts' bases out); its protein (since 2026-11-07)
        for sub, ext, seqs in (("cds", "fna", by_locus["cds"].get(locus, {})),
                               ("proteins", "faa", by_locus["proteins"].get(locus, {}))):
            with util.bgzf_open(os.path.join(part, set_of[locus], sub, f"{locus}.{ext}.gz")) as fh:
                fh.writelines(f">{acc}\n{seqs[acc].upper()}\n" for acc in sorted(seqs) if seqs[acc])
    # each copy whole and the features: read again from each done genome in turn into every locus' file at once
    # (since 2026-11-08: all held in memory from collect on); a Duplicated locus' copy taken by lineage (collect's
    # swaps); the genomes read `threads` at once (since 2026-11-10: genome_files), written here in their order, by
    # workers forked before the files' compression threads start (since 2026-11-10: after)
    done = sorted(collected.get("done", ()))
    items = [(collected["work"], acc, collected.get("swaps", {}).get(acc)) for acc in done]
    with util.forked_pool(max(1, min(threads, len(items)))) as pool:
        copy_files = {item["locus"]: util.bgzf_open(os.path.join(part, set_of[item["locus"]], "copies",
                                                                 f"{item['locus']}.fna.gz"))
                      for item in loci}
        feature_files = {s: open(util.part_of(os.path.join(part, s, "features.tsv")), "w", newline="\n")
                         for s in sets}
        try:
            for fh in feature_files.values():
                fh.write("\t".join(FEATURE_COLUMNS) + "\n")
            for got, feats in util.ordered(genome_files, items, pool=pool):
                for locus, text in got:
                    if locus in copy_files:
                        copy_files[locus].write(text)
                for locus, lines in feats:
                    s = set_of.get(locus, "busco")
                    if s in feature_files:
                        feature_files[s].write(lines)
        finally:
            for fh in list(copy_files.values()) + list(feature_files.values()):
                fh.close()
    for s in sets:
        path = os.path.join(part, s, "features.tsv")
        os.replace(util.part_of(path), path)

    def of_set(rows, s):
        return [r for r in rows if set_of.get(r["locus"], "busco") == s]
    for s in sets:
        _tsv(os.path.join(part, s, "copies.tsv"), COPY_COLUMNS,
             sorted(of_set(collected["copy_rows"], s), key=lambda r: (r["accession"], order.get(r["locus"], 1 << 20))))
    if "introns" in collected:  # package.genetics (since 2026-11-07)
        from . import codecheck, introns as introns_mod  # noqa: PLC0415
        for s in sets:
            _tsv(os.path.join(part, s, "introns.tsv"), introns_mod.INTRON_COLUMNS,
                 sorted(of_set(collected["introns"], s), key=lambda r: (r["accession"], order.get(r["locus"], 1 << 20),
                                                                         r["start"])))
            _tsv(os.path.join(part, s, "intron_sites.tsv"),
                 introns_mod.SITE_COLUMNS + tuple(collected["intron_site_columns"]),
                 sorted(of_set(collected["intron_sites"], s), key=lambda r: (order.get(r["locus"], 1 << 20),
                                                                              int(r["hmm_position"]), int(r["phase"]))))
        os.makedirs(os.path.join(part, "metadata"), exist_ok=True)
        path = os.path.join(part, "metadata", "codons.tsv.gz")
        with util.bgzf_open(util.part_of(path)) as fh:
            fh.write("\t".join(codecheck.CODON_COLUMNS) + "\n")
            fh.writelines("\t".join(str(r[c]) for c in codecheck.CODON_COLUMNS) + "\n" for r in collected["codons"])
        os.replace(util.part_of(path), path)
    n_records = collections.Counter()
    if rrna:
        from . import rrna as rrna_mod  # noqa: PLC0415
        rrna_mod.write_records(os.path.join(part, "rrna"), rrna)
        for _name, _seq, row in rrna:
            n_records[int(row["species_taxid"])] += 1
    from . import organelle as organelle_mod  # noqa: PLC0415
    organelle = organelle or []
    organelle_mod.write_genes(os.path.join(part, "organelle"), organelle)
    n_organelle = collections.Counter(int(row["species_taxid"]) for _n, _s, row in organelle)
    # every species' genome's copy numbers of its loci, rRNA and organelles in one table (since 2026-11-08)
    copy_number_rows = copy_numbers(genomes, collected["copy_rows"], loci, set_of, organelle)
    write_gz(os.path.join(part, COPY_NUMBERS), COPY_NUMBER_COLUMNS, copy_number_rows)
    for rel, columns, table_rows in tables:
        os.makedirs(os.path.dirname(os.path.join(part, rel)) or part, exist_ok=True)
        _tsv(os.path.join(part, rel), columns, table_rows)
    rows = []
    for acc, r in sorted(genomes.items(), key=lambda kv: int(kv[1]["species_taxid"])):
        if r["representative"] != "t":
            continue
        sp = int(r["species_taxid"])
        rows.append(dict(taxid=sp, name=ncbi.name.get(sp, r["organism"]), representative=acc, source=r["source"],
                         loci=r["loci"], rrna_records=n_records.get(sp, 0), organelle_genes=n_organelle.get(sp, 0),
                         type_material=r.get("type_material", "-")))
    held = {int(r["taxid"]) for r in rows}
    rows += [dict(taxid=sp, name=ncbi.name.get(sp, str(sp)), representative="-", source="rrna", loci=0,
                  rrna_records=n, organelle_genes=n_organelle.get(sp, 0))
             for sp, n in sorted(n_records.items()) if sp not in held]
    from . import labels as labels_mod  # noqa: PLC0415
    named = labels_mod.label_species(ncbi, [int(r["taxid"]) for r in rows], references)
    for r in rows:
        sp = int(r["taxid"])
        r.update(supergroup=labels_mod.supergroup(ncbi, sp), group=labels_mod.group(ncbi, sp),
                 ncbi_lineage=labels_mod.ncbi_lineage(ncbi, sp),
                 ncbi_lineage_filled=labels_mod.ncbi_lineage_filled(ncbi, sp))
        for ref in REFERENCES:
            lineage, how = named[sp].get(ref, ("-", "-"))
            r[ref.lower()], r[f"{ref.lower()}_match"] = lineage, how
    _tsv(os.path.join(part, "species.tsv"), SPECIES_COLUMNS, rows)
    taxa = ncbi.write_subset(os.path.join(part, "taxonomy"),
                             set(collected["taxids"]) | {int(r["species_taxid"]) for r in genomes.values()} |
                             set(n_records) | ({int(row["taxid"]) for _n, _s, row in rrna} if rrna else set()) |
                             {int(row["taxid"]) for _n, _s, row in organelle})
    for writer in writers:
        writer(part)
    described = metadata(part, genomes) if metadata else None
    # the package's tables are written: its benchmarks (benchmarks.run) read them, before the manifest lists them all
    measured = before_manifest(part) if before_manifest else None
    files = {}
    for dirpath, _dirs, names in os.walk(part):
        for name in sorted(names):
            path = os.path.join(dirpath, name)
            files[os.path.relpath(path, part).replace(os.sep, "/")] = [os.path.getsize(path), util.md5(path)]
    statuses = collections.Counter(r["status"] for r in collected["copy_rows"])
    inputs = dict(inputs or {}, references={ref.name: ref.stamp() for ref in references})
    manifest = {"format": FORMAT, "version": version, "eukore": __version__,
                "built": datetime.datetime.now().isoformat(timespec="seconds"),
                "inputs": inputs,
                "counts": {"species": sum(1 for r in rows if r["representative"] != "-"),
                           "species_rrna_only": sum(1 for r in rows if r["representative"] == "-"),
                           "genomes": len(genomes), "copies": collected.get("copy_count", 0),
                           "copy_status": dict(statuses), "rrna_records": sum(n_records.values()),
                           "organelle_genes": {k: dict(collections.Counter(row["gene"] for _n, _s, row in organelle
                                                                           if row["organelle"] == k))
                                               for k in organelle_mod.KINDS},
                           "species_organelle_genes": len(n_organelle), "taxa": taxa,
                           "copy_numbers": dict(collections.Counter(r["set"] for r in copy_number_rows)),
                           "busco_loci": sum(1 for s in set_of.values() if s == "busco"),
                           "marker_loci": sum(1 for s in set_of.values() if s == "markers"),
                           "copy_confidence": dict(collections.Counter(r["confidence"] for r in collected["copy_rows"]
                                                                       if r.get("confidence"))),
                           # of the unverified, the lone lineages' (since 2026-11-10)
                           "copy_unverified": unverified_alone(collected["copy_rows"], genomes, ncbi, set_of),
                           "rrna_confidence": dict(collections.Counter(row.get("confidence", "-")
                                                                       for _n, _s, row in rrna or ())),
                           "labels": {ref.lower(): dict(collections.Counter(r[f"{ref.lower()}_match"] for r in rows))
                                      for ref in REFERENCES}},
                "files": dict(sorted(files.items()))}
    if described is not None:
        manifest["counts"]["metadata"] = described
    if measured:
        manifest["benchmarks"] = measured
    with open(os.path.join(part, MANIFEST), "w", newline="\n") as fh:
        json.dump(manifest, fh, indent=1)
    return put_in_place(part, out)


def put_in_place(part, out):
    """A finished package folder (part) renamed to out, an earlier package there moved aside first and deleted after
    -> out. A killed run's <out>.old is deleted first (until 2026-11-06 the rename onto it failed, after hours)."""
    old = out + ".old"
    shutil.rmtree(old, ignore_errors=True)
    if os.path.isdir(out):
        os.replace(out, old)
    os.replace(part, out)
    shutil.rmtree(old, ignore_errors=True)
    return out


def read(folder):
    """A package's manifest, or None where there is none (unfinished)."""
    try:
        with open(os.path.join(folder, MANIFEST)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def check(folder, deep=False):
    """Why a package is not whole (no manifest, another format, a file missing or of another size; deep: or another
    MD5, or, since 2026-11-10, a BGZF file without its end block: cut short before its MD5 was taken, a node that
    went down after its rename), or None. A .gz of a package before BGZF (2026-11-09) is not BGZF: not judged so."""
    manifest = read(folder)
    if manifest is None:
        return f"{folder}: no Eukore package (no {MANIFEST})"
    if manifest.get("format") not in READABLE:
        return f"{folder}: format {manifest.get('format')}, this reader takes {', '.join(map(str, READABLE))}"
    for rel, (size, digest) in manifest["files"].items():
        path = os.path.join(folder, rel)
        if not os.path.isfile(path):
            return f"{folder}: {rel} is missing"
        if os.path.getsize(path) != size or (deep and util.md5(path) != digest):
            return f"{folder}: {rel} is not the file its manifest names"
        if deep and rel.endswith(".gz") and util.is_bgzf(path) and not util.bgzf_whole(path):
            return f"{folder}: {rel} is cut short (no BGZF end block)"
    return None


def read_table(folder, name):
    """A package table (genomes.tsv, species.tsv, busco/loci.tsv, busco/copies.tsv, busco/features.tsv,
    rrna/records.tsv, organelle/genes.tsv) -> [row dicts]."""
    with util.open_text(os.path.join(folder, name)) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_copies(folder, locus, kind="copies", locus_set="busco"):
    """A locus's copies (kind "copies": whole; "cds": their coding sequences; "proteins": their proteins; the last two
    since format 5; locus_set: busco, or markers for Eukore's own loci) -> {accession: sequence}."""
    out, name = {}, None
    ext = "faa" if kind == "proteins" else "fna"
    with gzip.open(os.path.join(folder, locus_set, kind, f"{locus}.{ext}.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            elif name:
                out[name] = line.strip()
    return out
