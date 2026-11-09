"""Eukore's genome set (since 2026-11-05: NCBI's own assemblies, no EukDetect): which eukaryote genome stands for each
species, one genome's run (fetched, its BUSCO genes called, its rRNA extracted, its results kept in the work folder),
and the prune of the kept genomes to a budget.

The species are NCBI's species with a GenBank or RefSeq assembly in the groups asked for (sources.ASSEMBLY_GROUPS:
fungi, protozoa, plant, invertebrate, vertebrate_mammalian, vertebrate_other), but NCBI's placeholders (uncultured
...) and interspecies hybrids (HYBRID). A species' genome is, by rank():
- not one NCBI flags as contaminated, of an unverified source organism, misassembled, chimeric or of low quality
  (EXCLUDE: never taken); one with another flag (NCBI's: too large or small for its species, a metagenome's, a
  single cell's, ...; "partial genome": NCBI's genome_rep; SMALL: under a quarter of its clade's typical genome) only
  where the species has nothing better; NCBI's own FCS-GX and ANI reports add their flags (quality.py: mostly another
  organism, a failed species check, never; contaminated, an inconclusive check, demoted);
- then a RefSeq reference, then a representative genome; then the better assembly level; then one FCS-GX finds no
  contaminant in; then the fewer contigs per megabase; then the newer.
A model species (data/model_species.tsv) takes its reference assembly whatever the ranking (why "model"). Where a
genome has a RefSeq copy (gbrs_paired_asm), that copy is fetched: RefSeq screens eukaryote genomes for contamination
(FCS-GX) before it takes them. A species with no assembly takes its records in data/no_assembly.fna.gz (why
"no_assembly": EukDetect2's genomes from EukProt and its Dientamoeba fragilis transcriptome, its only remnant).

A genome is done when WORK/done/<accession>/done.json is written (last); a run calls and extracts only the genomes not
done, so that one killed goes on where it stopped.
"""

import collections
import concurrent.futures
import contextlib
import gzip
import json
import os
import re
import shutil
import socket
import subprocess
import tempfile
import threading
import time

from . import annotation as annotation_mod
from . import loci as loci_mod
from . import metadata as metadata_mod
from . import odb as odb_mod
from . import placement, quality, species, util

EXCLUDE = {"contaminated", "unverified source organism", "misassembled", "chimeric", "low quality sequence"} \
    | quality.EXCLUDE
# NCBI's interspecies hybrids ("Saccharomyces cerevisiae x Saccharomyces eubayanus") have species rank: their genomes
# hold both parents' copies of nearly every locus (the 2026-11-05 trial: 24-30 of 30 Duplicated), and their reads are
# their parents'. A named hybrid species (Saccharomyces pastorianus; NCBI's "Genus x epithet": Fragaria x ananassa,
# Mentha x piperita) stays, its duplicated loci flagged: a cross names two species before its x (until 2026-11-06 any
# " x " did, and the strawberry was left out).
HYBRID = re.compile(r"^\S+\s+\S+.*\s[x×]\s")
CATEGORY = {"reference genome": 0, "representative genome": 1}
LEVEL = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
LOCI_COLUMNS = ("locus", "status", "contig", "contig_length", "start", "end", "strand", "length", "cds_length",
                "exons", "identity", "coverage", "hmm_score", "cutoff", "hmm_coverage", "copies", "query",
                "protein_length", "window_class", "rescued",  # rescued: loci.py's rescue pass (since 2026-11-06)
                "model",  # model: miniprot's, or the annotation's source (loci.annotated, since 2026-11-06)
                "dup_identity",  # a Duplicated locus' nearest other copy (loci.dup_identity, since 2026-11-07)
                "cds_check")  # its coding sequence against the protein confirmed (loci.cds_state, since 2026-11-09)
# A Duplicated locus' other complete copies kept whole (loci.ALTERNATES; since 2026-11-09, genome record 2): the
# package takes one where the chosen copy's lineage does not fit (package.collect). Their rows are LOCI_COLUMNS' with
# `alt` (1, 2, ...) first; their copies >locus|alt; their features with `alt` after the locus.
ALTERNATES = ("alternates.tsv", "alternates.fna", "alternates_features.tsv")
# The genome record (since 2026-11-09): what made a done genome's calls, in its done.json's `genome_record`. A change to
# what a genome's run makes (its calls, copies, features) that a package should not mix with older runs' is a new
# record: a done genome of an older one is run again (fetched again where pruned), unless the build's --reversion
# stamps it with this one (its calls then kept as they are). A done.json without one is of record 1.
# 1: the runs until 2026-11-09.
# 2: the coding sequence checked against the protein confirmed and recovered (loci.cds_state, cds_check); a region's
#    clean alignment first; a Duplicated locus' other complete copies kept whole (ALTERNATES); an annotated model
#    with a stop inside, or a fusion, not taken.
GENOME_RECORD = 2


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def genome_file(folder, accession):
    return os.path.join(folder, f"{accession}.fna.gz")


def rank(row):
    """A species' assemblies' order (module docstring); lower is better."""
    flags = set(row.get("flags") or ())
    size_mb = (row.get("size") or 0) / 1e6
    per_mb = (row.get("contigs") or 1 << 20) / max(size_mb, 1.0)
    return (bool(flags), row.get("category") == "derived from metagenome", CATEGORY.get(row.get("refseq_category"), 2),
            LEVEL.get(row.get("level"), 4), bool(row.get("fcs_bp")), per_mb,
            -int(re.sub(r"\D", "", row.get("date") or "") or 0), row["accession"])


def representatives(assemblies, ncbi, inside=None):
    """{NCBI species taxid: (best assembly row, the species' other rows)} of the assemblies (sources.read_assemblies
    rows of GenBank's tables; a RefSeq row paired with a GenBank one is that one's copy), EXCLUDE'd ones and NCBI's
    placeholder species left out; inside: a test of a species taxid (--species)."""
    by_species = collections.defaultdict(list)
    for row in assemblies:
        if set(row.get("flags") or ()) & EXCLUDE:
            continue
        sp = ncbi.species_of(row["species_taxid"])
        if ncbi.rank.get(sp) != "species" or ncbi.is_placeholder(sp) or HYBRID.search(ncbi.name.get(sp, "")) \
                or (inside and not inside(sp)):
            continue
        by_species[sp].append(row)
    by_species = flag_small(by_species, ncbi)
    out = {}
    for sp, rows in by_species.items():
        rows = sorted(rows, key=rank)
        out[sp] = (rows[0], rows[1:])
    return out


# A species' assembly under SMALL_SHARE of its clade's typical genome (the median of the other species' median sizes in
# its genus, else family, else order, of CLADE_MIN species or more) is a part of its genome, whatever NCBI calls it: the
# 2026-11-05 panel's "Full" vertebrate assemblies of 2-50 Mb (an 11 Mb Macaca silenus, a 241 Mb snake) held 0-25 of
# the 129 BUSCOs. Flagged SMALL (demoted by rank(): taken only where the species has nothing better).
SMALL = "small for its clade"
SMALL_SHARE = 0.25
CLADE_MIN = 3


def flag_small(by_species, ncbi):
    """{species: rows} -> the same, rows under SMALL_SHARE of their clade's median genome size flagged SMALL."""
    def median(xs):
        xs = sorted(xs)
        return xs[len(xs) // 2]
    own = {sp: median([r["size"] for r in rows if r.get("size")]) for sp, rows in by_species.items()
           if any(r.get("size") for r in rows)}
    clade_of = {}
    members = collections.defaultdict(list)
    for sp in own:
        at = {ncbi.rank.get(t): t for t in ncbi.ancestors(sp)}
        clade_of[sp] = [at.get(r) for r in ("genus", "family", "order")]
        for t in clade_of[sp]:
            if t:
                members[t].append(sp)
    out = {}
    for sp, rows in by_species.items():
        typical = None
        for t in clade_of.get(sp, ()):
            others = [own[o] for o in members.get(t, ()) if o != sp]
            if t and len(others) >= CLADE_MIN:
                typical = median(others)
                break
        if typical:
            rows = [dict(r, flags=sorted(set(r.get("flags") or ()) | {SMALL}))
                    if r.get("size") and r["size"] < SMALL_SHARE * typical else r for r in rows]
        out[sp] = rows
    return out


def download_of(row, refseq):
    """The accession and URL a genome is fetched from: its RefSeq copy's where it has one (refseq: {GenBank accession:
    RefSeq row}), else its own."""
    copy = refseq.get(row["accession"])
    use = copy or row
    return use["accession"], placement.genome_url(use)


def job_of(ncbi, row, why, sp, refseq, others=(), unit="-", path=None):
    """One assembly's job (jobs()'s), or None where it cannot be fetched."""
    from . import rrna  # noqa: PLC0415
    fetch, url = download_of(row, refseq) if path is None else (row["accession"], None)
    if url is None and path is None:
        return None
    sp = int(sp)
    genes = (refseq.get(row["accession"]) or row).get("genes")  # the fetched assembly's annotation
    return dict(
        accession=row["accession"], fetch=fetch, url=url, size=row.get("size"), species_taxid=sp,
        taxid=int(row.get("taxid") or sp), organism=row.get("organism", ""), why=[why],
        lineage=list(ncbi.ancestors(ncbi.merged.get(sp, sp))), itsx=rrna.itsx_groups(ncbi, sp), rdna_unit=unit,
        flags=sorted(row.get("flags") or ()), others=[o["accession"] for o in others], path=path,
        level=row.get("level", ""), category=row.get("category", ""),
        refseq_category=(refseq.get(row["accession"]) or row).get("refseq_category", ""),
        gencode=ncbi.code_of(int(row.get("taxid") or sp)), annotated=None if genes is None else genes > 0)


def jobs(ncbi, assemblies, refseq_rows=(), model_rows=(), inside=None, no_assembly=None):
    """The genome set -> {accession: job}: job = dict(accession (GenBank's), fetch (the accession fetched), url, size,
    species_taxid, taxid, organism, why [ncbi | model | no_assembly], lineage (ancestor taxids), itsx (ITSx's sets),
    rdna_unit, flags, others (the species' other assemblies' accessions), path (no_assembly: its records' FASTA), gencode
    (its NCBI genetic code), annotated (the fetched assembly has NCBI gene models, annotation.py; None: not known))."""
    refseq = {r["paired"]: r for r in refseq_rows if r.get("paired", "na").startswith("GCA_")}
    out = {}

    def add(row, why, sp, others=(), unit="-", path=None):
        job = job_of(ncbi, row, why, sp, refseq, others, unit, path)
        if job:
            out[row["accession"]] = job

    reps = representatives(assemblies, ncbi, inside)
    model_species = {}
    for m in model_rows:
        sp = ncbi.species_of(int(m["ncbi_taxid"]))
        if inside and not inside(sp):
            continue
        model_species[sp] = m
    for sp, (best, others) in sorted(reps.items()):
        if sp in model_species:
            continue
        # the best that can be fetched (until 2026-11-06 a best row without an FTP path lost the species)
        ranked = [best, *others]
        chosen = next((r for r in ranked if download_of(r, refseq)[1]), None)
        if chosen is not None:
            add(chosen, "ncbi", sp, [r for r in ranked if r is not chosen])
    for sp, m in sorted(model_species.items()):
        row = species.assembly_row(m)
        others = reps.get(sp, (None, []))
        rest = [r for r in ([others[0]] if others[0] else []) + list(others[1]) if r["accession"] != row["accession"]]
        add(row, "model", sp, rest, unit=m.get("rdna_unit", "-"))
    for acc, (sp, path, organism) in sorted((no_assembly or {}).items()):
        sp = ncbi.species_of(sp)
        if sp in reps or sp in model_species or (inside and not inside(sp)):
            continue
        add(dict(accession=acc, taxid=sp, organism=organism), "no_assembly", sp, path=path)
    return out


# NCBI's Transcriptome Shotgun Assemblies (TSA) as ENA mirrors them (since 2026-11-07, format 5): INSDC's TSA sets, their
# NCBI taxids, sizes and FASTA files, listed by ENA's portal (NCBI's own TSA lists name no taxids)
TSA_FIELDS = ("accession", "tax_id", "scientific_name", "base_count", "first_public", "strain", "environmental_sample",
              "set_fasta_ftp")
TSA_URL = ("https://www.ebi.ac.uk/ena/portal/api/search?result=tsa_set&query=tax_tree(2759)&format=tsv&limit=0"
           f"&fields={','.join(TSA_FIELDS)}")
TSA_FILE = "tsa_sets.tsv"


def tsa_path(src_dir):
    return os.path.join(src_dir, "tsa", TSA_FILE)


def fetch_tsa(src_dir, download, log=print, refresh=False):
    """ENA's list of eukaryote TSA sets into src_dir/tsa (kept when there, unless refresh)."""
    path = tsa_path(src_dir)
    if os.path.isfile(path) and os.path.getsize(path) > 0 and not refresh:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    log("TSA: fetching ENA's list of eukaryote transcriptome sets")
    download(TSA_URL, path)
    return path


def tsa_state(src_dir):
    """ENA's TSA list a build read, for the package manifest's inputs.tsa (None without one): its file, MD5, date
    fetched and sets."""
    path = tsa_path(src_dir)
    if not os.path.isfile(path):
        return None
    import datetime  # noqa: PLC0415
    fetched = datetime.datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds")
    return dict(file=TSA_FILE, md5=util.md5(path), fetched=fetched, sets=len(read_tsa(path)))


def read_tsa(path):
    """ENA's TSA list -> [row dicts] of its sets with a FASTA and a taxid, not of an environmental sample."""
    if not path or not os.path.isfile(path):
        return []
    with util.open_text(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    return [r for r in rows if r.get("set_fasta_ftp") and r.get("tax_id", "").isdigit()
            and r.get("environmental_sample", "false").lower() != "true"]


def tsa_jobs(ncbi, rows, job_map, inside=None):
    """The species with no genome (no NCBI, model or TOGA2 job; EukDetect2's records alone, or nothing) and a TSA set ->
    {set accession: job} (why "tsa", kind "transcriptome"): of each species, its set of the most bases, then the
    newest. A species' no_assembly job is then dropped from job_map: its transcriptome stands for it."""
    held = {int(j["species_taxid"]) for j in job_map.values() if j["why"][0] != "no_assembly"}
    best = {}
    for r in rows:
        sp = ncbi.species_of(int(r["tax_id"]))
        if ncbi.rank.get(sp) != "species" or sp in held or ncbi.is_placeholder(sp) or (inside and not inside(sp)):
            continue
        rank = (int(r["base_count"]) if r.get("base_count", "").isdigit() else 0, r.get("first_public", ""))
        if sp not in best or rank > best[sp][0]:
            best[sp] = (rank, r)
    out = {}
    from . import rrna  # noqa: PLC0415
    for sp, (_rank, r) in sorted(best.items()):
        url = "https://" + r["set_fasta_ftp"].split(";")[0].removeprefix("ftp://")
        out[r["accession"]] = dict(
            accession=r["accession"], fetch=r["accession"], url=url, size=int(r["base_count"] or 0) or None,
            species_taxid=sp, taxid=int(r["tax_id"]), organism=r.get("scientific_name", ""), why=["tsa"],
            kind="transcriptome", lineage=list(ncbi.ancestors(sp)), itsx=rrna.itsx_groups(ncbi, sp), rdna_unit="-",
            flags=[], others=[], path=None, level="transcriptome", category="transcriptome", refseq_category="",
            gencode=ncbi.code_of(int(r["tax_id"])), annotated=False, strain=r.get("strain", ""))
    for acc in [a for a, j in job_map.items() if j["why"][0] == "no_assembly"
                and int(j["species_taxid"]) in {j2["species_taxid"] for j2 in out.values()}]:
        del job_map[acc]
    return out


PANEL = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "intraspecific_panel.tsv")


def read_panel(path=PANEL):
    """The intraspecific panel (data/intraspecific_panel.tsv) -> {species taxid: extra assemblies wanted, or None for
    --intraspecific's number}."""
    out = {}
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    head = lines[0]
    for f in lines[1:]:
        r = dict(zip(head, f))
        out[int(r["ncbi_taxid"])] = int(r["extras"]) if r.get("extras", "-").isdigit() else None
    return out


def intraspecific_jobs(ncbi, assemblies, refseq_rows, job_map, panel, extras=3, inside=None):
    """The panel species' further assemblies (since 2026-11-07; the barcode gap, intraspecific.py) -> {accession:
    job} (why "intraspecific"): of each panel species with a genome in job_map, up to its number (panel's, else extras)
    of its other assemblies, by Eukore's ranking (NCBI's flagged ones last; never one it excludes), each of another
    BioSample, strain and isolate than the species' genome and those taken. They are runs of their own: job_map keeps
    one genome per species."""
    refseq = {r["paired"]: r for r in refseq_rows if r.get("paired", "na").startswith("GCA_")}
    reps = representatives(assemblies, ncbi, inside)
    by_species = {int(j["species_taxid"]): j for j in job_map.values()}
    by_accession = {r["accession"]: r for r in assemblies}
    out = {}
    for sp, n in sorted(panel.items()):
        sp = ncbi.species_of(sp)
        rep = by_species.get(sp)
        n = extras if n is None else min(n, extras) if extras else 0
        if rep is None or sp not in reps or not n:
            continue
        seen = [by_accession.get(rep["accession"], {})]
        best, rest = reps[sp]
        taken = 0
        for row in [best, *rest]:
            if taken >= n:
                break
            if row["accession"] == rep["accession"] or any(same_sample(row, s) for s in seen):
                continue
            job = job_of(ncbi, row, "intraspecific", sp, refseq)
            if job:
                out[row["accession"]] = job
                seen.append(row)
                taken += 1
    return out


def same_sample(a, b):
    """Whether two assemblies are of one sample: the same BioSample, strain or isolate (where named)."""
    for k in ("biosample", "strain", "isolate"):
        x, y = (a.get(k) or "").strip().lower(), (b.get(k) or "").strip().lower()
        if x and x not in ("na", "-") and x == y:
            return True
    return False


MAMMALIA = 40674


def toga2_jobs(ncbi, rows, job_map, inside=None, groups=("vertebrate_mammalian", "vertebrate_other")):
    """TOGA2's assemblies outside NCBI's set (annotation.read_toga2: DNA Zoo's, Bat1K's, the Hiller lab's; since
    2026-11-06) of the vertebrate species job_map has no genome of -> {accession TOGA2_<assembly>: job} (why "toga2"):
    of a species' the one of the longest scaffold N50, its genome TOGA2's .2bit (fetch_genome), its gene models TOGA2's
    (job["annotations"]). groups: NCBI's assembly groups wanted (a mammal's vertebrate_mammalian, else
    vertebrate_other)."""
    have = {int(j["species_taxid"]) for j in job_map.values()}
    best = {}
    for r in rows:
        if r.get("genome", "-") == "-":
            continue
        sp = ncbi.species_of(ncbi.merged.get(int(r["species_taxid"]), int(r["species_taxid"])))
        if sp in have or sp not in ncbi.parent or ncbi.rank.get(sp) != "species" or (inside and not inside(sp)):
            continue
        if ("vertebrate_mammalian" if MAMMALIA in ncbi.ancestors(sp) else "vertebrate_other") not in groups:
            continue
        if sp not in best or int(r["scaffold_n50"]) > int(best[sp]["scaffold_n50"]):
            best[sp] = r
    return {r["accession"]: toga2_job(ncbi, sp, r) for sp, r in sorted(best.items())}


def toga2_job(ncbi, sp, r, others=()):
    """A job of TOGA2's assembly r (annotation.read_toga2's row) of species sp."""
    from . import rrna  # noqa: PLC0415
    return dict(accession=r["accession"], fetch=r["accession"], url=r["genome"], size=None, species_taxid=sp, taxid=sp,
                organism=r["organism"], why=["toga2"], lineage=list(ncbi.ancestors(sp)), itsx=rrna.itsx_groups(ncbi, sp),
                rdna_unit="-", flags=[], others=list(others), path=None, level="", category="isolate",
                gencode=ncbi.code_of(sp), annotated=False, annotations=[("toga2", r["annotation"])],
                assembly_source=r["assembly_source"])


TOGA2_BETTER = 2.0  # TOGA2's assembly takes the place of NCBI's chosen one at this many times its contig N50 or more,
TOGA2_LENGTH = 0.9  # and at this share of its total length or more (not a part of the genome, however contiguous)


def prefer_toga2(ncbi, rows, job_map, ncbi_stats, toga2_length):
    """The NCBI genomes to replace by TOGA2's assemblies outside NCBI (since 2026-11-06) -> {NCBI accession: TOGA2
    job}: a species' TOGA2 assembly of the longest contig N50 where it is TOGA2_BETTER times that of the NCBI genome
    chosen for it and its total length TOGA2_LENGTH of that genome's or more (ncbi_stats(job): NCBI's assembly stats,
    dict(contig_n50, total_length); toga2_length(row): the sum of its chrom.sizes; either None where not known: NCBI's
    kept). A model species' genome stays. The swap gives up NCBI's contamination screens of the genome and its NCBI
    gene models (TOGA2's come with it); the NCBI genome is listed among the species' others."""
    by_species = {int(j["species_taxid"]): (a, j) for a, j in job_map.items()}
    best = {}
    for r in rows:
        if r.get("genome", "-") == "-" or r.get("ncbi_accession", "-") != "-" or not r.get("contig_n50", "").isdigit():
            continue
        sp = ncbi.species_of(ncbi.merged.get(int(r["species_taxid"]), int(r["species_taxid"])))
        acc, job = by_species.get(sp, (None, None))
        if not job or job["why"] != ["ncbi"]:
            continue
        if sp not in best or int(r["contig_n50"]) > int(best[sp]["contig_n50"]):
            best[sp] = r
    out = {}
    for sp, r in sorted(best.items()):
        acc, job = by_species[sp]
        theirs = ncbi_stats(job) or {}
        if not theirs.get("contig_n50") or not theirs.get("total_length"):
            continue
        if int(r["contig_n50"]) < TOGA2_BETTER * theirs["contig_n50"]:
            continue
        length = toga2_length(r)
        if length and length >= TOGA2_LENGTH * theirs["total_length"]:
            out[acc] = toga2_job(ncbi, sp, r, [acc] + list(job.get("others", ())))
    return out


def stats_url(job):
    """NCBI's assembly statistics beside a job's genome (its contig N50, its total length)."""
    url = job.get("url") or ""
    return url[:-len("_genomic.fna.gz")] + "_assembly_stats.txt" if url.endswith("_genomic.fna.gz") else None


def assembly_stats(stats):
    """NCBI's assembly statistics file -> dict(contig_n50, total_length) (None each where not given), or None."""
    if not stats or not os.path.isfile(stats):
        return None
    out = dict(contig_n50=None, total_length=None)
    with open(stats) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[:4] == ["all"] * 4 and len(f) > 5 and f[5].isdigit() and f[4] in ("contig-N50", "total-length"):
                out[f[4].replace("-N50", "_n50").replace("-", "_")] = int(f[5])
    return out


def chrom_sizes_url(row):
    """TOGA2's chrom.sizes beside its genome (annotation.read_toga2's row)."""
    return row["genome"].rsplit("/", 1)[0] + "/chrom.sizes" if row.get("genome", "-") != "-" else None


def total_length(chrom_sizes):
    """A chrom.sizes file (name, length) -> the lengths' sum, or None."""
    if not chrom_sizes or not os.path.isfile(chrom_sizes):
        return None
    with open(chrom_sizes) as fh:
        sizes = [line.split() for line in fh if line.strip()]
    return sum(int(f[1]) for f in sizes if len(f) > 1 and f[1].isdigit()) or None


def fetch_genome(job, path, **how):
    """A job's genome into path (gzipped FASTA): NCBI's as it is, TOGA2's .2bit converted (twobit.py) -> "" if it
    arrived, else why not (util.fetch_file's). how: util.fetch_file's scratch and checker."""
    url = job["url"]
    if not url.endswith(".2bit"):
        return util.fetch_file(url, path, **how)
    from . import twobit  # noqa: PLC0415
    # the .2bit in --local_tmp where given (since 2026-11-08: beside the genome, a killed conversion's left there for
    # good; the fetch stage now deletes leftovers of a day)
    raw = os.path.join(how.get("scratch") or os.path.dirname(path), os.path.basename(path) + ".2bit")
    why = util.fetch_file(url, raw, gzip_check=False, **{k: v for k, v in how.items() if k == "times"})
    if why:
        return why
    try:
        twobit.to_fasta(raw, path)
    except (OSError, ValueError, IndexError) as e:
        return f"not a readable .2bit ({e})"
    finally:
        os.remove(raw)
    return ""


FATES = ("chosen", "chosen, flagged", "model species", "no assembly (records)", "TOGA2's genome",
         "every assembly excluded", "placeholder", "hybrid", "not a species", "not in the taxdump", "no download")


def fates(ncbi, assemblies, job_map, inside=None):
    """What became of each NCBI taxon with a GenBank assembly (the benchmarks' selection funnel) -> {taxid: fate}:
    chosen (its genome unflagged, or flagged: the best it has), a model species' reference, no assembly (EukProt's
    records), every assembly excluded (EXCLUDE: NCBI's or its FCS-GX and ANI reports'), or left out as a placeholder,
    an interspecies hybrid, a taxon above species (assemblies labelled by a genus) or one the taxdump does not hold (a
    taxdump older than the assembly tables: until 2026-11-06 called "not a species")."""
    out = {}
    for job in job_map.values():
        why = job["why"][0]
        out[int(job["species_taxid"])] = ("model species" if why == "model" else "no assembly (records)"
                                          if why == "no_assembly" else "TOGA2's genome" if why == "toga2"
                                          else "chosen, flagged" if job.get("flags") else "chosen")
    excluded = collections.defaultdict(list)
    for row in assemblies:
        sp = ncbi.species_of(row["species_taxid"])
        if sp in out or (inside and not inside(sp)):
            continue
        excluded[sp].append(bool(set(row.get("flags") or ()) & EXCLUDE))
    for sp, flagged in excluded.items():
        out[sp] = ("not in the taxdump" if sp not in ncbi.parent else "not a species" if ncbi.rank.get(sp) != "species"
                   else "placeholder" if ncbi.is_placeholder(sp)
                   else "hybrid" if HYBRID.search(ncbi.name.get(sp, "")) else "every assembly excluded"
                   if all(flagged) else "no download")  # its genomes without an FTP path
    return out


def no_assembly_genomes(fasta, folder, ncbi):
    """data/no_assembly.fna.gz (records >taxid|organism|locus|source) as one FASTA per species in folder ->
    {accession EUK_<taxid>: (taxid, path, organism)}."""
    if not fasta or not os.path.isfile(fasta):
        return {}
    os.makedirs(folder, exist_ok=True)
    by_taxid = collections.defaultdict(list)
    names = {}
    for name, seq in odb_mod.read_fasta(fasta):
        f = name.split("|")
        if f and f[0].isdigit():
            by_taxid[int(f[0])].append((name, seq))
            names[int(f[0])] = f[1].replace("_", " ") if len(f) > 1 else ""
    out = {}
    for taxid, records in by_taxid.items():
        path = os.path.join(folder, f"EUK_{taxid}.fna")
        if not os.path.isfile(path):
            with open(util.part_of(path), "w", newline="\n") as fh:
                fh.writelines(f">{n}\n{s}\n" for n, s in records)
            os.replace(util.part_of(path), path)
        out[f"EUK_{taxid}"] = (taxid, path, names.get(taxid) or ncbi.name.get(taxid, ""))
    return out


def read_rows(path):
    """A table (plain or zstd; header first) -> its rows, split."""
    if not path or not os.path.isfile(path):
        return []
    with util.open_text(path) as fh:
        fh.readline()
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


HELD_OUT = "held_out.tsv"  # a dataset species' BUSCOs called without its own proteins (the OrthoDB benchmark's)
HELD_OUT_COLUMNS = ("locus", "status", "query", "contig", "start", "end", "rescued", "cds")  # place since 2026-11-06


def done_dir(work, accession):
    return os.path.join(work, "done", safe(accession))


def searched(work, accession):
    """The BUSCOs a done genome was searched for (done.json's `searched`; in a run before format 4, which searched
    euk30's 30, its loci.tsv's loci, Missing ones too), or None where it is not done."""
    try:
        with open(os.path.join(done_dir(work, accession), "done.json")) as fh:
            summary = json.load(fh)
    except (OSError, ValueError):
        return None
    if "searched" in summary:
        return set(summary["searched"])
    return {r[0] for r in read_rows(os.path.join(done_dir(work, accession), "loci.tsv"))}


def is_done(work, accession, loci=None):
    """Whether a genome's run is done; with loci, done for all of them (format 4 searches every BUSCO: a genome
    searched for euk30's 30 alone is run again)."""
    if loci is None:
        return os.path.isfile(os.path.join(done_dir(work, accession), "done.json"))
    got = searched(work, accession)
    return got is not None and set(loci) <= got


def record_of(summary):
    """A done genome's genome record (done.json's `genome_record`; 1 for a run before 2026-11-09), or None for none."""
    return None if summary is None else int(summary.get("genome_record", 1))


def record_stale(work, accession):
    """Whether a done genome's calls are of an older genome record than this Eukore's (GENOME_RECORD): run again."""
    got = record_of(done_summary(work, accession))
    return got is not None and got < GENOME_RECORD


def reversion(work, accession):
    """A done genome's done.json stamped with this Eukore's genome record, its calls kept as they are (--reversion;
    since 2026-11-09): `reversioned` notes the record it had and when -> that record, or None where it was not done or
    of this record already (or newer)."""
    path = os.path.join(done_dir(work, accession), "done.json")
    summary = done_summary(work, accession)
    old = record_of(summary)
    if old is None or old >= GENOME_RECORD:
        return None
    import datetime  # noqa: PLC0415
    from . import __version__  # noqa: PLC0415
    summary["genome_record"] = GENOME_RECORD
    summary.setdefault("reversioned", []).append(dict(
        record=old, to=GENOME_RECORD, eukore=__version__, at=datetime.datetime.now().isoformat(timespec="seconds")))
    with open(util.part_of(path), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    os.replace(util.part_of(path), path)
    return old


DONE_ELSEWHERE = "done by another build"  # finished by another build on the WORK since this build's plan (done_stamp)


def done_stamp(work, accession):
    """A genome's done.json's modification time (ns), or None where it is not done: what a plan saw, so that its run
    leaves a genome another build finished meanwhile alone (genomes.run)."""
    try:
        return os.stat(os.path.join(done_dir(work, accession), "done.json")).st_mtime_ns
    except OSError:
        return None


_ODB = {}


def odb_of(settings):
    """The worker's dataset (read once per process)."""
    key = (settings["odb"], settings.get("markers"), tuple(settings["loci"]))
    if key not in _ODB:
        _ODB.clear()
        _ODB[key] = odb_mod.Odb(settings["odb"], settings["loci"],
                                [settings["markers"]] if settings.get("markers") else [])
    return _ODB[key]


def tiara_classes(windows_fasta, tiara, threads=1, folder=None):
    """Tiara's class of each window (eukarya, bacteria, archaea, organelle, unknown, ...) -> {window: class}; {} without
    tiara or windows, None where tiara failed (since 2026-11-08: {}, the genome then done unclassified for good)."""
    if not tiara or not os.path.isfile(windows_fasta):
        return {}
    out_path = os.path.join(folder or os.path.dirname(windows_fasta), "tiara.tsv")
    run = subprocess.run([tiara, "-i", windows_fasta, "-o", out_path, "-t", str(max(1, threads)), "-m", "1000"],
                         capture_output=True, text=True, check=False)
    if run.returncode != 0 or not os.path.isfile(out_path):
        return None
    out = {}
    with open(out_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(head, line.rstrip("\n").split("\t")))
            name = f.get("sequence_id") or line.split("\t")[0]
            out[name.split()[0]] = f.get("class_fst_stage", "-")
    return out


FCS_COLUMNS = ("window", "start", "end", "action", "division", "top_taxon")


def fcs_screen(work, accession, taxid, gx_db, run_gx="run_gx", threads=1):
    """NCBI's FCS-GX on a done genome's windows (since 2026-11-05; run_gx of bioconda's ncbi-fcs-gx, its database of
    ~470 GB best in /dev/shm): each range it calls another organism's -> WORK/done/<accession>/fcs.tsv (FCS_COLUMNS:
    the window, the range on it 1-based, FCS-GX's action EXCLUDE/TRIM/FIX/REVIEW, its division and top taxon) -> the
    ranges, or None where it could not run (no windows, run_gx failed)."""
    folder = done_dir(work, accession)
    windows = os.path.join(folder, "windows.fna.gz")
    if not os.path.isfile(windows):
        return None
    import tempfile  # noqa: PLC0415
    with tempfile.TemporaryDirectory(dir=folder) as tmp:
        names, fasta = {}, os.path.join(tmp, "w.fna")
        with gzip.open(windows, "rt") as src, open(fasta, "w") as out:
            for line in src:  # short ids: FCS-GX reads NCBI-style ids with | as fields
                if line.startswith(">"):
                    names[f"w{len(names)}"] = line[1:].strip()
                    out.write(f">w{len(names) - 1}\n")
                else:
                    out.write(line)
        run = subprocess.run([run_gx, "--fasta", fasta, "--tax-id", str(taxid), "--gx-db", gx_db, "--out-dir", tmp],
                             capture_output=True, text=True, check=False,
                             env=dict(os.environ, GX_NUM_CORES=str(max(1, threads))))
        report = next((os.path.join(tmp, f) for f in os.listdir(tmp) if f.endswith("fcs_gx_report.txt")), None)
        if run.returncode != 0 or report is None:
            return None
        ranges = []
        with open(report) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5 and f[0] in names:
                    ranges.append(dict(window=names[f[0]], start=int(f[1]), end=int(f[2]), action=f[4],
                                       division=f[5] if len(f) > 5 else "-", top_taxon=f[7] if len(f) > 7 else "-"))
    path = os.path.join(folder, "fcs.tsv")
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(FCS_COLUMNS) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in FCS_COLUMNS) + "\n" for r in ranges)
    os.replace(util.part_of(path), path)
    return ranges


def fcs_of_copy(ranges, locus, contig, start, end):
    """The FCS-GX call over a copy (its contig span, 1-based): ACTION:division:taxon of the first range on its window
    overlapping it, or "-"."""
    for r in ranges:
        # locus|contig|w1-w2, the contig's name may hold | itself (a TSA's ENA|GKLH01011938|GKLH01011938.1: until
        # 2026-11-08 no transcriptome's copy met its FCS-GX call)
        head, _bar, span = r["window"].rpartition("|")
        name, _bar, on = head.partition("|")
        if name != locus or on != contig or "-" not in span:
            continue
        w1 = int(span.split("-")[0])
        a, b = w1 + r["start"] - 1, w1 + r["end"] - 1
        if a <= int(end) and b >= int(start):
            return f"{r['action']}:{r['division']}:{r['top_taxon']}"
    return "-"


RUNNING_ELSEWHERE = "running elsewhere"  # another build's worker holds the genome (claim)
# a claim of another host untouched this long (its heartbeat stopped: a job killed, e.g. at its time limit, or a
# node that died) is taken over: 5 beats, less than a restart's minutes of planning before its first genome (until
# 2026-11-08 2 days, the killed job's genomes, its largest, skipped by its restart as run by another build)
CLAIM_STALE = 300
CLAIM_BEAT = 60           # seconds between a running genome's claim touches


def claim_path(work, accession):
    return os.path.join(work, "done", f".{safe(accession)}.claim")


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except (PermissionError, OSError):
        return True
    return True


def claim_live(path):
    """Whether a genome's claim file is held: by a process alive on this host, or by another host's touched within
    CLAIM_STALE (a run touches its claim every CLAIM_BEAT seconds)."""
    try:
        with open(path) as fh:
            host, pid = fh.read().split()
        age = time.time() - os.path.getmtime(path)
    except FileNotFoundError:
        return False
    except (OSError, ValueError):  # being written this instant, or damaged: held while young
        try:
            return time.time() - os.path.getmtime(path) < 60
        except OSError:
            return False
    if host == socket.gethostname():
        return _alive(int(pid))
    return age < CLAIM_STALE


def claim(work, accession):
    """The genome claimed for this process (since 2026-11-08: two builds on one WORK, e.g. several nodes, never run the
    same genome; until then nothing stopped them) -> the claim's path, or None where another holds it. A claim of a
    dead process (claim_live) is taken over."""
    path = claim_path(work, accession)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for _ in range(2):
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            if claim_live(path):
                return None
            with contextlib.suppress(FileNotFoundError):
                os.remove(path)
            continue
        with os.fdopen(fd, "w") as fh:
            fh.write(f"{socket.gethostname()}\t{os.getpid()}\n")
        return path
    return None


@contextlib.contextmanager
def claimed(work, accession):
    """A genome claimed while the block runs (yields its claim's path, or None where another holds it), its claim
    touched every CLAIM_BEAT seconds and removed after."""
    path = claim(work, accession)
    if path is None:
        yield None
        return
    stop = threading.Event()

    def beat():
        while not stop.wait(CLAIM_BEAT):
            with contextlib.suppress(OSError):
                os.utime(path)
    heart = threading.Thread(target=beat, daemon=True)
    heart.start()
    try:
        yield path
    finally:
        stop.set()
        heart.join()
        with contextlib.suppress(OSError):
            os.remove(path)


def stale_partials(work):
    """The killed runs' partial done folders under WORK/done (util.part_of's .<host>-<pid>.<accession>.partial), not
    those of a run going on (since 2026-11-08: every one was deleted, another node's running genomes' too): a partial
    is kept while a live claim (claim_live) names its genome."""
    import glob  # noqa: PLC0415
    live = []
    for c in glob.glob(os.path.join(work, "done", ".*.claim")):
        if claim_live(c):
            live.append(os.path.basename(c)[1:-len(".claim")])
    return [p for p in glob.glob(os.path.join(work, "done", ".*.partial"))
            if not any(os.path.basename(p).endswith(f".{a}.partial") for a in live)]


def run(args):
    """One genome (a job of jobs(), with settings): fetched into settings["genomes"] unless there (or settings
    "staged_only": not fetched, left for a later run), claimed (claim: another build's worker may hold it), its BUSCO
    genes called (loci.call with odb's queries for its lineage), its rRNA extracted (settings["rrna"]: a function of
    (genome path, job, out folder) -> its summary), the windows around its copies classified (settings["tiara"]), a
    dataset species' BUSCOs called again without its own proteins (settings["odb_own"], HELD_OUT), its assembly's
    metadata measured (metadata.describe: stats.json, organelle.fna.gz, and sketch.sig.gz by settings["sourmash"]),
    its results written into WORK/done/<accession>/ (copies.fna, loci.tsv, features.tsv, windows.fna.gz, then
    done.json) -> (accession, summary or None, why it is not done or None).
    The rRNA and the metadata, which read the genome alone on one core each, run beside the locus calls (since
    2026-11-08: after them, the job's other threads idle). A run again for better gene models alone (job["rerun"]
    "annotation") keeps the done folder's rRNA and held-out calls, which do not depend on them."""
    job, settings = args
    acc = job["accession"]
    path = job.get("path") or genome_file(settings["genomes"], acc)
    if not os.path.isfile(path) and settings.get("staged_only"):
        return acc, None, "not fetched yet"
    with claimed(settings["work"], acc) as held:
        if held is None:
            return acc, None, RUNNING_ELSEWHERE
        # another build finished it since this one's plan noted its done.json (since 2026-11-08: run again whole,
        # that build's done folder replaced)
        if "done_stamp" in job and done_stamp(settings["work"], acc) != job["done_stamp"]:
            return acc, None, DONE_ELSEWHERE
        if not os.path.isfile(path):
            os.makedirs(settings["genomes"], exist_ok=True)
            why = fetch_genome(job, path)
            if why:
                return acc, None, f"fetch failed: {why}"
        folder = done_dir(settings["work"], acc)
        tmp = util.part_of(folder)
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        # the run's temporary files (miniprot's pieces of a large genome, the rescue's ORFs, tiara's and vsearch's
        # files) and a copy of the genome, read there by every step, in settings["scratch"] where given (--local_tmp: a
        # node's local disk; since 2026-11-08: all in WORK, the genome read from OUT some six times), else in tmp
        scratch = scratch_dir(settings, acc) or tmp
        side = concurrent.futures.ThreadPoolExecutor(2)
        try:
            if scratch != tmp and not job.get("path"):
                local = os.path.join(scratch, os.path.basename(path))
                shutil.copyfile(path, local)
            else:
                local = path
            got = _run(job, settings, local, folder, tmp, side, scratch)
        except BaseException:
            side.shutdown(wait=True, cancel_futures=True)
            shutil.rmtree(tmp, ignore_errors=True)
            raise
        finally:
            if scratch != tmp:
                shutil.rmtree(scratch, ignore_errors=True)
        side.shutdown(wait=True, cancel_futures=True)
        if got[2]:
            shutil.rmtree(tmp, ignore_errors=True)
            # its three alone (until 2026-11-09 _run's four: the build stopped at the first genome that failed)
            return got[:3]
        if os.path.isdir(folder):
            shutil.rmtree(folder)
        os.replace(tmp, folder)
        gff, report = got[3]
        for used in (gff, report):  # annotation.none_file stays: a run again knows none has one
            if used and os.path.isfile(used):  # compared and taken: annotation.tsv holds what the benchmark needs
                os.remove(used)
        return got[:3]


SCRATCH_PREFIX = ".eukore-"


def scratch_dir(settings, accession):
    """A genome run's own folder in settings["scratch"] (--local_tmp), or None."""
    if not settings.get("scratch"):
        return None
    os.makedirs(settings["scratch"], exist_ok=True)
    return tempfile.mkdtemp(prefix=f"{SCRATCH_PREFIX}{safe(accession)}-", dir=settings["scratch"])


def stale_scratch(work, scratch):
    """The killed runs' folders in scratch (scratch_dir's), not those of a run going on (its claim live)."""
    import glob  # noqa: PLC0415
    out = []
    for path in glob.glob(os.path.join(scratch, SCRATCH_PREFIX + "*")):
        acc = os.path.basename(path)[len(SCRATCH_PREFIX):].rpartition("-")[0]
        if not claim_live(claim_path(work, acc)):
            out.append(path)
    return out


def _run(job, settings, path, folder, tmp, side, scratch=None):
    """genomes.run's work in its partial folder tmp, side: a pool of two threads (the rRNA and the metadata), its
    temporary files in scratch (default tmp) -> (accession, summary or None, why not or None, (GFF, report) taken)."""
    scratch = scratch or tmp
    acc = job["accession"]
    started = time.time()  # the run's own time, the fetch left out (benchmarks.py's runtime.tsv)
    before = done_summary(settings["work"], acc) if job.get("rerun") == "annotation" else None
    odb = odb_of(settings)
    rrna_job = None
    # records of a species without an assembly, a transcriptome (since 2026-11-07): no genome's rRNA
    if settings.get("rrna") and not job.get("path") and job.get("kind") != "transcriptome":
        kept = os.path.join(folder, "rrna.tsv")
        if before is not None and os.path.isfile(kept):  # a run for its gene models alone: its rRNA as it was
            from . import rrna as rrna_mod  # noqa: PLC0415
            shutil.copyfile(kept, os.path.join(tmp, "rrna.tsv"))
            rrna_job = rrna_mod.counts(rrna_mod.load(tmp))
        else:
            rrna_job = side.submit(settings["rrna"], path, job, tmp)
    meta_job = None
    if not job.get("path"):  # the assembly's metadata while it is here (metadata.py); a failure the pass tries again
        for name in (metadata_mod.STATS, metadata_mod.ORGANELLE, metadata_mod.SKETCH):  # a run again: its own kept
            if os.path.isfile(os.path.join(folder, name)):
                shutil.copyfile(os.path.join(folder, name), os.path.join(tmp, name))
        want = metadata_mod.missing(tmp, settings.get("sourmash"))
        if want:
            meta_job = side.submit(metadata_mod.describe, path, tmp, acc, job.get("fetch", acc),
                                   settings.get("sourmash"), want)
    queries = odb.queries(job["lineage"], settings["odb_lineages"])
    gff, report, source = annotation_of(job, settings)
    taken = (gff, report)
    try:
        got = loci_mod.call(path, queries, odb, settings.get("miniprot", "miniprot"),
                            settings.get("hmmsearch", "hmmsearch"), settings.get("threads", 1), scratch,
                            settings.get("hmm_db"), job.get("gencode", 1), **rescue_of(settings),
                            annotation=(gff, report, source) if gff else None)
    except (OSError, RuntimeError) as e:
        return acc, None, f"locus calling failed: {e}", taken
    if got["failed"]:
        return acc, None, got["failed"], taken
    if gff is not False:  # miniprot's models against the annotation's (annotation.py's benchmark; no rows: none)
        annotation_mod.annotate(tmp, gff, report, {
            locus: dict(r["miniprot"], status=r["status"], contig=r["contig"], strand=r["strand"])
            for locus, r in got["loci"].items() if r["status"] != "Missing"}, source=source)
    windows = os.path.join(tmp, "windows.fna.gz")
    with util.bgzf_open(windows, level=util.WORK_LEVEL) as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in got["windows"].items())
    classes, tiara = {}, "-"
    if job.get("kind") != "transcriptome":  # a transcript's flanks are no window
        tiara = "none"
        if settings.get("tiara"):
            classes = classify(got["windows"], settings["tiara"], settings.get("threads", 1), scratch)
            tiara = "failed" if classes is None else "done"
            classes = classes or {}
    window_of = {n.split("|", 1)[0]: n for n in got["windows"]}
    dup = loci_mod.dup_identity(got["loci"], settings.get("vsearch"), settings.get("threads", 1), scratch)
    taken_loci = 0

    def row_of(locus, r, window_class, dup_identity):
        exons = sum(1 for f in r["features"] if f[0] == "exon")
        return [locus, r["status"], r["contig"], r["contig_length"], r["start"], r["end"], r["strand"], len(r["copy"]),
                len(r["spliced"]), exons, f"{r['identity']:.4f}", f"{r['coverage']:.3f}", f"{r['hmm_score']:.1f}",
                f"{r['cutoff']:.1f}", f"{r['hmm_coverage']:.3f}", r["copies"], r["query"], len(r["protein"]),
                window_class, "t" if r.get("rescued") else "f", r.get("model", "miniprot"), dup_identity,
                r.get("cds_check", "-")]
    alternates = [(locus, k, dict(alt, status=r["status"], copies=r["copies"]))
                  for locus in odb.loci for r in [got["loci"].get(locus, {})]
                  for k, alt in enumerate(r.get("alternates") or (), 1)]
    with open(os.path.join(tmp, "copies.fna"), "w", newline="\n") as fa, \
            open(os.path.join(tmp, "loci.tsv"), "w", newline="\n") as tb, \
            open(os.path.join(tmp, "features.tsv"), "w", newline="\n") as ft:
        tb.write("\t".join(LOCI_COLUMNS) + "\n")
        ft.write("locus\tkind\tstart\tend\n")
        for locus in odb.loci:
            r = got["loci"].get(locus, {"status": "Missing"})
            if r["status"] == "Missing":
                tb.write("\t".join([locus, "Missing"] + ["-"] * (len(LOCI_COLUMNS) - 2)) + "\n")
                continue
            taken_loci += 1
            fa.write(f">{locus}\n{r['copy']}\n")
            tb.write("\t".join(map(str, row_of(locus, r, classes.get(window_of.get(locus, ""), "-"),
                                               dup.get(locus, "-")))) + "\n")
            ft.writelines(f"{locus}\t{k}\t{s}\t{e}\n" for k, s, e in r["features"])
    if alternates:  # a Duplicated locus' other complete copies whole (ALTERNATES; since 2026-11-09)
        table, fasta, features = (os.path.join(tmp, n) for n in ALTERNATES)
        with open(table, "w", newline="\n") as tb, open(fasta, "w", newline="\n") as fa, \
                open(features, "w", newline="\n") as ft:
            tb.write("\t".join(("alt",) + LOCI_COLUMNS) + "\n")
            ft.write("locus\talt\tkind\tstart\tend\n")
            for locus, k, r in alternates:
                tb.write("\t".join(map(str, [k] + row_of(locus, r, "-", "-"))) + "\n")
                fa.write(f">{locus}|{k}\n{r['copy']}\n")
                ft.writelines(f"{locus}\t{k}\t{kind}\t{s}\t{e}\n" for kind, s, e in r["features"])
    rrna = rrna_job
    if isinstance(rrna_job, concurrent.futures.Future):
        try:
            rrna = rrna_job.result()
        except (OSError, ValueError, EOFError, RuntimeError) as e:  # RuntimeError: barrnap's, ITSx's
            return acc, None, f"rRNA extraction failed: {e}", taken
    seconds = round(time.time() - started, 1)
    own = own_of(job, settings) if job.get("kind") != "transcriptome" else []
    held = None
    if own:  # its species is in BUSCO's dataset: its BUSCOs called again without its own proteins (HELD_OUT)
        kept = os.path.join(folder, HELD_OUT)
        if before is not None and before.get("held_out") and {r[0] for r in read_rows(kept)} == set(odb.loci):
            shutil.copyfile(kept, os.path.join(tmp, HELD_OUT))  # a run for its gene models alone: as it was
            held = before["held_out"]
        else:
            why = held_out(path, job, settings, odb, own, tmp, scratch)
            if why:
                return acc, None, why, taken
            held = {"left_out": own, "seconds": round(time.time() - started - seconds, 1)}
    from . import __version__  # noqa: PLC0415
    summary = {"accession": acc, "fetched": job.get("fetch", acc), "why": job["why"], "loci": taken_loci,
               "genome_record": GENOME_RECORD, "eukore": __version__,  # what made its calls (since 2026-11-09)
               "cds_check": dict(collections.Counter(r.get("cds_check", "-") for r in got["loci"].values()
                                                     if r["status"] != "Missing")),
               "searched": list(odb.loci),
               "status": dict(collections.Counter(r["status"] for r in got["loci"].values())), "rrna": rrna,
               "bytes": os.path.getsize(path), "seconds": seconds, "threads": settings.get("threads", 1),
               "rescue": rescue_of(settings)["rescue_missing"],
               "rescued": sum(1 for r in got["loci"].values() if r.get("rescued")),
               "annotated": bool(gff), "annotation": source if gff else "-",
               # the sources asked before the one taken (all, where none had one): none (better_annotation)
               "annotation_none": [] if gff is False else annotation_mod.none_before(job, source if gff else None),
               "annotated_models": sum(1 for r in got["loci"].values() if r.get("model", "miniprot") != "miniprot"),
               "tiara": tiara}  # done, failed (backfill tries again), none (no tiara), - (a transcriptome)
    if held is not None:
        summary["held_out"] = held
    if meta_job is not None:
        why = meta_job.result()
        if why:
            summary["metadata_failed"] = why
    with open(os.path.join(tmp, "done.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    return acc, summary, None, taken


def classify(windows, tiara, threads, folder):
    """Tiara's classes of a genome's windows {name: sequence} (tiara_classes) -> {window: class}, or None where it
    failed."""
    plain = os.path.join(folder, "windows.fna")
    with open(plain, "w") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in windows.items())
    try:
        return tiara_classes(plain, tiara, threads, folder)
    finally:
        os.remove(plain)


TIARA_COLUMN = LOCI_COLUMNS.index("window_class")


def backfill(args):
    """A done genome's steps that a tool missing at its run left out, made now the tool is here (since 2026-11-08; a
    pool's task: (job, settings)) -> (accession, [steps made], why one failed or None):
    - tiara's classes of its windows, where its run had no tiara or tiara failed (done.json's `tiara`; a run before
      2026-11-08: its loci.tsv's window_class all `-`), into loci.tsv;
    - its rRNA (settings["rrna"]), where its run had no barrnap (no rrna.tsv) and its genome is here.
    Its done.json says what was made. A genome another build holds is left alone."""
    job, settings = args
    acc = job["accession"]
    folder = done_dir(settings["work"], acc)
    summary = done_summary(settings["work"], acc)
    if summary is None:
        return acc, [], None
    todo = backfill_needs(job, settings, summary)
    if not todo:
        return acc, [], None
    with claimed(settings["work"], acc) as held:
        if held is None:
            return acc, [], None
        made, why = [], None
        if "tiara" in todo:
            with tempfile.TemporaryDirectory(dir=folder) as tmp:
                with gzip.open(os.path.join(folder, "windows.fna.gz"), "rt") as fh:
                    names = [line[1:].strip() for line in fh if line.startswith(">")]
                windows = os.path.join(tmp, "windows.fna")
                with gzip.open(os.path.join(folder, "windows.fna.gz"), "rb") as src, open(windows, "wb") as out:
                    shutil.copyfileobj(src, out)
                classes = tiara_classes(windows, settings["tiara"], settings.get("threads", 1), tmp)
            if classes is None:
                why = "tiara failed"
                summary["tiara"] = "failed"
            else:
                window_of = {n.split("|", 1)[0]: n for n in names}
                table = os.path.join(folder, "loci.tsv")
                with open(table) as fh:
                    lines = fh.readlines()
                with open(util.part_of(table), "w", newline="\n") as fh:
                    fh.write(lines[0])
                    for line in lines[1:]:
                        f = line.rstrip("\n").split("\t")
                        if f[1] != "Missing":
                            f[TIARA_COLUMN] = classes.get(window_of.get(f[0], ""), "-")
                        fh.write("\t".join(f) + "\n")
                os.replace(util.part_of(table), table)
                summary["tiara"] = "done"
                made.append("tiara")
        if "rrna" in todo:
            path = genome_file(settings["genomes"], acc)
            with tempfile.TemporaryDirectory(dir=folder) as tmp:
                try:
                    summary["rrna"] = settings["rrna"](path, job, tmp)
                except (OSError, ValueError, EOFError, RuntimeError) as e:
                    why = why or f"rRNA extraction failed: {e}"
                else:
                    os.replace(os.path.join(tmp, "rrna.tsv"), os.path.join(folder, "rrna.tsv"))
                    made.append("rrna")
        path = os.path.join(folder, "done.json")
        with open(util.part_of(path), "w", newline="\n") as fh:
            json.dump(summary, fh, indent=1)
        os.replace(util.part_of(path), path)
        return acc, made, why


def backfill_needs(job, settings, summary=None):
    """The steps backfill would make of a done genome now: tiara, rrna (backfill)."""
    folder = done_dir(settings["work"], job["accession"])
    summary = summary or done_summary(settings["work"], job["accession"]) or {}
    out = []
    if settings.get("tiara") and job.get("kind") != "transcriptome" \
            and os.path.isfile(os.path.join(folder, "windows.fna.gz")) and summary.get("tiara") != "done":
        called = [r for r in read_rows(os.path.join(folder, "loci.tsv")) if len(r) > TIARA_COLUMN and r[1] != "Missing"]
        state = summary.get("tiara")
        if called and (state in ("none", "failed") or (state is None and all(r[TIARA_COLUMN] == "-" for r in called))):
            out.append("tiara")
    if settings.get("rrna") and not job.get("path") and job.get("kind") != "transcriptome" \
            and not os.path.isfile(os.path.join(folder, "rrna.tsv")) \
            and os.path.isfile(genome_file(settings["genomes"], job["accession"])):
        out.append("rrna")
    return out


def rescue_missing(work, accession, settings):
    """Whether a done genome was run without the rescue pass (done.json's `rescue` false: --no_rescue, or no
    esl-translate) where this run makes it: its run again (since 2026-11-08)."""
    summary = done_summary(work, accession) or {}
    return summary.get("rescue") is False and bool(rescue_of(settings)["rescue_missing"])


def annotation_of(job, settings):
    """A genome's annotation for its run (annotation.py) -> (GFF path, NCBI assembly report path or None, its source):
    its GFF beside the genome (the fetch stage's), else fetched where the run may fetch; (None, None, "-") where no
    source has one (annotation.tsv without rows); (False, None, "-") where the genome has none to try or it is not
    here yet (a staged-only run: compared later, annotate_pass)."""
    if not annotation_mod.wanted(job):
        return False, None, "-"
    gff, source = annotation_mod.found(settings["genomes"], job["accession"])
    report = annotation_mod.report_file(settings["genomes"], job["accession"])
    if not gff:
        if settings.get("staged_only"):  # the fetch stage found none (its note), or it is still to come
            none = annotation_mod.none_noted(settings["genomes"], job)
            return (None, report if os.path.isfile(report) else None, "-") if none else (False, None, "-")
        try:
            gff, source = annotation_mod.fetch(job, settings["genomes"], settings.get("veupathdb_key"))
        except OSError:
            return False, None, "-"
    return gff, report if os.path.isfile(report) else None, source or "-"


def done_summary(work, accession):
    """A done genome's done.json, or None."""
    try:
        with open(os.path.join(done_dir(work, accession), "done.json")) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def annotation_record(work, job):
    """A done genome's gene models -> (the source its run took them from, SOURCES' or "-", {sources found to have
    none}), or None where it is not done. A run before 2026-11-07 says only whether it took NCBI's (done.json's
    `annotated`); an empty annotation.tsv beside it (NCBI's 404, or a GFF without coding sequences) is NCBI's none."""
    summary = done_summary(work, job["accession"])
    if summary is None:
        return None
    if "annotation" in summary:
        return summary["annotation"], set(summary.get("annotation_none") or ())
    table = os.path.join(done_dir(work, job["accession"]), annotation_mod.TABLE)
    empty = os.path.isfile(table) and not read_rows(table)
    if summary.get("annotated") and not (empty and not summary.get("ncbi_models")):
        return "ncbi", set()
    return "-", {"ncbi"} if empty and annotation_mod.gff_url(job) else set()


def source_rank(source):
    return annotation_mod.SOURCES.index(source) if source in annotation_mod.SOURCES else len(annotation_mod.SOURCES)


def better_annotation(work, job):
    """The source of gene models better than those a done genome's run took (since 2026-11-07) -> its name, or None: the
    first of the job's annotation.candidates (NCBI's, then VEuPathDB's, Ensembl's, TOGA2's: annotation.assign's, by
    the build's indexes) not found to have none, where it comes before the source taken. A genome run before an
    index named its assembly, before a VEuPathDB key was given, or whose annotation was not fetched yet, is run again
    with it (eukore_build.py)."""
    got = annotation_record(work, job)
    if got is None:
        return None
    used, none = got
    for source, _url in annotation_mod.candidates(job):
        if source not in none:
            return source if source_rank(source) < source_rank(used) else None
    return None


def note_annotation_none(work, job, sources):
    """Sources found to have no annotation of a done genome added to its done.json's `annotation_none` (the fetch
    stage's, annotate_genome's), so that better_annotation does not ask them again."""
    summary = done_summary(work, job["accession"])
    if summary is None or not sources:
        return
    if "annotation" not in summary:  # a run before 2026-11-07: its record as annotation_record reads it
        from_old = annotation_record(work, job)
        summary["annotation"], summary["annotation_none"] = from_old[0], sorted(from_old[1])
    summary["annotation_none"] = sorted(set(summary.get("annotation_none") or ()) | set(sources))
    path = os.path.join(done_dir(work, job["accession"]), "done.json")
    with open(util.part_of(path), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    os.replace(util.part_of(path), path)


def rerun_for_annotation(work, job, genomes_folder, may_fetch):
    """Whether a done genome is to be run again for better gene models (better_annotation): its genome on disk and
    the better annotation beside it, or one the run may fetch."""
    better = better_annotation(work, job)
    if not better or not os.path.isfile(job.get("path") or genome_file(genomes_folder, job["accession"])):
        return False
    if may_fetch:
        return True
    _path, source = annotation_mod.found(genomes_folder, job["accession"])
    return source is not None and source_rank(source) < source_rank(annotation_record(work, job)[0])


# --- genomes superseded: a species' genome another assembly took the place of ------------------------------------

GENOME_SET = "genome_set.tsv"    # WORK/: each species' genome as the last plan chose it
SUPERSEDED = "superseded.tsv"    # WORK/: the genomes another took the place of (old, new), kept till deleted


def supersede(work, job_map):
    """The plan's genomes against the last plan's (WORK/GENOME_SET; since 2026-11-07) -> [(old accession, new)] of the
    species whose genome changed (a TOGA2 assembly in place of NCBI's; a newer NCBI assembly or version, after
    --refresh), added to WORK/SUPERSEDED. Species outside this plan (a --species build) keep their row."""
    path, chosen = os.path.join(work, GENOME_SET), {}
    if os.path.isfile(path):
        for f in read_rows(path):
            if len(f) >= 2:
                chosen[f[0]] = f[1]
    now = {str(j["species_taxid"]): a for a, j in job_map.items()}
    changed = [(chosen[sp], acc) for sp, acc in sorted(now.items()) if chosen.get(sp, acc) != acc]
    chosen.update(now)
    os.makedirs(work, exist_ok=True)
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("species_taxid\taccession\n" + "".join(f"{sp}\t{a}\n" for sp, a in sorted(chosen.items())))
    os.replace(util.part_of(path), path)
    if changed:
        table = os.path.join(work, SUPERSEDED)
        new = not os.path.isfile(table)
        with open(table, "a", newline="\n") as fh:
            fh.write("old\tnew\n" if new else "")
            fh.writelines(f"{old}\t{acc}\n" for old, acc in changed)
    return changed


def drop_superseded(genomes_folder, work, job_map):
    """The superseded genomes' files deleted (the genome, its annotation and assembly report), each once the genome
    that took its place is done and it is not in the plan again -> [accessions deleted]. Its done folder stays."""
    out = []
    for f in read_rows(os.path.join(work, SUPERSEDED)):
        old, new = f[0], f[1] if len(f) > 1 else None
        if old in job_map or not new or not is_done(work, new):
            continue
        paths = [genome_file(genomes_folder, old), annotation_mod.report_file(genomes_folder, old),
                 annotation_mod.none_file(genomes_folder, old)] + \
            [annotation_mod.gff_file(genomes_folder, old, s) for s in annotation_mod.SOURCES]
        if any(os.path.isfile(p) for p in paths):
            for p in paths:
                if os.path.isfile(p):
                    os.remove(p)
            out.append(old)
    return out


def rescue_of(settings):
    """loci.call's rescue options from the settings: the rescue pass unless settings["rescue"] is False (--no_rescue),
    with settings["esl_translate"]."""
    return dict(rescue_missing=settings.get("rescue", True),
                translate=settings.get("esl_translate") or "esl-translate")


def own_of(job, settings):
    """The dataset taxids of a job's species (settings["odb_own"], odb.own_taxids), [] where the dataset has none."""
    return (settings.get("odb_own") or {}).get(job["species_taxid"], [])


def held_out_missing(work, job, settings):
    """Whether a done genome of a dataset species lacks its held-out calls (a run before 2026-11-06)."""
    return bool(own_of(job, settings)) and not os.path.isfile(os.path.join(done_dir(work, job["accession"]),
                                                                           HELD_OUT))


def held_out(path, job, settings, odb, own, folder, scratch=None):
    """A dataset species' genome's BUSCOs called again with its own proteins (own: their dataset taxids) left out of
    the queries, into folder/HELD_OUT (HELD_OUT_COLUMNS; cds: the call's coding sequence): the OrthoDB benchmark
    compares these with OrthoDB's gene models, which the build's own calls were queried with, and their places with
    the build's own calls
    -> why it failed, or None."""
    queries = odb.queries(job["lineage"], settings["odb_lineages"], leave_out=own)
    try:
        got = loci_mod.call(path, queries, odb, settings.get("miniprot", "miniprot"),
                            settings.get("hmmsearch", "hmmsearch"), settings.get("threads", 1), scratch or folder,
                            settings.get("hmm_db"), job.get("gencode", 1), **rescue_of(settings))
    except (OSError, RuntimeError) as e:
        return f"held-out locus calling failed: {e}"
    if got["failed"]:
        return f"held-out locus calling failed: {got['failed']}"
    with open(os.path.join(folder, HELD_OUT), "w", newline="\n") as fh:
        fh.write("\t".join(HELD_OUT_COLUMNS) + "\n")
        for locus in odb.loci:
            r = got["loci"].get(locus, {"status": "Missing"})
            fh.write("\t".join(str(x) for x in (
                locus, r["status"], r.get("query", "-"), r.get("contig", "-"), r.get("start", "-"), r.get("end", "-"),
                ("t" if r.get("rescued") else "f") if r["status"] != "Missing" else "-", r.get("spliced", "-"))) + "\n")
    return None


def prune(genomes_folder, work, keep_bytes, protected=(), first=()):
    """The kept genomes cut back to keep_bytes once done (since 2026-11-05): `protected` accessions (a host's) kept
    whatever their size, then `first` (the model species'), then the smallest; a genome not done yet kept (pending);
    the rest deleted -> [(accession, bytes, status)], written as genomes_folder/manifest.tsv."""
    files = [f for f in os.listdir(genomes_folder) if f.endswith(".fna.gz")] if os.path.isdir(genomes_folder) else []
    order = []
    for f in files:
        acc = f[:-len(".fna.gz")]
        rank_ = 0 if acc in protected else 1 if acc in first else 2
        order.append((rank_, os.path.getsize(os.path.join(genomes_folder, f)), acc))
    kept, rows = 0, []
    for rank_, size, acc in sorted(order):
        if rank_ == 0 or not is_done(work, acc) or kept + size <= keep_bytes:
            kept += size
            rows.append((acc, size, "kept" if is_done(work, acc) or rank_ == 0 else "pending"))
        else:
            os.remove(genome_file(genomes_folder, acc))
            rows.append((acc, size, "deleted"))
    path = os.path.join(genomes_folder, "manifest.tsv")
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("accession\tbytes\tstatus\n" + "".join(f"{a}\t{s}\t{st}\n" for a, s, st in rows))
    os.replace(util.part_of(path), path)
    return rows


def why_counts(job_list):
    return collections.Counter(w for j in job_list for w in j["why"])
