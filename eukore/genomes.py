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
import gzip
import json
import os
import re
import shutil
import subprocess
import time

from . import loci as loci_mod
from . import odb as odb_mod
from . import placement, quality, species, util

EXCLUDE = {"contaminated", "unverified source organism", "misassembled", "chimeric", "low quality sequence"} \
    | quality.EXCLUDE
# NCBI's interspecies hybrids ("Saccharomyces cerevisiae x Saccharomyces eubayanus") have species rank: their genomes
# hold both parents' copies of nearly every locus (the 2026-11-05 trial: 24-30 of 30 Duplicated), and their reads are
# their parents'. A named hybrid species (Saccharomyces pastorianus) stays, its duplicated loci flagged.
HYBRID = re.compile(r"\s[x×]\s")
CATEGORY = {"reference genome": 0, "representative genome": 1}
LEVEL = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
LOCI_COLUMNS = ("locus", "status", "contig", "contig_length", "start", "end", "strand", "length", "cds_length",
                "exons", "identity", "coverage", "hmm_score", "cutoff", "hmm_coverage", "copies", "query",
                "protein_length", "window_class")


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


def jobs(ncbi, assemblies, refseq_rows=(), model_rows=(), inside=None, no_assembly=None):
    """The genome set -> {accession: job}: job = dict(accession (GenBank's), fetch (the accession fetched), url, size,
    species_taxid, taxid, organism, why [ncbi | model | no_assembly], lineage (ancestor taxids), itsx (ITSx's sets),
    rdna_unit, flags, others (the species' other assemblies' accessions), path (no_assembly: its records' FASTA), gencode
    (its NCBI genetic code))."""
    refseq = {r["paired"]: r for r in refseq_rows if r.get("paired", "na").startswith("GCA_")}
    from . import rrna  # noqa: PLC0415
    out = {}

    def add(row, why, sp, others=(), unit="-", path=None):
        fetch, url = download_of(row, refseq) if path is None else (row["accession"], None)
        if url is None and path is None:
            return
        sp = int(sp)
        out[row["accession"]] = dict(
            accession=row["accession"], fetch=fetch, url=url, size=row.get("size"), species_taxid=sp,
            taxid=int(row.get("taxid") or sp), organism=row.get("organism", ""), why=[why],
            lineage=list(ncbi.ancestors(ncbi.merged.get(sp, sp))), itsx=rrna.itsx_groups(ncbi, sp), rdna_unit=unit,
            flags=sorted(row.get("flags") or ()), others=[o["accession"] for o in others], path=path,
            level=row.get("level", ""), category=row.get("category", ""),
            gencode=ncbi.code_of(int(row.get("taxid") or sp)))

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
        add(best, "ncbi", sp, others)
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


RESULT_FILES = ("copies.fna", "loci.tsv", "features.tsv", "windows.fna.gz")
HELD_OUT = "held_out.tsv"  # a dataset species' BUSCOs called without its own proteins (the OrthoDB benchmark's)


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


_ODB = {}


def odb_of(settings):
    """The worker's dataset (read once per process)."""
    key = settings["odb"]
    if key not in _ODB:
        _ODB.clear()
        _ODB[key] = odb_mod.Odb(settings["odb"], settings["loci"])
    return _ODB[key]


def tiara_classes(windows_fasta, tiara, threads=1, folder=None):
    """Tiara's class of each window (eukarya, bacteria, archaea, organelle, unknown, ...) -> {window: class}."""
    if not tiara or not os.path.isfile(windows_fasta):
        return {}
    out_path = os.path.join(folder or os.path.dirname(windows_fasta), "tiara.tsv")
    run = subprocess.run([tiara, "-i", windows_fasta, "-o", out_path, "-t", str(max(1, threads)), "-m", "1000"],
                         capture_output=True, text=True, check=False)
    if run.returncode != 0 or not os.path.isfile(out_path):
        return {}
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
        name = r["window"].split("|")
        if len(name) < 3 or name[0] != locus or name[1] != contig:
            continue
        w1 = int(name[2].split("-")[0])
        a, b = w1 + r["start"] - 1, w1 + r["end"] - 1
        if a <= int(end) and b >= int(start):
            return f"{r['action']}:{r['division']}:{r['top_taxon']}"
    return "-"


def run(args):
    """One genome (a job of jobs(), with settings): fetched into settings["genomes"] unless there (or settings
    "staged_only": not fetched, left for a later run), its BUSCO genes called (loci.call with odb's queries for its
    lineage), its rRNA extracted (settings["rrna"]: a function of (genome path, job, out folder) -> its summary),
    the windows around its copies classified (settings["tiara"]), a dataset species' BUSCOs called again without its
    own proteins (settings["odb_own"], HELD_OUT), its results written into WORK/done/<accession>/
    (RESULT_FILES, then done.json) -> (accession, summary or None, why it is not done or None)."""
    job, settings = args
    acc = job["accession"]
    folder = done_dir(settings["work"], acc)
    path = job.get("path") or genome_file(settings["genomes"], acc)
    if not os.path.isfile(path):
        if settings.get("staged_only"):
            return acc, None, "not fetched yet"
        os.makedirs(settings["genomes"], exist_ok=True)
        why = util.fetch_file(job["url"], path)
        if why:
            return acc, None, f"fetch failed: {why}"
    started = time.time()  # the run's own time, the fetch left out (benchmarks.py's runtime.tsv)
    tmp = util.part_of(folder)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    odb = odb_of(settings)
    queries = odb.queries(job["lineage"], settings["odb_lineages"])
    try:
        got = loci_mod.call(path, queries, odb, settings.get("miniprot", "miniprot"),
                            settings.get("hmmsearch", "hmmsearch"), settings.get("threads", 1), tmp,
                            settings.get("hmm_db"), job.get("gencode", 1))
    except (OSError, RuntimeError) as e:
        shutil.rmtree(tmp, ignore_errors=True)
        return acc, None, f"locus calling failed: {e}"
    if got["failed"]:
        shutil.rmtree(tmp, ignore_errors=True)
        return acc, None, got["failed"]
    windows = os.path.join(tmp, "windows.fna.gz")
    with gzip.open(windows, "wt", compresslevel=6) as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in got["windows"].items())
    classes = {}
    if settings.get("tiara"):
        plain = os.path.join(tmp, "windows.fna")
        with open(plain, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in got["windows"].items())
        classes = tiara_classes(plain, settings["tiara"], settings.get("threads", 1), tmp)
        os.remove(plain)
    window_of = {n.split("|", 1)[0]: n for n in got["windows"]}
    taken = 0
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
            taken += 1
            fa.write(f">{locus}\n{r['copy']}\n")
            exons = sum(1 for f in r["features"] if f[0] == "exon")
            tb.write("\t".join(map(str, [
                locus, r["status"], r["contig"], r["contig_length"], r["start"], r["end"], r["strand"],
                len(r["copy"]), len(r["spliced"]), exons, f"{r['identity']:.4f}", f"{r['coverage']:.3f}",
                f"{r['hmm_score']:.1f}", f"{r['cutoff']:.1f}", f"{r['hmm_coverage']:.3f}", r["copies"], r["query"],
                len(r["protein"]), classes.get(window_of.get(locus, ""), "-")])) + "\n")
            ft.writelines(f"{locus}\t{k}\t{s}\t{e}\n" for k, s, e in r["features"])
    rrna = None
    if settings.get("rrna") and not job.get("path"):  # records of a species without an assembly: no genome's rRNA
        try:
            rrna = settings["rrna"](path, job, tmp)
        except (OSError, ValueError, EOFError, RuntimeError) as e:  # RuntimeError: barrnap's, ITSx's
            shutil.rmtree(tmp, ignore_errors=True)
            return acc, None, f"rRNA extraction failed: {e}"
    seconds = round(time.time() - started, 1)
    own = own_of(job, settings)
    if own:  # its species is in BUSCO's dataset: its BUSCOs called again without its own proteins (HELD_OUT)
        why = held_out(path, job, settings, odb, own, tmp)
        if why:
            shutil.rmtree(tmp, ignore_errors=True)
            return acc, None, why
    summary = {"accession": acc, "fetched": job.get("fetch", acc), "why": job["why"], "loci": taken,
               "searched": list(odb.loci),
               "status": dict(collections.Counter(r["status"] for r in got["loci"].values())), "rrna": rrna,
               "bytes": os.path.getsize(path), "seconds": seconds, "threads": settings.get("threads", 1)}
    if own:
        summary["held_out"] = {"left_out": own, "seconds": round(time.time() - started - seconds, 1)}
    with open(os.path.join(tmp, "done.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.replace(tmp, folder)
    return acc, summary, None


def own_of(job, settings):
    """The dataset taxids of a job's species (settings["odb_own"], odb.own_taxids), [] where the dataset has none."""
    return (settings.get("odb_own") or {}).get(job["species_taxid"], [])


def held_out_missing(work, job, settings):
    """Whether a done genome of a dataset species lacks its held-out calls (a run before 2026-11-06)."""
    return bool(own_of(job, settings)) and not os.path.isfile(os.path.join(done_dir(work, job["accession"]),
                                                                           HELD_OUT))


def held_out(path, job, settings, odb, own, folder):
    """A dataset species' genome's BUSCOs called again with its own proteins (own: their dataset taxids) left out of
    the queries, into folder/HELD_OUT (locus, status, query, cds: the call's coding sequence): the OrthoDB benchmark
    compares these with OrthoDB's gene models, which the build's own calls were queried with
    -> why it failed, or None."""
    queries = odb.queries(job["lineage"], settings["odb_lineages"], leave_out=own)
    try:
        got = loci_mod.call(path, queries, odb, settings.get("miniprot", "miniprot"),
                            settings.get("hmmsearch", "hmmsearch"), settings.get("threads", 1), folder,
                            settings.get("hmm_db"), job.get("gencode", 1))
    except (OSError, RuntimeError) as e:
        return f"held-out locus calling failed: {e}"
    if got["failed"]:
        return f"held-out locus calling failed: {got['failed']}"
    with open(os.path.join(folder, HELD_OUT), "w", newline="\n") as fh:
        fh.write("locus\tstatus\tquery\tcds\n")
        for locus in odb.loci:
            r = got["loci"].get(locus, {"status": "Missing"})
            fh.write(f"{locus}\t{r['status']}\t{r.get('query', '-')}\t{r.get('spliced', '-')}\n")
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
