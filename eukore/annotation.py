"""A genome's own gene models (since 2026-11-06): the BUSCO copies' models where they are good, and a truth set
independent of BUSCO's dataset.

A genome's annotation comes from the first of SOURCES that has one for its assembly (the same accession and version):
- `ncbi`: NCBI's GFF beside the genome on NCBI's server (a RefSeq copy, GCF_, mostly NCBI's Eukaryotic Genome
  Annotation Pipeline's; a GenBank assembly its submitter's: the assembly summary's protein-coding genes,
  sources.read_assemblies' `genes`). A 404, or a GFF without coding sequences, is none.
- `veupathdb`: VEuPathDB's (its organism table gives each genome's accession; its curated annotations of protists'
  and fungi's genomes, newer than Ensembl's copies of some), with a VEuPathDB user's API key alone (since 2026-08 its
  files need one; eukore_build.py --veupathdb_key). Its own sequence names are matched by length (aliases()).
- `ensembl`: Ensembl's (ftp.ebi.ac.uk/pub/ensemblorganisms, its species.json the index): Ensembl's own pipeline's, its
  BRAKER and Helixer runs, VEuPathDB's and community annotations; the biodiversity projects' genomes (Darwin Tree of
  Life, ERGA) RefSeq has not annotated. Of an assembly's builds the newest, Ensembl's own first on a tie.
- `toga2`: TOGA2's (Malovichko et al. 2026; genome.senckenberg.de/download/TOGA2), ~2,500 vertebrate assemblies'
  genes projected from a reference species' (the first of TOGA2_REFERENCES listing the assembly).
On 2026-11-06 VEuPathDB annotates 71 species of NCBI's assemblies that NCBI annotates none of, Ensembl 2,037 (+32%; mostly
invertebrates), TOGA2 907 more vertebrates.

The process stage gives it to loci.call: a copy takes the annotated transcript sharing the most coding bases with
miniprot's model (best_transcript) where the locus' HMM confirms it (loci.annotated; copies.tsv's `model` the
source). And each copy's miniprot model is laid on the annotation's coding sequences (CDS) on its contig
(WORK/done/<accession>/annotation.tsv, ANNOTATION_COLUMNS): its exons and introns, the transcript's, those exact on
both, and the coding bases each holds and shares. The introns' exact boundaries (splice sites) are the measure least
touched by start and stop codons (miniprot's CDS and an annotation's end differently at a stop).

The annotation is fetched by the fetch stage beside the genome (OUT/genomes/<accession>.gff.gz, or
<accession>.<source>.gff.gz), with NCBI's assembly report (its sequences' names: RefSeq's, GenBank's, the chromosomes'
that Ensembl uses, the unversioned accessions TOGA2 uses), or by the process stage where it may fetch, and deleted once
the genome is done. A genome without one gets an annotation.tsv without rows, so that it is not asked for again. A
genome done before its annotation was there is compared later (annotate_genome; its copies keep miniprot's models).
benchmarks.py gathers the tables (benchmarks/annotation.tsv).
"""

import collections
import json
import os
import re
import urllib.parse

from . import util

TABLE = "annotation.tsv"
ANNOTATION_COLUMNS = ("locus", "status", "contig", "source", "transcript", "gene", "exons", "ann_exons", "exact_exons",
                      "introns", "ann_introns", "exact_introns", "bases", "ann_bases", "shared_bases")
SOURCES = ("ncbi", "veupathdb", "ensembl", "toga2")  # the order a genome's annotation is taken in
INDEX_DIR = "annotation_sources"  # in the build's sources: the other sources' indexes (fetch_indexes)
# VEuPathDB's genomes (open: each one's GenBank or RefSeq accession and GFF); its files need a user's API key
VEUPATHDB_ORGANISMS = ("https://veupathdb.org/veupathdb/service/record-types/organism/searches/GenomeDataTypes/reports/"
                       "standard?reportConfig=%7B%22attributes%22%3A%5B%22genome_version%22%2C%22URLgff%22%2C"
                       "%22annotation_source%22%2C%22annotation_version%22%5D%7D")
VEUPATHDB_KEY = "VEUPATHDB_API_KEY"  # the environment variable eukore_build.py --veupathdb_key defaults to
ENSEMBL_FTP = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"
ENSEMBL_SKIP = {"prokka", "refseq", "genbank"}  # bacteria's; NCBI's own (the ncbi source)
TOGA2_URL = "https://genome.senckenberg.de/download/TOGA2"
# TOGA2's references, the first listing an assembly gives its annotation (a class's references list its assemblies)
TOGA2_REFERENCES = ("human", "mouse", "cow", "elephant", "chicken", "zebrafinch", "crow", "kittiwake", "emu",
                    "pondturtle", "riverturtle", "stickleback", "killifish", "mandarinfish", "rainbowfish", "halibut")
INDEX_COLUMNS = ("accession", "source", "build", "url")


def gff_url(job):
    """The NCBI annotation of a job's fetched genome, by default: None only where the assembly summary says it has no
    protein-coding genes (job["annotated"] False) or the genome is not NCBI's. Where that is not known (a job planned
    before 2026-11-06, a model species') it is asked for: NCBI's 404, or a GFF without coding sequences, is no
    annotation."""
    url = job.get("url") or ""
    if job.get("annotated") is False or not url.endswith("_genomic.fna.gz"):
        return None
    return url[:-len("_genomic.fna.gz")] + "_genomic.gff.gz"


def candidates(job):
    """A job's annotations to try, in SOURCES' order -> [(source, url)]: NCBI's (gff_url) and those assign() gave it
    (job["annotations"])."""
    out = [("ncbi", gff_url(job))] if gff_url(job) else []
    return out + [tuple(x) for x in job.get("annotations", ())]


def wanted(job):
    """Whether a job has an annotation to try."""
    return bool(candidates(job))


def gff_file(genomes_folder, accession, source="ncbi"):
    """Where a genome's annotation is kept beside it (VEuPathDB's GFF is served plain)."""
    return os.path.join(genomes_folder, f"{accession}.gff.gz" if source == "ncbi" else f"{accession}.{source}.gff"
                        if source == "veupathdb" else f"{accession}.{source}.gff.gz")


def found(genomes_folder, accession):
    """The annotation beside a genome (the fetch stage's) -> (path, source), or (None, None)."""
    for source in SOURCES:
        path = gff_file(genomes_folder, accession, source)
        if os.path.isfile(path):
            return path, source
    return None, None


def missing(work, job):
    """Whether a done genome with an annotation to try has not been compared with it yet."""
    from .genomes import done_dir, is_done  # noqa: PLC0415
    return wanted(job) and is_done(work, job["accession"]) and \
        not os.path.isfile(os.path.join(done_dir(work, job["accession"]), TABLE))


def ensembl_index(species_json):
    """Ensembl's species.json -> {assembly accession: (provider, genebuild, URL of its genes.gff3.gz)}: of an assembly's
    builds (ENSEMBL_SKIP's providers left out) the newest, Ensembl's own on a tie."""
    with util.open_text(species_json) as fh:
        species = json.load(fh).get("species", {})
    out = {}
    for s in species.values():
        for acc, a in s.get("assemblies", {}).items():
            best = None
            for provider, builds in a.get("genebuild_providers", {}).items():
                if provider in ENSEMBL_SKIP:
                    continue
                for date, b in builds.items():
                    path = (b.get("paths", {}).get("genebuild", {}).get("files", {}).get("annotations", {})
                            .get("genes.gff3.gz"))
                    key = (date, provider == "ensembl")
                    if path and (best is None or key > best[0]):
                        best = (key, provider, date, f"{ENSEMBL_FTP}/{path}")
            if best:
                out[acc] = best[1:]
    return out


def veupathdb_index(organisms_json):
    """VEuPathDB's organism table (VEUPATHDB_ORGANISMS) -> {assembly accession: (annotation source and version, URL
    of its GFF)}: the URL its file service serves (/a/service/raw-files/, with a user's API key), not the table's."""
    with util.open_text(organisms_json) as fh:
        records = json.load(fh).get("records", [])
    out = {}
    for r in records:
        a = r.get("attributes", {})
        acc, url = a.get("genome_version") or "", a.get("URLgff") or ""
        if acc.startswith("GC") and url:
            url = re.sub(r"^https?://([^/]+)/common/downloads/", r"https://\1/a/service/raw-files/", url)
            out[acc] = (f"{a.get('annotation_source') or '-'} {a.get('annotation_version') or ''}".strip(), url)
    return out


def toga2_rows(assemblies_tsv, listings, classes=None):
    """TOGA2's assemblies_and_species.tsv, its references' listings ({reference folder: [assembly folders]}) and its
    genomes' ({class folder: [assembly folders]}) -> [dict(TOGA2_COLUMNS)] of the assemblies it annotated: the
    reference the first of TOGA2_REFERENCES that lists it; genome the URL of its .2bit (TOGA2's own copy, for
    assemblies outside NCBI: genomes.toga2_jobs), "-" where not listed."""
    order = sorted(listings, key=lambda r: next((i for i, k in enumerate(TOGA2_REFERENCES) if f"_{k}_" in r),
                                                len(TOGA2_REFERENCES)))
    where, cls = {}, {}
    for ref in order:
        for folder in listings[ref]:
            where.setdefault(folder, ref)
    for c, folders in (classes or {}).items():
        for folder in folders:
            cls.setdefault(folder, c)
    out = []
    with util.open_text(assemblies_tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            folder, name = r.get("Directory Name", ""), r.get("Assembly name", "")
            if folder not in where or not name or not r.get("Species Taxonomy ID", "").isdigit():
                continue
            acc = r.get("NCBI accession", "")
            n50 = r.get("scaffold N50 (bp)", "")
            out.append(dict(
                accession=f"TOGA2_{name}", ncbi_accession=acc if acc.startswith("GC") else "-",
                species_taxid=r["Species Taxonomy ID"], organism=r.get("Species", "-") or "-",
                scaffold_n50=n50 if n50.isdigit() else "0", reference=where[folder],
                annotation=f"{TOGA2_URL}/TOGA2/{where[folder]}/{folder}/query_annotation.gtf.gz",
                genome=f"{TOGA2_URL}/Genomes/{cls[folder]}/{folder}/{name}.2bit" if folder in cls else "-",
                assembly_source=r.get("Source of assembly", "-") or "-"))
    return out


def toga2_index(assemblies_tsv, listings):
    """TOGA2's annotations of NCBI's assemblies (toga2_rows) -> {assembly accession: (reference folder, URL of its
    query_annotation.gtf.gz)}."""
    return {r["ncbi_accession"]: (r["reference"], r["annotation"])
            for r in toga2_rows(assemblies_tsv, listings) if r["ncbi_accession"] != "-"}


TOGA2_COLUMNS = ("accession", "ncbi_accession", "species_taxid", "organism", "scaffold_n50", "reference", "annotation",
                 "genome", "assembly_source")


def toga2_file(src_dir):
    return os.path.join(src_dir, INDEX_DIR, "toga2_genomes.tsv")


def read_toga2(src_dir):
    """TOGA2's assemblies (fetch_indexes' toga2_genomes.tsv) -> [dict(TOGA2_COLUMNS)], [] without."""
    path = toga2_file(src_dir)
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def index_file(src_dir):
    return os.path.join(src_dir, INDEX_DIR, "annotations.tsv")


def fetch_indexes(src_dir, download, log=print):
    """Ensembl's and TOGA2's indexes into src_dir/INDEX_DIR/annotations.tsv (INDEX_COLUMNS: accession, source, build,
    url), kept where there."""
    path = index_file(src_dir)
    if os.path.isfile(path) and os.path.isfile(toga2_file(src_dir)):
        return path
    folder = os.path.dirname(path)
    os.makedirs(folder, exist_ok=True)
    rows = []
    species = os.path.join(folder, "ensembl_species.json")
    download(f"{ENSEMBL_FTP}/species.json", species)
    rows += [(acc, "ensembl", f"{provider}/{date}", url)
             for acc, (provider, date, url) in sorted(ensembl_index(species).items())]
    os.remove(species)
    organisms = os.path.join(folder, "veupathdb_organisms.json")
    download(VEUPATHDB_ORGANISMS, organisms)
    rows += [(acc, "veupathdb", build.replace("\t", " "), url)
             for acc, (build, url) in sorted(veupathdb_index(organisms).items())]
    os.remove(organisms)
    table = os.path.join(folder, "toga2_assemblies.tsv")
    download(f"{TOGA2_URL}/assemblies_and_species.tsv", table)
    listing = os.path.join(folder, "toga2.html")
    listings, classes = {}, {}
    for kind, into in (("TOGA2", listings), ("Genomes", classes)):
        download(f"{TOGA2_URL}/{kind}/", listing)
        for sub in folders_of(listing):
            if kind == "Genomes" or sub.startswith("reference_"):
                download(f"{TOGA2_URL}/{kind}/{sub}/", listing)
                into[sub] = folders_of(listing)
    os.remove(listing)
    toga2 = toga2_rows(table, listings, classes)
    os.remove(table)
    rows += [(r["ncbi_accession"], "toga2", r["reference"], r["annotation"])
             for r in sorted(toga2, key=lambda r: r["ncbi_accession"]) if r["ncbi_accession"] != "-"]
    with open(util.part_of(toga2_file(src_dir)), "w", newline="\n") as fh:
        fh.write("\t".join(TOGA2_COLUMNS) + "\n")
        fh.writelines("\t".join(r[c] for c in TOGA2_COLUMNS) + "\n" for r in toga2)
    os.replace(util.part_of(toga2_file(src_dir)), toga2_file(src_dir))
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(INDEX_COLUMNS) + "\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    os.replace(util.part_of(path), path)
    log(f"annotations: Ensembl's of {sum(r[1] == 'ensembl' for r in rows):,} assemblies, VEuPathDB's of "
        f"{sum(r[1] == 'veupathdb' for r in rows):,}, TOGA2's of {sum(r[1] == 'toga2' for r in rows):,}")
    return path


def folders_of(listing):
    """An HTML directory listing -> its subfolders' names."""
    with open(listing) as fh:
        return sorted(set(re.findall(r'href="([^"/?]+)/"', fh.read())))


def read_index(src_dir):
    """The other sources' annotations (fetch_indexes) -> {accession: [(source, url)] in SOURCES' order}, {} without."""
    path = index_file(src_dir)
    out = collections.defaultdict(list)
    if os.path.isfile(path):
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                out[r["accession"]].append((r["source"], r["url"]))
    for got in out.values():
        got.sort(key=lambda x: SOURCES.index(x[0]))
    return dict(out)


def assign(job_map, index, veupathdb=False):
    """Each job's other sources' annotations of its own assembly (its GenBank accession, or the RefSeq copy it fetches)
    into job["annotations"] -> the jobs given one. veupathdb: whether there is a VEuPathDB API key (else its
    annotations are left out: its files need one)."""
    n = 0
    for job in job_map.values():
        got = index.get(job["accession"]) or index.get(job.get("fetch", "")) or []
        got = [x for x in got if veupathdb or x[0] != "veupathdb"]
        if got:
            job["annotations"] = got
            n += 1
    return n


def calls_of(folder):
    """A done genome's calls (loci.tsv, features.tsv) -> {locus: dict(status, contig, strand, exons [(start, end)] on
    the contig, 1-based inclusive, in order)}."""
    def rows(name):
        path = os.path.join(folder, name)
        if not os.path.isfile(path):
            return []
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    on_copy = collections.defaultdict(list)
    for f in rows("features.tsv"):
        if f["kind"] == "exon":
            on_copy[f["locus"]].append((int(f["start"]), int(f["end"])))
    out = {}
    for r in rows("loci.tsv"):
        if r["status"] == "Missing" or not r["start"].isdigit():
            continue
        start, end, plus = int(r["start"]), int(r["end"]), r["strand"] == "+"
        # features are on the copy (0-based, end excluded, the copy's strand): back onto the contig
        exons = sorted((start + s, start + e - 1) if plus else (end - e + 1, end - s) for s, e in on_copy[r["locus"]])
        out[r["locus"]] = dict(status=r["status"], contig=r["contig"], strand=r["strand"], exons=exons, start=start,
                               end=end)
    return out


def report_url(job):
    """The NCBI assembly report beside a job's genome: its sequences' GenBank and RefSeq names, and the chromosomes'."""
    url = job.get("url") or ""
    return url[:-len("_genomic.fna.gz")] + "_assembly_report.txt" if url.endswith("_genomic.fna.gz") else None


def report_file(genomes_folder, accession):
    return os.path.join(genomes_folder, f"{accession}.assembly_report.txt")


def aliases(report, gff=None):
    """An NCBI assembly report -> {another name of a sequence: a name the genome may give it}: its RefSeq name's
    GenBank name and the reverse (the GFF names the RefSeq copy's sequences, a genome fetched from GenBank the GenBank
    ones), a chromosome's name (Ensembl's: `1`, `X`) its GenBank name, and an accession without its version (TOGA2's)
    the accession. gff: the annotation, whose ##sequence-region lines name a sequence the report does not (VEuPathDB's
    own names: `Pf3D7_01_v3`) its GenBank name where its length is that of one sequence alone in both."""
    out = {}
    if not report or not os.path.isfile(report):
        return out
    by_length = collections.defaultdict(list)
    with open(report) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            gb, rs = (x if x not in ("na", "") else None for x in (f[4], f[6]))
            if gb and rs:
                out[rs], out[gb] = gb, rs
            for acc in (gb, rs):
                if acc and "." in acc:
                    out.setdefault(acc.rsplit(".", 1)[0], acc)
            if gb and f[1] == "assembled-molecule" and f[2] not in ("na", ""):
                out.setdefault(f[2], gb)
            if (gb or rs) and len(f) > 8 and f[8].isdigit():
                by_length[int(f[8])].append(gb or rs)
    if gff and os.path.isfile(gff):
        theirs = collections.defaultdict(list)
        with util.open_text(gff) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                f = line.split()
                if f[0] == "##sequence-region" and len(f) == 4 and f[3].isdigit():
                    theirs[int(f[3])].append(f[1])
        known = set(out) | {n for names in by_length.values() for n in names}
        for length, names in theirs.items():
            if len(names) == 1 and len(by_length.get(length, ())) == 1 and names[0] not in known:
                out[names[0]] = by_length[length][0]
    return out


def attributes(text):
    """A GFF3 line's attributes (key=value;...) or a GTF line's (key "value"; ...) -> dict."""
    if "=" in text.split(";", 1)[0]:
        return dict(kv.split("=", 1) for kv in text.split(";") if "=" in kv)
    out = {}
    for kv in text.split(";"):
        k, _s, v = kv.strip().partition(" ")
        if k:
            out.setdefault(k, v.strip().strip('"'))
    return out


def read_cds(gff, spans, alias=None):
    """The GFF's (or GTF's) CDS features overlapping spans ({contig: [(start, end)]}) -> {transcript: dict(contig,
    strand, gene, cds [(start, end)] sorted, phase {(start, end): its phase})}; a transcript is the CDS's Parent (a GTF's
    transcript_id). alias: {annotation's sequence name: the calls'} where the two name the sequences differently
    (aliases(); two steps: a chromosome's name to its GenBank name to its RefSeq name)."""
    out = {}
    alias = alias or {}
    with util.open_text(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            if f[0] not in spans:
                name = alias.get(f[0])
                f[0] = name if name in spans else alias.get(name, f[0])
            if f[0] not in spans:
                continue
            a, b = int(f[3]), int(f[4])
            if not any(a <= e and b >= s for s, e in spans[f[0]]):
                continue
            attrs = attributes(f[8])
            parent = attrs.get("Parent") or attrs.get("transcript_id") or attrs.get("ID") or f"{f[0]}:{a}"
            t = out.setdefault(parent, dict(contig=f[0], strand=f[6], gene=urllib.parse.unquote(
                attrs.get("gene") or attrs.get("gene_name") or attrs.get("locus_tag") or attrs.get("gene_id") or "-"),
                cds=[]))
            t["cds"].append((a, b))
            t.setdefault("phase", {})[(a, b)] = int(f[7]) if f[7].isdigit() else 0
    for t in out.values():
        t["cds"] = sorted(set(t["cds"]))
    return out


def has_cds(gff):
    """Whether a GFF holds a coding sequence at all (NCBI's GFF of an assembly without gene models: its regions)."""
    with util.open_text(gff) as fh:
        return any(len(f) > 2 and f[2] == "CDS" for f in (line.split("\t", 3) for line in fh if line[0] != "#"))


def introns(exons):
    return {(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:])}


def shared_bases(x, y):
    return sum(max(0, min(b, d) - max(a, c) + 1) for a, b in x for c, d in y)


def best_transcript(call, transcripts):
    """The annotated transcript on a call's contig and strand sharing the most coding bases with its exons (on a tie,
    the one of more exons), or None where none shares any."""
    best, shared = None, 0
    for name, t in transcripts.items():
        if t["contig"] != call["contig"] or t["strand"] != call["strand"]:
            continue
        s = shared_bases(call["exons"], t["cds"])
        if s > shared or (s == shared and s and best and len(t["cds"]) > len(transcripts[best]["cds"])):
            best, shared = name, s
    return best


def compare(call, transcripts):
    """One call against the annotated transcripts on its contig and strand -> ANNOTATION_COLUMNS' row (no locus)."""
    exons = call["exons"]
    best = best_transcript(call, transcripts)
    shared = shared_bases(exons, transcripts[best]["cds"]) if best else 0
    row = dict(status=call["status"], contig=call["contig"], exons=len(exons), introns=max(0, len(exons) - 1),
               bases=sum(b - a + 1 for a, b in exons))
    if best is None:
        return dict(row, transcript="-", gene="-", ann_exons=0, exact_exons=0, ann_introns=0, exact_introns=0,
                    ann_bases=0, shared_bases=0)
    ann = transcripts[best]["cds"]
    return dict(row, transcript=best, gene=transcripts[best]["gene"], ann_exons=len(ann),
                exact_exons=len(set(exons) & set(ann)), ann_introns=max(0, len(ann) - 1),
                exact_introns=len(introns(exons) & introns(ann)), ann_bases=sum(b - a + 1 for a, b in ann),
                shared_bases=shared)


def annotate(folder, gff, report=None, calls=None, source="ncbi"):
    """A genome's calls (calls_of(folder) where not given: a done genome's, folder WORK/done/<accession>) against its
    annotation (gff, or None: none; report: NCBI's assembly report, the sequences' names; source: SOURCES') into
    folder/TABLE -> the copies compared."""
    calls = calls_of(folder) if calls is None else calls
    rows = []
    if gff and calls:
        spans = collections.defaultdict(list)
        for c in calls.values():
            spans[c["contig"]].append((c["start"], c["end"]))
        transcripts = read_cds(gff, spans, aliases(report, gff))
        if transcripts or has_cds(gff):  # a GFF without coding sequences (an assembly's regions alone): no annotation
            rows = [dict(compare(c, transcripts), locus=locus, source=source) for locus, c in sorted(calls.items())]
    path = os.path.join(folder, TABLE)
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(ANNOTATION_COLUMNS) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in ANNOTATION_COLUMNS) + "\n" for r in rows)
    os.replace(util.part_of(path), path)
    return len(rows)


def fetch(job, genomes_folder, veupathdb_key=None, **how):
    """A job's annotation into genomes_folder (kept where there), from the first of its candidates() that has one
    (a 404, or a GFF without coding sequences, is none), with NCBI's assembly report -> (path, source), (None, None)
    where none has one, or raises OSError where one did not arrive (a later run asks again). veupathdb_key: a
    VEuPathDB user's API key (without, its annotation is not asked for); how: util.fetch_file's scratch and checker."""
    path, source = found(genomes_folder, job["accession"])
    if path:
        return path, source
    report = report_file(genomes_folder, job["accession"])
    if report_url(job) and not os.path.isfile(report):
        util.fetch_file(report_url(job), report, gzip_check=False,  # without it, the names may not meet
                        **{k: v for k, v in how.items() if k == "scratch"})
    for source, url in candidates(job):
        path = gff_file(genomes_folder, job["accession"], source)
        if source == "veupathdb":
            if not veupathdb_key:
                continue
            why = util.fetch_file(url, path, gzip_check=False, headers={"Authorization": f"Bearer {veupathdb_key}"},
                                  **{k: v for k, v in how.items() if k == "scratch"})
        else:
            why = util.fetch_file(url, path, **how)
        if why == "not found":
            continue
        if why:
            raise OSError(f"{url}: {why}")
        if has_cds(path):
            return path, source
        os.remove(path)  # NCBI's GFF of an unannotated assembly: its regions alone
    return None, None


def annotate_genome(work, job, genomes_folder, may_fetch=True, veupathdb_key=None):
    """A done genome compared with its annotation (fetched where may_fetch and not beside the genome), the GFF
    deleted after -> copies compared, or None where its GFF is not there (left for a stage that may fetch)."""
    from .genomes import done_dir  # noqa: PLC0415
    path, source = found(genomes_folder, job["accession"])
    if not path:
        if not may_fetch:
            return None
        path, source = fetch(job, genomes_folder, veupathdb_key)
    report = report_file(genomes_folder, job["accession"])
    n = annotate(done_dir(work, job["accession"]), path, report, source=source or "-")
    for done in (path, report):
        if done and os.path.isfile(done):
            os.remove(done)
    return n
