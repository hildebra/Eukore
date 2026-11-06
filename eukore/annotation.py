"""NCBI's own gene models of a genome (since 2026-11-06): the BUSCO copies' models where they are good, and a truth set
independent of BUSCO's dataset.

A genome with an annotation at NCBI (a RefSeq copy, GCF_, mostly NCBI's Eukaryotic Genome Annotation Pipeline's; a
GenBank assembly its submitter's: the assembly summary's protein-coding genes, sources.read_assemblies' `genes`) has a
GFF beside the genome on NCBI's server. The process stage gives it to loci.call: a copy takes the annotated transcript
sharing the most coding bases with miniprot's model (best_transcript) where the locus' HMM confirms it (loci.annotated,
copies.tsv's `model` ncbi). And each copy's miniprot model is laid on the GFF's coding sequences (CDS) on its contig
(WORK/done/<accession>/annotation.tsv, ANNOTATION_COLUMNS): its exons and introns, the transcript's, those exact on
both, and the coding bases each holds and shares. The introns' exact boundaries (splice sites) are the measure least
touched by start and stop codons (miniprot's CDS and NCBI's end differently at a stop).

The GFF is fetched by the fetch stage beside the genome (OUT/genomes/<accession>.gff.gz), or by the process stage
where it may fetch, and deleted once the genome is done. A genome without an annotation (NCBI's 404) gets an
annotation.tsv without rows, so that it is not asked for again. A genome done before its GFF was there is compared
later (annotate_genome; its copies keep miniprot's models). benchmarks.py gathers the tables (benchmarks/annotation.tsv).
"""

import collections
import os
import urllib.parse

from . import util

TABLE = "annotation.tsv"
ANNOTATION_COLUMNS = ("locus", "status", "contig", "transcript", "gene", "exons", "ann_exons", "exact_exons",
                      "introns", "ann_introns", "exact_introns", "bases", "ann_bases", "shared_bases")


def gff_url(job):
    """The NCBI annotation of a job's fetched genome, by default: None only where the assembly summary says it has no
    protein-coding genes (job["annotated"] False) or the genome is not NCBI's. Where that is not known (a job planned
    before 2026-11-06, a model species') it is asked for: NCBI's 404, or a GFF without coding sequences, is no
    annotation."""
    url = job.get("url") or ""
    if job.get("annotated") is False or not url.endswith("_genomic.fna.gz"):
        return None
    return url[:-len("_genomic.fna.gz")] + "_genomic.gff.gz"


def gff_file(genomes_folder, accession):
    return os.path.join(genomes_folder, f"{accession}.gff.gz")


def missing(work, job):
    """Whether a done genome with an NCBI annotation has not been compared with it yet."""
    from .genomes import done_dir, is_done  # noqa: PLC0415
    return bool(gff_url(job)) and is_done(work, job["accession"]) and \
        not os.path.isfile(os.path.join(done_dir(work, job["accession"]), TABLE))


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
    """The NCBI assembly report beside a job's annotation: its sequences' GenBank and RefSeq names."""
    url = gff_url(job)
    return url[:-len("_genomic.gff.gz")] + "_assembly_report.txt" if url else None


def report_file(genomes_folder, accession):
    return os.path.join(genomes_folder, f"{accession}.assembly_report.txt")


def aliases(report):
    """An NCBI assembly report -> {a sequence's RefSeq name: its GenBank name, and the reverse}: the GFF names the
    RefSeq copy's sequences, a genome fetched from GenBank (a model species', an older run's) the GenBank ones."""
    out = {}
    if not report or not os.path.isfile(report):
        return out
    with open(report) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7 and f[4] not in ("na", "") and f[6] not in ("na", ""):
                out[f[6]], out[f[4]] = f[4], f[6]
    return out


def read_cds(gff, spans, alias=None):
    """The GFF's CDS features overlapping spans ({contig: [(start, end)]}) -> {transcript: dict(contig, strand, gene,
    cds [(start, end)] sorted, phase {(start, end): its phase})}; a transcript is the CDS's Parent. alias: {GFF
    sequence name: the calls'} where the two name the sequences differently."""
    out = {}
    alias = alias or {}
    with util.open_text(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[0] not in spans and alias.get(f[0]) in spans:
                f[0] = alias[f[0]]
            if len(f) < 9 or f[2] != "CDS" or f[0] not in spans:
                continue
            a, b = int(f[3]), int(f[4])
            if not any(a <= e and b >= s for s, e in spans[f[0]]):
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            parent = attrs.get("Parent") or attrs.get("ID") or f"{f[0]}:{a}"
            t = out.setdefault(parent, dict(contig=f[0], strand=f[6], gene=urllib.parse.unquote(
                attrs.get("gene") or attrs.get("locus_tag") or "-"), cds=[]))
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


def annotate(folder, gff, report=None, calls=None):
    """A genome's calls (calls_of(folder) where not given: a done genome's, folder WORK/done/<accession>) against its
    NCBI annotation (gff, or None: none; report: its assembly report, the sequences' names) into folder/TABLE -> the
    copies compared."""
    calls = calls_of(folder) if calls is None else calls
    rows = []
    if gff and calls:
        spans = collections.defaultdict(list)
        for c in calls.values():
            spans[c["contig"]].append((c["start"], c["end"]))
        transcripts = read_cds(gff, spans, aliases(report))
        if transcripts or has_cds(gff):  # a GFF without coding sequences (an assembly's regions alone): no annotation
            rows = [dict(compare(c, transcripts), locus=locus) for locus, c in sorted(calls.items())]
    path = os.path.join(folder, TABLE)
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(ANNOTATION_COLUMNS) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in ANNOTATION_COLUMNS) + "\n" for r in rows)
    os.replace(util.part_of(path), path)
    return len(rows)


def fetch(job, genomes_folder, **how):
    """A job's NCBI annotation into genomes_folder (kept where there) -> its path, None where NCBI has none, or raises
    OSError where it did not arrive (a later run asks again). how: util.fetch_file's scratch and checker."""
    url = gff_url(job)
    if not url:
        return None
    path = gff_file(genomes_folder, job["accession"])
    report = report_file(genomes_folder, job["accession"])
    if not os.path.isfile(report):
        util.fetch_file(report_url(job), report, gzip_check=False,  # without it, the names may not meet
                        scratch=how.get("scratch"))
    if os.path.isfile(path):
        return path
    why = util.fetch_file(url, path, **how)
    if why == "not found":
        return None
    if why:
        raise OSError(f"{url}: {why}")
    return path


def annotate_genome(work, job, genomes_folder, may_fetch=True):
    """A done genome compared with its NCBI annotation (fetched where may_fetch and not beside the genome), the GFF
    deleted after -> copies compared, or None where its GFF is not there (left for a stage that may fetch)."""
    from .genomes import done_dir  # noqa: PLC0415
    path = gff_file(genomes_folder, job["accession"])
    if not os.path.isfile(path):
        if not may_fetch:
            return None
        path = fetch(job, genomes_folder)
    report = report_file(genomes_folder, job["accession"])
    n = annotate(done_dir(work, job["accession"]), path, report)
    for done in (path, report):
        if done and os.path.isfile(done):
            os.remove(done)
    return n
