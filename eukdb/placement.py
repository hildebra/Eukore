"""Placing euk30 copies in a genome with miniprot (since 2026-11-05; copied from protal's
scripts/mini_db/eukdetect1_genomic.py, whose per-genome run is genomes.py's here): a protein's best alignment in the
genome is its copy, the genome from the alignment's first codon to its last, introns and all, complete (EUKDB keeps
every copy whole; protal fits it into 16 kb when it converts the package), with its exons and introns (copies.fit's
features). A protein is the copy's own coding sequence (EukDetect1's), judged by MIN_COVERAGE, MIN_IDENTITY and
MIN_KMER_SHARE, or a relative's (species.best, species.taken: the model and RefSeq species).
"""

import collections
import gzip
import os
import re
import subprocess
import tempfile

from . import copies as euk30_copies
from . import sources as eukaryotes
from . import species as model_species

UNCAPPED = 10 ** 12  # a copy's length kept whole


MIN_COVERAGE = 0.95


MIN_IDENTITY = 0.9


MIN_KMER_SHARE = 0.5


MASK_MIN = 100


PLACED = ("genomic", "poor placement")


K = 21


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


CODE = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(cds):
    """A coding sequence's protein for miniprot: its stop codon dropped, an internal stop or unknown codon as X."""
    protein = "".join(CODE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))
    return (protein[:-1] if protein.endswith("*") else protein).replace("*", "X")


def kmers(seq, k=K):
    return {seq[i:i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i + k]}


def candidates(eukdetect1, eukdetect, locus_map, ncbi, groups=None):
    """EukDetect1's genomes whose copies build_release adds: {genome: (taxid, EukDetect1 group, {record: coding
    sequence})}, the genomes of species EukDetect2's link names (by their NCBI species) left out."""
    genomes1, records1 = eukaryotes.read_link1(os.path.join(eukdetect1, "busco_taxid_link.txt"), locus_map)
    copies1 = eukaryotes.read_loci(os.path.join(eukdetect1, "ncbi_eukprot_met_arch_markers.fna"), records1)
    genomes2, _records2 = eukaryotes.read_link(os.path.join(eukdetect, "busco_taxid_genome_link.txt"))
    known = {ncbi.species_of(taxid) for taxid, _group in genomes2.values()}
    group_of = {}  # EukDetect1's own group (worm, insect, ...): read_link1 keeps the broader one
    record_of = {}
    for record, (genome, locus) in records1.items():
        group_of[genome] = record.split("-", 1)[0]
        record_of[(genome, locus)] = record
    result = {}
    for genome, by_locus in copies1.items():
        taxid, _group = genomes1[genome]
        sp = ncbi.species_of(taxid)
        if sp in known or genome in genomes2 or ncbi.rank.get(sp) != "species" or eukaryotes.is_placeholder(ncbi, sp):
            continue
        if len(by_locus) < eukaryotes.MIN_EUKDETECT1_LOCI or (groups and group_of[genome] not in groups):
            continue
        result[genome] = (taxid, group_of[genome], {record_of[(genome, locus)]: seq for locus, seq in by_locus.items()})
    return result


def model_candidates(table, eukdetect1, eukdetect, locus_map, ncbi, loci, placed=None):
    """--models (since 2026-11-05): the eukaryotes of a model species table (model_species.tsv, or the download's
    table of the RefSeq species euk30 lacks) whose euk30 sources (EukDetect2's genomes of the species, EukDetect1's
    single copies of the loci locus_map maps) lack loci -> {genome name: (taxid, GROUP, {query record: coding
    sequence}, assembly row)}: each lacking locus's proteins of the species' nearest relatives (euk30_exons.Relatives:
    EukDetect1's coding sequences, EukDetect2's single-frame copies, the placed model species' coding sequences,
    placed: their copies' FASTA files), the species' reference assembly (model_species.assembly_row)."""
    rows = model_species.eukaryotes(model_species.read_table(table))
    loci = set(loci)
    have, (genomes2, records2) = eukaryotes.species_loci(eukdetect, eukdetect1, locus_map, ncbi, loci)
    wanted = {}
    for row in rows:
        lacking = loci - have[ncbi.species_of(int(row["ncbi_taxid"]))]
        if lacking:
            wanted[model_species.genome_name(row)] = (row, lacking)
    if not wanted:
        return {}
    from . import exons as euk30_exons  # noqa: PLC0415  (Relatives: each locus's copies by ancestor taxon)
    relatives, _link, _copies = euk30_exons.relatives_of(ncbi, eukdetect, eukdetect1, loci, locus_map, placed or (),
                                                         link2=(genomes2, records2))
    return {name: (int(row["ncbi_taxid"]), model_species.GROUP,
                   relatives.queries(name, int(row["ncbi_taxid"]), dict.fromkeys(lacking),
                                     model_species.record(row, "x").rsplit("-", 2)[0]),
                   model_species.assembly_row(row)) for name, (row, lacking) in wanted.items()}


def genome_url(row):
    """An assembly's genome on NCBI's FTP server (its ftp_path's folder and file)."""
    ftp = row["ftp"].rstrip("/")
    if not ftp or ftp == "na":
        return None
    return f"{ftp.replace('ftp://', 'https://')}/{os.path.basename(ftp)}_genomic.fna.gz"


def parse_gff(text):
    """miniprot's GFF -> {protein: (contig, start, end, strand, identity, (query start, query end), [(start, end)])}
    of each protein's first (best) alignment; positions 1-based inclusive, the exons in genome order."""
    mrna, exons = {}, collections.defaultdict(list)
    for line in text.splitlines():
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        attrs = dict(a.split("=", 1) for a in f[8].split(";") if "=" in a)
        if f[2] == "mRNA":
            target = attrs["Target"].split()
            if target[0] not in {v[0] for v in mrna.values()}:
                mrna[attrs["ID"]] = (target[0], f[0], int(f[3]), int(f[4]), f[6], float(attrs.get("Identity", 0)),
                                     (int(target[1]), int(target[2])))
        elif f[2] == "CDS":
            exons[attrs["Parent"]].append((int(f[3]), int(f[4])))
    return {protein: (contig, start, end, strand, identity, target, sorted(exons[mid]))
            for mid, (protein, contig, start, end, strand, identity, target) in mrna.items()}


def read_contigs(path, wanted):
    """The sequences (case kept) of the contigs named in `wanted` of a FASTA file."""
    out, name, parts = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    out[name] = "".join(parts)
                name, parts = (line[1:].split() or [""])[0], []
            elif name in wanted:
                parts.append(line.strip())
    if name in wanted:
        out[name] = "".join(parts)
    return out


def status(identity, coverage, share):
    """A placed copy's status: taken (genomic) or not."""
    return PLACED[0] if coverage >= MIN_COVERAGE and identity >= MIN_IDENTITY and share >= MIN_KMER_SHARE else PLACED[1]


def final(text):
    """Whether a copy's status is a final outcome (since 2026-11-03): placed (taken or not), not placed, no assembly, a
    genome over --max_genome_mb, or a fetch the server refused for good (not found, HTTP 4xx but 429:
    rdna_store.given_up's rule). Not final, and placed again by a later run: a fetch a retry may mend (a timeout, an
    HTTP 503 of a burst at -j 16, five tries in ~30 s), a miniprot run that failed (killed for memory), an rRNA
    extraction that failed. Before, those were written as done, merged into the checkout's placed table, and their
    genomes never placed again."""
    if text.startswith("fetch failed: "):
        return bool(given_up({"status": "failed: " + text[len("fetch failed: "):]}))
    return not text.startswith(("miniprot failed", "rRNA extraction failed"))


CHUNK_BASES = 4_000_000_000  # a genome over this many bases (gzipped: its file size / 0.3) is placed in pieces of it


GZIP_BYTES_PER_BASE = 0.3


def chunks(genome_path, folder, max_bases=None):
    """A genome as miniprot takes it: itself, or where it is over max_bases (CHUNK_BASES; since 2026-11-05: a 3 Gb
    genome's index takes ~22 GB, a lungfish's 40 Gb ~300) its contigs in pieces of max_bases at most, plain FASTA in
    folder (a contig longer than that alone in its own) -> [paths]."""
    max_bases = max_bases or CHUNK_BASES
    size = os.path.getsize(genome_path) / (GZIP_BYTES_PER_BASE if genome_path.endswith(".gz") else 1)
    if size <= max_bases:
        return [genome_path]
    paths, out, held = [], None, 0
    with (gzip.open(genome_path, "rt") if genome_path.endswith(".gz") else open(genome_path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if out is None or held >= max_bases:
                    if out:
                        out.close()
                    paths.append(os.path.join(folder, f"piece{len(paths)}.fna"))
                    out, held = open(paths[-1], "w"), 0
            else:
                held += len(line) - 1
            out.write(line)
    if out:
        out.close()
    return paths


def place(cds_by_record, genome_path, miniprot="miniprot", threads=1, mask_min=MASK_MIN, max_span=None, model=False):
    """Each coding sequence's genomic copy in a genome. -> {record: (copy or None, table fields from status on,
    features or None)}, the copy of every placed record, taken or not. Since 2026-11-05 a copy is fitted into max_span
    bases (default eukaryotes.MAX_LOCUS; euk30_copies.fit: every exon kept, the introns cut from their middles), its
    features (euk30_copies.fit's) where its exons and introns are and what was cut. A model species' copy (model: a
    relative's protein, none of whose 21-mers is the genome's) has as its kmer_share field the share of the alignment's
    coding bases the copy holds (1 but where its exons alone exceed max_span)."""
    hits = {}
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(genome_path))) as tmp:
        faa = os.path.join(tmp, "copies.faa")
        with open(faa, "w") as fh:
            for record, cds in cds_by_record.items():
                fh.write(f">{record}\n{translate(cds)}\n")
        for part in chunks(genome_path, tmp):  # a genome over CHUNK_BASES in pieces: miniprot's index of ~7.4 B a base
            run = subprocess.run([miniprot, "-t", str(threads), "-I", "--gff", "--outn=1", part, faa],
                                 capture_output=True, text=True, check=False)
            if run.returncode != 0:
                return {r: (None, ["miniprot failed: " + (run.stderr.strip().splitlines() or ["?"])[-1][:200]], None)
                        for r in cds_by_record}
            for protein, hit in parse_gff(run.stdout).items():  # each protein's best over the pieces
                score = hit[4] * (hit[5][1] - hit[5][0] + 1)
                if protein not in hits or score > hits[protein][4] * (hits[protein][5][1] - hits[protein][5][0] + 1):
                    hits[protein] = hit
    contigs = read_contigs(genome_path, {h[0] for h in hits.values()})
    result = {}
    for record, cds in cds_by_record.items():
        hit = hits.get(record)
        if hit is None:
            result[record] = (None, ["not placed"], None)
            continue
        contig, start, end, strand, identity, (q1, q2), exons = hit
        protein_length = len(translate(cds))
        coverage = (q2 - q1 + 1) / max(protein_length, 1)
        fitted = euk30_copies.fit(contigs[contig], start, end, strand, exons, max_span or UNCAPPED,
                                  mask_min=mask_min)
        copy, spliced = fitted["copy"], fitted["spliced"]
        if model:
            coding = sum(b - a + 1 for a, b in exons)
            share = len(spliced) / max(coding, 1)
        else:
            share = len(kmers(cds) & kmers(spliced)) / max(len(kmers(cds)), 1)
        result[record] = (copy, [status(identity, coverage, share), contig, start, end, strand, len(exons), len(cds),
                                 len(copy), f"{identity:.4f}", f"{coverage:.3f}", f"{share:.3f}", fitted["masked"]],
                          feature_rows_of(record, contig, strand, fitted))
    return result


def feature_rows_of(record, contig, strand, fitted):
    """A placed copy's features as euk30_copies.FEATURES rows, its record's locus (EukDetect1's for its records, euk30's
    for a model species') and no accession (the genome's: do_genome sets it)."""
    m = eukaryotes.RECORD_RE.match(record.split("~", 1)[0])
    return euk30_copies.feature_rows(record.split("~", 1)[0], "", m.group("locus") if m else "", contig, strand, fitted)


def given_up(row):
    """Why a genome is not asked for again (the server has it not, or refuses: not found, HTTP 4xx but 429), from its
    status row ({"status": "failed: <why>"}); None where a retry may mend it (protal's rdna_store.given_up rule)."""
    status = (row or {}).get("status", "")
    why = status[len("failed: "):] if status.startswith("failed: ") else None
    if why and (why == "not found" or (re.fullmatch(r"HTTP 4\d\d", why) and why != "HTTP 429")):
        return why
    return None
