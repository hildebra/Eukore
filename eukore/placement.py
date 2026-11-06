"""A genome as the locus caller reads it (since 2026-11-05; the rest of protal's eukdetect1_genomic.py, whose
EukDetect1 placements loci.py replaced): its URL on NCBI's server, its contigs, its pieces for miniprot, and which
failed fetches are final.
"""

import gzip
import os
import re


def genome_url(row):
    """An assembly's genome on NCBI's FTP server (its ftp_path's folder and file)."""
    ftp = (row.get("ftp") or "").rstrip("/")
    if not ftp or ftp == "na":
        return None
    return f"{ftp.replace('ftp://', 'https://')}/{os.path.basename(ftp)}_genomic.fna.gz"


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


CHUNK_BASES = 4_000_000_000  # a genome over this many bases (gzipped: its file size / 0.3) is placed in pieces of it
GZIP_BYTES_PER_BASE = 0.3


def chunks(genome_path, folder, max_bases=None):
    """A genome as miniprot takes it: itself, or where it is over max_bases (CHUNK_BASES: a 3 Gb genome's index takes
    ~22 GB, a lungfish's 40 Gb ~300) its contigs in pieces of max_bases at most, plain FASTA in folder (a contig longer
    than that alone in its own) -> [paths]."""
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


def given_up(row):
    """Why a genome is not asked for again (the server has it not, or refuses: not found, HTTP 4xx but 429), from its
    status row ({"status": "failed: <why>"}); None where a retry may mend it."""
    status = (row or {}).get("status", "")
    why = status[len("failed: "):] if status.startswith("failed: ") else None
    if why and (why == "not found" or (re.fullmatch(r"HTTP 4\d\d", why) and why != "HTTP 429")):
        return why
    return None
