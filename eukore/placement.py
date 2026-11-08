"""A genome as the locus caller reads it (since 2026-11-05; the rest of protal's eukdetect1_genomic.py, whose
EukDetect1 placements loci.py replaced): its URL on NCBI's server, its contigs, its pieces for miniprot, and which
failed fetches are final.
"""

import bisect
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


class Spans:
    """A contig's stretches as read_spans read them, sliced as the contig itself (seq[a:b], 0-based, end excluded) where
    a stretch holds the slice; len() is the contig's whole length. A slice outside what was read is an IndexError."""

    def __init__(self, length, pieces):
        self.length, self.pieces = length, pieces  # pieces: [(0-based start, text)] in order
        self.starts = [s for s, _t in pieces]

    def __len__(self):
        return self.length

    def __getitem__(self, key):
        if not isinstance(key, slice) or key.step not in (None, 1):
            raise TypeError("a Spans is sliced [a:b] alone")
        a = 0 if key.start is None else key.start
        b = self.length if key.stop is None else min(key.stop, self.length)
        if a < 0 or b < 0:
            raise IndexError("a Spans takes no negative index")
        if b <= a:
            return ""
        i = bisect.bisect_right(self.starts, a) - 1
        if i >= 0 and b <= self.starts[i] + len(self.pieces[i][1]):
            s, text = self.pieces[i]
            return text[a - s:b - s]
        raise IndexError(f"{a}-{b} was not read")


def merged(intervals):
    """[(start, end)] 1-based inclusive -> the same merged and sorted, starts at 1 at least."""
    out = []
    for a, b in sorted((max(1, int(a)), int(b)) for a, b in intervals):
        if b < a:
            continue
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


def read_spans(path, wanted):
    """The stretches wanted ({contig: [(start, end)]}, 1-based inclusive; merged, clipped to the contig) of a FASTA file,
    case kept -> {contig: Spans}. Only these stretches are held (since 2026-11-08: read_contigs held every contig with a
    call, a chromosome-level assembly whole, ~1 byte a base)."""
    spans = {c: merged(v) for c, v in wanted.items() if v}
    out, name, held, pos, k = {}, None, None, 0, 0

    def close():
        if held is not None:
            out[name] = Spans(pos, [(a - 1, "".join(parts)) for (a, _b), parts in zip(spans[name], held)])

    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                close()
                name = (line[1:].split() or [""])[0]
                held = [[] for _ in spans[name]] if name in spans else None
                pos = k = 0
            elif held is not None:
                s = line.strip()
                n, mine = len(s), spans[name]
                while k < len(mine) and mine[k][1] <= pos:  # stretches ended before this line
                    k += 1
                j = k
                while j < len(mine) and mine[j][0] - 1 < pos + n:
                    lo, hi = max(mine[j][0] - 1, pos), min(mine[j][1], pos + n)
                    if hi > lo:
                        held[j].append(s[lo - pos:hi - pos])
                    j += 1
                pos += n
    close()
    return out


CHUNK_BASES = 4_000_000_000  # a genome over this many bases (gzipped: its file size / 0.3) is placed in pieces of it
GZIP_BYTES_PER_BASE = 0.3


def genome_bases(path):
    """A genome's bases, estimated from its file's size (GZIP_BYTES_PER_BASE for a gzipped one)."""
    return os.path.getsize(path) / (GZIP_BYTES_PER_BASE if path.endswith(".gz") else 1)


def chunks(genome_path, folder, max_bases=None):
    """A genome as miniprot takes it: itself, or where it is over max_bases (CHUNK_BASES: a 3 Gb genome's index takes
    ~22 GB, a lungfish's 40 Gb ~300) its contigs in pieces of max_bases at most, plain FASTA in folder (a contig longer
    than that alone in its own) -> [paths]."""
    max_bases = max_bases or CHUNK_BASES
    if genome_bases(genome_path) <= max_bases:
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
