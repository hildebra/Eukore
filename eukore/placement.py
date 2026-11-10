"""A genome as the locus caller reads it (since 2026-11-05; the rest of protal's eukdetect1_genomic.py, whose
EukDetect1 placements loci.py replaced): its URL on NCBI's server, its contigs, its pieces for miniprot, miniprot's
saved indexes of them (since 2026-11-10: saved_indexes), and which failed fetches are final.
"""

import bisect
import contextlib
import gzip
import hashlib
import os
import re
import tempfile
import threading


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


def contig_lengths(path):
    """A FASTA's sequences (gzipped or plain) -> [(name, bases)] in its order."""
    out, name, n = [], None, 0
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, n))
                name, n = (line[1:].split() or [""])[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        out.append((name, n))
    return out


def chunks(genome_path, folder, max_bases=None):
    """A genome as miniprot takes it: itself, or where it is over max_bases (CHUNK_BASES: a 3 Gb genome's index takes
    ~22 GB, a lungfish's 40 Gb ~300) its contigs in order in pieces of max_bases at most, plain FASTA in folder (a
    contig longer than that alone in its own) -> [paths]. A piece is closed before a contig that would take it past
    max_bases (its length from a first pass; until 2026-11-08 a piece took a whole further contig past it: chromosomes
    of 3 Gb made pieces of 6, ~45 GB of index where the build reserved ~33)."""
    max_bases = max_bases or CHUNK_BASES
    if genome_bases(genome_path) <= max_bases:
        return [genome_path]
    lengths = iter(contig_lengths(genome_path))
    paths, out, held = [], None, 0
    with (gzip.open(genome_path, "rt") if genome_path.endswith(".gz") else open(genome_path)) as fh:
        for line in fh:
            if line.startswith(">"):
                n = next(lengths)[1]
                if out is None or (held and held + n > max_bases):
                    if out:
                        out.close()
                    paths.append(os.path.join(folder, f"piece{len(paths)}.fna"))
                    out, held = open(paths[-1], "w"), 0
                held += n
            if out is not None:  # (bases before any title: none of a contig)
                out.write(line)
    if out:
        out.close()
    return paths


INDEX_SUFFIX = ".mpi"
_INDEX_LOCKS, _INDEX_GUARD = {}, threading.Lock()


def index_file(index_dir, genome_path, key, piece=None, max_bases=None):
    """The file in index_dir of miniprot's saved index of a genome (piece None) or of one of its pieces (its number in
    chunks()' order, of max_bases at most). Its name holds all that shapes the index: the genome (its path, size and
    time: another file under that name is another genome), the piece, and key (the tool and its indexing options:
    loci.split_options')."""
    stat = os.stat(genome_path)
    what = [os.path.abspath(genome_path), stat.st_size, stat.st_mtime_ns, *key]
    stem = os.path.basename(genome_path)
    if piece is not None:
        what += ["piece", max_bases or CHUNK_BASES, piece]
        stem += f".piece{piece}"
    digest = hashlib.sha1("\n".join(str(x) for x in what).encode()).hexdigest()[:16]
    return os.path.join(index_dir, f"{stem}.{digest}{INDEX_SUFFIX}")


def saved_indexes(genome_path, folder, index_dir, key, build, max_bases=None):
    """A genome as miniprot takes it where it is aligned more than once (since 2026-11-10: a dataset species' genome is
    called up to five times a run, the main call and its hold-outs, and each call made miniprot's index, most of a
    large genome's time, again): the saved indexes of its pieces (chunks()) in index_dir, each made once (index_file
    names it) -> ([their paths] in the pieces' order, why one could not be made or None).
    build(source, dest): writes the index of a FASTA to dest -> why it failed, or None. It writes a name of its own
    that is renamed once it is whole, so a failed or killed build leaves no file a later call takes for an index; one
    build at a time a genome in this process (threads); two processes at worst build twice, neither reads half a
    file. The pieces are written into folder and each deleted once its index is there; their number is noted beside
    the indexes, so a later call does not cut the genome again. Nothing is kept outside index_dir, which is the
    caller's to delete (~7.5 bytes a base)."""
    max_bases = max_bases or CHUNK_BASES
    os.makedirs(index_dir, exist_ok=True)

    def names(n):
        return [index_file(index_dir, genome_path, key, i if n > 1 else None, max_bases) for i in range(n)]

    noted = index_file(index_dir, genome_path, ["pieces", max_bases])[:-len(INDEX_SUFFIX)] + ".pieces"
    with _INDEX_GUARD:
        lock = _INDEX_LOCKS.setdefault((os.path.abspath(index_dir), os.path.abspath(genome_path)), threading.Lock())
    with lock:
        try:
            with open(noted) as fh:
                paths = names(int(fh.read()))
        except (OSError, ValueError):
            paths = []
        if paths and all(os.path.isfile(p) for p in paths):
            return paths, None
        parts = chunks(genome_path, folder, max_bases)
        paths = names(len(parts))
        for part, dest in zip(parts, paths):
            if not os.path.isfile(dest):
                fd, partial = tempfile.mkstemp(prefix=os.path.basename(dest) + ".", suffix=".partial", dir=index_dir)
                os.close(fd)
                try:
                    why = build(part, partial)
                    if why:
                        return [], why
                    os.replace(partial, dest)
                finally:
                    with contextlib.suppress(OSError):
                        os.remove(partial)
            if part != genome_path:
                os.remove(part)  # read from its index from here on
        fd, partial = tempfile.mkstemp(prefix=os.path.basename(noted) + ".", suffix=".partial", dir=index_dir)
        with os.fdopen(fd, "w") as fh:
            fh.write(f"{len(parts)}\n")
        os.replace(partial, noted)
    return paths, None


def given_up(row):
    """Why a genome is not asked for again (the server has it not, or refuses: not found, HTTP 4xx but 429), from its
    status row ({"status": "failed: <why>"}); None where a retry may mend it."""
    status = (row or {}).get("status", "")
    why = status[len("failed: "):] if status.startswith("failed: ") else None
    if why and (why == "not found" or (re.fullmatch(r"HTTP 4\d\d", why) and why != "HTTP 429")):
        return why
    return None
