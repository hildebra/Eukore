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
    """A genome's bases, estimated from its file's size (GZIP_BYTES_PER_BASE for a gzipped one). Where nothing
    better is known: against the true length the estimate was 0.98 to 1.07 times it for the 2026-11 trial's 34
    genomes, and falls with the share of N (an assembly's gaps compress to nothing: the yeast's genome with 29% of
    it made N, 0.79 times; with 55%, 0.51); loci.call takes the true size where its caller has it (since
    2026-11-10)."""
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


# A contig over CHUNK_BASES cut into windows (since 2026-11-10, the audit's C8: it was alone in its piece, whole, and
# its index over the memory the build reserves, for CHUNK_BASES at most: eukore_build.run_memory; lungfish
# chromosomes are several Gb). A window reaches CHUNK_OVERLAP bases into the one before, so that a gene up to that
# long is whole in one of the two: miniprot's -I allows an intron of 300 kb at most (its manual), 5 Mb are 16 such.
CHUNK_OVERLAP = 5_000_000
WINDOW_MARK = "|eukore_window@"  # a window's name: its contig's, the mark, its first base's 0-based place
# miniprot's index holds a sequence's k-mers by blocks of 256 bases (-b 8), counted from its first base on the
# forward strand and from its last on the reverse: a window cut between two of the contig's blocks seeds a little
# otherwise far from its ends too (a made-up contig of 2.9 Mb with 48 C. elegans genes: its first 100, 1,000 or
# 5,000 bases left off, 3 to 8 of 971 alignments other than the whole contig's; its last 100, 128 or 300 off, 2 to 3,
# all on the reverse strand; 256, 512, 768, 1,024 or 65,536 off either end: none). So a window begins a multiple
# of 256 bases from the contig's start and ends one from its end.
WINDOW_GRID = 256


def windows(length, max_bases, overlap):
    """A contig's windows -> [(start, end)] 0-based, end excluded: itself where it is max_bases or shorter (or no
    overlap is given: not cut), else windows of max_bases or up to WINDOW_GRID - 1 bases less, each `overlap` bases
    or a little more into the one before, the last to the contig's end; each begins and ends on the contig's
    multiples of WINDOW_GRID from either end (where the sizes allow: not a made-up tiny one)."""
    if not overlap or length <= max_bases:
        return [(0, length)]
    size = max_bases - (max_bases - length) % WINDOW_GRID
    size = size if size > overlap else max_bases
    step = (size - overlap) // WINDOW_GRID * WINDOW_GRID or size - overlap
    out, start = [], 0
    while True:
        out.append((start, min(length, start + size)))
        if out[-1][1] == length:
            return out
        start += step


def window_of(name):
    """A sequence's name in a piece (chunks) -> (its contig, the window's first base's 0-based place on it); (name, 0)
    for a contig that was not cut."""
    contig, mark, start = name.rpartition(WINDOW_MARK)
    return (contig, int(start)) if mark and start.isdigit() else (name, 0)


def chunks(genome_path, folder, max_bases=None, overlap=0, bases=None):
    """A genome as miniprot takes it: itself, or where it is over max_bases (CHUNK_BASES: a 3 Gb genome's index takes
    ~22 GB, a lungfish's 40 Gb ~300) its contigs in order in pieces of max_bases at most, plain FASTA in folder
    -> [paths]. A piece is closed before a contig that would take it past max_bases (its length from a first pass;
    until 2026-11-08 a piece took a whole further contig past it: chromosomes of 3 Gb made pieces of 6, ~45 GB of
    index where the build reserved ~33).
    A contig longer than max_bases is alone in its own piece (overlap 0: the six-frame ORFs' pieces, which no memory
    bound hangs on), or (overlap given; since 2026-11-10) cut into windows of max_bases, each `overlap` bases (half
    a piece at most) into the one before and each in a piece of its own, the last with the contigs after it; a
    window is named by window_of's rule, and loci.align puts its alignments back on the contig.
    bases: the genome's bases where the caller knows them (since 2026-11-10; NCBI's total length), else estimated
    from the file's size (genome_bases)."""
    max_bases = max_bases or CHUNK_BASES
    if (bases or genome_bases(genome_path)) <= max_bases:
        return [genome_path]
    overlap = min(overlap, max_bases // 2)
    lengths = iter(contig_lengths(genome_path))
    paths, out, held = [], None, 0
    cut, k, pos = [], 0, 0  # a cut contig's windows [(start, end, its piece's file)], the first not yet passed
    try:
        with (gzip.open(genome_path, "rt") if genome_path.endswith(".gz") else open(genome_path)) as fh:
            for line in fh:
                if line.startswith(">"):
                    for _a, _b, piece in cut:
                        if piece is not out:
                            piece.close()
                    name, n = next(lengths)
                    wins = windows(n, max_bases, overlap)
                    cut, k, pos = [], 0, 0
                    for start, end in wins:
                        if out is None or (held and held + end - start > max_bases):
                            if out is not None and not cut:  # (a window's piece is closed once it is written)
                                out.close()
                            paths.append(os.path.join(folder, f"piece{len(paths)}.fna"))
                            out, held = open(paths[-1], "w"), 0
                        held += end - start
                        out.write(line if len(wins) == 1 else f">{name}{WINDOW_MARK}{start}\n")
                        cut.append((start, end, out))
                    if len(wins) == 1:
                        cut = []  # whole: its lines as they are
                elif cut:
                    s = line.strip()
                    a, b = pos, pos + len(s)
                    while k < len(cut) and cut[k][1] <= a:
                        if cut[k][2] is not out:
                            cut[k][2].close()
                        k += 1
                    for start, end, piece in cut[k:]:
                        if start >= b:
                            break
                        piece.write(s[max(start, a) - a:min(end, b) - a] + "\n")
                    pos = b
                elif out is not None:  # (bases before any title: none of a contig)
                    out.write(line)
    finally:
        for _a, _b, piece in cut:
            piece.close()
        if out:
            out.close()
    return paths


INDEX_SUFFIX = ".mpi"
_INDEX_LOCKS, _INDEX_GUARD = {}, threading.Lock()


def index_file(index_dir, genome_path, key, piece=None, max_bases=None, overlap=0):
    """The file in index_dir of miniprot's saved index of a genome (piece None) or of one of its pieces (its number in
    chunks()' order, of max_bases at most, a longer contig's windows `overlap` into each other). Its name holds all
    that shapes the index: the genome (its path, size and time: another file under that name is another genome), the
    piece, and key (the tool and its indexing options: loci.split_options')."""
    stat = os.stat(genome_path)
    what = [os.path.abspath(genome_path), stat.st_size, stat.st_mtime_ns, *key]
    stem = os.path.basename(genome_path)
    if piece is not None:
        what += ["piece", max_bases or CHUNK_BASES, piece, *(["overlap", overlap] if overlap else [])]
        stem += f".piece{piece}"
    digest = hashlib.sha1("\n".join(str(x) for x in what).encode()).hexdigest()[:16]
    return os.path.join(index_dir, f"{stem}.{digest}{INDEX_SUFFIX}")


def saved_indexes(genome_path, folder, index_dir, key, build, max_bases=None, overlap=0, bases=None):
    """A genome as miniprot takes it where it is aligned more than once (since 2026-11-10: a dataset species' genome is
    called up to five times a run, the main call and its hold-outs, and each call made miniprot's index, most of a
    large genome's time, again): the saved indexes of its pieces (chunks()) in index_dir, each made once (index_file
    names it) -> ([their paths] in the pieces' order, why one could not be made or None).
    build(source, dest): writes the index of a FASTA to dest -> why it failed, or None. It writes a name of its own
    that is renamed once it is whole, so a failed or killed build leaves no file a later call takes for an index; one
    build at a time a genome in this process (threads); two processes at worst build twice, neither reads half a
    file. The pieces are written into folder and each deleted once its index is there; their number is noted beside
    the indexes, so a later call does not cut the genome again. Nothing is kept outside index_dir, which is the
    caller's to delete (~7.5 bytes a base). overlap, bases: chunks()' (since 2026-11-10: a contig over max_bases in
    windows, each window's piece its own index; the genome's bases where the caller knows them)."""
    max_bases = max_bases or CHUNK_BASES
    os.makedirs(index_dir, exist_ok=True)

    def names(n):
        return [index_file(index_dir, genome_path, key, i if n > 1 else None, max_bases, overlap) for i in range(n)]

    noted = index_file(index_dir, genome_path, ["pieces", max_bases, overlap])[:-len(INDEX_SUFFIX)] + ".pieces"
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
        parts = chunks(genome_path, folder, max_bases, overlap, bases)
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
