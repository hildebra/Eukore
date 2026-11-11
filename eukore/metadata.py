"""Eukore's assembly metadata (since 2026-11-07): what the build measures of each genome while its file is on disk,
since the prune (genomes.prune) deletes most of them, and NCBI's own summaries of them.

Of each genome Eukore reads (its run, genomes.run, or a later pass over the done genomes still on disk: describe_task),
into WORK/done/<accession>/:
- STATS (stats.json): one pass over the FASTA (scan): the file's MD5; its sequences, length, ACGT, N and other IUPAC
  bases; scaffold N50/L50/N90/L90 and the longest; contigs (split at runs of GAP_MIN Ns or more, NCBI's gaps) with
  their N50/L50, gaps and gap bases; the share of bases in sequences under 1 kb and 10 kb; GC (of ACGT), soft-masked
  (lower-case) share of ACGT, CpG observed/expected; GC of WINDOW-base windows (their spread, quantiles, and the share
  GC_OUTLIER or more from their median: another organism's contigs, isochores); sequences identical to an earlier one;
  telomeric repeats at the sequence ends (TELOMERES); the organelle molecules (FASTA titles naming a mitochondrion,
  chloroplast, plastid, apicoplast, kinetoplast, chromatophore, ...: organelle.TITLE); the counts of the 136 canonical
  tetranucleotides;
- ORGANELLE (organelle.fna.gz): those organelle molecules, whole (up to organelle.MAX_MOLECULE each);
- SKETCH (sketch.sig.gz): sourmash's FracMinHash sketch of the genome (SKETCH_PARAMS), where sourmash is there.
NCBI's assembly_summary rows of each genome are its tables', read at package time (ncbi_rows).

The package's metadata/ (write): assembly.tsv.gz (STATS, ASSEMBLY_COLUMNS), tetranucleotides.tsv.gz, ncbi.tsv.gz
(NCBI's columns as NCBI names them) and sketches.zip (a sourmash collection); organelle/molecules.fna.gz.
"""

import collections
import contextlib
import gzip
import hashlib
import itertools
import json
import os
import re
import socket
import subprocess
import threading
import time
import zipfile
import zlib

import numpy as np

from . import rrna, util

STATS, ORGANELLE, SKETCH = "stats.json", "organelle.fna.gz", "sketch.sig.gz"
FOLDER = "metadata"
CHUNK = 1 << 23          # bases of a sequence measured at a time
GAP_MIN = 10             # a run of Ns this long: an assembly gap (NCBI's, and package.GAP_MIN)
WINDOW = 10_000          # GC windows
GC_OUTLIER = 0.10        # a window's GC this far from the windows' median: an outlier
END = 5_000              # a sequence end searched for telomeric repeats (sequences of 2 * END or more)
TELOMERE_MIN = 500       # bases of one motif's copies (either strand) in an end: a telomere (by chance: ~15)
# the telomeric repeats known across eukaryotes (an end takes the one whose copies cover most of it: TTAGG's copies in
# a TTAGGG array cover 5/6 of what TTAGGG's do); yeasts' irregular TG(1-3) repeats are none of them
TELOMERES = ("TTAGGG",     # vertebrates, many fungi and protists, some invertebrates
             "TTTAGGG",    # most plants
             "TTTTAGGG",   # Chlamydomonas and other green algae
             "TTAGG",      # insects and other arthropods
             "TTAGGC",     # nematodes
             "TTGGGG",     # Tetrahymena and other ciliates
             "TTTTGGGG")   # hypotrich ciliates (Oxytricha, Euplotes)
SKETCH_PARAMS = "k=21,k=31,scaled=10000,noabund"
SKETCH_PART = ".partial.sig.gz"  # a sketch being written: .<host>-<pid>.sketch.partial.sig.gz (sourmash gzips a .gz)
ASSEMBLY_COLUMNS = (
    "accession", "fetched", "genome_source", "genome_url", "annotation_source", "md5", "bytes", "sequences", "length",
    "acgt", "n_bases", "other_bases", "longest",
    "n50", "l50", "n90", "l90", "contigs", "contig_n50", "contig_l50", "gaps", "gap_bases", "short_1kb", "short_10kb",
    "gc", "gc_windows", "gc_window_sd", "gc_window_q05", "gc_window_q95", "gc_outlier_windows", "masked", "cpg_oe",
    "duplicate_sequences", "duplicate_bases", "telomere_motif", "telomeric_ends", "telomere_both",
    "mitochondrion_molecules", "mitochondrion_length", "mitochondrion_gc", "plastid_molecules", "plastid_length",
    "plastid_gc", "sketch")
NCBI_EXTRA = ("genome", "fetched")  # ncbi.tsv.gz's own columns before NCBI's

COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
KMERS = ["".join(p) for p in itertools.product("ACGT", repeat=4)]  # index: the bases' codes, A0 C1 G2 T3, first high


def revcomp(kmer):
    return kmer.translate(str.maketrans("ACGT", "TGCA"))[::-1]


CANONICAL = sorted({min(k, revcomp(k)) for k in KMERS})  # 136
FOLD = np.array([CANONICAL.index(min(k, revcomp(k))) for k in KMERS])
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
MOTIFS = [(m.encode(), m.encode().translate(COMPLEMENT)[::-1]) for m in TELOMERES]
GAP = re.compile(rb"N{%d,}" % GAP_MIN)
HEADER = re.compile(rb"^>([^\n]*)\n", re.MULTILINE)


# --- one pass over a genome ---------------------------------------------------------------------------------------

def _pieces(fh, size=None):
    """A FASTA's binary stream -> (title, None) for each record's title line, (None, bytes) for its sequence in pieces
    (line ends left in). A sequence line need not end within a block: only a title line is held whole."""
    rest, size = b"", size or CHUNK
    while True:
        block = fh.read(size)
        if not block:
            break
        block = rest + block
        cut = block.rfind(b"\n") + 1
        rest = block[cut:]
        if rest and not rest.startswith(b">"):  # a sequence's line going on in the next block
            cut, rest = len(block), b""
        yield from _split(block[:cut])
    if rest:
        yield from _split(rest + b"\n")


def _split(text):
    pos = 0
    for m in HEADER.finditer(text):
        if m.start() > pos:
            yield None, text[pos:m.start()]
        yield m.group(1), None
        pos = m.end()
    if pos < len(text):
        yield None, text[pos:]


class _Hashed:
    """A binary file read through to its end, its bytes into a digest as they are read: the file's MD5 from the pass
    that scans it (since 2026-11-10: describe read each genome a second time for it)."""

    def __init__(self, fh, digest):
        self.fh, self.digest = fh, digest

    def read(self, size=-1):
        data = self.fh.read(size)
        self.digest.update(data)
        return data

    def rest(self):
        """What the reader left unread (nothing, where it read to the file's end), into the digest."""
        for block in iter(lambda: self.fh.read(1 << 22), b""):
            self.digest.update(block)


class _Feed(threading.Thread):
    """A file's bytes into a decompressor's stdin and into a digest, by a thread beside the scan (the tool reads a
    pipe then, not the file itself: Python sees the bytes it hashes). A tool that ends before the file does (a damaged
    file: its error is the scan's) is fed no more; `given_up` (the scan ended early) stops the reading."""

    def __init__(self, src, sink, digest):
        super().__init__(daemon=True)
        self.src, self.sink, self.digest = src, sink, digest
        self.given_up, self.error = threading.Event(), None

    def run(self):
        piped = True
        try:
            for block in iter(lambda: self.src.read(1 << 20), b""):
                if self.given_up.is_set():
                    break
                self.digest.update(block)
                if piped:
                    try:
                        self.sink.write(block)
                    except OSError:  # a broken pipe: the tool ended
                        piped = False
        except OSError as e:  # the file's read failed: the digest is not the file's (raised by _pieces_of)
            self.error = e
        finally:
            with contextlib.suppress(OSError):
                self.sink.close()


def _pieces_of(path, digest=None):
    """A FASTA's pieces (_pieces), gzipped or plain (decompressed by util.GZIP_TEST's igzip or pigz where there: a
    core of its own, ~3x Python's zlib) -> raises OSError where the decompressor fails. digest: a hashlib object that
    is given the file's own bytes (the compressed ones) as they are read, all of them once the pieces are through."""
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    if not gz or not util.GZIP_TEST:
        with open(path, "rb") as raw:
            src = raw if digest is None else _Hashed(raw, digest)
            if gz:
                with gzip.GzipFile(fileobj=src, mode="rb") as fh:
                    yield from _pieces(fh)
            else:
                yield from _pieces(src)
            if digest is not None:
                src.rest()
        return
    with open(path, "rb") as src:
        proc = subprocess.Popen([util.GZIP_TEST, "-dc"], stdin=src if digest is None else subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        feed = None
        if digest is not None:
            feed = _Feed(src, proc.stdin, digest)
            feed.start()
        through = False
        try:
            yield from _pieces(proc.stdout)
            through = True
        finally:
            if feed is not None and not through:
                feed.given_up.set()
            proc.stdout.close()
            err = proc.stderr.read().decode("utf-8", "replace")
            proc.stderr.close()
            proc.wait()
            if feed is not None:
                feed.join()
    if proc.returncode:
        raise OSError(f"{os.path.basename(util.GZIP_TEST)} -dc {path}: {err.strip()[-200:] or proc.returncode}")
    if feed is not None and feed.error is not None:
        raise OSError(f"{path}: {feed.error}")


def _nx(lengths, total, share):
    """(Nx, Lx) of lengths sorted longest first."""
    run = 0
    for i, n in enumerate(lengths, 1):
        run += n
        if run >= share * total:
            return n, i
    return 0, 0


def _telomere(end):
    """The motif whose copies (either strand) cover TELOMERE_MIN or more of an end's bases, most of them; or None."""
    best, cover = None, TELOMERE_MIN - 1
    for (motif, rc), name in zip(MOTIFS, TELOMERES):
        c = (end.count(motif) + end.count(rc)) * len(motif)
        if c > cover:
            best, cover = name, c
    return best


class Scan:
    """A genome's measures, its sequences fed in turn (start, feed, then finish)."""

    def __init__(self, title_of_organelle=None, max_molecule=None):
        self.t = collections.Counter()
        self.lengths, self.contigs, self.windows = [], [], []
        self.tetra = np.zeros(256, dtype=np.int64)
        self.ends, self.pairs = collections.Counter(), collections.Counter()
        self.seen, self.molecules = set(), []
        self.organelle_title, self.max_molecule = title_of_organelle, max_molecule
        self.name = None

    def start(self, title):
        self.close()
        self.t["sequences"] += 1
        self.name = (title.split() or [b""])[0].decode("utf-8", "replace")
        self.length = self.contig = self.wpos = self.wgc = self.wat = 0
        self.pending, self.prev, self.head, self.tail = b"", b"", b"", b""
        self.digest = hashlib.md5()
        # its organelle by rrna.organelle_of (since 2026-11-10; one reading of the title's word with organelle.py's)
        self.organelle = rrna.organelle_of(title.decode("utf-8", "replace"), self.organelle_title) \
            if self.organelle_title else None
        self.parts = [] if self.organelle else None

    def feed(self, data):
        if self.name is None:  # bases before any title
            return
        self.pending += data.translate(None, b"\n\r \t")
        if len(self.pending) >= CHUNK:
            up = self.pending.upper()
            cut = len(up.rstrip(b"N"))  # a run of Ns whole in one chunk: the next begins with it
            if cut:
                self._chunk(self.pending[:cut], up[:cut])
                self.pending = self.pending[cut:]

    def _chunk(self, raw, up):
        t, n = self.t, len(up)
        self.length += n
        byte = np.bincount(np.frombuffer(raw, dtype=np.uint8), minlength=256)  # one pass: every base, either case
        a, c, g, tt, nn = (int(byte[ord(b)] + byte[ord(b.lower())]) for b in "ACGTN")
        t["A"], t["C"], t["G"], t["T"], t["N"] = t["A"] + a, t["C"] + c, t["G"] + g, t["T"] + tt, t["N"] + nn
        t["other"] += n - a - c - g - tt - nn
        t["lower"] += int(sum(byte[ord(b)] for b in "acgt"))
        t["CG"] += up.count(b"CG") + (self.prev[-1:] == b"C" and up[:1] == b"G")
        prev = 0
        for m in GAP.finditer(up):
            s, e = m.span()
            self.contig += s - prev
            self._close_contig()
            t["gaps"] += 1
            t["gap_bases"] += e - s
            prev = e
        self.contig += n - prev
        codes = CODE[np.frombuffer(self.prev + up, dtype=np.uint8)]
        own = codes[len(self.prev):]
        # the GC windows: the chunk cut where the open window fills, then every WINDOW bases
        starts = np.concatenate(([0], np.arange(WINDOW - self.wpos, n, WINDOW)))
        gcs = np.add.reduceat((own == 1) | (own == 2), starts, dtype=np.int64)
        ats = np.add.reduceat((own == 0) | (own == 3), starts, dtype=np.int64)
        for gc, at, size in zip(gcs.tolist(), ats.tolist(), np.diff(np.append(starts, n)).tolist()):
            self.wgc, self.wat, self.wpos = self.wgc + gc, self.wat + at, self.wpos + size
            if self.wpos == WINDOW:
                self._close_window()
        if len(codes) >= 4:
            ok = codes < 4
            ok = ok[:-3] & ok[1:-2] & ok[2:-1] & ok[3:]
            k = (codes[:-3] << 6) | (codes[1:-2] << 4) | (codes[2:-1] << 2) | codes[3:]
            self.tetra += np.bincount(k[ok], minlength=256)
        self.prev = (self.prev + up[-3:])[-3:]
        if len(self.head) < END:
            self.head += up[:END - len(self.head)]
        self.tail = (self.tail + up[-END:])[-END:]
        self.digest.update(up)
        if self.parts is not None:
            self.parts.append(up)
            if self.length > self.max_molecule:
                self.parts = None

    def _close_contig(self):
        if self.contig:
            self.contigs.append(self.contig)
        self.contig = 0

    def _close_window(self, least=WINDOW // 2):
        if self.wgc + self.wat >= least:
            self.windows.append(self.wgc / (self.wgc + self.wat))
        self.wpos = self.wgc = self.wat = 0

    def close(self):
        """The sequence read so far ended."""
        if self.name is None:
            return
        if self.pending:
            self._chunk(self.pending, self.pending.upper())
            self.pending = b""
        self._close_contig()
        if self.wpos:
            self._close_window()
        if self.length:
            self.lengths.append(self.length)
            digest = self.digest.digest()
            if digest in self.seen:
                self.t["duplicate_sequences"] += 1
                self.t["duplicate_bases"] += self.length
            self.seen.add(digest)
        if self.length >= 2 * END:
            first, last = _telomere(self.head), _telomere(self.tail)
            for m in (first, last):
                if m:
                    self.ends[m] += 1
            self.pairs[(first, last)] += 1
        if self.parts is not None and self.length:
            self.molecules.append((self.name, self.organelle, b"".join(self.parts).decode("ascii", "replace")))
        self.name = None

    def finish(self):
        """-> (STATS's measures, the 136 canonical tetranucleotides' counts, organelle molecules [(name, organelle,
        sequence)])."""
        self.close()
        t = self.t
        total = sum(self.lengths)
        acgt = t["A"] + t["C"] + t["G"] + t["T"]
        lengths = sorted(self.lengths, reverse=True)
        contigs = sorted(self.contigs, reverse=True)
        n50, l50 = _nx(lengths, total, 0.5)
        n90, l90 = _nx(lengths, total, 0.9)
        c50, cl50 = _nx(contigs, sum(contigs), 0.5)
        windows = sorted(self.windows)

        def q(share):
            return round(windows[min(len(windows) - 1, int(share * len(windows)))], 4) if windows else "-"
        median = windows[len(windows) // 2] if windows else 0
        mean = sum(windows) / len(windows) if windows else 0
        motif = max(TELOMERES, key=lambda m: (self.ends[m], -TELOMERES.index(m))) if self.ends else None
        out = dict(sequences=t["sequences"], length=total, acgt=acgt, n_bases=t["N"], other_bases=t["other"],
                   longest=lengths[0] if lengths else 0, n50=n50, l50=l50, n90=n90, l90=l90, contigs=len(contigs),
                   contig_n50=c50, contig_l50=cl50, gaps=t["gaps"], gap_bases=t["gap_bases"],
                   short_1kb=round(sum(n for n in lengths if n < 1000) / total, 4) if total else 0,
                   short_10kb=round(sum(n for n in lengths if n < 10_000) / total, 4) if total else 0,
                   gc=round((t["G"] + t["C"]) / acgt, 4) if acgt else "-", gc_windows=len(windows),
                   gc_window_sd=round((sum((w - mean) ** 2 for w in windows) / len(windows)) ** 0.5, 4)
                   if windows else "-", gc_window_q05=q(0.05), gc_window_q95=q(0.95),
                   gc_outlier_windows=round(sum(abs(w - median) >= GC_OUTLIER for w in windows) / len(windows), 4)
                   if windows else "-",
                   masked=round(t["lower"] / acgt, 4) if acgt else "-",
                   cpg_oe=round(t["CG"] * acgt / (t["C"] * t["G"]), 4) if t["C"] and t["G"] else "-",
                   duplicate_sequences=t["duplicate_sequences"], duplicate_bases=t["duplicate_bases"],
                   telomere_motif=motif or "-", telomeric_ends=self.ends[motif] if motif else 0,
                   telomere_both=self.pairs[(motif, motif)] if motif else 0)
        for kind in ("mitochondrion", "plastid"):
            seqs = [s for _n, k, s in self.molecules if k == kind]
            bases = sum(len(s) for s in seqs)
            gc = sum(s.count("G") + s.count("C") for s in seqs)
            at = sum(s.count("A") + s.count("T") for s in seqs)
            out.update({f"{kind}_molecules": len(seqs), f"{kind}_length": bases,
                        f"{kind}_gc": round(gc / (gc + at), 4) if gc + at else "-"})
        tetra = np.bincount(FOLD, weights=self.tetra, minlength=len(CANONICAL)).astype(np.int64)
        return out, [int(x) for x in tetra], self.molecules


def scan(path, digest=None):
    """A genome's FASTA (gzipped or plain) -> Scan.finish(). digest: a hashlib object given the file's bytes in the
    same pass (_pieces_of): its MD5."""
    from . import organelle  # noqa: PLC0415
    s = Scan(organelle.TITLE, organelle.MAX_MOLECULE)
    for title, data in _pieces_of(path, digest):
        if title is not None:
            s.start(title)
        else:
            s.feed(data)
    return s.finish()


# --- a genome's metadata into its done folder -------------------------------------------------------------------

def missing(folder, sourmash=None, names=None):
    """What a done genome's results (its folder, or the packed genome's: store.py) lack: STATS (and ORGANELLE),
    SKETCH (where sourmash is there) -> a set. names: its files' names where read already (genomes.names_of)."""
    from . import store  # noqa: PLC0415

    def has(name):
        return name in names if names is not None else store.has(folder, name)
    out = set()
    if not has(STATS):
        out.add(STATS)
    if sourmash and not has(SKETCH):
        out.add(SKETCH)
    return out


def describe(path, folder, accession, fetched=None, sourmash=None, want=None):
    """A genome's metadata into folder (want: of STATS and SKETCH, default both; SKETCH only where sourmash is
    given), sourmash's sketch made beside the scan -> why it failed, or None. The file's MD5 comes from the scan's own
    pass (since 2026-11-10: a second read of the genome). The sketch is written under this process' own name (host
    and pid; since 2026-11-10: one name, .sketch.partial.sig.gz, which two processes describing a genome both wrote)
    and, whatever fails, sourmash is not left running nor its partial left there."""
    want = {STATS, SKETCH} if want is None else set(want)
    started = time.time()
    part = proc = None
    if SKETCH in want and sourmash:
        for dead in util.dead_partials(folder, "sketch", SKETCH_PART):  # those of runs that were killed
            with contextlib.suppress(OSError):
                os.remove(dead)
        part = os.path.join(folder, f".{socket.gethostname()}-{os.getpid()}.sketch{SKETCH_PART}")
        proc = subprocess.Popen([sourmash, "sketch", "dna", "-p", SKETCH_PARAMS, "--name", accession, "-o", part, path],
                                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    why = None
    try:
        if STATS in want:
            digest = hashlib.md5()
            try:
                stats, tetra, molecules = scan(path, digest)
            # zlib.error: a damaged gzip stream read by Python's gzip (until 2026-11-10 not caught: it left the
            # function with sourmash running and its partial there)
            except (OSError, EOFError, ValueError, zlib.error) as e:
                why = f"metadata: {type(e).__name__}: {e}"
            else:
                with util.bgzf_open(util.part_of(os.path.join(folder, ORGANELLE)), level=util.WORK_LEVEL) as fh:
                    fh.writelines(f">{n} {k}\n{s}\n" for n, k, s in molecules)
                os.replace(util.part_of(os.path.join(folder, ORGANELLE)), os.path.join(folder, ORGANELLE))
                stats = dict(accession=accession, fetched=fetched or accession, md5=digest.hexdigest(),
                             bytes=os.path.getsize(path), **stats, tetranucleotides=tetra,
                             seconds=round(time.time() - started, 1))
                with open(util.part_of(os.path.join(folder, STATS)), "w", newline="\n") as fh:
                    json.dump(stats, fh)
                os.replace(util.part_of(os.path.join(folder, STATS)), os.path.join(folder, STATS))
        if proc:
            err = proc.communicate()[1]
            if proc.returncode == 0 and os.path.isfile(part):
                os.replace(part, os.path.join(folder, SKETCH))
            else:
                why = why or f"sourmash sketch failed (exit {proc.returncode}): {err.strip()[-200:]}"
    finally:
        if proc:
            if proc.poll() is None:  # an error left the block
                proc.kill()
                proc.communicate()
            if os.path.isfile(part):
                os.remove(part)
    return why


def describe_task(args):
    """A pool's task: (accession, genome path, done folder, fetched accession, sourmash) -> (accession, why)."""
    acc, path, folder, fetched, sourmash = args
    try:
        return acc, describe(path, folder, acc, fetched, sourmash, missing(folder, sourmash))
    except Exception as e:  # noqa: BLE001  a genome's failure, never the pool's
        return acc, f"metadata: {type(e).__name__}: {e}"


# --- the package's metadata/ ---------------------------------------------------------------------------------------

def ncbi_rows(paths, wanted):
    """NCBI's assembly_summary tables (GenBank's and RefSeq's) -> ({accession: row of all NCBI's columns}, the columns
    in NCBI's order) of the accessions wanted."""
    rows, columns = {}, []
    for path in paths:
        with util.open_text(path) as fh:
            head = None
            for line in fh:
                if line.startswith("#"):
                    head = line.lstrip("#").strip().split("\t")
                    continue
                if head is None:
                    continue
                f = line.rstrip("\n").split("\t")
                if f[0] in wanted:
                    rows[f[0]] = dict(zip(head, f))
                    columns += [c for c in head if c not in columns]
    return rows, columns


def _gz_table(path, columns, rows):
    from .package import write_gz  # noqa: PLC0415 (the package's gzipped tables, one writer)
    write_gz(path, columns, rows)


ANNOTATION_NAMES = {"veupathdb": "VEuPathDB", "ensembl": "Ensembl", "toga2": "TOGA2"}


def genome_source(row, job=None):
    """The database a genome was downloaded from (since 2026-11-07): `NCBI RefSeq` or `NCBI GenBank` (by the
    accession fetched), `TOGA2 (<who assembled it>)`, `EukProt` (a species without an assembly: EukDetect2's records
    of EukProt's genomes and transcriptomes, data/no_assembly.fna.gz)."""
    why = (job or {}).get("why", [row.get("source", "")])[0]
    if why == "toga2":
        who = (job or {}).get("assembly_source", "-")
        return f"TOGA2 ({who})" if who not in ("-", "") else "TOGA2"
    if why == "tsa":  # since 2026-11-07: INSDC's transcriptome, as ENA serves it
        return "ENA TSA"
    if why == "no_assembly":
        return "EukProt"
    fetched = row.get("fetched") or row.get("accession", "")
    return "NCBI RefSeq" if fetched.startswith("GCF_") else "NCBI GenBank" if fetched.startswith("GCA_") else "-"


def annotation_source(work, row, job=None):
    """The database of the gene models a genome's run took (genomes.annotation_record): `NCBI RefSeq` or `NCBI
    GenBank` (NCBI's GFF of the assembly fetched), `VEuPathDB`, `Ensembl`, `TOGA2`; `-` for none (miniprot's alone)."""
    from .genomes import annotation_record  # noqa: PLC0415
    got = annotation_record(work, job or dict(accession=row["accession"]))
    used = got[0] if got else "-"
    if used == "ncbi":
        return "NCBI RefSeq" if (row.get("fetched") or "").startswith("GCF_") else "NCBI GenBank"
    return ANNOTATION_NAMES.get(used, "-")


def write(part, genomes, work, tables=(), jobs=None):
    """The package's metadata (in its unfinished folder part) of its genomes (package.collect's {accession: row};
    jobs: the plan's {accession: job}, where the sources' details are): metadata/ncbi.tsv.gz (NCBI's rows of every
    genome and of the RefSeq copy fetched); metadata/assembly.tsv.gz, a row per species' genome (where it came from,
    and the measures of those Eukore described); of those, metadata/tetranucleotides.tsv.gz, metadata/sketches.zip and
    organelle/molecules.fna.gz -> counts for the manifest."""
    from . import store  # noqa: PLC0415
    from .genomes import done_dir  # noqa: PLC0415
    jobs = jobs or {}
    folder = os.path.join(part, FOLDER)
    os.makedirs(folder, exist_ok=True)
    os.makedirs(os.path.join(part, "organelle"), exist_ok=True)
    fetched = {a: g.get("fetched") or a for a, g in genomes.items()}
    rows, columns = ncbi_rows(tables, set(fetched) | set(fetched.values()))
    out = []
    for acc in sorted(genomes):
        for name in dict.fromkeys((acc, fetched[acc])):
            if name in rows:
                out.append(dict(rows[name], genome=acc, fetched="t" if name == fetched[acc] else "f"))
    _gz_table(os.path.join(folder, "ncbi.tsv.gz"), NCBI_EXTRA + tuple(c for c in columns if c not in NCBI_EXTRA), out)
    described, listed, sketches = [], [], []
    molecules = os.path.join(part, "organelle", "molecules.fna.gz")
    n_molecules = 0
    with util.bgzf_open(util.part_of(molecules)) as mol:
        for acc in sorted(genomes):
            if genomes[acc].get("representative") != "t":
                continue
            d = done_dir(work, acc)  # its folder, or its packed file (store.py; since 2026-11-10)
            try:
                stats = json.loads(store.text(d, STATS) or "")
            except ValueError:
                stats = None
            sketched = store.has(d, SKETCH)
            if sketched:
                sketches.append((acc, d))
            g, job = genomes[acc], jobs.get(acc)
            row = dict(stats or {}, accession=acc, fetched=g.get("fetched") or acc,
                       genome_source=genome_source(g, job), genome_url=(job or {}).get("url") or "-",
                       annotation_source=annotation_source(work, g, job),
                       sketch="t" if sketched else "f")
            listed.append(row)
            if stats is None:
                continue
            described.append(row)
            for line in store.open_text(d, ORGANELLE) or ():
                if line.startswith(">"):
                    n_molecules += 1
                    line = f">{acc}|{line[1:]}"
                mol.write(line)
    os.replace(util.part_of(molecules), molecules)
    _gz_table(os.path.join(folder, "assembly.tsv.gz"), ASSEMBLY_COLUMNS, listed)
    _gz_table(os.path.join(folder, "tetranucleotides.tsv.gz"), ("accession",) + tuple(CANONICAL),
              [dict(accession=r["accession"], **dict(zip(CANONICAL, r["tetranucleotides"]))) for r in described])
    if sketches:  # a zip of signatures sourmash reads as one collection (stored: each is gzipped)
        path = os.path.join(folder, "sketches.zip")
        with zipfile.ZipFile(util.part_of(path), "w", zipfile.ZIP_STORED) as z:
            for acc, d in sketches:
                z.writestr(f"signatures/{acc}.sig.gz", store.read(d, SKETCH))
        os.replace(util.part_of(path), path)
    return {"ncbi_rows": len(out), "genomes": len(listed), "assemblies_described": len(described),
            "sketches": len(sketches), "organelle_molecules": n_molecules,
            "genome_sources": dict(collections.Counter(r["genome_source"] for r in listed)),
            "annotation_sources": dict(collections.Counter(r["annotation_source"] for r in listed))}
