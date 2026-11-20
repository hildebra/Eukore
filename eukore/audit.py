"""The audit of the published rRNA records (since 2026-11-14, format 7): KSGP 4.1-4.5's record checks (KSGP 5.2's
README and KSGP 4.5's: 4.1's record audit, 4.4's 18S gene check, 4.5's flank screen; the chimera screens are in
chimeras.py) rebuilt for Eukore's published reads (EUKARYOME's long reads and SSU set, ROD's operons), before the checks
every record has (rrna.finish). What it finds beside a record's rRNA it cuts off; what it finds in its rRNA it flags
(records.tsv's `trimmed`, `flags`). A read's rRNA is what barrnap (and ITSx) placed in it (published_regions): its unit
(18S to 28S) or, without one, its 18S; the read is on that rRNA's strand.

ends
  - vectors, adapters, linkers and sequencing primers: NCBI's UniVec_Core (KSGP 4.1's pGEM/pUC polylinkers, T7, M13)
    and the long-read platforms' adapters (ADAPTERS), exact matches of K bases or more (K-mers of low complexity left
    out: oligo-dT primers match any poly-A run). Beside the rRNA or within END bases of a read's end, the read is cut
    there (and everything beyond); inside the rRNA a match of INSIDE bases or more is flagged `vector` (low).
  - a barcode or an adapter before an rDNA amplification primer (data/primers.tsv's 18S, ITS and 28S primers) at a
    read's start, or after one at its end: the bases between the primer (MAX_PRIMER mismatches or fewer, more than
    PRIMER_SLACK and at most MAX_SLACK bases in) and the read's end, where they are no rRNA, cut. The primer is kept:
    a gene-end primer's site is the 18S's own end, and most records end there (EUKARYOME's SSU set: 364,417 end 3
    bases past EukB's site, the 18S's last bases), primer-made or not.
  - N runs and a homopolymer of HOMOPOLYMER_END or more (a poly-A tail) at an end, beside the rRNA: cut.
  - a telomere repeat (rrna.telomere's units) beside the rRNA: the read cut there (KSGP 4.1 masked the two it found).
flanks (KSGP 4.4's 18S gene check, 4.5's flank screen)
  - a read's bases beside its rRNA, more than FLANK at an end, are searched again with the read whole (barrnap's
    models, nhmmer -E NHMMER_E: KSGP 4.5's second pass). Where the rRNA reaches its model's end on that side (within
    MODEL_END positions) the flank lies beyond the gene (an ETS, an IGS, an ITS after an 18S record's whole 18S, a
    leader) and is cut, keeping KEEP bases beside the rRNA and the model positions it did not reach; where the flank
    holds another rRNA gene, a reversed rRNA piece or a second unit (two parts joined) it is cut as well. Else (the
    rRNA stops short and other DNA follows: a clone's end, or a divergent stretch the model does not align) it is
    kept and flagged `flank` (medium): KSGP 4.5 kept such a flank where a relative's rRNA goes on into it and cut the
    others, a test Eukore does not make. After an operon read's partial 28S it is kept unflagged: the 28S going on to
    the primer site the read ends at (EUKARYOME 2.1: 255,487 long reads end at one 28S site), its divergent
    stretches aligned weakly by the model. An 18S read (no operon) whose flank holds its ITS's 5.8S or its 28S is cut
    to its 18S (KSGP 4.4's rule; a long read whose operon barrnap did not make whole, its 5.8S not found, so too).
    Since 2026-11-20, before these: where nhmmer finds the gene going on into the flank past barrnap's end, its
    region grows to it (a stretch the model aligns weakly), and where the gene goes on only after INSERTION bases or
    more than the model skips (Drosophila's R1/R2 retrotransposons in the 28S), the read is cut before the insertion
    (`flank_insertion`; its 28S partial; nhmmer's one hit across the insertion cut at barrnap's end). The rest of
    the flank, beyond the gene, is then judged as above; one kept where the rRNA stops short is cut where it is
    longer than twice the model positions not reached and INSERTION (`flank_long`: more than the rest of the gene
    could be).
  - a read in which barrnap finds no rRNA: flagged `no_rrna` (low; KSGP 4.4 dropped such 18S).
A read of which fewer than MIN_LEFT bases would stay is not cut, flagged `trim_short` (low).
"""

import collections
import os
import re
import shutil
import subprocess
import tempfile

import numpy as np

from . import primers as primers_mod, rrna, util, versions

K = 25              # a vector's match: this many bases or more (NCBI's VecScreen: a terminal match of 24-25 bases)
LOW_SHARE = 0.6     # a K-mer of low complexity: this share of it one base
INSIDE = 50         # a vector inside a read's rRNA: a match this long at least (VecScreen's internal match: 30)
END = 60            # a vector within this many bases of a read's end, inside its rRNA, still an end's
FLANK = 50          # a flank of more bases searched again (KSGP 4.5's)
KEEP = 10           # bases kept beside the rRNA of a flank cut off (KSGP 4.5's)
MODEL_END = 60      # an rRNA hit this close to its model's end reaches it (KSGP 4.5's)
NHMMER_E = 1e-3     # the second pass' E-value (KSGP 4.5's; barrnap's own 1e-6)
MAX_PRIMER = 2
PRIMER_SLACK = 3
HOMOPOLYMER_END = 12
MIN_LEFT = 300      # KSGP 4.5's flank screen: a record of fewer bases left is dropped
# a read's bases between two pieces of one rRNA gene beyond the model positions between them: this many more make an
# insertion (since 2026-11-20; flank_cut). The 2026.11-013 package's ROD and EUKARYOME reads: a 28S barrnap stopped
# short of went on within 5 bases of its model; Drosophila's R1/R2 insertions in the 28S held 1,388-8,319 more
INSERTION = 500
VECTOR_BATCH = 1 << 23  # bases encoded at a time (a uint64 a base)
ADAPTERS = {  # the long-read platforms' adapters, which UniVec_Core (2017) lacks
    "PacBio_SMRTbell_hairpin": "ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
    "ONT_ligation_Y_top": "AATGTACTTCGTTCAGTTACGTATTGCT",
    "ONT_rapid_adapter": "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA",
}
RDNA_TARGETS = ("18S", "ITS", "28S")
TRIM_COLUMNS = ("record", "source", "length_before", "start", "end", "why")


# --- vectors ---------------------------------------------------------------------------------------------------------

def _codes(data, k=K):
    """bytes -> (each position's k-mer code (uint64; 2 bits a base), whether that k-mer is all ACGT)."""
    c = rrna._CODE[np.frombuffer(bytes(data), dtype=np.uint8)]
    n = len(c) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=bool)
    bits = (c & 3).astype(np.uint64)
    code = bits[:n].copy()
    for j in range(1, k):
        code <<= np.uint64(2)
        code |= bits[j:j + n]
    bad = np.concatenate(([0], np.cumsum(c > 3, dtype=np.int64)))
    return code, bad[k:k + n] == bad[:n]


def low_complexity(data, k=K):
    """bytes -> whether each position's k-mer is of low complexity (LOW_SHARE of it one base, or two bases alone):
    UniVec's oligo-dT primers (SMART CDS III, the I.M.A.G.E. libraries' NotI oligo-dT) and poly-A vectors (pSP64)
    match any poly-A run (EUKARYOME's SSU set: most of its 328 records with a UniVec 25-mer, 2026-11-14)."""
    c = rrna._CODE[np.frombuffer(bytes(data), dtype=np.uint8)]
    n = len(c) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=bool)
    counts = np.stack([np.concatenate(([0], np.cumsum(c == b))) for b in range(4)])
    window = counts[:, k:k + n] - counts[:, :n]
    return (window.max(axis=0) >= LOW_SHARE * k) | ((window > 0).sum(axis=0) <= 2)


class Vectors:
    """UniVec_Core's sequences and ADAPTERS as a table of their K-mers (both strands), but those of low
    complexity."""

    def __init__(self, entries):
        codes, names = [], []
        for name, seq in list(entries) + list(ADAPTERS.items()):
            for s in (seq.upper(), rrna.revcomp(seq.upper())):
                c, ok = _codes(s.encode())
                ok &= ~low_complexity(s.encode())[:len(ok)]
                codes.append(c[ok])
                names.append(np.full(int(ok.sum()), len(names), dtype=np.int32))
        self.labels = [n for n, _s in list(entries) + list(ADAPTERS.items()) for _ in (0, 1)]
        code = np.concatenate(codes) if codes else np.zeros(0, dtype=np.uint64)
        label = np.concatenate(names) if names else np.zeros(0, dtype=np.int32)
        order = np.argsort(code, kind="stable")
        self.code, self.label = code[order], label[order]
        self.entries = len(self.labels) // 2

    def hits(self, seqs):
        """[sequence] -> {index: [(start, end, entry name)]}: each sequence's stretches of matched K-mers, merged."""
        out = {}
        batch, at = [], 0
        for i, s in enumerate(seqs):
            batch.append((i, s))
            at += len(s) + 1
            if at >= VECTOR_BATCH:
                out.update(self._hits(batch))
                batch, at = [], 0
        if batch:
            out.update(self._hits(batch))
        return out

    def _hits(self, batch):
        data = b"N".join(s.encode() for _i, s in batch)
        starts = np.cumsum([0] + [len(s) + 1 for _i, s in batch[:-1]])
        code, ok = _codes(data)
        if not len(self.code) or not len(code):
            return {}
        pos = np.searchsorted(self.code, code)
        pos[pos >= len(self.code)] = len(self.code) - 1
        found = np.nonzero(ok & (self.code[pos] == code))[0]
        out = collections.defaultdict(list)
        for p in found:
            j = int(np.searchsorted(starts, p, side="right")) - 1
            i, a = batch[j][0], int(p - starts[j])
            name = self.labels[int(self.label[pos[p]])]
            spans = out[i]
            if spans and a <= spans[-1][1]:
                spans[-1] = (spans[-1][0], a + K, spans[-1][2])
            else:
                spans.append((a, a + K, name))
        return out


def rdna_primers(path=primers_mod.PRIMERS):
    """data/primers.tsv's rDNA primers -> [(name, sequence as the read's 5' end holds it, end: 5 or 3)]: a forward
    primer at a read's start, a reverse one's reverse complement at its end."""
    out = []
    for r in primers_mod.read_tsv(path):
        if r["target"] in RDNA_TARGETS:
            fwd = r["strand"] == "forward"
            out.append((r["name"], r["sequence"].upper() if fwd else primers_mod.revcomp(r["sequence"]),
                        5 if fwd else 3))
    return out


MAX_SLACK = 60  # a primer this far in at most (behind a barcode and an adapter)


def primer_ends(seqs, primers, slacks=None):
    """[sequence], rdna_primers, {index: (bases before its rRNA at its start, after it at its end)} -> {index: {5 or 3:
    (primer, offset from that end, its length, mismatches)}}: the best primer (fewest mismatches, then nearest the end)
    more than PRIMER_SLACK bases in from an end and no further than its rRNA begins there (slacks; MAX_SLACK at
    most), MAX_PRIMER mismatches or fewer."""
    slacks = slacks or {}
    longest = max((len(p) for _n, p, _e in primers), default=0)
    out = collections.defaultdict(dict)
    for end in (5, 3):
        mine = [(name, np.array([primers_mod.BITS.get(c, 15) for c in seq], dtype=np.uint8))
                for name, seq, e in primers if e == end]  # as the read's end holds it, read from that end inwards
        if not mine:
            continue
        best = {}

        def take(i, key, hit):
            if i not in best or key < best[i][0]:
                best[i] = (key, hit)
        for i, sl in slacks.items():  # only where the bases before it are no rRNA
            reach = min(MAX_SLACK, sl[0] if end == 5 else sl[1])
            if reach <= PRIMER_SLACK:
                continue
            s = seqs[i]
            bits = primers_mod.encode(s[:reach + longest] if end == 5 else s[-(reach + longest):][::-1])
            for name, p in mine:
                q = p if end == 5 else p[::-1]  # a 3' end read backwards: the primer's reverse complement, reversed
                if len(bits) < len(q):
                    continue
                win = np.lib.stride_tricks.sliding_window_view(bits, len(q))[PRIMER_SLACK + 1:reach + 1]
                if not len(win):
                    continue
                miss = ((win & q) == 0).sum(axis=1)
                k = int(miss.argmin())
                if miss[k] <= MAX_PRIMER:
                    off = PRIMER_SLACK + 1 + k
                    take(i, (int(miss[k]), off), (name, off, len(q), int(miss[k])))
        for i, (_k, hit) in best.items():
            out[i][end] = hit
    return out


# --- flanks ----------------------------------------------------------------------------------------------------------

def barrnap_hmm(barrnap, kingdom="euk"):
    """barrnap's HMM file of a kingdom (its db/<kingdom>.hmm), or None."""
    if not barrnap:
        return None
    real = os.path.realpath(barrnap)
    for d in (os.path.join(os.path.dirname(real), "..", "db"), os.path.join(os.path.dirname(real), "..", "lib",
                                                                              "barrnap", "db"),
              os.path.join(os.path.dirname(barrnap), "..", "lib", "barrnap", "db")):
        path = os.path.normpath(os.path.join(d, f"{kingdom}.hmm"))
        if os.path.isfile(path):
            return path
    return None


def model_lengths(hmm):
    """An HMMER3 file's models -> {name: length}."""
    out, name = {}, None
    with open(hmm) as fh:
        for line in fh:
            if line.startswith("NAME"):
                name = line.split()[1]
            elif line.startswith("LENG") and name:
                out[name] = int(line.split()[1])
    return out


def nhmmer_hits(reads, nhmmer, hmm, threads=1, tmp=None, evalue=NHMMER_E):
    """reads [(name, sequence)] -> {name: [(model, hmm_from, hmm_to, start, end, strand)]} (start/end 0-based, end
    exclusive, on the read as given; hmm positions 1-based)."""
    out = collections.defaultdict(list)
    if not reads:
        return out
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        fasta, tbl = os.path.join(d, "reads.fna"), os.path.join(d, "hits.tbl")
        rrna.barrnap_fasta(fasta, reads)
        subprocess.run([nhmmer, "--noali", "--cpu", str(max(1, threads)), "-E", str(evalue), "--tblout",
                        tbl, "-o", os.devnull, hmm, fasta], check=True)
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                a, b = int(f[6]), int(f[7])
                strand = f[11]
                out[f[0]].append((f[2], int(f[4]), int(f[5]), min(a, b) - 1, max(a, b), strand))
    return out


def flank_cut(length, unit, hits, lengths, side):
    """One side of a read (side 5 or 3; the read on its rRNA's strand): (how many bases of it to cut, why, how many
    bases its rRNA grows by on that side) or (0, `flank` where kept and flagged, or None where nothing beyond FLANK lies
    there, ...). unit: its rRNA [start, end) (barrnap's); hits: nhmmer_hits' of the read; lengths: the models'
    lengths. The bases cut are counted from nhmmer's rRNA on that side, which may reach further into the flank than
    barrnap's (the second pass' E-value).
    Since 2026-11-20 the gene going on into the flank (nhmmer's pieces of it, model positions after barrnap's) is
    the rRNA's, piece by piece, while the bases between two pieces are about the model positions between them (a
    stretch the model aligns weakly: its region grows to it); where they are INSERTION more, the gene holds an
    insertion (Drosophila's R1/R2 retrotransposons in the 28S: 3-5 kb) and the read is cut before it
    (`flank_insertion`; so too where nhmmer's one hit spans it, INSERTION bases longer on the read than on its model
    and as many past barrnap's end: cut at barrnap's end). Until then such a flank was kept whole (Drosophila's ROD
    units ended 3-5 kb past their region, the insertion and the 28S's 3' part beyond it). A flank kept where the
    rRNA stops short is cut where it is longer than twice the model positions not reached and INSERTION
    (`flank_long`: more than the rest of the gene could be; an 18S read with 4 kb of no rRNA after it)."""
    u0, u1 = unit
    flank = u0 if side == 5 else length - u1
    if flank <= FLANK:
        return 0, None, 0
    own = [h for h in hits if h[5] == "+" and h[3] < u1 and h[4] > u0]
    edge = (min(own, key=lambda h: h[3]) if side == 5 else max(own, key=lambda h: h[4])) if own else None
    grow = 0
    if edge is not None:
        grow = max(0, u0 - edge[3] if side == 5 else edge[4] - u1)  # nhmmer's rRNA reaching further than barrnap's
        if grow > INSERTION and (edge[4] - edge[3]) - (edge[2] - edge[1] + 1) > INSERTION:
            # nhmmer's one hit across an insertion (its alignment's gap; Drosophila's ROD units): the read cut at
            # barrnap's end, the insertion's start
            return (u0 if side == 5 else length - u1) - KEEP, "flank_insertion", 0

        def goes_on(h, e):  # h the same gene's, further on the model and on the read than e
            return h[5] == "+" and h[0] == e[0] and (h[2] > e[2] and h[4] > e[4] if side == 3 else
                                                     h[1] < e[1] and h[3] < e[3])
        for h in sorted(hits, key=lambda h: h[3] if side == 3 else -h[4]):
            if not goes_on(h, edge):
                continue
            bases, positions = (h[3] - edge[4], h[1] - edge[2] - 1) if side == 3 else (edge[3] - h[4],
                                                                                        edge[1] - h[2] - 1)
            if bases - positions > INSERTION:
                return (edge[3] if side == 5 else length - edge[4]) - KEEP, "flank_insertion", grow
            edge = h
            grow = max(0, u0 - edge[3] if side == 5 else edge[4] - u1)
        # the bases beside nhmmer's rRNA
        flank = min(flank, edge[3] if side == 5 else length - edge[4])
        if flank <= FLANK:
            return 0, None, grow

    def kept():  # the flank where the rRNA stops short, kept: unless more than the rest of its gene could be
        if edge is not None:
            unreached = edge[1] - 1 if side == 5 else lengths.get(edge[0], edge[2]) - edge[2]
            if flank > 2 * unreached + INSERTION:
                return flank - KEEP, "flank_long", grow
        return 0, "flank", grow
    lo, hi = (0, flank) if side == 5 else (length - flank, length)
    inside = [h for h in hits if h[3] < hi and h[4] > lo and min(h[4], hi) - max(h[3], lo) >= 20]
    if any(h[5] == "-" for h in inside):
        return flank - KEEP, "flank_joined", grow  # a reversed rRNA piece beside it: two parts joined
    if edge is not None:
        if any(h[0] != edge[0] for h in inside):  # another gene beside it (an 18S record's ITS's 5.8S, a second unit)
            return flank - KEEP, "flank_gene", grow
        same = [h for h in inside if h[0] == edge[0]]
        if same and (side == 5 and min(h[1] for h in same) < edge[1] or side == 3 and max(h[2] for h in same) >
                     edge[2]):
            return kept()  # the same gene goes on into the flank: the rRNA stops short
        reached = edge[1] <= MODEL_END if side == 5 else edge[2] >= lengths.get(edge[0], 10 ** 9) - MODEL_END
        if reached:
            unreached = edge[1] - 1 if side == 5 else lengths.get(edge[0], edge[2]) - edge[2]
            return max(0, flank - KEEP - unreached), "flank_beyond", grow
    return kept()


# --- the audit -------------------------------------------------------------------------------------------------------

def unit_of(row, length):
    """A record's rRNA [start, end) (its regions' extent), or None."""
    regions = row.get("regions_list") or ()
    if not regions:
        return None
    return max(0, min(a for _n, a, _b in regions)), min(length, max(b for _n, _a, b in regions))


def cut_row(row, a, b):
    """A record's row with its read cut to [a, b): its regions shifted and clipped."""
    regions = [(n, max(x, a) - a, min(y, b) - a) for n, x, y in row.get("regions_list") or () if min(y, b) > max(x, a)]
    return dict(row, regions_list=regions or None)


TELOMERE_RUN = re.compile("|".join(f"(?:{u}){{{rrna.TELOMERE},}}|(?:{rrna.revcomp(u)}){{{rrna.TELOMERE},}}"
                                   for u in rrna.TELOMERE_UNITS))
END_RUN = re.compile(rf"(N+|A{{{HOMOPOLYMER_END},}}|C{{{HOMOPOLYMER_END},}}|G{{{HOMOPOLYMER_END},}}|"
                     rf"T{{{HOMOPOLYMER_END},}})")


def end_cuts(seq, unit, vector_spans, primer_hits):
    """The ends' cuts of one read: (bases off its start, bases off its end, [why], vector inside its rRNA (name or
    None))."""
    length = len(seq)
    u0, u1 = unit if unit else (0, length)
    a, b, why, inside = 0, length, [], None
    for s, e, name in vector_spans:
        if e <= u0 + END or e <= END:  # before the rRNA, or at the read's start
            if e > a:
                a = e
                why.append(f"5':vector {name}")
        elif s >= u1 - END or s >= length - END:
            if s < b:
                b = s
                why.append(f"3':vector {name}")
        elif e - s >= INSIDE:
            inside = inside or name
    # a leader before a primer at the start (after one at the end) that is no rRNA: a barcode, an adapter. The
    # primer itself is kept: a gene-end primer's site is the 18S's own end (EUKARYOME's SSU set: 364,417 of its
    # records end 3 bases past EukB's site, the 18S's last bases; 2026-11-14)
    for end, (name, off, _plen, _mm) in sorted(primer_hits.items()):
        if end == 5 and off > PRIMER_SLACK and u0 >= off - PRIMER_SLACK and off > a:
            a = off
            why.append(f"5':leader before primer {name}")
        elif end == 3 and off > PRIMER_SLACK and length - u1 >= off - PRIMER_SLACK and length - off < b:
            b = length - off
            why.append(f"3':leader after primer {name}")
    for m in TELOMERE_RUN.finditer(seq):
        if m.end() <= u0 and m.end() > a:
            a = m.end()
            why.append("5':telomere")
        elif m.start() >= u1 and m.start() < b:
            b = m.start()
            why.append("3':telomere")
    # Ns at an end, and a homopolymer there beside the rRNA (a poly-A tail), one after the other
    for _ in range(3):
        m = END_RUN.match(seq, a, b)
        if not m or not (m.group(1)[0] == "N" or m.end() <= u0):
            break
        a = m.end()
        why.append("5':" + ("Ns" if m.group(1)[0] == "N" else "homopolymer"))
    for _ in range(3):
        m = END_RUN.match(seq[a:b][::-1])
        if not m or not (m.group(1)[0] == "N" or b - m.end() >= u1):
            break
        b -= m.end()
        why.append("3':" + ("Ns" if m.group(1)[0] == "N" else "homopolymer"))
    return a, length - b, why, inside




AUDIT_RECORD = versions.major("rrna_audit")  # what makes a read's cached audit (audit_text): versions.STEPS
MEDIUM = {"flank", "chimera_cut", "chimera_weak", "chimera_recurrent"}  # the audit's flags of tier medium; else low
PUBLISHED = ("EUKARYOME", "ROD")  # the sources audited
_VECTORS = {}  # {UniVec folder: Vectors}, built once a process


def vectors_of(folder):
    """A UniVec folder's Vectors, built once a process (published reads are audited in pools of processes)."""
    if folder and folder not in _VECTORS:
        from . import sources  # noqa: PLC0415
        _VECTORS[folder] = Vectors(list(sources.univec(folder)))
    return _VECTORS.get(folder)


def audit_reads(items, vectors=None, primers=None, nhmmer=None, hmm=None, threads=1, tmp=None):
    """The ends' and flanks' audit of reads [(name, sequence on its rRNA's strand, regions [(name, start, end)] or
    None)] -> {name: (start, end, [why], [flags], (grow5, grow3))}: the read's part kept [start, end) (0-based), each
    cut's why (`5':vector NGB00848.1`, `3':flank_beyond`, ...; none: kept whole), its flags (vector, flank, no_rrna,
    trim_short) and the bases its first and last regions grow by (since 2026-11-20: the gene going on past barrnap's
    end, flank_cut)."""
    seqs = [s for _n, s, _r in items]
    units = [unit_of({"regions_list": r}, len(s)) for _n, s, r in items]
    spans = vectors.hits(seqs) if vectors is not None else {}
    slacks = {j: (u[0], len(seqs[j]) - u[1]) for j, u in enumerate(units)
              if u and (u[0] > PRIMER_SLACK or len(seqs[j]) - u[1] > PRIMER_SLACK)}
    ends = primer_ends(seqs, primers, slacks) if primers else {}
    out = {}
    for j, (name, seq, _r) in enumerate(items):
        unit, flags = units[j], []
        if unit is None:
            flags.append("no_rrna")
        a, cut_b, why, inside = end_cuts(seq, unit, spans.get(j, ()), ends.get(j, {}))
        if inside:
            flags.append("vector")
        out[name] = [a, len(seq) - cut_b, why, flags, (0, 0)]
    # the flanks: more than FLANK bases beside the rRNA after the ends' cuts, the read searched again whole
    flanked = [j for j, (name, seq, _r) in enumerate(items) if units[j] and
               (units[j][0] - out[name][0] > FLANK or out[name][1] - units[j][1] > FLANK)]
    hits = nhmmer_hits([items[j][:2] for j in flanked], nhmmer, hmm, threads, tmp) if flanked and nhmmer and hmm \
        else {}
    lengths = model_lengths(hmm) if hmm and os.path.isfile(hmm) else {}
    for j in flanked:
        name, seq, regions = items[j]
        a, b, why, flags, grow = out[name]
        # after an operon's partial 28S the read's bases are its 28S going on (EUKARYOME's reads end at a 28S
        # primer's site, the model aligns its divergent stretches weakly): kept, not flagged
        lsu = any(n == "28S" for n, _x, _y in regions or ())
        if not (nhmmer and hmm):
            flags.append("flank")
            continue
        for side in (5, 3):
            sub_hits = [(m, f, t, s - a, e - a, st) for m, f, t, s, e, st in hits.get(name, ()) if e > a and s < b]
            n, kind, more = flank_cut(b - a, (units[j][0] - a, units[j][1] - a), sub_hits, lengths, side)
            grow = (more, grow[1]) if side == 5 else (grow[0], more)
            if n > 0:
                if side == 5:
                    a += n
                else:
                    b -= n
                why.append(f"{side}':{kind}")
            elif kind == "flank" and not (side == 3 and lsu):
                flags.append("flank")
        out[name] = [a, b, why, flags, grow]
    for name, seq, _r in items:
        a, b, why, flags, grow = out[name]
        if why and b - a < MIN_LEFT:  # not cut: too little would stay
            out[name] = [0, len(seq), [], flags + ["trim_short"], grow]
    return {n: tuple(v) for n, v in out.items()}


def audit_text(result):
    """audit_reads' result of a read as the regions cache's field:
    `<AUDIT_RECORD>|start|end|why;..|flag,..|grow5,grow3` (the last since AUDIT_RECORD 2)."""
    a, b, why, flags, grow = result
    return f"{AUDIT_RECORD}|{a}|{b}|{';'.join(why)}|{','.join(sorted(set(flags)))}|{grow[0]},{grow[1]}"


def parse_audit(text):
    """audit_text's field -> audit_reads' result, or None where missing or of another AUDIT_RECORD."""
    f = (text or "").split("|")
    if len(f) != 6 or f[0] != str(AUDIT_RECORD):
        return None
    g5, g3 = f[5].split(",")
    return (int(f[1]), int(f[2]), [w for w in f[3].split(";") if w], [x for x in f[4].split(",") if x],
            (int(g5), int(g3)))


def grown(regions, grow):
    """A read's regions with its first one's start and its last one's end moved by grow (grow5, grow3)."""
    if not regions or grow == (0, 0):
        return regions
    first = min(range(len(regions)), key=lambda i: regions[i][1])
    last = max(range(len(regions)), key=lambda i: regions[i][2])
    out = [list(r) for r in regions]
    out[first][1] -= grow[0]
    out[last][2] += grow[1]
    return [tuple(r) for r in out]


def apply(name, seq, row, result):
    """A record with its audit applied -> (sequence, row, its TRIM_COLUMNS row or None): its regions grown where its
    gene goes on past barrnap's ends (since 2026-11-20), the read cut to its part, its regions shifted, `trimmed`
    written and the flags added to row["flags_extra"]. A 28S cut at an insertion is partial (lsu_complete false)."""
    a, b, why, flags, grow = result
    flags = sorted(set(list(row.get("flags_extra") or ()) + list(flags)))
    row = dict(row, regions_list=grown(row.get("regions_list"), grow))
    if any(w == "3':flank_insertion" for w in why) and any(n == "28S" for n, _x, _y in row["regions_list"] or ()):
        row["lsu_complete"] = False
    if not why or (a, b) == (0, len(seq)):
        return seq, dict(row, flags_extra=flags, trimmed=row.get("trimmed") or "-"), None
    trim = dict(record=name, source=row.get("source", "-"), length_before=len(seq), start=a + 1, end=b,
                why="; ".join(why))
    row = dict(cut_row(row, a, b), flags_extra=flags, trimmed=f"{a + 1}-{b} of {len(seq)}: " + "; ".join(why))
    return seq[a:b], row, trim


def clean(records, vectors=None, primers=None, nhmmer=None, hmm=None, threads=1, tmp=None, log=None, cached=None):
    """The ends' and flanks' audit of the published records [(name, sequence, row)] (a genome's or a curated unit,
    passed as it is) -> (records, trims [TRIM_COLUMNS rows], counts). cached: {name: audit_reads' result} of the reads
    the published batches audited (eukore_build.published: across the array's nodes); the others audited here."""
    cached = cached or {}
    todo = [(n, s, r.get("regions_list")) for n, s, r in records if r.get("source") in PUBLISHED and n not in cached]
    found = dict(cached)
    if todo:
        found.update(audit_reads(todo, vectors, primers, nhmmer, hmm, threads, tmp))
        if log:
            log(f"  audit: {len(todo):,} published reads' ends and flanks audited here "
                f"({len(cached):,} by the published batches)")
    out, trims, counts = [], [], collections.Counter()
    for name, seq, row in records:
        if name not in found:
            out.append((name, seq, row))
            continue
        _a, _b, why, flags, grow = found[name]
        for f in flags:
            counts[f"flag: {f}"] += 1
        if grow != (0, 0):
            counts["rRNA grown past barrnap's end"] += 1
        for w in why:
            counts["cut: " + w.split(":", 1)[1].split()[0].split("+")[0]] += 1
        seq, row, trim = apply(name, seq, row, found[name])
        if trim:
            trims.append(trim)
            counts["reads cut"] += 1
        out.append((name, seq, row))
    if log:
        log("  audit (ends, flanks): " + ", ".join(f"{k} {v:,}" for k, v in sorted(counts.items())))
    return out, trims, dict(counts)


def tools(src_dir, sources_mod):
    """The audit's inputs of a build: (UniVec's folder or None, the rDNA primers, nhmmer, barrnap's eukaryote HMM)."""
    return (sources_mod.ready(src_dir, "UniVec") if src_dir else None, rdna_primers(), shutil.which("nhmmer"),
            barrnap_hmm(shutil.which("barrnap")))


def trim_rows(trims):
    """clean's trims as rrna/audit.tsv's rows (chimeras.AUDIT_COLUMNS): check `ends`, each cut's why its call."""
    return [dict(record=t["record"], check="ends", call=t["why"], action=f"cut to {t['start']}-{t['end']} of "
                 f"{t['length_before']}", evidence="-") for t in trims]


def stamp(univec, primers, nhmmer, hmm, counts, screen_counts):
    """The manifest's inputs.audit: what the audit ran with and what it found."""
    from . import sources  # noqa: PLC0415
    return {"record": AUDIT_RECORD, "univec": sources.SOURCES["UniVec"]["version"] if univec else None,
            "primers": len(primers or ()), "nhmmer": bool(nhmmer), "barrnap_hmm": os.path.basename(hmm) if hmm else None,
            "ends_and_flanks": counts, "chimeras": screen_counts}


def write_trims(path, trims):
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(TRIM_COLUMNS) + "\n")
        for t in trims:
            fh.write("\t".join(str(t[c]) for c in TRIM_COLUMNS) + "\n")
    os.replace(util.part_of(path), path)


# --- the chimera screens ---------------------------------------------------------------------------------------------

KINGDOM, FAMILY, GENUS = 1, 7, 8  # eukaryome.RANKS' places
# the window screen's clade: EUKARYOME's kingdom (Fungi, Metazoa, Alveolata, Amoebozoa, ...), the 18S's match of
# KSGP's 16S phylum: in the 18S, phyla of one kingdom are as alike as a bacterial phylum's families in the 16S
# (EUKARYOME 2.1 screened by phylum: 4,987 "mislabels" and 1,388 chimeras, most of them Endogonomycetes (Mucoromycota)
# reads 97-99% like Mortierellomycota genomes' 18S in windows, a lineage without a genome of its own; 2026-11-14)
CLADE = KINGDOM


def _named(name):
    from . import eukaryome  # noqa: PLC0415
    return name if name and name != "unclassified" and not eukaryome.is_placeholder(name) else None


def lineages_of(records, tree, ncbi):
    """{record: its EUKARYOME lineage, the ten ranks' names ("" unnamed)}: a EUKARYOME record's label's (the label
    the screens test), another record's (a genome's, ROD's) by its NCBI taxon (eukaryome.record_node)."""
    from . import eukaryome  # noqa: PLC0415
    cache, out = {}, {}
    for name, _seq, row in records:
        label = row.get("source_label")
        if row.get("source") == "EUKARYOME" and label not in (None, "-", ""):
            key = ("label", label)
            if key not in cache:
                cache[key] = tree.text(tree.place_label(label))
        else:
            key = ("taxid", int(row["taxid"]))
            if key not in cache:
                cache[key] = tree.text(eukaryome.record_node(tree, ncbi, int(row["taxid"]))[0])
        names = cache[key].split(";") if cache[key] else []
        out[name] = names + [""] * (len(eukaryome.RANKS) - len(names))
    return out


def is_reference(name, row):
    """A genome's 18S: Eukore's own, ROD's, EUKARYOME's GCA_ records."""
    return row.get("source") in ("genome", "ROD") or (row.get("source") == "EUKARYOME" and
                                                     name.startswith("EUK_GCA_"))


def genome_id(name, row):
    """A reference's genome: one id for an assembly's every copy (screen.assembly_of: GCA_x.1, GCF_x.2, EUKARYOME's
    EUK_GCA_x and its _dup records; since 2026-11-15 they were three genomes, and the window screen's "a second genome"
    was met by one assembly), else its accession, else its name."""
    if row.get("source") == "EUKARYOME":
        acc = name[4:] if name.startswith("EUK_") else name
    else:
        acc = row.get("accession") if row.get("accession") not in (None, "-", "") else name
    from .screen import assembly_of  # noqa: PLC0415
    return assembly_of(acc)


# a genome's references in the window screen, at most (since 2026-11-15; the audit's R1): every genome 18S and operon
# was one, 4.43 M in the 2026-11 build, 8,140 of 22,740 genomes holding 31 or more; minimap2 skips a seed seen in over
# 500 references, so a window like a genome of 600 copies had no hit, and at 40 copies a window's 30 secondary hits
# were all one genome's: the second genome the screen asks for never showed
REFS_PER_GENOME = 3


def references(records, ssu, lin):
    """The window screen's references {name: (18S, its EUKARYOME clade, its genome)}: of each genome its distinct 18S,
    the REFS_PER_GENOME longest (then by name)."""
    by_genome = collections.defaultdict(dict)
    for name, _s, row in sorted(records, key=lambda r: r[0]):
        clade = _named(lin[name][CLADE])
        if is_reference(name, row) and clade and ssu[name]:
            by_genome[genome_id(name, row)].setdefault(ssu[name], (name, clade))
    out = {}
    for genome, seqs in by_genome.items():
        for seq, (name, clade) in sorted(seqs.items(), key=lambda kv: (-len(kv[0]), kv[1][0]))[:REFS_PER_GENOME]:
            out[name] = (seq, clade, genome)
    return out


def screen(records, tree, ncbi, minimap2=None, vsearch=None, searcher=None, threads=1, tmp=None, log=None):
    """The chimera screens (chimeras.py) of the published records' 18S -> (records, audit rows (chimeras.AUDIT_COLUMNS),
    counts). A chimera of an 18S record cut to its own part (`chimera_cut`); else the flags of chimeras.py's calls in
    row["flags_extra"]; each call a row of rrna/audit.tsv."""
    from . import chimeras  # noqa: PLC0415
    lin = lineages_of(records, tree, ncbi)
    ssu = {name: rrna.ssu_part(seq, row) for name, seq, row in records}
    screened = {name for name, _s, row in records if row.get("source") == "EUKARYOME" and not is_reference(name, row)}
    rows, counts = [], collections.Counter()
    found = {}  # name -> [(check, call, keep, evidence)]
    if minimap2:
        refs = references(records, ssu, lin)
        reps, members = {}, collections.defaultdict(list)
        for name in sorted(screened):
            key = (ssu[name], _named(lin[name][CLADE]))
            rep = reps.setdefault(key, f"w{len(reps)}")
            members[rep].append(name)
        calls = chimeras.window_screen({rep: key for key, rep in reps.items()}, refs, minimap2, threads, tmp,
                                       searcher, log)
        for rep, (call, keep, evidence) in calls.items():
            for name in members[rep]:
                found.setdefault(name, []).append(("window", call, keep, evidence))
    elif log:
        log("WARNING: no minimap2: the 18S not screened for pieces of another phylum (KSGP 4.5's window screen)")
    if vsearch:
        calls = chimeras.uchime_screen({n: ssu[n] for n in screened}, {n: ssu[n] for n, _s, r in records
                                                                      if n not in screened},
                                       {n: (_named(lin[n][GENUS]), _named(lin[n][FAMILY])) for n in lin},
                                       searcher, vsearch, threads, tmp, log)
        for name, (call, keep, evidence) in calls.items():
            found.setdefault(name, []).append(("uchime", call, keep, evidence))
    flag_of = {"chimera": "chimera", "mislabel": "mislabel", "weak": "chimera_weak", "recurrent": "chimera_recurrent"}
    out = []
    for name, seq, row in records:
        if name not in found:
            out.append((name, seq, row))
            continue
        flags = list(row.get("flags_extra") or ())
        for check, call, keep, evidence in found[name]:
            action = "-"
            if call == "chimera" and keep and row.get("kind") == "18S" and "chimera_cut" not in flags:
                a, b = keep
                seq, row = seq[a:b], cut_row(row, a, b)
                cut = f"{check} chimera: kept {a + 1}-{b}"
                row["trimmed"] = cut if row.get("trimmed") in (None, "-") else f"{row['trimmed']}; then {cut}"
                flags.append("chimera_cut")
                action = f"cut to {a + 1}-{b}"
            elif call in flag_of:
                flags.append(flag_of[call])
                action = f"flagged {flag_of[call]}"
            counts[f"{check}: {call}" + (" (own part kept)" if action.startswith("cut") else "")] += 1
            rows.append(dict(record=name, check=check, call=call, action=action, evidence=evidence))
        out.append((name, seq, dict(row, flags_extra=sorted(set(flags)))))
    if log:
        log("  chimera screens: " + (", ".join(f"{k} {v:,}" for k, v in sorted(counts.items())) or "none called"))
    return out, rows, dict(counts)
