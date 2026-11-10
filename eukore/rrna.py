"""Eukore's rRNA (since 2026-11-05; from protal's genome_18s.py and rdna_freeze.py, made complete and KSGP-free): every
18S span and complete rDNA operon of the build's genomes and of the published long reads (EUKARYOME 2.1, ROD v1.2),
labelled by NCBI taxid, checked without any KSGP.

- A genome's units: windows around conserved 18S sites (SEEDS) reaching DOWN bases past them on the 18S's strand;
  barrnap (--kingdom euk) finds the 18S, 5.8S, 28S and 5S; an 18S cut by introns is one span. A unit is the 18S span's
  start to its 28S's end, uncut (protal cut its records at 8,191 bases), with its 5.8S; ITSx gives the ITS1/5.8S/ITS2
  borders and the 28S's start where installed. A 5S in the same block (within GAP_5S of the 28S's end or of the 18S's
  start, on either strand, as in some fungi's and protists' repeats) extends the unit and is one of its regions; a 5S
  elsewhere in the genome is not taken. `lsu_complete`: barrnap's 28S hit is not partial.
- Published records: ROD's operon variants by their own NCBI taxid; EUKARYOME's long reads by their species name
  matched to NCBI's names and synonyms (else their genus's), its compromised-quality ids left out. Their regions as a
  genome's units'.
- Checks, each a column (a record failing one is kept and flagged): ambiguous bases, the longest homopolymer, the 99%
  18S type within its species, and own_species: a genome's 18S types against every other source's 18S (vsearch),
  agreeing to the family or deeper with its NCBI lineage (t), sharing no more than the phylum at FOREIGN_ID or more
  (f, a contaminant's rDNA), else undecided (-).
"""

import collections
import os
import re
import subprocess
import tempfile

import numpy as np

from . import util
from .copies import revcomp  # IUPAC's codes too (the copies' and the records' one table)

MOTIFS = ("GTGCCAGCAGCCGCGGTAA", "AAACTTAAAGGAATTGACGG", "GGTGGTGCATGGCCGTTCTTAGTT")  # 18S ~570, ~1,140, ~1,200
SEED = 12
UP = 3000           # a window's bases before a seed on its 18S's strand
DOWN = 15000        # and after it: the whole 28S (human's ends ~9.5 kb past the 18S's start) and an inline 5S
MAX_INTRON = 1500   # an 18S cut by introns: pieces this close are one span
MIN_LEN = 1200      # an 18S span this long at least
MAX_ITS = 3000      # the 18S's end to the 5.8S, the 5.8S's end to the 28S: at most
GAP_5S = 3000       # an inline 5S: this close to the unit's ends
SAME_TYPE = 0.99
NEAR = 0.5
FOREIGN_ID = 97.0
PREFIX = 1000       # a published read's first bases: its 18S, for the own-species pool
ITSX_FIELD = re.compile(r"(SSU|ITS1|5\.8S|ITS2|LSU): (\d+)-(\d+)")
# ITSx 1.1.3 without both of the 5.8S's ends still writes ITS1 and ITS2, their inner borders guessed (an ITS2 45
# bases after the 5.8S start where it found no end: Arabidopsis with set T, Blastocystis, 2026-11-06)
ITSX_58_PROBLEM = re.compile(r"5\.8S: (No start|No end|Not found|Overlap SSU|Overlap LSU)")
RECORD_COLUMNS = ("record", "kind", "taxid", "species_taxid", "accession", "source", "regions", "length",
                  "lsu_complete", "ambiguous", "homopolymer", "type", "own_species", "screen", "label_from",
                  "contig", "start", "end", "strand", "contig_edge", "introns", "gaps", "confidence", "flags",
                  # since 2026-11-09 (format 6): the other clade's species its 18S is shared with (qc.verdict); the
                  # SSU model its 18S scores best by (regions.ssu_models: eukarya, else bacteria or archaea: not a
                  # nucleus' 18S; - unscored)
                  "cross", "ssu_model")
GAP_MIN = 10        # a run of Ns this long: an assembly gap (or a read's)
EDGE = 1000         # a genome's record this close to its contig's end: maybe cut, or joined to another sequence
MAX_AMBIGUOUS = 0.005
HOMOPOLYMER = 36

SEEDS = sorted({s[i:i + SEED] for m in MOTIFS for s in (m, revcomp(m)) for i in range(len(m) - SEED + 1)})
FORWARD = {m[i:i + SEED] for m in MOTIFS for i in range(len(m) - SEED + 1)}


class BarrnapError(RuntimeError):
    pass


class ItsxError(RuntimeError):
    pass


def check_itsx(itsx, threads=2):
    """Why ITSx cannot be used (it fails on a short sequence, every profile set searched), or None."""
    try:
        run_itsx(["ACGTTGCA" * 40], itsx, threads)
    except (ItsxError, OSError) as e:
        return str(e)
    return None


def barrnap_gff(barrnap, fasta, threads=1, kingdom="euk"):
    command = [barrnap, "--kingdom", kingdom, "--threads", str(threads), "--quiet", fasta]
    run = subprocess.run(command, capture_output=True, text=True, check=False)
    if run.returncode != 0:
        raise BarrnapError(f"barrnap failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
    return run.stdout


# an organelle molecule's FASTA title (organelle.py's too, and metadata.py's: both read its word as the organelle,
# mitochondrion or else a plastid, so that a word of another organelle is not added here: NOT_NUCLEAR's are)
ORGANELLE_TITLE = re.compile(r"\b(mitochondrion|chloroplast|plastid|apicoplast)\b", re.IGNORECASE)
# not the nucleus either: a cryptophyte's or chlorarachniophyte's nucleomorph (Guillardia theta, Bigelowiella natans),
# the endosymbiont's red or green algal nucleus with rDNA operons of its own (until 2026-11-06 labelled the host's);
# since 2026-11-10 the organelles' other titles too (until then read as the nucleus', their SSU an 18S): a
# kinetoplastid's mitochondrial DNA (kinetoplast, maxicircle, minicircle), Paulinella's chromatophore, a glaucophyte's
# cyanelle, and DDBJ's "mitochondrial DNA, complete genome". "mitochondrial" alone is not taken: a title may name a
# nuclear gene's product by it (a clone "... gene for mitochondrial ribosomal protein", "mitochondrial DNA polymerase
# gamma"), so only "mitochondrial DNA" or "mitochondrial genome" where the title's part ends
NOT_NUCLEAR = re.compile(r"\b(mitochondrion|chloroplast|plastid|apicoplast|nucleomorph|kinetoplast|maxicircle|"
                         r"minicircle|chromatophore|cyanelle)\b|\bmitochondrial\s+(?:DNA|genome)\s*(?:[,;]|$)",
                         re.IGNORECASE)


_CODE = np.full(256, 4, dtype=np.uint8)  # a base's 2 bits (either case); 4: not ACGT, no seed across it
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _CODE[_b + 32] = _i
_SEED_KIND = None  # every 12-mer's code -> 0 no seed, 1 a forward seed, 2 a reverse one (built once a process)
SEED_BATCH = 1 << 24  # bases scanned at a time


def _seed_kinds():
    global _SEED_KIND
    if _SEED_KIND is None:
        table = np.zeros(1 << (2 * SEED), dtype=np.uint8)
        for s in SEEDS:
            code = 0
            for c in s:
                code = code * 4 + "ACGT".index(c)
            table[code] = 1 if s in FORWARD else 2
        _SEED_KIND = table
    return _SEED_KIND


def seed_hits(data):
    """The seeds in a stretch of bases (bytes, no line ends) -> (positions of their first base, kinds 1 forward, 2
    reverse), of every 12-mer wholly inside it. One numpy pass (since 2026-11-08: a str.find per seed, 60 passes over
    each contig, ~66 s a Gb)."""
    c = _CODE[np.frombuffer(bytes(data), dtype=np.uint8)]
    n = len(c) - SEED + 1
    if n <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.uint8)
    bits = (c & 3).astype(np.uint32)
    k = bits[:n].copy()
    for j in range(1, SEED):
        k <<= 2
        k |= bits[j:j + n]
    bad = np.concatenate(([0], np.cumsum(c > 3, dtype=np.int64)))
    kind = _seed_kinds()[k]
    kind[bad[SEED:SEED + n] != bad[:n]] = 0
    pos = np.nonzero(kind)[0]
    return pos, kind[pos]


def _spans_of(hits, length):
    """Seed hits [(position, kind)] on a sequence of `length` -> merged [start, end) windows around them."""
    sites = sorted((max(0, i - (UP if k == 1 else DOWN)), min(length, i + (DOWN if k == 1 else UP)))
                   for i, k in ((int(p), int(k)) for p, k in hits))
    spans = []
    for a, b in sites:
        if spans and a <= spans[-1][1]:
            spans[-1][1] = max(spans[-1][1], b)
        else:
            spans.append([a, b])
    return spans


def windows(seq):
    """Merged [start, end) windows around the seeds in a sequence."""
    pos, kind = seed_hits(seq.encode("ascii", "replace") if isinstance(seq, str) else seq)
    return _spans_of(zip(pos.tolist(), kind.tolist()), len(seq))


class SeedScan:
    """The seeds of a genome's sequences fed in turn (start, feed, finish), scanned in batches of SEED_BATCH bases
    across sequences (an N between two: no seed across them), so that a fragmented assembly's millions of contigs
    are not millions of numpy calls -> finish(): {sequence: merged [start, end) windows}, and lengths {sequence: its
    length}."""

    def __init__(self):
        self.buf, self.marks, self.hits, self.lengths = bytearray(), [], collections.defaultdict(list), {}
        self.name = None

    def start(self, name):
        self.name = name
        self.lengths[name] = 0
        self.buf += b"N"
        self.marks.append((len(self.buf), name, 0))

    def feed(self, data):
        self.buf += data
        self.lengths[self.name] += len(data)
        if len(self.buf) >= SEED_BATCH:
            self._scan(keep=True)

    def _scan(self, keep=False):
        pos, kind = seed_hits(self.buf)
        if len(pos):
            at = np.array([m[0] for m in self.marks])
            which = np.searchsorted(at, pos, side="right") - 1
            for p, k, w in zip(pos.tolist(), kind.tolist(), which.tolist()):
                start, name, offset = self.marks[w]
                self.hits[name].append((offset + p - start, k))
        if keep and self.name is not None:  # the current sequence's last bases: a seed may start in them
            start, name, offset = self.marks[-1]
            carry = min(SEED - 1, len(self.buf) - start)
            self.buf = self.buf[len(self.buf) - carry:]
            self.marks = [(0, name, self.lengths[name] - carry)]
        else:
            self.buf, self.marks = bytearray(), []

    def finish(self):
        self._scan()
        return {n: _spans_of(h, self.lengths[n]) for n, h in self.hits.items()}, self.lengths


ALIGNED = re.compile(r"aligned only (\d+) percent")
WHOLE = 0.8         # barrnap's own: a hit of this share of its model or more is not partial


def parse_gff(gff):
    """barrnap's GFF -> {window: {feature name: [(start, end, strand, partial, share of the model aligned)]}} (1-based
    inclusive)."""
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in gff.splitlines():
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        name = f[8].split("Name=", 1)[1].split(";", 1)[0] if "Name=" in f[8] else ""
        partial = "aligned only" in f[8] or "partial" in f[8]
        share = ALIGNED.search(f[8])
        out[f[0]][name].append((int(f[3]), int(f[4]), f[6], partial,
                                int(share.group(1)) / 100 if share else 0.5 if partial else 1.0))
    return out


def joined(features, gap=MAX_INTRON):
    """Features (parse_gff's) on a strand this close one span -> [[start, end, strand, share of the model]]: an rRNA cut
    by introns (barrnap's pieces each partial; their shares summed)."""
    out = []
    for strand in "+-":
        for f in sorted(f for f in features if f[2] == strand):
            start, end = f[0], f[1]
            share = f[4] if len(f) > 4 else 0.5 if f[3] else 1.0  # a feature of (start, end, strand, partial)
            if out and out[-1][2] == strand and start - out[-1][1] <= gap:
                out[-1][1] = max(out[-1][1], end)
                out[-1][3] += share
            else:
                out.append([start, end, strand, share])
    return out


def units_of(seq, feats):
    """A window's (or a read's) 18S spans and their units -> [dict(span (start, end, strand), sequence (the span,
    its strand), exons (the 18S pieces joined), unit (sequence or None), regions [(name, start, end)] on the unit,
    lsu_complete, its stretch for ITSx and where that starts, unit_span (a, b) in `seq`)]. Coordinates 1-based in
    `seq`. A 28S cut by introns (barrnap's pieces, each partial) is one, whole where its pieces' shares make WHOLE
    (until 2026-11-06 the unit ended at its first piece). A 5S is one unit's at most (until 2026-11-06 a tandem
    array's spacer 5S was both its neighbours'): the unit whose 28S it follows within GAP_5S, else whose 18S it
    precedes; a unit keeps the nearest of its own."""
    out, cores = [], []
    s18 = feats.get("18S_rRNA", [])
    lsu_spans = joined(feats.get("28S_rRNA", []))
    for start, end, strand, _share in joined(s18):
        pieces = [f for f in s18 if f[2] == strand and start <= f[0] and f[1] <= end]
        span = seq[start - 1:end]
        exons = "".join(seq[f[0] - 1:f[1]] for f in sorted(pieces))
        if strand == "-":
            span, exons = revcomp(span), revcomp(exons)
        # the 18S's introns: the gaps between its pieces, on the span ([start, end) on its strand)
        ordered = sorted(pieces)
        cuts = [(x[1] + 1, y[0] - 1) for x, y in zip(ordered, ordered[1:]) if y[0] > x[1] + 1]
        on_span = sorted((a - start, b - start + 1) if strand == "+" else (end - b, end - a + 1) for a, b in cuts)
        item = dict(span=(start, end, strand), sequence=span, exons=exons, unit=None, regions=None,
                    lsu_complete=False, introns=[("18S",) + x for x in on_span])
        out.append(item)
        fwd = strand == "+"
        s58 = [f for f in feats.get("5_8S_rRNA", []) if f[2] == strand and
               (0 < f[0] - end <= MAX_ITS if fwd else 0 < start - f[1] <= MAX_ITS)]
        if not s58:
            continue
        f58 = min(s58, key=lambda f: f[0]) if fwd else max(s58, key=lambda f: f[1])
        lsu = [f for f in lsu_spans if f[2] == strand and
               ((f[1] > f58[1] + 100 and f[0] <= f58[1] + MAX_ITS) if fwd else
                (f[0] < f58[0] - 100 and f[1] >= f58[0] - MAX_ITS))]
        if lsu:
            f28 = min(lsu, key=lambda f: f[0]) if fwd else max(lsu, key=lambda f: f[1])
            cores.append((item, start, end, fwd, f58, f28, cuts) + ((start, f28[1]) if fwd else (f28[0], end)))
    five_of = {}  # unit -> ((rank, gap), its 5S)
    for f in feats.get("5S_rRNA", []):  # inline: in the same block, either strand
        best = None
        for i, (_item, _s, _e, fwd, _f58, _f28, _cuts, a, b) in enumerate(cores):
            after, before = (f[0] - b, a - f[1]) if fwd else (a - f[1], f[0] - b)  # past the 28S, before the 18S
            for rank, gap in ((0, after), (1, before)):
                if 0 <= gap <= GAP_5S and (best is None or (rank, gap) < best[0]):
                    best = ((rank, gap), i)
        if best and (best[1] not in five_of or best[0][1] < five_of[best[1]][0][1]):
            five_of[best[1]] = (best[0], f)
    for i, (item, start, end, fwd, f58, f28, cuts, a, b) in enumerate(cores):
        five = [five_of[i][1]] if i in five_of else []
        a = min([a] + [f[0] for f in five])
        b = max([b] + [f[1] for f in five])
        unit = seq[a - 1:b]
        if not fwd:
            unit = revcomp(unit)

        def at(x1, x2, a=a, b=b, fwd=fwd):  # a feature's [start, end) on the unit
            return (x1 - a, x2 - a + 1) if fwd else (b - x2, b - x1 + 1)

        r18, r58 = ("18S",) + at(start, end), ("5.8S",) + at(f58[0], f58[1])
        r28 = ("28S",) + at(f28[0], f28[1])
        # barrnap's 28S model can reach back over the 5.8S (Drosophila's units): the 28S starts after it
        r28 = ("28S", max(r28[1], r58[2]), r28[2])
        regions = [r18, ("ITS1", r18[2], r58[1]), r58, r28] + [("5S",) + at(f[0], f[1]) for f in five]
        if r28[1] > r58[2]:
            regions.insert(3, ("ITS2", r58[2], r28[1]))
        item.update(unit=unit, regions=sorted(regions, key=lambda r: r[1]), lsu_complete=f28[3] >= WHOLE,
                    unit_span=(a, b), unit_introns=sorted(("18S",) + at(x1, x2) for x1, x2 in cuts),
                    its=unit[max(0, r18[2] - 300):min(len(unit), r58[2] + 3300)],
                    its_at=max(0, r18[2] - 300))
    return out


# ITSx 1.1.3's profile sets by the NCBI taxon each stands for (checked against NCBI's taxdump of 2026-11-05). ITSx
# searches every set with -t all, in turn and without HMMER's filters (--max): ~50 core-minutes per 2,000 published
# reads (the 2026-11-05 trial). A record gets its lineage's sets; one in none of them (Blastocystis, Cryptista, ...)
# all of them.
ITSX_SETS = {33630: "A", 3208: "B", 2836: "C", 554915: "D", 33682: "E", 4751: "F", 3041: "G", 2763: "H", 2870: "I",
             3195: "L", 33208: "M", 6029: "N", 4762: "O", 2830: "P", 38410: "Q", 543769: "R", 33859: "S", 58023: "T",
             5747: "U", 554296: "X", 5719: "Y"}


def itsx_groups(ncbi, taxid):
    """ITSx's -t for a taxon: the profile sets of its lineage (comma-separated), else "all"."""
    if taxid in (None, "", "-"):
        return "all"
    taxid = int(taxid)
    sets = sorted({ITSX_SETS[a] for a in ncbi.ancestors(ncbi.merged.get(taxid, taxid)) if a in ITSX_SETS})
    return ",".join(sets) or "all"


def run_itsx(stretches, itsx, threads=1, groups=None):
    """{stretch: ITSx's positions {region: (begin, end) 1-based}} for each distinct stretch. groups: {stretch: ITSx's
    -t (itsx_groups)}, one ITSx run per distinct -t; none given: all sets."""
    if not stretches or not itsx:
        return {}
    groups = groups or {}
    by_set = collections.defaultdict(list)
    for s in dict.fromkeys(stretches):
        by_set[groups.get(s, "all")].append(s)
    out = {}
    for sets, part in sorted(by_set.items()):
        out.update(_itsx(part, itsx, threads, sets))
    return out


def _itsx(stretches, itsx, threads, sets):
    with tempfile.TemporaryDirectory() as tmp:
        fasta, prefix = os.path.join(tmp, "its.fasta"), os.path.join(tmp, "itsx")
        names = {}
        with open(fasta, "w") as fh:
            for s in stretches:
                if s not in names:
                    names[s] = f"u{len(names)}"
                    fh.write(f">{names[s]}\n{s}\n")
        run = subprocess.run([itsx, "-i", fasta, "-o", prefix, "-t", sets, "--cpu", str(max(1, threads)),
                              "--save_regions", "none", "--preserve", "T", "--complement", "F", "--graphical", "F",
                              "--summary", "F"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        # ITSx with --cpu over 1 exits 0 when HMMER fails on a profile set (a damaged ITSx_db): its stderr tells
        if run.returncode or "ERROR" in run.stderr:
            lines = [x for x in run.stderr.splitlines() if "rror" in x or "not found" in x]
            raise ItsxError(f"ITSx failed (exit {run.returncode}): {' | '.join(lines[:4]) or run.stderr[-300:]}")
        found = {}
        if os.path.isfile(prefix + ".positions.txt"):
            with open(prefix + ".positions.txt") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    got = {m.group(1): (int(m.group(2)), int(m.group(3)))
                           for m in (ITSX_FIELD.match(x) for x in f[2:]) if m}
                    problem = [m.group(1) for m in (ITSX_58_PROBLEM.match(x) for x in f[2:]) if m]
                    if problem:
                        got["5.8S problem"] = problem[0]
                    found[f[0]] = got
        return {s: found.get(n, {}) for s, n in names.items()}


def refine(item, found):
    """A unit's ITS1/5.8S/ITS2/28S borders from ITSx where it found ITS1 and ITS2 (barrnap's otherwise); the 5.8S's
    start (ITS1's end) and end (ITS2's start) ITSx's only where it found them (ITSX_58_PROBLEM), barrnap's otherwise."""
    if not item.get("unit") or not all(k in found for k in ("ITS1", "ITS2")):
        return item
    regions = {r[0]: r for r in item["regions"]}
    problem = found.get("5.8S problem")
    at = item["its_at"]
    its1 = (at + found["ITS1"][0] - 1, at + found["ITS1"][1] if problem in (None, "No end") else regions["5.8S"][1])
    its2 = (at + found["ITS2"][0] - 1 if problem in (None, "No start") else regions["5.8S"][2], at + found["ITS2"][1])
    if not its1[0] <= its1[1] < its2[0] <= its2[1]:
        return item
    end28 = regions["28S"][2]
    new = [r for r in item["regions"] if r[0] not in ("18S", "ITS1", "5.8S", "ITS2", "28S")]
    new += [("18S", regions["18S"][1], its1[0]), ("ITS1",) + its1, ("5.8S", its1[1], its2[0]), ("ITS2",) + its2,
            ("28S", its2[1], end28)]
    item["regions"] = sorted(new, key=lambda r: r[1])
    return item


def regions_text(regions):
    return ",".join(f"{n}:{a}-{b}" for n, a, b in regions) if regions else "-"


def extract_genome(path, barrnap, itsx=None, threads=1, itsx_sets="all"):
    """A genome's rRNA -> [dict(kind 18S or rDNA, contig, start, end, strand, sequence, exons, regions, lsu_complete)]:
    each 18S span (MIN_LEN at least), and each unit with its 5.8S and 28S as an rDNA record, uncut. itsx_sets: ITSx's
    -t for the genome's taxon (itsx_groups). Organelle molecules are left out: barrnap's eukaryote models call a
    mitochondrial SSU an 18S (Arabidopsis' BK010421.1, 2026-11-06), whose bacterial kind then takes bacterial reads;
    organelle.py takes their rRNA."""
    from . import metadata, placement  # noqa: PLC0415
    # two passes, neither holding a sequence whole (since 2026-11-08: each contig was held 2-3 times, a 5 Gb
    # chromosome ~15 GB): the seeds (SeedScan, over igzip's output where there), then the windows' bases alone
    scan, keep = SeedScan(), False
    for title, data in metadata._pieces_of(path):  # noqa: SLF001 (the package's own FASTA reader)
        if title is not None:
            text = title.decode("utf-8", "replace")
            keep = not NOT_NUCLEAR.search(text)
            if keep:
                scan.start((text.split() or [""])[0])
        elif keep:
            scan.feed(data.translate(None, b"\n\r \t"))
    found, lengths = scan.finish()
    spans = placement.read_spans(path, {n: [(a + 1, b) for a, b in w] for n, w in found.items()}) if found else {}
    seqs, origin = {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        win = os.path.join(tmp, "windows.fna")
        with open(win, "w") as fh:
            order = {n: i for i, n in enumerate(lengths)}
            for name in sorted(found, key=order.get):
                for a, b in found[name]:
                    wid = f"w{len(seqs)}"
                    seqs[wid], origin[wid] = spans[name][a:b].upper(), (name, a)
                    fh.write(f">{wid}\n{seqs[wid]}\n")
        gff = barrnap_gff(barrnap, win, threads) if seqs else ""
    feats = parse_gff(gff)
    items = []
    for wid, f in feats.items():
        for item in units_of(seqs[wid], f):
            contig, offset = origin[wid]
            item.update(contig=contig, start=offset + item["span"][0], end=offset + item["span"][1],
                        strand=item["span"][2], contig_length=lengths.get(contig))
            if item.get("unit_span"):
                item["unit_at"] = (offset + item["unit_span"][0], offset + item["unit_span"][1])
            items.append(item)
    stretches = [i["its"] for i in items if i.get("unit")]
    found = run_itsx(stretches, itsx, threads, dict.fromkeys(stretches, itsx_sets))
    rows = []
    for item in items:
        if len(item["sequence"]) >= MIN_LEN:
            rows.append(dict(kind="18S", contig=item["contig"], start=item["start"], end=item["end"],
                             strand=item["strand"], sequence=item["sequence"], exons=item["exons"], regions=None,
                             lsu_complete=False, contig_length=item["contig_length"], introns=item["introns"]))
        if item.get("unit"):
            refine(item, found.get(item["its"], {}))
            rows.append(dict(kind="rDNA", contig=item["contig"], start=item["unit_at"][0], end=item["unit_at"][1],
                             strand=item["strand"], sequence=item["unit"], exons=item["exons"],
                             regions=item["regions"], lsu_complete=item["lsu_complete"],
                             contig_length=item["contig_length"], introns=item["unit_introns"]))
    return rows


def published_regions(reads, barrnap, itsx=None, threads=1, sets=None):
    """Published long reads [(name, sequence)] -> {name: (regions, lsu_complete, 18S exons, strand)}: each read as a
    window, its regions on the record (the read itself on +, its reverse complement on -: on_record). sets: {name:
    ITSx's -t for its taxon (itsx_groups)}, else all sets."""
    if not reads:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "reads.fna")
        with open(fasta, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in reads)
        gff = barrnap_gff(barrnap, fasta, threads)
    feats = parse_gff(gff)
    seqs = dict(reads)
    items = {}
    for name, f in feats.items():
        units = [u for u in units_of(seqs[name], f) if u.get("unit")]
        if units:
            items[name] = max(units, key=lambda u: len(u["unit"]))
    sets = sets or {}
    found = run_itsx([u["its"] for u in items.values()], itsx, threads,
                     {u["its"]: sets.get(n, "all") for n, u in items.items()})
    return {n: (on_record(refine(u, found.get(u["its"], {})), len(seqs[n])), u["lsu_complete"], u["exons"],
                u["span"][2]) for n, u in items.items()}


def on_record(unit, length):
    """A read's unit's regions (on the unit: units_of) on its record, the read on the unit's strand (its reverse
    complement on -) of `length` bases (until 2026-11-06 the regions were the unit's, the record the read: a read whose
    18S began at base 150 had them 149 bases early, one on - mirrored)."""
    a, b = unit["unit_span"]
    shift = a - 1 if unit["span"][2] == "+" else length - b
    return [(name, s + shift, e + shift) for name, s, e in unit["regions"]]


# --- labels by NCBI taxid ----------------------------------------------------------------------------------------

def name_index(ncbi):
    """{lower-case name: taxid} of NCBI's scientific names and synonyms. A name several taxa share is the one's among
    them that holds all the others (since 2026-11-10: NCBI names a genus and its nominate subgenus alike, Aedes, Culex,
    Anopheles and ~160 more, a phylum and its class; until then every shared name was dropped, and EUKARYOME's records
    of such a genus had no label). Unrelated taxa's (homonyms: Morus the plant and the bird): none, as before."""
    index, shared = {}, {}
    for taxid, name in ncbi.name.items():
        for n in [name] + ncbi.synonyms.get(taxid, []):
            key = n.lower()
            if index.setdefault(key, taxid) != taxid:
                shared.setdefault(key, {index[key]}).add(taxid)
    for key, taxids in shared.items():
        holder = set.intersection(*(set(ncbi.ancestors(t)) for t in taxids)) & taxids  # one at most
        if holder:
            index[key] = holder.pop()
        else:
            del index[key]
    return index


def species_name(genus, species):
    """A EUKARYOME species' whole name: its genus and its s__ epithet ("" for unclassified); an s__ that already starts
    with its genus as it is."""
    species = species.replace("_", " ").strip()
    genus = re.sub(r"\(.*?\)", "", genus).strip()  # EUKARYOME's homonyms: Drosophila(Metazoa)
    if not species or species.lower() in ("unclassified", "unidentified") or species.endswith(" sp."):
        return ""
    if genus and not species.lower().startswith(genus.lower() + " ") and " " not in species:
        return f"{genus} {species}"
    return species


def eukaryome_taxid(label, index, ncbi):
    """A EUKARYOME label (k__...;g__Genus;s__epithet) -> (taxid, how) by its species name, else its genus's;
    (None, None) where neither is NCBI's. EUKARYOME 2.1's s__ is the epithet alone (g__Geotrupes;s__spiniger), or
    unclassified: until format 4 it was read as the whole name, and every record was labelled by its genus."""
    parts = {p[:1]: p[3:].replace("_", " ") for p in label.split(";") if "__" in p}
    genus = parts.get("g", "")
    sp = species_name(genus, parts.get("s", ""))
    if sp and (t := index.get(sp.lower())) and ncbi.rank.get(ncbi.species_of(t)) == "species":
        return ncbi.species_of(t), "eukaryome_name"
    if genus and (t := index.get(genus.lower())) and ncbi.rank.get(t) == "genus":
        return t, "eukaryome_genus"
    return None, None


# --- checks ------------------------------------------------------------------------------------------------------

def homopolymer(seq):
    """A sequence's longest run of one character (by numpy since 2026-11-08: ~40x a regex's, ~0.08 s a Mb, minutes
    on one core over every record)."""
    if not seq:
        return 0
    a = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    edges = np.concatenate(([0], np.flatnonzero(a[1:] != a[:-1]) + 1, [len(a)]))
    return int(np.diff(edges).max())


def ambiguous(seq):
    """A sequence's share of characters other than ACGT."""
    return (len(seq) - sum(seq.count(c) for c in "ACGT")) / max(1, len(seq))


def identity(a, b):
    """1 - edits / the shorter's length, the shorter aligned to its best stretch of the longer (Myers' bit-vector)."""
    if len(a) > len(b):
        a, b = b, a
    if not a:
        return 0.0
    m = len(a)
    peq = {}
    for i, c in enumerate(a):
        peq[c] = peq.get(c, 0) | (1 << i)
    mask, high = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score, best = mask, 0, m, m
    for c in b:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (mask & ~(xh | pv))
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph, mh = (ph << 1) & mask, (mh << 1) & mask
        pv, mv = mh | (mask & ~(xv | ph)), ph & xv
        best = min(best, score)
    return 1 - best / m


def types_of(exons):
    """Sequences -> {sequence: type number}, each the first type whose first sequence it is SAME_TYPE identical to,
    the commonest first."""
    reps, out = [], {}
    for seq, _n in sorted(collections.Counter(exons).items(), key=lambda kv: (-kv[1], -len(kv[0]), kv[0])):
        for i, rep in enumerate(reps):
            if identity(seq, rep) >= SAME_TYPE:
                out[seq] = i + 1
                break
        else:
            reps.append(seq)
            out[seq] = len(reps)
    return out


def ssu_part(seq, row):
    """A record's 18S: an 18S record whole, an operon's 18S region, else (a read without one found) its first PREFIX."""
    if row["kind"] == "18S":
        return seq
    return next((seq[a:b] for n, a, b in row.get("regions_list") or () if n == "18S"), seq[:PREFIX])


def own_species(records, ncbi, vsearch, threads=1):
    """The own-species check (module docstring) -> {record: (t, f or -, the other clade's species its 18S is shared
    with or None)}: each genome's 18S types against the 18S of every other genome and source (vsearch
    --usearch_global at 80%, the first 32 hits), their nearest (within NEAR points of the best) agreeing with the
    genome's NCBI lineage to the family or deeper: t; at FOREIGN_ID or more and known to part from it at the class or
    above, its owner told (qc.verdict): f; else -. records: [(name, sequence, row)] (row: accession, source,
    species_taxid, regions); a published record is not checked (its label is its source's). A genome's (or a
    source's) 18S are searched as one per 99% type (screen.representatives), so that its own copies do not fill the
    first hits."""
    from . import qc, screen  # noqa: PLC0415
    if not vsearch or not records:
        return {}
    ssu = {name: ssu_part(seq, row)[:2500] for name, seq, row in records}
    # a record's genome, its RefSeq copy's and every version's one (screen.genome_of; since 2026-11-09: ROD's copy of
    # it was another genome's, a contaminant's operon confirmed by its own copy)
    owner = {name: screen.genome_of(name, row) for name, _s, row in records}
    group = {name: (row["source"], screen.genome_of(name, row)) for name, _s, row in records}
    reps = screen.representatives({n: s for n, s in ssu.items() if s}, group)
    genome_reps = {k for k in reps if group[k][0] == "genome"}
    if not genome_reps:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        found = screen.search({k: ssu[k] for k in sorted(genome_reps)}, {k: ssu[k] for k in sorted(reps)}, vsearch,
                              threads, tmp, ident=0.8, accepts=32, rejects=256)
    taxid_of = {name: int(row["species_taxid"]) for name, _s, row in records}
    cache, out = {}, {}
    # since 2026-11-09 qc.verdict's: another clade's only where known to part at the class or above (apart), and of
    # a sequence two clades' genomes hold, the clade whose other species hold it nearer its owner (until then both
    # genomes' records were f, the true owner's too)
    members = qc.families(ncbi, {taxid_of[k] for k in reps}, cache)
    for rep in genome_reps:
        hits = [(i, taxid_of[t]) for i, t in found.get(rep, []) if owner[t] != owner[rep]]
        verdict, cross, _top, _agree = qc.verdict(taxid_of[rep], hits, ncbi, cache, members, NEAR, FOREIGN_ID)
        for member in reps[rep]:
            out[member] = (verdict, cross)
    return out


def write_records(folder, records):
    """records [(name, sequence, row)] -> folder/records.fna.gz and records.tsv (RECORD_COLUMNS)."""
    os.makedirs(folder, exist_ok=True)
    fasta, table = os.path.join(folder, "records.fna.gz"), os.path.join(folder, "records.tsv")
    with util.bgzf_open(util.part_of(fasta)) as fa, open(util.part_of(table), "w", newline="\n") as tb:
        tb.write("\t".join(RECORD_COLUMNS) + "\n")
        for name, seq, row in records:
            fa.write(f">{name}\n{seq}\n")
            tb.write("\t".join(str(row.get(c, "-")) if row.get(c) not in (None, "") else "-" for c in RECORD_COLUMNS)
                     + "\n")
    os.replace(util.part_of(fasta), fasta)
    os.replace(util.part_of(table), table)


def finish(records, ncbi, vsearch=None, threads=1, log=None, ssu_models=None):
    """The checks over all records [(name, sequence, row)] (row: kind, accession, source, species_taxid, taxid,
    regions_list, lsu_complete, label_from) -> the same with their RECORD_COLUMNS filled, sorted by species and name.
    The 99% types are within a species: records labelled by a genus alone (EUKARYOME's) get none (until 2026-11-06 a
    genus of thousands of environmental reads was typed in pure Python for an hour). log: the screen's counts.
    ssu_models: a function of {record: its 18S} -> {record: the SSU model it scores best by} (regions.ssu_models;
    since 2026-11-09: a bacterium's or an organelle's SSU, barrnap's 18S by its eukaryote models alone, is told)."""
    by_species = collections.defaultdict(list)
    for name, seq, row in records:
        if ncbi is None or ncbi.rank.get(int(row["species_taxid"])) == "species":
            by_species[row["species_taxid"]].append((name, ssu_part(seq, row)))
    type_of = {}
    for sp, items in by_species.items():
        types = types_of([p for _n, p in items])
        for name, part in items:
            type_of[name] = types[part]
    own = own_species(records, ncbi, vsearch, threads)
    from . import screen as screen_mod  # noqa: PLC0415
    screened, counts_ = screen_mod.screen(records, ncbi, vsearch, threads)
    if log and counts_:
        log("rRNA part screen: " + ", ".join(f"{k} {v}" for k, v in sorted(counts_.items())))
    models = {}
    if ssu_models:
        models = ssu_models({name: ssu_part(seq, row) for name, seq, row in records})
        if log:
            log("rRNA SSU models: " + ", ".join(f"{k} {v}" for k, v in sorted(collections.Counter(
                models.values()).items())))
    out = []
    for name, seq, row in records:
        verdict, cross = own.get(name, ("-", None))
        row = dict(row, record=name, length=len(seq), regions=regions_text(row.get("regions_list")),
                   lsu_complete="t" if row.get("lsu_complete") else "f", ambiguous=f"{ambiguous(seq):.4f}",
                   homopolymer=homopolymer(seq), type=type_of.get(name, "-"), own_species=verdict,
                   cross=cross or "-", ssu_model=models.get(name, "-"),
                   screen=screened.get(name, "-"), introns=regions_text(row.get("introns_list")), gaps=gaps(seq))
        edge = contig_edge(row)
        row["contig_edge"] = "-" if edge is None else edge
        row["confidence"], row["flags"] = confidence(row)
        out.append((name, seq, row))
    return sorted(out, key=lambda r: (int(r[2]["species_taxid"]), r[0]))


def genome_records(accession, rows, taxid, species_taxid, source="genome"):
    """A genome's extraction (extract_genome's rows) -> records [(name, sequence, row)]."""
    out, n = [], collections.Counter()
    for r in rows:
        n[r["kind"]] += 1
        name = f"{accession}_{r['kind']}_{n[r['kind']]}"
        out.append((name, r["sequence"], dict(kind=r["kind"], accession=accession, source=source, taxid=taxid,
                                              species_taxid=species_taxid, regions_list=r["regions"],
                                              lsu_complete=r["lsu_complete"], label_from="assembly",
                                              contig=r["contig"], start=r["start"], end=r["end"], strand=r["strand"],
                                              contig_length=r.get("contig_length"),
                                              introns_list=r.get("introns"))))
    return out


SAVED = ("kind", "contig", "start", "end", "strand", "regions", "lsu_complete", "contig_length", "introns", "sequence")


def parse_regions(text):
    """regions_text's "name:start-end,..." -> [(name, start, end)], None for "-"."""
    if text in ("-", "", None):
        return None
    return [(p.rpartition(":")[0], int(p.rpartition(":")[2].split("-")[0]), int(p.rpartition(":")[2].split("-")[1]))
            for p in text.split(",")]


def save(folder, rows):
    """A genome's extraction into its done folder (rrna.tsv: SAVED's columns; contig_length and introns since format
    4) -> its counts."""
    with open(os.path.join(folder, "rrna.tsv"), "w", newline="\n") as fh:
        fh.write("\t".join(SAVED) + "\n")
        for r in rows:
            fh.write(f"{r['kind']}\t{r['contig']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{regions_text(r['regions'])}"
                     f"\t{'t' if r['lsu_complete'] else 'f'}\t{r.get('contig_length') or '-'}"
                     f"\t{regions_text(r.get('introns'))}\t{r['sequence']}\n")
    return counts(rows)


def counts(rows):
    """A genome's rRNA as genomes.tsv counts it: 18S spans, rDNA units (those with a whole 28S), the contigs holding
    units."""
    units = [r for r in rows if r["kind"] == "rDNA"]
    return {"18S": sum(r["kind"] == "18S" for r in rows), "rDNA": len(units),
            "rDNA_complete": sum(1 for r in units if r["lsu_complete"]),
            "rDNA_contigs": len({r["contig"] for r in units})}


COPY_ELEMENTS = ("rDNA_unit", "18S", "5.8S", "28S", "5S")


def element_counts(rows):
    """A genome's rRNA as copy_numbers.tsv counts it (since 2026-11-08) -> {element: (copies, partial)}: rDNA units
    (partial: without a whole 28S), 18S (its spans of MIN_LEN or more, in units or alone, and partial: a unit's 18S
    shorter, at a contig's end: 22 of 1,851 units on the trial of 2026-11-07), the units' 5.8S, 28S (partial: not
    whole) and 5S (an inline 5S alone: 5S arrays elsewhere are not extracted). Copies in the assembly, a lower bound
    (arrays collapse)."""
    units = [r for r in rows if r["kind"] == "rDNA"]
    names = [collections.Counter(n for n, _a, _b in (r.get("regions") or ())) for r in units]
    spans = collections.defaultdict(list)
    for r in rows:
        if r["kind"] == "18S":
            spans[(r.get("contig"), r.get("strand"))].append((r.get("start", 0), r.get("end", 0)))
    short = sum(1 for u, n in zip(units, names) if n["18S"] and not any(
        a <= u.get("end", 0) and u.get("start", 0) <= b for a, b in spans[(u.get("contig"), u.get("strand"))]))
    return {"rDNA_unit": (len(units), sum(1 for r in units if not r["lsu_complete"])),
            "18S": (sum(r["kind"] == "18S" for r in rows) + short, short),
            "5.8S": (sum(n["5.8S"] for n in names), 0),
            "28S": (sum(n["28S"] for n in names),
                    sum(n["28S"] for r, n in zip(units, names) if not r["lsu_complete"])),
            "5S": (sum(n["5S"] for n in names), 0)}


def load(folder):
    """A genome's saved extraction (folder: its done folder, or the packed genome's, store.py) -> extract_genome's
    rows ([] where there is none; contig_length and introns None in one saved before format 4)."""
    from . import store  # noqa: PLC0415
    rows = []
    for f in store.table(folder, "rrna.tsv"):
        length = f.get("contig_length", "-")
        rows.append(dict(kind=f["kind"], contig=f["contig"], start=int(f["start"]), end=int(f["end"]),
                         strand=f["strand"], regions=parse_regions(f["regions"]),
                         lsu_complete=f["lsu_complete"] == "t", sequence=f["sequence"],
                         contig_length=int(length) if length.isdigit() else None,
                         introns=parse_regions(f.get("introns", "-")) if "introns" in f else None))
    return rows


def gaps(seq, min_run=GAP_MIN):
    """The runs of min_run Ns or more in a sequence (assembly gaps)."""
    return len(re.findall(f"N{{{min_run},}}", seq.upper()))


def contig_edge(row):
    """A genome's record's bases to its contig's nearer end (None where not known)."""
    length = row.get("contig_length")
    if not length or row.get("start") in (None, "-"):
        return None
    return min(int(row["start"]) - 1, int(length) - int(row["end"]))


def confidence(row):
    """A record's tier and why (rrna/records.tsv's confidence and flags; package.tier's four since 2026-11-09): low
    where a check calls it another organism's or of bad sequence (own_species f, the part screen, ambiguous bases over
    MAX_AMBIGUOUS, a homopolymer of HOMOPOLYMER; since 2026-11-09 an 18S a bacterial or archaeal SSU model scores best:
    `ssu_prokaryote`, a bacterium's or an organelle's); medium where something may be wrong (an assembly gap in it, at
    a contig's end, an operon without a whole 28S, labelled by its genus alone, an 18S shared with another clade's
    genome, its owner undecided); unverified where nothing is wrong but its own-species check could not decide or was
    not made (a published read's: its label is its source's; until 2026-11-09 medium); else high."""
    from .package import tier  # noqa: PLC0415
    low, medium, unverified = [], [], []
    if row.get("own_species") == "f":
        low.append("own_species")
    if row.get("screen", "-") not in ("-", ""):
        low.append("screen")
    if row.get("ssu_model") in ("bacteria", "archaea"):
        low.append("ssu_prokaryote")
    if float(row.get("ambiguous") or 0) > MAX_AMBIGUOUS:
        low.append("ambiguous")
    if int(row.get("homopolymer") or 0) >= HOMOPOLYMER:
        low.append("homopolymer")
    if int(row.get("gaps") or 0):
        medium.append("gap")
    edge = row.get("contig_edge")
    if edge not in (None, "-") and int(edge) < EDGE:
        medium.append("contig_edge")
    if row.get("kind") == "rDNA" and row.get("lsu_complete") != "t":
        medium.append("lsu_partial")
    if row.get("label_from") == "eukaryome_genus":
        medium.append("genus_label")
    if row.get("own_species") in ("-", None):
        if row.get("cross") not in (None, "-", ""):
            medium.append("cross_shared")
        else:
            unverified.append("own_species_undecided")
    return tier(low, medium, unverified)
