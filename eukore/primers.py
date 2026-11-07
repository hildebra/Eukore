"""In-silico PCR of the barcode and amplicon-survey primers (since 2026-11-07; data/primers.tsv, data/primer_pairs.tsv)
on the package's rRNA records and organelle genes: rrna/amplicons.tsv, organelle/amplicons.tsv.

A primer's site on a record is where it matches with the fewest mismatches (IUPAC: a base the primer's code allows
matches; an N in the record matches nothing), the forward primer on the record's strand, the reverse primer's reverse
complement on it after the forward's site. A pair amplifies a record where each primer has MAX_MISMATCHES or fewer and
none in its 3'-most THREE_PRIME bases, and the product (primers included) is within the pair's length range. The
mismatch counts are given so that a user applies a rule of its own.
"""

import os

import numpy as np

from . import util

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIMERS = os.path.join(HERE, "data", "primers.tsv")
PAIRS = os.path.join(HERE, "data", "primer_pairs.tsv")
MAX_MISMATCHES = 3
THREE_PRIME = 2  # no mismatch in a primer's 3'-most two bases (ecoPCR's rule)
BITS = {"A": 1, "C": 2, "G": 4, "T": 8, "U": 8, "R": 5, "Y": 10, "S": 6, "W": 9, "K": 12, "M": 3, "B": 14, "D": 13,
        "H": 11, "V": 7, "N": 15, "I": 15}
COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNI", "TGCAAYRSWMKVHDBNN")
AMPLICON_COLUMNS = ("record", "pair", "forward_start", "forward_mismatches", "forward_3prime", "reverse_start",
                    "reverse_mismatches", "reverse_3prime", "start", "end", "length", "amplifies")


def read_tsv(path):
    with util.open_text(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    return [dict(zip(lines[0], f)) for f in lines[1:]]


def load(primers=PRIMERS, pairs=PAIRS):
    """The primers {name: row} and the pairs [row]."""
    return {r["name"]: r for r in read_tsv(primers)}, read_tsv(pairs)


def encode(seq):
    """A sequence as IUPAC bit masks (numpy uint8; an unknown character 0: matches nothing)."""
    table = np.zeros(256, dtype=np.uint8)
    for c, b in BITS.items():
        table[ord(c)] = table[ord(c.lower())] = b
    table[ord("N")] = table[ord("n")] = 0  # a record's N matches nothing (a primer's N is its own mask)
    return table[np.frombuffer(seq.encode(), dtype=np.uint8)]


def revcomp(seq):
    return seq.upper().translate(COMPLEMENT)[::-1]


def best_site(record_bits, primer, after=0):
    """The place (0-based, on the record) a primer's sequence (as written on the record's strand) matches with the
    fewest mismatches at or after `after` -> (start, mismatches, mismatches in its 3'-most THREE_PRIME bases, as the
    primer reads), or None where the record is shorter than the primer."""
    p = np.array([BITS.get(c, 15) for c in primer.upper()], dtype=np.uint8)
    if len(record_bits) - after < len(p):
        return None
    windows = np.lib.stride_tricks.sliding_window_view(record_bits[after:], len(p))
    miss = (windows & p) == 0
    counts = miss.sum(axis=1)
    i = int(counts.argmin())
    return after + i, int(counts[i]), miss[i]


def amplify(seq, pair, primers):
    """One pair on one record (its sequence, as the package holds it) -> an AMPLICON_COLUMNS row without `record`, or
    None where a primer is longer than what is left of the record."""
    bits = encode(seq)
    fwd, rev = primers[pair["forward"]]["sequence"], primers[pair["reverse"]]["sequence"]
    f = best_site(bits, fwd)
    if f is None:
        return None
    f_start, f_mm, f_miss = f
    rev_on_record = revcomp(rev)  # the reverse primer reads the other strand: its 3' end is the site's first base
    r = best_site(bits, rev_on_record, f_start + len(fwd))
    if r is None:
        return dict(pair=pair["pair"], forward_start=f_start, forward_mismatches=f_mm,
                    forward_3prime=int(f_miss[-THREE_PRIME:].sum()), amplifies="f")
    r_start, r_mm, r_miss = r
    f3, r3 = int(f_miss[-THREE_PRIME:].sum()), int(r_miss[:THREE_PRIME].sum())
    end = r_start + len(rev)
    length = end - f_start
    ok = (f_mm <= MAX_MISMATCHES and r_mm <= MAX_MISMATCHES and not f3 and not r3
          and int(pair["min_length"]) <= length <= int(pair["max_length"]))
    return dict(pair=pair["pair"], forward_start=f_start, forward_mismatches=f_mm, forward_3prime=f3,
                reverse_start=r_start, reverse_mismatches=r_mm, reverse_3prime=r3, start=f_start, end=end,
                length=length, amplifies="t" if ok else "f")


def runs_on(pair, row, kind):
    """Whether a pair runs on a record: an rRNA pair on rRNA records; an organelle pair on its gene (an animal
    mitochondrion's 12S and 16S: SSU and LSU of a mitochondrion)."""
    if kind == "rrna":
        return pair["records"] == "rrna"
    if pair["records"] != row.get("gene"):
        return False
    return row.get("organelle") == "mitochondrion" if row.get("gene") in ("SSU", "LSU") else True


def run(records, kind, primers=None, pairs=None):
    """Every pair on every record it runs on (records [(name, sequence, row)]; kind "rrna" or "organelle") -> rows
    (AMPLICON_COLUMNS) of the pairs whose forward primer fits the record (a pair's every record: none left out)."""
    if primers is None or pairs is None:
        primers, pairs = load()
    out = []
    for name, seq, row in records:
        for pair in pairs:
            if runs_on(pair, row, kind):
                got = amplify(seq, pair, primers)
                if got:
                    out.append(dict(got, record=name))
    return out


def write(path, rows):
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(AMPLICON_COLUMNS) + "\n")
        fh.writelines("\t".join(str(r.get(c, "-")) for c in AMPLICON_COLUMNS) + "\n" for r in rows)
    os.replace(util.part_of(path), path)
