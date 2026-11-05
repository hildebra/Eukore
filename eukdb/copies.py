#!/usr/bin/env python3
"""euk30_copies.py - a euk30 gene copy fitted into the index's MAX_LOCUS bases, every exon kept (since 2026-11-05;
docs/databases.md).

A euk30 copy is the genome from a gene's first codon to its last, introns and all (eukdetect1_genomic.py, EukDetect2's
copies): a read of the genome holds the introns, and protal aligns a read end to end within one gene. The index holds
a gene of MAX_LOCUS bases at most (eukaryotes.MAX_LOCUS, 16,000: its 14-bit position field). Many an animal's or a
plant's genes are longer with their introns (10 of the 30 loci in human and zebrafish, FNTB 74 kb): until 2026-11-05 such
a copy was left out, or cut to the run of its exons with the most coding bases (model species). Now it is fitted:

- every exon is kept whole;
- the bases over MAX_LOCUS come out of the introns, in proportion to each intron's length, from its middle (the
  inside): each intron keeps FLANK bases at either end at least, so that a read across a splice site still aligns (half
  of FLANK, then none, where the introns cannot give enough otherwise); the rounding remainder comes off the longest
  intron that has room;
- a gene whose exons alone are longer than MAX_LOCUS (none of euk30's loci in the genomes checked) keeps the run of its
  exons with the most coding bases (model_species.window), marked "exon window".

What was kept and cut is a copy's annotation (fit's features): its exons and introns on the copy as the reference
holds it (0-based, end excluded, the copy's own strand), each with where it lies in the genome, and each cut as a
zero-width row with the bases it removed. EUKDB fits nothing (no cap: loci.call); protal's converter does.
"""

import os
import re
import sys


FLANK = 50          # an intron's bases kept at either end at least
MAX_LOCUS = 16000   # eukaryotes.MAX_LOCUS (the caller passes it: a test may lower it)
MASK_MIN = 100      # eukdetect1_genomic.MASK_MIN: soft-masked runs this long outside the exons become N
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def masked_span(contig_seq, start, end, exons, mask_min=MASK_MIN):
    """contig_seq[start-1:end] (1-based inclusive) upper-case, its soft-masked runs of mask_min bases or more outside
    the exons N (eukdetect1_genomic.cut's rule) -> (span, masked bases)."""
    span = list(contig_seq[start - 1:end])
    inside = [False] * len(span)
    for a, b in exons:
        for i in range(max(a, start) - start, min(b, end) - start + 1):
            inside[i] = True
    masked, run = 0, []
    for i, c in enumerate(span + ["."]):
        if i < len(span) and c.islower() and not inside[i]:
            run.append(i)
            continue
        if len(run) >= mask_min:
            for j in run:
                span[j] = "N"
            masked += len(run)
        run = []
    return "".join(span).upper(), masked


def allocate(lengths, excess, flank):
    """Bases to cut from each intron (lengths) to remove `excess`: in proportion to the lengths, each keeping 2*flank
    bases at least, the rounding remainder off the longest that has room -> [cut per intron], or None where they cannot
    give enough."""
    room = [max(0, n - 2 * flank) for n in lengths]
    if sum(room) < excess:
        return None
    cuts = [0] * len(lengths)
    remaining = excess
    active = [i for i, r in enumerate(room) if r > 0]
    while remaining > 0 and active:
        total = sum(lengths[i] for i in active)
        give = {i: min(room[i] - cuts[i], remaining * lengths[i] // total) for i in active}
        if not any(give.values()):  # the remainder: off the longest with room
            i = max(active, key=lambda i: (lengths[i], -i))
            give = {i: min(room[i] - cuts[i], remaining)}
        for i, g in give.items():
            cuts[i] += g
            remaining -= g
        active = [i for i in active if cuts[i] < room[i]]
    return cuts if remaining == 0 else None


def fit(contig_seq, start, end, strand, exons, max_len=MAX_LOCUS, flank=FLANK, mask_min=MASK_MIN):
    """A placed gene's copy fitted into max_len bases (module docstring). start, end: the gene's span on the contig,
    1-based inclusive; exons [(start, end)] in genome order, 1-based inclusive, within the span; strand + or -.
    -> dict(copy (its strand, upper case, soft-masked runs N), spliced (its exons joined: the coding sequence), masked
    (bases masked), features [(kind, start, end, removed, genome_start, genome_end)] on the copy, how ("whole", "introns
    cut" or "exon window"))"""
    exons = sorted((max(a, start), min(b, end)) for a, b in exons if b >= start and a <= end)
    how = "whole"
    if sum(b - a + 1 for a, b in exons) > max_len:  # the exons alone over it: the run of them with the most coding
        from . import species as model_species  # noqa: PLC0415
        exons = model_species.window(exons, max_len)
        start, end, how = exons[0][0], exons[-1][1], "exon window"
    span, masked = masked_span(contig_seq, start, end, exons, mask_min)
    # the span's segments in genome order: (kind, 0-based start, end) of exons and the introns between them
    segments, at = [], 0
    for a, b in exons:
        a0, b0 = a - start, b - start + 1
        if a0 > at:
            segments.append(("intron", at, a0))
        segments.append(("exon", a0, b0))
        at = b0
    if at < len(span):
        segments.append(("intron", at, len(span)))
    introns = [i for i, s in enumerate(segments) if s[0] == "intron"]
    cuts = {}
    excess = len(span) - max_len
    if excess > 0:
        lengths = [segments[i][2] - segments[i][1] for i in introns]
        for f in (flank, flank // 2, 0):
            got = allocate(lengths, excess, f)
            if got is not None:
                cuts = {i: c for i, c in zip(introns, got) if c}
                break
        how = "introns cut" if how == "whole" else how
    pieces, features, pos = [], [], 0
    for i, (kind, a0, b0) in enumerate(segments):
        n = b0 - a0
        c = cuts.get(i, 0)
        g1, g2 = start + a0, start + b0 - 1  # genome, 1-based inclusive
        if not c:
            pieces.append(span[a0:b0])
            features.append((kind, pos, pos + n, 0, g1, g2))
            pos += n
            continue
        left = (n - c) // 2
        pieces += [span[a0:a0 + left], span[b0 - (n - c - left):b0]]
        features.append((kind, pos, pos + n - c, c, g1, g2))
        features.append(("cut", pos + left, pos + left, c, g1 + left, g1 + left + c - 1))
        pos += n - c
    copy = "".join(pieces)
    spliced = "".join(span[a0:b0] for kind, a0, b0 in segments if kind == "exon")
    if strand == "-":
        length = len(copy)
        copy, spliced = revcomp(copy), revcomp(spliced)
        features = sorted(((k, length - e, length - s, r, g1, g2) for k, s, e, r, g1, g2 in features),
                          key=lambda f: (f[1], f[0] == "cut"))
    return dict(copy=copy, spliced=spliced, masked=masked, features=features, how=how)
