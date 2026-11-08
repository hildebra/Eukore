"""A BUSCO gene's copy as Eukore keeps it (loci.call): the genome from the gene's first codon to its last on its own
strand, introns and all, uncut, unmasked and in upper case, with its exons and introns on it.

Until 2026-11-06 this was protal's euk30_copies.py, which fitted a copy into protal's 16 kb index by cutting its
introns' middles and masking soft-masked runs; Eukore fits nothing (its callers turned both off), and protal's own
database build does that fitting.
"""

# IUPAC's codes too, either case (since 2026-11-08 one table for the copies and rrna.py: the copies' left an R, Y,
# K, M, B, V, D or H of a minus-strand copy uncomplemented)
COMPLEMENT = str.maketrans("ACGTNRYKMBVDHSWacgtnrykmbvdhsw", "TGCANYRMKVBHDSWtgcanyrmkvbhdsw")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def whole(contig_seq, start, end, strand, exons):
    """A placed gene's copy. start, end: the gene's span on the contig, 1-based inclusive; exons [(start, end)] in
    genome order, 1-based inclusive, within the span; strand + or -. -> dict(copy (its strand, upper case), spliced
    (its exons joined: the coding sequence), features [(kind exon or intron, start, end)] on the copy (0-based, end
    excluded, the copy's own strand))."""
    exons = sorted((max(a, start), min(b, end)) for a, b in exons if b >= start and a <= end)
    span = contig_seq[start - 1:end].upper()
    segments, at = [], 0  # the span's exons and the introns between them, in genome order (0-based, end excluded)
    for a, b in exons:
        a0, b0 = a - start, b - start + 1
        if a0 > at:
            segments.append(("intron", at, a0))
        segments.append(("exon", a0, b0))
        at = b0
    if at < len(span):
        segments.append(("intron", at, len(span)))
    spliced = "".join(span[a0:b0] for kind, a0, b0 in segments if kind == "exon")
    if strand == "-":
        n = len(span)
        return dict(copy=revcomp(span), spliced=revcomp(spliced),
                    features=sorted(((k, n - e, n - s) for k, s, e in segments), key=lambda f: f[1]))
    return dict(copy=span, spliced=spliced, features=segments)
