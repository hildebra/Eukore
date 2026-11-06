#!/usr/bin/env python3
"""Tests of eukore/copies.py (a gene's copy whole, introns and all): the copy is the genome's span, upper case; its
features say where every exon and intron is, on either strand, and give the genome back; the exons joined are the
coding sequence.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import copies as ec  # noqa: E402


def gene(rng, exon_lengths, intron_lengths, before=100, after=100):
    """A contig with a gene of those exons and introns -> (contig, start, end, exons (1-based inclusive))."""
    seq = [rng.choice("ACGT") for _ in range(before)]
    exons = []
    for i, n in enumerate(exon_lengths):
        a = len(seq) + 1
        seq += [rng.choice("ACGT") for _ in range(n)]
        exons.append((a, len(seq)))
        if i < len(intron_lengths):
            seq += [rng.choice("acgt") for _ in range(intron_lengths[i])]  # soft-masked: kept, upper-cased
    seq += [rng.choice("ACGT") for _ in range(after)]
    return "".join(seq), exons[0][0], exons[-1][1], exons


class WholeTest(unittest.TestCase):
    def test_the_forward_strand(self):
        contig, s, e, exons = gene(random.Random(1), [100, 200], [300])
        got = ec.whole(contig, s, e, "+", exons)
        self.assertEqual(got["copy"], contig[s - 1:e].upper())
        self.assertEqual(got["features"], [("exon", 0, 100), ("intron", 100, 400), ("exon", 400, 600)])
        self.assertEqual(got["spliced"], "".join(contig[a - 1:b] for a, b in exons))

    def test_the_reverse_strand(self):
        contig, s, e, exons = gene(random.Random(5), [100, 100, 100], [1000, 3000])
        fwd = ec.whole(contig, s, e, "+", exons)
        rev = ec.whole(contig, s, e, "-", exons)
        self.assertEqual((rev["copy"], rev["spliced"]), (ec.revcomp(fwd["copy"]), ec.revcomp(fwd["spliced"])))
        self.assertEqual([f[0] for f in rev["features"]], ["exon", "intron", "exon", "intron", "exon"])
        self.assertEqual(rev["features"][0], ("exon", 0, 100))
        self.assertEqual(rev["features"][1], ("intron", 100, 3100))  # the last intron first
        n = e - s + 1
        for (kind, a, b), (a1, b1) in zip(rev["features"][::2], exons[::-1]):
            self.assertEqual(rev["copy"][a:b], ec.revcomp(contig[a1 - 1:b1]))
            self.assertEqual((a, b), (n - (b1 - s + 1), n - (a1 - s)))


if __name__ == "__main__":
    unittest.main()
