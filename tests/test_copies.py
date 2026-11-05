#!/usr/bin/env python3
"""Tests of eukdb/copies.py (copied from protal's test_euk30_copies.py) (a gene copy fitted into MAX_LOCUS bases, every exon kept): a gene that fits stays whole;
one over it loses bases from its introns' middles in proportion to their lengths, each intron keeping FLANK bases at
either end; the features say where every exon, intron and cut is, on either strand, and give the genome back; the
exons joined are the coding sequence; soft-masked intron runs are N; exons alone over the cap keep their best window.

    python3 -m unittest discover EUKDBbuild/tests
"""

import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukdb import copies as ec  # noqa: E402


def gene(rng, exon_lengths, intron_lengths, before=100, after=100):
    """A contig with a gene of those exons and introns -> (contig, start, end, exons (1-based inclusive))."""
    seq = [rng.choice("ACGT") for _ in range(before)]
    exons = []
    for i, n in enumerate(exon_lengths):
        a = len(seq) + 1
        seq += [rng.choice("ACGT") for _ in range(n)]
        exons.append((a, len(seq)))
        if i < len(intron_lengths):
            seq += [rng.choice("ACGT") for _ in range(intron_lengths[i])]
    seq += [rng.choice("ACGT") for _ in range(after)]
    return "".join(seq), exons[0][0], exons[-1][1], exons


class FitTest(unittest.TestCase):
    def test_a_gene_that_fits_stays_whole(self):
        contig, s, e, exons = gene(random.Random(1), [100, 200], [300])
        got = ec.fit(contig, s, e, "+", exons, max_len=1000)
        self.assertEqual((got["copy"], got["how"]), (contig[s - 1:e], "whole"))
        self.assertEqual([f[:4] for f in got["features"]], [("exon", 0, 100, 0), ("intron", 100, 400, 0),
                                                            ("exon", 400, 600, 0)])

    def test_introns_cut_in_proportion_from_their_middle(self):
        contig, s, e, exons = gene(random.Random(2), [100, 100, 100], [1000, 3000])
        got = ec.fit(contig, s, e, "+", exons, max_len=1300, flank=50)
        self.assertEqual((len(got["copy"]), got["how"]), (1300, "introns cut"))
        f = got["features"]
        introns = [x for x in f if x[0] == "intron"]
        self.assertEqual([x[3] for x in introns], [750, 2250])  # 3,000 over: 1/4 and 3/4
        self.assertEqual([x[2] - x[1] for x in introns], [250, 750])
        # each exon whole, the coding sequence the genome's exons joined
        self.assertEqual(got["spliced"], "".join(contig[a - 1:b] for a, b in exons))
        for kind, a, b, _r, g1, g2 in f:
            if kind == "exon":
                self.assertEqual(got["copy"][a:b], contig[g1 - 1:g2])
        # the cut in the first intron's middle: its 125 first and 125 last bases kept
        i1 = exons[0][1] + 1  # the first intron's first base (1-based)
        self.assertEqual(got["copy"][100:225], contig[i1 - 1:i1 - 1 + 125])
        self.assertEqual(got["copy"][225:350], contig[i1 - 1 + 875:i1 - 1 + 1000])
        cuts = [x for x in f if x[0] == "cut"]
        self.assertEqual([(x[1], x[2], x[3], x[4], x[5]) for x in cuts],
                         [(225, 225, 750, i1 + 125, i1 + 874), (350 + 100 + 375, 350 + 100 + 375, 2250,
                                                                 exons[1][1] + 1 + 375, exons[1][1] + 375 + 2250)])

    def test_flanks_kept_and_the_remainder_off_the_longest(self):
        contig, s, e, exons = gene(random.Random(3), [10, 10, 10, 10], [101, 300, 1000])
        got = ec.fit(contig, s, e, "+", exons, max_len=40 + 101 + 300 + 1000 - 7)  # 7 over: rounding
        cut = [x[3] for x in got["features"] if x[0] == "intron"]
        self.assertEqual(sum(cut), 7)
        self.assertEqual(cut[0], 0)  # 101 bases: 1 over its 2 x 50 flanks, its share (0.5) rounds to none
        self.assertGreaterEqual(cut[2], cut[1])
        # introns too short for full flanks: half of them, then none
        contig, s, e, exons = gene(random.Random(4), [100, 100], [120])
        got = ec.fit(contig, s, e, "+", exons, max_len=240)
        self.assertEqual(len(got["copy"]), 240)
        self.assertEqual([x[2] - x[1] for x in got["features"] if x[0] == "intron"], [40])  # 2 x 25 would be 50

    def test_the_reverse_strand(self):
        contig, s, e, exons = gene(random.Random(5), [100, 100, 100], [1000, 3000])
        fwd = ec.fit(contig, s, e, "+", exons, max_len=1300)
        rev = ec.fit(contig, s, e, "-", exons, max_len=1300)
        self.assertEqual(rev["copy"], ec.revcomp(fwd["copy"]))
        self.assertEqual(rev["spliced"], ec.revcomp(fwd["spliced"]))
        kinds = [x[0] for x in rev["features"]]
        self.assertEqual(kinds, ["exon", "intron", "cut", "exon", "intron", "cut", "exon"])
        for kind, a, b, _r, g1, g2 in rev["features"]:
            if kind == "exon":
                self.assertEqual(rev["copy"][a:b], ec.revcomp(contig[g1 - 1:g2]))
        self.assertEqual(rev["features"][0][1:3], (0, 100))
        self.assertEqual(rev["features"][-1][1:3], (1200, 1300))

    def test_soft_masked_intron_runs_are_n(self):
        contig, s, e, exons = gene(random.Random(6), [100, 100], [500])
        a = exons[0][1] + 100
        contig = contig[:a] + contig[a:a + 200].lower() + contig[a + 200:]
        got = ec.fit(contig, s, e, "+", exons, max_len=10000, mask_min=100)
        self.assertEqual((got["masked"], got["copy"].count("N")), (200, 200))

    def test_exons_alone_over_the_cap(self):
        contig, s, e, exons = gene(random.Random(7), [600, 300, 500], [100, 100])
        got = ec.fit(contig, s, e, "+", exons, max_len=1000)
        self.assertEqual(got["how"], "exon window")
        self.assertLessEqual(len(got["copy"]), 1000)
        self.assertEqual([x[2] - x[1] for x in got["features"] if x[0] == "exon"], [600, 300])  # 900 coding bases in 1,000


if __name__ == "__main__":
    unittest.main()
