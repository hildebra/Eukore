#!/usr/bin/env python3
"""Tests of the copies' genetics (since 2026-11-07): proteins placed on their HMMs (align.py), introns' phase, class,
U12 look and shared sites (introns.py), a genome's genetic code checked (codecheck.py).

    python3 -m unittest discover -s tests -t tests
"""

import hashlib
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import align, codecheck, gencode, introns, util  # noqa: E402
import fixtures as fx  # noqa: E402

AMINO = "ACDEFGHIKLMNPQRSTVWY"


def favouring(aa, p=0.9):
    return {a: p if a == aa else (1 - p) / 19 for a in AMINO}


class AlignTest(unittest.TestCase):
    def test_match_states(self):
        # a residue, a deletion, an insert (lower case), a residue: states 1, (2 skipped), None, 3
        self.assertEqual(align.states_of("A-gC", "xx.x"), [1, None, 3])
        self.assertEqual(align.states_of("..AC", "..xx"), [1, 2])
        # since 2026-11-10 a row kept on the profile alone (A2M): the same states, the other rows' inserts left out
        for row, rf in (("A-gC", "xx.x"), ("..AC", "..xx"), ("A-.C", "xx.x"), ("Mk-l", "x.xx")):
            self.assertEqual(align.states_of_a2m(align.a2m_of(row, rf)), align.states_of(row, rf))
        self.assertEqual(align.a2m_of("A-.gC", "xx..x"), "A-gC")

    def test_emissions_and_hmmalign(self):
        with tempfile.TemporaryDirectory() as tmp:
            hmm = os.path.join(tmp, "1at2759.hmm")
            with open(hmm, "w") as fh:
                fh.write(fx.hmm_text("1", [favouring("L"), favouring("S", 0.6)]))
            em = align.read_emissions(hmm)
            self.assertEqual(len(em), 3)
            self.assertAlmostEqual(em[1]["L"], 0.9, places=3)
            self.assertAlmostEqual(em[2]["S"], 0.6, places=3)
            got = align.align_locus(hmm, {"a": "MKL", "b": "MK*LV"}, fx.hmmalign(tmp), os.path.join(tmp, "cache"))
            self.assertEqual(got, {"a": [1, 2, 3], "b": [1, 2, 3, 4, 5]})
            self.assertEqual(len(os.listdir(os.path.join(tmp, "cache"))), 1)  # kept: a second call reads it
            self.assertEqual(align.align_locus(hmm, {"a": "MKL", "b": "MK*LV"}, "no-hmmalign",
                                               os.path.join(tmp, "cache")), got)
            # since 2026-11-10 kept by protein: a locus' protein more is aligned alone, a protein twice once; the
            # rest read whatever hmmalign is there
            cache = os.path.join(tmp, "cache")
            more = align.align_locus(hmm, {"a": "MKL", "b": "MK*LV", "c": "MKLVS", "d": "MKL"}, fx.hmmalign(tmp),
                                     cache)
            self.assertEqual(more, dict(got, c=[1, 2, 3, 4, 5], d=[1, 2, 3]))
            self.assertEqual(sorted(util.cache_read(os.path.join(cache, f"1at2759.{util.md5(hmm)[:16]}"))), sorted(
                hashlib.md5(p.encode()).hexdigest() for p in ("MKL", "MKXLV", "MKLVS")))
            self.assertEqual(align.align_locus(hmm, {"c": "MKLVS", "e": "MKL"}, "no-hmmalign", cache),
                             {"c": [1, 2, 3, 4, 5], "e": [1, 2, 3]})
            self.assertEqual(align.align_locus(hmm, {"f": "WWW"}, "no-hmmalign", cache), {})


class IntronTest(unittest.TestCase):
    U12 = "ATATCCTT" + "G" * 20 + "TCCTTAAC" + "C" * 10 + "AC"  # a U12 donor, its branch point 12 from the end
    U2 = "GTAAGT" + "T" * 40 + "CAG"

    def test_an_intron_described(self):
        copy = "ATGAAAC" + self.U12 + "AAGTAA"
        a, b = 7, 7 + len(self.U12)
        feats = [("exon", 0, 7), ("intron", a, b), ("exon", b, len(copy))]
        (row,) = introns.of_copy(copy, feats, [11, 12, 13, 14])
        self.assertEqual((row["start"], row["end"], row["phase"], row["donor"], row["acceptor"], row["class"],
                          row["u12"], row["aa_position"], row["hmm_position"]), (a, b, 1, "AT", "AC", "AT-AC", "t", 2, 13))
        copy = "ATGAAAC" + self.U2 + "AAGTAA"
        b = 7 + len(self.U2)
        (row,) = introns.of_copy(copy, [("exon", 0, 7), ("intron", 7, b), ("exon", b, len(copy))])
        self.assertEqual((row["class"], row["u12"], row["hmm_position"]), ("GT-AG", "f", "-"))
        # a frameshift's base in the first exon is not coding: the intron falls between codons
        (row,) = introns.of_copy(copy, [("exon", 0, 7), ("frameshift", 2, 3), ("intron", 7, b),
                                        ("exon", b, len(copy))])
        self.assertEqual((row["phase"], row["aa_position"]), (0, 2))

    def test_sites(self):
        rows = [dict(hmm_position=10, phase=0)] * 3 + [dict(hmm_position=11, phase=0), dict(hmm_position=10, phase=1),
                                                       dict(hmm_position=30, phase=2), dict(hmm_position="-", phase=0)]
        self.assertEqual(introns.sites(rows), {(10, 0): "10.0", (11, 0): "10.0", (10, 1): "10.1", (30, 2): "30.2"})

    def test_a_genome_without_coding_bases(self):
        # since 2026-11-10: a genome whose copies' coding sequences are all empty gets `-`, not a division by zero
        row = dict(accession="B", locus="L1", intron=1, start=7, end=57, length=50, phase=0, donor="GT", acceptor="AG",
                   u12="f", aa_position=2, **{"class": "GT-AG"})
        collected = dict(intron_rows=[row], cds={("A", "L1"): "", ("A", "L2"): "", ("B", "L1"): "ATGAAA" * 100})
        _rows, _sites, _shares, genome = introns.annotate(collected, {})
        self.assertEqual(genome, {"A": dict(introns=0, introns_per_kb_cds="-", median_intron_length="-"),
                                  "B": dict(introns=1, introns_per_kb_cds="1.67", median_intron_length=50)})


class CodeTest(unittest.TestCase):
    def genome(self, table, times=12):
        """A CTG-clade yeast's gene: each sense codon of `table` at a position favouring its amino acid there."""
        code = gencode.code(table)
        codons = [c for c in gencode.CODONS if code[c] != "*"]
        cds = "".join(c * 1 for c in codons) * times
        emissions = [None] + [favouring(code[c]) for c in codons] * times
        return cds, list(range(1, len(codons) * times + 1)), {"L1": emissions}

    def test_a_ctg_clade_genome(self):
        cds, states, emissions = self.genome(12)
        tallies = codecheck.tally([(cds, states, "L1")], emissions)
        rows, verdict, best = codecheck.check(tallies, 1)  # NCBI's table wrong: the standard code
        self.assertEqual((verdict, best), ("differs:CTG>S", 12))
        ctg = next(r for r in rows if r["codon"] == "CTG")
        self.assertEqual((ctg["table_aa"], ctg["inferred_aa"], ctg["n"]), ("L", "S", 12))
        self.assertEqual(codecheck.check(tallies, 12)[1:], ("agrees", 12))

    def test_a_ciliate_read_through_its_stops(self):
        cds, states, emissions = self.genome(6)  # TAA and TAG are Q
        tallies = codecheck.tally([(cds, states, "L1")], emissions)
        self.assertEqual(codecheck.check(tallies, 1)[1:], ("differs:TAA>Q,TAG>Q", 6))
        self.assertEqual(codecheck.check(tallies, 6)[1:], ("agrees", 6))

    def test_only_the_reassignable_codons_judged(self):
        # ATG at positions favouring L (Met rarely a conserved position's best): reported, not judged
        cds, states, emissions = self.genome(1)
        for i in range(0, len(cds), 3):
            if cds[i:i + 3] == "ATG":
                emissions["L1"][i // 3 + 1] = favouring("L")
        tallies = codecheck.tally([(cds, states, "L1")], emissions)
        rows, verdict, best = codecheck.check(tallies, 1)
        self.assertEqual(next(r["inferred_aa"] for r in rows if r["codon"] == "ATG"), "L")
        self.assertEqual((verdict, best, codecheck.judged(tallies)), ("agrees", 1, "CTG"))  # no stop seen

    def test_a_stop_the_table_reads_but_the_genome_does_not_use(self):
        # since 2026-11-11: Fabrea salina's case, NCBI's table 10 (TGA Cys), no TGA inside its copies: TGA a stop
        cds, states, emissions = self.genome(1, times=100)  # 200 TGT and TGC, no TGA
        tallies = codecheck.tally([(cds, states, "L1")], emissions)
        within = codecheck.inside([cds])
        self.assertEqual((within["TGA"], within["TGT"] + within["TGC"]), (0, 200))
        self.assertEqual(codecheck.check(tallies, 10, within)[1:], ("differs:TGA>*", 1))
        self.assertEqual(codecheck.check(tallies, 10)[1:], ("agrees", 10))  # not judged without the copies' codons
        self.assertEqual(codecheck.check(tallies, 1, within)[1:], ("agrees", 1))  # a stop read as one: not judged
        # 1% of the other codons is unused; more is a codon used rarely; fewer than 200 other codons say nothing
        self.assertEqual(codecheck.unused(dict(TGA=2, TGT=100, TGC=100), gencode.code(10)), ["TGA"])
        self.assertEqual(codecheck.unused(dict(TGA=3, TGT=100, TGC=100), gencode.code(10)), [])
        self.assertEqual(codecheck.unused(dict(TGT=100, TGC=99), gencode.code(10)), [])
        # Strombidium's TAA, read as Q at 7% of its Q codons, Oxytricha's TAG at 26%: used
        self.assertEqual(codecheck.unused(dict(TAA=98, TAG=351, CAA=700, CAG=677), gencode.code(6)), [])
        self.assertEqual(codecheck.unused(dict(TAA=98, CAA=700, CAG=677), gencode.code(6)), ["TAG"])

    def test_a_ciliate_using_its_stops_inside(self):
        cds, states, emissions = self.genome(6, times=100)  # TAA and TAG Q, as often as CAA and CAG
        tallies = codecheck.tally([(cds, states, "L1")], emissions)
        self.assertEqual(codecheck.check(tallies, 6, codecheck.inside([cds]))[1:], ("agrees", 6))
        self.assertEqual(codecheck.inside(["ATGTAATAGTGA"]), {"ATG": 1, "TAA": 1, "TAG": 1})  # the last left out

    def test_too_little_to_say(self):
        cds, states, emissions = self.genome(1, times=2)
        self.assertEqual(codecheck.check(codecheck.tally([(cds, states, "L1")], emissions), 1)[1:], ("undecided", "-"))


if __name__ == "__main__":
    unittest.main()
