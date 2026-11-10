#!/usr/bin/env python3
"""Tests of the audit's rRNA findings (2026-11-10: B1, L6, B7). The 18S' own-species rule (rrna.verdict): a lone
lineage's 18S is undecided, not another clade's; another clade's only where the copy's own clade holds the 18S a
margin farther and the other clade's other species hold it nearer; a hit a species and the clades searched again, so
that a species of many reads does not hide the others; the coding sequences' rule (qc.verdict) as it was, but for
its tie at the top identity. And the units: the check on an 18S' exons, an ITS1 over 3 kb, a unit's 5S.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import unittest

import tempfile
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi as ncbi_mod, qc, rrna  # noqa: E402
import fixtures as fx  # noqa: E402

# the fixtures' taxonomy and: a second class of Basidiomycota with a family of two species, and a second family of
# Cryptococcus' order
MORE = [(9910001, 5204, "class", "Agaricomycetes"), (9910002, 9910001, "order", "Agaricales"),
        (9910003, 9910002, "family", "Agaricaceae"), (9910004, 9910003, "genus", "Agaricus"),
        (9910005, 9910004, "species", "Agaricus bisporus"), (9910006, 9910004, "species", "Agaricus campestris"),
        (9910010, 5234, "family", "Tremellaceae"), (9910011, 9910010, "genus", "Tremella"),
        (9910012, 9910011, "species", "Tremella mesenterica")]
YEAST, PARADOXUS, LACTIS, CRYPTO, ASCARIS = 4932, 27291, 28985, 40410, 6253
AGARIC, AGARIC2, TREMELLA = 9910005, 9910006, 9910012


def taxonomy(tmp):
    with mock.patch.object(fx, "NODES", fx.NODES + MORE):
        return ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))


def mutate(seq, rate, rng):
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in seq)


class VerdictTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory() as tmp:
            cls.ncbi = taxonomy(tmp)

    def verdict(self, sp, found, holders, own_too=False):
        """rrna.verdict of a record of sp among the holders' records (own_too: another genome holds sp's 18S)."""
        cache = {}
        clades = rrna.clades_of(self.ncbi, holders, cache)
        own = rrna.clade_of(sp, clades, self.ncbi, cache, () if own_too else (sp,))
        return rrna.verdict(sp, found, own, clades, self.ncbi, cache)

    def coding(self, sp, found, holders):
        """The rule until 2026-11-10 (qc.verdict by the 18S' thresholds), found as its list."""
        cache = {}
        return qc.verdict(sp, [(i, s) for s, i in found.items()], self.ncbi, cache,
                          qc.families(self.ncbi, holders, cache), rrna.NEAR, rrna.FOREIGN_ID, rrna.RELATED)[:2]

    def test_a_lone_lineages_18s_is_undecided(self):
        # B1: Cryptococcus alone in its family, order and class among the records, its 18S 98.6% like another class'
        # (18S is that alike across classes): until 2026-11-10 another clade's (f, low), that family's second
        # species "holding it nearer" than the none of its own
        holders = {CRYPTO, AGARIC, AGARIC2, YEAST, PARADOXUS}
        found = {AGARIC: 98.6, AGARIC2: 98.2}
        self.assertEqual(self.coding(CRYPTO, found, holders), ("f", AGARIC))
        self.assertEqual(self.verdict(CRYPTO, found, holders), ("-", None))
        # one sequence in two clades (SHARED_ID): of its own phylum undecided and shared; another phylum's, that
        # clade's other species holding it too: another clade's; no other species of theirs holding it: shared
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 99.8, AGARIC2: 99.6}, holders), ("-", AGARIC))
        self.assertEqual(self.verdict(CRYPTO, {YEAST: 99.8, PARADOXUS: 99.0}, holders), ("f", YEAST))
        self.assertEqual(self.verdict(CRYPTO, {YEAST: 99.8, PARADOXUS: 89.0}, holders), ("-", YEAST))
        self.assertEqual(self.verdict(CRYPTO, {YEAST: 99.4, PARADOXUS: 99.0}, holders), ("-", None))
        # no hit, a hit under FOREIGN_ID, a near hit of its own class: nothing against it
        self.assertEqual(self.verdict(CRYPTO, {}, holders), ("-", None))
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 96.9}, holders), ("-", None))

    def test_another_clades_by_its_own_clades_margin(self):
        # a species of its order holds the 18S (its own clade: the deepest of family, order, class another holder is
        # of): another clade's only where that is MARGIN under the top and the partner's clade holds it nearer
        holders = {CRYPTO, TREMELLA, AGARIC, AGARIC2}
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, AGARIC2: 98.2, TREMELLA: 96.7}, holders), ("-", None))
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, AGARIC2: 98.2, TREMELLA: 96.5}, holders), ("f", AGARIC))
        # the partner's clade's other species not holding it nearer than its own: shared, undecided
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, AGARIC2: 92.0, TREMELLA: 93.0}, holders), ("-", AGARIC))
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, AGARIC2: 89.9, TREMELLA: 85.0}, holders), ("-", AGARIC))
        # the partner alone in its clade (none to tell): another clade's where its own clade holds nothing like it
        # (under RELATED), else shared
        alone = holders - {AGARIC2}
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, TREMELLA: 85.0}, alone), ("f", AGARIC))
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, TREMELLA: 93.0}, alone), ("-", AGARIC))
        # its nearest of its family, or of its class within NEAR of the top: its own, or nothing against it
        self.assertEqual(self.verdict(YEAST, {PARADOXUS: 99.0, AGARIC: 99.4}, holders | {YEAST, PARADOXUS}),
                         ("t", None))
        self.assertEqual(self.verdict(CRYPTO, {AGARIC: 98.6, TREMELLA: 98.2}, holders), ("-", None))

    def test_the_owner_of_a_sequence_two_clades_hold(self):
        # Ascaris' assembly holds the yeast's 18S: the yeast's is its own (its family within MARGIN of the top and
        # nearer than Ascaris' clade), Ascaris' the yeast's (alone in its clade, 100% like another phylum's whose
        # family holds it too)
        holders = {YEAST, PARADOXUS, LACTIS, ASCARIS}
        self.assertEqual(self.verdict(YEAST, {ASCARIS: 100.0, PARADOXUS: 99.0, LACTIS: 97.0}, holders), ("t", None))
        self.assertEqual(self.verdict(ASCARIS, {YEAST: 100.0, PARADOXUS: 99.0, LACTIS: 97.0}, holders), ("f", YEAST))
        # its family far from the copy and the other clade's nearer: not the owner, whatever the family's 91%
        # (qc.verdict without `related`, until b84d54f, called such a copy t)
        found = {AGARIC: 99.0, AGARIC2: 98.0, PARADOXUS: 91.0}
        self.assertEqual(self.verdict(YEAST, found, holders | {AGARIC, AGARIC2}), ("f", AGARIC))
        # another genome of its own species is of its clade: far from this copy, the copy is not its species'
        self.assertEqual(self.verdict(YEAST, {AGARIC: 99.0, AGARIC2: 98.0, YEAST: 91.0}, {YEAST, AGARIC, AGARIC2},
                                      own_too=True), ("f", AGARIC))
        self.assertEqual(self.verdict(YEAST, {AGARIC: 99.0, AGARIC2: 98.0}, {YEAST, AGARIC, AGARIC2}), ("-", None))

    def test_a_tie_at_the_top_is_not_told_by_the_taxid(self):
        # L6: two species at the top identity, one whose family's second species holds the sequence, one alone in
        # its family: the partner is the first, whichever has the higher taxid (qc.verdict takes the higher taxid:
        # its verdict turns with NCBI's numbering)
        holders = {YEAST, PARADOXUS, AGARIC, AGARIC2, TREMELLA, CRYPTO}
        for lone in (TREMELLA, CRYPTO):  # a higher taxid than Agaricus', and a lower
            found = {AGARIC: 99.0, lone: 99.0, AGARIC2: 98.0, PARADOXUS: 91.0}
            self.assertEqual(self.verdict(YEAST, found, holders), ("f", AGARIC))
        self.assertTrue(CRYPTO < AGARIC < TREMELLA)
        found = {AGARIC: 99.0, AGARIC2: 98.0, PARADOXUS: 91.0}
        self.assertEqual({self.coding(YEAST, found | {lone: 99.0}, holders) for lone in (TREMELLA, CRYPTO)},
                         {("t", None), ("f", AGARIC)})

    def test_the_coding_sequences_rule_as_it_was(self):
        # qc.verdict, the BUSCO copies', markers' and organelle genes' rule, is not the 18S': a lone lineage's
        # coding sequence 90% or more like another class' whose family has a second species is another clade's
        cache = {}
        holders = {CRYPTO, AGARIC, AGARIC2, YEAST, PARADOXUS, LACTIS, ASCARIS}
        members = qc.families(self.ncbi, holders, cache)
        self.assertEqual(qc.verdict(CRYPTO, [(93.0, AGARIC), (92.5, AGARIC2)], self.ncbi, cache, members),
                         ("f", AGARIC, 93.0, 1))
        self.assertEqual(qc.verdict(CRYPTO, [(89.0, AGARIC), (88.0, AGARIC2)], self.ncbi, cache, members)[:2],
                         ("-", None))
        # test_genomes' cases of the owner told (2026-11-09), as they were
        relatives = [(93.0, PARADOXUS), (85.0, LACTIS)]
        self.assertEqual(qc.verdict(YEAST, [(100.0, ASCARIS)] + relatives, self.ncbi, cache, members)[:2],
                         ("t", None))
        self.assertEqual(qc.verdict(ASCARIS, [(100.0, YEAST)] + relatives, self.ncbi, cache, members)[:2],
                         ("f", YEAST))
        self.assertEqual(qc.verdict(YEAST, [(100.0, CRYPTO)], self.ncbi, cache, members)[:2], ("f", CRYPTO))
        self.assertEqual(qc.verdict(CRYPTO, [(100.0, YEAST)], self.ncbi, cache, members)[:2], ("-", None))
        self.assertEqual(qc.verdict(YEAST, [(99.0, PARADOXUS), (99.5, ASCARIS)], self.ncbi, cache, members)[0], "t")
        self.assertEqual((qc.NEAR, qc.FOREIGN_ID), (1.0, 90.0))


def genome(n, taxid, seq, **row):
    """A genome's 18S record (assembly GCA_10000000n.1)."""
    return (f"GCA_10000000{n}.1_18S_1", seq, dict(dict(kind="18S", source="genome", accession=f"GCA_10000000{n}.1",
                                                        taxid=taxid, species_taxid=taxid, regions_list=None), **row))


def read(name, taxid, seq):
    """A EUKARYOME read's record."""
    return (name, seq, dict(kind="18S", source="EUKARYOME", accession="-", taxid=taxid, species_taxid=taxid,
                            regions_list=None))


class OwnSpeciesTest(unittest.TestCase):
    """rrna.own_species with vsearch's stand-in (fixtures.VSEARCH: a query's first 32 hits by identity)."""

    def setUp(self):
        self.rng = random.Random(31)
        self.yeast = "".join(self.rng.choice("ACGT") for _ in range(1200))

    def test_a_species_of_many_reads_does_not_hide_the_clades(self):
        # L6: a yeast genome holding an agaric's 18S, the agaric a species of 40 reads: vsearch's 32 hits are all
        # its reads. Until 2026-11-10 neither its family's second species (98% like the copy) nor the yeast's own
        # family (85%) was seen, and the copy was undecided ("neither clade tells"); the two clades searched alone
        # tell: the agaric's family holds it, the yeast's nothing like it
        rng = self.rng
        agaric = mutate(self.yeast, 0.15, rng)
        records = [genome(1, YEAST, mutate(agaric, 0.002, rng)), genome(2, PARADOXUS, mutate(self.yeast, 0.01, rng)),
                   genome(3, LACTIS, mutate(self.yeast, 0.03, rng))]
        records += [read(f"EUK_a{i}", AGARIC, mutate(agaric, 0.002, rng)) for i in range(40)]
        records += [read("EUK_b", AGARIC2, mutate(agaric, 0.02, rng))]
        with tempfile.TemporaryDirectory() as tmp:
            own = rrna.own_species(records, taxonomy(tmp), fx.vsearch(tmp))
        self.assertEqual(own["GCA_100000001.1_18S_1"], ("f", AGARIC))
        self.assertEqual((own["GCA_100000002.1_18S_1"], own["GCA_100000003.1_18S_1"]), (("t", None), ("t", None)))
        # the agaric's reads, set beside the genomes alone: like a yeast genome's record, and no genome of their
        # class, order or family to tell whose it is: undecided, not the yeast's
        self.assertEqual({own[f"EUK_a{i}"][0] for i in range(40)} | {own["EUK_b"][0]}, {"-"})

    def test_an_18s_with_an_intron_is_searched_as_its_exons_and_its_span(self):
        # B7: a 300-base intron in a yeast's 18S (vsearch counts its bases as a gap's): its relative's 18S is 99%
        # like its exons; until 2026-11-10 the span was searched, no relative found
        rng = self.rng
        intron = "".join(rng.choice("ACGT") for _ in range(300))
        relative = genome(2, PARADOXUS, mutate(self.yeast, 0.01, rng))
        cut = genome(1, YEAST, self.yeast[:500] + intron + self.yeast[500:], introns_list=[("18S", 500, 800)])
        with tempfile.TemporaryDirectory() as tmp:
            ncbi, vsearch = taxonomy(tmp), fx.vsearch(tmp)
            own = rrna.own_species([cut, relative], ncbi, vsearch)
            self.assertEqual((own[cut[0]], own[relative[0]]), (("t", None), ("t", None)))
            # the record without its introns' list (one saved before format 4): its span alone, as it was
            bare = (cut[0], cut[1], dict(cut[2], introns_list=None))
            self.assertEqual(rrna.own_species([bare, relative], ncbi, vsearch)[cut[0]], ("-", None))
            # an "intron" that is 18S (barrnap's pieces parted by an assembly gap, the trial's GCA_000151665.1): the
            # span is the better of the two, the record its family's still
            gap = genome(1, YEAST, self.yeast, introns_list=[("18S", 500, 800)])
            self.assertEqual(rrna.own_species([gap, relative], ncbi, vsearch)[gap[0]], ("t", None))


class UnitTest(unittest.TestCase):
    def setUp(self):
        rng = random.Random(21)
        self.seq = "".join(rng.choice("ACGT") for _ in range(30000))

    def test_a_units_5s_after_its_28s_before_one_before_its_18s(self):
        # one unit, a 5S 580 bases before its 18S and one 2,000 after its 28S: the one after the 28S, as the units
        # further on in an array take theirs (until 2026-11-10 the nearer whatever its side: the one before)
        feats = {"18S_rRNA": [(5001, 6800, "+", False)], "5_8S_rRNA": [(7201, 7350, "+", False)],
                 "28S_rRNA": [(7601, 16100, "+", False)]}
        before, after = (4301, 4420, "+", False), (18101, 18220, "+", False)
        for five in ([before, after], [after, before]):  # in either order of barrnap's lines
            unit = rrna.units_of(self.seq, dict(feats, **{"5S_rRNA": five}))[0]
            self.assertEqual(unit["unit_span"], (5001, 18220))
            self.assertEqual([r[0] for r in unit["regions"]], ["18S", "ITS1", "5.8S", "ITS2", "28S", "5S"])
        # of two after its 28S the nearer; one before its 18S alone is still taken
        near = (16601, 16720, "-", False)
        unit = rrna.units_of(self.seq, dict(feats, **{"5S_rRNA": [after, near, before]}))[0]
        self.assertEqual(unit["unit_span"], (5001, 16720))
        unit = rrna.units_of(self.seq, dict(feats, **{"5S_rRNA": [before]}))[0]
        self.assertEqual((unit["unit_span"], unit["regions"][0][0]), ((4301, 16100), "5S"))

    def test_an_its1_over_3_kb_is_a_unit(self):
        # an ITS1 of 3,400 bases (Trinematidae's, an opossum's 3.1 kb; until 2026-11-10 MAX_ITS 3000 for both
        # spacers: no unit, the read an "18S" record of its whole length); one over MAX_ITS1 is none still, and
        # the ITS2's bound is as it was
        def feats(its1, its2=250):
            s58 = 2800 + its1 + 1
            s28 = s58 + 150 + its2
            return {"18S_rRNA": [(1001, 2800, "+", False)], "5_8S_rRNA": [(s58, s58 + 149, "+", False)],
                    "28S_rRNA": [(s28, s28 + 3399, "+", False)]}
        unit = rrna.units_of(self.seq, feats(3400))[0]
        self.assertEqual([r for r in unit["regions"] if r[0] == "ITS1"], [("ITS1", 1800, 5200)])
        self.assertEqual(unit["unit_span"], (1001, 2800 + 3400 + 150 + 250 + 3400))
        self.assertIsNone(rrna.units_of(self.seq, feats(rrna.MAX_ITS1))[0]["unit"])
        self.assertIsNotNone(rrna.units_of(self.seq, feats(rrna.MAX_ITS1 - 1))[0]["unit"])
        self.assertIsNotNone(rrna.units_of(self.seq, feats(300, rrna.MAX_ITS - 150))[0]["unit"])
        self.assertIsNone(rrna.units_of(self.seq, feats(300, rrna.MAX_ITS + 1))[0]["unit"])

    def test_an_18s_without_its_introns(self):
        seq = self.seq[:6000]
        # an 18S record: the span, its introns on it
        row = dict(kind="18S", introns_list=[("18S", 700, 900), ("18S", 1500, 1600)])
        self.assertEqual(rrna.ssu_exons(seq[:2000], row), seq[:700] + seq[900:1500] + seq[1600:2000])
        self.assertEqual(rrna.ssu_exons(seq[:2000], dict(kind="18S")), seq[:2000])  # none, or none saved
        # an operon: its 18S region less the introns on the record; an intron of another region is not the 18S'
        row = dict(kind="rDNA", regions_list=[("5S", 0, 120), ("18S", 500, 2500), ("28S", 3000, 6000)],
                   introns_list=[("18S", 1200, 1400), ("28S", 4000, 4100)])
        self.assertEqual(rrna.ssu_exons(seq, row), seq[500:1200] + seq[1400:2500])
        # a read without regions: its first bases, as ssu_part
        self.assertEqual(rrna.ssu_exons(seq, dict(kind="rDNA", regions_list=None)), seq[:rrna.PREFIX])


if __name__ == "__main__":
    unittest.main()
