#!/usr/bin/env python3
"""Tests of the clades' expected BUSCOs and the genomes' tiers (eukore/tiers.py; since 2026-11-14).

    python3 -m unittest discover -s tests -t tests
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import package, tiers  # noqa: E402

LOCI = [f"L{i}" for i in range(10)]


def genome(level="Chromosome", size=10_000_000, contigs=10, source="ncbi", flags="-"):
    return dict(level=level, genome_size=size, contigs=contigs, source=source, flags=flags, representative="t")


class TiersTest(unittest.TestCase):
    def test_a_clade_that_lost_loci(self):
        # a parasite clade: its three chromosome-level genomes hold L0-L3 alone; a scaffold genome of it holding three
        # of them is a genome, not poor; a fragment of an animal (whose expected set is every locus) is poor
        genomes = {f"P{i}": genome() for i in range(3)}
        genomes["P3"] = genome("Scaffold", contigs=100)
        genomes["A0"] = genome("Scaffold", contigs=1000)
        statuses = {f"P{i}": {l: "Single" for l in LOCI[:4]} for i in range(3)}
        statuses["P3"] = {l: "Single" for l in LOCI[:3]}
        statuses["A0"] = {l: "Single" for l in LOCI[:3]}
        group = {"P0": "Metamonada", "P1": "Metamonada", "P2": "Metamonada", "P3": "Metamonada", "A0": "Metazoa"}
        sets = tiers.assign(genomes, statuses, group, set(LOCI))
        self.assertEqual(sets["Metamonada"], frozenset(LOCI[:4]))
        self.assertNotIn("Metazoa", sets)  # no reference genome: every BUSCO expected
        self.assertEqual((genomes["P3"]["complete_expected"], genomes["P3"]["tier"]), ("0.750", "genome"))
        self.assertEqual((genomes["A0"]["expected_loci"], genomes["A0"]["tier"]), (10, "poor"))

    def test_tiers(self):
        whole = {l: "Single" for l in LOCI}
        self.assertEqual(tiers.tier_of(genome(flags="partial genome"), 0.2), "partial")
        self.assertEqual(tiers.tier_of(genome(flags="small for its clade,derived from metagenome"), 0.3), "partial")
        self.assertEqual(tiers.tier_of(genome(flags="small for its clade"), None), "partial")
        # a "small" genome holding its expected BUSCOs is a genome (Setaria italica, 2026-11: 138 of 138 loci)
        self.assertEqual(tiers.tier_of(genome(flags="small for its clade"), 1.0), "genome")
        # every assembly of its species excluded (genomes.ONLY_EXCLUDED), as genomes.tsv joins the flags
        from eukore import genomes as genomes_mod  # noqa: PLC0415
        excluded = ",".join(sorted({"contaminated", genomes_mod.ONLY_EXCLUDED}))
        self.assertEqual(tiers.tier_of(genome(flags=excluded), 1.0), "contaminated")
        self.assertFalse(tiers.is_reference(genome(flags=excluded)))
        self.assertEqual(tiers.tier_of(genome(flags="contaminated," + tiers.OLD_EXCLUDED_FLAG), 1.0), "contaminated")
        self.assertEqual(tiers.tier_of(genome("Contig", size=3_000_000, contigs=2000), 0.9), "draft")
        self.assertEqual(tiers.tier_of(genome(level="transcriptome", source="tsa"), 0.9), "transcriptome")
        self.assertEqual(tiers.tier_of(genome(level="transcriptome", source="tsa"), 0.1), "poor")
        self.assertEqual(tiers.complete_expected(whole, frozenset(LOCI[:5])), 1.0)
        self.assertEqual(tiers.complete_expected({"L0": "Fragmented"}, frozenset(["L0", "L1"])), 0.0)
        self.assertFalse(tiers.is_reference(genome(flags="partial genome")))
        self.assertTrue(tiers.is_reference(genome("Scaffold", size=50_000_000, contigs=100)))

    def test_isoforms(self):
        # a transcriptome's Duplicated locus of copies 97% alike or more: its transcripts (package.collect: Single)
        self.assertTrue(package.isoforms(dict(dup_identity="99.2")))
        self.assertFalse(package.isoforms(dict(dup_identity="91.0")))  # a paralog
        self.assertFalse(package.isoforms(dict(dup_identity="<50")))
        self.assertFalse(package.isoforms(dict(dup_identity="-")))

    def test_eukore_queries(self):
        # an Odb with one dataset protein of a far species (taxid 1) and two Eukore proteins of the genome's genus (10:
        # another species; 11: the genome's own species)
        from eukore import odb as odb_mod  # noqa: PLC0415
        o = odb_mod.Odb.__new__(odb_mod.Odb)
        o.loci, o.markers = ["L0"], set()
        o.proteins = {"L0": [(1, "L0_1_x", "MAAA")]}
        o.ancestral = {}
        o.extra = {"L0": [(10, "L0_10_eukore_G1", "MCCC"), (11, "L0_11_eukore_G2", "MDDD")]}
        lineages = {1: [1, 100, 1000], 10: [10, 200, 300, 1000], 11: [11, 200, 300, 1000]}
        got = o.queries([11, 200, 300, 1000], lineages)
        self.assertEqual(sorted(got), ["L0~L0_10_eukore_G1", "L0~L0_1_x"])  # its own species' never
        self.assertEqual(sorted(o.queries([11, 200, 300, 1000], lineages, leave_out=[1])), [])  # a hold-out: none
        far = o.queries([1, 100, 1000], lineages)  # the dataset's own species is nearer than any Eukore protein
        self.assertEqual(sorted(far), ["L0~L0_1_x"])
        self.assertEqual(o.all_species(), [1, 10, 11])
        # the proteins taken: a clean Single copy of a genome of a served group, one per locus and genus
        rows = [dict(accession=a, locus="L0", status="Single", confidence="high", own_species="t", rescued="f",
                     internal_stops=0) for a in ("G1", "G2", "G3")]
        genomes = {"G1": dict(tier="genome", species_taxid=10), "G2": dict(tier="genome", species_taxid=11),
                   "G3": dict(tier="draft", species_taxid=12)}
        out = tiers.query_proteins(rows, {("G1", "L0"): "MCCC*", ("G2", "L0"): "MDDD", ("G3", "L0"): "MEEE"},
                                   genomes, {"Discoba"}, {"G1": "Discoba", "G2": "Discoba", "G3": "Discoba"},
                                   lambda sp: 200)
        self.assertEqual(out, [("L0_10_eukore_G1", "MCCC")])

    def test_resolution_regions(self):
        # since 2026-11-14 (benchmarks.region_span): the ITS from ITS1's start to ITS2's end; the 28S' first 900 bases
        from eukore import benchmarks  # noqa: PLC0415
        r = dict(regions="18S:0-1800,ITS1:1800-2000,5.8S:2000-2160,ITS2:2160-2400,28S:2400-5800")
        self.assertEqual(benchmarks.region_span(r, "ITS"), (1800, 2400))
        self.assertEqual(benchmarks.region_span(r, "28S D1-D2"), (2400, 3300))
        self.assertIsNone(benchmarks.region_span(dict(regions="18S:0-1800,28S:1900-2500"), "28S D1-D2"))
        self.assertIsNone(benchmarks.region_span(dict(regions="-"), "ITS"))

    def test_dropped_genome_searched_for_none(self):
        # package.collect gives a dropped genome's `_searched` () (since 2026-11-15): copy_numbers.tsv.gz holds no
        # zeros of it, busco/loci.tsv's call rates do not count it
        g = dict(representative="t", source="ncbi", species_taxid=4932, _searched=())
        kept = dict(representative="t", source="ncbi", species_taxid=4932, _searched=("L0",))
        loci = [dict(locus="L0")]
        self.assertEqual(package.copy_numbers({"D": g}, [], loci, {}), [])
        self.assertEqual(len(package.copy_numbers({"K": kept}, [], loci, {})), 1)

    def test_copy_flags(self):
        self.assertEqual(package.copy_confidence(dict(genome_tier="poor", own_species="t"))[0], "low")
        self.assertIn("poor_genome", package.copy_confidence(dict(genome_tier="partial", own_species="t"))[1])
        self.assertEqual(package.copy_confidence(dict(genome_tier="draft", own_species="t")),
                         ("medium", "draft_genome"))
        self.assertEqual(package.copy_confidence(dict(genome_tier="genome", own_species="t")), ("high", "-"))


if __name__ == "__main__":
    unittest.main()
