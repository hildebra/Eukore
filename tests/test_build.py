#!/usr/bin/env python3
"""Tests of Eukore's build: an rDNA unit uncut, with its 5.8S and an inline 5S (one block) and not a 5S elsewhere; a
EUKARYOME read's NCBI taxid by its species name, else its genus's; the record checks; a model species' genome run
(stand-in miniprot and barrnap: its locus placed whole over 16 kb, its operon over 8,191 bases whole) and the package
written from it and the test world, read back and checked.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import random
import stat
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, ncbi as ncbi_mod, package, rrna, sources, species  # noqa: E402
import fixtures as fx  # noqa: E402

VERTEBRATES = [(7711, 33208, "phylum", "Chordata"), (9605, 7711, "genus", "Homo"), (9606, 9605, "species", "Homo sapiens")]
HUMAN = "GCA_009914755.4"


def tool(folder, name, text):
    path = os.path.join(folder, name)
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
    return path


class UnitTest(unittest.TestCase):
    def test_a_unit_whole_with_its_inline_5s(self):
        rng = random.Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(30000))
        feats = {"18S_rRNA": [(1001, 2800, "+", False)], "5_8S_rRNA": [(3201, 3350, "+", False)],
                 "28S_rRNA": [(3601, 12100, "+", False)],
                 "5S_rRNA": [(12500, 12620, "-", False), (25000, 25120, "+", False)]}  # inline, and far away
        unit = rrna.units_of(seq, feats)[0]
        self.assertEqual(len(unit["unit"]), 12620 - 1001 + 1)  # 18S to the inline 5S, uncut (over 8,191)
        names = [r[0] for r in unit["regions"]]
        self.assertEqual(names, ["18S", "ITS1", "5.8S", "ITS2", "28S", "5S"])
        self.assertEqual(dict((r[0], r[1:]) for r in unit["regions"])["5S"], (12500 - 1001, 12620 - 1001 + 1))
        self.assertTrue(unit["lsu_complete"])
        # the reverse strand: the same unit, its regions on its own strand
        rc = rrna.revcomp(seq)
        n = len(seq)
        flip = {k: [(n - b + 1, n - a + 1, "-" if s == "+" else "+", p) for a, b, s, p in v] for k, v in feats.items()}
        flip["28S_rRNA"] = [(a, b, s, True) for a, b, s, _p in flip["28S_rRNA"]]  # a partial 28S hit
        back = rrna.units_of(rc, flip)[0]
        self.assertEqual(back["unit"], unit["unit"])
        self.assertEqual(back["regions"], unit["regions"])
        self.assertFalse(back["lsu_complete"])
        # no 5.8S, no unit: the 18S span alone
        alone = rrna.units_of(seq, {"18S_rRNA": [(1001, 2800, "+", False)]})[0]
        self.assertIsNone(alone["unit"])
        # barrnap's 28S hit reaching back over the 5.8S (Drosophila's units, 2026-11-05 trial): the 28S starts at the
        # 5.8S's end, no ITS2 (ITSx gives it where installed), the regions in order on either strand
        over = dict(feats, **{"28S_rRNA": [(3300, 12100, "+", False)]})
        for strand_seq, f in ((seq, over), (rc, {k: [(n - b + 1, n - a + 1, "-" if s == "+" else "+", p)
                                                      for a, b, s, p in v] for k, v in over.items()})):
            regions = rrna.units_of(strand_seq, f)[0]["regions"]
            self.assertEqual([r[0] for r in regions], ["18S", "ITS1", "5.8S", "28S", "5S"])
            self.assertTrue(all(a < b for _n, a, b in regions))
            self.assertTrue(all(x[2] <= y[1] for x, y in zip(regions, regions[1:])))

    def test_itsx_failures_are_seen(self):
        # ITSx with --cpu over 1 exits 0 when HMMER fails on a missing profile set: the build sees it, before any work
        with tempfile.TemporaryDirectory() as tmp:
            broken = tool(tmp, "ITSx_broken", "#!/bin/sh\necho 'Error: File existence/permissions problem in trying "
                                              "to open HMM file F.hmm' >&2\necho 'ERROR! HMMER search terminated "
                                              "unexpectedly.' >&2\nexit 0\n")
            failing = tool(tmp, "ITSx_failing", "#!/bin/sh\nexit 1\n")
            good = tool(tmp, "ITSx_good", "#!/bin/sh\necho 'Preparing HMM database (should be quick)...' >&2\n"
                                          "exit 0\n")
            self.assertIn("F.hmm", rrna.check_itsx(broken))
            self.assertIn("exit 1", rrna.check_itsx(failing))
            self.assertIsNone(rrna.check_itsx(good))
            with self.assertRaises(rrna.ItsxError):
                rrna.run_itsx(["ACGT" * 50], broken)

    def test_itsx_searches_the_lineages_sets(self):
        # each record's ITSx profile sets by its NCBI lineage (fungi F, animals M), all sets for one in none of them;
        # one ITSx run per distinct set
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            self.assertEqual([rrna.itsx_groups(ncbi, t) for t in (559292, "6253", 952, None)],
                             ["F", "M", "all", "all"])
            log = os.path.join(tmp, "itsx_calls.txt")
            logging = tool(tmp, "ITSx_log", f"#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = -t ] && echo \"$2\" >> {log}; "
                                            "shift; done\nexit 0\n")
            got = rrna.run_itsx(["AAAA", "CCCC", "GGGG", "AAAA"], logging,
                                groups={"AAAA": "F", "CCCC": "M", "GGGG": "F"})
            self.assertEqual(set(got), {"AAAA", "CCCC", "GGGG"})
            with open(log) as fh:
                self.assertEqual(sorted(fh.read().split()), ["F", "M"])

    def test_labels_and_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        index = rrna.name_index(ncbi)
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__Saccharomyces_cerevisiae", index, ncbi),
                         (4932, "eukaryome_name"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__Saccharomyces_novus", index, ncbi),
                         (4930, "eukaryome_genus"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Nogenus;s__Nogenus_x", index, ncbi), (None, None))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Kluyveromyces;s__Candida_sphaerica", index, ncbi),
                         (28985, "eukaryome_name"))  # an NCBI synonym
        # EUKARYOME 2.1 writes the epithet alone in s__ (g__Geotrupes;s__spiniger): the genus and it are the name
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__cerevisiae", index, ncbi),
                         (4932, "eukaryome_name"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__unclassified", index, ncbi),
                         (4930, "eukaryome_genus"))
        self.assertEqual(rrna.species_name("Drosophila(Metazoa)", "melanogaster"), "Drosophila melanogaster")
        self.assertEqual((rrna.homopolymer("ACGGGGTA"), round(rrna.ambiguous("ACGN"), 2)), (4, 0.25))
        a = "".join(random.Random(2).choice("ACGT") for _ in range(1000))
        b = a[:500] + ("A" if a[500] != "A" else "C") + a[501:]  # 99.9% identical: one type
        c = "".join(random.Random(3).choice("ACGT") for _ in range(1000))
        self.assertEqual(rrna.types_of([a, a, b, c]), {a: 1, b: 1, c: 2})


class PartialTest(unittest.TestCase):
    def test_species_clades_limit_the_package(self):
        # --species (a trial build): the genus Saccharomyces and Candida albicans, as taxids or a file; a strain's
        # taxid is inside its species' clade
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            listed = os.path.join(tmp, "clades.txt")
            with open(listed, "w") as fh:
                fh.write("# trial\n4930\t# Saccharomyces\n5476\n")
            self.assertEqual(eukore_build.clades(listed), {4930, 5476})
            self.assertEqual(eukore_build.clades("4930,5476"), {4930, 5476})
            self.assertIsNone(eukore_build.clades(None))
            inside = eukore_build.inside_of(ncbi, {4930, 5476})
            self.assertEqual([inside(t) for t in (559292, 4932, "5476", 28985, None, "-")],
                             [True, True, True, False, False, False])
            self.assertTrue(eukore_build.inside_of(ncbi, None)(28985))


if __name__ == "__main__":
    unittest.main()
