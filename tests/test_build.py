#!/usr/bin/env python3
"""Tests of Eukore's build: an rDNA unit uncut, with its 5.8S and an inline 5S (one block) and not a 5S elsewhere, a
tandem array's spacer 5S one unit's, a 28S cut by an intron one; a published read's regions on its record; a
EUKARYOME read's NCBI taxid by its species name, else its genus's; the record checks; ITSx's failures seen; the
build's resume (a cut cache line, a killed run's old package) and a genome's unexpected error its own failure.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi as ncbi_mod, rrna  # noqa: E402
import fixtures as fx  # noqa: E402

tool = fx.tool


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

    def test_a_tandem_arrays_spacer_5s_is_one_units(self):
        # two units 9 kb apart, a 5S in the spacer: within 3 kb of the first's 28S end and of the second's 18S start;
        # it is the first's (the 28S it follows), the second has none (until 2026-11-06 both took it)
        rng = random.Random(4)
        seq = "".join(rng.choice("ACGT") for _ in range(40000))
        feats = {"18S_rRNA": [(1001, 2800, "+", False), (15001, 16800, "+", False)],
                 "5_8S_rRNA": [(3201, 3350, "+", False), (17201, 17350, "+", False)],
                 "28S_rRNA": [(3601, 12100, "+", False), (17601, 26100, "+", False)],
                 "5S_rRNA": [(13000, 13120, "+", False)]}
        first, second = rrna.units_of(seq, feats)
        self.assertEqual([r[0] for r in first["regions"]][-1], "5S")
        self.assertNotIn("5S", [r[0] for r in second["regions"]])
        self.assertEqual(first["unit_span"], (1001, 13120))

    def test_a_28s_cut_by_an_intron_is_one(self):
        # barrnap's 28S in two partial pieces (a group I intron, Tetrahymena's): the unit reaches the second's end,
        # whole where their shares make 80% (until 2026-11-06 it ended at the first piece)
        rng = random.Random(5)
        seq = "".join(rng.choice("ACGT") for _ in range(20000))
        feats = {"18S_rRNA": [(1001, 2800, "+", False, 1.0)], "5_8S_rRNA": [(3201, 3350, "+", False, 1.0)],
                 "28S_rRNA": [(3601, 6000, "+", True, 0.45), (6400, 9000, "+", True, 0.5)]}
        unit = rrna.units_of(seq, feats)[0]
        self.assertEqual((unit["unit_span"], unit["lsu_complete"]), ((1001, 9000), True))
        short = dict(feats, **{"28S_rRNA": [(3601, 6000, "+", True, 0.45)]})
        self.assertFalse(rrna.units_of(seq, short)[0]["lsu_complete"])
        gff = ("w\tbarrnap:0.9\trRNA\t1\t900\t1e-9\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA (partial);"
               "note=aligned only 39 percent of the 28S ribosomal RNA\n")
        self.assertEqual(rrna.parse_gff(gff)["w"]["28S_rRNA"], [(1, 900, "+", True, 0.39)])

    def test_a_published_reads_regions_on_its_record(self):
        # a read whose unit starts at its base 150 (+), or lies on its reverse strand: the regions on the record, the
        # read itself or its reverse complement (until 2026-11-06 the unit's, 149 bases early or mirrored)
        unit = dict(unit_span=(150, 1000), span=(150, 900, "+"), regions=[("18S", 0, 751), ("ITS1", 751, 851)])
        self.assertEqual(rrna.on_record(unit, 1200), [("18S", 149, 900), ("ITS1", 900, 1000)])
        rev = dict(unit_span=(201, 1051), span=(300, 1051, "-"), regions=[("18S", 0, 752)])
        self.assertEqual(rrna.on_record(rev, 1200), [("18S", 149, 901)])  # the record: the read's reverse complement
        rng = random.Random(6)
        read = "".join(rng.choice("ACGT") for _ in range(1200))
        record = rrna.revcomp(read)
        self.assertEqual(record[149:901], rrna.revcomp(read[300 - 1:1051]))  # the 18S, read on its strand


class HomonymTest(unittest.TestCase):
    def test_a_genus_of_two_kingdoms(self):
        from eukore import labels  # noqa: PLC0415
        # EUKARYOME's classification names Drosophila twice (a fly, a fungus); a long read of the fly's takes the
        # fly's lineage by its kingdom-specific name, and a plain-name read none of either (its own k__ .. g__)
        rows = [("superkingdom",) * 10,
                ("Eukaryota", "Fungi", "", "Ascomycota", "", "Sordariomycetes", "Hypocreales", "Fam", "Drosophila",
                 "Drosophila(Fungi)"),
                ("Eukaryota", "Metazoa", "", "Arthropoda", "", "Insecta", "Diptera", "Drosophilidae", "Drosophila",
                 "Drosophila(Metazoa)")]
        ref = labels.read_eukaryome(rows, [("EUK1", "k__Metazoa;p__Arthropoda;g__Drosophila(Metazoa);s__melanogaster"),
                                           ("EUK2", "k__Metazoa;p__Arthropoda;g__Drosophila;s__simulans")])
        fly = next(iter(ref.species["drosophila melanogaster"]))
        self.assertEqual(fly[1], "Metazoa")
        self.assertEqual(next(iter(ref.species["drosophila simulans"]))[1], "Metazoa")  # its own k__
        # a genus match of another kingdom's lineage alone: none
        self.assertTrue(labels.other_kingdom(("Eukaryota", "Archaeplastida", "Streptophyta", "Morus"),
                                             {"eukaryota", "metazoa", "chordata", "morus"}))
        self.assertFalse(labels.other_kingdom(("Eukaryota", "Obazoa", "Metazoa", "Morus"),
                                              {"eukaryota", "metazoa", "chordata", "morus"}))


class ResumeTest(unittest.TestCase):
    def test_a_cut_cache_line_is_annotated_again(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            cache = os.path.join(tmp, "regions.tsv")
            with open(cache, "w", newline="\n") as fh:
                fh.write("EUK_1\t18S:0-1800\tt\t+\nEUK_2\t18S:0-1800,ITS1:1800-2\tf\t-\nEUK_3\t18S:1-")
            self.assertEqual(eukore_build.read_regions(cache), {"EUK_1": ("18S:0-1800", "t", "+"),
                                                                "EUK_2": ("18S:0-1800,ITS1:1800-2", "f", "-")})
            self.assertFalse(eukore_build.ends_in_newline(cache))

    def test_a_genomes_unexpected_error_is_its_failure(self):
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with mock.patch.object(eukore_build.genomes, "run", side_effect=UnicodeDecodeError("utf-8", b"\xff", 0, 1,
                                                                                         "invalid start byte")):
            acc, summary, why = eukore_build.run_one(({"accession": "GCA_1.1"}, {}))
        self.assertEqual((acc, summary), ("GCA_1.1", None))
        self.assertTrue(why.startswith("UnicodeDecodeError"))

    def test_a_package_written_over_a_killed_runs_old(self):
        from eukore import package  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "Eukore_2026.11")
            for folder in (out, out + ".old"):  # a package, and a killed run's .old left beside it
                os.makedirs(os.path.join(folder, "busco"))
                with open(os.path.join(folder, "busco", "x.tsv"), "w") as fh:
                    fh.write("x\n")
            part = os.path.join(tmp, "new")
            os.makedirs(part)
            with open(os.path.join(part, "manifest.json"), "w") as fh:
                fh.write("{}\n")
            self.assertEqual(package.put_in_place(part, out), out)
            self.assertEqual(os.listdir(out), ["manifest.json"])
            self.assertEqual(sorted(os.listdir(tmp)), ["Eukore_2026.11"])

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

    def test_itsx_guessed_borders_are_not_taken(self):
        # ITSx without the 5.8S's end still writes an ITS2, 45 bases after the 5.8S start (Arabidopsis with set T,
        # 2026-11-06): ITS1 and the 28S start are ITSx's, the 5.8S's end barrnap's
        line = "u0\t3755 bp.\tSSU: 1-300\tITS1: 301-590\t5.8S: No end\tITS2: 636-1206\tLSU: 1207-3755\t"
        with tempfile.TemporaryDirectory() as tmp:
            itsx = tool(tmp, "ITSx_noend", "#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = -o ] && o=$2; shift; done\n"
                                           f"printf '{line}\\n' > $o.positions.txt\nexit 0\n")
            found = rrna.run_itsx(["ACGT" * 50], itsx)["ACGT" * 50]
        self.assertEqual(found["5.8S problem"], "No end")
        self.assertEqual((found["ITS1"], found["ITS2"]), ((301, 590), (636, 1206)))

        def item():  # barrnap's regions on a unit, the ITSx stretch from 1500
            return dict(unit="A" * 5000, its_at=1500, regions=[("18S", 0, 1800), ("ITS1", 1800, 2100),
                                                               ("5.8S", 2100, 2260), ("ITS2", 2260, 2500),
                                                               ("28S", 2500, 5000)])
        got = {r[0]: r[1:] for r in rrna.refine(item(), found)["regions"]}
        self.assertEqual(got, {"18S": (0, 1800), "ITS1": (1800, 2090), "5.8S": (2090, 2260), "ITS2": (2260, 2706),
                               "28S": (2706, 5000)})
        whole = {k: v for k, v in found.items() if k != "5.8S problem"}  # what was taken before: a 45-base 5.8S
        self.assertEqual({r[0]: r[1:] for r in rrna.refine(item(), whole)["regions"]}["5.8S"], (2090, 2135))
        neither = dict(found, **{"5.8S problem": "Not found"})
        self.assertEqual({r[0]: r[1:] for r in rrna.refine(item(), neither)["regions"]}["5.8S"], (2100, 2260))
        start = dict(found, **{"5.8S problem": "No start"})
        self.assertEqual({r[0]: r[1:] for r in rrna.refine(item(), start)["regions"]}["5.8S"], (2100, 2135))

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
