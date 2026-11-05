#!/usr/bin/env python3
"""Tests of EUKDB's build: an rDNA unit uncut, with its 5.8S and an inline 5S (one block) and not a 5S elsewhere; a
EUKARYOME read's NCBI taxid by its species name, else its genus's; the record checks; a model species' genome run
(stand-in miniprot and barrnap: its locus placed whole over 16 kb, its operon over 8,191 bases whole) and the package
written from it and the test world, read back and checked.

    python3 -m unittest discover -s EUKDBbuild/tests -t EUKDBbuild/tests
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
from eukdb import genomes, ncbi as ncbi_mod, package, rrna, sources, species  # noqa: E402
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
            eukdetect = fx.world(tmp)[0]
            fx.eukdetect1_world(tmp, eukdetect)
            ncbi = ncbi_mod.Ncbi(os.path.join(eukdetect, "taxdump.tar.gz"))
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
            eukdetect, *_ = fx.world(tmp)
            ncbi = ncbi_mod.Ncbi(os.path.join(eukdetect, "taxdump.tar.gz"))
        index = rrna.name_index(ncbi)
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__Saccharomyces_cerevisiae", index, ncbi),
                         (4932, "eukaryome_name"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__Saccharomyces_novus", index, ncbi),
                         (4930, "eukaryome_genus"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Nogenus;s__Nogenus_x", index, ncbi), (None, None))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Kluyveromyces;s__Candida_sphaerica", index, ncbi),
                         (28985, "eukaryome_name"))  # an NCBI synonym
        self.assertEqual((rrna.homopolymer("ACGGGGTA"), round(rrna.ambiguous("ACGN"), 2)), (4, 0.25))
        a = "".join(random.Random(2).choice("ACGT") for _ in range(1000))
        b = a[:500] + ("A" if a[500] != "A" else "C") + a[501:]  # 99.9% identical: one type
        c = "".join(random.Random(3).choice("ACGT") for _ in range(1000))
        self.assertEqual(rrna.types_of([a, a, b, c]), {a: 1, b: 1, c: 2})


class PartialTest(unittest.TestCase):
    def test_species_clades_limit_the_package(self):
        # --species (a trial build): the genus Saccharomyces and Candida albicans, as taxids or a file; a strain's
        # taxid is inside its species' clade; the collected genomes, copies, features and taxids of the rest left out
        import eukdb_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            eukdetect = fx.world(tmp)[0]
            ncbi = ncbi_mod.Ncbi(os.path.join(eukdetect, "taxdump.tar.gz"))
            listed = os.path.join(tmp, "clades.txt")
            with open(listed, "w") as fh:
                fh.write("# trial\n4930\t# Saccharomyces\n5476\n")
            self.assertEqual(eukdb_build.clades(listed), {4930, 5476})
            self.assertEqual(eukdb_build.clades("4930,5476"), {4930, 5476})
            self.assertIsNone(eukdb_build.clades(None))
            inside = eukdb_build.inside_of(ncbi, {4930, 5476})
            self.assertEqual([inside(t) for t in (559292, 4932, "5476", 28985, None, "-")],
                             [True, True, True, False, False, False])
            self.assertTrue(eukdb_build.inside_of(ncbi, None)(28985))
            collected = dict(genomes={"A": dict(species_taxid=4932), "B": dict(species_taxid=28985)},
                             copies={"A": {"1at2759": "ACGT"}, "B": {"1at2759": "ACGT"}},
                             features=[dict(accession="A"), dict(accession="B")], taxids={559292, 28985})
            got = eukdb_build.subset(collected, inside)
            self.assertEqual((set(got["genomes"]), set(got["copies"]), [f["accession"] for f in got["features"]],
                              got["taxids"]), ({"A"}, {"A"}, ["A"], {559292}))


class RunTest(unittest.TestCase):
    def test_a_model_species_genome_into_the_package(self):
        rng = random.Random(4)
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(fx, "NODES", fx.NODES + VERTEBRATES):
            eukdetect, _ksgp, assemblies_path, _gtdb = fx.world(tmp)
            eukdetect1 = fx.eukdetect1_world(tmp, eukdetect)
            ncbi = ncbi_mod.Ncbi(os.path.join(eukdetect, "taxdump.tar.gz"))
            # the human genome: an rDNA unit (18S motif at 1,500) and a gene of 18,001 bases with three exons
            c1 = list("".join(rng.choice("ACGT") for _ in range(60000)))
            motif = "GTGCCAGCAGCCGCGGTAA"
            c1[1500:1500 + len(motif)] = motif
            c1 = "".join(c1)
            genomes_dir = os.path.join(tmp, "out", "genomes")
            os.makedirs(genomes_dir)
            with gzip.open(genomes.genome_file(genomes_dir, HUMAN), "wt") as fh:
                fh.write(f">chr1\n{c1}\n")
            q = "model-Homo_sapiens-1at2759-S1~0~Ascaris"
            gff = ["##gff-version 3", f"chr1\tminiprot\tmRNA\t30001\t48001\t900\t+\t.\tID=M1;Identity=0.6000;Target={q} 1 99",
                   f"chr1\tminiprot\tCDS\t30001\t30100\t9\t+\t0\tParent=M1;Target={q} 1 33",
                   f"chr1\tminiprot\tCDS\t40001\t40100\t9\t+\t0\tParent=M1;Target={q} 34 66",
                   f"chr1\tminiprot\tCDS\t47902\t48001\t9\t+\t0\tParent=M1;Target={q} 67 99"]
            with open(os.path.join(tmp, "mp.gff"), "w") as fh:
                fh.write("\n".join(gff) + "\n")
            miniprot = tool(tmp, "miniprot", f"#!/bin/sh\ncat '{tmp}/mp.gff'\n")
            bar = ["##gff-version 3", "w0\tbarrnap\trRNA\t1001\t2800\t0\t+\t.\tName=18S_rRNA;product=18S",
                   "w0\tbarrnap\trRNA\t3201\t3350\t0\t+\t.\tName=5_8S_rRNA;product=5.8S",
                   "w0\tbarrnap\trRNA\t3601\t12100\t0\t+\t.\tName=28S_rRNA;product=28S",
                   "w0\tbarrnap\trRNA\t12500\t12620\t0\t-\t.\tName=5S_rRNA;product=5S"]
            with open(os.path.join(tmp, "bar.gff"), "w") as fh:
                fh.write("\n".join(bar) + "\n")
            barrnap = tool(tmp, "barrnap", f"#!/bin/sh\ncat '{tmp}/bar.gff'\n")
            row = dict(next(r for r in species.read_table() if r["species"] == "Homo sapiens"), why="model")
            job = dict(accession=HUMAN, url="https://x/y", size=3e9, species_taxid=9606, why=["model"],
                       queries={q: "ATG" + "GCT" * 98}, own={q: False}, rdna_unit="-", name="model_Homo_sapiens")
            work = os.path.join(tmp, "work")
            import eukdb_build  # noqa: PLC0415
            import functools  # noqa: PLC0415
            settings = dict(work=work, genomes=genomes_dir, miniprot=miniprot, threads=1, staged_only=True,
                            rrna=functools.partial(eukdb_build.rrna_hook, barrnap=barrnap, itsx=None, threads=1))
            acc, summary, why = genomes.run((job, settings))
            self.assertEqual((acc, why, summary["copies_taken"], summary["rrna"]), (HUMAN, None, 1, {"18S": 1, "rDNA": 1}))
            self.assertTrue(genomes.is_done(work, HUMAN))
            # not fetched and staged only: waits
            self.assertEqual(genomes.run((dict(job, accession="GCA_000000001.1"), settings))[2], "not fetched yet")
            # the package
            assemblies = sources.read_assemblies([assemblies_path])
            src = {"eukdetect2": eukdetect, "eukdetect1": eukdetect1}
            loci = ["1at2759", "2at2759", "3at2759"]
            loci_table = os.path.join(tmp, "loci.tsv")
            with open(loci_table, "w") as fh:
                fh.write("locus\tfungi_share\tprotist_share\tmedian_length\teukdetect1_locus\n")
                fh.writelines(f"{a}\t1\t1\t400\t{b}\n" for a, b in fx.MAP1.items())
            collected = package.collect(ncbi, src, assemblies, loci, fx.MAP1, [row], work, [], [])
            records = rrna.finish(rrna.genome_records(HUMAN, rrna.load(genomes.done_dir(work, HUMAN)), 9606, 9606),
                                  ncbi)
            out = package.write(os.path.join(tmp, "out", "EUKDB_2026.11"), "2026.11", ncbi, collected, loci_table,
                                records, genomes_dir, {"test": True})
            self.assertIsNone(package.check(out, deep=True))
            g = {r["accession"]: r for r in package.read_table(out, "genomes.tsv")}
            self.assertEqual((g[HUMAN]["source"], g[HUMAN]["representative"], g[HUMAN]["species_taxid"],
                              g[HUMAN]["local"]), ("model", "t", "9606", f"genomes/{HUMAN}.fna.gz"))
            self.assertEqual(g["GCA_000146045.2"]["source"], "eukdetect2")  # S288C's assembly
            human = package.read_copies(out, "1at2759")[HUMAN]
            self.assertEqual(human, c1[30000:48001])  # whole: 18,001 bases (protal fits it into 16 kb, not EUKDB)
            feats = [f for f in package.read_table(out, "euk30/features.tsv") if f["accession"] == HUMAN]
            self.assertEqual([(f["kind"], f["start"], f["end"], f["method"]) for f in feats],
                             [("exon", "0", "100", "placed_relative"), ("intron", "100", "10000", "placed_relative"),
                              ("exon", "10000", "10100", "placed_relative"),
                              ("intron", "10100", "17901", "placed_relative"),
                              ("exon", "17901", "18001", "placed_relative")])
            coding = [f for f in package.read_table(out, "euk30/features.tsv") if f["method"] == "coding"]
            self.assertTrue(coding)  # EukDetect1's Ascaris and Mysterium rarum: coding copies, one exon each
            recs = {r["record"]: r for r in package.read_table(out, "rrna/records.tsv")}
            unit = recs[f"{HUMAN}_rDNA_1"]
            self.assertEqual((unit["kind"], unit["length"], unit["species_taxid"], unit["label_from"]),
                             ("rDNA", str(12620 - 1001 + 1), "9606", "assembly"))
            self.assertIn("5S:", unit["regions"])
            sp = {r["taxid"]: r for r in package.read_table(out, "species.tsv")}
            self.assertEqual((sp["9606"]["representative"], sp["9606"]["rrna_records"]), (HUMAN, "2"))
            tax = ncbi_mod.Ncbi(os.path.join(out, "taxonomy"))
            self.assertEqual(tax.ranks(9606)["phylum"], "Chordata")
            with open(os.path.join(out, "manifest.json")) as fh:
                manifest = json.load(fh)
            self.assertEqual((manifest["format"], manifest["counts"]["rrna_records"]), (1, 2))
            # the prune: the model species' genome kept first
            size = os.path.getsize(genomes.genome_file(genomes_dir, HUMAN))
            self.assertEqual(genomes.prune(genomes_dir, work, size, set(), {HUMAN}), [(HUMAN, size, "kept")])
            self.assertEqual(genomes.prune(genomes_dir, work, 10, set(), {HUMAN}), [(HUMAN, size, "deleted")])
            self.assertFalse(os.path.exists(genomes.genome_file(genomes_dir, HUMAN)))


if __name__ == "__main__":
    unittest.main()
