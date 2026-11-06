#!/usr/bin/env python3
"""Tests of Eukore's core modules: the NCBI taxonomy and its subset (written into a package and read back the same for
the package's taxa), NCBI's assembly tables read, the model species table.

    python3 -m unittest discover -s tests -t tests
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi, sources, species, util  # noqa: E402
import fixtures as fx  # noqa: E402


class CoreTest(unittest.TestCase):
    def test_taxonomy_subset(self):
        with tempfile.TemporaryDirectory() as tmp:
            full = ncbi.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            self.assertEqual(full.species_of(559292), 4932)
            n = full.write_subset(os.path.join(tmp, "taxonomy"), [559292, 28985])
            sub = ncbi.Ncbi(os.path.join(tmp, "taxonomy"))
            self.assertEqual(n, len(set(full.ancestors(559292)) | set(full.ancestors(28985))))
            for taxid in (559292, 4932, 28985):
                self.assertEqual((sub.species_of(taxid), sub.ranks(taxid)), (full.species_of(taxid), full.ranks(taxid)))
            self.assertEqual(sub.synonyms[28985], ["Candida sphaerica"])
            self.assertNotIn(40410, sub.name)  # not on the lineages asked for
            # the genetic code, NCBI's nodes.dmp column 7, kept in the subset
            self.assertEqual((full.code_of(28985), sub.code_of(28985), sub.code_of(4932)), (6, 6, 1))

    def test_assembly_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            by = {r["accession"]: r for r in rows}
            self.assertNotIn("GCA_000000560.1", by)  # a gene's record as a "complete genome": no genome
            self.assertEqual(by["GCA_000000011.1"]["flags"], ["contaminated"])
            self.assertEqual(by["GCA_000000012.1"]["category"], "derived from metagenome")
            self.assertEqual((by["GCA_000146045.2"]["refseq_category"], by["GCA_000146045.2"]["paired"],
                              by["GCA_000146045.2"]["group"]), ("reference genome", "GCF_000146045.2", "fungi"))

    def test_the_model_species_table(self):
        rows = species.read_table()
        self.assertIn("Homo sapiens", {r["species"] for r in species.eukaryotes(rows)})
        self.assertTrue(os.path.isfile(species.UNITS))
        self.assertFalse(os.path.exists(os.path.join(sources.DATA, "euk30_loci.tsv")))  # the user picks loci


class ProgressTest(unittest.TestCase):
    def test_reports_every_interval_with_the_time_left(self):
        # four genomes of 1, 1, 2 and 4 Gb; a report every 600 s; a Gb takes 100 s
        now, said = [0.0], []
        p = util.Progress("genomes fetched", [1e9, 1e9, 2e9, 4e9], said.append, every=600, clock=lambda: now[0])
        now[0] = 100.0
        p.step(1e9)
        self.assertEqual(said, [])  # not 600 s yet
        now[0] = 700.0  # a slow one
        p.step(1e9, failed=True)
        self.assertEqual(len(said), 1)
        self.assertIn("2/4 genomes fetched (1 failed), 2 left (6.0 Gb)", said[0])
        self.assertIn("11m 40s so far, ~35m 00s left", said[0])  # 2 Gb in 700 s: 6 Gb in 2,100 s, by size not count
        now[0] = 900.0
        p.step(2e9, skipped=True)  # not fetched yet: no work, out of the estimate
        self.assertEqual(len(said), 1)
        now[0] = 1000.0
        p.step(4e9, note="; done")
        self.assertEqual(said[1], "  4/4 genomes fetched (1 failed), 0 left; took 16m 40s; done")
        self.assertEqual([util.duration(s) for s in (45, 725, 11220, 187200)], ["45s", "12m 05s", "3h 07m", "2d 04h"])
        # no sizes: each item weighs the same
        q = util.Progress("genomes run", [None, None], said.append, every=0, clock=lambda: now[0])
        now[0] += 10
        q.step(None)
        self.assertIn("~10s left", said[-1])


if __name__ == "__main__":
    unittest.main()
