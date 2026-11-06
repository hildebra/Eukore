#!/usr/bin/env python3
"""Tests of EUKDB's core modules: the NCBI taxonomy and its subset (written into a package and read back the same for
the package's taxa), NCBI's assembly tables read, the model species table.

    python3 -m unittest discover -s tests -t tests
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukdb import ncbi, sources, species  # noqa: E402
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


if __name__ == "__main__":
    unittest.main()
