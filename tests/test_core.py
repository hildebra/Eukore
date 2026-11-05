#!/usr/bin/env python3
"""Tests of EUKDB's core modules: the NCBI taxonomy and its subset (written into a package and read back the same for
the package's taxa), NCBI's assembly tables read, the model species table; and that EUKDBbuild imports nothing
from protal's scripts/.

    python3 -m unittest discover EUKDBbuild/tests
"""

import ast
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
        self.assertEqual(len(sources.read_loci_table()), 30)

    def test_self_contained(self):
        # EUKDBbuild imports nothing from protal's scripts/ (its modules import each other relatively); outside
        # protal's repository there is nothing to check
        if not os.path.isdir(os.path.join(os.path.dirname(os.path.dirname(HERE)), "scripts", "mini_db")):
            self.skipTest("not inside protal's repository")
        protal ={os.path.splitext(f)[0] for d in ("scripts", os.path.join("scripts", "mini_db"))
                  for f in os.listdir(os.path.join(os.path.dirname(os.path.dirname(HERE)), d)) if f.endswith(".py")}
        folder = os.path.join(os.path.dirname(HERE), "eukdb")
        for name in os.listdir(folder):
            if not name.endswith(".py"):
                continue
            with open(os.path.join(folder, name)) as fh:
                tree = ast.parse(fh.read())
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    bad = [a.name for a in node.names if a.name.split(".")[0] in protal]
                elif isinstance(node, ast.ImportFrom) and node.level == 0:
                    bad = [node.module] if (node.module or "").split(".")[0] in protal else []
                else:
                    continue
                self.assertEqual(bad, [], f"eukdb/{name} imports {bad} from protal's scripts/")


if __name__ == "__main__":
    unittest.main()
