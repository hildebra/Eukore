#!/usr/bin/env python3
"""Tests of EUKDB's core modules: the NCBI taxonomy and its subset (written into a package and read back the same for
the package's taxa), the EukDetect readers on the test world, the model species table; and that EUKDBbuild imports
nothing from protal's scripts/.

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
            eukdetect, *_ = fx.world(tmp)
            full = ncbi.Ncbi(os.path.join(eukdetect, "taxdump.tar.gz"))
            self.assertEqual(full.species_of(559292), 4932)
            n = full.write_subset(os.path.join(tmp, "taxonomy"), [559292, 28985])
            sub = ncbi.Ncbi(os.path.join(tmp, "taxonomy"))
            self.assertEqual(n, len(set(full.ancestors(559292)) | set(full.ancestors(28985))))
            for taxid in (559292, 4932, 28985):
                self.assertEqual((sub.species_of(taxid), sub.ranks(taxid)), (full.species_of(taxid), full.ranks(taxid)))
            self.assertEqual(sub.synonyms[28985], ["Candida sphaerica"])
            self.assertNotIn(40410, sub.name)  # not on the lineages asked for

    def test_eukdetect_readers(self):
        with tempfile.TemporaryDirectory() as tmp:
            eukdetect, _ksgp, assemblies, _gtdb = fx.world(tmp)
            eukdetect1 = fx.eukdetect1_world(tmp, eukdetect)
            genomes, records = sources.read_link(os.path.join(eukdetect, "busco_taxid_genome_link.txt"))
            copies = sources.read_loci(os.path.join(eukdetect, "eukdb.fasta"), records)
            self.assertEqual(sorted(copies["Saccharomyces_cerevisiae_S288C"]), ["1at2759", "2at2759", "3at2759"])
            self.assertEqual(len(copies["Candida_albicans_SC5314"]["1at2759"]), 400)  # the longer of its two pieces
            genomes1, records1 = sources.read_link1(os.path.join(eukdetect1, "busco_taxid_link.txt"), fx.MAP1)
            self.assertIn(("Ascaris_lumbricoides", "1at2759"), set(records1.values()))
            self.assertNotIn(("Ascaris_lumbricoides", "3at2759"), set(records1.values()))  # duplicated: left out
            rows = sources.read_assemblies([assemblies])
            self.assertNotIn("GCA_000000560.1", {r["accession"] for r in rows})  # a gene's record, no genome
            taxonomy = ncbi.Ncbi(os.path.join(eukdetect, "taxdump.tar.gz"))
            have, _link = sources.species_loci(eukdetect, eukdetect1, fx.MAP1, taxonomy, {"1at2759", "2at2759"})
            self.assertEqual(have[6253], {"1at2759", "2at2759"})

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
