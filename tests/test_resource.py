#!/usr/bin/env python3
"""Tests of the package described as a resource (benchmarks/eukore_bench/resource.py; since 2026-11-12): the taxa at
every rank NCBI names, holds a genome of and Eukore holds; the markers a species holds from one assembly; the named
lists' organisms found by their names.

    python3 -m unittest discover -s tests -t tests
"""

import importlib.util
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore import benchmarks, ncbi as ncbi_mod  # noqa: E402
from eukore_bench import resource  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST = "GCA_000146045.2"


def species_row(representative="-", group="Fungi", loci=0, rrna=0, organelle=0):
    return dict(representative=representative, supergroup=group, loci=str(loci), rrna_records=str(rrna),
                organelle_genes=str(organelle))


class ResourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(cls.tmp.name, "taxdump.tar.gz")))
        # the yeast with a genome, S. paradoxus with rRNA records alone, a published read labelled by its genus
        cls.species = {"4932": species_row(YEAST, loci=3, rrna=2, organelle=1), "27291": species_row(rrna=1),
                       "951": species_row(group="other eukaryotes", rrna=1)}
        cls.taxa = [dict(taxid="4932", supergroup="Fungi", fate="chosen"),
                    dict(taxid="28985", supergroup="Fungi", fate="chosen, flagged"),
                    dict(taxid="5476", supergroup="Fungi", fate="every assembly excluded"),
                    dict(taxid="1080349", supergroup="Fungi", fate="hybrid"),
                    dict(taxid="4930", supergroup="Fungi", fate="not a species")]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_lineages(self):
        lin = resource.Lineages(self.ncbi)
        self.assertEqual(lin.of(4932), ("Fungi", (4930, 4893, 4892, 4891, 4890)))
        self.assertEqual(lin.of(559292), lin.of(4932))  # a strain: its species' lineage
        self.assertEqual(lin.of(952), ("other eukaryotes", (951, None, None, None, 950)))
        self.assertEqual(lin.of(9606)[0], "Metazoa")
        self.assertIsNone(lin.of(1)[0])  # the root: no eukaryote
        self.assertIsNone(lin.of(123456789)[0])  # not in the taxonomy
        self.assertNotIn(4932, lin.cache)  # the taxa above the species alone are kept
        self.assertIn(4930, lin.cache)

    def test_named_species(self):
        self.assertTrue(resource.named_species(self.ncbi, 4932))
        self.assertFalse(resource.named_species(self.ncbi, 953))  # uncultured Mysterium
        self.assertFalse(resource.named_species(self.ncbi, 1080349))  # a hybrid

    def test_breadth_at_every_rank(self):
        lin = resource.Lineages(self.ncbi)
        rows = {(r["supergroup"], r["rank"]): r for r in resource.rank_rows(self.species, self.taxa, self.ncbi, lin)}
        # NCBI's named species: every binomial of the fixtures; with a genome: three (the hybrid and the genus none)
        self.assertEqual(rows[("all", "species")],
                         dict(supergroup="all", rank="species", ncbi_named=8, ncbi_genome=3, eukore_genome=1,
                              eukore=2))
        # genera: the yeast's with a genome; Mysterium by its published read
        self.assertEqual([rows[("all", "genus")][k] for k in resource.RANK_COLUMNS[2:]], [7, 3, 1, 2])
        self.assertEqual([rows[("Fungi", "family")][k] for k in resource.RANK_COLUMNS[2:]], [3, 2, 1, 1])
        self.assertEqual(rows[("other eukaryotes", "phylum")]["eukore"], 1)
        self.assertEqual(rows[("other eukaryotes", "species")]["eukore"], 0)  # the read names a genus
        # a package's own taxonomy, a package before selection_taxa.tsv: neither of NCBI's counted
        rows = resource.rank_rows(self.species, [], self.ncbi, lin, whole=False)
        self.assertEqual({(r["ncbi_named"], r["ncbi_genome"]) for r in rows}, {("-", "-")})

    def test_markers_of_one_assembly(self):
        folder = os.path.join(self.tmp.name, "package")
        os.makedirs(os.path.join(folder, "rrna"))
        os.makedirs(os.path.join(folder, "organelle"))
        with open(os.path.join(folder, "rrna", "records.tsv"), "w") as fh:
            fh.write("record\tkind\tspecies_taxid\taccession\tsource\tregions\n"
                     f"a\trDNA\t4932\t{YEAST}\tgenome\t18S:0-1800,ITS1:1800-2100,5.8S:2100-2260,ITS2:2260-2500,"
                     "28S:2500-5900\n"
                     "b\t18S\t27291\t-\tEUKARYOME\t-\n"                       # a published 18S: no assembly's
                     "c\trDNA\t27291\t-\tROD\t18S:0-1800,ITS1:1800-2100,5.8S:2100-2260\n"   # no ITS2: no ITS
                     "d\t18S\t999\t-\tROD\t-\n")                              # no species of the package
        with open(os.path.join(folder, "organelle", "genes.tsv"), "w") as fh:
            fh.write("record\tgene\torganelle\tspecies_taxid\tsource\taccession\n"
                     "m\tCOX1\tmitochondrion\t4932\trefseq\t-\n")             # RefSeq's organelle genome
        at = {(r["supergroup"], r["marker"]): (r["of"], r["species"], r["same_assembly"])
              for r in resource.content_rows(folder, self.species)}
        self.assertEqual(at[("all", "BUSCO")], (3, 1, 1))
        self.assertEqual(at[("all", "18S")], (3, 2, 1))
        self.assertEqual(at[("all", "ITS")], (3, 1, 1))
        self.assertEqual(at[("all", "operon")], (3, 2, 1))  # a unit's record, whatever its regions
        self.assertEqual(at[("all", "mitochondrion")], (3, 1, 0))
        self.assertEqual(at[("all", "BUSCO + operon")], (3, 1, 1))
        self.assertEqual(at[("all", "BUSCO + operon + organelle")], (3, 1, 0))  # its organelle genes RefSeq's
        self.assertEqual(at[("other eukaryotes", "18S")], (1, 0, 0))

    def test_named_lists(self):
        lin = resource.Lineages(self.ncbi)
        self.assertEqual(resource.resolve(["Saccharomyces cerevisiae", "candida_sphaerica", "Nonexistus"], self.ncbi,
                                          lin), {"Saccharomyces cerevisiae": 4932, "candida_sphaerica": 28985})
        self.assertEqual(resource.name_key("[Candida] auris"), "candida auris")
        entries = [("p", "Saccharomyces cerevisiae", None, ""), ("p", "Saccharomyces", None, "a genus"),
                   ("p", "Candida sphaerica", None, ""), ("p", "Saccharomyces paradoxus", None, ""),
                   ("q", "Nonexistus", None, ""), ("q", "by its taxid", 5476, "")]
        rows = {r["entry"]: r for r in resource.panel_rows(entries, self.species, self.taxa, self.ncbi, lin)}
        self.assertEqual([rows[e]["status"] for _p, e, _t, _n in entries],
                         ["genome", "genome", "absent", "records", "not found", "absent"])
        genus = rows["Saccharomyces"]
        self.assertEqual((genus["rank"], genus["eukore_genome"], genus["eukore"], genus["ncbi_genome"], genus["loci"],
                          genus["rrna_records"], genus["note"]), ("genus", 1, 2, 1, 3, 3, "a genus"))
        self.assertEqual((rows["Candida sphaerica"]["ncbi_name"], rows["Candida sphaerica"]["ncbi_genome"]),
                         ("Kluyveromyces lactis", 1))
        self.assertEqual(rows["by its taxid"]["ncbi_name"], "Candida albicans")
        self.assertEqual(resource.panel_summary(rows.values())["p"],
                         {"entries": 4, "genome": 2, "records": 1, "absent": 1, "not found": 0})

    def test_the_lists_file(self):
        entries = resource.read_panels()
        self.assertEqual(len([e for e in entries if e[0] == "WHO fungal priority pathogens"]), 19)
        self.assertTrue(all(name and panel for panel, name, _t, _n in entries))
        self.assertEqual(len({(p, n) for p, n, _t, _n in entries}), len(entries))  # no entry twice in a list

    def test_run_and_figures(self):
        package = os.path.join(self.tmp.name, "run_package")
        out = os.path.join(self.tmp.name, "results")
        os.makedirs(out)
        benchmarks.write(out, "selection_taxa.tsv", benchmarks.TAXA_COLUMNS,
                         [dict(r, name="-", rank="species") for r in self.taxa])
        summary = {"species": 1, "busco_loci": 3, "left_out": {}}
        resource.run(package, out, self.ncbi, self.species, summary, log=lambda t: None)
        self.assertEqual(summary["breadth_ranks"]["genus"],
                         {"ncbi_named": 7, "ncbi_genome": 3, "eukore_genome": 1, "eukore": 2})
        self.assertEqual(summary["content"]["BUSCO"], {"species": 1, "same_assembly": 1})
        self.assertEqual(summary["panels"]["mock communities"]["genome"], 1)  # the yeast
        for name in ("breadth_ranks.tsv", "content.tsv", "panels.tsv"):
            self.assertTrue(os.path.isfile(os.path.join(out, name)), name)
        md = benchmarks.summary_md(summary, "Eukore_x")
        self.assertIn("| genera with a genome; with any record | 1 of NCBI's 3 with a genome (7 named); 2 |", md)
        self.assertIn("| phyla with a genome; with any record |", md)
        self.assertIn("BUSCO 1 (1)", md)
        self.assertIn("| mock communities: entries with a genome | 1 of ", md)
        if importlib.util.find_spec("matplotlib") is None:
            return
        from eukore_bench import figures  # noqa: PLC0415
        figs = os.path.join(self.tmp.name, "figures")
        os.makedirs(figs)
        for fn in (figures.breadth_ranks, figures.content, figures.panels):
            self.assertTrue(os.path.isfile(fn(out, package, figs, "png")))

    def test_taxa_of_the_selection(self):
        rows = benchmarks.taxa_rows({4932: "chosen", 1080349: "hybrid"}, self.ncbi)
        self.assertEqual(rows[0], dict(taxid=4932, name="Saccharomyces cerevisiae", rank="species",
                                       supergroup="Fungi", fate="chosen", assemblies="-", reasons="-"))
        self.assertEqual(rows[1]["fate"], "hybrid")


if __name__ == "__main__":
    unittest.main()
