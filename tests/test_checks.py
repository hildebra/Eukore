#!/usr/bin/env python3
"""Tests of the records' quality as headline numbers (benchmarks/eukore_bench/checks.py and tree.py; since
2026-11-12): each check's share per group and level, the genetic codes, the rDNA copies beside published numbers;
the tree of the BUSCO proteins against NCBI's taxonomy (the genomes and loci chosen, the proteins on their profiles,
a taxon as one branch of the unrooted tree, the genomes that stray).

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import importlib.util
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore import benchmarks, ncbi as ncbi_mod  # noqa: E402
from eukore_bench import checks, tree  # noqa: E402
import fixtures as fx  # noqa: E402

# the tree's genomes: two Saccharomyces, Kluyveromyces (their family), Candida (their order), a basidiomycete, a
# nematode, human; M of too few loci
SPECIES = {"Y": 4932, "P": 27291, "K": 28985, "C": 5476, "D": 40410, "A": 6253, "H": 9606, "M": 952}
RIGHT = "((((Y:1,P:1):1,K:1):1,C:1):1,D:1,(A:1,H:1):1);"
WRONG = "(((Y:1,P:1)0.99:1,C:1):1,D:1,((A:1,H:1):1,K:1):1);"   # Kluyveromyces beside the animals
FASTTREE = "#!/bin/sh\necho '%s'\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with (gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")) as fh:
        fh.write(text)


class ChecksTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.package = os.path.join(cls.tmp.name, "Eukore_x")
        head = "accession\tspecies_taxid\trepresentative\tsource\tflags\tforeign_copies\tduplication\tgencode_check" \
               "\tbest_table\trdna_units\trdna_copies_reads\n"
        write(os.path.join(cls.package, "genomes.tsv"), head
              + "Y\t4932\tt\tncbi\t-\t0\tlow\tagrees\t1\t2\t120.4\n"
              + "P\t27291\tt\tncbi\tFCS-GX: contaminated,excluded_from_refseq\t3\thaplotypes\tdiffers:CTG>S\t12\t1\t-\n"
              + "A\t6253\tt\tncbi\t-\t0\t-\tundecided\t-\t-\t-\n"
              + "N\t28985\tt\tno_assembly\t-\t0\t-\t-\t-\t-\t-\n"   # records of a species without an assembly
              + "Y2\t4932\tf\tassembly\tcontaminated\t0\t-\t-\t-\t-\t-\n")  # another assembly of the yeast
        write(os.path.join(cls.package, "busco", "copies.tsv"), "accession\tlocus\tconfidence\tflags\n"
              + "Y\t1at2759\thigh\t-\nY\t2at2759\tmedium\tfragmented,contig_edge\n"
              + "P\t1at2759\tlow\town_species\nA\t1at2759\tunverified\town_species_undecided\n"
              + "Y2\t1at2759\thigh\t-\n")
        write(os.path.join(cls.package, "rrna", "records.tsv"), "record\tspecies_taxid\tconfidence\tflags\n"
              + "r1\t4932\thigh\t-\nr2\t6253\tmedium\tlsu_partial\nr3\t999\tlow\tambiguous\n")
        cls.species = {"4932": dict(name="Saccharomyces cerevisiae", supergroup="Fungi", representative="Y"),
                       "27291": dict(name="Saccharomyces paradoxus", supergroup="Fungi", representative="P"),
                       "6253": dict(name="Ascaris lumbricoides", supergroup="Metazoa", representative="A"),
                       "28985": dict(name="Kluyveromyces lactis", supergroup="Fungi", representative="N")}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_each_checks_share(self):
        at = {(r["supergroup"], r["level"], r["check"]): (r["flagged"], r["of"])
              for r in checks.quality_rows(self.package, self.species)}
        # the genomes: the species' own, those of an assembly (no records' species, no further assembly)
        self.assertEqual(at[("all", "genome", "any flag")], (1, 3))
        self.assertEqual(at[("Fungi", "genome", "FCS-GX: contaminated")], (1, 2))
        self.assertEqual(at[("Fungi", "genome", "copies another clade's")], (1, 2))
        self.assertEqual(at[("Fungi", "genome", "duplication: haplotypes")], (1, 2))
        self.assertEqual(at[("all", "genome", "genetic code differs from NCBI's")], (1, 3))
        self.assertNotIn(("Metazoa", "genome", "any flag"), at)
        # the copies of the species' genomes, each flag and tier
        self.assertEqual(at[("all", "BUSCO copy", "confidence: low")], (1, 4))
        self.assertEqual(at[("Fungi", "BUSCO copy", "fragmented")], (1, 3))
        self.assertEqual(at[("Metazoa", "BUSCO copy", "own_species_undecided")], (1, 1))
        # the rRNA records of the package's species
        self.assertEqual(at[("all", "rRNA record", "lsu_partial")], (1, 2))
        self.assertNotIn(("all", "organelle gene", "confidence: high"), at)  # no table: no level

    def test_genetic_codes_and_rdna(self):
        rows = {r["supergroup"]: r for r in checks.gencode_rows(self.package, self.species)}
        self.assertEqual(rows["all"], dict(supergroup="all", genomes=3, agrees=1, differs=1, undecided=1, unchecked=0,
                                           best_tables="12:1"))
        self.assertEqual((rows["Fungi"]["differs"], rows["Metazoa"]["undecided"]), (1, 1))
        rows = checks.rdna_rows(self.package, self.species)  # the yeast: 100-200 published, 120 by its reads
        self.assertEqual([(r["species"], r["accession"], r["copies_reads"], r["units_assembled"], r["within"])
                          for r in rows], [("Saccharomyces cerevisiae", "Y", "120.4", "2", "t")])

    def test_run_and_figures(self):
        out = os.path.join(self.tmp.name, "results")
        os.makedirs(out, exist_ok=True)
        summary = {"species": 3, "busco_loci": 3, "left_out": {}}
        checks.run(self.package, out, self.species, summary, log=lambda t: None)
        self.assertEqual(summary["quality"]["genome"]["checks"]["any flag"], 1)
        self.assertEqual(summary["gencode"]["differs"], 1)
        self.assertEqual(summary["rdna_published"], {"species": 1, "measured": 1, "within": 1})
        md = benchmarks.summary_md(summary, "Eukore_x")
        self.assertIn("| genomes NCBI or its reports flag; with copies another clade's | 33.3%; 33.3% of 3 |", md)
        self.assertIn("agree with NCBI's; differ; undecided | 1; 1; 1 of 3 genomes |", md)
        if importlib.util.find_spec("matplotlib") is None:
            return
        from eukore_bench import figures  # noqa: PLC0415
        for fn in (figures.quality, figures.gencode):
            self.assertTrue(os.path.isfile(fn(out, self.package, out, "png")), fn.__name__)


class TreeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = cls.tmp.name
        cls.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        cls.odb = fx.odb_world(tmp)
        cls.package = os.path.join(tmp, "Eukore_x")
        write(os.path.join(cls.package, "genomes.tsv"), "accession\tspecies_taxid\trepresentative\tsource\n"
              + "".join(f"{a}\t{t}\tt\tncbi\n" for a, t in SPECIES.items()) + "Y2\t4932\tf\tassembly\n")
        rows = []
        for a in SPECIES:
            for i, locus in enumerate(fx.LOCI):
                status = "Single" if not (a == "M" and i) and not (a == "P" and i == 2) else "Fragmented"
                rows.append(f"{a}\t{locus}\t{status}\t{'low' if (a, i) == ('D', 2) else 'high'}\n")
        write(os.path.join(cls.package, "busco", "copies.tsv"), "accession\tlocus\tstatus\tconfidence\n" + "".join(rows))
        for locus in fx.LOCI:
            write(os.path.join(cls.package, "busco", "proteins", f"{locus}.faa.gz"),
                  "".join(f">{a}\n{fx.protein(hash((a, locus)) % 1000)}\n" for a in SPECIES))
        cls.species = {str(t): dict(supergroup="Metazoa" if t in (6253, 9606) else "Fungi") for t in SPECIES.values()}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_the_genomes_and_loci(self):
        tips, loci = tree.choose(self.package, self.species, self.ncbi)
        self.assertEqual(sorted(tips), ["A", "C", "D", "H", "K", "P", "Y"])  # M: one locus of three
        self.assertEqual(loci, sorted(fx.LOCI))
        self.assertEqual(tips["Y"][:2], (4932, dict(genus=4930, family=4893, order=4892, phylum=4890,
                                                    **{"class": 4891})))
        self.assertEqual((len(tips["P"][2]), len(tips["D"][2])), (2, 2))  # a fragment; a copy of confidence low
        tips, loci = tree.choose(self.package, self.species, self.ncbi, per_genus=1, n_loci=2)
        self.assertNotIn("P", tips)  # the genus' genome of the most loci
        self.assertEqual(len(loci), 2)

    def test_the_proteins_on_their_profiles(self):
        tips, loci = tree.choose(self.package, self.species, self.ncbi)
        rows, kept, total = tree.alignment(self.package, tips, loci, self.odb, fx.hmmalign(self.tmp.name),
                                           log=lambda t: None)
        self.assertEqual((kept, total), (360, 360))
        self.assertEqual({len(r) for r in rows.values()}, {360})
        self.assertEqual(rows["P"][240:], "-" * 120)  # its third locus not Single: no residue of it
        self.assertEqual(rows["Y"][:120], fx.protein(hash(("Y", fx.LOCI[0])) % 1000))

    def test_newick(self):
        t = tree.Tree(WRONG)
        self.assertEqual([t.name[n] for n in t.tips], ["Y", "P", "C", "D", "A", "H", "K"])
        self.assertEqual(len(t.children[0]), 3)
        masks, bit = t.masks()
        self.assertEqual(masks[0], 127)
        self.assertIn(bit["Y"] | bit["P"], masks)
        self.assertEqual(tree.Tree("(a,b);").name[1:], ["a", "b"])

    def evaluate(self, newick):
        tips, _loci = tree.choose(self.package, self.species, self.ncbi)
        group = {a: self.species[str(t[0])]["supergroup"] for a, t in tips.items()}
        return tree.evaluate(tree.Tree(newick), tips, self.ncbi, group)

    def test_taxa_as_branches(self):
        taxa, strays = self.evaluate(RIGHT)
        self.assertEqual([(r["rank"], r["name"], r["genomes"], r["monophyletic"]) for r in taxa],
                         [("genus", "Saccharomyces", 2, "t"), ("family", "Saccharomycetaceae", 3, "t"),
                          ("order", "Saccharomycetales", 4, "t"), ("class", "Saccharomycetes", 4, "t"),
                          ("phylum", "Ascomycota", 4, "t")])
        self.assertEqual(strays, [])
        self.assertEqual(tree.rank_summary(taxa)["family"], {"taxa": 1, "monophyletic": 1, "share": 1.0})
        # the tree is unrooted: the animals, the other side of the fungi's edge, are a branch however it is written
        t = tree.Tree("((A:1,H:1):1,D:1,(C:1,(K:1,(Y:1,P:1):1):1):1);")
        masks, bit = t.masks()
        self.assertIn(masks[0] ^ (bit["A"] | bit["H"]), set(masks) | {masks[0] ^ m for m in masks})

    def test_a_genome_apart_from_its_taxon(self):
        taxa, strays = self.evaluate(WRONG)
        at = {r["rank"]: r for r in taxa}
        self.assertEqual((at["genus"]["monophyletic"], at["family"]["monophyletic"]), ("t", "f"))
        self.assertEqual((at["family"]["f1"], at["family"]["strays"], at["family"]["intruders"]), ("0.800", 1, 0))
        # the order: all but the animals are a branch of the unrooted tree, which holds it whole with the
        # basidiomycete in it
        self.assertEqual((at["order"]["f1"], at["order"]["strays"], at["order"]["intruders"]), ("0.889", 0, 1))
        self.assertEqual(strays, [dict(accession="K", species_taxid=28985, name="Kluyveromyces lactis",
                                       supergroup="Fungi", rank="family", taxon="Saccharomycetaceae", genomes=3,
                                       ranks="family", beside="Ascaris lumbricoides, Homo sapiens")])

    def test_run_and_figure(self):
        out = os.path.join(self.tmp.name, "results")
        os.makedirs(out, exist_ok=True)
        summary = {"species": 8, "busco_loci": 3, "left_out": {}}
        tree.run(self.package, out, self.ncbi, self.species, summary, None, log=lambda t: None)
        self.assertIn("no BUSCO dataset", summary["left_out"]["tree"])
        tools = os.path.join(self.tmp.name, "tools")
        os.makedirs(tools, exist_ok=True)
        tree.run(self.package, out, self.ncbi, self.species, summary, self.odb, fx.hmmalign(tools),
                 fx.tool(tools, "FastTree", FASTTREE % WRONG), log=lambda t: None)
        self.assertNotIn("tree", summary["left_out"])
        self.assertEqual((summary["tree"]["genomes"], summary["tree"]["loci"], summary["tree"]["strays"]), (7, 3, 1))
        self.assertEqual(summary["tree"]["ranks"]["family"]["share"], 0.0)
        self.assertEqual([r["accession"] for r in fx.read_tsv(os.path.join(out, "tree_strays.tsv"))], ["K"])
        self.assertEqual(len(fx.read_tsv(os.path.join(out, "tree_tips.tsv"))), 7)
        self.assertFalse([n for n in os.listdir(out) if n.startswith("tree.") and n != "tree.nwk"])  # no work left
        md = benchmarks.summary_md(summary, "Eukore_x")
        self.assertIn("| genomes placed apart from their taxon in that tree | 1 |", md)
        self.assertIn("genus 100.0% of 1; family 0.0% of 1", md)
        if importlib.util.find_spec("matplotlib") is None:
            return
        from eukore_bench import figures  # noqa: PLC0415
        self.assertTrue(os.path.isfile(figures.tree(out, self.package, out, "png")))


if __name__ == "__main__":
    unittest.main()
