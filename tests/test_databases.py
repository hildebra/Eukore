#!/usr/bin/env python3
"""Tests of the package beside the other databases (benchmarks/eukore_bench/databases.py; since 2026-11-12): the
species by the databases that hold them, the species a genome is the first sequence of a marker for, the same
species' marker in two databases (the build's own checks, and another 18S reference's records aligned to Eukore's),
the comparison table.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import importlib.util
import json
import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore import benchmarks, ncbi as ncbi_mod  # noqa: E402
from eukore_bench import databases  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST, PARADOXUS = "GCA_000146045.2", "GCA_000000013.1"

MINIMAP2 = '''#!/usr/bin/env python3
"""minimap2's stand-in: each query (the last argument's FASTA) against each target (the one before) of its length,
their bases compared place by place; an alignment where half or more are alike."""
import sys


def fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


targets, queries = fasta(sys.argv[-2]), fasta(sys.argv[-1])
for q, s in queries.items():
    for t, u in targets.items():
        if len(s) == len(u):
            same = sum(a == b for a, b in zip(s, u))
            if same >= len(s) / 2:
                print(q, len(s), 0, len(s), "+", t, len(u), 0, len(u), same, len(s), 60, sep="\\t")
'''


def species_row(representative="-", **matches):
    row = dict(representative=representative, supergroup="Fungi", loci="3", rrna_records="1", organelle_genes="0")
    row.update({f"{ref.lower()}_match": "-" for ref in benchmarks.REFERENCES})
    row.update(matches)
    return row


class DatabasesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = cls.tmp.name
        cls.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        rng = random.Random(1)
        cls.ssu = "".join(rng.choice("ACGT") for _ in range(600))                   # the yeast's 18S
        cls.other = "".join(rng.choice("ACGT".replace(c, "")) if i % 10 == 0 else c  # S. paradoxus': 90% alike
                            for i, c in enumerate(cls.ssu))
        # the yeast named by PR2 and UNITE; S. paradoxus by no reference (EUKARYOME names its genus); K. lactis known
        # by published reads alone
        cls.species = {"4932": species_row(YEAST, pr2_match="species", unite_match="synonym"),
                       "27291": species_row(PARADOXUS, eukaryome_match="genus"),
                       "28985": species_row()}
        cls.package = os.path.join(tmp, "Eukore_x")
        os.makedirs(os.path.join(cls.package, "rrna"))
        os.makedirs(os.path.join(cls.package, "organelle"))
        columns = ("record", "kind", "species_taxid", "accession", "source", "regions", "length", "type",
                   "own_species", "confidence")
        records = [("y1", "rDNA", "4932", YEAST, "genome", "18S:0-600,ITS1:600-700,5.8S:700-800,ITS2:800-900,"
                    "28S:900-1500", "1500", "1", "t", "high"),
                   ("p1", "18S", "27291", PARADOXUS, "genome", "-", "600", "1", "t", "high"),
                   ("p2", "18S", "27291", PARADOXUS, "genome", "-", "600", "1", "t", "high"),   # its type again
                   ("p3", "18S", "27291", PARADOXUS, "genome", "-", "600", "2", "f", "low"),    # a contaminant's
                   ("e1", "18S", "4932", "-", "EUKARYOME", "-", "600", "-", "t", "high"),
                   ("e2", "18S", "28985", "-", "EUKARYOME", "-", "600", "-", "f", "low"),
                   ("r1", "rDNA", "28985", "-", "ROD", "18S:0-600,28S:900-1500", "1500", "-", "-", "unverified")]
        with open(os.path.join(cls.package, "rrna", "records.tsv"), "w") as fh:
            fh.write("\t".join(columns) + "\n" + "".join("\t".join(r) + "\n" for r in records))
        with gzip.open(os.path.join(cls.package, "rrna", "records.fna.gz"), "wt") as fh:
            fh.write(f">y1\n{cls.ssu}{'A' * 900}\n>p1\n{cls.other}\n>p2\n{cls.other}\n>p3\n{'C' * 600}\n")
        with open(os.path.join(cls.package, "organelle", "genes.tsv"), "w") as fh:
            fh.write("record\tgene\torganelle\tspecies_taxid\tsource\taccession\tlength\n")
        with open(os.path.join(cls.package, "manifest.json"), "w") as fh:
            json.dump({"version": "x", "eukore": "0.2.4", "inputs": {
                "odb": {"lineage": "eukaryota_odb12", "date": "2026-05-22"}, "sources": {"ROD": "v1.2"},
                "references": {"PR2": {"version": "5.1.1"}, "MIDORI2": {"version": "GB273", "screen": {
                    "records": 120, "screened": 100, "flagged": 7, "misnamed": 3, "chimera": 1, "other_gene": 1,
                    "weak": 2}}}}}, fh)
        cls.out = os.path.join(tmp, "results")
        os.makedirs(cls.out)
        with open(os.path.join(cls.out, "database_species.tsv"), "w") as fh:
            fh.write("database\tspecies_taxid\nPR2\t4932\nPR2\t28985\nPR2\t5476\nROD\t28985\n"
                     f"{benchmarks.ODB_NAME}\t4932\n")
        cls.theirs = {"Other": {4932, 6253}}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_species_by_the_databases_that_hold_them(self):
        sets = databases.species_sets(self.out, self.species, self.theirs, self.ncbi)
        self.assertEqual(sets[databases.EUKORE], {4932, 27291, 28985})
        self.assertEqual((sets["PR2"], sets["Other"]), ({4932, 28985, 5476}, {4932, 6253}))
        rows = {(r["supergroup"], r["databases"]): r for r in databases.upset_rows(sets, self.species, self.ncbi)}
        everywhere = f"Eukore + {benchmarks.ODB_NAME} + Other + PR2"
        self.assertEqual({d: r["species"] for (g, d), r in rows.items() if g == "all"},
                         {everywhere: 1, "Eukore": 1, "Eukore + PR2 + ROD": 1, "PR2": 1, "Other": 1})
        self.assertEqual((rows[("all", "Eukore")]["eukore"], rows[("all", "PR2")]["eukore"],
                          rows[("all", everywhere)]["n"]), ("t", "f", 4))
        self.assertIn(("Metazoa", "Other"), rows)  # Ascaris, which Eukore does not hold: by NCBI's lineage
        self.assertEqual(rows[("Fungi", "PR2")]["species"], 1)

    def test_species_a_genome_is_the_first_sequence_for(self):
        sets = databases.species_sets(self.out, self.species, self.theirs, self.ncbi)
        rows = {(r["supergroup"], r["marker"]): r for r in databases.novel_rows(self.package, self.species, sets)}
        # the 18S of two genomes: the yeast's PR2 names, S. paradoxus' no reference (EUKARYOME names its genus)
        self.assertEqual(rows[("all", "18S")], dict(supergroup="all", marker="18S", references="PR2, EUKARYOME",
                                                    genome=2, named=1, new=1))
        self.assertEqual([rows[("all", "ITS")][k] for k in ("references", "genome", "new")],
                         ["UNITE, EUKARYOME", 1, 0])
        # the yeast's operon: ROD (its species of database_species.tsv) holds K. lactis' alone
        self.assertEqual([rows[("all", "operon")][k] for k in ("references", "genome", "new")],
                         ["ROD, EUKARYOME", 1, 1])
        # no MIDORI2 read: nothing said of the mitochondrial genes
        self.assertEqual([rows[("Fungi", "mitochondrial genes")][k] for k in ("references", "named", "new")],
                         ["-", "-", "-"])

    def test_the_builds_own_checks(self):
        rows = {r["database"]: r for r in databases.published_agreement(self.package)}
        self.assertEqual([rows["EUKARYOME"][k] for k in databases.AGREEMENT_COLUMNS[2:10]], [2, 1, 1, 0, 0, 2, 1, 1])
        self.assertEqual([rows["ROD"][k] for k in ("records", "undecided")], [1, 1])
        self.assertEqual([rows["MIDORI2"][k] for k in ("records", "agree", "disagree", "undecided")], [100, 93, 5, 2])

    def test_a_records_verdict(self):
        taxa = {"y": 4932, "p": 27291}
        self.assertEqual(databases.verdict_of([(0.99, "y"), (1.0, "p")], taxa, 4932), ("agree", 0.99, (1.0, 27291)))
        self.assertEqual(databases.verdict_of([(0.90, "y"), (0.99, "p")], taxa, 4932)[0], "disagree")
        self.assertEqual(databases.verdict_of([(0.99, "p")], taxa, 4932), ("disagree", None, (0.99, 27291)))
        self.assertEqual(databases.verdict_of([(0.90, "y"), (0.95, "p")], taxa, 4932)[0], "unlike")  # no one's gene
        self.assertEqual(databases.verdict_of([(0.96, "y"), (0.975, "p")], taxa, 4932)[0], "unlike")  # no margin
        self.assertEqual(databases.verdict_of([(0.80, "p")], taxa, 4932)[0], "unlike")  # aligned to others alone
        self.assertEqual(databases.verdict_of([], taxa, 4932), ("undecided", None, None))
        self.assertEqual(databases.shared_rank(self.ncbi, 4932, 27291), "genus")
        self.assertEqual(databases.shared_rank(self.ncbi, 4932, 6253), "none")

    def reference(self):
        """Another 18S reference: the yeast's 18S under its name, S. paradoxus' under the yeast's name (a strain's
        taxid), one 95% like the yeast's under its name, a sequence like nothing under S. paradoxus', and records of
        a species without a genome's 18S and of a genus."""
        fasta, table = os.path.join(self.tmp.name, "PR2.fna.gz"), os.path.join(self.tmp.name, "PR2_taxids.tsv")
        far = "".join("ACGT"["ACGT".index(c) - 1] if i % 20 == 5 else c for i, c in enumerate(self.ssu))
        with gzip.open(fasta, "wt") as fh:
            fh.write(f">a\n{self.ssu}\n>b\n{self.other}\n>u\n{far}\n>c\n{'G' * 600}\n>d\n{self.ssu}\n"
                     f">g\n{self.ssu}\n>short\n{self.ssu[:100]}\n")
        with open(table, "w") as fh:
            fh.write("a\t4932\nb\t559292\nu\t4932\nc\t27291\nd\t28985\ng\t4930\nshort\t4932\n")
        return fasta, table

    def test_another_reference_aligned(self):
        fasta, table = self.reference()
        work = os.path.join(self.tmp.name, "work")
        os.makedirs(work, exist_ok=True)
        taxa = databases.eukore_ssu(self.package, self.species, os.path.join(work, "ref.fna"))
        self.assertEqual(taxa, {"y1": 4932, "p1": 27291})  # one of a species' type, none of confidence low
        with open(os.path.join(work, "ref.fna")) as fh:
            self.assertIn(f">y1\n{self.ssu}\n", fh.read())  # the unit's 18S alone
        row, records = databases.ssu_agreement(self.package, self.species, "PR2", fasta, table, self.ncbi,
                                               fx.tool(self.tmp.name, "minimap2", MINIMAP2), 1, work)
        # agree a; disagree b; unlike u; c aligned to nothing; the yeast agrees by one record, S. paradoxus is
        # undecided: no species differs
        self.assertEqual([row[k] for k in databases.AGREEMENT_COLUMNS[2:10]], [4, 1, 1, 1, 1, 2, 1, 0])
        self.assertEqual(records[0], dict(database="PR2", record="b", verdict="disagree", species_taxid=4932,
                                          name="Saccharomyces cerevisiae", own_identity="90.0", other_taxid=27291,
                                          other_name="Saccharomyces paradoxus", other_identity="100.0",
                                          shared_rank="genus"))
        self.assertEqual([(r["record"], r["verdict"], r["own_identity"], r["other_identity"]) for r in records[1:]],
                         [("u", "unlike", "95.0", "85.0")])
        self.assertEqual(databases.differing({1: {"unlike"}, 2: {"agree", "disagree"}, 3: {"undecided"}}), 1)

    def test_run_table_and_figures(self):
        fasta, table = self.reference()
        out = os.path.join(self.tmp.name, "run")
        os.makedirs(out)
        with open(os.path.join(self.out, "database_species.tsv")) as src, \
                open(os.path.join(out, "database_species.tsv"), "w") as fh:
            fh.write(src.read())
        summary = {"species": 2, "busco_loci": 3, "left_out": {}}
        databases.run(self.package, None, out, self.ncbi, self.species, summary, {"PR2": {4932, 27291}},
                      [("PR2", fasta, table)], fx.tool(self.tmp.name, "minimap2", MINIMAP2), log=lambda t: None)
        rows = {r["database"]: r for r in fx.read_tsv(os.path.join(out, "databases.tsv"))}
        self.assertEqual(list(rows)[0], "Eukore")
        # Eukore's records: its rRNA records' (the table's lengths); PR2's: the reference's records of its table
        self.assertEqual((rows["Eukore"]["version"], rows["Eukore"]["species"], rows["Eukore"]["records"],
                          rows["Eukore"]["bases"]), ("x (Eukore 0.2.4)", "3", "7", "6000"))
        self.assertEqual((rows["PR2"]["version"], rows["PR2"]["species"], rows["PR2"]["both"],
                          rows["PR2"]["eukore_only"], rows["PR2"]["records"], rows["PR2"]["bases"]),
                         ("5.1.1", "2", "2", "1", "7", "3700"))  # compare.py's species of it, not the build's
        self.assertTrue(rows["PR2"]["markers"].startswith("18S"))
        self.assertEqual((rows["ROD"]["version"], rows[benchmarks.ODB_NAME]["version"]),
                         ("v1.2", "eukaryota_odb12 2026-05-22"))
        self.assertEqual(rows["SILVA"]["species"], "-")  # listed, not read
        self.assertEqual(summary["upset"]["eukore_alone"], 0)
        self.assertEqual(summary["novel"]["18S"]["new"], 1)
        self.assertEqual(summary["agreement"]["PR2 (18S)"]["disagree"], 1)
        self.assertEqual(len(fx.read_tsv(os.path.join(out, "agreement_records.tsv"))), 2)
        self.assertFalse([n for n in os.listdir(out) if n.startswith("agreement.") and n != "agreement.tsv"])
        md = benchmarks.summary_md(summary, "Eukore_x")
        self.assertIn("| species in Eukore and in none of 3 other databases | 0 |", md)
        self.assertIn("| species with a genome's 18S named in none of PR2, EUKARYOME | 1 of 2 |", md)
        self.assertIn("| PR2 (18S): records agreeing with Eukore's genomes; another species'; unlike | 1; 1; 1 of 4 |", md)
        if importlib.util.find_spec("matplotlib") is None:
            return
        from eukore_bench import figures  # noqa: PLC0415
        figs = os.path.join(self.tmp.name, "figures")
        os.makedirs(figs, exist_ok=True)
        for fn in (figures.upset, figures.novel, figures.agreement):
            self.assertTrue(os.path.isfile(fn(out, self.package, figs, "png")), fn.__name__)
        self.assertIn("| PR2 | 5.1.1 | 2 | 2 | 1 | 0 | 7 | 3,700 |", "\n".join(figures.databases_md(out)))

    def test_the_tables_words(self):
        words = databases.read_words()
        self.assertEqual(list(words)[0], "Eukore")
        self.assertTrue({"EukDetect2", "MetaPhlAn4", "PR2", "SILVA", "UNITE", "EUKARYOME", "KSGP", "MIDORI2", "ROD",
                         benchmarks.ODB_NAME} <= set(words))
        self.assertTrue(all(set(w) == set(databases.WORDS) and all(w.values()) for w in words.values()))


if __name__ == "__main__":
    unittest.main()
