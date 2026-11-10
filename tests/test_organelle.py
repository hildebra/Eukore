#!/usr/bin/env python3
"""Tests of Eukore's organelle genes (eukore/organelle.py): RefSeq's GenBank records read (a COX1 across a circular
genome's origin on the minus strand, a qualifier over two lines, an animal's 16S as its LSU, a plastid's rRNA twice,
its 4.5S left out), a species' record chosen, an assembly's mitochondrion searched (stand-in barrnap and nhmmer: its
rRNA by bacterial models, its COX1 of two exons), a mitochondrion's CYTB and a plastid's rbcL, matK and tufA (since
2026-11-07; an assembly searched before them searched again, one searched for them all not), and the genes written into
the package; since 2026-11-10 (AuditTest) a gene's best feature, /pseudo, /codon_start, /transl_table, the spacer
across the origin and the caches' keys.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import gencode, genomes, ncbi as ncbi_mod, organelle, package, rrna  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST = "GCA_000146045.2"


def feature(key, location, *quals):
    lines = [f"     {key:<16}{location}\n"]
    for q in quals:
        for i, part in enumerate(q.split("\n")):
            lines.append(" " * 21 + part + "\n")
    return "".join(lines)


def gbff_record(accession, seq, taxid, organelle_name, features):
    origin = "".join(f"{i + 1:>9} " + " ".join(seq[j:j + 10].lower() for j in range(i, min(i + 60, len(seq)), 10))
                     + "\n" for i in range(0, len(seq), 60))
    return (f"LOCUS       {accession.split('.')[0]}  {len(seq)} bp    DNA     circular MAM 01-JAN-2026\n"
            f"DEFINITION  test {organelle_name}, complete genome.\n"
            f"ACCESSION   {accession.split('.')[0]}\nVERSION     {accession}\n"
            f"SOURCE      test\n  ORGANISM  test\n            Eukaryota.\n"
            "FEATURES             Location/Qualifiers\n"
            + feature("source", f"1..{len(seq)}", '/organism="test"', f'/organelle="{organelle_name}"',
                      f'/db_xref="taxon:{taxid}"')
            + "".join(features) + "ORIGIN      \n" + origin + "//\n")


def random_seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


HUMAN_MT = random_seq(3000, 1)
HUMAN_MT2 = random_seq(2500, 2)
PLASTID = random_seq(7000, 3)  # rbcL and matK beyond its rRNA
YEAST_MT = random_seq(4000, 4)


def release(src_dir):
    folder = organelle.folder_of(src_dir)
    os.makedirs(folder)
    with gzip.open(os.path.join(folder, "mitochondrion.1.genomic.gbff.gz"), "wt") as fh:
        fh.write(gbff_record("NC_000001.1", HUMAN_MT, 9606, "mitochondrion", [
            feature("rRNA", "101..500", '/product="12S ribosomal RNA"'),
            feature("rRNA", "601..1500", '/gene="RNR2"', '/product="16S ribosomal RNA"'),
            feature("CDS", "complement(join(2901..3000,1..200))", '/gene="COX1"',
                    '/product="cytochrome c oxidase\nsubunit I"'),
            feature("CDS", "1601..2000", '/product="cytochrome c oxidase subunit II"'),
            feature("CDS", "2101..2800", '/gene="CYTB"', '/product="cytochrome b"'),
            feature("rRNA", "J01415.2:1..10", '/product="12S ribosomal RNA"')]))
        fh.write(gbff_record("NC_000002.1", HUMAN_MT2, 9606, "mitochondrion", [
            feature("rRNA", "101..500", '/product="12S ribosomal RNA"')]))
    with gzip.open(os.path.join(folder, "plastid.1.genomic.gbff.gz"), "wt") as fh:
        fh.write(gbff_record("NC_000003.1", PLASTID, 28985, "plastid:chloroplast", [
            feature("rRNA", "101..1600", '/gene="rrn16"'),
            feature("rRNA", "complement(2001..2400)", '/product="16S ribosomal RNA"'),
            feature("rRNA", "2501..2600", '/product="4.5S ribosomal RNA"'),
            feature("rRNA", "2701..2820", '/product="5S ribosomal RNA"'),
            feature("CDS", "3001..4428", '/gene="rbcL"',
                    '/product="ribulose-1,5-bisphosphate carboxylase/oxygenase\nlarge subunit"'),
            feature("CDS", "complement(4501..6000)", '/product="maturase K"'),
            feature("CDS", "6101..6400", '/gene="rbcS"',
                    '/product="ribulose-1,5-bisphosphate carboxylase/oxygenase small subunit"'),
            feature("CDS", "6401..6900", '/product="translation elongation factor Tu"'),
            feature("CDS", "6901..6990", '/product="cytochrome b6"')]))
    with open(os.path.join(folder, "release.txt"), "w") as fh:
        fh.write("mitochondrion.1.genomic.gbff.gz\nplastid.1.genomic.gbff.gz\n")


BARRNAP = r"""#!/usr/bin/env python3
import sys
kingdom = sys.argv[sys.argv.index("--kingdom") + 1]
names = [l[1:].split()[0] for l in open(sys.argv[-1]) if l.startswith(">")]
models = {"mito": ("12S_rRNA", "16S_rRNA"), "bac": ("16S_rRNA", "23S_rRNA")}[kingdom]
for n in names:
    print(f"{n}\tbarrnap\trRNA\t11\t60\t1e-9\t+\t.\tName={models[0]};product=x")
    print(f"{n}\tbarrnap\trRNA\t101\t300\t1e-9\t-\t.\tName={models[1]};product=x (partial)")
    print(f"{n}\tbarrnap\trRNA\t3101\t3500\t1e-9\t+\t.\tName={models[1]};product=x")
"""

NHMMER = r"""#!/usr/bin/env python3
import sys
args = sys.argv[1:]
tbl = args[args.index("--tblout") + 1]
q, t = args[-2], args[-1]
seqs, name = {}, None
for l in open(q):
    if l.startswith(">"):
        name = l[1:].split()[0]; seqs[name] = ""
    else:
        seqs[name] += l.strip()
target = [l[1:].split()[0] for l in open(t) if l.startswith(">")][0]
first = next((n for n in seqs if n.startswith("COX1|")), None)
with open(tbl, "w") as fh:
    fh.write("# target query\n")
    if first is None:
        sys.exit(0)
    half = len(seqs[first]) // 2
    fh.write(f"{target} - {first} - 3 {half} 1003 {1000 + half} 1003 {1000 + half} 4000 + 1e-50 500.0 0 -\n")
    fh.write(f"{target} - {first} - {half + 1} {2 * half} 1501 {1500 + half} 1501 {1500 + half} "
             "4000 + 1e-40 300.0 0 -\n")
"""


class ParseTest(unittest.TestCase):
    def test_genes_named(self):
        self.assertEqual(organelle.rrna_gene("16s ribosomal rna", True), "LSU")
        self.assertEqual(organelle.rrna_gene("16s ribosomal rna", False), "SSU")
        self.assertEqual(organelle.rrna_gene("rrn23"), "LSU")
        self.assertEqual(organelle.rrna_gene("rrn5 5s ribosomal rna"), "5S")
        self.assertEqual(organelle.rrna_gene("l-rrna 16s ribosomal rna", False), "LSU")
        self.assertIsNone(organelle.rrna_gene("4.5s ribosomal rna"))
        self.assertEqual(organelle.cds_gene([("product", "cytochrome c oxidase subunit I")], "mitochondrion"), "COX1")
        self.assertNotEqual(organelle.cds_gene([("product", "cytochrome c oxidase subunit II")], "mitochondrion"),
                            "COX1")
        self.assertEqual(organelle.cds_gene([("gene", "COI")], "mitochondrion"), "COX1")
        self.assertEqual(organelle.cds_gene([("gene", "rbcL")], "plastid"), "rbcL")
        self.assertEqual(organelle.cds_gene([("product", "RuBisCO large subunit")], "plastid"), "rbcL")
        self.assertEqual(organelle.cds_gene([("gene", "matK"), ("product", "maturase K")], "plastid"), "matK")
        self.assertIsNone(organelle.cds_gene([("gene", "rbcS")], "plastid"))
        self.assertIsNone(organelle.cds_gene([("gene", "rbcL")], "mitochondrion"))
        self.assertIsNone(organelle.cds_gene([("gene", "COX1")], "plastid"))
        for quals in ([("gene", "cob"), ("product", "apocytochrome b")], [("product", "cytochrome b")],
                      [("gene", "MT-CYB")]):
            self.assertEqual(organelle.cds_gene(quals, "mitochondrion"), "CYTB", quals)
        for product in ("cytochrome b6", "cytochrome b559 alpha subunit", "cytochrome b-245", "cytochrome b5"):
            self.assertIsNone(organelle.cds_gene([("product", product)], "mitochondrion"), product)
        self.assertEqual(organelle.cds_gene([("gene", "tufA")], "plastid"), "tufA")
        self.assertEqual(organelle.cds_gene([("product", "elongation factor Tu")], "plastid"), "tufA")
        self.assertIsNone(organelle.cds_gene([("gene", "tufA")], "mitochondrion"))

    def test_refseq_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            release(os.path.join(tmp, "sources"))
            rows = organelle.refseq_genes(os.path.join(tmp, "sources"), os.path.join(tmp, "work"), ncbi, lambda _t: 0)
            got = {(r["seqid"], r["gene"]): r for r in rows}
            self.assertEqual(sorted(got), [("NC_000001.1", "COX1"), ("NC_000001.1", "CYTB"), ("NC_000001.1", "LSU"),
                                           ("NC_000001.1", "SSU"), ("NC_000002.1", "SSU"), ("NC_000003.1", "5S"),
                                           ("NC_000003.1", "SSU"), ("NC_000003.1", "matK"), ("NC_000003.1", "rbcL"),
                                           ("NC_000003.1", "tufA")])
            self.assertEqual(got[("NC_000001.1", "CYTB")]["sequence"], HUMAN_MT[2100:2800])
            self.assertEqual(got[("NC_000003.1", "tufA")]["sequence"], PLASTID[6400:6900])
            rbcl, matk = got[("NC_000003.1", "rbcL")], got[("NC_000003.1", "matK")]
            self.assertEqual(rbcl["sequence"], PLASTID[3000:4428])
            self.assertEqual((matk["strand"], matk["sequence"]), ("-", rrna.revcomp(PLASTID[4500:6000])))
            cox1 = got[("NC_000001.1", "COX1")]
            self.assertEqual((cox1["start"], cox1["end"], cox1["strand"], cox1["exons"]), ("2901", "200", "-", "2"))
            self.assertEqual(cox1["sequence"], rrna.revcomp(HUMAN_MT[2900:] + HUMAN_MT[:200]))
            self.assertEqual(cox1["spliced"], cox1["sequence"])
            self.assertEqual(got[("NC_000001.1", "LSU")]["sequence"], HUMAN_MT[600:1500])  # an animal's 16S
            ssu = got[("NC_000003.1", "SSU")]
            self.assertEqual((ssu["copies"], ssu["sequence"]), ("2", PLASTID[100:1600]))  # the first of two
            # parsed once: the cache is read again
            self.assertEqual(len(organelle.refseq_genes(os.path.join(tmp, "sources"), os.path.join(tmp, "work"),
                                                        ncbi)), len(rows))
            chosen = organelle.choose(rows, ncbi)
            self.assertEqual({r["seqid"] for r in chosen[(9606, "mitochondrion")]}, {"NC_000001.1"})  # most genes
            self.assertEqual(sorted(chosen), [(9606, "mitochondrion"), (28985, "plastid")])
            near = organelle.Nearest(chosen, ncbi)
            self.assertEqual([s for s, _q in near.queries(4932, "mitochondrion")],
                             ["SSU|NC_000001.1", "LSU|NC_000001.1", "COX1|NC_000001.1",
                              "CYTB|NC_000001.1"])  # the nearest with each gene
            self.assertEqual([s for s, _q in near.queries(952, "plastid")],
                             ["SSU|NC_000003.1", "5S|NC_000003.1", "rbcL|NC_000003.1", "matK|NC_000003.1",
                              "tufA|NC_000003.1"])


class JoinTest(unittest.TestCase):
    def test_hits_joined_in_the_querys_order_alone(self):
        # a 16 kb mitochondrion opened inside COX1: its first half at the molecule's end, its second at its start;
        # neither is joined to the other (until 2026-11-06 the gene was the molecule, 1-16,000)
        a = dict(start=15200, end=16000, q1=1, q2=800, strand="+")
        b = dict(start=1, end=800, q1=801, q2=1600, strand="+")
        self.assertEqual(organelle.collinear(a, [a, b]), [a])
        # one turn on (across_origin), the second half follows the first: the gene joined across the origin, its
        # start after its end, its sequence the molecule's end then its start
        chain = organelle.collinear(a, [a, b] + organelle.across_origin([a, b], 16000))
        self.assertEqual([(h["start"], h["q1"]) for h in chain], [(15200, 1), (16001, 801)])
        self.assertEqual(organelle.gene_span(chain, "+", 1600, 16000), (15200, 800))
        rng = random.Random(8)
        molecule = "".join(rng.choice("ACGT") for _ in range(16000))
        gene = molecule[15199:] + molecule[:800]
        self.assertEqual(organelle.subsequence(molecule, 15200, 800, "+"), gene)
        # on the minus strand: the gene's first half at the molecule's start, its second at its end
        m_a = dict(start=1, end=800, q1=1, q2=800, strand="-")
        m_b = dict(start=15200, end=16000, q1=801, q2=1600, strand="-")
        chain = organelle.collinear(m_a, [m_a, m_b] + organelle.across_origin([m_a, m_b], 16000))
        self.assertEqual(organelle.gene_span(chain, "-", 1600, 16000), (15200, 800))
        self.assertEqual(organelle.subsequence(molecule, 15200, 800, "-"), rrna.revcomp(gene))
        # a gene inside the molecule: as before, its ends clipped to the molecule
        inside = dict(start=100, end=900, q1=11, q2=800, strand="+")
        self.assertEqual(organelle.gene_span([inside], "+", 800, 16000), (90, 900))
        self.assertEqual(organelle.gene_span([dict(inside, start=5, end=805)], "+", 800, 16000), (1, 805))
        # a gene's exons in order on either strand, a concatemer's next copy (the query again) left out
        c = dict(start=17000, end=17800, q1=801, q2=1600, strand="+")
        again = dict(start=18000, end=18800, q1=1, q2=800, strand="+")
        self.assertEqual(organelle.collinear(a, [a, c, again]), [a, c])
        m1, m2 = dict(start=5000, end=5800, q1=1, q2=800, strand="-"), dict(start=3000, end=3800, q1=801, q2=1600,
                                                                          strand="-")
        self.assertEqual(organelle.collinear(m2, [m1, m2]), [m2, m1])

    def test_a_mitogenome_of_several_records(self):
        # lice's minicircles: rrnS on one record, rrnL and COX1 on another; the species takes them all
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        rows = [dict(seqid="NC_1.1", organelle="mitochondrion", taxid="6253", record_length="3000", gene="SSU"),
                dict(seqid="NC_2.1", organelle="mitochondrion", taxid="6253", record_length="4000", gene="LSU"),
                dict(seqid="NC_2.1", organelle="mitochondrion", taxid="6253", record_length="4000", gene="COX1"),
                dict(seqid="NC_3.1", organelle="mitochondrion", taxid="6253", record_length="5000", gene="COX1")]
        got = organelle.choose(rows, ncbi)[(6253, "mitochondrion")]
        self.assertEqual(sorted((r["gene"], r["seqid"]) for r in got),
                         [("COX1", "NC_2.1"), ("LSU", "NC_2.1"), ("SSU", "NC_1.1")])
        # an animal or a fungus holds no plastid: not searched for one
        self.assertEqual(organelle.organelles_of(ncbi, 6253), ("mitochondrion",))
        self.assertEqual(organelle.organelles_of(ncbi, 952), organelle.KINDS)


class AssemblyTest(unittest.TestCase):
    def test_an_assembly_into_the_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            src = os.path.join(tmp, "sources")
            release(src)
            work = os.path.join(tmp, "work", "organelle")
            chosen = organelle.choose(organelle.refseq_genes(src, work, ncbi, lambda _t: 0), ncbi)
            genomes_dir = os.path.join(tmp, "genomes")
            os.makedirs(genomes_dir)
            with gzip.open(genomes.genome_file(genomes_dir, YEAST), "wt") as fh:
                fh.write(f">c1 Saccharomyces cerevisiae chromosome I\n{random_seq(5000, 5)}\n"
                         f">m1 Saccharomyces cerevisiae mitochondrion, complete genome\n{YEAST_MT}\n")
            jobs = {YEAST: dict(accession=YEAST, species_taxid=4932, taxid=559292),
                    "EUK_6253": dict(accession="EUK_6253", species_taxid=6253, path="x.fna"),  # records: not searched
                    "GCA_000000099.1": dict(accession="GCA_000000099.1", species_taxid=9606)}  # not on disk
            barrnap, nhmmer = fx.tool(tmp, "barrnap", BARRNAP), fx.tool(tmp, "nhmmer", NHMMER)
            ran = organelle.assembly_pass(work, jobs, lambda a: genomes.genome_file(genomes_dir, a), chosen, ncbi,
                                          barrnap, nhmmer, 1, 1, lambda _t: 0)
            self.assertEqual(ran, (1, 0))
            self.assertEqual(organelle.assembly_pass(work, jobs, lambda a: genomes.genome_file(genomes_dir, a),
                                                     chosen, ncbi, barrnap, nhmmer, 1, 1, lambda _t: 0), (0, 0))
            records = organelle.collect(ncbi, jobs, {4932, 9606, 6253}, chosen, work)
            rows = {(r["species_taxid"], r["gene"]): r for _n, _s, r in records}
            self.assertEqual(sorted(rows), [(4932, "COX1"), (4932, "LSU"), (4932, "SSU"),
                                            (9606, "COX1"), (9606, "CYTB"), (9606, "LSU"), (9606, "SSU")])
            lsu = rows[(4932, "LSU")]  # bacterial models (a yeast's mitochondrion), the whole hit of two
            self.assertEqual((lsu["start"], lsu["end"], lsu["copies"], lsu["partial"], lsu["source"]),
                             ("3101", "3500", "2", "f", "assembly"))
            cox1 = rows[(4932, "COX1")]  # two hits joined, the query's first two bases reached
            seqs = {r["record"]: s for _n, s, r in records}
            self.assertEqual((cox1["start"], cox1["end"], cox1["exons"], cox1["partial"]), ("1001", "1650", "2", "f"))
            self.assertEqual(seqs[cox1["record"]], YEAST_MT[1000:1650])
            self.assertEqual((rows[(9606, "SSU")]["source"], rows[(9606, "SSU")]["accession"]), ("refseq", "-"))
            # into the package
            collected = dict(genomes={YEAST: dict(accession=YEAST, species_taxid=4932, taxid=559292, organism="y",
                                                  source="ncbi", representative="t", loci=0)},
                             copies={YEAST: {}}, features=[], copy_rows=[], taxids=set())
            out = package.write(os.path.join(tmp, "out", "Eukore_2026.11"), "2026.11", ncbi, collected, [],
                                None, None, {}, records)
            self.assertIsNone(package.check(out, deep=True))
            table = package.read_table(out, "organelle/genes.tsv")
            self.assertEqual((len(table), set(table[0])), (7, set(organelle.GENE_COLUMNS)))
            self.assertEqual(package.read_table(out, "species.tsv")[0]["organelle_genes"], "3")
            with open(os.path.join(out, "manifest.json")) as fh:
                counts = json.load(fh)["counts"]
            self.assertEqual(counts["organelle_genes"]["mitochondrion"], {"COX1": 2, "CYTB": 1, "LSU": 2, "SSU": 2})
            with gzip.open(os.path.join(out, "organelle", "genes.fna.gz"), "rt") as fh:
                self.assertEqual(fh.read().count(">"), 7)

    def test_an_assembly_searched_before_its_genes_were_added(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            src = os.path.join(tmp, "sources")
            release(src)
            work = os.path.join(tmp, "work", "organelle")
            chosen = organelle.choose(organelle.refseq_genes(src, work, ncbi, lambda _t: 0), ncbi)
            genomes_dir = os.path.join(tmp, "genomes")
            os.makedirs(genomes_dir)
            plant, yeast, worm = "GCA_000000952.1", YEAST, "GCA_000006253.1"
            for acc in (plant, yeast, worm):
                with gzip.open(genomes.genome_file(genomes_dir, acc), "wt") as fh:
                    fh.write(f">m1 mitochondrion\n{YEAST_MT}\n>p1 chloroplast\n{PLASTID}\n")
            jobs = {plant: dict(accession=plant, species_taxid=952), yeast: dict(accession=yeast, species_taxid=4932),
                    worm: dict(accession=worm, species_taxid=6253)}
            head = "\t".join(organelle.ASSEMBLY_COLUMNS) + "\n"
            for acc, searched in ((plant, ""),  # a file of before 2026-11-07: no #searched line
                                  (yeast, "#searched\tmitochondrion:SSU,LSU,5S,COX1\n"),  # before CYTB
                                  # every gene of its mitochondrion, but not what it was searched with (a file of
                                  # before 2026-11-10: no `inputs`): searched again, once
                                  (worm, "#searched\tmitochondrion:" + ",".join(organelle.GENES_OF["mitochondrion"])
                                   + "\n")):
                os.makedirs(os.path.dirname(organelle.assembly_file(work, acc)), exist_ok=True)
                with open(organelle.assembly_file(work, acc), "w") as fh:
                    fh.write(searched + head)
            self.assertEqual(organelle.searched_of(organelle.assembly_file(work, plant))["plastid"],
                             set(organelle.OLD_GENES))
            barrnap, nhmmer = fx.tool(tmp, "barrnap", BARRNAP), fx.tool(tmp, "nhmmer", NHMMER)

            def run():
                return organelle.assembly_pass(work, jobs, lambda a: genomes.genome_file(genomes_dir, a), chosen,
                                               ncbi, barrnap, nhmmer, 1, 1, lambda _t: 0)

            def searched(acc):  # its genes and models; its inputs' organelles (their keys: test_caches_keyed)
                got = organelle.searched_of(organelle.assembly_file(work, acc))
                return {k: v for k, v in got.items() if k != "inputs"}, sorted(i.split("=")[0] for i in got["inputs"])
            self.assertEqual(run(), (3, 0))
            self.assertEqual(searched(plant), (dict({k: set(organelle.GENES_OF[k]) for k in organelle.KINDS},
                                                    models={"nhmmer"}), ["mitochondrion", "plastid"]))
            self.assertEqual(searched(yeast), ({"mitochondrion": set(organelle.GENES_OF["mitochondrion"]),
                                                "models": {"nhmmer"}}, ["mitochondrion"]))
            self.assertEqual(run(), (0, 0))  # searched for every gene, with these inputs: not again
            records = organelle.collect(ncbi, jobs, {952}, chosen, work)
            self.assertIn(("mitochondrion", "COX1"), {(r["organelle"], r["gene"]) for _n, _s, r in records})


class Format5Test(unittest.TestCase):
    """Since 2026-11-07 (format 5): GenBank locations read whole (a trans-spliced gene), more genes (COX3, nad1, nad5,
    atp6, psbA, ycf1, the trnH-psbA spacer), the genes' common names, coding sequences, proteins and checks."""

    def test_locations(self):
        self.assertEqual(organelle.parts_of("complement(join(2901..3000,1..200))"),
                         ([(1, 200, "-"), (2901, 3000, "-")], False))
        self.assertEqual(organelle.parts_of("join(complement(5000..5100),<200..300)"),
                         ([(5000, 5100, "-"), (200, 300, "+")], True))
        self.assertEqual(organelle.parts_of("join(1..10,20..30)"), ([(1, 10, "+"), (20, 30, "+")], False))
        self.assertIsNone(organelle.parts_of("join(J01415.2:1..10,20..30)"))
        self.assertTrue(organelle.trans_spliced([(5000, 5100, "-"), (200, 300, "+")]))
        self.assertFalse(organelle.trans_spliced([(1, 200, "-"), (2901, 3000, "-")]))  # across the origin
        self.assertTrue(organelle.trans_spliced([(500, 600, "+"), (100, 200, "+"), (900, 950, "+"), (10, 20, "+")]))
        self.assertEqual(organelle.parts_read(organelle.parts_text([(1, 5, "+"), (9, 12, "-")])),
                         [(1, 5, "+"), (9, 12, "-")])

    def test_genes_named(self):
        for quals, kind, gene in (([("gene", "nad1")], "mitochondrion", "nad1"),
                                  ([("product", "NADH dehydrogenase subunit 5")], "mitochondrion", "nad5"),
                                  ([("product", "ATP synthase F0 subunit 6")], "mitochondrion", "atp6"),
                                  ([("gene", "COX3")], "mitochondrion", "COX3"),
                                  ([("product", "cytochrome c oxidase subunit III")], "mitochondrion", "COX3"),
                                  ([("product", "photosystem II protein D1")], "plastid", "psbA"),
                                  ([("gene", "ycf1")], "plastid", "ycf1")):
            self.assertEqual(organelle.cds_gene([(n, v) for n, v in quals], kind), gene, quals)
        for quals, kind in (([("product", "NADH dehydrogenase subunit 4L")], "mitochondrion"),
                            ([("product", "NADH dehydrogenase subunit 1")], "plastid"),
                            ([("product", "cytochrome c oxidase subunit II")], "mitochondrion")):
            self.assertIsNone(organelle.cds_gene(quals, kind), quals)
        self.assertTrue(organelle.is_his_trna([("gene", "trnH-GUG"), ("product", "tRNA-His")]))
        self.assertFalse(organelle.is_his_trna([("product", "tRNA-Lys")]))
        self.assertEqual([organelle.name_of("mitochondrion", "SSU", animal=True), organelle.name_of("plastid", "LSU"),
                          organelle.name_of("mitochondrion", "LSU", plant=True), organelle.name_of("mitochondrion", "SSU"),
                          organelle.name_of("mitochondrion", "COX1", animal=True)], ["12S", "23S", "26S", "SSU", "COX1"])

    def test_spacer_and_a_trans_spliced_gene(self):
        seq = random_seq(4000, 9)
        rec = dict(seqid="NC_9.1", taxid=952, length=len(seq), circular=True, sequence=seq, features=[
            ("tRNA", ["complement(1..75)"], [("gene", "trnH-GUG"), ("product", "tRNA-His")]),
            ("CDS", ["complement(400..1459)"], [("gene", "psbA")])])
        got = organelle.record_genes(rec, "plastid")
        spacer = got["trnH-psbA"]
        self.assertEqual((spacer["start"], spacer["end"], spacer["strand"], spacer["regions"]),
                         (1, 599, "-", "psbA3':0-200,spacer:200-524,trnH:524-599"))
        self.assertEqual(spacer["sequence"], rrna.revcomp(seq[:599]))
        mito = dict(rec, features=[("CDS", ["join(complement(3001..3100),100..200)"], [("gene", "nad1")])])
        nad1 = organelle.record_genes(mito, "mitochondrion")["nad1"]
        self.assertEqual((nad1["trans_spliced"], nad1["exons"], nad1["sequence"]),
                         (True, 2, rrna.revcomp(seq[3000:3100]) + seq[99:200]))

    def test_models_combined(self):
        # nhmmer's single hit kept; its several hits (introns) replaced by miniprot's whole model; a gene nhmmer missed
        # taken; miniprot's partial model never replaces nhmmer's
        found = {"CYTB": dict(gene="CYTB", exons=1, model="nhmmer"), "COX1": dict(gene="COX1", exons=3, model="nhmmer"),
                 "nad5": dict(gene="nad5", exons=2, model="nhmmer"), "SSU": dict(gene="SSU", exons=1, model="barrnap")}
        models = {"CYTB": dict(gene="CYTB", partial=False, model="miniprot"),
                  "COX1": dict(gene="COX1", partial=False, model="miniprot"),
                  "nad5": dict(gene="nad5", partial=True, model="miniprot"),
                  "atp6": dict(gene="atp6", partial=True, model="miniprot")}
        got = organelle.with_miniprot(found, models)
        self.assertEqual({g: r["model"] for g, r in got.items()},
                         {"CYTB": "nhmmer", "COX1": "miniprot", "nad5": "nhmmer", "SSU": "barrnap", "atp6": "miniprot"})

    def test_coding_sequences_proteins_and_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        # a yeast's COX1 by miniprot: two exons on the minus strand of its record, read with its code (3: TGA is W)
        cds = "ATG" + "TGA" * 20 + "AAA" * 20 + "TAA"
        exon1, exon2 = cds[:60], cds[60:]
        molecule = "C" * 50 + rrna.revcomp(exon2) + "G" * 30 + rrna.revcomp(exon1) + "C" * 50
        start, end = 51, 50 + len(cds) + 30
        record = organelle.subsequence(molecule, start, end, "-")
        parts = [(end - 59, end, "-"), (start, start + len(exon2) - 1, "-")]
        row = dict(gene="COX1", start=start, end=end, parts=organelle.parts_text(parts), exons=2)
        self.assertEqual(organelle.coding_of(row, record), cds)
        self.assertIsNone(organelle.coding_of(dict(row, parts="-"), record))  # nhmmer's hits joined: no exons
        self.assertIsNone(organelle.coding_of(dict(gene="SSU"), "ACGT"))
        records = [(f"m{i}_COX1", record, dict(record=f"m{i}_COX1", gene="COX1", organelle="mitochondrion",
                                                taxid=t, species_taxid=t, partial="f", ambiguous="0.0000",
                                                model="miniprot", exons=2, cds=cds))
                   for i, t in enumerate((4932, 27291))]
        organelle.check(records, ncbi)  # no vsearch: no own-species check
        y = records[0][2]
        self.assertEqual((y["protein"], y["internal_stops"]), ("M" + "W" * 20 + "K" * 20, 0))  # yeast's code: W
        # no vsearch: undecided, unverified (since 2026-11-09; medium until then)
        self.assertEqual((y["confidence"], y["flags"]), ("unverified", "own_species_undecided"))
        self.assertEqual(organelle.confidence(dict(y, internal_stops=2, partial="t"))[1],
                         "internal_stops,partial,own_species_undecided")


class AuditTest(unittest.TestCase):
    """Since 2026-11-10 (the audit's C19, C22, L10): a gene's best feature and no pseudogene, a CDS read by its
    /codon_start and /transl_table, the spacer across a circular record's origin, the caches' keys."""

    def record(self, features, seq=None, circular=True):
        seq = seq or random_seq(6000, 11)
        return dict(seqid="NC_8.1", taxid=952, length=len(seq), circular=circular, sequence=seq, features=features)

    def test_a_genes_best_feature_and_no_pseudogene(self):
        rec = self.record([
            # a truncated ycf1 at the inverted repeat's border first, the gene later, a longer pseudogene between
            ("CDS", ["complement(100..400)"], [("gene", "ycf1")]),
            ("CDS", ["500..2999"], [("gene", "ycf1"), ("pseudo", "")]),
            ("CDS", ["3000..4499"], [("gene", "ycf1")]),
            # a partial 16S (longer) before a whole one; a 16S pseudogene
            ("rRNA", ["<4500..5600"], [("gene", "rrn16")]),
            ("rRNA", ["5601..5900"], [("gene", "rrn16")]),
            ("rRNA", ["1..5999"], [("gene", "rrn16"), ("pseudogene", "unknown")]),
            # two of one length (an inverted repeat's): the first, as before
            ("rRNA", ["complement(410..490)"], [("product", "5S ribosomal RNA")]),
            ("rRNA", ["5901..5981"], [("product", "5S ribosomal RNA")]),
            ("CDS", ["10..90"], [("gene", "psbA"), ("pseudo", "")])])  # a pseudogene alone: no gene
        seq = rec["sequence"]
        got = organelle.record_genes(rec, "plastid")
        self.assertEqual(sorted(got), ["5S", "SSU", "ycf1"])
        self.assertEqual({g: (r["start"], r["end"], r["strand"], r["copies"], r["partial"]) for g, r in got.items()},
                         {"ycf1": (3000, 4499, "+", 2, False), "SSU": (5601, 5900, "+", 2, False),
                          "5S": (410, 490, "-", 2, False)})
        self.assertEqual(got["ycf1"]["sequence"], seq[2999:4499])
        # partial features alone: the longest of them
        rec = self.record([("rRNA", ["<1..300"], [("gene", "rrn23")]), ("rRNA", ["400..>1500"], [("gene", "rrn23")])])
        lsu = organelle.record_genes(rec, "plastid")["LSU"]
        self.assertEqual((lsu["start"], lsu["end"], lsu["partial"], lsu["copies"]), (400, 1500, True, 2))

    def test_a_cds_by_its_codon_start(self):
        coding = "AAA" * 30 + "TAA"
        seq = random_seq(100, 12) + "CA" + coding + random_seq(100, 13) + rrna.revcomp("G" + coding) + random_seq(50, 4)
        a, b = 101, 102 + len(coding)  # the plus strand's feature: two bases of a codon, then the coding sequence
        c, d = b + 101, b + 101 + len(coding)  # the minus strand's: one base, on the molecule its last
        rec = self.record([
            ("CDS", [f"<{a}..{b}"], [("gene", "nad1"), ("codon_start", "3"), ("transl_table", "4")]),
            ("CDS", [f"complement({c}..>{d})"], [("gene", "nad5"), ("codon_start", "2")]),
            ("CDS", [f"join(<{a}..{a + 10},{a + 20}..{a + 30})"], [("gene", "atp6"), ("codon_start", "2")]),
            ("CDS", [f"{a}..{b}"], [("gene", "COX3"), ("codon_start", "1")])], seq)
        got = organelle.record_genes(rec, "mitochondrion")
        for gene, start, end in (("nad1", a + 2, b), ("nad5", c, d - 1)):
            row = got[gene]
            self.assertEqual((row["start"], row["end"], row["spliced"], row["sequence"]), (start, end, coding, coding))
            self.assertEqual(gencode.protein_of(row["spliced"], 1), ("K" * 30, 0))  # in frame
        self.assertEqual((got["nad1"]["parts"], got["nad5"]["parts"]), (f"{a + 2}-{b}:+", f"{c}-{d - 1}:-"))
        self.assertEqual(got["atp6"]["parts"], f"{a + 1}-{a + 10}:+,{a + 20}-{a + 30}:+")  # off its first exon alone
        self.assertEqual((got["COX3"]["start"], got["COX3"]["parts"]), (a, f"{a}-{b}:+"))  # 1: as it is
        self.assertEqual([got[g]["transl_table"] for g in ("nad1", "nad5")], ["4", "-"])
        self.assertEqual(organelle.in_frame([(5, 6, "+")], [("codon_start", "3")]), [(5, 6, "+")])  # nothing left: kept

    def test_a_cds_by_its_own_transl_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        row = dict(taxid=27291, organelle="mitochondrion")  # NCBI's mitochondrial code of the species: 1
        self.assertEqual([organelle.table_of(dict(row, transl_table=t), ncbi) for t in ("4", "-", "", "99", "x")],
                         [4, 1, 1, 1, 1])
        self.assertEqual(organelle.table_of(dict(taxid=4932, organelle="mitochondrion"), ncbi), 3)  # an assembly's
        cds = "ATG" + "TGA" * 5 + "AAA" * 5 + "TAA"  # TGA: W by table 4, a stop by table 1
        records = [(f"m{i}_COX1", cds, dict(record=f"m{i}_COX1", gene="COX1", organelle="mitochondrion", taxid=27291,
                                            species_taxid=27291, partial="f", ambiguous="0.0000", model="refseq",
                                            exons=1, cds=cds, transl_table=t)) for i, t in enumerate(("4", "-"))]
        organelle.check(records, ncbi)
        self.assertEqual([(r["protein"], r["internal_stops"]) for _n, _s, r in records],
                         [("M" + "W" * 5 + "K" * 5, 0), ("M" + "*" * 5 + "K" * 5, 5)])
        # the assemblies' protein queries: a RefSeq gene read by its own table is whole, by the species' it is not
        chosen = {(27291, "mitochondrion"): [dict(seqid=f"NC_{i}.1", gene="COX1", taxid="27291", partial="f",
                                                  trans_spliced="f", spliced=cds, organelle="mitochondrion",
                                                  transl_table=t) for i, t in enumerate(("4", "-"))]}
        self.assertEqual(sorted(organelle.Nearest(chosen, ncbi).protein), [("NC_0.1", "COX1")])

    def test_the_spacer_across_the_origin(self):
        # test_spacer_and_a_trans_spliced_gene's record opened 300 bases on: between trnH and psbA's 3' end
        seq = random_seq(4000, 9)
        turned = seq[300:] + seq[:300]
        features = [("tRNA", ["complement(3701..3775)"], [("gene", "trnH-GUG"), ("product", "tRNA-His")]),
                    ("CDS", ["complement(100..1159)"], [("gene", "psbA")])]
        spacer = organelle.record_genes(self.record(features, turned), "plastid")["trnH-psbA"]
        self.assertEqual((spacer["start"], spacer["end"], spacer["strand"], spacer["regions"], spacer["parts"]),
                         (3701, 299, "-", "psbA3':0-200,spacer:200-524,trnH:524-599", "1-299:-,3701-4000:-"))
        self.assertEqual(spacer["sequence"], rrna.revcomp(seq[:599]))  # the record's, as before it was opened there
        self.assertNotIn("trnH-psbA", organelle.record_genes(self.record(features, turned, circular=False), "plastid"))
        # on the plus strand: psbA at the record's end, trnH at its start
        features = [("CDS", ["3000..3990"], [("gene", "psbA")]),
                    ("tRNA", ["100..174"], [("gene", "trnH-GUG"), ("product", "tRNA-His")])]
        spacer = organelle.record_genes(self.record(features, seq), "plastid")["trnH-psbA"]
        self.assertEqual((spacer["start"], spacer["end"], spacer["regions"], spacer["parts"]),
                         (3791, 174, "psbA3':0-200,spacer:200-309,trnH:309-384", "3791-4000:+,1-174:+"))
        self.assertEqual(spacer["sequence"], seq[3790:] + seq[:174])
        # psbA itself across the origin (until 2026-11-10: no spacer)
        features = [("CDS", ["join(3900..4000,1..900)"], [("gene", "psbA")]),
                    ("tRNA", ["1000..1074"], [("gene", "trnH-GUG"), ("product", "tRNA-His")])]
        spacer = organelle.record_genes(self.record(features, seq), "plastid")["trnH-psbA"]
        self.assertEqual((spacer["start"], spacer["end"], spacer["sequence"]), (701, 1074, seq[700:1074]))
        # a trnH on the other side of psbA, or too far beyond it: no neighbour, on a circular record neither
        for trnh in ("2000..2074", "complement(100..174)"):
            features = [("CDS", ["3000..3990"], [("gene", "psbA")]), ("tRNA", [trnh], [("gene", "trnH-GUG")])]
            self.assertNotIn("trnH-psbA", organelle.record_genes(self.record(features, seq), "plastid"))

    def test_caches_keyed_by_their_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            src, work = os.path.join(tmp, "sources"), os.path.join(tmp, "work", "organelle")
            release(src)

            def caches():
                return sorted(f for f in os.listdir(work) if f.startswith("refseq_"))
            rows = organelle.refseq_genes(src, work, ncbi, lambda _t: 0)
            first = caches()
            self.assertEqual((len(first), {r["transl_table"] for r in rows}), (1, {"-"}))
            genomes_dir = os.path.join(tmp, "genomes")
            os.makedirs(genomes_dir)
            with gzip.open(genomes.genome_file(genomes_dir, YEAST), "wt") as fh:
                fh.write(f">m1 Saccharomyces cerevisiae mitochondrion, complete genome\n{YEAST_MT}\n")
            jobs = {YEAST: dict(accession=YEAST, species_taxid=4932, taxid=559292)}
            barrnap, nhmmer = fx.tool(tmp, "barrnap", BARRNAP), fx.tool(tmp, "nhmmer", NHMMER)

            def run(rows):
                return organelle.assembly_pass(work, jobs, lambda a: genomes.genome_file(genomes_dir, a),
                                               organelle.choose(rows, ncbi), ncbi, barrnap, nhmmer, 1, 1, lambda _t: 0)
            self.assertEqual((run(rows), run(rows)), ((1, 0), (0, 0)))
            # another release under the same names (release.txt as it was): parsed again, not read from the cache;
            # the yeast's queries are another species' now (a nearer one's COX1), its assembly searched again
            other = random_seq(2000, 21)
            with gzip.open(os.path.join(organelle.folder_of(src), "mitochondrion.1.genomic.gbff.gz"), "at") as fh:
                fh.write(gbff_record("NC_000005.1", other, 27291, "mitochondrion", [
                    feature("CDS", "101..1600", '/gene="COX1"', "/transl_table=3")]))
            rows = organelle.refseq_genes(src, work, ncbi, lambda _t: 0)
            self.assertEqual(len(caches()), 2)
            self.assertEqual({r["seqid"]: r["transl_table"] for r in rows if r["seqid"] == "NC_000005.1"},
                             {"NC_000005.1": "3"})
            self.assertEqual((run(rows), run(rows)), ((1, 0), (0, 0)))
            # the cache's version raised (the parsing or the search changed): both made again
            organelle.CACHE_VERSION += 1
            try:
                self.assertEqual(len(organelle.refseq_genes(src, work, ncbi, lambda _t: 0)), len(rows))
                self.assertEqual((len(caches()), run(rows), run(rows)), (3, (1, 0), (0, 0)))
            finally:
                organelle.CACHE_VERSION -= 1
            self.assertEqual(len(organelle.refseq_genes(src, work, ncbi, lambda _t: 0)), len(rows))
            self.assertEqual(len(caches()), 3)  # the version as it was: its cache read again, none made
            # a rule changed: another key without the version raised
            before = organelle.rules_text()
            organelle.SPACER_MAX += 1
            try:
                self.assertNotEqual(organelle.rules_text(), before)
            finally:
                organelle.SPACER_MAX -= 1


if __name__ == "__main__":
    unittest.main()
