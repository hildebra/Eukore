#!/usr/bin/env python3
"""Tests of Eukore's organelle genes (eukore/organelle.py): RefSeq's GenBank records read (a COX1 across a circular
genome's origin on the minus strand, a qualifier over two lines, an animal's 16S as its LSU, a plastid's rRNA twice,
its 4.5S left out), a species' record chosen, an assembly's mitochondrion searched (stand-in barrnap and nhmmer: its
rRNA by bacterial models, its COX1 of two exons), a mitochondrion's CYTB and a plastid's rbcL, matK and tufA (since
2026-11-07; an assembly searched before them searched again, one searched for them all not), and the genes written into
the package.

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
from eukore import genomes, ncbi as ncbi_mod, organelle, package, rrna  # noqa: E402
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
        self.assertTrue(organelle.is_cox1([("product", "cytochrome c oxidase subunit I")]))
        self.assertFalse(organelle.is_cox1([("product", "cytochrome c oxidase subunit II")]))
        self.assertTrue(organelle.is_cox1([("gene", "COI")]))
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
                                  (worm, "#searched\tmitochondrion:" + ",".join(organelle.GENES_OF["mitochondrion"])
                                   + "\n")):  # every gene of its mitochondrion: not searched again
                os.makedirs(os.path.dirname(organelle.assembly_file(work, acc)), exist_ok=True)
                with open(organelle.assembly_file(work, acc), "w") as fh:
                    fh.write(searched + head)
            self.assertEqual(organelle.searched_of(organelle.assembly_file(work, plant))["plastid"],
                             set(organelle.OLD_GENES))
            barrnap, nhmmer = fx.tool(tmp, "barrnap", BARRNAP), fx.tool(tmp, "nhmmer", NHMMER)

            def run():
                return organelle.assembly_pass(work, jobs, lambda a: genomes.genome_file(genomes_dir, a), chosen,
                                               ncbi, barrnap, nhmmer, 1, 1, lambda _t: 0)
            self.assertEqual(run(), (2, 0))  # the plant and the yeast, not the worm
            self.assertEqual(organelle.searched_of(organelle.assembly_file(work, plant)),
                             {k: set(organelle.GENES_OF[k]) for k in organelle.KINDS})
            self.assertEqual(organelle.searched_of(organelle.assembly_file(work, yeast)),
                             {"mitochondrion": set(organelle.GENES_OF["mitochondrion"])})
            self.assertEqual(run(), (0, 0))
            records = organelle.collect(ncbi, jobs, {952}, chosen, work)
            self.assertIn(("mitochondrion", "COX1"), {(r["organelle"], r["gene"]) for _n, _s, r in records})


if __name__ == "__main__":
    unittest.main()
