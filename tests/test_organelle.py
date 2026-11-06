#!/usr/bin/env python3
"""Tests of Eukore's organelle genes (eukore/organelle.py): RefSeq's GenBank records read (a COX1 across a circular
genome's origin on the minus strand, a qualifier over two lines, an animal's 16S as its LSU, a plastid's rRNA twice,
its 4.5S left out), a species' record chosen, an assembly's mitochondrion searched (stand-in barrnap and nhmmer: its
rRNA by bacterial models, its COX1 of two exons), and the genes written into the package.

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
PLASTID = random_seq(3000, 3)
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
            feature("rRNA", "J01415.2:1..10", '/product="12S ribosomal RNA"')]))
        fh.write(gbff_record("NC_000002.1", HUMAN_MT2, 9606, "mitochondrion", [
            feature("rRNA", "101..500", '/product="12S ribosomal RNA"')]))
    with gzip.open(os.path.join(folder, "plastid.1.genomic.gbff.gz"), "wt") as fh:
        fh.write(gbff_record("NC_000003.1", PLASTID, 28985, "plastid:chloroplast", [
            feature("rRNA", "101..1600", '/gene="rrn16"'),
            feature("rRNA", "complement(2001..2400)", '/product="16S ribosomal RNA"'),
            feature("rRNA", "2501..2600", '/product="4.5S ribosomal RNA"'),
            feature("rRNA", "2701..2820", '/product="5S ribosomal RNA"')]))
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
first = next(n for n in seqs if n.startswith("COX1|")); half = len(seqs[first]) // 2
with open(tbl, "w") as fh:
    fh.write("# target query\n")
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

    def test_refseq_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            release(os.path.join(tmp, "sources"))
            rows = organelle.refseq_genes(os.path.join(tmp, "sources"), os.path.join(tmp, "work"), ncbi, lambda _t: 0)
            got = {(r["seqid"], r["gene"]): r for r in rows}
            self.assertEqual(sorted(got), [("NC_000001.1", "COX1"), ("NC_000001.1", "LSU"), ("NC_000001.1", "SSU"),
                                           ("NC_000002.1", "SSU"), ("NC_000003.1", "5S"), ("NC_000003.1", "SSU")])
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
                             ["SSU|NC_000001.1", "LSU|NC_000001.1", "COX1|NC_000001.1"])  # the nearest with each gene


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
                                            (9606, "COX1"), (9606, "LSU"), (9606, "SSU")])
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
            self.assertEqual((len(table), set(table[0])), (6, set(organelle.GENE_COLUMNS)))
            self.assertEqual(package.read_table(out, "species.tsv")[0]["organelle_genes"], "3")
            with open(os.path.join(out, "manifest.json")) as fh:
                counts = json.load(fh)["counts"]
            self.assertEqual(counts["organelle_genes"]["mitochondrion"], {"COX1": 2, "LSU": 2, "SSU": 2})
            with gzip.open(os.path.join(out, "organelle", "genes.fna.gz"), "rt") as fh:
                self.assertEqual(fh.read().count(">"), 6)


if __name__ == "__main__":
    unittest.main()
