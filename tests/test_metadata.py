#!/usr/bin/env python3
"""Tests of Eukore's assembly metadata (eukore/metadata.py, since 2026-11-07): one pass over a genome's FASTA (its
lengths, N50s, contigs split at gaps, GC and its windows, soft-masking, CpG, IUPAC bases, duplicate sequences,
telomeres at both ends, organelle molecules, canonical tetranucleotides) checked against counts made the slow way and
the same in any chunking and decompressor; a genome's metadata into its done folder with a stand-in sourmash; the
package's metadata/ and organelle/molecules.fna.gz.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import random
import shutil
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, metadata, ncbi as ncbi_mod, package, util  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST, YEAST_GCF, OTHER = "GCA_000146045.2", "GCF_000146045.2", "GCA_000000003.1"


def rand(n, seed, alphabet="ACGT"):
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(n))


CHR2 = rand(15_000, 2)[:7000] + "RY" + rand(15_000, 2)[7002:]
RECORDS = [
    ("chr1 Testus primus chromosome 1",
     "CCCTAA" * 100 + rand(9000, 1) + "N" * 50 + rand(4000, 3).lower() + "NNNNN" + rand(8000, 4) + "TTAGGG" * 100),
    ("chr2 Testus primus chromosome 2", CHR2),
    ("chr3 Testus primus chromosome 3", CHR2),  # the same as chr2: a duplicate
    ("s1 unplaced", "GGCC" * 125),  # 500 bases
    ("mt1 Testus primus mitochondrion, complete genome", rand(3000, 5, "AATTGC")),
]
SKETCHER = r"""#!/usr/bin/env python3
import gzip, sys
a = sys.argv
with gzip.open(a[a.index("-o") + 1], "wt") as fh:
    fh.write('[{"name": "%s"}]' % a[a.index("--name") + 1])
"""


def fasta(path, records=RECORDS, width=60):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as fh:
        for title, seq in records:
            fh.write(f">{title}\n" + "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width)))
    return path


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def slow_tetranucleotides(seqs):
    counts = dict.fromkeys(metadata.CANONICAL, 0)
    for s in seqs:
        s = s.upper()
        for i in range(len(s) - 3):
            k = s[i:i + 4]
            if set(k) <= set("ACGT"):
                counts[min(k, revcomp(k))] += 1
    return [counts[k] for k in metadata.CANONICAL]


class ScanTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = fasta(os.path.join(self.tmp.name, "g.fna.gz"))

    def tearDown(self):
        self.tmp.cleanup()

    def test_canonical_tetranucleotides(self):
        self.assertEqual(len(metadata.CANONICAL), 136)
        self.assertIn("AAAA", metadata.CANONICAL)
        self.assertNotIn("TTTT", metadata.CANONICAL)

    def test_measures_against_slow_counts(self):
        stats, tetra, molecules = metadata.scan(self.path)
        seqs = [s for _t, s in RECORDS]
        up = "".join(seqs).upper()
        acgt = sum(up.count(b) for b in "ACGT")
        lengths = sorted(map(len, seqs), reverse=True)
        self.assertEqual((stats["sequences"], stats["length"], stats["longest"]), (5, sum(lengths), lengths[0]))
        self.assertEqual((stats["acgt"], stats["n_bases"], stats["other_bases"]), (acgt, 55, 4))
        n50 = next(n for i, n in enumerate(lengths) if sum(lengths[:i + 1]) >= sum(lengths) / 2)
        self.assertEqual((stats["n50"], stats["l50"]), (n50, lengths.index(n50) + 1))
        # contigs: chr1 split at its 50 Ns, not at its 5
        self.assertEqual((stats["contigs"], stats["gaps"], stats["gap_bases"]), (6, 1, 50))
        self.assertEqual(stats["contig_n50"], 15_000)
        self.assertEqual(stats["gc"], round((up.count("G") + up.count("C")) / acgt, 4))
        self.assertEqual(stats["masked"], round(4000 / acgt, 4))
        cg = sum(s.upper().count("CG") for s in seqs)
        self.assertEqual(stats["cpg_oe"], round(cg * acgt / (up.count("C") * up.count("G")), 4))
        self.assertEqual((stats["duplicate_sequences"], stats["duplicate_bases"]), (1, 15_000))
        self.assertEqual(stats["short_1kb"], round(500 / sum(lengths), 4))
        # chr1: telomeric repeats at both ends (CCCTAA: TTAGGG's other strand); chr2, chr3 none
        self.assertEqual((stats["telomere_motif"], stats["telomeric_ends"], stats["telomere_both"]),
                         ("TTAGGG", 2, 1))
        self.assertEqual((stats["mitochondrion_molecules"], stats["mitochondrion_length"], stats["plastid_molecules"]),
                         (1, 3000, 0))
        self.assertEqual(molecules, [("mt1", "mitochondrion", RECORDS[4][1])])
        self.assertEqual(tetra, slow_tetranucleotides(seqs))
        # 10 kb windows of half ACGT or more: chr1 2 (its last 2,255 bases too few), chr2 and chr3 2 each (their last
        # 5,000 enough), mt1 and s1 none
        self.assertEqual(stats["gc_windows"], 6)

    def test_the_same_in_any_chunks_and_decompressor(self):
        got = metadata.scan(self.path)
        chunk, tool = metadata.CHUNK, util.GZIP_TEST
        try:
            metadata.CHUNK = 7  # pieces and chunks cut everywhere: in titles, gaps, windows, tetranucleotides
            self.assertEqual(metadata.scan(self.path), got)
            metadata.CHUNK = chunk
            util.GZIP_TEST = shutil.which("gzip")  # gzip -dc as igzip's stand-in
            self.assertEqual(metadata.scan(self.path), got)
            util.GZIP_TEST = None  # Python's gzip
            self.assertEqual(metadata.scan(self.path), got)
        finally:
            metadata.CHUNK, util.GZIP_TEST = chunk, tool
        plain = fasta(os.path.join(self.tmp.name, "g.fna"), width=100_000)  # plain, one line a sequence
        self.assertEqual(metadata.scan(plain), got)

    def test_a_cut_gzip_fails(self):
        with open(self.path, "rb") as fh:
            data = fh.read()
        cut = os.path.join(self.tmp.name, "cut.fna.gz")
        with open(cut, "wb") as fh:
            fh.write(data[:len(data) // 2])
        tool = util.GZIP_TEST
        try:
            for util.GZIP_TEST in (shutil.which("gzip"), None):
                with self.assertRaises((OSError, EOFError)):
                    metadata.scan(cut)
        finally:
            util.GZIP_TEST = tool


class PackageTest(unittest.TestCase):
    def test_a_genomes_metadata_into_the_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            path = fasta(os.path.join(tmp, f"{YEAST}.fna.gz"))
            work = os.path.join(tmp, "work")
            folder = genomes.done_dir(work, YEAST)
            os.makedirs(folder)
            sourmash = fx.tool(tmp, "sourmash", SKETCHER)
            self.assertEqual(metadata.missing(folder, sourmash), {metadata.STATS, metadata.SKETCH})
            self.assertEqual(metadata.describe_task((YEAST, path, folder, YEAST_GCF, sourmash)), (YEAST, None))
            self.assertEqual(metadata.missing(folder, sourmash), set())
            self.assertEqual(sorted(os.listdir(folder)), sorted([metadata.STATS, metadata.ORGANELLE, metadata.SKETCH]))
            with open(os.path.join(folder, metadata.STATS)) as fh:
                stats = json.load(fh)
            self.assertEqual((stats["fetched"], stats["md5"]), (YEAST_GCF, util.md5(path)))
            # a sketch that fails: said, the measures kept, nothing half-written left
            broken = fx.tool(tmp, "broken", "#!/bin/sh\necho no >&2\nexit 3\n")
            os.remove(os.path.join(folder, metadata.SKETCH))
            why = metadata.describe(path, folder, YEAST, YEAST_GCF, broken, {metadata.SKETCH})
            self.assertTrue(why.startswith("sourmash sketch failed (exit 3)"), why)
            self.assertEqual(metadata.missing(folder, broken), {metadata.SKETCH})
            self.assertFalse([f for f in os.listdir(folder) if "partial" in f])
            metadata.describe(path, folder, YEAST, YEAST_GCF, sourmash, {metadata.SKETCH})
            tables = [fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK),
                      fx.assembly_table(os.path.join(tmp, "refseq_fungi.txt"), fx.REFSEQ)]
            # the runs' gene models (done.json): the yeast's Ensembl's, the TOGA2 genome's TOGA2's
            toga2, eukprot = "TOGA2_HLtesSpe1", "EUK_6253"
            for acc, source in ((YEAST, "ensembl"), (toga2, "toga2"), (eukprot, "-")):
                os.makedirs(genomes.done_dir(work, acc), exist_ok=True)
                with open(os.path.join(genomes.done_dir(work, acc), "done.json"), "w") as fh:
                    json.dump({"accession": acc, "annotation": source, "annotation_none": ["ncbi"]}, fh)
            jobs = {YEAST: dict(accession=YEAST, fetch=YEAST_GCF, why=["ncbi"], url="https://ftp.example/y.fna.gz"),
                    toga2: dict(accession=toga2, fetch=toga2, why=["toga2"], assembly_source="DNA Zoo",
                                url="https://toga2.example/HLtesSpe1.2bit"),
                    eukprot: dict(accession=eukprot, fetch=eukprot, why=["no_assembly"], url=None, path="x.fna")}
            collected = dict(genomes={
                YEAST: dict(accession=YEAST, fetched=YEAST_GCF, species_taxid=4932, taxid=559292, organism="y",
                            source="ncbi", representative="t", loci=0),
                OTHER: dict(accession=OTHER, fetched=OTHER, species_taxid=4932, taxid=4932, organism="y",
                            source="assembly", representative="f", loci=0),
                toga2: dict(accession=toga2, fetched=toga2, species_taxid=9606, taxid=9606, organism="h",
                            source="toga2", representative="t", loci=0),  # not described: no stats.json
                eukprot: dict(accession=eukprot, fetched=eukprot, species_taxid=6253, taxid=6253, organism="w",
                              source="no_assembly", representative="t", loci=0)},
                copies={YEAST: {}, toga2: {}, eukprot: {}}, features=[], copy_rows=[], taxids=set())
            hook = lambda part, g: metadata.write(part, g, work, tables, jobs)  # noqa: E731
            out = package.write(os.path.join(tmp, "out", "Eukore_2026.11"), "2026.11", ncbi, collected, [],
                                metadata=hook)
            self.assertIsNone(package.check(out, deep=True))
            with open(os.path.join(out, "manifest.json")) as fh:
                counts = json.load(fh)["counts"]
            self.assertEqual(counts["metadata"], {
                "ncbi_rows": 3, "genomes": 3, "assemblies_described": 1, "sketches": 1, "organelle_molecules": 1,
                "genome_sources": {"EukProt": 1, "NCBI RefSeq": 1, "TOGA2 (DNA Zoo)": 1},
                "annotation_sources": {"-": 1, "Ensembl": 1, "TOGA2": 1}})
            ncbi_rows = package.read_table(out, "metadata/ncbi.tsv.gz")
            self.assertEqual([(r["genome"], r["assembly_accession"], r["fetched"]) for r in ncbi_rows],
                             [(OTHER, OTHER, "t"), (YEAST, YEAST, "f"), (YEAST, YEAST_GCF, "t")])
            self.assertEqual(list(ncbi_rows[0])[:3], ["genome", "fetched", "assembly_accession"])
            self.assertEqual(ncbi_rows[2]["refseq_category"], "reference genome")
            rows = {r["accession"]: r for r in package.read_table(out, "metadata/assembly.tsv.gz")}
            self.assertEqual((len(rows), list(rows[YEAST])), (3, list(metadata.ASSEMBLY_COLUMNS)))
            y = rows[YEAST]
            self.assertEqual((y["fetched"], y["genome_source"], y["genome_url"], y["annotation_source"], y["n50"],
                              y["sketch"]),
                             (YEAST_GCF, "NCBI RefSeq", "https://ftp.example/y.fna.gz", "Ensembl", str(stats["n50"]),
                              "t"))
            t = rows[toga2]  # where it came from, though not measured
            self.assertEqual((t["genome_source"], t["genome_url"], t["annotation_source"], t["md5"], t["n50"]),
                             ("TOGA2 (DNA Zoo)", "https://toga2.example/HLtesSpe1.2bit", "TOGA2", "-", "-"))
            self.assertEqual((rows[eukprot]["genome_source"], rows[eukprot]["genome_url"]), ("EukProt", "-"))
            self.assertEqual(metadata.genome_source(dict(accession="GCA_1.1", fetched="GCA_1.1")), "NCBI GenBank")
            tetra = package.read_table(out, "metadata/tetranucleotides.tsv.gz")
            self.assertEqual(len(tetra[0]), 137)
            self.assertEqual([int(tetra[0][k]) for k in metadata.CANONICAL], stats["tetranucleotides"])
            with zipfile.ZipFile(os.path.join(out, "metadata", "sketches.zip")) as z:
                self.assertEqual(z.namelist(), [f"signatures/{YEAST}.sig.gz"])
            with gzip.open(os.path.join(out, "organelle", "molecules.fna.gz"), "rt") as fh:
                lines = fh.read().split("\n")
            self.assertEqual((lines[0], lines[1]), (f">{YEAST}|mt1 mitochondrion", RECORDS[4][1]))


if __name__ == "__main__":
    unittest.main()
