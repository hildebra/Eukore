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


# sourmash's stand-in that does not end: its partial half-written, its pid and its partial's name noted beside it
STALLED = r"""#!/usr/bin/env python3
import os, sys, time
a = sys.argv
part = a[a.index("-o") + 1]
with open(part, "w") as fh:
    fh.write("half")
note = os.path.join(os.path.dirname(os.path.abspath(a[0])), "noted")
with open(note + ".tmp", "w") as fh:
    fh.write("%d\t%s" % (os.getpid(), part))
os.replace(note + ".tmp", note)
time.sleep(60)
"""


class DescribeTest(unittest.TestCase):
    """metadata.describe since 2026-11-10: the MD5 from the scan's own pass, a damaged gzip a genome's failure, the
    sketch's partial this process' own and never left."""

    def test_the_md5_comes_from_the_scan(self):
        import hashlib  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            gz = fasta(os.path.join(tmp, "g.fna.gz"))
            plain = fasta(os.path.join(tmp, "g.fna"))
            bgzf = os.path.join(tmp, "b.fna.gz")  # BGZF: a gzip member a block, and the empty one at its end
            util.bgzf_copy(plain, bgzf, util.WORK_LEVEL)
            # 4 Mb: a compressed file of more than one of the feeder's blocks and many of gzip's reads
            big = os.path.join(tmp, "big.fna.gz")
            with gzip.open(big, "wb", compresslevel=1) as fh:
                fh.write(b">big1\n" + os.urandom(4_000_000).translate(bytes(b"ACGT"[i & 3] for i in range(256))))
            self.assertGreater(os.path.getsize(big), 1 << 20)

            def md5(path):
                with open(path, "rb") as fh:
                    return hashlib.md5(fh.read()).hexdigest()
            tools = [t for t in (shutil.which("gzip"), shutil.which("pigz"), shutil.which("igzip")) if t] + [None]
            tool, expected = util.GZIP_TEST, metadata.scan(gz)
            try:
                for path in (gz, bgzf, plain):
                    for util.GZIP_TEST in tools if path != plain else [None]:  # a plain file: no decompressor
                        digest = hashlib.md5()
                        self.assertEqual(metadata.scan(path, digest), expected)  # the scan the same with a digest
                        self.assertEqual(digest.hexdigest(), md5(path), (path, util.GZIP_TEST))
                for util.GZIP_TEST in tools:
                    digest = hashlib.md5()
                    for _piece in metadata._pieces_of(big, digest):  # noqa: SLF001
                        pass
                    self.assertEqual(digest.hexdigest(), md5(big), util.GZIP_TEST)
                # describe takes it from there: the file is not read a second time
                for util.GZIP_TEST in tools:
                    folder = os.path.join(tmp, f"done-{tools.index(util.GZIP_TEST)}")
                    os.makedirs(folder)
                    with mock.patch.object(util, "md5", side_effect=AssertionError("a second read")):
                        self.assertIsNone(metadata.describe(gz, folder, YEAST, YEAST_GCF))
                    with open(os.path.join(folder, metadata.STATS)) as fh:
                        stats = json.load(fh)
                    self.assertEqual((stats["md5"], stats["bytes"]), (md5(gz), os.path.getsize(gz)))
                # a scan given up (its reader's error) ends the tool and its feeder: nothing waits
                util.GZIP_TEST = shutil.which("gzip")
                pieces = metadata._pieces_of(big, hashlib.md5())  # noqa: SLF001
                next(pieces)
                pieces.close()
            finally:
                util.GZIP_TEST = tool

    def test_a_damaged_gzip_is_the_genomes_failure(self):
        import zlib  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            with open(fasta(os.path.join(tmp, "g.fna.gz")), "rb") as fh:
                data = bytearray(fh.read())
            data[len(data) // 2:len(data) // 2 + 64] = b"\xff" * 64  # the deflate stream damaged, not cut
            bad = os.path.join(tmp, "bad.fna.gz")
            with open(bad, "wb") as fh:
                fh.write(data)
            with self.assertRaises(zlib.error), gzip.open(bad, "rb") as fh:
                fh.read()
            folder = os.path.join(tmp, "done")
            os.makedirs(folder)
            sourmash = fx.tool(tmp, "sourmash", SKETCHER)
            tool = util.GZIP_TEST
            try:
                util.GZIP_TEST = None  # Python's gzip: zlib.error, where a tool's failure is an OSError
                why = metadata.describe(bad, folder, YEAST, YEAST_GCF, sourmash)
            finally:
                util.GZIP_TEST = tool
            self.assertTrue(why.startswith("metadata: error: "), why)
            self.assertEqual(os.listdir(folder), [metadata.SKETCH])  # no stats of it, nothing half-written

    def test_the_sketchs_partial_is_this_process_own_and_never_left(self):
        import socket  # noqa: PLC0415
        import subprocess  # noqa: PLC0415
        import time  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            path = fasta(os.path.join(tmp, f"{YEAST}.fna.gz"))
            folder = os.path.join(tmp, "done")
            os.makedirs(folder)
            host = socket.gethostname()
            gone = subprocess.Popen(["true"])
            gone.wait()
            # another process' partials: the name every process wrote until 2026-11-10, another host's of now, and
            # one of a process of this host that is gone
            old, theirs, dead = (os.path.join(folder, n) for n in (
                ".sketch.partial.sig.gz", ".elsewhere-7.sketch.partial.sig.gz",
                f".{host}-{gone.pid}.sketch.partial.sig.gz"))
            for other in (old, theirs, dead):
                with open(other, "w") as fh:
                    fh.write("another's")
            sourmash = fx.tool(tmp, "sourmash", SKETCHER)
            self.assertIsNone(metadata.describe(path, folder, YEAST, YEAST_GCF, sourmash, {metadata.SKETCH}))
            self.assertEqual(sorted(os.listdir(folder)), sorted([os.path.basename(old), os.path.basename(theirs),
                                                                 metadata.SKETCH]))
            for other in (old, theirs):  # not written over, not renamed away
                with open(other) as fh:
                    self.assertEqual(fh.read(), "another's")
            os.remove(os.path.join(folder, metadata.SKETCH))
            # an error nobody foresaw while sourmash runs: sourmash ended, its partial removed
            stalled = fx.tool(tmp, "stalled", STALLED)
            note = os.path.join(tmp, "noted")

            def failing(*_args):
                for _ in range(400):
                    if os.path.isfile(note):
                        break
                    time.sleep(0.05)
                raise RuntimeError("the scan broke")
            with mock.patch.object(metadata, "scan", failing), self.assertRaises(RuntimeError):
                metadata.describe(path, folder, YEAST, YEAST_GCF, stalled)
            with open(note) as fh:
                pid, part = fh.read().split("\t")
            self.assertEqual(part, os.path.join(folder, f".{host}-{os.getpid()}.sketch.partial.sig.gz"))
            self.assertFalse(util.pid_alive(int(pid)))
            self.assertEqual(sorted(os.listdir(folder)), sorted([os.path.basename(old), os.path.basename(theirs)]))


if __name__ == "__main__":
    unittest.main()
