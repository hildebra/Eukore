#!/usr/bin/env python3
"""Tests of eukore/depth.py and tools/rdna_depth.py (read-depth copy numbers, since 2026-11-06): which alignments
count; per-base depth over a target's interior; the run taken of a genome's runs (curated, its BioSample's, its
strain's; WGS of genomic DNA, Illumina first); the rDNA units and single-copy BUSCOs mapped to; the first reads of a
stream; the copy numbers and their flags; and, where minimap2 is on PATH, a small genome with 40 rDNA units and 30
mitochondria per haploid genome measured from its simulated reads end to end (staged reads mapped, depth.tsv, the
package's genomes.tsv and manifest).

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import io
import json
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import depth as dp, util  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))
import rdna_depth  # noqa: E402


def paf(q, qlen, qs, qe, strand, t, tlen, ts, te, matches=None, block=None, tp="P"):
    block = block if block is not None else qe - qs
    matches = matches if matches is not None else block
    return "\t".join(map(str, (q, qlen, qs, qe, strand, t, tlen, ts, te, matches, block, 60, f"tp:A:{tp}")))


class AcceptTest(unittest.TestCase):
    def test_a_read_inside_a_target_counts(self):
        self.assertTrue(dp.accept(paf("r", 150, 0, 150, "+", "t", 1000, 100, 250).split("\t"), 0.95))

    def test_too_divergent_does_not(self):
        self.assertFalse(dp.accept(paf("r", 150, 0, 150, "+", "t", 1000, 100, 250, 140, 150).split("\t"), 0.95))

    def test_a_read_over_a_targets_end_counts(self):
        # its first 60 bases lie before the target's start, on either strand
        self.assertTrue(dp.accept(paf("r", 150, 60, 150, "+", "t", 1000, 0, 90).split("\t"), 0.95))
        self.assertTrue(dp.accept(paf("r", 150, 0, 90, "-", "t", 1000, 0, 90).split("\t"), 0.95))
        self.assertTrue(dp.accept(paf("r", 150, 0, 90, "+", "t", 1000, 910, 1000).split("\t"), 0.95))

    def test_a_read_clipped_inside_a_target_does_not(self):
        # a paralog's or a shared domain's read: 60 bases unaligned where the target goes on
        self.assertFalse(dp.accept(paf("r", 150, 60, 150, "+", "t", 1000, 400, 490).split("\t"), 0.95))
        self.assertFalse(dp.accept(paf("r", 150, 60, 150, "-", "t", 1000, 400, 490).split("\t"), 0.95))

    def test_a_long_read_over_a_whole_target_counts(self):
        self.assertTrue(dp.accept(paf("r", 15000, 4000, 7000, "+", "t", 3000, 0, 3000).split("\t"), 0.95))


class DepthTest(unittest.TestCase):
    def test_per_base_depth_and_secondary_alignments_left_out(self):
        lines = [paf("a", 10, 0, 10, "+", "t", 30, 0, 10), paf("b", 10, 0, 10, "+", "t", 30, 5, 15),
                 paf("c", 10, 0, 10, "+", "t", 30, 5, 15, tp="S"), paf("d", 10, 0, 10, "+", "u", 30, 0, 10)]
        cov, counted, mean = dp.depths(lines, {"t": 30}, 0.95)
        self.assertEqual(cov["t"][:16], [1] * 5 + [2] * 5 + [1] * 5 + [0])
        self.assertEqual((counted, mean), (2, 10.0))

    def test_interior_depth(self):
        cov = [0] * 100 + [10] * 300 + [4] * 300 + [0] * 100
        self.assertAlmostEqual(dp.interior_depth(cov, [(0, 800)], 100), 7.0)
        self.assertEqual(dp.interior_depth(cov, [(0, 800)], 100, median=True), 7.0)  # 300 of 10, 300 of 4
        self.assertEqual(dp.interior_depth(cov, [(100, 400)], 100, median=True), 10.0)
        self.assertIsNone(dp.interior_depth(cov, [(0, 250)], 100))  # 150 bases left: too few


def run_row(acc, platform="ILLUMINA", bases=1000, selection="RANDOM", strategy="WGS", source="GENOMIC", model="x"):
    return {"run_accession": acc, "instrument_platform": platform, "instrument_model": model,
            "library_strategy": strategy, "library_source": source, "library_selection": selection,
            "base_count": str(bases), "fastq_ftp": f"ftp.sra.ebi.ac.uk/{acc}_1.fastq.gz;ftp.sra.ebi.ac.uk/{acc}_2.fastq.gz"}


class RunTest(unittest.TestCase):
    def test_pick_run(self):
        rows = [run_row("cdna", selection="cDNA", bases=9e9), run_row("hic", strategy="Hi-C", bases=9e9),
                run_row("rna", source="TRANSCRIPTOMIC", bases=9e9), run_row("pb", "PACBIO_SMRT", 9e9),
                run_row("small", bases=10), run_row("big", bases=5000)]
        self.assertEqual(dp.pick_run(rows)["run_accession"], "big")
        self.assertEqual(dp.pick_run(rows[:4])["run_accession"], "pb")
        self.assertIsNone(dp.pick_run(rows[:3]))

    def test_presets(self):
        self.assertEqual(dp.preset_of("ILLUMINA", "NovaSeq"), "sr")
        self.assertEqual(dp.preset_of("PACBIO_SMRT", "Sequel IIe"), "map-hifi")
        self.assertEqual(dp.preset_of("PACBIO_SMRT", "PacBio RS II"), "map-pb")
        self.assertEqual(dp.preset_of("OXFORD_NANOPORE", "PromethION"), "map-ont")

    def test_first_file(self):
        self.assertEqual(dp.first_file("h/x.fastq.gz;h/x_1.fastq.gz;h/x_2.fastq.gz"), "https://h/x_1.fastq.gz")
        self.assertEqual(dp.first_file("h/x.fastq.gz"), "https://h/x.fastq.gz")

    def fake_get(self, biosample_runs=(), strain_runs=(), biosample="SAMN1", project_runs=()):
        asked = []

        def get(url):
            asked.append(url)
            header = "\t".join(dp.RUN_FIELDS)

            def table(rows):
                return "\n".join([header] + ["\t".join(r.get(f, "") for f in dp.RUN_FIELDS) for r in rows])
            if "datasets" in url:
                info = {"bioproject_lineage": [{"bioprojects": [{"accession": "PRJNA1"}]}],
                        "bioproject_accession": "PRJNA1"}
                if biosample:
                    info["biosample"] = {"accession": biosample}
                return json.dumps({"reports": [{"assembly_info": info}]})
            if "filereport" in url and "SAMN1" in url:
                return table(biosample_runs)
            if "filereport" in url and "PRJNA1" in url:
                return table(project_runs)
            if "filereport" in url:  # a curated run
                return table([run_row(url.split("accession=")[1].split("&")[0])])
            if "search" in url:
                return table(strain_runs)
            return ""
        return get, asked

    def test_find_runs_in_turn(self):
        genome = {"accession": "GCA_1.1", "strain": "S288C", "taxid": "559292", "species_taxid": "4932"}
        get, _ = self.fake_get([run_row("bs")], [run_row("st")])
        self.assertEqual(dp.find_runs(genome, {"GCA_1.1": ["SRR9"]}, get)["run"]["run_accession"], "SRR9")
        found = dp.find_runs(genome, {}, get)
        self.assertEqual((found["how"], found["run"]["run_accession"]), ("biosample", "bs"))
        # a BioProject's runs of the genome's taxa only
        other, own = dict(run_row("other", bases=9e9), tax_id="9606"), dict(run_row("pr"), tax_id="4932")
        get, _ = self.fake_get([], [run_row("st")], project_runs=[other, own])
        found = dp.find_runs(genome, {}, get)
        self.assertEqual((found["how"], found["run"]["run_accession"]), ("bioproject", "pr"))
        get, asked = self.fake_get([], [run_row("st")], project_runs=[other])
        found = dp.find_runs(genome, {}, get)
        self.assertEqual((found["how"], found["run"]["run_accession"]), ("strain", "st"))
        self.assertIn("tax_eq%28559292%29", asked[-1])
        get, _ = self.fake_get([], [run_row("st")], biosample="")
        self.assertEqual(dp.find_runs(dict(genome, strain="na"), {}, get)["how"], "none")


def record(name, acc, typ, confidence="high", lsu="t", length=7000, taxid="4932", kind="rDNA",
           regions="18S:0-1800,ITS1:1800-2100,5.8S:2100-2260,ITS2:2260-2500,28S:2500-6000"):
    return {"record": name, "kind": kind, "accession": acc, "species_taxid": taxid, "type": typ,
            "confidence": confidence, "lsu_complete": lsu, "length": str(length), "regions": regions}


class ReferenceTest(unittest.TestCase):
    def test_one_unit_per_type_the_best(self):
        rows = [record("a1", "G", "1", "low", length=9000), record("a2", "G", "1", "high"),
                record("a3", "G", "1", "high", lsu="f", length=8000), record("b1", "G", "2", "low"),
                record("s", "G", "1", kind="18S", regions="18S:0-1800"), record("p", "-", "1")]
        units, source = dp.rdna_units(rows, "G", "4932")
        self.assertEqual(([r["record"] for r in units], source), (["a2", "b1"], "genome"))

    def test_the_species_published_units_without_the_genomes(self):
        units, source = dp.rdna_units([record("p", "-", "1"), record("q", "-", "1", taxid="5")], "G", "4932")
        self.assertEqual(([r["record"] for r in units], source), (["p"], "published"))
        self.assertEqual(dp.rdna_units([], "G", "4932"), ([], "-"))

    def test_single_copies(self):
        rows = [{"accession": "G", "locus": "a", "status": "Single", "confidence": "medium", "flags": "x"},
                {"accession": "G", "locus": "b", "status": "Duplicated"},
                {"accession": "G", "locus": "c", "status": "Single", "confidence": "low"},
                {"accession": "G", "locus": "d", "status": "Single", "confidence": "medium", "flags": "gap"},
                {"accession": "H", "locus": "e", "status": "Single"},
                {"accession": "G", "locus": "f", "status": "Single"}]  # format 3: no confidence
        self.assertEqual(dp.single_copies(rows, "G"), ["a", "f"])


class ReadsTest(unittest.TestCase):
    def test_the_first_reads_of_a_stream(self):
        records = [b"@r%d\n%s\n+\n%s\n" % (i, b"A" * 100, b"I" * 100) for i in range(50)]
        data = b"".join(records)
        stream = io.BytesIO(gzip.compress(data))
        out = io.BytesIO()
        self.assertEqual(dp.take_reads(stream, out, 10, 10 ** 9), (10, 1000, 1))
        self.assertEqual(out.getvalue(), b"".join(records[:10]))
        self.assertEqual(dp.take_reads(io.BytesIO(gzip.compress(data)), io.BytesIO(), 99, 450), (5, 500, 1))
        self.assertEqual(dp.take_reads(io.BytesIO(gzip.compress(data)), io.BytesIO(), 99, 10 ** 9), (50, 5000, 1))

    def test_a_small_file_thinned_over_its_whole_length(self):
        # 50 reads of 100 bases, 1,000 bases wanted: every 5th, from the first to the last fifth
        records = [b"@r%d\n%s\n+\n%s\n" % (i, b"A" * 100, b"I" * 100) for i in range(50)]
        out = io.BytesIO()
        got = dp.take_reads(io.BytesIO(gzip.compress(b"".join(records))), out, 99, 1000, spots=50)
        self.assertEqual(got, (10, 1000, 5))
        self.assertEqual(out.getvalue(), b"".join(records[::5]))
        # a file over thin_max_bases: its first reads
        self.assertEqual(dp.take_reads(io.BytesIO(gzip.compress(b"".join(records))), io.BytesIO(), 99, 1000,
                                       spots=50, thin_max_bases=4000), (10, 1000, 1))
        self.assertEqual(dp.thinning(0, 100, 1000), 1)
        self.assertEqual(dp.thinning(10 ** 6, 150, 1.5e7), 10)

    def test_stream_reads_writes_whole_or_nothing(self):
        data = gzip.compress(b"@r\nACGT\n+\nIIII\n" * 3)
        with tempfile.TemporaryDirectory() as tmp:
            dest = os.path.join(tmp, "x.fq.gz")
            self.assertEqual(dp.stream_reads("u", dest, 2, 10 ** 9, opener=lambda u: io.BytesIO(data)), (2, 8, 1))
            with gzip.open(dest) as fh:
                self.assertEqual(fh.read(), b"@r\nACGT\n+\nIIII\n" * 2)

            def broken(u):
                raise OSError("refused")
            other = os.path.join(tmp, "y.fq.gz")
            with self.assertRaises(OSError):
                dp.stream_reads("u", other, 2, 10 ** 9, attempts=2, opener=broken)
            self.assertEqual(os.listdir(tmp), ["x.fq.gz"])


class EstimateTest(unittest.TestCase):
    def test_copies_and_flags(self):
        meta = {f"busco|{i}": ("busco", [], []) for i in range(12)}
        meta.update({"rdna|a": ("rdna", [], []), "rdna|b": ("rdna", [], []), "rdna|a#18S": ("rdna18s", [], []),
                     "rdna|b#18S": ("rdna18s", [], []), "mitochondrion|m": ("mitochondrion", [], [])})
        depths = {f"busco|{i}": 10.0 + (i % 3) for i in range(12)}
        depths.update({"rdna|a": 1500.0, "rdna|b": 500.0, "rdna|a#18S": 1400.0, "rdna|b#18S": 600.0,
                       "mitochondrion|m": 550.0})
        got = dp.estimate(depths, meta)
        self.assertEqual((got["sc_depth"], got["rdna_copies_reads"], got["rdna_18s_copies"]), (11.0, 2000 / 11,
                                                                                             2000 / 11))
        self.assertEqual((got["mt_genes"], got["mt_copies_reads"], got["pt_copies_reads"]), (1, 50.0, None))
        self.assertEqual(got["flags"], "-")
        few = {k: (v / 5 if k.startswith("busco") else v) for k, v in depths.items() if k not in ("busco|0",
                                                                                             "busco|1", "busco|2")}
        self.assertEqual(dp.estimate(few, meta)["flags"], "few_loci,low_depth")
        # reads of a few regions alone (a run sorted by position): most copies near nothing
        clumped = dict(depths, **{f"busco|{i}": 0.1 for i in range(8)})
        self.assertIn("uneven", dp.estimate(clumped, meta)["flags"].split(","))
        self.assertEqual(dp.estimate({"rdna|a": 3.0}, meta)["flags"], "no_single_copy")

    def test_gc(self):
        self.assertEqual(dp.gc("ATGCnnAT"), 1 / 3)
        self.assertEqual(dp.gc("NNN"), 0.0)


MINIMAP2 = os.environ.get("EUKORE_MINIMAP2") or shutil.which("minimap2")


def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


@unittest.skipUnless(MINIMAP2, "minimap2 is not on PATH (or EUKORE_MINIMAP2)")
class EndToEndTest(unittest.TestCase):
    """A haploid genome: 25 single-copy genes of 2 kb among unique sequence, an array of 40 rDNA units of 7 kb
    (18S 1.8 kb, 5.8S, 28S 3.5 kb) with 2 kb spacers, and 30 mitochondria of 8 kb; reads of 150 bases at 15x."""

    UNITS, MITO, DEPTH = 40, 30, 15

    def setUp(self):
        rng = random.Random(7)
        self.tmp = tempfile.mkdtemp()
        pkg = self.pkg = os.path.join(self.tmp, "Eukore_test")
        for sub in ("rrna", "busco/copies", "organelle"):
            os.makedirs(os.path.join(pkg, sub))
        genes = {f"{i}at2759": rand_seq(rng, 2000) for i in range(25)}
        unit, igs, mito = rand_seq(rng, 7000), rand_seq(rng, 2000), rand_seq(rng, 8000)
        nuclear = rand_seq(rng, 5000) + "".join(g + rand_seq(rng, 5000) for g in genes.values())
        nuclear += (unit + igs) * self.UNITS + rand_seq(rng, 5000)
        genome = [nuclear] + [mito] * self.MITO
        regions = "18S:0-1800,ITS1:1800-2100,5.8S:2100-2260,ITS2:2260-2500,28S:2500-6000"
        with open(os.path.join(pkg, "rrna", "records.tsv"), "w") as fh:
            fh.write("record\tkind\ttaxid\tspecies_taxid\taccession\tsource\tregions\tlength\tlsu_complete\ttype"
                     "\tconfidence\n")
            fh.write(f"GCA_9.1_rDNA_1\trDNA\t4932\t4932\tGCA_9.1\tgenome\t{regions}\t7000\tt\t1\thigh\n")
        with gzip.open(os.path.join(pkg, "rrna", "records.fna.gz"), "wt") as fh:
            fh.write(f">GCA_9.1_rDNA_1\n{unit}\n")
        with open(os.path.join(pkg, "busco", "copies.tsv"), "w") as fh:
            fh.write("accession\tlocus\tstatus\tconfidence\tflags\n")
            for locus in genes:
                fh.write(f"GCA_9.1\t{locus}\tSingle\thigh\t-\n")
        for locus, seq in genes.items():
            with gzip.open(os.path.join(pkg, "busco", "copies", f"{locus}.fna.gz"), "wt") as fh:
                fh.write(f">GCA_9.1\n{seq}\n>GCA_8.1\n{rand_seq(rng, 2000)}\n")
        with open(os.path.join(pkg, "organelle", "genes.tsv"), "w") as fh:
            fh.write("record\tgene\torganelle\ttaxid\tspecies_taxid\tpartial\n")
            fh.write("NC_1_COX1\tCOX1\tmitochondrion\t4932\t4932\tf\n")
        with gzip.open(os.path.join(pkg, "organelle", "genes.fna.gz"), "wt") as fh:
            fh.write(f">NC_1_COX1\n{mito[2000:3600]}\n")
        size = sum(map(len, genome))
        with open(os.path.join(pkg, "genomes.tsv"), "w") as fh:
            fh.write("accession\tspecies_taxid\ttaxid\tstrain\tgenome_size\trepresentative\tlocal\n")
            fh.write(f"GCA_9.1\t4932\t4932\tna\t{len(nuclear)}\tt\t-\n")
        with open(os.path.join(pkg, "manifest.json"), "w") as fh:
            json.dump({"format": 3, "inputs": {}, "files": {}}, fh)
        # reads: uniform starts over every molecule, either strand
        work = self.work = os.path.join(self.tmp, "work")
        os.makedirs(os.path.join(work, "reads"))
        comp = str.maketrans("ACGT", "TGCA")
        n = int(self.DEPTH * size / 150)
        with gzip.open(os.path.join(work, "reads", "GCA_9.1.fq.gz"), "wt", compresslevel=1) as fh:
            weights = [len(m) for m in genome]
            for i in range(n):
                m = rng.choices(genome, weights)[0]
                s = rng.randrange(len(m) - 150)
                r = m[s:s + 150]
                if rng.random() < 0.5:
                    r = r.translate(comp)[::-1]
                fh.write(f"@r{i}\n{r}\n+\n{'I' * 150}\n")
        with open(os.path.join(work, "reads", "GCA_9.1.json"), "w") as fh:
            json.dump({"run": "SRR0", "how": "manual", "platform": "ILLUMINA", "preset": "sr", "reads": n,
                       "bases": n * 150}, fh)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_copies_measured_and_written_into_the_package(self):
        args = ["--package", self.pkg, "--work", self.work, "--minimap2", MINIMAP2, "-t", "2"]
        rdna_depth.main(["measure"] + args)
        with open(os.path.join(self.work, "depth", "GCA_9.1.json")) as fh:
            got = json.load(fh)
        self.assertEqual(got["sc_loci"], 25)
        self.assertAlmostEqual(got["sc_depth"], self.DEPTH, delta=self.DEPTH * 0.15)
        self.assertAlmostEqual(got["rdna_copies_reads"], self.UNITS, delta=self.UNITS * 0.1)
        self.assertAlmostEqual(got["rdna_18s_copies"], self.UNITS, delta=self.UNITS * 0.1)
        self.assertAlmostEqual(got["mt_copies_reads"], self.MITO, delta=self.MITO * 0.1)
        from eukore import package as package_mod  # noqa: PLC0415
        copies = os.path.join(self.pkg, package_mod.COPY_NUMBERS)  # a format-5 package's copy numbers (2026-11-08)
        package_mod.write_gz(copies, package_mod.COPY_NUMBER_COLUMNS, [dict(
            accession="GCA_9.1", species_taxid=4932, set="rrna", element="rDNA_unit", copies=1, partial=0,
            basis="assembly")])
        rdna_depth.main(["table", "--into_package"] + args)
        cn = {(r["set"], r["element"], r["basis"]): r["copies"] for r in package_mod.read_gz(copies)}
        self.assertEqual(cn[("rrna", "rDNA_unit", "assembly")], "1")  # the assembly's kept, the reads' added
        self.assertAlmostEqual(float(cn[("rrna", "rDNA_unit", "reads")]), self.UNITS, delta=self.UNITS * 0.1)
        self.assertIn(("mitochondrion", "genome", "reads"), cn)
        self.assertNotIn(("plastid", "genome", "reads"), cn)  # not measured
        rows = dp.read_dicts(os.path.join(self.pkg, "genomes.tsv"))
        self.assertEqual(list(rows[0])[-4:], ["rdna_copies_reads", "mt_copies_reads", "pt_copies_reads", "local"])
        self.assertAlmostEqual(int(rows[0]["rdna_copies_reads"]), self.UNITS, delta=self.UNITS * 0.1)
        self.assertEqual(rows[0]["pt_copies_reads"], "-")
        with open(os.path.join(self.pkg, "manifest.json")) as fh:
            manifest = json.load(fh)
        self.assertEqual(manifest["files"]["genomes.tsv"][1], util.md5(os.path.join(self.pkg, "genomes.tsv")))
        self.assertEqual(manifest["inputs"]["read_depth"]["measured"], 1)
        rdna_depth.main(["table", "--into_package"] + args)  # again: the columns replaced, not added twice
        self.assertEqual(len(dp.read_dicts(os.path.join(self.pkg, "genomes.tsv"))[0]), 10)
        self.assertEqual(len(package_mod.read_gz(copies)), len(cn))  # the reads' rows replaced, not added twice


if __name__ == "__main__":
    unittest.main()
