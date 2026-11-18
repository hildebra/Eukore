#!/usr/bin/env python3
"""Tests of the benchmark step's evaluation bundle (benchmarks/eukore_bench/bundle.py; since 2026-11-18): the results'
tables a step wrote, not the package's own copied ones, the records that attract background reads with their rows of
the package, the logs; a table too large cut short.

    python3 -m unittest discover -s tests -t tests
"""

import json
import os
import shutil
import sys
import tarfile
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore_bench import bundle  # noqa: E402


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


class BundleTest(unittest.TestCase):
    def test_a_benchmark_run_bundled(self):
        with tempfile.TemporaryDirectory() as tmp:
            package, results, bench = (os.path.join(tmp, n) for n in ("Eukore_2026.11", "Eukore_2026.11_benchmarks",
                                                                      "bench"))
            write(os.path.join(package, "benchmarks", "completeness.tsv"), "accession\tcomplete\nA\t0.9\n")
            write(os.path.join(package, "benchmarks", "summary.json"), "{}")
            write(os.path.join(package, "genomes.tsv"),
                  "accession\ttier\tflags\nGCA_1.1\tgenome\t-\nGCA_2.1\tcontaminated\tFCS-GX\nGCA_3.1\tgenome\t-\n")
            write(os.path.join(package, "rrna", "records.tsv"),
                  "record\tkind\taccession\tspecies_taxid\nr1\t18S\tGCA_2.1\t7\nr2\t18S\tGCA_3.1\t8\nr3\trDNA\t-\t9\n")
            # the results: the package's own table copied (copy2: same size and time), the tables the steps wrote
            os.makedirs(results)
            shutil.copy2(os.path.join(package, "benchmarks", "completeness.tsv"),
                         os.path.join(results, "completeness.tsv"))
            write(os.path.join(results, "summary.json"), json.dumps({"reads": {}}))
            write(os.path.join(results, "reads.tsv"), "test\tn\nY\t3\n")
            write(os.path.join(results, "big.tsv"), "a\tb\n" + "1\t2\n" * 500)
            write(os.path.join(results, ".host-1.reads.tsv.partial"), "half")
            write(os.path.join(results, "background_hits.tsv"),
                  "test\tsource\treference\trecord\ttaxid\tpart\treads\tidentity\n"
                  "Bug\tbackground\tEukore\tr1\t7\trrna:18S\t30\t95.0\n"
                  "Bug\tbackground\tEukore\tGCA_1.1|88at2759\t7\tbusco:88at2759\t5\t99.0\n"
                  "Bug\tbackground\tPR2\tx9\t7\t-\t2\t92.0\n")
            write(os.path.join(results, "figures", "Fig1.png"), "png")
            write(os.path.join(results, "figures", "notes.txt"), "not a figure")
            write(os.path.join(bench, "compare", "references.tsv"), "name\tfasta\ttable\nOther\ta\tb\n")
            write(os.path.join(bench, "read_genomes", "tests.tsv"), "accession\tspecies_taxid\nGCA_1.1\t7\n")
            log = os.path.join(tmp, "benchmark.1.etxt")
            write(log, "x" * 50 + "the end\n")
            logged = []
            with mock.patch.object(bundle, "BUNDLE_FILE_BYTES", 100), mock.patch.object(bundle, "BUNDLE_LINES", 10), \
                    mock.patch.object(bundle, "LOG_BYTES", 20):
                path = bundle.make(package, results, bench, logs=[log, os.path.join(tmp, "missing.log")],
                                   say=logged.append)
            self.assertEqual(path, os.path.join(tmp, "Eukore_2026.11_benchmarks_evaluation.tar.gz"))
            self.assertFalse([n for n in os.listdir(tmp) if n.startswith(".") or n.startswith("bundle.")])  # no leftovers
            top = "Eukore_2026.11_benchmarks_evaluation/"
            with tarfile.open(path) as tar:
                names = set(tar.getnames())
                text = lambda n: tar.extractfile(top + n).read().decode()  # noqa: E731
                self.assertTrue({top + n for n in (
                    "BUNDLE.md", "results/summary.json", "results/reads.tsv", "results/big.tsv",
                    "results/background_hits.tsv", "results/figures/Fig1.png", "bench/references.tsv",
                    "bench/tests.tsv", "attracting/records.tsv", "attracting/genomes.tsv",
                    "logs/benchmark.1.etxt")} <= names)
                self.assertNotIn(top + "results/completeness.tsv", names)  # the package's own, copied
                self.assertNotIn(top + "results/figures/notes.txt", names)
                self.assertFalse([n for n in names if "partial" in n])
                self.assertEqual(len(text("results/big.tsv").splitlines()), 10)  # cut to its first lines
                self.assertEqual(text("logs/benchmark.1.etxt"), "xxxxxxxxxxxx" + "the end\n"[-8:])
                # the package's rows of what attracts the background: the rRNA record r1 (30 reads) and its genome, the
                # BUSCO copy's genome; not r2 (not named), nor the record of another reference (x9: not in the package)
                records = [l.split("\t") for l in text("attracting/records.tsv").splitlines()]
                self.assertEqual([(r[0], r[-1]) for r in records], [("record", "reads_in_background"), ("r1", "30")])
                self.assertEqual([l.split("\t")[0] for l in text("attracting/genomes.tsv").splitlines()],
                                 ["accession", "GCA_1.1", "GCA_2.1"])
                md = text("BUNDLE.md")
                self.assertIn("completeness.tsv", md)  # said to be left out
                self.assertIn("big.tsv", md)
                self.assertIn("Other\ta\tb", md)
            self.assertTrue(any("benchmark bundle" in t for t in logged))

    def test_a_folder_that_cannot_be_written_is_said_not_raised(self):
        with tempfile.TemporaryDirectory() as tmp:
            results = os.path.join(tmp, "Eukore_2026.11_benchmarks")
            os.makedirs(results)
            logged = []
            self.assertIsNone(bundle.make(tmp, results, None, out=os.path.join(tmp, "no", "such", "x.tar.gz"),
                                          say=logged.append))
            self.assertTrue(any(t.startswith("WARNING: no benchmark bundle") for t in logged))


if __name__ == "__main__":
    unittest.main()
