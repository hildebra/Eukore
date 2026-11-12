#!/usr/bin/env python3
"""Tests of tools/run_report.py (since 2026-11-11): a work folder's genomes told done, failed, killed, running and
never tried, the failures' reasons counted with their paths left out, the done genomes' checks and completeness; a
package's benchmarks turned into findings; nothing in either changed.

    python3 -m unittest discover -s tests -t tests
"""

import json
import os
import socket
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))
import run_report  # noqa: E402

PLAN = {"GCA_000001.1": "11", "GCA_000002.1": "12", "GCA_000003.1": "13", "GCA_000004.1": "14",
        "GAJY01000000": "15", "EUK_99": "16", "GCA_000007.1": "17"}


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


def tsv(path, head, rows):
    write(path, "\t".join(head) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


def done(work, acc, single, duplicated=0, missing=0, **extra):
    loci = single + duplicated + missing
    summary = dict(dict(accession=acc, searched=[f"L{i}" for i in range(loci)], genome_record=2, seconds=3600,
                        threads=8, bytes=1e8, status={"Single": single, "Duplicated": duplicated, "Missing": missing}),
                   **extra)
    write(os.path.join(work, "done", acc, "done.json"), json.dumps(summary))


class RunReport(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.work = os.path.join(self.tmp.name, "WORK")
        tsv(os.path.join(self.work, "genome_set.tsv"), ("species_taxid", "accession"),
            [(t, a) for a, t in PLAN.items()])
        done(self.work, "GCA_000001.1", 95, 3, 2)
        done(self.work, "GCA_000002.1", 20, 0, 80, rrna_failed="nhmmer: Unable to guess alphabet", genome_record=1)
        for acc, path in (("GCA_000003.1", "/scratch/a/x.fna"), ("GAJY01000000", "/scratch/b/y.fna")):
            write(os.path.join(self.work, "failed", acc + ".txt"),
                  f"{acc}\t2026-11-11 10:00:00\tq512n19\t24855336\tOSError: [Errno 28] No space left on device: "
                  f"'{path}'\n")
        stale = os.path.join(self.work, "done", ".GCA_000004.1.claim")
        write(stale, "otherhost.cluster 123\n")
        os.utime(stale, (time.time() - 3600, time.time() - 3600))
        write(os.path.join(self.work, "done", ".q512n19.a-77.GCA_000004.1.partial", "x"), "")
        write(os.path.join(self.work, "done", ".EUK_99.claim"), f"{socket.gethostname()} {os.getpid()}\n")
        self.before = self.listing()

    def tearDown(self):
        self.tmp.cleanup()

    def listing(self):
        return sorted((os.path.relpath(os.path.join(d, f), self.tmp.name), os.path.getmtime(os.path.join(d, f)))
                      for d, _s, fs in os.walk(self.tmp.name) for f in fs)

    def test_fates(self):
        facts = run_report.run_facts(self.work, threads=2)
        self.assertEqual(facts["fate"], {
            "GCA_000001.1": "done", "GCA_000002.1": "done", "GCA_000003.1": "failed", "GAJY01000000": "failed",
            "GCA_000004.1": "killed", "EUK_99": "running", "GCA_000007.1": "never tried"})
        self.assertEqual(facts["stale_claims"], 1)
        genomes_folder = os.path.join(self.tmp.name, "genomes")
        os.makedirs(genomes_folder)
        self.assertEqual(run_report.run_facts(self.work, genomes_folder, 2)["fate"]["GCA_000007.1"], "not fetched")

    def test_kinds_and_reasons(self):
        self.assertEqual([run_report.kind_of(a) for a in ("GCF_1.1", "GAJY01000000", "EUK_9", "XYZ")],
                         ["assembly", "transcriptome", "records alone", "other"])
        self.assertEqual(run_report.reason_of("OSError: no space: '/a/b'"), run_report.reason_of("OSError: no space: "
                                                                                                  "'/c/d'"))

    def test_report(self):
        rep = run_report.Report()
        run_report.report_run(rep, self.work, run_report.run_facts(self.work, threads=2))
        text = rep.text("t")
        self.assertIn("5 of 7 planned genomes not done", text)
        self.assertIn("the commonest failure (2 of 2): `OSError: [Errno 28] No space left on device: ...`", text)
        self.assertIn("1 done genomes with rrna_failed", text)
        self.assertIn("1 done genomes with older genome record", text)
        self.assertIn("| GCA_000002.1 | 0.200 |", text)  # the poorest
        self.assertEqual(self.listing(), self.before)  # read only

    def test_package(self):
        pkg = os.path.join(self.tmp.name, "Eukore_x")
        bench = os.path.join(pkg, "benchmarks")
        write(os.path.join(bench, "summary.json"), json.dumps(dict(
            species=2, busco_loci=100, left_out={"reads": "no minimap2"},
            completeness=dict(median=0.6, at_least_90=0.5, by_supergroup={"Obazoa": 0.97, "Discoba": 0.41}),
            orthodb=dict(called=0.95, concordant_of_called=0.97, median_recall=0.99, verdicts={},
                         held_out=dict(called=0.8, concordant_of_called=0.7, median_recall=0.9, verdicts={})),
            barcode_gap={"18S": {"resolved": 0.3}, "BUSCO": {"resolved": 0.9}})))
        tsv(os.path.join(bench, "composition.tsv"),
            ("supergroup", "species", "with_genome", "rrna_only", "complete_median", "complete_90", "with_rrna",
             "with_operon", "with_organelle"),
            [("Rhizaria", 100, 10, 90, "0.8", "0.5", 95, 3, 1), ("all", 100, 10, 90, "0.8", "0.5", 95, 3, 1)])
        tsv(os.path.join(bench, "flags.tsv"), ("supergroup", "confidence", "flags", "copies"),
            [("Discoba", "high", "-", 100), ("Discoba", "low", "own_species", 50)])
        tsv(os.path.join(bench, "completeness.tsv"),
            ("accession", "name", "supergroup", "source", "level", "loci", "duplicated", "complete", "missing"),
            [("GCA_1", "A b", "Obazoa", "ncbi", "Contig", 100, 40, "0.95", 5)])
        tsv(os.path.join(pkg, "species.tsv"), ("taxid", "group"), [("11", "Obazoa"), ("13", "Discoba")])
        rep = run_report.Report()
        run_report.report_run(rep, self.work, run_report.run_facts(self.work, threads=2),
                              run_report.supergroups(pkg))
        run_report.report_package(rep, pkg)
        text = rep.text("t")
        for want in ("1 benchmarks left out: reads", "Discoba 0.41", "the held-out calls concordant 0.7",
                     "the held-out calls called 0.8", "18S 0.30", "Rhizaria (10 of 100 with a genome)",
                     "Discoba: 33.3% of its BUSCO copies low or unverified", "| own_species | 50 |",
                     "1 genomes with 20% or more of their BUSCOs duplicated", "| Discoba | 1 | 1 | 100.0% |"):
            self.assertIn(want, text)
        self.assertNotIn("the build's calls concordant", text)
        findings = text.split("## Findings")[1].split("\n## ")[0]
        self.assertLess(findings.index("**high**"), findings.index("**medium**"))
        self.assertLess(findings.index("**medium**"), findings.index("**low**"))

    def test_no_package_benchmarks(self):
        rep = run_report.Report()
        run_report.report_package(rep, os.path.join(self.tmp.name, "none"))
        self.assertIn("no benchmarks in", rep.text("t"))

    def test_main(self):
        out = os.path.join(self.tmp.name, "report.md")
        self.assertEqual(run_report.main(["--work", self.work, "--out", out, "-t", "2"]), 0)
        with open(out) as fh:
            self.assertTrue(fh.read().startswith("# Eukore run report"))


if __name__ == "__main__":
    unittest.main()
