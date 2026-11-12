#!/usr/bin/env python3
"""Tests of the package stage's reports (since 2026-11-11: a step of hours, the copies' vsearch checks, said nothing):
each step said when it begins and ends and while it runs (eukore_build.timed), and qc.check's loci searched reported
by a progress (util.Progress), its checks the same with or without one.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import eukore_build  # noqa: E402
from eukore import ncbi as ncbi_mod, qc, util  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST, PARADOXUS, LACTIS = 4932, 27291, 28985


class TimedTest(unittest.TestCase):
    def test_a_step_is_said_begun_running_and_ended(self):
        said, stages = [], []
        usage = types.SimpleNamespace(line=lambda: "; 3 cores busy")
        opts = types.SimpleNamespace(report_minutes=1 / 60, usage=usage)  # a report every second
        with mock.patch.object(eukore_build, "say", said.append):
            with eukore_build.timed(stages, "BUSCO copies and their checks", opts):
                time.sleep(1.5)
        self.assertEqual(said[0], "package: BUSCO copies and their checks ..")
        self.assertEqual(said[1], "  BUSCO copies and their checks: 1s so far; 3 cores busy")
        self.assertRegex(said[-1], r"^package: BUSCO copies and their checks: took 1s; 3 cores busy$")
        self.assertEqual([s for s, _t in stages], ["BUSCO copies and their checks"])
        # without opts (a caller of before): said begun and ended, no reports while it runs
        said.clear()
        with mock.patch.object(eukore_build, "say", said.append):
            with eukore_build.timed(stages, "FCS-GX"):
                pass
        self.assertEqual(said, ["package: FCS-GX ..", "package: FCS-GX: took 0s"])
        # a step that fails: its reports stop, nothing more said
        said.clear()
        with mock.patch.object(eukore_build, "say", said.append), self.assertRaises(ValueError):
            with eukore_build.timed(stages, "organelle genes", opts):
                raise ValueError("x")
        self.assertEqual(said, ["package: organelle genes .."])


class CheckProgressTest(unittest.TestCase):
    def test_the_loci_searched_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            vsearch = fx.vsearch(tmp)
            gene = fx.cds_of(fx.protein(3))
            other = fx.cds_of(fx.protein(4))
            copies = {("A", "L1"): gene, ("B", "L1"): gene, ("C", "L1"): gene,
                      ("A", "L2"): other, ("B", "L2"): other,
                      ("A", "L3"): gene}  # one species alone: not searched
            species = {"A": YEAST, "B": PARADOXUS, "C": LACTIS}
            said, made = [], []

            def progress(what, weights):
                made.append(weights)
                return util.Progress(what, weights, said.append, every=0, unit="")
            got = qc.check(copies, species, ncbi, vsearch, 1, progress=progress)
            self.assertEqual(got, qc.check(copies, species, ncbi, vsearch, 1))
            self.assertEqual(made, [[9, 4]])  # each locus by its copies squared
            self.assertEqual(len(said), 2)
            self.assertRegex(said[-1], r"^  2/2 loci searched \(vsearch\) \(0 failed\), 0 left; took \d+s$")
            # no vsearch: no progress made
            made.clear()
            qc.check(copies, species, ncbi, None, 1, progress=progress)
            self.assertEqual(made, [])


class SearchesTest(unittest.TestCase):
    """qc.check's searches kept in a folder and shared with helper nodes (qc.searched_loci; since 2026-11-11)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = self.tmp.name
        self.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        self.vsearch = fx.vsearch(tmp)
        self.root = os.path.join(tmp, "WORK", "qc")
        self.folder = os.path.join(self.root, "busco")
        genes = [fx.cds_of(fx.protein(s)) for s in (3, 4, 5)]
        self.copies = {(a, f"L{i}"): g for i, g in enumerate(genes) for a in "ABC"}
        self.copies[("C", "L2")] = genes[2][:-30] + "A" * 30  # a copy apart
        self.species = {"A": YEAST, "B": PARADOXUS, "C": LACTIS}
        self.alone = qc.check(self.copies, self.species, self.ncbi, self.vsearch, 1)

    def tearDown(self):
        self.tmp.cleanup()

    def check(self, helpers=None, vsearch=None, copies=None):
        return qc.check(copies or self.copies, self.species, self.ncbi, vsearch or self.vsearch, 1,
                        searches={"folder": self.folder, "helpers": helpers})

    def hits(self):
        return sorted(n for n in os.listdir(self.folder) if n.endswith(".hits.gz"))

    def test_kept_and_read_again(self):
        self.assertEqual(self.check(), self.alone)
        kept = self.hits()
        self.assertEqual(len(kept), 3)
        self.assertEqual(sorted(os.listdir(self.folder)), kept)  # no inputs, claims or listing left
        # read again, vsearch not run (a stand-in that fails)
        broken = fx.tool(self.tmp.name, "broken", "#!/bin/sh\nexit 3\n")
        self.assertEqual(self.check(vsearch=broken), self.alone)
        # a locus whose copies changed is searched again, its old search removed
        changed = {**self.copies, ("B", "L1"): self.copies[("A", "L2")]}
        self.assertEqual(self.check(copies=changed), qc.check(changed, self.species, self.ncbi, self.vsearch, 1))
        now = self.hits()
        self.assertEqual(len(now), 3)
        self.assertEqual(len(set(now) & set(kept)), 2)

    def test_helpers_search_the_loci(self):
        searched, stopped = [], []

        def helpers(n):  # another node, through before this one claims any locus
            self.assertEqual(n, 3)
            listing = os.path.join(self.folder, "todo.tsv")
            with open(listing) as fh:
                weights = [int(line.split("\t")[1]) for line in fh]
            self.assertEqual(weights, sorted(weights, reverse=True))  # the largest first
            searched.append(qc.help_search(self.root, self.vsearch, 1, log=lambda _t: None))
            return lambda: stopped.append(True)
        self.assertEqual(self.check(helpers), self.alone)
        self.assertEqual((searched, stopped), ([3], [True]))
        self.assertEqual(len(self.hits()), 3)
        # nothing left: a helper started late ends at once
        self.assertEqual(qc.help_search(self.root, self.vsearch, 1, log=lambda _t: None), 0)
        # helpers not started (None): this node searches them all
        shutil.rmtree(self.folder)
        self.assertEqual(self.check(lambda n: None), self.alone)
        # the helper's command, another process: its claims its own
        shutil.rmtree(self.folder)
        runs = []

        def command(_n):
            runs.append(subprocess.run([sys.executable, "-m", "eukore.qc", "helper", self.root, "--vsearch",
                                        self.vsearch, "-t", "1"], cwd=os.path.dirname(HERE), capture_output=True,
                                       text=True, check=False))
            return None
        self.assertEqual(self.check(command), self.alone)
        self.assertEqual(runs[0].returncode, 0, runs[0].stderr)
        self.assertIn("3 loci searched here; none left", runs[0].stdout)

    def test_another_nodes_locus_awaited_and_a_dead_ones_taken_over(self):
        todo = [(f"L{i}", {a: self.copies[(a, f"L{i}")] for a in "ABC"}, None) for i in range(3)]
        os.makedirs(self.folder)
        bases = [qc.search_base(l, q) for l, q, _d in todo]
        live, dead = (qc.search_files(self.folder, b)[3] for b in bases[:2])
        for path in (live, dead):
            with open(path, "w") as fh:
                fh.write("othernode\t123\n")
        os.utime(dead, (time.time() - 3600,) * 2)  # its node dead: taken over

        def other_node():  # the live claim's node ends its search
            time.sleep(0.3)
            q, _d, hits, claim = qc.search_files(self.folder, bases[0])
            qc.write_inputs(self.folder, bases[0], todo[0][1])
            out = hits[:-3]
            util.vsearch_files(self.vsearch, q, q, out, qc.MIN_ID, 1, *qc.SEARCH_LIMITS, qc.SEARCH_OPTIONS)
            with open(out, "rb") as src, gzip.open(hits, "wb") as dst:
                dst.write(src.read())
            os.remove(claim)
        t = threading.Thread(target=other_node)
        t.start()
        said = []
        got = qc.searched_loci(todo, self.vsearch, 1, self.folder, poll=0.05,
                               progress=lambda w, ws: util.Progress(w, ws, said.append, every=0, unit=""))
        t.join()
        self.assertEqual(sorted(got), ["L0", "L1", "L2"])
        self.assertTrue(all(got.values()))
        self.assertIn("searched by other nodes", "".join(said))
        self.assertFalse(any(n.endswith(".claim") for n in os.listdir(self.folder)))

    def test_the_build_starts_and_stops_its_helpers(self):
        said = []
        opts = types.SimpleNamespace(work=os.path.join(self.tmp.name, "WORK"), qc_helpers="false")
        with mock.patch.object(eukore_build, "say", said.append):
            searches = eukore_build.qc_searches(opts, "busco")
            self.assertEqual(searches["folder"], self.folder)
            self.assertIsNone(searches["helpers"](3))  # refused: this node alone
            self.assertIn("WARNING: the helpers were not started (exit 1", said[-1])
            opts.qc_helpers = "echo 4711"
            with mock.patch.object(eukore_build.shutil, "which", return_value="/bin/scancel"), \
                    mock.patch.object(eukore_build.subprocess, "run", wraps=subprocess.run) as run:
                stop = eukore_build.qc_searches(opts, "busco")["helpers"](3)
                self.assertIn("helpers: job 4711", said[-1])
                with mock.patch.object(eukore_build.subprocess, "run") as cancel:
                    stop()
                cancel.assert_called_once_with(["scancel", "4711"], capture_output=True, check=False)
                self.assertEqual(run.call_args_list[0].args[0], ["echo", "4711"])
        self.assertIsNone(eukore_build.qc_searches(types.SimpleNamespace(work="W", qc_helpers=None), "x")["helpers"])


if __name__ == "__main__":
    unittest.main()
