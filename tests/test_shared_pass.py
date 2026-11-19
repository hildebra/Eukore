#!/usr/bin/env python3
"""Tests of the final job's passes shared with the helper nodes (since 2026-11-19: qc.write_pass, qc.help_pass,
qc.shared_results): the organelle genes of the assemblies and the gene models' comparison, run by the helpers too,
each task under its item's claim; the final job awaits the items other nodes hold.

    python3 -m unittest discover -s tests -t tests
"""

import concurrent.futures
import contextlib
import os
import sys
import tempfile
import threading
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import eukore_build  # noqa: E402
from eukore import genomes, qc  # noqa: E402


def task(folder, item):
    """A pass' task as the passes' own: its item made under its claim, ELSEWHERE where another holds it or it is
    made."""
    with genomes.claimed_file(os.path.join(folder, f".{item}.claim")) as held:
        out = os.path.join(folder, f"{item}.done")
        if not held or os.path.isfile(out):
            return item, qc.ELSEWHERE, None
        if item == "bad":
            return item, None, "it failed"
        with open(out, "w") as fh:
            fh.write(f"{os.getpid()}\n")
        return item, 1, None


@contextlib.contextmanager
def thread_pool(workers):
    """util.process_pool's stand-in: this test module's task found in a thread as in a worker."""
    with concurrent.futures.ThreadPoolExecutor(max(1, workers)) as pool:
        pool.size = max(1, workers)
        yield pool


class SharedResultsTest(unittest.TestCase):
    def test_an_item_another_node_holds_is_awaited(self):
        with tempfile.TemporaryDirectory() as tmp:
            held = genomes.claim_file(os.path.join(tmp, ".b.claim"))  # a helper's

            def helper():  # its task ends: made, its claim let go
                with open(os.path.join(tmp, "b.done"), "w") as fh:
                    fh.write("helper\n")
                os.remove(held)
            threading.Timer(0.3, helper).start()
            with thread_pool(2) as pool:
                got = sorted(qc.shared_results(pool, lambda i: task(tmp, i), ["a", "b"], lambda i: i,
                                               lambda i: os.path.isfile(os.path.join(tmp, f"{i}.done")), poll=0.1))
            self.assertEqual(got, [("a", 1, None), ("b", qc.ELSEWHERE, None)])

    def test_a_dead_helpers_item_is_made_here(self):
        with tempfile.TemporaryDirectory() as tmp:
            held = genomes.claim_file(os.path.join(tmp, ".b.claim"))
            threading.Timer(0.3, os.remove, (held,)).start()  # the helper gone, its item not made
            with thread_pool(2) as pool:
                got = sorted(qc.shared_results(pool, lambda i: task(tmp, i), ["a", "b"], lambda i: i,
                                               lambda i: os.path.isfile(os.path.join(tmp, f"{i}.done")), poll=0.1))
            self.assertEqual(got, [("a", 1, None), ("b", 1, None)])

    def test_without_wait_another_nodes_item_is_left(self):
        # a SLURM array's tasks: another task's genome is left to it, as before
        with tempfile.TemporaryDirectory() as tmp:
            genomes.claim_file(os.path.join(tmp, ".b.claim"))
            with thread_pool(1) as pool:
                got = sorted(qc.shared_results(pool, lambda i: task(tmp, i), ["a", "b"], lambda i: i,
                                               lambda _i: False, wait=False))
            self.assertEqual(got, [("a", 1, None), ("b", qc.ELSEWHERE, None)])


class HelpPassTest(unittest.TestCase):
    def test_a_helper_runs_a_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            items = os.path.join(tmp, "items")
            os.makedirs(items)
            genomes.claim_file(os.path.join(items, ".c.claim"))  # the final job's own
            path = qc.write_pass(os.path.join(tmp, "qc", "organelle_pass"), "test_shared_pass:task", (items,),
                                 ["a", "bad", "c", "d"])
            said = []
            with mock.patch.object(qc.util, "process_pool", thread_pool):
                self.assertEqual(qc.help_pass(path, 2, said.append), (2, True))
            self.assertEqual(sorted(n for n in os.listdir(items) if n.endswith(".done")), ["a.done", "d.done"])
            self.assertIn("organelle_pass: bad failed: it failed", said)
            self.assertTrue(any("2 tasks made here, 1 failed" in s for s in said))
            # taken back: nothing run
            qc.end_pass(os.path.dirname(path))
            with mock.patch.object(qc.util, "process_pool", thread_pool):
                self.assertEqual(qc.help_pass(path, 2, said.append), (0, False))

    def test_a_pass_gives_way_to_searches(self):
        with tempfile.TemporaryDirectory() as tmp:
            items = os.path.join(tmp, "items")
            os.makedirs(items)
            path = qc.write_pass(os.path.join(tmp, "qc", "p"), "test_shared_pass:task", (items,), ["a", "b", "c"])
            with mock.patch.object(qc.util, "process_pool", thread_pool):
                made, tried = qc.help_pass(path, 1, lambda _t: None, busy=lambda: True)
            self.assertEqual((made, tried), (0, False))  # nothing begun while a search waits

    def test_the_helper_loop_takes_the_passes(self):
        # help_search: no search left, the pass' tasks run; each version of a pass gone through once
        with tempfile.TemporaryDirectory() as tmp:
            items = os.path.join(tmp, "items")
            os.makedirs(items)
            root = os.path.join(tmp, "qc")
            qc.write_pass(os.path.join(root, "annotation_pass"), "test_shared_pass:task", (items,), ["a", "b"])
            said = []
            with mock.patch.object(qc.util, "process_pool", thread_pool):
                qc.help_search(root, "vsearch", 2, said.append)
            self.assertEqual(sorted(os.listdir(items)), ["a.done", "b.done"])
            self.assertIn("0 loci searched here, 2 tasks of the passes made; none left", said)


class SharedPassTest(unittest.TestCase):
    def test_the_build_shares_a_pass_with_its_helpers(self):
        with tempfile.TemporaryDirectory() as tmp:
            opts = types.SimpleNamespace(work=tmp, qc_helpers="echo 4713")
            self.assertIsNone(eukore_build.shared_pass(types.SimpleNamespace(work=tmp, qc_helpers=None), "x"))
            with mock.patch.object(eukore_build, "say"), mock.patch.object(eukore_build.atexit, "register"):
                share = eukore_build.shared_pass(opts, "organelle_pass")
                end = share("eukore.organelle:run_claimed", (), [1, 2])
            folder = os.path.join(tmp, "qc", "organelle_pass")
            self.assertTrue(os.path.isfile(os.path.join(folder, qc.PASS)))
            end()
            self.assertFalse(os.path.exists(os.path.join(folder, qc.PASS)))
            # the pass' files where the helpers read them: under the WORK, not --local_tmp
            with eukore_build.pass_files(opts, share) as files:
                self.assertEqual(os.path.dirname(files), os.path.join(tmp, "qc"))


if __name__ == "__main__":
    unittest.main()
