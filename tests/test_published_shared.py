#!/usr/bin/env python3
"""Tests of the package stage's single-node steps spread out (since 2026-11-12): the published reads' annotation in
fixed batches every run on the WORK shares by claims (eukore_build.annotate_published: an array's tasks, then the
final job), its done batches a line each in one file a process; the 18S types of each species by a pool of
processes (rrna.species_types), the pairs the q-gram lemma rules out not aligned (rrna.may_reach), the types the same.

    python3 -m unittest discover -s tests -t tests
"""

import collections
import os
import random
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
from eukore import rrna  # noqa: E402
import fixtures as fx  # noqa: E402


def mutate(seq, rate, rng):
    out = []
    for c in seq:
        r = rng.random()
        if r < rate / 3:
            continue  # a deletion
        if r < 2 * rate / 3:
            out.append(rng.choice("ACGT"))  # an insertion
        out.append(rng.choice("ACGT".replace(c, "")) if rate / 3 * 2 <= r < rate else c)
    return "".join(out)


def types_before(exons):
    """types_of as it was until 2026-11-12: every pair aligned."""
    reps, out = [], {}
    for seq, _n in sorted(collections.Counter(exons).items(), key=lambda kv: (-kv[1], -len(kv[0]), kv[0])):
        for i, rep in enumerate(reps):
            if rrna.identity(seq, rep) >= rrna.SAME_TYPE:
                out[seq] = i + 1
                break
        else:
            reps.append(seq)
            out[seq] = len(reps)
    return out


class TypesTest(unittest.TestCase):
    def setUp(self):
        rng = random.Random(7)
        bases = ["".join(rng.choice("ACGT") for _ in range(rng.randint(300, 600))) for _ in range(4)]
        # around the 99% line: copies 0-3% apart, some cut short (identity over the shorter)
        self.seqs = []
        for b in bases:
            for _ in range(25):
                s = mutate(b, rng.choice([0, 0.003, 0.008, 0.01, 0.012, 0.02, 0.03]), rng)
                self.seqs.append(s[rng.randint(0, 30):] if rng.random() < 0.2 else s)

    def test_the_filter_never_rules_out_a_pair_that_reaches_it(self):
        got = collections.Counter()
        for a in self.seqs[:40]:
            for b in self.seqs[:40]:
                same = rrna.identity(a, b) >= rrna.SAME_TYPE
                may = rrna.may_reach(a, rrna.grams(a), b, rrna.grams(b), rrna.SAME_TYPE)
                self.assertTrue(may or not same, (a, b))
                got[(same, may)] += 1
        self.assertGreater(got[(False, False)], got[(False, True)])  # it rules out most pairs of other types
        self.assertEqual([rrna.max_edits(m, 0.99) for m in (50, 99, 100, 101, 1800)], [0, 0, 1, 1, 18])

    def test_the_types_are_the_same(self):
        self.assertEqual(rrna.types_of(self.seqs), types_before(self.seqs))
        a, b, c = "ACGT" * 30, "ACGT" * 29 + "ACGA", "TTTT" * 30
        self.assertEqual(rrna.types_of([a, a, b, c]), {a: 1, b: 1, c: 2})

    def test_the_first_filter_against_every_type_at_once(self):
        # since 2026-11-12 (later): the 4-gram lemma against every type at once (short_grams) rules a pair out before
        # may_reach and identity; an assembly's error-rich copies (each its own type), Ns, lower case and fragments
        # typed as every pair aligned
        rng = random.Random(11)
        base = "".join(rng.choice("ACGT") for _ in range(500))
        rich = [mutate(base, 0.02, rng) for _ in range(60)]
        odd = [s.lower() if i % 7 == 0 else s[:200] + "NN" + s[202:] if i % 5 == 0 else s[50:250] if i % 4 == 0
               else s for i, s in enumerate(self.seqs)] + ["", "ACG"]
        for seqs in (rich, odd, rich + self.seqs):
            self.assertEqual(rrna.types_of(seqs), types_before(seqs))
        self.assertEqual(len(set(rrna.types_of(rich).values())), 60)
        counts, other = rrna.short_grams("ACGTnACGTA")
        self.assertEqual((int(counts.sum()), other), (3, 4))  # ACGT, CGTA... the four holding n apart
        self.assertEqual(int(counts[0b00011011]), 2)  # ACGT twice
        self.assertEqual((int(rrna.short_grams("ACG")[0].sum()), rrna.short_grams("ACG")[1]), (0, 0))

    def test_species_typed_by_a_pool(self):
        parts = {1: self.seqs[:50], 2: self.seqs[50:], 3: [self.seqs[0]], 4: [self.seqs[1]] * 3}
        one = rrna.species_types(parts, 1)
        self.assertEqual(one, rrna.species_types(parts, 3))
        self.assertEqual(one[1], types_before(self.seqs[:50]))
        self.assertEqual(one[4], {self.seqs[1]: 1})


BARRNAP = """#!/bin/sh
echo "$$" >> "$(dirname "$0")/calls.txt"
echo '##gff-version 3'
"""


class SharedBatchesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.work = os.path.join(self.tmp.name, "WORK")
        self.barrnap = fx.tool(self.tmp.name, "barrnap", BARRNAP)
        rng = random.Random(3)
        self.reads = [(f"EUK_{i}", "".join(rng.choice("ACGT") for _ in range(200)), 4932, "EUKARYOME", "species",
                       "-") for i in range(5)]
        self.opts = types.SimpleNamespace(work=self.work, threads=1, report_minutes=10,
                                          usage=types.SimpleNamespace(line=lambda: ""))
        self.patches = [mock.patch.object(eukore_build, "PUBLISHED_BATCH", 2),
                        mock.patch.object(eukore_build, "PUBLISHED_POLL", 0.05),
                        mock.patch.object(eukore_build.rrna, "itsx_groups", return_value="all"),
                        mock.patch.object(eukore_build.shutil, "which",
                                          side_effect=lambda t: self.barrnap if t == "barrnap" else None),
                        mock.patch.object(eukore_build, "say", lambda _t: None)]
        for p in self.patches:
            p.start()
        self.folder = os.path.join(self.work, "published")
        self.batches = eukore_build.published_batches(self.reads, None)

    def tearDown(self):
        for p in self.patches:
            p.stop()
        self.tmp.cleanup()

    def calls(self):
        path = os.path.join(self.tmp.name, "calls.txt")
        if not os.path.exists(path):
            return 0
        with open(path) as fh:
            return len(fh.read().split())

    def claim(self, key, age=0):
        os.makedirs(self.folder, exist_ok=True)
        path = os.path.join(self.folder, f".{key}.claim")
        with open(path, "w") as fh:
            fh.write("othernode\t123\n")
        os.utime(path, (time.time() - age,) * 2)
        return path

    def test_batches_fixed_for_every_run(self):
        self.assertEqual([len(b) for _k, b in self.batches], [2, 2, 1])
        self.assertEqual(self.batches, eukore_build.published_batches(list(reversed(self.reads)), None))

    def test_an_array_task_leaves_the_batches_others_hold_and_the_final_job_awaits_them(self):
        first = self.batches[0][0]
        other = self.claim(first)  # another live run's
        known, unread = eukore_build.annotate_published(self.opts, self.reads, None, wait=False)
        self.assertEqual((len(known), unread, self.calls()), (3, set(), 2))
        self.assertTrue(os.path.exists(other))

        def other_run():  # its batch annotated, its lines and its done mark written, its claim removed
            time.sleep(0.3)
            with open(os.path.join(self.folder, "regions.3.othernode-123.tsv"), "w") as fh:
                fh.writelines(f"{n}\t-\tf\t+\t{eukore_build.read_digest(s)}\n" for n, s, _t in self.batches[0][1])
            with open(os.path.join(self.folder, "done.othernode-123.txt"), "w") as fh:
                fh.write(first + "\n")
            os.remove(other)
        t = threading.Thread(target=other_run)
        t.start()
        known, unread = eukore_build.annotate_published(self.opts, self.reads, None, wait=True)
        t.join()
        self.assertEqual((len(known), unread, self.calls()), (5, set(), 2))  # none annotated twice
        names = sorted(os.listdir(self.folder))
        self.assertFalse([n for n in names if n.endswith(".claim")])
        self.assertEqual(len(names), 4)  # a regions file and a done file a process: no file a batch
        # all known: a rerun annotates nothing
        self.assertEqual(eukore_build.annotate_published(self.opts, self.reads, None)[0].keys(), known.keys())
        self.assertEqual(self.calls(), 2)

    def test_a_dead_runs_batch_taken_over(self):
        self.claim(self.batches[1][0], age=3600)
        known, unread = eukore_build.annotate_published(self.opts, self.reads, None, wait=True)
        self.assertEqual((len(known), unread, self.calls()), (5, set(), 3))
        self.assertEqual(eukore_build.published_done(self.folder), {k for k, _b in self.batches})

    def test_a_failed_batch_left_out(self):
        fx.tool(self.tmp.name, "barrnap", "#!/bin/sh\nexit 3\n")
        known, unread = eukore_build.annotate_published(self.opts, self.reads, None, wait=True)
        self.assertEqual((len(known), len(unread)), (0, 5))  # every batch failed, none awaited for ever


if __name__ == "__main__":
    unittest.main()
