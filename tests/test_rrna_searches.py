#!/usr/bin/env python3
"""Tests of the rRNA checks' searches (since 2026-11-13; screen.Searcher): a large search kept in WORK/qc/rrna_<phase>/
in parts of its queries against one shared database (qc.searched_loci with its own settings, qc.Params), searched by
the helper nodes too, its inputs on the node's disk where no helper reads them; the follow-up searches several at
once; the own-species check's and the screen's results the same as with the searches of before.

    python3 -m unittest discover -s tests -t tests
"""

import collections
import hashlib
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import qc, rrna, screen, util  # noqa: E402
import fixtures as fx  # noqa: E402
import test_foreign as tf  # noqa: E402


def normal(hits):
    return {q: sorted(v) for q, v in hits.items() if v}


class KeysTest(unittest.TestCase):
    def test_the_loci_searches_keep_their_names(self):
        # the BUSCO and organelle loci's kept searches (DEFAULT) are named as before 2026-11-13: read, not run again
        q, db = {"a": "ACGT", "b": "ACGA"}, {"c": "TTTT"}
        for d in (None, db):
            h = hashlib.sha1(f"{qc.SEARCH_RECORD}\t{qc.MIN_ID}\t{qc.SEARCH_LIMITS}\t{qc.SEARCH_OPTIONS}\n".encode())
            for part in (q,) + (() if d is None else (d,)):
                h.update(b">\n")
                for name in sorted(part):
                    h.update(f"{name}\t{part[name]}\n".encode())
            self.assertEqual(qc.search_key(q, d), h.hexdigest()[:20])
        other = qc.Params(0.8, 32, 256, ())
        self.assertNotEqual(qc.search_key(q, db, other), qc.search_key(q, db))
        self.assertEqual(qc.search_key(q, qc.Shared(db), other), qc.search_key(q, qc.Shared(dict(db)), other))

    def test_todo_lines(self):
        self.assertEqual(qc.todo_line("L1.x", 9), "L1.x\t9\n")  # the loci's: as before
        self.assertEqual(qc.todo_entry("L1.x\t9\n"), ("L1.x", qc.DEFAULT, None))
        line = qc.todo_line("own-1.y", 16, qc.Params(0.8, 32, 256, ()), "db.abc.fna")
        self.assertEqual(line, "own-1.y\t16\t0.8\t32\t256\t-\tdb.abc.fna\n")
        self.assertEqual(qc.todo_entry(line), ("own-1.y", qc.Params(0.8, 32, 256, ()), "db.abc.fna"))
        line = qc.todo_line("L2.z", 4, shared="db.d.fna")  # a shared database with the loci's settings
        self.assertEqual(qc.todo_entry(line), ("L2.z", qc.DEFAULT, "db.d.fna"))


class SearcherTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = self.tmp.name
        self.vsearch = fx.vsearch(tmp)
        self.root = os.path.join(tmp, "WORK", "qc")
        self.scratch = os.path.join(tmp, "local")
        os.makedirs(self.scratch)
        rng = random.Random(5)
        bases = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(2)]
        self.seqs = {f"s{i:02d}": tf.mutate(bases[i % 2], 0.05, rng) for i in range(11)}
        self.db = dict(self.seqs, extra="".join(rng.choice("ACGT") for _ in range(300)))
        self.plain = screen.search(self.seqs, self.db, self.vsearch, 1, tmp, 0.75, 32, 64)
        self.patch = mock.patch.multiple(screen, SEARCH_MIN_PART=3, HELPED=2)  # 11 queries: 3 parts for one node
        self.patch.start()

    def tearDown(self):
        self.patch.stop()
        self.tmp.cleanup()

    def searcher(self, helpers=None, vsearch=None, log=None, nodes=1):
        return screen.Searcher(vsearch or self.vsearch, 1, self.scratch, folder=os.path.join(self.root, "rrna"),
                               helpers=helpers, log=log, nodes=nodes)

    def test_parts_interleaved_and_sized_by_the_nodes(self):
        names = [f"n{i:03d}" for i in range(100)]
        with mock.patch.object(screen, "SEARCH_MIN_PART", 10):
            parts = screen.search_parts(names, 1)
            self.assertEqual(len(parts), screen.PER_NODE)  # one node: PER_NODE parts
            self.assertEqual(parts[1][:3], ["n001", "n004", "n007"])  # every n-th name: a clade spread over them
            self.assertEqual(sorted(n for p in parts for n in p), names)  # each name in one part
            self.assertEqual(len(screen.search_parts(names, 2)), 2 * screen.PER_NODE)
            self.assertEqual(len(screen.search_parts(names, 9)), 10)  # no part under SEARCH_MIN_PART
            self.assertEqual(screen.search_parts(names[:5], 9), [names[:5]])
            self.assertEqual(screen.search_parts([], 9), [])
        sizes = {len(p) for p in screen.search_parts(names * 1000, 9)}
        self.assertLessEqual(max(sizes) - min(sizes), 1)  # alike in size

    def test_more_nodes_more_parts_alike(self):
        with mock.patch.object(screen, "SEARCH_MIN_PART", 1):
            many = self.searcher(nodes=3).one("own_genomes", self.seqs, self.db, 0.75, 32, 64)
        self.assertEqual(normal(many), normal(self.plain))
        self.assertEqual(len(os.listdir(os.path.join(self.root, "rrna_own_genomes"))), 3 * screen.PER_NODE)

    def test_each_part_searched_here_logged_with_its_database_time(self):
        said = []
        self.searcher(log=said.append).one("screen 18S", self.seqs, self.seqs, 0.75, 32, 64)
        parts = [s for s in said if "searched in" in s]
        self.assertEqual(len(parts), 3)
        self.assertTrue(all("database read and indexed in" in s for s in parts), parts)
        self.assertTrue(parts[0].startswith("  screen_18S-"))

    def test_vsearch_timed(self):
        q = os.path.join(self.scratch, "q.fna")
        with open(q, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in self.seqs.items())
        out = os.path.join(self.scratch, "o.tsv")
        setup, took = util.vsearch_files(self.vsearch, q, q, out, 0.75, timed=True)
        self.assertTrue(os.path.isfile(out) and 0 <= setup <= took)
        silent = fx.tool(self.tmp.name, "silent", f"#!/bin/sh\ntouch {out}\n")  # no progress lines: no setup time
        self.assertIsNone(util.vsearch_files(silent, q, q, out, 0.75, timed=True)[0])
        broken = fx.tool(self.tmp.name, "broken", "#!/bin/sh\necho 'Fatal error: no memory' >&2\nexit 3\n")
        with self.assertRaisesRegex(OSError, "exit 3.*no memory"):
            util.vsearch_files(broken, q, q, out + "2", 0.75, timed=True)

    def test_kept_in_parts_against_one_database(self):
        got = self.searcher().one("screen 28S 5'", self.seqs, self.db, 0.75, 32, 64)
        self.assertEqual(normal(got), normal(self.plain))
        folder = os.path.join(self.root, "rrna_screen_28S_5_")
        kept = sorted(os.listdir(folder))
        self.assertEqual(len(kept), 3)  # 11 queries in 3 parts (one node)
        self.assertTrue(all(n.endswith(".hits.gz") for n in kept))  # no inputs, database or claim left
        self.assertEqual(os.listdir(self.scratch), [])  # the node's inputs gone
        # read again: vsearch not run
        broken = fx.tool(self.tmp.name, "broken", "#!/bin/sh\nexit 3\n")
        again = self.searcher(vsearch=broken).one("screen 28S 5'", self.seqs, self.db, 0.75, 32, 64)
        self.assertEqual(normal(again), normal(self.plain))
        # the same queries against itself (the screen's): self hits left out, as search's
        alone = screen.search(self.seqs, self.seqs, self.vsearch, 1, self.tmp.name, 0.75, 32, 64)
        self.assertEqual(normal(self.searcher().one("screen 18S", self.seqs, self.seqs, 0.75, 32, 64)), normal(alone))

    def test_inputs_on_the_nodes_disk_where_no_helper_reads_them(self):
        seen = []
        real = util.vsearch_files

        def spy(vsearch, q, db, out, *args, **kw):
            seen.append((q, db, out))
            return real(vsearch, q, db, out, *args, **kw)
        with mock.patch.object(util, "vsearch_files", spy):
            self.searcher().one("own_genomes", self.seqs, self.db, 0.8, 32, 256)
        self.assertEqual(len(seen), 3)
        self.assertTrue(all(p.startswith(self.scratch) for paths in seen for p in paths))
        self.assertEqual(len({db for _q, db, _o in seen}), 1)  # one database file for the three parts

    def test_helpers_search_the_parts_with_their_settings(self):
        started = []

        def helpers(n):  # another node's run of the helper's command: its settings read from todo.tsv
            started.append(n)
            run = subprocess.run([sys.executable, "-m", "eukore.qc", "helper", self.root, "--vsearch", self.vsearch,
                                  "-t", "1"], cwd=os.path.dirname(HERE), capture_output=True, text=True, check=False)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("3 loci searched here", run.stdout)
            self.assertEqual(run.stdout.count("(database read and indexed in "), 3, run.stdout)
            return None
        got = self.searcher(helpers).one("own_published", self.seqs, self.db, 0.75, 32, 64)
        self.assertEqual((started, normal(got)), ([3], normal(self.plain)))
        # a search of fewer than HELPED parts: this node alone
        started.clear()
        few = dict(list(self.seqs.items())[:4])
        self.searcher(helpers).one("own_genomes", few, self.db, 0.75, 32, 64)
        self.assertEqual(started, [])

    def test_follow_up_searches_several_at_once_in_their_order(self):
        items = [({k: self.seqs[k]}, self.db) for k in sorted(self.seqs)]
        got = screen.Searcher(self.vsearch, 4, self.scratch).many("x", items, 0.75, 32, 64)
        self.assertEqual([normal(g) for g in got],
                         [normal(screen.search(q, d, self.vsearch, 1, self.scratch, 0.75, 32, 64)) for q, d in items])
        self.assertEqual(screen.Searcher(self.vsearch).many("x", [], 0.75, 32, 64), [])


class RepresentativesTest(unittest.TestCase):
    def test_owners_typed_by_a_pool(self):
        # since 2026-11-13 each genome's copies typed by a pool (rrna.species_types), each distinct copy sent once with
        # its count: the representatives as one owner after another gave them
        rng = random.Random(9)
        base = "".join(rng.choice("ACGT") for _ in range(400))
        seqs, owner = {}, {}
        for g in range(6):
            unit = tf.mutate(base, 0.03, rng)
            for i in range(rng.randint(1, 25)):
                name = f"g{g}_{i:02d}"
                seqs[name] = unit if i % 3 else tf.mutate(unit, rng.choice((0.002, 0.02)), rng)
                owner[name] = f"G{g}"
        seqs["lone"], owner["lone"] = base, "L"

        def before():  # representatives as it was: types_of of each owner's copies in turn
            by_owner = {}
            for name in sorted(seqs):
                by_owner.setdefault(owner[name], []).append(name)
            out = {}
            for names in by_owner.values():
                types = rrna.types_of([seqs[k] for k in names]) if len(names) > 1 else {seqs[names[0]]: 1}
                rep = {}
                for k in names:
                    t = types[seqs[k]]
                    if t not in rep or len(seqs[k]) > len(seqs[rep[t]]):
                        rep[t] = k
                for k in names:
                    out.setdefault(rep[types[seqs[k]]], []).append(k)
            return out
        self.assertEqual(screen.representatives(seqs, owner, 3), before())
        self.assertEqual(screen.representatives(seqs, owner), before())
        copies = [seqs[k] for k in sorted(seqs)]
        self.assertEqual(rrna.types_of(collections.Counter(copies)), rrna.types_of(copies))


class ChecksTest(unittest.TestCase):
    """The own-species check and the screen give what they gave with the searches of before."""

    def test_the_same_verdicts(self):
        rng = random.Random(31)
        yeast = "".join(rng.choice("ACGT") for _ in range(1200))
        agaric = tf.mutate(yeast, 0.15, rng)
        records = [tf.genome(1, tf.YEAST, tf.mutate(agaric, 0.002, rng)),
                   tf.genome(2, tf.PARADOXUS, tf.mutate(yeast, 0.01, rng)),
                   tf.genome(3, tf.LACTIS, tf.mutate(yeast, 0.03, rng))]
        records += [tf.read(f"EUK_a{i}", tf.AGARIC, tf.mutate(agaric, 0.002, rng)) for i in range(12)]
        records += [tf.read("EUK_b", tf.AGARIC2, tf.mutate(agaric, 0.02, rng))]
        with tempfile.TemporaryDirectory() as tmp, mock.patch.multiple(screen, SEARCH_MIN_PART=3, HELPED=1000):
            ncbi, vsearch = tf.taxonomy(tmp), fx.vsearch(tmp)
            before = rrna.own_species(records, ncbi, vsearch)
            said = []
            searcher = screen.Searcher(vsearch, 2, tmp, folder=os.path.join(tmp, "WORK", "qc", "rrna"),
                                       log=said.append)
            self.assertEqual(rrna.own_species(records, ncbi, vsearch, 2, searcher, said.append), before)
            self.assertEqual(before["GCA_100000001.1_18S_1"], ("f", tf.AGARIC))
            self.assertTrue(any("own_genomes" in s for s in said) and any("suspects" in s for s in said))
            kept = os.listdir(os.path.join(tmp, "WORK", "qc"))
            self.assertEqual(sorted(kept), ["rrna_own_genomes", "rrna_own_published"])
            flags = screen.screen(records, ncbi, vsearch, 2, searcher, said.append)
            self.assertEqual(flags, screen.screen(records, ncbi, vsearch))
            # finish with the searches kept, and again with them read
            rows = [(n, s, dict(r, lsu_complete=False)) for n, s, r in records]
            one = rrna.finish(rows, ncbi, vsearch, 2, log=said.append)
            searches = {"folder": os.path.join(tmp, "WORK", "qc", "rrna2"), "helpers": None}
            self.assertEqual(rrna.finish(rows, ncbi, vsearch, 2, log=said.append, searches=searches, tmp=tmp), one)
            shutil.rmtree(os.path.join(tmp, "WORK", "qc", "rrna"), ignore_errors=True)
            self.assertEqual(rrna.finish(rows, ncbi, vsearch, 2, searches=searches), one)
            self.assertTrue(any(s.startswith("  rRNA own-species check: took") for s in said))


if __name__ == "__main__":
    unittest.main()
