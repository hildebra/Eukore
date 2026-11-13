#!/usr/bin/env python3
"""Tests of the 18S variable regions' and the SSU models' steps (eukore/regions.py; 2026-11-10): the joined models' file
follows the models' list, a batch cmalign fails on is aligned in halves, the sequences it fails on alone are kept as
such for that cmalign, each sequence's SSU model is kept. (The regions' cache as it goes: test_resume.py; the SSU
models' calls: test_screen.py.)

    python3 -m unittest discover -s tests -t tests
"""

import glob
import os
import random
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import regions  # noqa: E402
import fixtures as fx  # noqa: E402

CMALIGN = """#!/bin/sh
# cmalign's stand-in (test_resume.py's, its runs counted in <itself>.calls): each base on the match state of its own
# position; fails on a sequence holding NNNN (as on one over its matrix limit)
while [ $# -gt 0 ]; do case $1 in -o) o=$2;; esac; last=$1; shift; done
grep -c '>' $last >> $0.calls
if grep -q NNNN $last; then exit 1; fi
awk 'BEGIN{n=0} /^>/{name=substr($1,2); next}
     {seq[name]=seq[name] $0; if (length(seq[name])>n) n=length(seq[name]); order[++k]=name}
     END{print "# STOCKHOLM 1.0";
         for (i=1;i<=k;i++) {s=seq[order[i]]; while (length(s)<n) s=s "-"; print order[i] " " s}
         rf=""; for (i=0;i<n;i++) rf=rf "x"; print "#=GC RF " rf; print "//"}' $last > $o
"""
FAILS = "if grep -q NNNN $last; then exit 1; fi"

CMSEARCH = """#!/bin/sh
# cmsearch --tblout's stand-in (its runs' sequences counted in <itself>.calls): a sequence of A first is the eukaryote
# model's, one of C first the bacterial model's by 40 bits, any other no model's
while [ "$#" -gt 0 ]; do [ "$1" = "--tblout" ] && out=$2; last=$1; shift; done
grep -c '>' "$last" >> "$0.calls"
awk '/^>/{name=substr($1,2); next}
     /^A/{print name " - SSU_rRNA_eukarya RF01960 hmm 1 1800 1 1800 + - 6 0.5 0 1500.0 0 ! -"}
     /^C/{print name " - SSU_rRNA_eukarya RF01960 hmm 1 1800 1 1500 + - 6 0.5 0 700.0 0 ! -"
          print name " - SSU_rRNA_bacteria RF00177 hmm 1 1500 1 1500 + - 6 0.5 0 740.0 0 ! -"}' "$last" > "$out"
"""


def calls(tool):
    """The sequences of each run of a stand-in so far."""
    try:
        with open(tool + ".calls") as fh:
            return [int(x) for x in fh.read().split()]
    except OSError:
        return []


def lines(pattern):
    out = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as fh:
            out += fh.read().splitlines()
    return out


class FetchTest(unittest.TestCase):
    def test_the_joined_models_follow_their_list(self):
        # until 2026-11-10 the joined file was kept where it was an Infernal model: a model added to SSU_MODELS, or one
        # fetched again, never reached cmsearch
        fetched = []

        def download(url, path):
            fetched.append(url.split("/")[-2])
            with open(path, "w", newline="\n") as fh:
                fh.write(f"INFERNAL1/a [stand-in]\nNAME {fetched[-1]}\n//\n")

        def names(path):
            with open(path) as fh:
                return [line.split()[1] for line in fh if line.startswith("NAME")]

        with tempfile.TemporaryDirectory() as src:
            said = []
            self.assertEqual(regions.fetch(src, download, said.append), regions.cm_path(src))
            joined = regions.models_path(src)
            listed = [acc for acc, _kind in regions.SSU_MODELS]
            self.assertEqual((fetched, names(joined)), (listed, listed))
            first = os.stat(joined)
            regions.fetch(src, download, said.append)  # nothing new: nothing fetched, the file as it is
            self.assertEqual((fetched, os.stat(joined).st_ino, os.stat(joined).st_mtime_ns),
                             (listed, first.st_ino, first.st_mtime_ns))
            with mock.patch.object(regions, "SSU_MODELS", regions.SSU_MODELS + (("RF99999", "archaea"),)):
                regions.fetch(src, download, said.append)
                self.assertEqual((fetched, names(joined)), (listed + ["RF99999"], listed + ["RF99999"]))
            regions.fetch(src, download, said.append)  # the list as before: the file too
            self.assertEqual(names(joined), listed)
            with open(os.path.join(src, "rfam", f"{listed[1]}.cm"), "w", newline="\n") as fh:  # a model's new release
                fh.write("INFERNAL1/a [stand-in]\nNAME new\n//\n")
            regions.fetch(src, download, said.append)
            self.assertEqual(names(joined), [listed[0], "new"] + listed[2:])
            self.assertTrue(regions.is_cm(joined))
            self.assertEqual(os.listdir(os.path.join(src, "rfam")).count(regions.SSU_MODELS_FILE), 1)
            self.assertFalse([f for f in os.listdir(os.path.join(src, "rfam")) if f.endswith(".partial")])


class CmalignTest(unittest.TestCase):
    def setUp(self):
        rng = random.Random(5)
        self.seqs = {f"s{i:02d}": "".join(rng.choice("ACGT") for _ in range(150 + i)) for i in range(40)}
        for name in ("s03", "s37"):  # one of the first batch of 16, one of the last of 8
            self.seqs[name] = self.seqs[name][:50] + "NNNN" + self.seqs[name][54:]
        self.good = sorted(set(self.seqs) - {"s03", "s37"})

    def test_a_failing_batch_in_halves(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(regions, "BATCH", 16):
            tool = fx.tool(tmp, "cmalign", CMALIGN)
            alone = {}
            for name in self.good:  # each alone: what the halves must give too
                alone.update(regions.cmalign("x.cm", {name: self.seqs[name]}, tool, tmp=tmp))
            self.assertEqual(sorted(alone), self.good)
            os.remove(tool + ".calls")
            kept, failed, said = [], {}, []
            got = regions.cmalign("x.cm", self.seqs, tool, tmp=tmp, log=said.append, each=kept.append,
                                  fail=failed.update)
            self.assertEqual(got, alone)
            self.assertEqual(got["s00"][:3] + got["s00"][-1:], [1, 2, 3, 150])
            # each sequence given to each once, as its run ends; the two that fail alone to fail
            self.assertEqual(sorted(n for part in kept for n in part), self.good)
            self.assertEqual({n: got[n] for part in kept for n in part}, alone)
            self.assertEqual(sorted(failed), ["s03", "s37"])
            self.assertIn("2 it could not align left out", said[0])
            # 16 in 1 + 2 * 4 runs, 16 in 1, 8 in 1 + 2 * 3 (one sequence at a time: 3 + 16 + 8)
            runs = calls(tool)  # a half that fails: its first half next
            self.assertEqual(runs, [16, 8, 4, 2, 2, 1, 1, 4, 8, 16, 8, 4, 4, 2, 1, 1, 2])

    def test_a_killed_cmalign_says_nothing_of_its_sequence(self):
        for end in ("kill -9 $$", "exit 137"):  # a signal's, a shell's word of one
            with tempfile.TemporaryDirectory() as tmp:
                tool = fx.tool(tmp, "cmalign", CMALIGN.replace(FAILS, f"if grep -q NNNN $last; then {end}; fi"))
                failed = {}
                got = regions.cmalign("x.cm", self.seqs, tool, tmp=tmp, fail=failed.update)
                self.assertEqual((sorted(got), failed), (self.good, {}))


class UnalignedTest(unittest.TestCase):
    def setUp(self):
        ref = regions.reference()
        self.records = [("r1", ref[:1700], dict(kind="18S")),
                        ("r2", ref[:800] + "NNNN" + ref[804:1700], dict(kind="18S")),
                        ("r3", ref[5:1700], dict(kind="18S"))]

    def model(self, tmp):
        work, cm = os.path.join(tmp, "work"), os.path.join(tmp, "ssu.cm")
        os.makedirs(work)
        with open(cm, "w") as fh:
            fh.write("stand-in\n")
        return work, cm

    def test_a_sequence_cmalign_fails_on_is_not_aligned_again(self):
        # until 2026-11-10 every run aligned it again (its batch, then the batch's sequences one at a time)
        with tempfile.TemporaryDirectory() as tmp:
            work, cm = self.model(tmp)
            tool = fx.tool(tmp, "cmalign", CMALIGN)
            rows = regions.run(self.records, cm, tool, 1, work, log=lambda *_: None)
            self.assertEqual({r["record"] for r in rows}, {"r1", "r3"})
            unaligned = os.path.join(work, "cmalign", "unaligned.*.tsv")
            [line] = lines(unaligned)
            self.assertEqual(line.split("\t")[0], regions._digest(self.records[1][1]))  # noqa: SLF001
            self.assertEqual(len(lines(os.path.join(work, "cmalign", "ssu.*.tsv"))), 2)  # not among the aligned
            os.remove(tool + ".calls")
            said = []
            again = regions.run(self.records, cm, tool, 1, work, log=said.append)
            self.assertEqual((again, calls(tool)), (rows, [1]))  # the reference alone
            self.assertIn("1 distinct 18S cmalign failed on before left out", said[0])
            # another record of the same 18S: not aligned either; a new 18S is
            more = self.records + [("r4", self.records[1][1].lower(), dict(kind="18S")),
                                   ("r5", self.records[0][1][:1600], dict(kind="18S"))]
            os.remove(tool + ".calls")
            rows = regions.run(more, cm, tool, 1, work, log=lambda *_: None)
            self.assertEqual(({r["record"] for r in rows}, calls(tool)), ({"r1", "r3", "r5"}, [1, 1]))
            # another cmalign (here the file's new size and time): tried again, and aligned
            fx.tool(tmp, "cmalign", CMALIGN.replace(FAILS, ""))
            rows = regions.run(more, cm, tool, 1, work, log=lambda *_: None)
            self.assertEqual({r["record"] for r in rows}, {"r1", "r2", "r3", "r4", "r5"})
            # the model's change: its own failures (and its own regions)
            with open(cm, "a") as fh:
                fh.write("changed\n")
            self.assertEqual(len(glob.glob(unaligned)), 1)
            fx.tool(tmp, "cmalign", CMALIGN)
            os.utime(tool, ns=(1, 1))  # the first cmalign's size, another time
            regions.run(self.records, cm, tool, 1, work, log=lambda *_: None)
            self.assertEqual(len(glob.glob(unaligned)), 2)

    def test_a_failure_that_is_not_the_sequences_is_not_kept(self):
        for name, fails in (("killed", "if grep -q NNNN $last; then kill -9 $$; fi"),
                            # a cmalign that fails on all once it failed (a disk that filled): the reference too
                            ("full", "if [ -e $0.full ] || grep -q NNNN $last; then touch $0.full; exit 1; fi")):
            with tempfile.TemporaryDirectory() as tmp:
                work, cm = self.model(tmp)
                tool = fx.tool(tmp, "cmalign", CMALIGN.replace(FAILS, fails))
                rows = regions.run(self.records, cm, tool, 1, work, log=lambda *_: None)
                self.assertEqual({r["record"] for r in rows}, {"r1", "r3"} if name == "killed" else set())
                self.assertEqual(glob.glob(os.path.join(work, "cmalign", "unaligned.*.tsv")), [], name)
                fixed = fx.tool(tmp, "cmalign", CMALIGN.replace(FAILS, ""))
                rows = regions.run(self.records, cm, fixed, 1, work, log=lambda *_: None)
                self.assertEqual({r["record"] for r in rows}, {"r1", "r2", "r3"})

    def test_no_work_folder_no_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            _work, cm = self.model(tmp)
            tool = fx.tool(tmp, "cmalign", CMALIGN)
            rows = regions.run(self.records, cm, tool, 1, None, log=lambda *_: None)
            self.assertEqual({r["record"] for r in rows}, {"r1", "r3"})
            self.assertEqual(sorted(os.listdir(tmp)), ["cmalign", "cmalign.calls", "ssu.cm", "work"])


class AlignedTogetherTest(unittest.TestCase):
    """The 18S alignments several cmalign runs at once and in kept parts shared with the helpers (since 2026-11-13)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = self.tmp.name
        ref = regions.reference()
        self.records = [(f"r{i:02d}", ref[i * 7:1700 - i], dict(kind="18S")) for i in range(14)]
        self.records.append(("bad", ref[:800] + "NNNN" + ref[804:1700], dict(kind="18S")))
        self.records.append(("long", ref * 3, dict(kind="18S")))
        self.cm = os.path.join(tmp, "ssu.cm")
        with open(self.cm, "w") as fh:
            fh.write("stand-in\n")
        self.tool = fx.tool(tmp, "cmalign", CMALIGN)
        self.plain = regions.run(self.records, self.cm, fx.tool(tmp, "plain", CMALIGN), 1, None, log=lambda *_: None)
        self.local = os.path.join(tmp, "local")
        os.makedirs(self.local)
        self.patch = mock.patch.multiple(regions, BATCH=2, SSU_HELPED=2)
        self.patch.start()

    def tearDown(self):
        self.patch.stop()
        self.tmp.cleanup()

    def work(self, name):
        path = os.path.join(self.tmp.name, name)
        os.makedirs(path)
        return path

    def searches(self, work, helpers=None):
        return {"folder": os.path.join(work, "qc", "cmalign"), "helpers": helpers, "nodes": 1, "scratch": self.local}

    def test_several_runs_at_once_alike(self):
        said = []
        rows = regions.run(self.records, self.cm, self.tool, 8, self.work("w"), log=said.append)  # 2 runs of 4
        self.assertEqual(rows, self.plain)
        self.assertEqual({r["record"] for r in rows}, {f"r{i:02d}" for i in range(14)})
        self.assertTrue(any("1 18S over 4000 bases and 1 it could not align left out" in s for s in said), said)
        self.assertTrue(any("15 18S aligned (cmalign): took" in s for s in said), said)
        self.assertEqual(len(lines(os.path.join(self.tmp.name, "w", "cmalign", "unaligned.*.tsv"))), 1)

    def test_kept_parts_alike_and_removed(self):
        work, said = self.work("w"), []
        rows = regions.run(self.records, self.cm, self.tool, 1, work, log=said.append, searches=self.searches(work))
        self.assertEqual(rows, self.plain)
        self.assertFalse(os.path.exists(os.path.join(work, "qc", "cmalign")))  # in the cache now
        self.assertEqual(os.listdir(self.local), [])
        self.assertTrue(any("(cmalign, 3 parts): took" in s for s in said), said)
        self.assertEqual(sum("cmalign-" in s and "searched in" in s for s in said), 3)
        unaligned = lines(os.path.join(work, "cmalign", "unaligned.*.tsv"))
        self.assertEqual([u.split("\t")[0] for u in unaligned], [regions._digest(self.records[14][1])])  # noqa: SLF001
        os.remove(self.tool + ".calls")
        again = regions.run(self.records, self.cm, self.tool, 1, work, log=lambda *_: None,
                            searches=self.searches(work))
        self.assertEqual((again, calls(self.tool)), (rows, [1]))  # the reference alone

    def test_the_helpers_align_the_parts(self):
        work, started = self.work("w"), []

        def helpers(n):  # another node's run of the helper's command: the model and its columns read from todo.tsv
            started.append(n)
            env = dict(os.environ, PYTHONPATH=os.path.dirname(HERE))
            code = ("import sys; from unittest import mock; from eukore import qc, regions\n"
                    "with mock.patch.object(regions, 'BATCH', 2):\n"
                    "    sys.exit(qc.main(sys.argv[1:]))\n")
            run = subprocess.run([sys.executable, "-c", code, "helper", os.path.join(work, "qc"), "--vsearch",
                                  self.tool, "--cmalign", self.tool, "-t", "8"], capture_output=True, text=True,
                                 check=False, env=env)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("3 loci searched here", run.stdout)
            return None
        rows = regions.run(self.records, self.cm, self.tool, 1, work, log=lambda *_: None,
                           searches=self.searches(work, helpers))
        self.assertEqual((started, rows), ([3], self.plain))

    def test_columns_text(self):
        cols = {"V4": (560, 1010), "V9": (1620, 1700)}
        self.assertEqual(regions.cols_of(regions.cols_text(cols)), cols)
        from eukore import qc
        params = qc.CmalignParams("/x/ssu.cm", "abc", regions.cols_text(cols))
        self.assertEqual(qc.todo_entry(qc.todo_line("cmalign-1.k", 4, params)), ("cmalign-1.k", params, None))


class ModelsCacheTest(unittest.TestCase):
    def test_a_sequences_model_is_kept(self):
        # until 2026-11-10 every package stage searched every 18S again
        seqs = {"euk": "ACGT" * 10, "bac": "CCGT" * 10, "none": "TCGT" * 10, "twin": "acgt" * 10, "empty": ""}
        calls_of = {"euk": "eukarya", "twin": "eukarya", "bac": "bacteria", "none": "-"}
        with tempfile.TemporaryDirectory() as tmp:
            work, cm = os.path.join(tmp, "work"), os.path.join(tmp, "models.cm")
            os.makedirs(work)
            with open(cm, "w") as fh:
                fh.write("stand-in\n")
            tool = fx.tool(tmp, "cmsearch", CMSEARCH)
            self.assertEqual(regions.ssu_models(seqs, cm, tool, work=work), calls_of)
            kept = os.path.join(work, "cmalign", "models.*.tsv")
            self.assertEqual((calls(tool), len(lines(kept))), ([3], 3))  # each distinct sequence, - too
            self.assertEqual(regions.ssu_models(seqs, cm, tool, work=work), calls_of)
            self.assertEqual(calls(tool), [3])  # none searched again
            more = dict(seqs, new="CCCC" * 10)
            self.assertEqual(regions.ssu_models(more, cm, tool, work=work), dict(calls_of, new="bacteria"))
            self.assertEqual((calls(tool), len(lines(kept))), ([3, 1], 4))  # the new one alone
            # the build gives its WORK as tmp: the cache there
            self.assertEqual(regions.ssu_models(more, cm, tool, tmp=work), dict(calls_of, new="bacteria"))
            self.assertEqual(calls(tool), [3, 1])
            # a line a killed run cut: that sequence searched again
            [path] = glob.glob(kept)
            with open(path) as fh:
                text = fh.read()
            with open(path, "w", newline="\n") as fh:
                fh.write(text[:-4])
            self.assertEqual(regions.ssu_models(more, cm, tool, work=work), dict(calls_of, new="bacteria"))
            self.assertEqual((calls(tool), len(lines(kept))), ([3, 1, 1], 5))
            # other models (the file, the list or the margin): a cache of their own
            with mock.patch.object(regions, "SSU_MARGIN", 50.0):
                self.assertEqual(regions.ssu_models(seqs, cm, tool, work=work), dict(calls_of, bac="eukarya"))
            self.assertEqual(regions.ssu_models(seqs, cm, tool, work=work), calls_of)
            with open(cm, "a") as fh:
                fh.write("another release\n")
            self.assertEqual(regions.ssu_models(seqs, cm, tool, work=work), calls_of)
            self.assertEqual((calls(tool), len(glob.glob(kept))), ([3, 1, 1, 3, 3], 3))
            # a cmsearch that fails: said, nothing kept
            broken = fx.tool(tmp, "cmsearch_fails", "#!/bin/sh\nexit 1\n")
            with self.assertRaises(RuntimeError):
                regions.ssu_models({"x": "ACGA" * 10}, cm, broken, work=work)
            self.assertEqual(len(lines(kept)), 5 + 3 + 3)
            # no folder: no cache
            os.remove(tool + ".calls")
            for _ in range(2):
                self.assertEqual(regions.ssu_models(seqs, cm, tool), calls_of)
            self.assertEqual(calls(tool), [3, 3])


class ModelsPartsTest(unittest.TestCase):
    """The SSU models' search in parts kept in WORK/qc/ssu_models/ and shared with the helpers (since 2026-11-13)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        tmp = self.tmp.name
        rng = random.Random(3)
        self.seqs = {f"r{i:02d}": rng.choice("AC") + "".join(rng.choice("ACGT") for _ in range(39)) for i in range(12)}
        self.seqs["none"] = "T" + "G" * 39
        self.work, self.cm = os.path.join(tmp, "WORK"), os.path.join(tmp, "models.cm")
        os.makedirs(self.work)
        with open(self.cm, "w") as fh:
            fh.write("stand-in\n")
        self.local = os.path.join(tmp, "local")
        os.makedirs(self.local)
        self.folder = os.path.join(self.work, "qc", "ssu_models")
        self.tool = fx.tool(tmp, "cmsearch", CMSEARCH)
        self.plain = regions.ssu_models(self.seqs, self.cm, fx.tool(tmp, "plain", CMSEARCH))  # no cache, one run
        self.patch = mock.patch.multiple(regions, SSU_PART_MIN=1, SSU_HELPED=2)
        self.patch.start()

    def tearDown(self):
        self.patch.stop()
        self.tmp.cleanup()

    def searches(self, helpers=None, nodes=1):
        return {"folder": self.folder, "helpers": helpers, "nodes": nodes, "scratch": self.local}

    def test_parts(self):
        names = [f"n{i:03d}" for i in range(100)]
        with mock.patch.multiple(regions, SSU_PART_MIN=10, SSU_PART_MAX=40):
            self.assertEqual(len(regions.ssu_parts(names, 1)), regions.SSU_PER_NODE)
            self.assertEqual(len(regions.ssu_parts(names, 3)), 9)
            self.assertEqual(len(regions.ssu_parts(names, 20)), 10)  # none under SSU_PART_MIN
            self.assertEqual(sorted(n for p in regions.ssu_parts(names, 2) for n in p), names)
            self.assertEqual(regions.ssu_parts(names, 1)[1][:2], ["n001", "n004"])  # interleaved
        with mock.patch.multiple(regions, SSU_PART_MIN=1, SSU_PART_MAX=10):
            self.assertEqual(len(regions.ssu_parts(names, 1)), 10)  # each part SSU_PART_MAX at most
        self.assertEqual(regions.ssu_parts([], 3), [])

    def test_kept_parts_give_the_same_calls(self):
        said = []
        got = regions.ssu_models(self.seqs, self.cm, self.tool, work=self.work, searches=self.searches(),
                                 log=said.append)
        self.assertEqual(got, self.plain)
        self.assertEqual((len(calls(self.tool)), sum(calls(self.tool))), (3, 13))  # 3 parts, each sequence once
        self.assertFalse(os.path.exists(self.folder))  # the calls in the cache: the parts removed
        self.assertEqual(os.listdir(self.local), [])
        self.assertEqual(sum("ssu-" in s and "searched in" in s for s in said), 3)
        self.assertIn("SSU models: 13 18S searched (cmsearch, 3 parts)", said[-1])
        # the cache read again: no search
        self.assertEqual(regions.ssu_models(self.seqs, self.cm, self.tool, work=self.work, searches=self.searches()),
                         self.plain)
        self.assertEqual(len(calls(self.tool)), 3)

    def test_a_killed_run_keeps_the_parts_done(self):
        # the third part's search fails (as a killed node): the two searched before it are kept
        failing = fx.tool(self.tmp.name, "failing", CMSEARCH.replace(
            'grep -c', '[ -f "$0.calls" ] && [ "$(wc -l < "$0.calls")" -ge 2 ] && exit 1\ngrep -c'))
        with self.assertRaises(RuntimeError):
            regions.ssu_models(self.seqs, self.cm, failing, work=self.work, searches=self.searches(nodes=3))
        done = sum(calls(failing))
        kept = [n for n in os.listdir(self.folder) if n.endswith(".hits.gz")]
        self.assertEqual((len(kept), len(calls(failing))), (2, 2))
        got = regions.ssu_models(self.seqs, self.cm, self.tool, work=self.work, searches=self.searches(nodes=3))
        self.assertEqual(got, self.plain)
        self.assertEqual(sum(calls(self.tool)), 13 - done)  # those of the parts kept not searched again

    def test_the_helpers_search_the_parts(self):
        started = []

        def helpers(n):  # another node's run of the helper's command: the models' file read from todo.tsv
            started.append(n)
            run = subprocess.run([sys.executable, "-m", "eukore.qc", "helper", os.path.dirname(self.folder),
                                  "--vsearch", self.tool, "--cmsearch", self.tool, "-t", "1"],
                                 cwd=os.path.dirname(HERE), capture_output=True, text=True, check=False)
            self.assertEqual(run.returncode, 0, run.stderr)
            self.assertIn("3 loci searched here", run.stdout)
            return None
        got = regions.ssu_models(self.seqs, self.cm, self.tool, work=self.work, searches=self.searches(helpers))
        self.assertEqual((started, got), ([3], self.plain))
        self.assertEqual(sum(calls(self.tool)), 13)

    def test_a_helper_without_cmsearch_leaves_the_parts(self):
        from eukore import qc
        part = {"x": "ACGT" * 10}
        params = qc.CmParams(self.cm, "abc")
        base = qc.search_base("ssu-1", part, None, params)
        os.makedirs(self.folder)
        qc.write_inputs(self.folder, base, part)
        with open(os.path.join(self.folder, "todo.tsv"), "w") as fh:
            fh.write(qc.todo_line(base, 1, params))
        self.assertEqual(qc.todo_entry(qc.todo_line(base, 1, params)), (base, params, None))
        self.assertEqual(qc.help_search(os.path.dirname(self.folder), self.tool, 1, log=lambda *_: None), 0)
        self.assertEqual(qc.help_search(os.path.dirname(self.folder), self.tool, 1, log=lambda *_: None,
                                        cmsearch=self.tool), 1)
        self.assertEqual(qc.read_hits(self.folder, base, params), [("x", "RF01960", 1500.0)])  # A first: eukaryote
        self.assertNotEqual(qc.search_key(part, None, params), qc.search_key(part, None, qc.CmParams(self.cm, "abd")))


if __name__ == "__main__":
    unittest.main()
