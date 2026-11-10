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


if __name__ == "__main__":
    unittest.main()
