#!/usr/bin/env python3
"""Tests of the process stage on a large node (2026-11-08): genomes begun as their memory fits, each with threads by
its share of it, the intraspecific panel's in the same pool; hmmsearch in runs of a share of the models; the rescue's
ORFs read in blocks and made in pieces; a run's temporary files in --scratch.

    python3 -m unittest discover -s tests -t tests
"""

import os
import shutil
import stat
import sys
import tempfile
import time
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import eukore_build  # noqa: E402
from eukore import genomes, loci, placement, util  # noqa: E402


def fake_run(item):
    """genomes.run's stand-in: its threads and when it ran."""
    job, settings = item
    began = time.time()
    time.sleep(0.3)
    return job["accession"], {"threads": settings["threads"], "began": began, "ended": time.time()}, None


def killed_run(item):
    """genomes.run's stand-in whose worker is killed (the OOM killer) on the genome "killed"."""
    if item[0]["accession"] == "killed":
        os.kill(os.getpid(), 9)
    return fake_run(item)


class Quiet:
    def tick(self, note=""):
        pass


@unittest.skipUnless(sys.platform.startswith("linux"), "fork's pools")
class ScheduleTest(unittest.TestCase):
    def run_all(self, opts, sizes, run=fake_run):
        old = eukore_build.run_one
        eukore_build.run_one = run  # submitted by name: a module-level stand-in pickles under fork and forkserver
        try:
            todo = [({"accession": a}, {"threads": opts.threads_per_genome, "rrna": None})
                    for a in sorted(sizes, key=lambda a: -sizes[a])]
            return {acc: s for acc, s, _why in eukore_build.schedule(opts, todo, sizes, Quiet(), lambda: "", poll=1)}
        finally:
            eukore_build.run_one = old

    def test_memory_and_threads(self):
        # 8 cores, 4 at once at 2 threads each; 30 GB: a 3 Gb genome (~25.5 GB) at a time, on more threads; the small
        # ones (~3.75 GB) beside it where threads are free, 4 at once after
        opts = types.SimpleNamespace(threads=8, threads_per_genome=2, jobs=4, memory=30e9)
        sizes = {"big1": 3e9, "big2": 3e9, **{f"small{i}": 1e8 for i in range(6)}}
        got = self.run_all(opts, sizes)
        self.assertEqual(set(got), set(sizes))
        self.assertGreater(got["big1"]["threads"], 2)  # memory's share: round(8 * 25.5 / 30) = 7

        def overlap(a, b):
            return got[a]["began"] < got[b]["ended"] and got[b]["began"] < got[a]["ended"]
        self.assertFalse(overlap("big1", "big2"))  # never both: 51 GB
        for t in sorted({s["began"] for s in got.values()} | {s["ended"] - 1e-6 for s in got.values()}):
            on = [a for a, s in got.items() if s["began"] <= t < s["ended"]]
            self.assertLessEqual(sum(eukore_build.run_memory(sizes[a]) for a in on), 30e9 + 1)
            self.assertLessEqual(sum(got[a]["threads"] for a in on), 8)
            self.assertLessEqual(len(on), 4)
        # no memory known: -j at once, -t / -j threads each
        opts.memory = None
        got = self.run_all(opts, sizes)
        self.assertEqual({s["threads"] for s in got.values()}, {2})
        # a genome over the memory: run alone
        opts.memory = 10e9
        got = self.run_all(opts, {"huge": 5e9, "small": 1e8})
        self.assertFalse(got["huge"]["began"] < got["small"]["ended"] and got["small"]["began"] < got["huge"]["ended"])
        # a worker killed: the build stops (until 2026-11-08 the scheduler waited for its genome forever)
        with self.assertRaises(SystemExit):
            self.run_all(opts, {"killed": 1e8, "small": 1e8}, killed_run)

    def test_memory_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            proc, cgroup = os.path.join(tmp, "proc"), os.path.join(tmp, "cgroup")
            os.makedirs(proc)
            os.makedirs(cgroup)
            with open(os.path.join(proc, "meminfo"), "w") as fh:
                fh.write("MemTotal: 263000000 kB\nMemAvailable: 250000000 kB\n")
            self.assertEqual(util.memory_budget({"SLURM_MEM_PER_NODE": "204800"}, proc, cgroup),
                             (204800 * 2 ** 20, "SLURM_MEM_PER_NODE"))
            self.assertEqual(util.memory_budget({"SLURM_MEM_PER_CPU": "2000", "SLURM_CPUS_ON_NODE": "80"}, proc,
                                                cgroup), (2000 * 2 ** 20 * 80, "SLURM_MEM_PER_CPU"))
            with open(os.path.join(cgroup, "memory.max"), "w") as fh:
                fh.write("100000000000\n")
            self.assertEqual(util.memory_budget({}, proc, cgroup), (100000000000, "the cgroup's limit"))
            with open(os.path.join(cgroup, "memory.max"), "w") as fh:
                fh.write("max\n")
            self.assertEqual(util.memory_budget({}, proc, cgroup), (250000000 * 1024, "the node's memory available"))
            # an allocation over what the node has free: the node's
            self.assertEqual(util.memory_budget({"SLURM_MEM_PER_NODE": "512000"}, proc, cgroup)[1],
                             "the node's memory available")
            self.assertEqual(util.memory_budget({}, os.path.join(tmp, "none"), cgroup), (None, None))

    def test_threads_by_memory(self):
        opts = types.SimpleNamespace(threads=120, threads_per_genome=4, memory=170e9)
        self.assertEqual(eukore_build.genome_threads(opts, eukore_build.run_memory(1e8)), 4)
        self.assertEqual(eukore_build.genome_threads(opts, eukore_build.run_memory(3e9)), 18)  # 25.5 of 170 GB
        self.assertEqual(eukore_build.genome_threads(opts, eukore_build.run_memory(30e9)), 23)  # pieces of 4 Gb
        opts.memory = None
        self.assertEqual(eukore_build.genome_threads(opts, 1e12), 4)


HMM = "HMMER3/f [3.4]\nNAME  {}\nLENG  10\nHMM ...\n//\n"
FAKE_HMMSEARCH = """#!/bin/sh
# hmmsearch's stand-in: a domtblout line for each model of the HMM file (the one but last argument) and each target
while [ $# -gt 2 ]; do case $1 in --domtblout) out=$2; shift;; esac; shift; done
hmm=$1; seqs=$2
echo "# target" > $out
for name in $(grep '^NAME' $hmm | awk '{print $2}'); do
  for t in $(grep '>' $seqs | sed 's/>//' | awk '{print $1}'); do
    echo "$t - 100 $name - 10 1e-9 42.0 0 1 1 1 1e-9 1e-9 40.0 2 8 1 90 1 90 0.9 -" >> $out
  done
done
"""


class HmmsearchRunsTest(unittest.TestCase):
    def test_runs_of_a_share_of_the_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            hmm = os.path.join(tmp, "db.hmm")
            names = [f"m{i}" for i in range(7)]
            with open(hmm, "w") as fh:
                fh.writelines(HMM.format(n) for n in names)
            shards = loci.hmm_shards(hmm, 3)
            self.assertEqual(len(shards), 3)
            with open(shards[1]) as fh:
                self.assertEqual([l.split()[1] for l in fh if l.startswith("NAME")], ["m1", "m4"])
            self.assertEqual(loci.hmm_shards(hmm, 3), shards)  # made once
            self.assertEqual(len(loci.hmm_shards(hmm, 20)), 7)  # no more than the models
            self.assertEqual(loci.hmm_shards(hmm, 1), [hmm])
            fake = os.path.join(tmp, "hmmsearch")
            with open(fake, "w") as fh:
                fh.write(FAKE_HMMSEARCH)
            os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
            proteins = {"p1": "MKV", "p2": "MAL*"}
            one = loci.hmmsearch(proteins, hmm, fake, 1, tmp)
            self.assertEqual(loci.hmmsearch(proteins, hmm, fake, 3, tmp), one)
            self.assertEqual(set(one["p1"]), set(names))
            target = os.path.join(tmp, "t.faa")
            with open(target, "w") as fh:
                fh.write(">a\nMK\n>b\nMK\n")
            got = loci.hmmsearch_runs(hmm, target, os.path.join(tmp, "dom"), 3, fake)
            self.assertEqual([l.split()[3] for l in loci.domtbl_lines(got)], [n for n in names for _t in "ab"])


class OrfTest(unittest.TestCase):
    TEXT = (">orf1 source=c1 coords=1..120 length=40 frame=1  desc=\nMKPGFKPGFKPGFKPGFKPGFKPGFK\nPGFKPGFKPGFKPG\n"
            ">orf2 source=c1 coords=900..781 length=40 frame=4  desc=\nMKPGFKPGFKPGFKPGFKPGFKPGXXXXXXXXXXXXXXXX"
            "XXXXXXXXXXXXXXXXXXXXXKPGFKPGFKPGFKPGFKPGFK\n>orf3 source=c2 coords=4..30 length=9 frame=2  desc=\nMKPGFKPGF\n")

    def test_blocks_and_pieces(self):
        import io
        whole = io.StringIO()
        loci.split_orfs(self.TEXT.splitlines(True), whole)
        for block in (1, 7, 64, 1 << 20):  # a record cut anywhere by the blocks
            out = io.StringIO()
            for rec in loci._records(io.StringIO(self.TEXT), block):
                loci.split_orfs([">" + rec.lstrip(">") + "\n"], out)
            self.assertEqual(out.getvalue(), whole.getvalue(), block)
        self.assertEqual(whole.getvalue().count(">"), 3)  # orf1 whole, orf2 in two pieces, orf3 too short
        self.assertIn(">orf1.0 source=c1 coords=1..120 length=40\n", whole.getvalue())
        out = io.StringIO()
        loci.split_orfs(io.StringIO(self.TEXT), out, "p2.")
        self.assertEqual(out.getvalue(), whole.getvalue().replace(">orf", ">p2.orf"))

    def test_a_genomes_pieces_stay_under_their_limit(self):
        # chromosomes of 3 units and a limit of 4: a piece each (until 2026-11-08 a piece took the next whole one,
        # 6 units: a giant genome's miniprot index over the memory the build reserved); a contig over it alone
        with tempfile.TemporaryDirectory() as tmp:
            genome = os.path.join(tmp, "g.fna")
            with open(genome, "w") as fh:
                fh.writelines(f">c{i}\n" + "ACGT" * (n // 4) + "\n" for i, n in enumerate((3000, 3000, 1000, 9000, 500)))
            pieces = placement.chunks(genome, tmp, max_bases=4000)
            self.assertEqual([[n for _c, n in placement.contig_lengths(p)] for p in pieces],
                             [[3000], [3000, 1000], [9000], [500]])

    @unittest.skipUnless(shutil.which("esl-translate"), "esl-translate")
    def test_pieces_as_one_run(self):
        import random
        rnd = random.Random(5)
        with tempfile.TemporaryDirectory() as tmp:
            genome = os.path.join(tmp, "g.fna")
            with open(genome, "w") as fh:
                for i in range(9):
                    fh.write(f">c{i} x\n" + "".join(rnd.choice("ACGT") for _ in range(rnd.randint(500, 20000)))
                             + "N" * 300 + "ATGAAA" * 100 + "\n")
            got = {}
            for threads in (1, 4):
                orfs = os.path.join(tmp, f"o{threads}.faa")
                loci.orf_file(genome, orfs, 1, threads, tmp)
                with open(orfs) as fh:
                    got[threads] = [line.split(None, 1)[1] if line.startswith(">") else line for line in fh]
            self.assertTrue(got[1])
            self.assertEqual(got[4], got[1])  # the same ORFs in the same order, the names aside
            self.assertEqual(os.listdir(os.path.join(tmp, "orf_pieces")), [])


class ScratchTest(unittest.TestCase):
    def test_a_run_folder_and_the_stale_ones(self):
        with tempfile.TemporaryDirectory() as tmp:
            work, scratch = os.path.join(tmp, "work"), os.path.join(tmp, "scratch")
            self.assertIsNone(genomes.scratch_dir({}, "GCA_1.1"))
            mine = genomes.scratch_dir({"scratch": scratch}, "GCA_1.1")
            other = genomes.scratch_dir({"scratch": scratch}, "GCA_2.1")
            self.assertTrue(os.path.isdir(mine) and os.path.dirname(mine) == scratch)
            with genomes.claimed(work, "GCA_1.1"):
                self.assertEqual(genomes.stale_scratch(work, scratch), [other])  # GCA_1.1's run goes on
            self.assertEqual(sorted(genomes.stale_scratch(work, scratch)), sorted([mine, other]))


if __name__ == "__main__":
    unittest.main()
