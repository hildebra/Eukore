#!/usr/bin/env python3
"""Tests of the process and package stages' speed and resume fixes (2026-11-08): stretches of a genome read alone, the
rRNA seeds in one pass, FCS-GX on a transcriptome's windows, genomes claimed by one build at a time, steps a missing tool
left out made later, the cmalign cache kept as it goes and without failures, copies held compressed, the organelle
molecules the metadata scan kept, the fetch stage's note of genomes without an annotation.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import random
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import annotation, genomes, intraspecific, organelle, package, placement, regions, rrna  # noqa: E402
import fixtures as fx  # noqa: E402


def old_windows(seq):
    """rrna.windows as it was until 2026-11-08 (a str.find per seed): the reference the one-pass scan must equal."""
    sites = []
    for seed in rrna.SEEDS:
        before, after = (rrna.UP, rrna.DOWN) if seed in rrna.FORWARD else (rrna.DOWN, rrna.UP)
        i = seq.find(seed)
        while i >= 0:
            sites.append((max(0, i - before), min(len(seq), i + after)))
            i = seq.find(seed, i + 1)
    spans = []
    for a, b in sorted(sites):
        if spans and a <= spans[-1][1]:
            spans[-1][1] = max(spans[-1][1], b)
        else:
            spans.append([a, b])
    return spans


class StreamTest(unittest.TestCase):
    def test_stretches_read_alone(self):
        rng = random.Random(2)
        c1, c2 = "".join(rng.choice("ACGTacgt") for _ in range(50)), "".join(rng.choice("ACGT") for _ in range(30))
        with tempfile.TemporaryDirectory() as tmp:
            for path in (os.path.join(tmp, "g.fna"), os.path.join(tmp, "g.fna.gz")):
                text = ">c1 x\n" + "\n".join(c1[i:i + 7] for i in range(0, 50, 7)) + "\n>c2\n" + c2 + "\n"
                with (gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")) as fh:
                    fh.write(text)
                got = placement.read_spans(path, {"c1": [(5, 12), (10, 20), (40, 60)], "c9": [(1, 5)]})
                self.assertEqual(set(got), {"c1"})
                s = got["c1"]
                self.assertEqual((len(s), s[4:20], s[39:50], s[10:15], s[45:]), (50, c1[4:20], c1[39:50], c1[10:15],
                                                                              c1[45:]))
                with self.assertRaises(IndexError):
                    s[0:3]  # not read
                with self.assertRaises(IndexError):
                    s[15:45]  # across what was not read

    def test_rrna_seeds_in_one_pass(self):
        rng = random.Random(4)
        for n in (5, 900, 40_000, 120_000):
            seq = "".join(rng.choice("ACGT") if rng.random() > 0.01 else "N" for _ in range(n))
            for _ in range(3):
                m = rng.choice(rrna.MOTIFS)
                m = m if rng.random() < 0.5 else rrna.revcomp(m)
                if n > 40:
                    i = rng.randrange(n - 30)
                    seq = seq[:i] + m + seq[i + len(m):]
            self.assertEqual(rrna.windows(seq.lower()), old_windows(seq))
            # fed in pieces, batches cut anywhere, a second sequence after it
            old, rrna.SEED_BATCH = rrna.SEED_BATCH, 333
            try:
                scan = rrna.SeedScan()
                scan.start("a")
                for i in range(0, len(seq), 101):
                    scan.feed(seq[i:i + 101].encode())
                scan.start("b")
                scan.feed(rrna.MOTIFS[0].encode())
                found, lengths = scan.finish()
            finally:
                rrna.SEED_BATCH = old
            self.assertEqual((found.get("a", []), lengths["a"], found["b"]),
                             (old_windows(seq), len(seq), old_windows(rrna.MOTIFS[0])))


class FcsTsaTest(unittest.TestCase):
    def test_a_transcriptomes_contig_name_holds_bars(self):
        ranges = [dict(window="88at2759|ENA|GKLH01011938|GKLH01011938.1|101-682", start=10, end=20, action="EXCLUDE",
                       division="prok:bacteria", top_taxon="E. coli")]
        self.assertEqual(genomes.fcs_of_copy(ranges, "88at2759", "ENA|GKLH01011938|GKLH01011938.1", 105, 120),
                         "EXCLUDE:prok:bacteria:E. coli")  # 110-120 on the contig
        self.assertEqual(genomes.fcs_of_copy(ranges, "88at2759", "ENA|GKLH01011938|GKLH01011938.1", 200, 300), "-")
        self.assertEqual(genomes.fcs_of_copy(ranges, "88at2759", "ENA", 105, 120), "-")


class ClaimTest(unittest.TestCase):
    def test_one_build_a_genome(self):
        with tempfile.TemporaryDirectory() as work:
            acc = "GCA_000000001.1"
            with genomes.claimed(work, acc) as held:
                self.assertTrue(held and os.path.isfile(held))
                self.assertIsNone(genomes.claim(work, acc))  # held by a live process
                partial = os.path.join(work, "done", f".otherhost-1.{acc}.partial")
                os.makedirs(partial)
                self.assertEqual(genomes.stale_partials(work), [])  # its genome's claim lives
            self.assertFalse(os.path.isfile(genomes.claim_path(work, acc)))  # released
            self.assertEqual(genomes.stale_partials(work), [partial])
            # another host's claim: held while touched, taken over once stale (its node died)
            path = genomes.claim_path(work, acc)
            with open(path, "w") as fh:
                fh.write("some-other-node\t12345\n")
            self.assertIsNone(genomes.claim(work, acc))
            old = time.time() - genomes.CLAIM_STALE - 10
            os.utime(path, (old, old))
            self.assertEqual(genomes.claim(work, acc), path)
            # this host's claim of a process that ended: taken over
            with open(path, "w") as fh:
                fh.write(f"{__import__('socket').gethostname()}\t999999999\n")
            self.assertFalse(genomes.claim_live(path))


TIARA = """#!/bin/sh
while [ $# -gt 0 ]; do case $1 in -i) i=$2;; -o) o=$2;; esac; shift; done
printf 'sequence_id\\tclass_fst_stage\\tclass_snd_stage\\n' > $o
grep '>' $i | sed 's/>//' | awk '{print $1 "\\teukarya\\tunknown"}' >> $o
"""


class BackfillTest(unittest.TestCase):
    def test_tiara_made_later(self):
        with tempfile.TemporaryDirectory() as tmp:
            work, acc = os.path.join(tmp, "work"), "GCA_000000001.1"
            folder = genomes.done_dir(work, acc)
            os.makedirs(folder)
            row = ["1at2759", "Single", "c1", "1000", "101", "400", "+"] + ["1"] * (len(genomes.LOCI_COLUMNS) - 7)
            row[genomes.TIARA_COLUMN] = "-"
            with open(os.path.join(folder, "loci.tsv"), "w") as fh:
                fh.write("\t".join(genomes.LOCI_COLUMNS) + "\n" + "\t".join(row) + "\n"
                         + "\t".join(["2at2759", "Missing"] + ["-"] * (len(genomes.LOCI_COLUMNS) - 2)) + "\n")
            with gzip.open(os.path.join(folder, "windows.fna.gz"), "wt") as fh:
                fh.write(">1at2759|c1|1-1000\n" + "ACGT" * 250 + "\n")
            with open(os.path.join(folder, "done.json"), "w") as fh:
                json.dump({"accession": acc, "tiara": "failed"}, fh)
            job = dict(accession=acc)
            failing = dict(work=work, genomes=tmp, tiara=fx.tool(tmp, "tiara_fails", "#!/bin/sh\nexit 1\n"), rrna=None)
            self.assertEqual(genomes.backfill_needs(job, failing), ["tiara"])
            self.assertEqual(genomes.backfill((job, failing)), (acc, [], "tiara failed"))
            settings = dict(failing, tiara=fx.tool(tmp, "tiara", TIARA))
            self.assertEqual(genomes.backfill((job, settings)), (acc, ["tiara"], None))
            rows = {r["locus"]: r for r in fx.read_tsv(os.path.join(folder, "loci.tsv"))}
            self.assertEqual((rows["1at2759"]["window_class"], rows["2at2759"]["window_class"]), ("eukarya", "-"))
            self.assertEqual(genomes.done_summary(work, acc)["tiara"], "done")
            self.assertEqual(genomes.backfill_needs(job, settings), [])
            self.assertEqual(genomes.backfill_needs(job, dict(settings, tiara=None)), [])  # none asked for


CMALIGN = """#!/bin/sh
# cmalign's stand-in: each base on the match state of its own position; fails on a sequence holding NNNN (as on one
# over its matrix limit)
while [ $# -gt 0 ]; do case $1 in -o) o=$2;; esac; last=$1; shift; done
if grep -q NNNN $last; then exit 1; fi
awk 'BEGIN{n=0} /^>/{name=substr($1,2); next} {seq[name]=seq[name] $0; if (length(seq[name])>n) n=length(seq[name]); order[++k]=name}
     END{print "# STOCKHOLM 1.0"; for (i=1;i<=k;i++) {s=seq[order[i]]; while (length(s)<n) s=s "-"; print order[i] " " s}
         rf=""; for (i=0;i<n;i++) rf=rf "x"; print "#=GC RF " rf; print "//"}' $last > $o
"""


class CmalignCacheTest(unittest.TestCase):
    def test_kept_as_it_goes_failures_tried_again(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = os.path.join(tmp, "work")
            os.makedirs(work)
            cm = os.path.join(tmp, "ssu.cm")
            with open(cm, "w") as fh:
                fh.write("stand-in\n")
            ref = regions.reference()
            good, bad = ref[:1700], ref[:800] + "NNNN" + ref[804:1700]
            records = [("r1", good, dict(kind="18S")), ("r2", bad, dict(kind="18S"))]
            tool = fx.tool(tmp, "cmalign", CMALIGN)
            rows = regions.run(records, cm, tool, 1, work, log=lambda *_: None)
            self.assertEqual({r["record"] for r in rows}, {"r1"})
            cache = os.path.join(work, "cmalign", f"ssu.{regions.util.md5(cm)[:10]}.tsv")
            with open(cache) as fh:
                self.assertEqual(len(fh.readlines()), 1)  # r2's failure not cached
            fixed = fx.tool(tmp, "cmalign2", CMALIGN.replace("if grep -q NNNN $last; then exit 1; fi", ""))
            rows = regions.run(records, cm, fixed, 1, work, log=lambda *_: None)
            self.assertEqual({r["record"] for r in rows}, {"r1", "r2"})  # aligned now
            with open(cache) as fh:
                self.assertEqual(len(fh.readlines()), 2)


class SmallTest(unittest.TestCase):
    def test_a_genomes_features_in_the_packages_order(self):
        with tempfile.TemporaryDirectory() as work:
            folder = genomes.done_dir(work, "G1")
            os.makedirs(folder)
            with open(os.path.join(folder, "features.tsv"), "w") as fh:
                fh.write("locus\tkind\tstart\tend\n2at2759\texon\t0\t30\n1at2759\tintron\t10\t20\n1at2759\texon\t0\t10\n")
            got = package.genome_features(work, "G1", {"1at2759": "A" * 5 + "N" * 12 + "A" * 13, "2at2759": "C" * 30})
            self.assertEqual([(f["locus"], f["kind"], f["start"], f["end"]) for f in got],
                             [("1at2759", "exon", 0, 10), ("1at2759", "gap", 5, 17), ("1at2759", "intron", 10, 20),
                              ("2at2759", "exon", 0, 30)])

    def test_organelle_molecules_the_scan_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, organelle.SAVED)
            with gzip.open(path, "wt") as fh:
                fh.write(">m1 mitochondrion\nacgt\nACGT\n>p1 plastid\nGG\n")
            self.assertEqual(organelle.saved_molecules(path), [("m1", "mitochondrion", "ACGTACGT"), ("p1", "plastid",
                                                                                                    "GG")])

    def test_none_noted_for_every_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            job = dict(accession="GCA_1.1", url="https://x/GCA_1.1_a_genomic.fna.gz", annotated=None)
            self.assertFalse(annotation.none_noted(tmp, job))
            annotation.note_none(tmp, "GCA_1.1", ["ncbi"])
            self.assertTrue(annotation.none_noted(tmp, job))
            # a source named since (an index, a VEuPathDB key): asked again
            self.assertFalse(annotation.none_noted(tmp, dict(job, annotations=[("ensembl", "https://e/x.gff3.gz")])))

    def test_barcode_gap_without_vsearch_is_not_measured(self):
        seqs = {"A1": {"18S": "ACGT" * 50}, "A2": {"18S": "ACGT" * 49 + "TTTT"}}
        rows = intraspecific.barcode_gap({1: "A1"}, {"A1": seqs["A1"]}, {"A2": seqs["A2"]}, {"A1": 1, "A2": 1},
                                         {(1, "18S"): 90.0}, {}, None, None)
        self.assertEqual(rows, [])  # until 2026-11-08: every pair read as 50% and the gap called unresolved


if __name__ == "__main__":
    unittest.main()
