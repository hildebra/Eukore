#!/usr/bin/env python3
"""Tests of Eukore's build: an rDNA unit uncut, with its 5.8S and an inline 5S (one block) and not a 5S elsewhere, a
tandem array's spacer 5S one unit's, a 28S cut by an intron one; a published read's regions on its record; a
EUKARYOME read's NCBI taxid by its species name, else its genus's; the record checks; ITSx's failures seen; the
build's resume (a cut cache line, a killed run's old package) and a genome's unexpected error its own failure.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi as ncbi_mod, rrna  # noqa: E402
import fixtures as fx  # noqa: E402

tool = fx.tool


class UnitTest(unittest.TestCase):
    def test_a_unit_whole_with_its_inline_5s(self):
        rng = random.Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(30000))
        feats = {"18S_rRNA": [(1001, 2800, "+", False)], "5_8S_rRNA": [(3201, 3350, "+", False)],
                 "28S_rRNA": [(3601, 12100, "+", False)],
                 "5S_rRNA": [(12500, 12620, "-", False), (25000, 25120, "+", False)]}  # inline, and far away
        unit = rrna.units_of(seq, feats)[0]
        self.assertEqual(len(unit["unit"]), 12620 - 1001 + 1)  # 18S to the inline 5S, uncut (over 8,191)
        names = [r[0] for r in unit["regions"]]
        self.assertEqual(names, ["18S", "ITS1", "5.8S", "ITS2", "28S", "5S"])
        self.assertEqual(dict((r[0], r[1:]) for r in unit["regions"])["5S"], (12500 - 1001, 12620 - 1001 + 1))
        self.assertTrue(unit["lsu_complete"])
        # the reverse strand: the same unit, its regions on its own strand
        rc = rrna.revcomp(seq)
        n = len(seq)
        flip = {k: [(n - b + 1, n - a + 1, "-" if s == "+" else "+", p) for a, b, s, p in v] for k, v in feats.items()}
        flip["28S_rRNA"] = [(a, b, s, True) for a, b, s, _p in flip["28S_rRNA"]]  # a partial 28S hit
        back = rrna.units_of(rc, flip)[0]
        self.assertEqual(back["unit"], unit["unit"])
        self.assertEqual(back["regions"], unit["regions"])
        self.assertFalse(back["lsu_complete"])
        # no 5.8S, no unit: the 18S span alone
        alone = rrna.units_of(seq, {"18S_rRNA": [(1001, 2800, "+", False)]})[0]
        self.assertIsNone(alone["unit"])
        # barrnap's 28S hit reaching back over the 5.8S (Drosophila's units, 2026-11-05 trial): the 28S starts at the
        # 5.8S's end, no ITS2 (ITSx gives it where installed), the regions in order on either strand
        over = dict(feats, **{"28S_rRNA": [(3300, 12100, "+", False)]})
        for strand_seq, f in ((seq, over), (rc, {k: [(n - b + 1, n - a + 1, "-" if s == "+" else "+", p)
                                                      for a, b, s, p in v] for k, v in over.items()})):
            regions = rrna.units_of(strand_seq, f)[0]["regions"]
            self.assertEqual([r[0] for r in regions], ["18S", "ITS1", "5.8S", "28S", "5S"])
            self.assertTrue(all(a < b for _n, a, b in regions))
            self.assertTrue(all(x[2] <= y[1] for x, y in zip(regions, regions[1:])))

    def test_a_tandem_arrays_spacer_5s_is_one_units(self):
        # two units 9 kb apart, a 5S in the spacer: within 3 kb of the first's 28S end and of the second's 18S start;
        # it is the first's (the 28S it follows), the second has none (until 2026-11-06 both took it)
        rng = random.Random(4)
        seq = "".join(rng.choice("ACGT") for _ in range(40000))
        feats = {"18S_rRNA": [(1001, 2800, "+", False), (15001, 16800, "+", False)],
                 "5_8S_rRNA": [(3201, 3350, "+", False), (17201, 17350, "+", False)],
                 "28S_rRNA": [(3601, 12100, "+", False), (17601, 26100, "+", False)],
                 "5S_rRNA": [(13000, 13120, "+", False)]}
        first, second = rrna.units_of(seq, feats)
        self.assertEqual([r[0] for r in first["regions"]][-1], "5S")
        self.assertNotIn("5S", [r[0] for r in second["regions"]])
        self.assertEqual(first["unit_span"], (1001, 13120))

    def test_a_28s_cut_by_an_intron_is_one(self):
        # barrnap's 28S in two partial pieces (a group I intron, Tetrahymena's): the unit reaches the second's end,
        # whole where their shares make 80% (until 2026-11-06 it ended at the first piece)
        rng = random.Random(5)
        seq = "".join(rng.choice("ACGT") for _ in range(20000))
        feats = {"18S_rRNA": [(1001, 2800, "+", False, 1.0)], "5_8S_rRNA": [(3201, 3350, "+", False, 1.0)],
                 "28S_rRNA": [(3601, 6000, "+", True, 0.45), (6400, 9000, "+", True, 0.5)]}
        unit = rrna.units_of(seq, feats)[0]
        self.assertEqual((unit["unit_span"], unit["lsu_complete"]), ((1001, 9000), True))
        short = dict(feats, **{"28S_rRNA": [(3601, 6000, "+", True, 0.45)]})
        self.assertFalse(rrna.units_of(seq, short)[0]["lsu_complete"])
        gff = ("w\tbarrnap:0.9\trRNA\t1\t900\t1e-9\t+\t.\tName=28S_rRNA;product=28S ribosomal RNA (partial);"
               "note=aligned only 39 percent of the 28S ribosomal RNA\n")
        self.assertEqual(rrna.parse_gff(gff)["w"]["28S_rRNA"], [(1, 900, "+", True, 0.39)])

    def test_a_published_reads_regions_on_its_record(self):
        # a read whose unit starts at its base 150 (+), or lies on its reverse strand: the regions on the record, the
        # read itself or its reverse complement (until 2026-11-06 the unit's, 149 bases early or mirrored)
        unit = dict(unit_span=(150, 1000), span=(150, 900, "+"), regions=[("18S", 0, 751), ("ITS1", 751, 851)])
        self.assertEqual(rrna.on_record(unit, 1200), [("18S", 149, 900), ("ITS1", 900, 1000)])
        rev = dict(unit_span=(201, 1051), span=(300, 1051, "-"), regions=[("18S", 0, 752)])
        self.assertEqual(rrna.on_record(rev, 1200), [("18S", 149, 901)])  # the record: the read's reverse complement
        rng = random.Random(6)
        read = "".join(rng.choice("ACGT") for _ in range(1200))
        record = rrna.revcomp(read)
        self.assertEqual(record[149:901], rrna.revcomp(read[300 - 1:1051]))  # the 18S, read on its strand


class HomonymTest(unittest.TestCase):
    def test_a_genus_of_two_kingdoms(self):
        from eukore import labels  # noqa: PLC0415
        # EUKARYOME's classification names Drosophila twice (a fly, a fungus); a long read of the fly's takes the
        # fly's lineage by its kingdom-specific name, and a plain-name read none of either (its own k__ .. g__)
        rows = [("superkingdom",) * 10,
                ("Eukaryota", "Fungi", "", "Ascomycota", "", "Sordariomycetes", "Hypocreales", "Fam", "Drosophila",
                 "Drosophila(Fungi)"),
                ("Eukaryota", "Metazoa", "", "Arthropoda", "", "Insecta", "Diptera", "Drosophilidae", "Drosophila",
                 "Drosophila(Metazoa)")]
        ref = labels.read_eukaryome(rows, [("EUK1", "k__Metazoa;p__Arthropoda;g__Drosophila(Metazoa);s__melanogaster"),
                                           ("EUK2", "k__Metazoa;p__Arthropoda;g__Drosophila;s__simulans")])
        fly = next(iter(ref.species["drosophila melanogaster"]))
        self.assertEqual(fly[1], "Metazoa")
        self.assertEqual(next(iter(ref.species["drosophila simulans"]))[1], "Metazoa")  # its own k__
        # a genus match of another kingdom's lineage alone: none
        self.assertTrue(labels.other_kingdom(("Eukaryota", "Archaeplastida", "Streptophyta", "Morus"),
                                             {"eukaryota", "metazoa", "chordata", "morus"}))
        self.assertFalse(labels.other_kingdom(("Eukaryota", "Obazoa", "Metazoa", "Morus"),
                                              {"eukaryota", "metazoa", "chordata", "morus"}))


class ResumeTest(unittest.TestCase):
    def test_a_cut_cache_line_is_annotated_again(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            cache = os.path.join(tmp, "regions.tsv")
            with open(cache, "w", newline="\n") as fh:
                fh.write("EUK_1\t18S:0-1800\tt\t+\nEUK_2\t18S:0-1800,ITS1:1800-2\tf\t-\tabc123abc123\nEUK_3\t18S:1-")
            # a line before 2026-11-08 has no digest of its read; a later one has (a changed read: annotated again)
            self.assertEqual(eukore_build.read_regions(cache), {"EUK_1": ("18S:0-1800", "t", "+", None),
                                                                "EUK_2": ("18S:0-1800,ITS1:1800-2", "f", "-",
                                                                          "abc123abc123")})
            self.assertFalse(eukore_build.util.ends_in_newline(cache))
            self.assertEqual(eukore_build.read_digest("ACGT"), "f1f8f4bf413b")

    def test_no_genome_deleted_without_keep_gb(self):
        # since 2026-11-08 deleting genome files is asked for (--keep_gb): a superseded assembly's kept without it
        import types  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        for keep_gb, calls in ((None, 0), (1000.0, 1)):
            opts = types.SimpleNamespace(keep_gb=keep_gb, genomes="OUT/genomes", work="WORK")
            with mock.patch.object(eukore_build.genomes, "drop_superseded", return_value=[]) as drop:
                eukore_build.drop_superseded(opts, {})
            self.assertEqual(drop.call_count, calls)

    def test_a_genomes_unexpected_error_is_its_failure(self):
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with mock.patch.object(eukore_build.genomes, "run", side_effect=UnicodeDecodeError("utf-8", b"\xff", 0, 1,
                                                                                         "invalid start byte")):
            acc, summary, why = eukore_build.run_one(({"accession": "GCA_1.1"}, {}))
        self.assertEqual((acc, summary), ("GCA_1.1", None))
        self.assertTrue(why.startswith("UnicodeDecodeError"))
        # a step's failure (_run's four values, its gene models' files the fourth): three, as a done genome's (until
        # 2026-11-09 the four, and the process stage stopped at the first genome that failed: too many values)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "g.fna.gz")
            open(path, "w").close()
            failed = ("GCA_1.1", None, "locus calling failed: x", (None, None))
            with mock.patch.object(eukore_build.genomes, "_run", return_value=failed):
                got = eukore_build.genomes.run(({"accession": "GCA_1.1", "path": path},
                                                {"work": os.path.join(d, "WORK"), "genomes": d}))
            self.assertEqual(got, failed[:3])
            self.assertFalse(os.path.isdir(eukore_build.genomes.done_dir(os.path.join(d, "WORK"), "GCA_1.1")))
        # a batch of published reads barrnap fails on, a GFF comparison's unexpected error: theirs, not the stage's
        with mock.patch.object(eukore_build.rrna, "published_regions", side_effect=RuntimeError("barrnap failed")):
            self.assertEqual(eukore_build._annotate_batch([("EUK_1", "ACGT", "all")], "barrnap"),
                             ({}, "RuntimeError: barrnap failed"))
        with mock.patch.object(eukore_build.annotation, "annotate_genome", side_effect=ValueError("bad GFF")):
            self.assertEqual(eukore_build.annotate_one("w", "g", False, None, {"accession": "GCA_1.1"}),
                             ("GCA_1.1", None, "ValueError: bad GFF"))

    def test_the_banner_says_the_version_and_the_run(self):
        # since 2026-11-08 a run opens with Eukore's version and what it will do, before any work; a key never said
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            args = ["--build_dir", os.path.join(tmp, "work"), "--out", os.path.join(tmp, "out"), "--species", "4751,5794",
                    "--groups", "fungi,protozoa", "--memory_gb", "100", "-t", "16", "--no_tsa", "--veupathdb_key",
                    "s3cret", "--veupathdb=s3cret2", "--version", "2026.11"]
            with mock.patch.object(eukore_build, "fetch_sources", side_effect=RuntimeError("fetched")), \
                    contextlib.redirect_stdout(io.StringIO()) as out, self.assertRaises(RuntimeError):
                eukore_build.main(args)
        text = out.getvalue()
        self.assertTrue(text.startswith("=" * 100 + f"\n Eukore {eukore_build.__version__}: "))
        for said in ("stage     all: fetch -> process -> package", "groups    fungi, protozoa",
                     "species   a trial: the clades of 2 taxids (4751, 5794)",
                     "threads   16; 4 genomes at a time, 4 threads each; 8 downloads at a time",
                     "memory    100.0 GB for the genomes' runs (--memory_gb)", "genomes   every one kept (no --keep_gb)",
                     "off       transcriptomes", "OUT/Eukore_2026.11/", "--veupathdb_key HIDDEN",
                     "--veupathdb=HIDDEN"):
            self.assertIn(said, text)
        self.assertNotIn("s3cret", text)
        self.assertTrue(text.isascii())  # a node's log in any locale
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()) as out, \
                mock.patch.object(eukore_build, "fetch_sources", side_effect=RuntimeError("fetched")), \
                mock.patch.object(eukore_build.ncbi_mod, "Ncbi", side_effect=RuntimeError("planned")), \
                self.assertRaises(RuntimeError):
            eukore_build.main(["--build_dir", tmp, "--out", tmp, "--stage", "package", "--keep_gb", "500"])
        text = out.getvalue()
        self.assertIn("stage     package\n   package ", text)
        self.assertIn("genomes   cut back to 500 GB (--keep_gb)", text)
        self.assertNotIn("memory ", text)  # the process stage's alone
        self.assertNotIn("downloads", text)  # the fetch stage's alone

    def test_a_scratch_that_cannot_be_made(self):
        # a SLURM job's $SLURM_LOCAL_SCRATCH named, its folder not made (no local disk asked for): the fetch stage
        # puts its .part files beside the genomes, the process stage stops at once (not after the planning's minutes)
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(eukore_build.scratch_problem(os.path.join(tmp, "job", "x")))
            self.assertTrue(os.path.isdir(os.path.join(tmp, "job", "x")))
            locked = os.path.join(tmp, "locked")
            os.makedirs(locked)
            os.chmod(locked, 0o500)
            try:
                if os.access(locked, os.W_OK):
                    self.skipTest("root writes anywhere")
                self.assertIn("Permission denied", eukore_build.scratch_problem(os.path.join(locked, "job", "x")))
                self.assertIn("not writable", eukore_build.scratch_problem(locked))
                args = ["--build_dir", os.path.join(tmp, "work"), "--out", os.path.join(tmp, "out"),
                        "--local_tmp", os.path.join(locked, "job", "x")]
                with contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit):
                    eukore_build.main(args + ["--stage", "process"])
                self.assertIn("--local_tmp", err.getvalue())
                self.assertFalse(os.path.exists(os.path.join(tmp, "work")))  # stopped before any work
                with mock.patch.object(eukore_build, "fetch_sources", side_effect=RuntimeError("fetched")), \
                        contextlib.redirect_stdout(io.StringIO()) as out, self.assertRaises(RuntimeError):
                    eukore_build.main(args + ["--stage", "fetch"])
                self.assertIn("WARNING: --local_tmp", out.getvalue())
            finally:
                os.chmod(locked, 0o700)

    def test_genomes_only_runs_the_genomes_alone(self):
        # several process runs on one WORK (a SLURM array, since 2026-11-09): with --genomes_only each runs its share
        # of the genomes and stops, the passes over all genomes left to one run without it once all have ended; the
        # passes that take a genome at a time under its claim (the backfill, the metadata) are every run's (since
        # 2026-11-10: the last run's alone, on one node)
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        passes = ("backfill_pass", "intraspecific_pass", "describe_pass", "drop_superseded", "organelle_genes",
                  "annotate_pass")
        with tempfile.TemporaryDirectory() as tmp:
            args = ["--build_dir", os.path.join(tmp, "work"), "--out", os.path.join(tmp, "out"), "--stage", "process"]
            with contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit):
                eukore_build.main(args[:-1] + ["package", "--genomes_only"])  # `all` would package without them
            self.assertIn("--genomes_only: with --stage process alone", err.getvalue())
            for only in (True, False):
                with contextlib.ExitStack() as stack:
                    out = stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
                    stack.enter_context(mock.patch.object(eukore_build.ncbi_mod, "Ncbi"))
                    for name, value in (("inside_of", None), ("plan", ({}, {})), ("extras_of", {}),
                                        ("settings_of", {}), ("run_genomes", ([], []))):
                        stack.enter_context(mock.patch.object(eukore_build, name, return_value=value))
                    stack.enter_context(mock.patch.object(eukore_build.shutil, "which", side_effect=lambda t: "/x/" + t))
                    stack.enter_context(mock.patch.object(eukore_build.rrna, "check_itsx", return_value=None))
                    called = {n: stack.enter_context(mock.patch.object(eukore_build, n)) for n in passes}
                    self.assertEqual(eukore_build.main(args + (["--genomes_only"] if only else [])), 0)
                self.assertEqual({n for n, m in called.items() if m.called},
                                 {"backfill_pass", "describe_pass"} if only else set(passes))
                self.assertEqual("left to a process run without it" in out.getvalue(), only)

    def test_the_folders_old_flags_still_taken(self):
        # --build_dir and --local_tmp since 2026-11-09; a running build's scripts' --work and --scratch still mean what
        # they did (never --scratch the build's folder: a node's local disk, emptied after its job)
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        self.assertEqual(eukore_build.renamed(["--work", "W", "--scratch=/l", "--out", "O", "--work_x"]),
                         ["--build_dir", "W", "--local_tmp=/l", "--out", "O", "--work_x"])
        with tempfile.TemporaryDirectory() as tmp:
            work, local = os.path.join(tmp, "work"), os.path.join(tmp, "local")
            with mock.patch.object(eukore_build, "fetch_sources", side_effect=RuntimeError("fetched")), \
                    contextlib.redirect_stdout(io.StringIO()) as said, self.assertRaises(RuntimeError):
                eukore_build.main(["--work", work, "--out", os.path.join(tmp, "out"), "--scratch", local])
            self.assertIn(f"build_dir {work}\n", said.getvalue())
            self.assertIn(f"local_tmp {local}\n", said.getvalue())
            self.assertTrue(os.path.isdir(os.path.join(work)))

    def test_the_genomes_folder_apart_from_out(self):
        # --genomes (since 2026-11-09): the downloads in a folder of their own, not OUT/genomes; the banner says it
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            out, apart = os.path.join(tmp, "out"), os.path.join(tmp, "downloads")
            for given in ([], ["--genomes", apart]):
                with mock.patch.object(eukore_build, "fetch_sources", side_effect=RuntimeError("fetched")), \
                        contextlib.redirect_stdout(io.StringIO()) as said, self.assertRaises(RuntimeError):
                    eukore_build.main(["--build_dir", os.path.join(tmp, "work"), "--out", out, "--stage", "fetch"] + given)
                self.assertEqual(f"the genomes in {apart} (--genomes)" in said.getvalue(), bool(given))
                self.assertEqual("into --genomes" in said.getvalue(), bool(given))
            self.assertTrue(os.path.isdir(apart))
            self.assertEqual(os.listdir(out), ["genomes"])  # the run without --genomes'

    def test_a_package_written_over_a_killed_runs_old(self):
        from eukore import package  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "Eukore_2026.11")
            for folder in (out, out + ".old"):  # a package, and a killed run's .old left beside it
                os.makedirs(os.path.join(folder, "busco"))
                with open(os.path.join(folder, "busco", "x.tsv"), "w") as fh:
                    fh.write("x\n")
            part = os.path.join(tmp, "new")
            os.makedirs(part)
            with open(os.path.join(part, "manifest.json"), "w") as fh:
                fh.write("{}\n")
            self.assertEqual(package.put_in_place(part, out), out)
            self.assertEqual(os.listdir(out), ["manifest.json"])
            self.assertEqual(sorted(os.listdir(tmp)), ["Eukore_2026.11"])

    def test_only_ended_runs_partials_are_leftovers(self):
        # a package stage going on keeps its unfinished package: this host's of a live process, another host's
        # touched lately; a killed run's (its process gone, or another host's untouched for days) is a leftover
        import socket  # noqa: PLC0415
        import time  # noqa: PLC0415
        from eukore import util  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            host = socket.gethostname()
            names = {"live": f".{host}-{os.getpid()}.Eukore_2026.11.partial",
                     "dead": f".{host}-999999999.Eukore_2026.11.partial",
                     "other": ".node-b-12.Eukore_2026.11.partial", "old": ".node-c-12.Eukore_2026.11.partial",
                     "another": f".{host}-999999999.Eukore_2026.12.partial"}
            for name in names.values():
                os.makedirs(os.path.join(tmp, name))
            old = time.time() - util.PARTIAL_STALE - 60
            os.utime(os.path.join(tmp, names["old"]), (old, old))
            self.assertEqual(sorted(util.dead_partials(tmp, "Eukore_2026.11")),
                             sorted(os.path.join(tmp, names[k]) for k in ("dead", "old")))

    def test_itsx_failures_are_seen(self):
        # ITSx with --cpu over 1 exits 0 when HMMER fails on a missing profile set: the build sees it, before any work
        with tempfile.TemporaryDirectory() as tmp:
            broken = tool(tmp, "ITSx_broken", "#!/bin/sh\necho 'Error: File existence/permissions problem in trying "
                                              "to open HMM file F.hmm' >&2\necho 'ERROR! HMMER search terminated "
                                              "unexpectedly.' >&2\nexit 0\n")
            failing = tool(tmp, "ITSx_failing", "#!/bin/sh\nexit 1\n")
            good = tool(tmp, "ITSx_good", "#!/bin/sh\necho 'Preparing HMM database (should be quick)...' >&2\n"
                                          "exit 0\n")
            self.assertIn("F.hmm", rrna.check_itsx(broken))
            self.assertIn("exit 1", rrna.check_itsx(failing))
            self.assertIsNone(rrna.check_itsx(good))
            with self.assertRaises(rrna.ItsxError):
                rrna.run_itsx(["ACGT" * 50], broken)

    def test_itsx_guessed_borders_are_not_taken(self):
        # ITSx without the 5.8S's end still writes an ITS2, 45 bases after the 5.8S start (Arabidopsis with set T,
        # 2026-11-06): ITS1 and the 28S start are ITSx's, the 5.8S's end barrnap's
        line = "u0\t3755 bp.\tSSU: 1-300\tITS1: 301-590\t5.8S: No end\tITS2: 636-1206\tLSU: 1207-3755\t"
        with tempfile.TemporaryDirectory() as tmp:
            itsx = tool(tmp, "ITSx_noend", "#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = -o ] && o=$2; shift; done\n"
                                           f"printf '{line}\\n' > $o.positions.txt\nexit 0\n")
            found = rrna.run_itsx(["ACGT" * 50], itsx)["ACGT" * 50]
        self.assertEqual(found["5.8S problem"], "No end")
        self.assertEqual((found["ITS1"], found["ITS2"]), ((301, 590), (636, 1206)))

        def item():  # barrnap's regions on a unit, the ITSx stretch from 1500
            return dict(unit="A" * 5000, its_at=1500, regions=[("18S", 0, 1800), ("ITS1", 1800, 2100),
                                                               ("5.8S", 2100, 2260), ("ITS2", 2260, 2500),
                                                               ("28S", 2500, 5000)])
        got = {r[0]: r[1:] for r in rrna.refine(item(), found)["regions"]}
        self.assertEqual(got, {"18S": (0, 1800), "ITS1": (1800, 2090), "5.8S": (2090, 2260), "ITS2": (2260, 2706),
                               "28S": (2706, 5000)})
        whole = {k: v for k, v in found.items() if k != "5.8S problem"}  # what was taken before: a 45-base 5.8S
        self.assertEqual({r[0]: r[1:] for r in rrna.refine(item(), whole)["regions"]}["5.8S"], (2090, 2135))
        neither = dict(found, **{"5.8S problem": "Not found"})
        self.assertEqual({r[0]: r[1:] for r in rrna.refine(item(), neither)["regions"]}["5.8S"], (2100, 2260))
        start = dict(found, **{"5.8S problem": "No start"})
        self.assertEqual({r[0]: r[1:] for r in rrna.refine(item(), start)["regions"]}["5.8S"], (2100, 2135))

    def test_itsx_searches_the_lineages_sets(self):
        # each record's ITSx profile sets by its NCBI lineage (fungi F, animals M), all sets for one in none of them;
        # one ITSx run per distinct set
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            self.assertEqual([rrna.itsx_groups(ncbi, t) for t in (559292, "6253", 952, None)],
                             ["F", "M", "all", "all"])
            log = os.path.join(tmp, "itsx_calls.txt")
            logging = tool(tmp, "ITSx_log", f"#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = -t ] && echo \"$2\" >> {log}; "
                                            "shift; done\nexit 0\n")
            got = rrna.run_itsx(["AAAA", "CCCC", "GGGG", "AAAA"], logging,
                                groups={"AAAA": "F", "CCCC": "M", "GGGG": "F"})
            self.assertEqual(set(got), {"AAAA", "CCCC", "GGGG"})
            with open(log) as fh:
                self.assertEqual(sorted(fh.read().split()), ["F", "M"])

    def test_labels_and_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        index = rrna.name_index(ncbi)
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__Saccharomyces_cerevisiae", index, ncbi),
                         (4932, "eukaryome_name"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__Saccharomyces_novus", index, ncbi),
                         (4930, "eukaryome_genus"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Nogenus;s__Nogenus_x", index, ncbi), (None, None))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Kluyveromyces;s__Candida_sphaerica", index, ncbi),
                         (28985, "eukaryome_name"))  # an NCBI synonym
        # EUKARYOME 2.1 writes the epithet alone in s__ (g__Geotrupes;s__spiniger): the genus and it are the name
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__cerevisiae", index, ncbi),
                         (4932, "eukaryome_name"))
        self.assertEqual(rrna.eukaryome_taxid("k__Fungi;g__Saccharomyces;s__unclassified", index, ncbi),
                         (4930, "eukaryome_genus"))
        self.assertEqual(rrna.species_name("Drosophila(Metazoa)", "melanogaster"), "Drosophila melanogaster")
        self.assertEqual((rrna.homopolymer("ACGGGGTA"), round(rrna.ambiguous("ACGN"), 2)), (4, 0.25))
        a = "".join(random.Random(2).choice("ACGT") for _ in range(1000))
        b = a[:500] + ("A" if a[500] != "A" else "C") + a[501:]  # 99.9% identical: one type
        c = "".join(random.Random(3).choice("ACGT") for _ in range(1000))
        self.assertEqual(rrna.types_of([a, a, b, c]), {a: 1, b: 1, c: 2})


class PartialTest(unittest.TestCase):
    def test_species_clades_limit_the_package(self):
        # --species (a trial build): the genus Saccharomyces and Candida albicans, as taxids or a file; a strain's
        # taxid is inside its species' clade
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            listed = os.path.join(tmp, "clades.txt")
            with open(listed, "w") as fh:
                fh.write("# trial\n4930\t# Saccharomyces\n5476\n")
            self.assertEqual(eukore_build.clades(listed), {4930, 5476})
            self.assertEqual(eukore_build.clades("4930,5476"), {4930, 5476})
            self.assertIsNone(eukore_build.clades(None))
            inside = eukore_build.inside_of(ncbi, {4930, 5476})
            self.assertEqual([inside(t) for t in (559292, 4932, "5476", 28985, None, "-")],
                             [True, True, True, False, False, False])
            self.assertTrue(eukore_build.inside_of(ncbi, None)(28985))


class RerunTest(unittest.TestCase):
    def test_which_genomes_are_run_again(self):
        # genomes_to_run (2026-11-08): each condition alone, on done folders as their runs left them
        import json  # noqa: PLC0415
        import types  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import genomes  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            opts = types.SimpleNamespace(work=os.path.join(tmp, "work"), genomes=os.path.join(tmp, "genomes"))
            os.makedirs(opts.genomes)
            loci = ["1at2759", "2at2759"]
            settings = dict(loci=loci, staged_only=True, odb_own={7: [7]}, rescue=True, esl_translate="esl-translate")
            ncbi_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/x/GCA_6.1_asm_genomic.fna.gz"

            def done(acc, here=True, **summary):
                os.makedirs(genomes.done_dir(opts.work, acc))
                with open(os.path.join(genomes.done_dir(opts.work, acc), "done.json"), "w") as fh:
                    json.dump(dict(dict(searched=loci, rescue=True, annotation="-", annotation_none=[],
                                        genome_record=genomes.GENOME_RECORD), **summary), fh)
                if here:
                    open(genomes.genome_file(opts.genomes, acc), "w").close()
            jobs = {a: dict(accession=a, species_taxid=1, url="") for a in (
                "new", "done", "short", "held_out", "held_out_gone", "rescue", "annotation", "old_record",
                "old_record_gone", "reversioned")}
            done("done")
            # since 2026-11-09: a genome of an older genome record is run again (its genome here), unless stamped
            done("old_record", genome_record=1)
            done("old_record_gone", here=False, genome_record=1)  # the fetch stage fetches it again
            done("reversioned", genome_record=1)
            self.assertEqual(genomes.reversion(opts.work, "reversioned"), 1)
            self.assertIsNone(genomes.reversion(opts.work, "reversioned"))  # of this record now: left alone
            self.assertEqual(genomes.done_summary(opts.work, "reversioned")["reversioned"][0]["record"], 1)
            with open(os.path.join(genomes.done_dir(opts.work, "old_record"), "done.json"), "w") as fh:
                json.dump(dict(searched=loci, rescue=True, annotation="-", annotation_none=[]), fh)  # no record: 1
            self.assertTrue(genomes.record_stale(opts.work, "old_record"))
            done("short", searched=loci[:1])  # searched for fewer loci than called now
            jobs["held_out"]["species_taxid"] = jobs["held_out_gone"]["species_taxid"] = 7  # a dataset species
            done("held_out")
            done("held_out_gone", here=False)  # its genome no longer here: not fetched again for the benchmark
            done("rescue", rescue=False)  # run without the rescue pass, which this run makes
            jobs["annotation"]["url"] = ncbi_url  # NCBI's gene models beside it now: run again for them alone
            done("annotation")
            open(os.path.join(opts.genomes, "annotation.gff.gz"), "w").close()
            todo, again = eukore_build.genomes_to_run(opts, jobs, settings)
            got = {j["accession"]: j.get("rerun", "whole") for j, _s in todo}
            self.assertEqual(got, {"new": "whole", "short": "whole", "held_out": "whole", "rescue": "whole",
                                   "annotation": "annotation", "old_record": "whole"})
            self.assertEqual(again, 1)
            # each with its done.json as the plan saw it: a genome another build finished meanwhile is left alone
            stamps = {j["accession"]: j["done_stamp"] for j, _s in todo}
            self.assertIsNone(stamps["new"])
            self.assertEqual(stamps["short"], genomes.done_stamp(opts.work, "short"))
            done("new")
            new = next(j for j, _s in todo if j["accession"] == "new")
            run = dict(settings, work=opts.work, genomes=opts.genomes)
            self.assertEqual(genomes.run((new, run)), ("new", None, genomes.DONE_ELSEWHERE))
            self.assertFalse(os.path.isfile(genomes.claim_path(opts.work, "new")))  # its claim released


class Steps:
    """A pass's progress, its reports unsaid."""
    def step(self, *_args, **_kwargs):
        pass

    tick = step


class ChainTest(unittest.TestCase):
    """A cluster's build as one chain after the download (since 2026-11-10): the process runs of a SLURM array
    (--genomes_only), then one run of --stage process,package (tools/eukore.sh)."""

    def test_stages_that_follow_one_another(self):
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        for text, stages in (("all", ("fetch", "process", "package")), ("process", ("process",)),
                             ("process,package", ("process", "package")), ("fetch,process", ("fetch", "process")),
                             ("fetch, process ,package", ("fetch", "process", "package")), ("release", ("release",)),
                             ("fetch,package", None), ("package,process", None), ("process,process", None),
                             ("process,release", None), ("all,package", None), ("", None), ("pack", None)):
            self.assertEqual(eukore_build.stages_of(text), stages, text)
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stderr(io.StringIO()) as err, \
                self.assertRaises(SystemExit):
            eukore_build.main(["--build_dir", tmp, "--out", tmp, "--stage", "fetch,package"])
        self.assertIn("stages that follow one another", err.getvalue())

    def run_main(self, tmp, stage, waiting=(), extra=(), prepare=None):
        """main() of --stage (and the options extra), its plan empty and its passes stand-ins -> (its exit, or
        "package" where it came to the package stage; the stand-ins by name, Ncbi's and the marker set's build
        (prepare: its side effect) among them; its STDOUT)."""
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        passes = ("backfill_pass", "intraspecific_pass", "describe_pass", "drop_superseded", "organelle_genes",
                  "annotate_pass")
        args = ["--build_dir", os.path.join(tmp, "work"), "--out", os.path.join(tmp, "out"), "--stage", stage,
                "--local_tmp", os.path.join(tmp, "local"), *extra]
        with contextlib.ExitStack() as stack:
            out = stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            stack.enter_context(mock.patch.dict(os.environ))  # main sets TMPDIR
            stack.enter_context(mock.patch.object(tempfile, "tempdir", None))
            called = {name: stack.enter_context(mock.patch.object(eukore_build, name, return_value=value))
                      for name, value in (("inside_of", None), ("plan", ({}, {})), ("extras_of", {}),
                                          ("settings_of", {}), ("run_genomes", (list(waiting), [])))}
            called["Ncbi"] = stack.enter_context(mock.patch.object(eukore_build.ncbi_mod, "Ncbi"))
            called["prepare"] = stack.enter_context(mock.patch.object(eukore_build.markers_mod, "prepare",
                                                                      side_effect=prepare))
            stack.enter_context(mock.patch.object(eukore_build.shutil, "which", side_effect=lambda t: "/x/" + t))
            stack.enter_context(mock.patch.object(eukore_build.rrna, "check_itsx", return_value=None))
            called.update({n: stack.enter_context(mock.patch.object(eukore_build, n)) for n in passes})
            # the package stage's first step after its leftovers
            stack.enter_context(mock.patch.object(eukore_build, "called_loci", side_effect=RuntimeError("package")))
            try:
                got = eukore_build.main(args)
            except RuntimeError as e:
                got = str(e)
            tmpdir = os.environ.get("TMPDIR")
        return got, called, out.getvalue(), tmpdir

    def test_process_and_package_in_one_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            got, called, text, tmpdir = self.run_main(tmp, "process,package")
            self.assertEqual(got, "package")  # on to the package stage in the same run
            self.assertIn("stage     process,package: process -> package", text)
            self.assertEqual(called["plan"].call_count, 1)  # one plan for both
            self.assertEqual(called["backfill_pass"].call_count, 1)  # the process stage's, not the package's own too
            self.assertTrue(called["settings_of"].call_args.kwargs["staged_only"])  # no fetch stage: nothing fetched
            self.assertFalse(called["annotate_pass"].call_args.kwargs["may_fetch"])
            # the package stage alone: its own backfill, and --local_tmp its tools' TMPDIR too
            got, called, text, tmpdir = self.run_main(tmp, "package")
            self.assertEqual((got, called["backfill_pass"].call_count, called["run_genomes"].call_count),
                             ("package", 1, 0))
            self.assertEqual(tmpdir, os.path.join(tmp, "local"))
            # genomes not fetched: a run that cannot fetch them packages without them; the process stage alone stops
            waiting = [("GCA_1.1", "not fetched yet")]
            got, _called, text, _tmpdir = self.run_main(tmp, "process,package", waiting)
            self.assertEqual(got, "package")
            self.assertIn("1 genomes not fetched yet, 0 downloads to try again", text)
            self.assertIn("the package is written without them", text)
            got, _called, text, _tmpdir = self.run_main(tmp, "process", waiting)
            self.assertEqual(got, 0)
            self.assertIn("run again (--stage fetch first on a cluster)", text)

    def test_the_process_stage_builds_the_marker_set(self):
        # since 2026-11-11 the fetch stage fetches the set's inputs alone (it built the set, its hmmsearch of ~230,000
        # candidates on the download node): the process stage builds it before its genomes, or says why not
        import argparse  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "work", "sources")
            before = []  # whether the genomes had run when it was built (run_genomes is run_main's stand-in then)
            got, called, text, _tmpdir = self.run_main(
                tmp, "process", prepare=lambda *a: before.append(eukore_build.run_genomes.called))
            self.assertEqual((got, before), (0, [False]))
            self.assertEqual(called["prepare"].call_args.args[:3], (src, "hmmsearch", 8))
            got, called, text, _tmpdir = self.run_main(tmp, "process", prepare=OSError("marker set not built: no"))
            self.assertEqual(got, 0)  # the genomes called all the same
            self.assertIn("WARNING: marker set not built: no: the genomes are called without Eukore's markers", text)
            self.assertEqual(called["run_genomes"].call_count, 1)
            got, called, _text, _tmpdir = self.run_main(tmp, "process", extra=["--no_markers"])
            self.assertFalse(called["prepare"].called)
            with mock.patch.object(eukore_build, "fetch_sources"):
                got, called, _text, _tmpdir = self.run_main(tmp, "fetch")
            self.assertEqual(got, "package")  # stopped at run_main's called_loci
            self.assertFalse(called["prepare"].called)  # the fetch stage builds nothing
        # the fetch stage counts a done genome's markers before their set is built (data/markers.tsv's)
        table = [m["locus"] for m in eukore_build.markers_mod.read_table()]
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(eukore_build, "loci_of", return_value=["B1"]):
            for no_markers, unbuilt, loci in ((False, False, ["B1"]), (False, True, ["B1"] + table),
                                              (True, True, ["B1"])):
                self.assertEqual(eukore_build.called_loci(argparse.Namespace(no_markers=no_markers), tmp, unbuilt),
                                 loci)

    def test_sources_of_their_own(self):
        # --sources (since 2026-11-11): the sources apart from WORK, e.g. kept on a file system that keeps its files
        # while WORK is on scratch: every stage reads them there, a fetch holds them under a lock of their own too
        import contextlib  # noqa: PLC0415
        import io  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "kept", "sources")
            got, called, text, _tmpdir = self.run_main(tmp, "process,package", extra=["--sources", src])
            self.assertEqual(got, "package")
            called["Ncbi"].assert_called_once_with(os.path.join(src, "taxdump.tar.gz"))
            self.assertEqual(called["plan"].call_args.args[1], src)
            self.assertIn(f"; the sources in {src} (--sources)", text)
            seen = []

            def fetched(src_dir, *_args):
                seen.append((src_dir, os.path.isfile(os.path.join(src_dir, ".fetch.lock")),
                             os.path.isfile(os.path.join(tmp, "work", ".fetch.lock"))))
            with mock.patch.object(eukore_build, "fetch_sources", side_effect=fetched):
                got, _called, _text, _tmpdir = self.run_main(tmp, "fetch", extra=["--sources", src])
            self.assertEqual((got, seen), ("package", [(src, True, True)]))  # stopped at run_main's called_loci
            self.assertFalse(os.path.exists(os.path.join(src, ".fetch.lock")))  # let go after the run
            # no sources there: a stage other than fetch says so (a traceback in every task of an array until then)
            with mock.patch.dict(os.environ), mock.patch.object(tempfile, "tempdir", None), \
                    contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as stop:
                eukore_build.main(["--build_dir", os.path.join(tmp, "new"), "--out", os.path.join(tmp, "out"),
                                   "--stage", "process", "--local_tmp", os.path.join(tmp, "local")])
            self.assertIn(f"no sources in {os.path.join(tmp, 'new', 'sources')}", str(stop.exception.code))

    def test_a_run_waits_for_another_builds_genomes(self):
        # the final job after an array: a genome another build's claim holds is waited for, then run here where that
        # build left it (until 2026-11-10 the run went on, and the package was written, without it)
        import contextlib  # noqa: PLC0415
        import types  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import genomes  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        jobs = {a: {"accession": a, "size": 1e6} for a in ("held", "done_there", "mine")}
        for only, rounds in ((False, 2), (True, 1)):
            with tempfile.TemporaryDirectory() as tmp, contextlib.ExitStack() as stack:
                opts = types.SimpleNamespace(work=tmp, genomes=tmp, scratch=None, genomes_only=only, report_minutes=1)
                seen, to_run = [], [[(j, {}) for j in jobs.values()], [(jobs["held"], {})]]

                def planned(_opts, job_map, _settings):
                    seen.append(sorted(job_map))
                    return (to_run[0] if len(job_map) == 3 else [it for it in to_run[1] if it[0]["accession"]
                                                                 in job_map]), 0

                def scheduled(_opts, todo, _sizes, _progress, _note):
                    if len(todo) == 3:
                        yield "held", None, genomes.RUNNING_ELSEWHERE
                        yield "done_there", None, genomes.DONE_ELSEWHERE
                        yield "mine", {}, None
                    else:
                        yield "held", {}, None
                stack.enter_context(mock.patch.object(eukore_build, "say"))
                stack.enter_context(mock.patch.object(eukore_build, "progress_of", return_value=Steps()))
                stack.enter_context(mock.patch.object(eukore_build, "genomes_to_run", side_effect=planned))
                schedule = stack.enter_context(mock.patch.object(eukore_build, "schedule", side_effect=scheduled))
                ended = stack.enter_context(mock.patch.object(eukore_build, "others_end"))
                self.assertEqual(eukore_build.run_genomes(opts, jobs, {"tiara": None}), ([], []))
                self.assertEqual(schedule.call_count, rounds)
                if only:  # one of the array's runs: the other's genome is the other's
                    ended.assert_not_called()
                else:  # waited for the held genome alone (not the one done there), then asked what is left of it
                    ended.assert_called_once_with(opts, ["held"])  # (and of the panel's assemblies: none held)
                    self.assertEqual(seen, [sorted(jobs), ["held"], []])

    def test_failed_genomes_noted_with_their_reasons(self):
        # each failed genome's why in WORK/failed, gone once it is done, and WORK/failed.tsv of those not done (since
        # 2026-11-11: the log named three of them); a full WORK notes none and the run goes on
        import contextlib  # noqa: PLC0415
        import types  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import genomes  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        jobs = {a: {"accession": a, "size": 1e6} for a in ("full", "bad", "ok", "later", "held")}
        nospace = "OSError: [Errno 28] No space left on device: '/w/done/.x.partial/rrna.tsv'"
        with tempfile.TemporaryDirectory() as tmp, contextlib.ExitStack() as stack:
            opts = types.SimpleNamespace(work=tmp, genomes=tmp, scratch=None, genomes_only=True, report_minutes=1)
            self.assertTrue(genomes.note_failed(tmp, "ok", "an earlier run's"))  # done in this run: its note goes

            def scheduled(_opts, _todo, _sizes, _progress, _note):
                yield "full", None, nospace
                yield "bad", None, "locus calling failed:\tminiprot\nkilled"
                yield "ok", {}, None
                yield "later", None, "not fetched yet"  # waiting: no failure
                yield "held", None, genomes.RUNNING_ELSEWHERE
            said = stack.enter_context(mock.patch.object(eukore_build, "say"))
            stack.enter_context(mock.patch.object(eukore_build, "progress_of", return_value=Steps()))
            stack.enter_context(mock.patch.object(eukore_build, "genomes_to_run",
                                                  return_value=([(j, {}) for j in jobs.values()], 0)))
            stack.enter_context(mock.patch.object(eukore_build, "schedule", side_effect=scheduled))
            waiting, failed = eukore_build.run_genomes(opts, jobs, {"tiara": None})
            self.assertEqual([a for a, _w in failed], ["full", "bad"])
            self.assertEqual(sorted(os.listdir(os.path.join(tmp, genomes.FAILED))), ["bad.txt", "full.txt"])
            self.assertFalse(any("could not be noted" in c.args[0] for c in said.call_args_list))  # every note made
            # "full" done since by another build: its note dropped, the table lists "bad" alone, its why on one line
            os.makedirs(os.path.join(tmp, "done", "full"))
            with open(os.path.join(tmp, "done", "full", "done.json"), "w") as fh:
                fh.write("{}")
            rows = genomes.failed_table(tmp)
            self.assertEqual([(r[0], r[-1]) for r in rows], [("bad", "locus calling failed: miniprot killed")])
            self.assertEqual(genomes.read_rows(os.path.join(tmp, genomes.FAILED_TABLE)), rows)
            self.assertEqual(os.listdir(os.path.join(tmp, genomes.FAILED)), ["bad.txt"])
            # reasons counted without their paths
            self.assertEqual(eukore_build.failure_kinds([("a", nospace), ("b", nospace.replace(".x.", ".y.")),
                                                         ("c", "tiara failed")]),
                             "2x OSError: [Errno 28] No space left on device: ...; 1x tiara failed")
            # a WORK that takes no note: the run goes on, and says so once
            with mock.patch.object(genomes, "note_failed", return_value=False):
                eukore_build.run_genomes(opts, jobs, {"tiara": None})
            self.assertIn("3 genomes' failures or ends could not be noted", said.call_args.args[0])
        with tempfile.TemporaryDirectory() as tmp:  # nothing failed, ever: no table
            self.assertEqual(genomes.failed_table(tmp), [])
            self.assertFalse(os.path.exists(os.path.join(tmp, genomes.FAILED_TABLE)))

    def test_another_builds_claims_end(self):
        import threading  # noqa: PLC0415
        import time  # noqa: PLC0415
        import types  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import genomes  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(eukore_build, "say") as said:
            opts = types.SimpleNamespace(work=tmp, report_minutes=0)
            os.makedirs(os.path.join(tmp, "done"))

            def claim(acc, age=0):
                path = genomes.claim_path(tmp, acc)
                with open(path, "w") as fh:
                    fh.write("another-node\t4242\n")
                os.utime(path, (time.time() - age,) * 2)
                return path
            claim("killed", age=genomes.CLAIM_STALE + 60)  # its heartbeat stopped: a job killed at its time limit
            began = time.time()
            eukore_build.others_end(opts, ["killed", "no_claim"], poll=0.05)
            self.assertLess(time.time() - began, 2)
            live = claim("running")  # another node's, touched now: waited for until that run ends (its claim removed)
            timer = threading.Timer(0.4, os.remove, [live])
            timer.start()
            began = time.time()
            eukore_build.others_end(opts, ["running", "killed"], poll=0.05)
            timer.join()
            self.assertGreaterEqual(time.time() - began, 0.3)
            self.assertFalse(os.path.exists(live))
            self.assertTrue(any("still run by another build" in str(c) for c in said.call_args_list))

    def test_the_backfill_runs_as_its_memory_fits(self):
        # a clade hold-out is a call of the whole genome: begun by the scheduler with a genome run's memory, the
        # other steps with a worker's alone (until 2026-11-10 -j at once whatever the genomes' sizes)
        import types  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        jobs = {a: {"accession": a, "size": 3e9} for a in ("hold_out", "tiara", "whole", "not_done")}
        needs = {"hold_out": ["held_out"], "tiara": ["tiara", "rrna"], "whole": []}

        def scheduled(_opts, todo, sizes, _progress, _note, task=None):
            self.assertIs(task, eukore_build.backfill_one)
            self.assertEqual(sorted(j["accession"] for j, _s in todo), ["hold_out", "tiara"])
            self.assertEqual(sizes, {"hold_out": 3e9, "tiara": 0})
            yield "hold_out", ["held_out_genus"], None
            yield "tiara", ["tiara"], "rRNA extraction failed: x"
        opts = types.SimpleNamespace(work="WORK", genomes="GENOMES")
        with mock.patch.object(eukore_build, "say") as said, \
                mock.patch.object(eukore_build, "progress_of", return_value=Steps()), \
                mock.patch.object(eukore_build.genomes, "is_done", side_effect=lambda _w, a: a != "not_done"), \
                mock.patch.object(eukore_build.genomes, "backfill_needs",
                                  side_effect=lambda j, _s: needs[j["accession"]]), \
                mock.patch.object(eukore_build, "schedule", side_effect=scheduled) as schedule:
            eukore_build.backfill_pass(opts, jobs, {})
            self.assertEqual(schedule.call_count, 1)
            self.assertIn("{'held_out_genus': 1, 'tiara': 1} made, 1 failed", said.call_args_list[-1].args[0])
            self.assertLess(eukore_build.run_memory(0), eukore_build.run_memory(3e9) / 5)
            schedule.reset_mock()
            needs["hold_out"] = needs["tiara"] = []
            eukore_build.backfill_pass(opts, jobs, {})  # nothing left out: no pool
            schedule.assert_not_called()

    def test_a_backfilled_hold_out_comes_from_the_local_disk(self):
        # --local_tmp is another file system than WORK: the hold-out's table is moved from it, not renamed (until
        # 2026-11-10 os.replace failed there after the whole call, and on every rerun)
        import errno  # noqa: PLC0415
        import json  # noqa: PLC0415
        from eukore import genomes  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        replace = os.replace

        def one_file_system(src, dst):
            if os.path.commonpath([src, local]) == local and os.path.commonpath([dst, local]) != local:
                raise OSError(errno.EXDEV, "Invalid cross-device link")
            replace(src, dst)

        def held_out(_genome, _job, _settings, _odb, folder, _scratch, rank, _earlier, _index_dir=None):
            with open(os.path.join(folder, genomes.held_out_file(rank)), "w") as fh:
                fh.write("its calls\n")
            return {"left_out": 3}, None
        with tempfile.TemporaryDirectory() as tmp:
            work, local = os.path.join(tmp, "work"), os.path.join(tmp, "local")
            os.makedirs(local)
            folder = genomes.done_dir(work, "GCA_1.1")
            os.makedirs(folder)
            with open(os.path.join(folder, "done.json"), "w") as fh:
                json.dump({"held_out": {"calls": 5}}, fh)
            job = {"accession": "GCA_1.1", "path": os.path.join(tmp, "genome.fna.gz")}
            settings = {"work": work, "genomes": tmp, "scratch": local, "held_out_ranks": ("genus",)}
            with mock.patch.object(genomes, "backfill_needs", return_value=["held_out"]), \
                    mock.patch.object(genomes, "clade_held_out", side_effect=held_out), \
                    mock.patch.object(genomes, "odb_of", return_value=None), \
                    mock.patch.object(genomes.os, "replace", side_effect=one_file_system):
                self.assertEqual(genomes.backfill((job, settings)), ("GCA_1.1", ["held_out_genus"], None))
            with open(os.path.join(folder, genomes.held_out_file("genus"))) as fh:
                self.assertEqual(fh.read(), "its calls\n")
            self.assertEqual(genomes.done_summary(work, "GCA_1.1")["held_out"]["ranks"], {"genus": {"left_out": 3}})
            self.assertEqual(os.listdir(local), [])  # its folder there removed

    @unittest.skipUnless(sys.platform.startswith("linux"), "bash, and stand-ins on the PATH")
    def test_the_slurm_chain(self):
        # tools/eukore.sh with stand-ins of sbatch, micromamba and python3: `build` submits the array and the final
        # job after it (afterany), each job's step makes one call of the build
        import stat  # noqa: PLC0415
        import subprocess  # noqa: PLC0415
        script = os.path.join(os.path.dirname(HERE), "tools", "eukore.sh")
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir, log = os.path.join(tmp, "bin"), os.path.join(tmp, "calls")
            os.makedirs(bin_dir)
            # sbatch --parsable says "<job id>;<cluster>"; it refuses the job named REFUSE as a QOS does
            sbatch = ('if [[ -n ${REFUSE:-} && " $* " == *" -J $REFUSE "* ]]; then echo "QOS policy" >&2; exit 1; fi\n'
                      'echo "sbatch $*" >> "$CALLS"; n=$(grep -c ^sbatch "$CALLS"); echo "$((100 + n));cluster"')
            python3 = 'echo "python3 $*" >> "$CALLS"'
            scancel = 'echo "scancel $*" >> "$CALLS"'
            for name, body in (("sbatch", sbatch), ("micromamba", "true"), ("python3", python3), ("scancel", scancel)):
                with open(os.path.join(bin_dir, name), "w", newline="\n") as fh:
                    fh.write("#!/bin/bash\n" + body + "\n")
                os.chmod(os.path.join(bin_dir, name), stat.S_IRWXU)
            conf = os.path.join(tmp, "eukore.conf")
            with open(conf, "w", newline="\n") as fh:
                fh.write(f"PDIR=/p\nWORK={tmp}/w\nOUT=/o\nGEN=/g\nVERSION=2026.11\nLOGS={tmp}/runs\nTASKS=3\n"
                         "SOURCES=${TEST_SOURCES:-}\n"
                         "SCR=${SLURM_LOCAL_SCRATCH:+$SLURM_LOCAL_SCRATCH/Eukore}\nBUILD_OPTS=(--groups fungi)\n"
                         "FINAL_OPTS=(--keep_gb 1000)\nFETCH_OPTS=(--downloads 12)\nSBATCH_ALL=(--time=1:00:00)\n"
                         "SBATCH_PROCESS=(--mem=200G)\nSBATCH_FINAL=(--mem=800G)\n")
            env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "EUKORE_")) and k != "AFTER"}
            env.update(PATH=bin_dir + os.pathsep + env["PATH"], CALLS=log, EUKORE_CONF=conf)

            def run(*args, **more):
                got = subprocess.run(["bash", script, *args], env=dict(env, **more), capture_output=True, text=True)
                calls = []
                if os.path.isfile(log):
                    with open(log) as fh:
                        calls = fh.read().splitlines()
                    os.remove(log)
                return got, calls
            # no chain on a WORK no fetch filled (2026-11-11: every task crashed on NCBI's taxonomy missing), unless
            # it waits for a fetch job
            got, calls = run("build")
            self.assertEqual((got.returncode, calls), (1, []))
            self.assertIn(f"no sources in {tmp}/w", got.stderr)
            got, calls = run("build", AFTER="77")
            self.assertEqual((got.returncode, len(calls)), (0, 2), got.stderr)
            os.makedirs(os.path.join(tmp, "w", "sources"))
            with open(os.path.join(tmp, "w", "sources", "taxdump.tar.gz"), "w") as fh:
                fh.write("its taxonomy")
            got, calls = run("build")
            self.assertEqual(got.returncode, 0, got.stderr)
            self.assertEqual(len(calls), 2)
            for said in ("-J Eukore-process", "--time=1:00:00", "--mem=200G", "--array=1-3", f"{script} _process",
                         f"-o {tmp}/runs/eukore.%j.otxt"):
                self.assertIn(said, calls[0])
            self.assertNotIn("--dependency", calls[0])
            for said in ("-J Eukore-final", "--mem=800G", "--dependency=afterany:101 ", f"{script} _final"):
                self.assertIn(said, calls[1])
            self.assertIn("array job 101, 3 tasks", got.stdout)
            self.assertIn("job 102", got.stdout)
            # after a fetch job: both wait for it to end well, and are cancelled where it did not
            got, calls = run("build", "5", AFTER="77")
            self.assertIn("--array=1-5 --dependency=afterok:77 --kill-on-invalid-dep=yes", calls[0])
            self.assertIn("--dependency=afterany:101,afterok:77 --kill-on-invalid-dep=yes", calls[1])
            # a job sbatch refuses ends the call (2026-11-10: the QOS refused the array, and `build` went on to a
            # final job with a dependency on no job); a final job refused cancels the array before it
            got, calls = run("build", REFUSE="Eukore-process")
            self.assertEqual((got.returncode, calls, got.stdout), (1, [], ""))
            self.assertIn("sbatch refused Eukore-process, submitted as:\n sbatch --parsable -J Eukore-process "
                          "--time=1:00:00 --mem=200G --array=1-3 -o", got.stderr)
            got, calls = run("build", REFUSE="Eukore-final")
            self.assertEqual((got.returncode, calls[1:], got.stdout), (1, ["scancel 101"], ""))
            self.assertIn("array job 101 cancelled", got.stderr)
            got, calls = run("fetch", REFUSE="Eukore-fetch")
            self.assertEqual((got.returncode, calls, got.stdout), (1, [], ""))
            # inside the jobs: one call of the build each
            job = dict(SLURM_JOB_ID="101", SLURM_CPUS_PER_TASK="64", SLURM_LOCAL_SCRATCH="/ssd", EUKORE_SH=script)
            common = f"python3 /p/eukore_build.py --build_dir {tmp}/w --out /o --version 2026.11 --stage "
            got, calls = run("_process", **job)
            self.assertEqual(calls, [common + "process --genomes /g --groups fungi -t 64 --genomes_only "
                                     "--local_tmp /ssd/Eukore"], got.stderr)
            got, calls = run("_final", REVERSION="1", **job)
            self.assertEqual(calls, [common + "process,package --genomes /g --reversion --groups fungi -t 64 "
                                     "--local_tmp /ssd/Eukore --keep_gb 1000"], got.stderr)
            got, calls = run("_fetch", **dict(job, SLURM_LOCAL_SCRATCH=""))  # no local disk there: beside the genomes
            self.assertEqual(calls, [common + "fetch --genomes /g --groups fungi -t 64 --downloads 12"], got.stderr)
            # SOURCES (since 2026-11-11): each job's --sources, and the folder `build` looks for the fetch's files in
            got, calls = run("_process", TEST_SOURCES="/s", **job)
            self.assertIn("--genomes /g --sources /s --groups fungi", calls[0])
            got, calls = run("build", TEST_SOURCES=f"{tmp}/kept")
            self.assertEqual((got.returncode, calls), (1, []))
            self.assertIn(f"no sources in {tmp}/kept", got.stderr)
            got, calls = run("_process", **dict(job, SLURM_LOCAL_SCRATCH=""))  # the process stage needs one
            self.assertEqual((got.returncode, calls), (1, []))
            self.assertIn("no local temporary folder", got.stderr)
            got, calls = run("build", EUKORE_CONF=os.path.join(tmp, "none.conf"))
            self.assertEqual((got.returncode, calls), (1, []))
            self.assertIn("no settings at", got.stderr)


if __name__ == "__main__":
    unittest.main()
