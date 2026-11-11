#!/usr/bin/env python3
"""Tests of a done genome's results as one file (store.py; 2026-11-10) and of that day's fixes of the build's core: a
genome read in both forms, changed and packed beside other builds; a packed build's package equal to a folder
build's; a run again that takes the old results' place; a failed rRNA extraction that keeps the genome's calls; a
dead claim taken over by one build; a plan that reads each done.json once; miniprot's index made once a genome.

    python3 -m unittest discover -s tests -t tests
"""

import contextlib
import gzip
import io
import json
import os
import sys
import tempfile
import threading
import time
import types
import unittest
import zipfile
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, ncbi as ncbi_mod, odb as odb_mod, package, sources, store  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST = "GCA_000146045.2"
LOCI = "locus\tstatus\n1at2759\tSingle\n2at2759\tMissing\n"
WINDOWS = ">1at2759|c1|1-8\nACGTACGT\n"


def folder_of(work, acc="G1", summary=None):
    """A done genome's folder: done.json, loci.tsv, windows.fna.gz."""
    folder = genomes.done_dir(work, acc)
    os.makedirs(folder)
    with open(os.path.join(folder, "loci.tsv"), "w") as fh:
        fh.write(LOCI)
    with gzip.open(os.path.join(folder, "windows.fna.gz"), "wt") as fh:
        fh.write(WINDOWS)
    with open(os.path.join(folder, "done.json"), "w") as fh:
        json.dump(summary or {"accession": acc}, fh)
    return folder


def read_all(folder):
    return (store.exists(folder), store.names(folder), store.has(folder, "loci.tsv"), store.has(folder, "none.tsv"),
            store.read(folder, "loci.tsv"), store.read(folder, "none.tsv"), store.text(folder, "windows.fna.gz"),
            store.open_text(folder, "loci.tsv").read(), store.open_text(folder, "none.tsv"),
            store.rows(folder, "loci.tsv"), store.rows(folder, "none.tsv"), store.table(folder, "loci.tsv"),
            store.stamp(folder), store.stamp(folder, "none.tsv"))


class StoreTest(unittest.TestCase):
    def tearDown(self):
        store.close()

    def test_a_genome_read_in_both_forms(self):
        with tempfile.TemporaryDirectory() as work:
            folder = folder_of(work)
            none = genomes.done_dir(work, "G0")  # no such genome
            self.assertEqual((store.exists(none), store.names(none), store.has(none, "loci.tsv"),
                              store.text(none, "loci.tsv"), store.rows(none, "loci.tsv"), store.stamp(none)),
                             (False, set(), False, None, [], None))
            self.assertEqual(store.unpacked(os.path.join(work, "done")), [folder])
            before = read_all(folder)
            self.assertEqual(before[:8], (True, {"done.json", "loci.tsv", "windows.fna.gz"}, True, False,
                                          LOCI.encode(), None, WINDOWS, LOCI))
            self.assertEqual(before[8:12], (None, [["1at2759", "Single"], ["2at2759", "Missing"]], [],
                                            [dict(locus="1at2759", status="Single"),
                                             dict(locus="2at2759", status="Missing")]))
            self.assertIsNotNone(before[12])
            self.assertIsNone(before[13])
            # several files in one look, as the threads of the package's passes read them (since 2026-11-11)
            wanted = ("loci.tsv", "windows.fna.gz", "none.tsv")
            looked = store.members(folder, wanted)
            self.assertEqual((looked["loci.tsv"], looked["none.tsv"]), (LOCI.encode(), None))
            self.assertEqual(store.decoded("windows.fna.gz", looked["windows.fna.gz"]), WINDOWS)
            self.assertEqual(store.table_of(store.decoded("loci.tsv", looked["loci.tsv"])), before[11])
            self.assertEqual((store.members(none, wanted), store.decoded("x.gz", None), store.table_of(None)),
                             (None, None, []))
            self.assertFalse(store.is_packed(folder))
            with open(os.path.join(folder, ".host-1.fcs.tsv.partial"), "w") as fh:  # a killed change's
                fh.write("half")
            self.assertEqual(store.pack(folder), folder + ".zip")
            # one file in place of the folder, read as it was; its stamp the same (a pack is no new result)
            self.assertEqual(os.listdir(os.path.join(work, "done")), ["G1.zip"])
            self.assertTrue(store.is_packed(folder))
            self.assertEqual(read_all(folder), before)
            self.assertEqual(store.members(folder, wanted), looked)
            self.assertEqual(store.unpacked(os.path.join(work, "done")), [])
            dest = os.path.join(work, "w.fna.gz")
            self.assertTrue(store.extract(folder, "windows.fna.gz", dest))
            with gzip.open(dest, "rt") as fh:
                self.assertEqual(fh.read(), WINDOWS)
            self.assertFalse(store.extract(folder, "none.tsv", dest + ".x"))
            self.assertFalse(os.path.exists(dest + ".x"))
            with zipfile.ZipFile(folder + ".zip") as z:  # text deflated, a compressed file as it is
                self.assertEqual({i.filename: i.compress_type for i in z.infolist()},
                                 {"done.json": zipfile.ZIP_DEFLATED, "loci.tsv": zipfile.ZIP_DEFLATED,
                                  "windows.fna.gz": zipfile.ZIP_STORED})

    def test_a_genome_changed_in_both_forms(self):
        with tempfile.TemporaryDirectory() as work:
            for packed in (False, True):
                folder = folder_of(work, f"G{packed:d}")
                if packed:
                    store.pack(folder)
                stamp = store.stamp(folder)
                source = os.path.join(work, "rrna.tsv")
                with open(source, "w") as fh:
                    fh.write("kind\n18S\n")
                store.put(folder, {"done.json": b'{"tiara": "done"}', "rrna.tsv": source, "loci.tsv": None,
                                   "never.tsv": None})
                self.assertEqual(store.names(folder), {"done.json", "rrna.tsv", "windows.fna.gz"})
                self.assertEqual((store.text(folder, "done.json"), store.rows(folder, "rrna.tsv"),
                                  store.text(folder, "windows.fna.gz")), ('{"tiara": "done"}', [["18S"]], WINDOWS))
                self.assertNotEqual(store.stamp(folder), stamp)  # a new done.json: another build's plan sees it
                self.assertTrue(os.path.isfile(source))  # read, left where it is
                self.assertEqual(store.is_packed(folder), packed)
                self.assertEqual(os.path.isdir(folder), not packed)
                # nothing left beside it: no lock, no half-written file
                self.assertEqual([n for n in os.listdir(os.path.join(work, "done")) if n.startswith(f"G{packed:d}")],
                                 [f"G{packed:d}" + (".zip" if packed else "")])
            new = genomes.done_dir(work, "G9")  # a genome with nothing yet: a folder
            store.put(new, {"fcs.tsv": b"window\n"})
            self.assertEqual((os.path.isdir(new), store.names(new)), (True, {"fcs.tsv"}))

    def test_a_killed_packs_folder_left_behind(self):
        with tempfile.TemporaryDirectory() as work:
            folder = folder_of(work)
            with mock.patch.object(store.shutil, "rmtree"):  # killed once the zip was in place
                store.pack(folder)
            self.assertTrue(os.path.isdir(folder) and store.is_packed(folder))
            with open(os.path.join(folder, "loci.tsv"), "w") as fh:
                fh.write("locus\tstatus\nstale\tSingle\n")
            self.assertEqual(store.rows(folder, "loci.tsv")[0][0], "1at2759")  # the zip is read
            store.put(folder, {"fcs.tsv": b"window\n"})  # a change is the zip's, the folder removed
            self.assertFalse(os.path.isdir(folder))
            self.assertEqual(store.names(folder), {"done.json", "loci.tsv", "windows.fna.gz", "fcs.tsv"})
            # the same met by a pack: the folder removed, the zip (changed since) kept
            os.makedirs(folder)
            with open(os.path.join(folder, "done.json"), "w") as fh:
                fh.write("{}")
            self.assertEqual(store.unpacked(os.path.join(work, "done")), [folder])
            store.pack(folder)
            self.assertFalse(os.path.isdir(folder))
            self.assertEqual((store.has(folder, "fcs.tsv"), json.loads(store.text(folder, "done.json"))),
                             (True, {"accession": "G1"}))
            # a folder another build's pack took meanwhile: nothing written over its zip
            before = store.stamp(folder)
            self.assertEqual(store.pack(os.path.join(work, "done", "gone"), folder), folder + ".zip")
            self.assertEqual(store.stamp(folder), before)
            self.assertIsNone(store.pack(os.path.join(work, "done", "gone")))

    def test_changes_of_a_genome_wait_for_each_other(self):
        with tempfile.TemporaryDirectory() as work:
            folder = folder_of(work)
            store.pack(folder)
            lock = folder + ".zip.lock"
            with open(lock, "w"):  # another build's change going on
                pass
            threading.Timer(0.3, os.remove, [lock]).start()
            began = time.time()
            store.put(folder, {"a.tsv": b"a\n"})
            self.assertGreaterEqual(time.time() - began, 0.25)
            self.assertFalse(os.path.exists(lock))
            with open(lock, "w"):  # a killed build's: taken over
                pass
            old = time.time() - store.LOCK_STALE - 10
            os.utime(lock, (old, old))
            began = time.time()
            store.put(folder, {"b.tsv": b"b\n"})
            self.assertLess(time.time() - began, 5)
            # this process' threads: each one's file in the zip, none lost to another's rewrite
            threads = [threading.Thread(target=store.put, args=(folder, {f"t{i}.tsv": b"x\n"})) for i in range(8)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            self.assertEqual(store.names(folder), {"done.json", "loci.tsv", "windows.fna.gz", "a.tsv", "b.tsv"}
                             | {f"t{i}.tsv" for i in range(8)})
            # folders packed by several threads (pack_pass), one of them met twice
            folders = [folder_of(work, f"P{i}") for i in range(6)]
            threads = [threading.Thread(target=store.pack, args=(f,)) for f in folders + folders[:2]]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            self.assertEqual([(os.path.isdir(f), store.rows(f, "loci.tsv")[0]) for f in folders],
                             [(False, ["1at2759", "Single"])] * 6)

    @unittest.skipUnless(hasattr(os, "fork"), "no fork")
    def test_a_pools_worker_opens_its_own_zip(self):
        with tempfile.TemporaryDirectory() as work:
            folder = folder_of(work)
            store.pack(folder)
            self.assertEqual(store.text(folder, "loci.tsv"), LOCI)  # open in this process
            pid = os.fork()
            if pid == 0:  # the worker: nothing of its parent's open files, its own read right
                ok = not store._OPEN and store.text(folder, "loci.tsv") == LOCI
                os._exit(0 if ok else 1)
            self.assertEqual(os.waitpid(pid, 0)[1], 0)
            self.assertEqual(store.text(folder, "windows.fna.gz"), WINDOWS)


def world(tmp, **more):
    """The fixtures' yeast genome and what a run of it needs -> (ncbi, jobs, assemblies, settings)."""
    ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
    root = fx.odb_world(tmp)
    odb = odb_mod.Odb(root, fx.LOCI)
    seq, gff, _spans = fx.genome()
    genomes_dir = os.path.join(tmp, "out", "genomes")
    os.makedirs(genomes_dir)
    with gzip.open(genomes.genome_file(genomes_dir, YEAST), "wt") as fh:
        fh.write(f">c1\n{seq}\n")
    settings = dict(work=os.path.join(tmp, "work"), genomes=genomes_dir, miniprot=fx.miniprot(tmp, gff),
                    hmmsearch=fx.hmmsearch(tmp), threads=1, staged_only=True, odb=root, loci=fx.LOCI,
                    hmm_db=odb.hmm_db(os.path.join(tmp, "loci.hmm")),
                    odb_lineages=odb_mod.species_lineages(odb, ncbi), tiara=None, rrna=None,
                    vsearch=fx.vsearch(tmp))
    settings.update(more)
    gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
    return ncbi, odb, genomes.jobs(ncbi, gb), gb, settings


def rrna_hook(path, job, out):
    with open(os.path.join(out, "rrna.tsv"), "w") as fh:
        fh.write("kind\n")
    return {"18S": 0}


def rrna_fails(path, job, out):
    with open(os.path.join(out, "rrna.tsv"), "w") as fh:
        fh.write("half")
    raise RuntimeError("barrnap failed (exit 1)")


class PackedBuildTest(unittest.TestCase):
    def tearDown(self):
        store.close()

    def test_a_packed_genome_into_the_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi, _odb, jobs, gb, folders = world(tmp)
            packed = dict(folders, work=os.path.join(tmp, "packed"), pack=True)
            jobs = {YEAST: jobs[YEAST]}
            collected = {}
            for settings in (folders, packed):
                work = settings["work"]
                acc, summary, why = genomes.run((jobs[YEAST], settings))
                self.assertEqual((acc, why, summary["loci"]), (YEAST, None, 2))
                # one entry a genome: its folder or its file, no claim, no partial left
                self.assertEqual(os.listdir(os.path.join(work, "done")),
                                 [YEAST + (".zip" if settings is packed else "")])
                self.assertTrue(genomes.is_done(work, YEAST, fx.LOCI))
                self.assertEqual(genomes.searched(work, YEAST), set(fx.LOCI))
                collected[work] = package.collect(ncbi, jobs, gb, work, fx.vsearch(tmp), 1)
            a, b = (genomes.done_dir(s["work"], YEAST) for s in (folders, packed))
            self.assertEqual(store.names(a), store.names(b))
            for name in ("loci.tsv", "copies.fna", "features.tsv", "windows.fna.gz"):
                self.assertEqual(store.text(a, name), store.text(b, name), name)
            for key in ("copy_rows", "cds", "intron_rows", "copy_count", "done"):
                self.assertEqual(collected[folders["work"]][key], collected[packed["work"]][key], key)
            # a run again of the packed genome: its plan's stamp met, the new results in the old one's place
            stamp = genomes.done_stamp(packed["work"], YEAST)
            self.assertEqual(genomes.run((dict(jobs[YEAST], done_stamp="another"), packed))[2], genomes.DONE_ELSEWHERE)
            store.put(b, {"fcs.tsv": b"window\n"})  # of the old results: gone with them
            self.assertIsNone(genomes.run((dict(jobs[YEAST], done_stamp=stamp), packed))[2])
            self.assertEqual(os.listdir(os.path.join(packed["work"], "done")), [YEAST + ".zip"])
            self.assertFalse(store.has(b, "fcs.tsv"))
            # a WORK of folders run by a build that packs: the genome packed by its run again; without (the library's
            # default, --done_folders), a packed genome stays packed
            self.assertIsNone(genomes.run((jobs[YEAST], dict(folders, pack=True)))[2])
            self.assertEqual(os.listdir(os.path.join(folders["work"], "done")), [YEAST + ".zip"])
            self.assertIsNone(genomes.run((jobs[YEAST], folders))[2])
            self.assertEqual(os.listdir(os.path.join(folders["work"], "done")), [YEAST + ".zip"])

    def test_a_run_again_takes_the_old_results_place(self):
        # until 2026-11-10 the old folder was deleted, file by file, before the new one was renamed in: a kill
        # meanwhile left half a genome that read as done (done.json is its last file by name only). Now each step is
        # one rename: the genome is whole, or (a folder's two renames apart) not there and run again; a packed
        # genome's results change in one rename, whole throughout
        for pack, expected in ((False, [True, None, None, True]), (True, [True, True])):
            with tempfile.TemporaryDirectory() as tmp:
                _ncbi, _odb, jobs, _gb, settings = world(tmp, pack=pack)
                work = settings["work"]
                self.assertIsNone(genomes.run((jobs[YEAST], settings))[2])
                folder = genomes.done_dir(work, YEAST)
                store.put(folder, {"fcs.tsv": b"window\n"})  # of the old results: gone with them
                seen, real = [], os.replace

                def state():
                    names = store.names(folder)
                    return None if not names else {"done.json", "loci.tsv", "copies.fna"} <= names

                def replace(source, dest):
                    ours = os.path.dirname(dest) == os.path.dirname(folder) and YEAST in os.path.basename(dest)
                    if ours:
                        seen.append(state())
                    real(source, dest)
                    if ours:
                        seen.append(state())
                with mock.patch.object(os, "replace", replace):
                    self.assertIsNone(genomes.run((jobs[YEAST], settings))[2])
                self.assertEqual(seen, expected)
                # the old results gone, nothing left beside the new
                self.assertEqual(os.listdir(os.path.join(work, "done")), [YEAST + (".zip" if pack else "")])
                self.assertFalse(store.has(folder, "fcs.tsv"))

    def test_a_failed_rrna_extraction_keeps_the_genomes_calls(self):
        # until 2026-11-10 the genome failed with it: its locus calls thrown away, and made again by every run
        for pack in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                _ncbi, _odb, jobs, _gb, settings = world(tmp, rrna=rrna_fails, pack=pack)
                acc, summary, why = genomes.run((jobs[YEAST], settings))
                self.assertEqual((why, summary["loci"]), (None, 2))
                self.assertIn("barrnap failed", summary["rrna_failed"])
                work, folder = settings["work"], genomes.done_dir(settings["work"], YEAST)
                self.assertTrue(genomes.is_done(work, YEAST, fx.LOCI))
                self.assertFalse(store.has(folder, "rrna.tsv"))  # no half-written table taken for its rRNA
                self.assertEqual(genomes.backfill_needs(jobs[YEAST], settings), ["rrna"])
                # the backfill tries again: failing still, the genome stays as it is; then made
                self.assertEqual(genomes.backfill((jobs[YEAST], settings))[1:],
                                 ([], "rRNA extraction failed: barrnap failed (exit 1)"))
                self.assertIn("rrna_failed", genomes.done_summary(work, YEAST))
                mended = dict(settings, rrna=rrna_hook)
                self.assertEqual(genomes.backfill((jobs[YEAST], mended)), (YEAST, ["rrna"], None))
                summary = genomes.done_summary(work, YEAST)
                self.assertEqual((summary.get("rrna_failed"), summary["rrna"], store.text(folder, "rrna.tsv")),
                                 (None, {"18S": 0}, "kind\n"))
                self.assertEqual(genomes.backfill_needs(jobs[YEAST], mended), [])
                self.assertEqual(os.listdir(os.path.join(work, "done")), [YEAST + (".zip" if pack else "")])
                self.assertEqual(store.is_packed(folder), pack)

    def test_a_backfill_looks_again_once_the_genome_is_its_own(self):
        # two builds' passes (a SLURM array's tasks) planned the same step: the second, the genome claimed once the
        # first is through, finds it made (until 2026-11-10 it made it again)
        with tempfile.TemporaryDirectory() as tmp:
            _ncbi, _odb, jobs, _gb, settings = world(tmp, rrna=rrna_fails)
            self.assertIsNone(genomes.run((jobs[YEAST], settings))[2])
            calls = []

            def counted(path, job, out):
                calls.append(job["accession"])
                return rrna_hook(path, job, out)
            mended = dict(settings, rrna=counted)
            real = genomes.claim

            def claim(work, accession):  # the other build's backfill, through before this one's claim
                if not calls:
                    with mock.patch.object(genomes, "claim", real):
                        self.assertEqual(genomes.backfill((jobs[YEAST], mended))[1], ["rrna"])
                return real(work, accession)
            with mock.patch.object(genomes, "claim", claim):
                self.assertEqual(genomes.backfill((jobs[YEAST], mended)), (YEAST, [], None))
            self.assertEqual(calls, [YEAST])
            # a genome another build holds is left to it
            with genomes.claimed(settings["work"], YEAST):
                with mock.patch.object(genomes, "claim", return_value=None):
                    self.assertEqual(genomes.backfill((jobs[YEAST], dict(mended, tiara="tiara"))), (YEAST, [], None))

    def test_the_index_made_once_for_a_genomes_calls(self):
        # a dataset species' genome is called for itself and again for each hold-out: with --local_tmp, miniprot's
        # index of it is made by the first call and read by the others (until 2026-11-10 made by each)
        with tempfile.TemporaryDirectory() as tmp:
            ncbi, odb, jobs, _gb, settings = world(tmp)
            log = os.path.join(tmp, "miniprot.log")
            logged = fx.tool(tmp, "miniprot_logged",
                             f'#!/bin/sh\necho "$@" >> {log}\nexec {settings["miniprot"]} "$@"\n')
            own = dict(settings, miniprot=logged, odb_own=odb_mod.own_taxids(odb, ncbi), held_out_ranks=[])
            for scratch, builds in ((None, 0), (os.path.join(tmp, "scratch"), 1)):
                if os.path.exists(log):
                    os.remove(log)
                work = os.path.join(tmp, f"work{builds}")
                self.assertIsNone(genomes.run((jobs[YEAST], dict(own, work=work, scratch=scratch)))[2])
                self.assertTrue(store.has(genomes.done_dir(work, YEAST), genomes.HELD_OUT))
                with open(log) as fh:
                    calls = [line.split() for line in fh]
                self.assertEqual(sum("-d" in c for c in calls), builds)
                # the genome's own call and its hold-out's: both on the index where it is kept
                self.assertEqual(sum(c[-2].endswith(".mpi") for c in calls if "-d" not in c), 2 * builds)
                self.assertEqual(os.listdir(os.path.join(work, "done")), [YEAST])  # no index among its results
            self.assertEqual(os.listdir(os.path.join(tmp, "scratch")), [])
            # not for a genome called once, nor for one whose index would fill the node's disk
            job = jobs[YEAST]
            path, scratch = genomes.genome_file(settings["genomes"], YEAST), os.path.join(tmp, "scratch")
            with_scratch = dict(own, scratch=scratch)
            self.assertEqual(genomes.index_dir_of(job, with_scratch, path, "/s/x"), os.path.join("/s/x", "mpi"))
            self.assertIsNone(genomes.index_dir_of(job, dict(with_scratch, odb_own={}), path, "/s/x"))
            self.assertIsNone(genomes.index_dir_of(dict(job, size=genomes.INDEX_MAX_BASES + 1), with_scratch, path,
                                                   "/s/x"))
            self.assertIsNone(genomes.index_dir_of(dict(job, kind="transcriptome"), with_scratch, path, "/s/x"))


class PlanTest(unittest.TestCase):
    def tearDown(self):
        store.close()

    def test_a_plan_reads_each_done_json_once(self):
        with tempfile.TemporaryDirectory() as work:
            for acc in ("G1", "G2"):
                folder_of(work, acc, {"accession": acc, "searched": ["1at2759"], "genome_record": 1})
            store.pack(genomes.done_dir(work, "G2"))
            real, read = store.text, []

            def text(folder, name):
                read.append((os.path.basename(folder), name))
                return real(folder, name)
            with mock.patch.object(store, "text", text):
                with genomes.cached_summaries():
                    for acc in ("G1", "G2", "G3"):
                        for _ in range(3):
                            self.assertEqual(genomes.is_done(work, acc), acc != "G3")
                            self.assertEqual(genomes.is_done(work, acc, ["1at2759"]), acc != "G3")
                            self.assertEqual(genomes.record_stale(work, acc), acc != "G3")
                    self.assertEqual(sorted(read), [("G1", "done.json"), ("G2", "done.json"), ("G3", "done.json")])
                    # this process' own change is seen at once
                    summary = genomes.done_summary(work, "G1")
                    summary["genome_record"] = genomes.GENOME_RECORD
                    genomes.write_summary(work, "G1", summary)
                    self.assertFalse(genomes.record_stale(work, "G1"))
                    # each caller its own copy: one's change is not another's
                    genomes.done_summary(work, "G2")["searched"].append("x")
                    self.assertEqual(genomes.done_summary(work, "G2")["searched"], ["1at2759"])
                del read[:]
                # outside a plan: every read the file's (another build may have written since)
                self.assertTrue(genomes.is_done(work, "G1") and genomes.is_done(work, "G1"))
                self.assertEqual(len(read), 2)

    def test_a_plan_reads_the_done_genomes_at_once_and_alone(self):
        # prefetch (since 2026-11-11): each done genome read once by threads, outside store's lock and its open zips;
        # a genome WORK/done does not list read no further; the plan's questions after it answered from what it read
        from eukore import metadata  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as work:
            for acc in ("G1", "G2"):
                folder_of(work, acc, {"accession": acc, "searched": ["1at2759"], "genome_record": 1})
            store.pack(genomes.done_dir(work, "G2"))
            os.makedirs(genomes.done_dir(work, "G4"))  # a folder without its done.json: not done
            stamps = {a: store.stamp(genomes.done_dir(work, a)) for a in ("G1", "G2")}
            self.assertEqual(store.snapshot(genomes.done_dir(work, "G3")), None)
            names, data, stamp = store.snapshot(genomes.done_dir(work, "G2"))
            self.assertEqual((json.loads(data)["accession"], stamp), ("G2", stamps["G2"]))
            self.assertIn("windows.fna.gz", names)
            looked, real = [], store.snapshot

            def snapshot(folder, name="done.json"):
                looked.append(os.path.basename(folder))
                return real(folder, name)
            with mock.patch.object(store, "snapshot", snapshot), \
                    mock.patch.object(store, "text", side_effect=AssertionError("read one at a time")), \
                    mock.patch.object(store, "has", side_effect=AssertionError("read one at a time")), \
                    mock.patch.object(store, "stamp", side_effect=AssertionError("read one at a time")):
                with genomes.cached_summaries():
                    genomes.prefetch(work, ["G1", "G2", "G3", "G4"], threads=3)
                    self.assertEqual(sorted(looked), ["G1", "G2", "G4"])  # G3: not listed, not read
                    self.assertEqual([genomes.is_done(work, a, ["1at2759"]) for a in ("G1", "G2", "G3", "G4")],
                                     [True, True, False, False])
                    self.assertEqual([genomes.done_stamp(work, a) for a in ("G1", "G2", "G3")],
                                     [stamps["G1"], stamps["G2"], None])
                    self.assertTrue(genomes.has(work, "G2", "loci.tsv"))
                    self.assertFalse(genomes.has(work, "G2", genomes.HELD_OUT))
                    self.assertEqual(metadata.missing(genomes.done_dir(work, "G1"), "sourmash",
                                                      genomes.names_of(work, "G1")),
                                     {metadata.STATS, metadata.SKETCH})
            # outside a plan: store's own reads again
            self.assertIsNone(genomes.names_of(work, "G1"))
            self.assertTrue(genomes.has(work, "G1", "loci.tsv"))
            self.assertEqual(genomes.done_stamp(work, "G1"), stamps["G1"])

    def test_a_cache_inside_another_is_the_outers(self):
        # since 2026-11-11 (the package's passes cache too): an inner block's end no longer ends the outer's cache
        from eukore import util  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as work:
            folder_of(work, "G1", {"accession": "G1"})
            real, read = store.text, []

            def text(folder, name):
                read.append(name)
                return real(folder, name)
            with mock.patch.object(store, "text", text), genomes.cached_summaries():
                self.assertTrue(genomes.is_done(work, "G1"))
                with genomes.cached_summaries():
                    self.assertTrue(genomes.is_done(work, "G1"))
                self.assertTrue(genomes.is_done(work, "G1"))
            self.assertEqual(read, ["done.json"])
            # threads read in the items' order, whatever order they finish in
            self.assertEqual(list(util.threaded(lambda i: (time.sleep(0.001 * (i % 3)), i)[1], range(50), 4)),
                             list(range(50)))

    def test_the_rrna_records_from_collects_reads(self):
        # since 2026-11-11 the package's rRNA records are the rows collect's pool read, not each genome read again
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import eukore_build  # noqa: PLC0415
        from eukore import rrna  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            work = os.path.join(tmp, "work")
            folder = folder_of(work, YEAST, {"accession": YEAST})
            with open(os.path.join(folder, "loci.tsv"), "w") as fh:
                fh.write("locus\tstatus\n2at2759\tMissing\n")
            with open(os.path.join(folder, "rrna.tsv"), "w") as fh:
                fh.write("kind\tcontig\tstart\tend\tstrand\tregions\tlsu_complete\tsequence\tcontig_length\tintrons\n"
                         "rDNA\tc1\t1\t40\t+\t18S:1-20,28S:21-40\tt\t" + "ACGT" * 10 + "\t5000\t-\n")
            store.pack(folder)
            part = package.genome_part((work, YEAST, 1))
            self.assertEqual(part["rrna_rows"], rrna.load(folder))
            self.assertEqual(len(part["rrna_rows"]), 1)
            jobs = {YEAST: dict(accession=YEAST, species_taxid=4932, taxid=559292),
                    "G0": dict(accession="G0", species_taxid=4932)}  # not done
            opts = types.SimpleNamespace(work=work)
            want = eukore_build.genome_rrna_records(opts, jobs, ncbi)
            collected = {"rrna_rows": {YEAST: part["rrna_rows"]}}
            with mock.patch.object(rrna, "load", side_effect=AssertionError("read again")):
                self.assertEqual(eukore_build.genome_rrna_records(opts, jobs, ncbi, collected), want)
            self.assertEqual((len(want), collected), (1, {}))  # taken from collect's: held once

    def test_stamping_an_older_record_keeps_a_packed_genome_packed(self):
        with tempfile.TemporaryDirectory() as work:
            folder = folder_of(work, "G1", {"accession": "G1", "genome_record": 1})
            store.pack(folder)
            stamp = genomes.done_stamp(work, "G1")
            self.assertEqual(genomes.reversion(work, "G1"), 1)
            self.assertEqual((genomes.record_stale(work, "G1"), store.is_packed(folder), os.path.isdir(folder)),
                             (False, True, False))
            self.assertNotEqual(genomes.done_stamp(work, "G1"), stamp)
            self.assertEqual(store.rows(folder, "loci.tsv")[0], ["1at2759", "Single"])  # its tables as they were


class BuildPassesTest(unittest.TestCase):
    """eukore_build.py's passes over the done genomes."""

    def tearDown(self):
        store.close()

    @staticmethod
    def quiet(function, *args):
        with contextlib.redirect_stdout(io.StringIO()) as out:
            return function(*args), out.getvalue()

    def test_the_done_folders_packed(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as work:
            folders = [folder_of(work, f"G{i}") for i in range(5)]
            done = os.path.join(work, "done")
            os.makedirs(os.path.join(done, ".node-1.G7.partial"))  # a run going on: not a done genome
            os.makedirs(os.path.join(done, "G8"))  # no done.json: no genome's results
            opts = types.SimpleNamespace(work=work, done_folders=True, threads=4)
            self.assertEqual(self.quiet(eukore_build.pack_pass, opts), (None, ""))  # --done_folders: as they are
            self.assertTrue(all(os.path.isdir(f) for f in folders))
            opts.done_folders = False
            said = self.quiet(eukore_build.pack_pass, opts)[1]
            self.assertIn("5 to pack", said)
            self.assertIn("5 packed", said)
            self.assertEqual(sorted(os.listdir(done)), [".node-1.G7.partial"] + [f"G{i}.zip" for i in range(5)] + ["G8"])
            self.assertEqual([genomes.is_done(work, f"G{i}") for i in range(5)], [True] * 5)
            self.assertEqual(self.quiet(eukore_build.pack_pass, opts)[1], "")  # nothing left to pack: nothing said

    def test_the_metadata_made_under_the_genomes_claim(self):
        import eukore_build  # noqa: PLC0415
        from eukore import metadata  # noqa: PLC0415
        for pack in (False, True):
            with tempfile.TemporaryDirectory() as tmp:
                work, acc = os.path.join(tmp, "work"), "G1"
                folder = folder_of(work, acc)
                if pack:
                    store.pack(folder)
                genome = os.path.join(tmp, "g.fna.gz")
                with gzip.open(genome, "wt") as fh:
                    fh.write(">c1\n" + "ACGT" * 300 + "\n")
                task = (acc, genome, work, acc, None)
                self.assertEqual(metadata.missing(folder), {metadata.STATS})
                with genomes.claimed(work, acc):  # another build's: left to it
                    with mock.patch.object(genomes, "claim", return_value=None):
                        self.assertEqual(eukore_build.describe_one(task), (acc, None))
                self.assertEqual(metadata.missing(folder), {metadata.STATS})
                self.assertEqual(eukore_build.describe_one(task), (acc, None))
                self.assertEqual(metadata.missing(folder), set())
                stats = json.loads(store.text(folder, metadata.STATS))
                self.assertEqual((stats["length"], stats["accession"]), (1200, acc))
                self.assertTrue(store.has(folder, metadata.ORGANELLE))
                self.assertEqual(store.is_packed(folder), pack)
                self.assertEqual(os.listdir(os.path.join(work, "done")), [acc + (".zip" if pack else "")])
                # made by another build since this one's plan: not made again
                with mock.patch.object(metadata, "describe", side_effect=AssertionError("made again")):
                    self.assertEqual(eukore_build.describe_one(task), (acc, None))
                # the organelle pass' molecules of a packed genome: written out for its tools, once
                opts = types.SimpleNamespace(work=work)
                saved = eukore_build.saved_molecules(opts, tmp, acc)
                self.assertEqual(saved, os.path.join(tmp if pack else folder, f"{acc}." * pack + metadata.ORGANELLE))
                self.assertTrue(os.path.isfile(saved))
                self.assertIsNone(eukore_build.saved_molecules(opts, tmp, "G2"))  # not measured: its genome read

    def test_a_published_read_whose_batch_failed_is_left_out(self):
        # until 2026-11-10 it went into the package as an 18S record of its whole length (a 4-6 kb operon read among
        # the 18S), on the strand it came in
        import eukore_build  # noqa: PLC0415
        reads = [dict(seqid=n, taxid="4932", sequence="ACGT" * 20, assembly_id="-") for n in ("good", "bad")]
        failing = [True]

        def annotate(batch, barrnap=None, itsx=None):
            if failing[0] and any(name == "ROD_bad" for name, _s, _t in batch):
                return {}, "RuntimeError: barrnap failed"
            return {name: ("18S:0-40", "f", "-") for name, _s, _t in batch}, None
        with tempfile.TemporaryDirectory() as work, contextlib.ExitStack() as stack:
            for target, name, value in (
                    (eukore_build.sources, "ready", lambda _src, name: "rod" if name == "ROD" else None),
                    (eukore_build.sources, "rod_operons", lambda _path: iter(reads)),
                    (eukore_build.rrna, "name_index", lambda _ncbi: {}),
                    (eukore_build.rrna, "itsx_groups", lambda _ncbi, _taxid: "all"),
                    (eukore_build.shutil, "which", lambda tool: "/x/" + tool),
                    (eukore_build.util, "unordered", lambda _pool, task, items, *_more: (task(i) for i in items)),
                    (eukore_build, "_annotate_batch", annotate), (eukore_build, "PUBLISHED_BATCH", 1)):
                stack.enter_context(mock.patch.object(target, name, value))
            opts = types.SimpleNamespace(work=work, threads=1, inside=lambda _taxid: True)
            ncbi = mock.Mock(species_of=lambda taxid: taxid)
            got, said = self.quiet(eukore_build.published, opts, "src", ncbi)
            self.assertEqual([(name, row["kind"]) for name, _seq, row in got], [("ROD_good", "rDNA")])
            self.assertIn("1 batches of 1 failed", said)
            self.assertIn("their 1 reads left out of this package", said)
            failing[0] = False  # the next run annotates it: both, the good one from the cache
            got, said = self.quiet(eukore_build.published, opts, "src", ncbi)
            self.assertEqual([name for name, _seq, _row in got], ["ROD_good", "ROD_bad"])
            self.assertIn("published reads: 1 to annotate", said)
            # no barrnap: no read has regions, each an 18S record as before
            stack.enter_context(mock.patch.object(eukore_build.shutil, "which", lambda _tool: None))
            stack.enter_context(mock.patch.object(eukore_build, "read_regions", lambda _path: {}))
            got = self.quiet(eukore_build.published, opts, "src", ncbi)[0]
            self.assertEqual([(name, row["kind"]) for name, _seq, row in got], [("ROD_good", "18S"), ("ROD_bad", "18S")])

    def test_one_fetch_at_a_time(self):
        # until 2026-11-10 nothing stopped two fetch stages on one WORK (a job submitted twice): the downloads have
        # no claims of their own
        import eukore_build  # noqa: PLC0415
        from eukore import util  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            work = os.path.join(tmp, "work")
            os.makedirs(work)
            args = ["--build_dir", work, "--out", os.path.join(tmp, "out"), "--stage", "fetch"]
            lock = util.run_lock(work, "fetch")
            with open(lock, "w") as fh:
                fh.write("another-node\t4242\n")  # a fetch going on there
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as stop, \
                    mock.patch.object(eukore_build, "fetch_sources", side_effect=AssertionError("fetched beside it")):
                eukore_build.main(args)
            self.assertIn("--stage fetch", str(stop.exception))
            self.assertIn("process 4242 on another-node", str(stop.exception))
            old = time.time() - 2 * util.RUN_STALE  # that node died: its lock taken over, let go when this run ends
            os.utime(lock, (old, old))
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(RuntimeError), \
                    mock.patch.object(eukore_build, "fetch_sources", side_effect=RuntimeError("the fetch ran")):
                eukore_build.main(args)
            self.assertFalse(os.path.exists(lock))

    def test_the_python_the_build_needs(self):
        import eukore_build  # noqa: PLC0415
        self.assertIsNone(eukore_build.python_problem())  # the one that runs this test
        self.assertIsNone(eukore_build.python_problem((3, 12, 3), True))
        self.assertIsNone(eukore_build.python_problem((3, 11, 4), True))
        self.assertIn("Python 3.11.2: the build needs tarfile's extraction filter",
                      eukore_build.python_problem((3, 11, 2), False))
        self.assertIn("Python 3.8.17", eukore_build.python_problem((3, 8, 17), True))
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(eukore_build, "python_problem",
                                                                      return_value="Python 3.11.2: too old"), \
                contextlib.redirect_stderr(io.StringIO()) as err, self.assertRaises(SystemExit):
            eukore_build.main(["--build_dir", os.path.join(tmp, "work"), "--out", os.path.join(tmp, "out")])
        self.assertIn("Python 3.11.2: too old", err.getvalue())  # said at the start, before any work


class ClaimTest(unittest.TestCase):
    def test_a_dead_claim_taken_over_by_one_build(self):
        # until 2026-11-10 two builds that both found a claim dead each removed what was there and made their own:
        # the second removed the first's new claim, and both ran the genome
        with tempfile.TemporaryDirectory() as work:
            acc = "GCA_000000001.1"
            path = genomes.claim_path(work, acc)
            os.makedirs(os.path.dirname(path))

            def dead():
                with open(path, "w") as fh:
                    fh.write("some-other-node\t12345\n")
                old = time.time() - 2 * genomes.CLAIM_STALE
                os.utime(path, (old, old))
            dead()
            with open(path + ".take", "w"):  # another build is taking it over right now
                pass
            self.assertIsNone(genomes.claim(work, acc))
            with open(path) as fh:
                self.assertEqual(fh.read(), "some-other-node\t12345\n")  # left to that build
            old = time.time() - 2 * genomes.TAKE_STALE  # that build was killed in the act: its lock removed, for the
            os.utime(path + ".take", (old, old))        # next try
            self.assertIsNone(genomes.claim(work, acc))
            self.assertEqual(genomes.claim(work, acc), path)
            self.assertTrue(genomes.claim_own(path))
            self.assertFalse(os.path.exists(path + ".take"))
            # many at once: one takes it
            dead()
            got = []
            threads = [threading.Thread(target=lambda: got.append(genomes.claim(work, acc))) for _ in range(8)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
            self.assertEqual(sorted(g is not None for g in got), [False] * 7 + [True])

    def test_a_claim_taken_over_is_not_its_old_owners(self):
        # a run stalled longer than CLAIM_STALE, its claim taken over by another build: its end leaves that build's
        # claim (until 2026-11-10 it removed it, a third build then ran the genome beside the second)
        with tempfile.TemporaryDirectory() as work:
            acc = "GCA_000000001.1"
            with genomes.claimed(work, acc) as held:
                with open(held, "w") as fh:
                    fh.write("some-other-node\t12345\n")
            self.assertTrue(os.path.isfile(genomes.claim_path(work, acc)))

    def test_superseded_written_whole(self):
        with tempfile.TemporaryDirectory() as work:
            def plan(**accession_of):
                return genomes.supersede(work, {a: dict(species_taxid=sp) for sp, a in accession_of.items()})
            self.assertEqual(plan(s1="A1", s2="B1"), [])
            self.assertEqual(plan(s1="A2", s2="B1"), [("A1", "A2")])
            table = os.path.join(work, genomes.SUPERSEDED)
            with open(table, "a") as fh:  # a build killed in its line, as the appending of before left it
                fh.write("B1")
            self.assertEqual(plan(s1="A2", s2="B2"), [("B1", "B2")])
            with open(table) as fh:
                self.assertEqual(fh.read(), "old\tnew\nA1\tA2\nB1\tB2\n")
            self.assertEqual([n for n in os.listdir(work) if n.endswith(".partial")], [])


if __name__ == "__main__":
    unittest.main()
