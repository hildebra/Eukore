#!/usr/bin/env python3
"""Tests of Eukore's core modules: the NCBI taxonomy and its subset (written into a package and read back the same for
the package's taxa), NCBI's assembly tables read, the model species table, a genome's nuclear molecules (its
rRNA's).

    python3 -m unittest discover -s tests -t tests
"""

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi, rrna, sources, species, util  # noqa: E402
import fixtures as fx  # noqa: E402


class CoreTest(unittest.TestCase):
    def test_taxonomy_subset(self):
        with tempfile.TemporaryDirectory() as tmp:
            full = ncbi.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            self.assertEqual(full.species_of(559292), 4932)
            n = full.write_subset(os.path.join(tmp, "taxonomy"), [559292, 28985])
            sub = ncbi.Ncbi(os.path.join(tmp, "taxonomy"))
            self.assertEqual(n, len(set(full.ancestors(559292)) | set(full.ancestors(28985))))
            for taxid in (559292, 4932, 28985):
                self.assertEqual((sub.species_of(taxid), sub.ranks(taxid)), (full.species_of(taxid), full.ranks(taxid)))
            self.assertEqual(sub.synonyms[28985], ["Candida sphaerica"])
            self.assertNotIn(40410, sub.name)  # not on the lineages asked for
            # the genetic code, NCBI's nodes.dmp column 7, kept in the subset
            self.assertEqual((full.code_of(28985), sub.code_of(28985), sub.code_of(4932)), (6, 6, 1))

    def test_assembly_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            by = {r["accession"]: r for r in rows}
            self.assertNotIn("GCA_000000560.1", by)  # a gene's record as a "complete genome": no genome
            self.assertEqual(by["GCA_000000011.1"]["flags"], ["contaminated"])
            self.assertEqual(by["GCA_000000012.1"]["category"], "derived from metagenome")
            self.assertEqual((by["GCA_000146045.2"]["refseq_category"], by["GCA_000146045.2"]["paired"],
                              by["GCA_000146045.2"]["group"]), ("reference genome", "GCF_000146045.2", "fungi"))

    def test_the_model_species_table(self):
        rows = species.read_table()
        self.assertIn("Homo sapiens", {r["species"] for r in species.eukaryotes(rows)})
        self.assertTrue(os.path.isfile(species.UNITS))
        self.assertFalse(os.path.exists(os.path.join(sources.DATA, "euk30_loci.tsv")))  # the user picks loci

    def test_the_genome_rrna_reads_nuclear_molecules_only(self):
        # Arabidopsis' assembly holds its mitochondrion and chloroplast, whose SSU barrnap calls an 18S (2026-11-06)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "g.fna")
            with open(path, "w") as fh:
                fh.write(">CP002684.1 Arabidopsis thaliana chromosome 1 sequence\nACGT\nAC\n"
                         ">BK010421.1 TPA_asm: Arabidopsis thaliana ecotype Col-0 mitochondrion, complete genome\nGGGG\n"
                         ">AP000423.1 Arabidopsis thaliana chloroplast DNA, complete genome\nTTTT\n"
                         ">x1 Plasmodium falciparum apicoplast\nCCCC\n>CP002685.1 chromosome 2\nAAAA\n"
                         ">CP000881.1 Guillardia theta nucleomorph chromosome 1\nGGCC\n")
            self.assertEqual(list(rrna.contigs(path, nuclear=True)), [("CP002684.1", "ACGTAC"), ("CP002685.1", "AAAA")])
            self.assertEqual(len(list(rrna.contigs(path))), 6)

    def test_lineages_agree_by_taxid(self):
        # two homonym genera of two kingdoms (Morus the plant, the bird) share a name, not a taxon
        from eukore import qc  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            tax = ncbi.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        tax.parent.update({7001: 4893, 7002: 7001, 7003: 6231, 7004: 7003})
        tax.rank.update({7001: "genus", 7002: "species", 7003: "genus", 7004: "species"})
        tax.name.update({7001: "Morus", 7002: "Morus alba", 7003: "Morus", 7004: "Morus bassanus"})
        self.assertEqual(qc.depth(tax, 7002, 7004), 0)
        self.assertEqual(qc.depth(tax, 4932, 28985), 4)
        # a taxid the taxonomy lacks is left out of an LCA, not making it none
        self.assertEqual(qc.lca(tax, {4932, 27291, 99999999}), 4930)
        self.assertIsNone(qc.lca(tax, {99999999}))


class FetchTest(unittest.TestCase):
    def test_scratch_and_a_checker_process(self):
        # a whole gzip file and one cut short (its announced length right): the check, in a process of its own,
        # keeps the cut one out; the .part files stay in the scratch folder, none is left anywhere
        import concurrent.futures, functools, gzip, http.server, multiprocessing, threading  # noqa: E401, PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            served, out, scratch = (os.path.join(tmp, d) for d in ("served", "out", "scratch"))
            for d in (served, out, scratch):
                os.makedirs(d)
            whole = gzip.compress(os.urandom(1 << 16))
            with open(os.path.join(served, "whole.fna.gz"), "wb") as fh:
                fh.write(whole)
            with open(os.path.join(served, "cut.fna.gz"), "wb") as fh:
                fh.write(whole[:len(whole) // 2])
            quiet = type("Quiet", (http.server.SimpleHTTPRequestHandler,), {"log_message": lambda *a: None})
            server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(quiet, directory=served))
            serving = threading.Thread(target=server.serve_forever, daemon=True)
            serving.start()
            url = f"http://127.0.0.1:{server.server_address[1]}/"
            curl = util.CURL
            try:
                for util.CURL in {curl, None}:  # by curl (where there) and by urllib
                    with concurrent.futures.ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("spawn")) \
                            as checker:
                        times, counter = {}, util.collections.Counter
                        got = {name: util.fetch_file(url + name, os.path.join(out, name), attempts=1, scratch=scratch,
                                                     checker=checker, times=times.setdefault(name, counter()))
                               for name in ("whole.fna.gz", "cut.fna.gz", "missing.fna.gz")}
                    self.assertEqual(times["whole.fna.gz"]["bytes"], len(whole))  # the file placed, its steps timed
                    self.assertTrue(all(times["whole.fna.gz"][k] > 0 for k in ("transfer", "check", "copy")))
                    self.assertEqual(times["cut.fna.gz"]["bytes"], 0)
                    self.assertEqual(got["whole.fna.gz"], "", util.CURL)
                    self.assertTrue(got["cut.fna.gz"].startswith(("EOFError", "OSError")), got)  # Python's; a tool's
                    self.assertEqual(got["missing.fna.gz"], "not found")
                    self.assertEqual(sorted(os.listdir(out)), ["whole.fna.gz"])
                    self.assertEqual(os.listdir(scratch), [])
                    with open(os.path.join(out, "whole.fna.gz"), "rb") as fh:
                        self.assertEqual(fh.read(), whole)
                    os.remove(os.path.join(out, "whole.fna.gz"))
            finally:
                util.CURL = curl
                server.shutdown()
                server.server_close()
                serving.join()


    def test_a_stop_ends_the_retries(self):
        # a stop (Ctrl-C: eukore_build.at_a_time sets util.STOP) while a download fails: not tried again
        from unittest import mock  # noqa: PLC0415
        tries = []

        def failing(url, path):
            tries.append(url)
            util.STOP.set()
            raise OSError("curl exit 2: interrupted")
        try:
            with mock.patch.object(util, "CURL", "curl"), mock.patch.object(util, "curl_get", failing), \
                    mock.patch.object(util, "RETRY_WAIT", 0):
                with tempfile.TemporaryDirectory() as tmp:
                    why = util.fetch_file("http://x/a.fna.gz", os.path.join(tmp, "a.fna.gz"))
        finally:
            util.STOP.clear()
        self.assertEqual((len(tries), why), (1, "OSError: curl exit 2: interrupted"))

    def test_gzip_checks(self):
        # igzip's, pigz's (where there) and Python's checks: a whole file passes, one cut short or not gzip fails
        import gzip, shutil  # noqa: E401, PLC0415
        tools = [t for t in (shutil.which("igzip"), shutil.which("pigz")) if t] + [None]
        with tempfile.TemporaryDirectory() as tmp:
            whole = gzip.compress(os.urandom(1 << 16))
            files = {}
            for name, data in (("whole", whole), ("cut", whole[:len(whole) // 2]), ("bad", b"not gzip")):
                files[name] = os.path.join(tmp, name + ".part")
                with open(files[name], "wb") as fh:
                    fh.write(data)
            for tool in tools:
                util.gzip_whole(files["whole"], tool)
                for name in ("cut", "bad"):
                    with self.assertRaises((OSError, EOFError), msg=f"{tool} {name}"):
                        util.gzip_whole(files[name], tool)


class ProgressTest(unittest.TestCase):
    def test_reports_every_interval_with_the_time_left(self):
        # four genomes of 1, 1, 2 and 4 Gb; a report every 600 s; a Gb takes 100 s
        now, said = [0.0], []
        p = util.Progress("genomes fetched", [1e9, 1e9, 2e9, 4e9], said.append, every=600, clock=lambda: now[0])
        now[0] = 100.0
        p.step(1e9)
        self.assertEqual(said, [])  # not 600 s yet
        now[0] = 700.0  # a slow one
        p.step(1e9, failed=True)
        self.assertEqual(len(said), 1)
        self.assertIn("2/4 genomes fetched (1 failed), 2 left (6.0 Gb)", said[0])
        self.assertIn("11m 40s so far, ~35m 00s left", said[0])  # 2 Gb in 700 s: 6 Gb in 2,100 s, by size not count
        now[0] = 900.0
        p.step(2e9, skipped=True)  # not fetched yet: no work, out of the estimate
        self.assertEqual(len(said), 1)
        now[0] = 1000.0
        p.step(4e9, note="; done")
        self.assertEqual(said[1], "  4/4 genomes fetched (1 failed), 0 left; took 16m 40s; done")
        self.assertEqual([util.duration(s) for s in (45, 725, 11220, 187200)], ["45s", "12m 05s", "3h 07m", "2d 04h"])
        # no sizes: each item weighs the same
        q = util.Progress("genomes run", [None, None], said.append, every=0, clock=lambda: now[0])
        now[0] += 10
        q.step(None)
        self.assertIn("~10s left", said[-1])

    def test_early_reports_then_every_interval(self):
        # the downloads' first reports at 1, 2, 3, 5 and 8 minutes, then every 10; an item finishing every 30 s
        now, said = [0.0], []
        p = util.Progress("genomes fetched", [1] * 100, said.append, every=600, clock=lambda: now[0],
                          early=util.REPORT_EARLY)
        at = []
        for _ in range(60):
            now[0] += 30
            n = len(said)
            p.step(1)
            if len(said) > n:
                at.append(now[0])
        self.assertEqual(at, [60, 120, 180, 300, 480, 1080, 1680])


if __name__ == "__main__":
    unittest.main()
