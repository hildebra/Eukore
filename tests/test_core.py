#!/usr/bin/env python3
"""Tests of Eukore's core modules: the NCBI taxonomy and its subset (written into a package and read back the same for
the package's taxa), NCBI's assembly tables read, the model species table, a genome's nuclear molecules (its
rRNA's).

    python3 -m unittest discover -s tests -t tests
"""

import os
import signal
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi, rrna, sources, species, util  # noqa: E402
import fixtures as fx  # noqa: E402


def _slow(x):  # a pool's task (util.unordered)
    time.sleep(0.2)
    return 2 * x


def _killed(x):  # a pool's task whose worker is killed, as by the OOM killer
    if x == 2:
        os.kill(os.getpid(), signal.SIGKILL)
    time.sleep(0.2)
    return x


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
            # the mitochondrial code (column 9) and the plastid code (11: this taxdump has no plastid column)
            for t in (full, sub):
                self.assertEqual((t.code_of(4932, "mito"), t.code_of(559292, "mito"), t.code_of(4932, "plastid")),
                                 (3, 1, 11))

    def test_genetic_codes(self):
        from eukore import gencode
        self.assertEqual(len(gencode.TABLES), 27)
        self.assertTrue(all(len(t) == 64 for t in gencode.TABLES.values()))
        self.assertEqual(gencode.translate("AGAAGGATATGA", 2), "**MW")  # vertebrate mitochondria
        self.assertEqual(gencode.translate("TGA", 4), "W")
        self.assertEqual(gencode.translate("TAATAG", 6), "QQ")  # ciliates
        self.assertEqual(gencode.translate("CTG", 12), "S")  # the CTG clade's yeasts
        self.assertEqual(gencode.translate("CTGTGA", 3), "TW")  # yeast mitochondria
        self.assertEqual(gencode.translate("ATGTAAAAA", 99), "M*K")  # an unknown table reads as the standard one
        self.assertEqual(gencode.protein_of("ATGTAAAAATAA"), ("M*K", 1))  # the last stop dropped, internal kept

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
        # (extract_genome's title test; until 2026-11-08 through rrna.contigs, which nothing else read)
        titles = ["CP002684.1 Arabidopsis thaliana chromosome 1 sequence",
                  "BK010421.1 TPA_asm: Arabidopsis thaliana ecotype Col-0 mitochondrion, complete genome",
                  "AP000423.1 Arabidopsis thaliana chloroplast DNA, complete genome",
                  "x1 Plasmodium falciparum apicoplast", "CP002685.1 chromosome 2",
                  "CP000881.1 Guillardia theta nucleomorph chromosome 1"]
        self.assertEqual([t.split()[0] for t in titles if not rrna.NOT_NUCLEAR.search(t)], ["CP002684.1", "CP002685.1"])

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

    def test_an_api_key_stays_on_its_host(self):
        # a redirect elsewhere, or off https, drops the request's Authorization (urllib's own passes it on); VEuPathDB's
        # key goes to its sites alone
        import urllib.request  # noqa: PLC0415
        from eukore import annotation  # noqa: PLC0415
        handler = util.KeepAuthorization()

        def redirected(to):
            req = urllib.request.Request("https://plasmodb.org/a/service/raw-files/x.gff",
                                         headers={"Authorization": "Bearer secret"})
            return handler.redirect_request(req, None, 302, "Found", {}, to).has_header("Authorization")
        self.assertTrue(redirected("https://plasmodb.org/elsewhere/x.gff"))
        self.assertFalse(redirected("https://cdn.example.com/x.gff"))
        self.assertFalse(redirected("http://plasmodb.org/x.gff"))
        self.assertTrue(annotation.veupathdb_host("https://toxodb.org/a/service/raw-files/x.gff"))
        self.assertFalse(annotation.veupathdb_host("http://toxodb.org/x.gff"))
        self.assertFalse(annotation.veupathdb_host("https://toxodb.org.example.com/x.gff"))

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

    def test_reports_while_nothing_finishes_with_the_resources(self):
        # a pool whose first result comes 25 minutes in: reports at 10 and 20 minutes all the same, each with usage
        now, said = [0.0], []
        p = util.Progress("genomes run", [3e9, 1e9], said.append, every=600, clock=lambda: now[0],
                          usage=lambda: "; 4.0 cores busy")

        for minute in range(1, 26):  # util.unordered ticks at each poll where nothing finished
            now[0] = 60.0 * minute
            p.tick(", 2 running")
        p.step(3e9)
        p.step(1e9)
        self.assertEqual(said[0], "  0/2 genomes run (0 failed), 2 left (4.0 Gb); 10m 00s so far, 2 running; 4.0 cores "
                                  "busy")
        self.assertTrue(said[1].startswith("  0/2 genomes run (0 failed), 2 left (4.0 Gb); 20m 00s so far"))
        self.assertTrue(said[-1].startswith("  2/2 genomes run (0 failed), 0 left; took"))
        # no unit: the items left alone
        q = util.Progress("genomes compared", [1, 1], said.append, every=0, unit="", scale=1, clock=lambda: now[0])
        q.step(1)
        self.assertRegex(said[-1], r"^  1/2 genomes compared \(0 failed\), 1 left; ")

    def test_a_pool_ticks_while_it_waits_and_stops_where_a_worker_dies(self):
        ticks = []

        class Ticks:
            def tick(self, note=""):
                ticks.append(note)
        with util.process_pool(2) as pool:
            got = sorted(util.unordered(pool, _slow, [3, 1, 2, 4, 5], Ticks(), lambda: ", 2 running", poll=0.05))
        self.assertEqual(got, [2, 4, 6, 8, 10])
        self.assertTrue(ticks and ticks[0] == ", 2 running")
        # a worker killed (the OOM killer): the build stops, where multiprocessing.Pool waited for it forever
        with self.assertRaises(SystemExit):
            with util.process_pool(2) as pool:
                list(util.unordered(pool, _killed, [1, 2, 3], poll=0.05))

    def test_usage_of_the_process_tree(self):
        tick, page = os.sysconf("SC_CLK_TCK"), os.sysconf("SC_PAGE_SIZE")
        me = os.getpid()
        with tempfile.TemporaryDirectory() as proc:
            def stat(pid, ppid, cpu, pages, name="python"):
                os.makedirs(os.path.join(proc, str(pid)), exist_ok=True)
                fields = (["S", str(ppid)] + ["0"] * 9 + [str(cpu), "0", "0", "0"] + ["0"] * 6 + [str(pages)]
                          + ["0"] * 20)
                with open(os.path.join(proc, str(pid), "stat"), "w") as fh:
                    fh.write(f"{pid} ({name}) {' '.join(fields)}\n")
            stat(me, 1, 10 * tick, 1000)
            stat(me + 1, me, 5 * tick, 500, name="mini prot) (x")  # a child; a name with spaces and brackets
            stat(me + 2, me + 1, 5 * tick, 500)  # a grandchild
            stat(me + 3, 1, 999 * tick, 99999)  # another user's
            with open(os.path.join(proc, "meminfo"), "w") as fh:
                fh.write("MemTotal:       1000000 kB\nMemFree:  10 kB\nMemAvailable:    250000 kB\n")
            now = [0.0]
            u = util.Usage([("WORK", proc), ("OUT", proc)], cores=8, clock=lambda: now[0], proc=proc)
            self.assertEqual(u.tree(), (20.0, 2000 * page))
            stat(me, 1, 40 * tick, 1000)  # 30 CPU seconds more in 10 s: 3 cores
            now[0] = 10.0
            line = u.line()
            self.assertIn("; 3.0 cores busy (of 8), ", line)
            self.assertIn(f"resident (peak {util.size_text(2000 * page)})", line)
            self.assertIn("; node: 256 MB of 1.0 GB free", line)
            self.assertIn("; disk free: WORK, OUT ", line)  # one file system: one figure
        self.assertEqual([util.size_text(n) for n in (5e8, 8.41e10, 2.1e12)], ["500 MB", "84.1 GB", "2.1 TB"])
        # no /proc: the disk alone
        self.assertEqual(util.Usage(proc="/no/such/proc").line(), "")


if __name__ == "__main__":
    unittest.main()
