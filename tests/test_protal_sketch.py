#!/usr/bin/env python3
"""Tests of protal's genome sketch of the package's species (eukore/protal_sketch.py, eukore_build.py --protal; since
2026-11-19) with a stand-in protal: the genome table (each species' representative first, keyed by its accession; its
further assemblies on disk after it; transcriptomes, EukProt's records and the genomes whose copies the package leaves
out not taken), the copies' spans, the command run (the entries scaled by genome size), the rerun cache (reused;
made again; a stale sketch kept where genomes of its table were pruned), a failure, the sketch in the package's
manifest, and its release download of its own.

    python3 -m unittest discover -s tests -t tests
"""

import contextlib
import gzip
import io
import json
import os
import sys
import tempfile
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, ncbi as ncbi_mod, package, protal_sketch as ps, release, store, util  # noqa: E402
import fixtures as fx  # noqa: E402

# protal's stand-in: --version; sketch: its arguments noted (one JSON line a call), a file of the table's species and
# --kmers written to --out, the summary lines protal prints
PROTAL = r'''#!/usr/bin/env python3
import json, sys
a = sys.argv[1:]
if a == ["--version"]:
    print("protal 9.9.9 (stand-in)")
    sys.exit(0)
with open(CALLS, "a") as fh:
    fh.write(json.dumps(a) + "\n")
table = [line.rstrip("\n").split("\t") for line in open(a[a.index("--genomes") + 1]) if line.strip()]
species = sorted({k for k, _p in table})
kmers = a[a.index("--kmers") + 1]
with open(a[a.index("--out") + 1], "wb") as fh:
    fh.write(("PSK " + kmers + " " + ",".join(species)).encode())
print("protal sketch: %d genomes of %d species; %s entries a species" % (len(table), len(species), kmers))
print("protal sketch: %d species, %d in the sketch (0 species with fewer than %s own windows)"
      % (len(species), 2 * len(species), kmers))
'''
BROKEN = "#!/bin/sh\necho 'protal sketch: cannot read the genome /x.fna.gz' >&2\nexit 3\n"
SILENT = "#!/bin/sh\nexit 0\n"  # ends well and writes no sketch

YEAST, OTHER, PANEL, GONE = "GCA_000146045.2", "GCA_000000003.1", "GCA_000000004.1", "GCA_000000005.1"
ASCARIS, TSA, POOR, PRUNED, EUKPROT = "GCA_000951015.1", "GAAA01000000", "GCA_000000009.1", "GCA_000000010.1", \
    "EUK_1"


def row(acc, sp, source, rep):
    return dict(accession=acc, species_taxid=sp, taxid=sp, organism="x", source=source,
                representative="t" if rep else "f", loci=0)


GENOMES = {YEAST: row(YEAST, 4932, "ncbi", True), OTHER: row(OTHER, 4932, "assembly", False),
           PANEL: row(PANEL, 4932, "intraspecific", False), GONE: row(GONE, 4932, "assembly", False),
           ASCARIS: row(ASCARIS, 6253, "model", True), TSA: row(TSA, 28985, "tsa", True),
           POOR: row(POOR, 5476, "ncbi", True), PRUNED: row(PRUNED, 9606, "ncbi", True),
           EUKPROT: row(EUKPROT, 952, "no_assembly", True)}
ON_DISK = (YEAST, OTHER, PANEL, ASCARIS, TSA, POOR, EUKPROT)  # GONE and PRUNED not


def loci_row(locus, status, contig="-", start="-", end="-", alt=None):
    values = dict.fromkeys(genomes.LOCI_COLUMNS, "-")
    values.update(locus=locus, status=status, contig=contig, start=str(start), end=str(end), strand="+")
    return ([] if alt is None else [str(alt)]) + [values[c] for c in genomes.LOCI_COLUMNS]


def done(work, acc, rows, alternates=(), pack=False):
    folder = genomes.done_dir(work, acc)
    os.makedirs(folder)
    with open(os.path.join(folder, "loci.tsv"), "w", newline="\n") as fh:
        fh.write("\t".join(genomes.LOCI_COLUMNS) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    if alternates:
        with open(os.path.join(folder, genomes.ALTERNATES[0]), "w", newline="\n") as fh:
            fh.write("\t".join(("alt",) + genomes.LOCI_COLUMNS) + "\n"
                     + "".join("\t".join(r) + "\n" for r in alternates))
    if pack:
        store.pack(folder)


class World:
    """A build's genomes folder, WORK with done genomes, and protal's stand-ins."""

    def __init__(self, tmp):
        self.tmp, self.work, self.folder = tmp, os.path.join(tmp, "work"), os.path.join(tmp, "genomes")
        os.makedirs(self.folder)
        for acc in ON_DISK:
            with gzip.open(genomes.genome_file(self.folder, acc), "wt") as fh:
                fh.write(f">{acc}_c1 a contig\nACGT\n")
        done(self.work, YEAST, [loci_row("1at2759", "Single", "chrI", 100, 900),
                                loci_row("2at2759", "Missing"),
                                loci_row("3at2759", "Duplicated", "chrII", 5000, 7200),
                                loci_row("RPB1", "Single", "chrIV", 20, 4020)],
             [loci_row("3at2759", "Duplicated", "chrX", 300, 2500, alt=1)])
        done(self.work, PANEL, [loci_row("1at2759", "Single", "scf7", 1, 800)])
        done(self.work, ASCARIS, [loci_row("1at2759", "Fragmented", "AL1", 10, 50)], pack=True)
        self.calls = os.path.join(tmp, "calls.jsonl")
        self.protal = fx.tool(tmp, "protal", PROTAL.replace("CALLS", repr(self.calls)))
        self.broken = fx.tool(tmp, "protal_broken", BROKEN)
        self.silent = fx.tool(tmp, "protal_silent", SILENT)
        self.said = []

    def path(self, acc):
        return genomes.genome_file(self.folder, acc)

    def build(self, protal=None, kmers=ps.KMERS, mb=ps.MB, exponent=ps.EXPONENT):
        return ps.build(protal or self.protal, GENOMES, self.path, self.work, kmers=kmers, mb=mb, threads=4,
                        left_out=[POOR], log=self.said.append, exponent=exponent)

    def ran(self):
        if not os.path.isfile(self.calls):
            return []
        with open(self.calls) as fh:
            return [json.loads(line) for line in fh]

    def sketch_folder(self, name=""):
        return os.path.join(self.work, ps.FOLDER, name)


class ProtalSketchTest(unittest.TestCase):
    def test_the_table_positions_and_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            path, about = w.build()
            folder = w.sketch_folder()
            self.assertEqual(path, os.path.join(folder, ps.SKETCH))
            # the representative first, keyed by its accession; its further assemblies on disk after it (the one
            # pruned left out); the Ascaris model species; not the transcriptome, EukProt's records, the poor genome
            # (its copies left out), nor the species whose genome is gone
            with open(os.path.join(folder, ps.TABLE)) as fh:
                table = [line.rstrip("\n").split("\t") for line in fh]
            self.assertEqual(table, [[YEAST, w.path(YEAST)], [YEAST, w.path(OTHER)], [YEAST, w.path(PANEL)],
                                     [ASCARIS, w.path(ASCARIS)]])
            # each copy's span (the chosen ones and a Duplicated locus' alternate; Missing left out), of the done
            # genomes in either form; the panel's assembly's too; OTHER has no done genome
            with open(os.path.join(folder, ps.POSITIONS)) as fh:
                positions = [line.rstrip("\n").split("\t") for line in fh]
            self.assertEqual(positions, [list(ps.POSITION_COLUMNS),
                                         ["chrI", "100", "900", YEAST, "1at2759"],
                                         ["chrII", "5000", "7200", YEAST, "3at2759"],
                                         ["chrIV", "20", "4020", YEAST, "RPB1"],
                                         ["chrX", "300", "2500", YEAST, "3at2759"],
                                         ["scf7", "1", "800", PANEL, "1at2759"],
                                         ["AL1", "10", "50", ASCARIS, "1at2759"]])
            # protal's command: exactly these arguments, the sketch written beside and renamed once whole
            self.assertEqual(w.ran(), [["sketch", "--genomes", os.path.join(folder, ps.TABLE),
                                        "--out", os.path.join(folder, ps.SKETCH + ps.PARTIAL),
                                        "--gene_positions", os.path.join(folder, ps.POSITIONS),
                                        "--kmers", "5000", "--kmers_mb", "4", "--kmers_exponent", "0.25",
                                        "--kmers_min", "5000", "--kmers_max", "0", "-t", "4"]])
            self.assertFalse(os.path.exists(os.path.join(folder, ps.SKETCH + ps.PARTIAL)))
            with open(path) as fh:
                self.assertEqual(fh.read(), f"PSK 5000 {YEAST},{ASCARIS}")
            self.assertEqual(about, dict(file=ps.IN_PACKAGE, protal=os.path.realpath(w.protal),
                                         version="protal 9.9.9 (stand-in)", kmers=5000, kmers_mb=4, kmers_exponent=0.25,
                                         kmers_min=5000, kmers_max=0, species=2, genomes=4,
                                         gene_spans=6, made=about["made"], state="made"))
            # the log: the species without their genome file, the summary line of protal's echoed
            said = "\n".join(w.said)
            self.assertIn(f"1 species of the package without their genome file (pruned?), not in it (e.g. "
                          f"['{PRUNED}'])", said)
            self.assertIn("  protal sketch: 2 species, 4 in the sketch (0 species with fewer than 5000 own windows)",
                          w.said)
            self.assertIn("6 BUSCO and marker copies' spans left out", said)
            self.assertIn("1 further assemblies without (none called)", said)

    def test_a_representative_without_copies_is_said(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            store.pack(genomes.done_dir(w.work, YEAST))  # either form read
            os.remove(store.packed(genomes.done_dir(w.work, ASCARIS)))
            _path, about = w.build()
            self.assertEqual(about["gene_spans"], 5)
            self.assertIn(f"no copies known (no done loci.tsv) of 1 representatives (e.g. ['{ASCARIS}']), the sketch "
                          "made without their spans", "\n".join(w.said))

    def test_the_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            first, made = w.build()
            # the same protal, --kmers, genomes and copies: kept, protal not run again
            again, about = w.build()
            self.assertEqual((again, len(w.ran()), about["state"], about["made"]), (first, 1, "reused", made["made"]))
            self.assertFalse([f for f in os.listdir(w.sketch_folder()) if f.endswith(ps.NEXT)])
            # another --kmers: made again
            _p, about = w.build(kmers=500)
            self.assertEqual((len(w.ran()), w.ran()[-1][w.ran()[-1].index("--kmers") + 1], about["state"]),
                             (2, "500", "made"))
            # another genome size and power for those entries: made again, at least --kmers, no upper bound
            _p, about = w.build(kmers=500, mb=100, exponent=0.5)
            self.assertEqual((len(w.ran()), w.ran()[-1][-12:], about["state"], about["kmers_exponent"]),
                             (3, ["--kmers", "500", "--kmers_mb", "100", "--kmers_exponent", "0.5", "--kmers_min", "500",
                                  "--kmers_max", "0", "-t", "4"], "made", 0.5))
            _p, about = w.build(kmers=500)  # and back: made again
            self.assertEqual(len(w.ran()), 4)
            # a new genome on disk (a further assembly fetched since): made again with it
            with gzip.open(w.path(GONE), "wt") as fh:
                fh.write(">g\nACGT\n")
            _p, about = w.build(kmers=500)
            self.assertEqual((len(w.ran()), about["genomes"]), (5, 5))
            # pruned since (genomes of its table gone): kept though stale, and said why; protal not run
            os.remove(w.path(OTHER))
            os.remove(w.path(YEAST))
            w.said.clear()
            path, about = w.build(kmers=500)
            self.assertEqual((path, len(w.ran()), about["state"], about["genomes"], about["kmers"]),
                             (os.path.join(w.sketch_folder(), ps.SKETCH), 5, "stale", 5, 500))
            self.assertEqual(about["stale"], "another genome table, other gene positions; 2 genomes of its table gone")
            self.assertTrue(any(s.startswith("WARNING: protal sketch: ") and "kept though stale" in s
                                for s in w.said), w.said)
            with open(os.path.join(w.sketch_folder(), ps.TABLE)) as fh:  # the kept sketch's table, as it was
                self.assertEqual(len(fh.readlines()), 5)
            # no sketch kept: made of the genomes there, the species without one said
            for name in (ps.SKETCH, ps.KEY):
                os.remove(w.sketch_folder(name))
            w.said.clear()
            _p, about = w.build(kmers=500)
            self.assertEqual((len(w.ran()), about["state"], about["species"], about["genomes"]), (6, "made", 1, 1))
            self.assertIn(f"2 species of the package without their genome file (pruned?), not in it (e.g. "
                          f"['{PRUNED}', '{YEAST}'])", "\n".join(w.said))

    def test_another_binary_is_made_again(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            w.build()
            other = fx.tool(tmp, "protal2", PROTAL.replace("CALLS", repr(w.calls)))
            _p, about = w.build(protal=other)
            self.assertEqual((len(w.ran()), about["state"], about["protal"]), (2, "made", os.path.realpath(other)))

    def test_a_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            w.build()
            # protal fails: said with its output, no sketch (the kept one, of other inputs, removed), no partial
            for bad in (w.broken, w.silent):
                w.said.clear()
                self.assertEqual(w.build(protal=bad), (None, None))
                folder = w.sketch_folder()
                self.assertFalse([f for f in os.listdir(folder) if f in (ps.SKETCH, ps.KEY) or ps.PARTIAL in f])
                said = [s for s in w.said if s.startswith("WARNING: protal sketch failed")]
                self.assertEqual(len(said), 1, w.said)
                self.assertIn("exit 3" if bad == w.broken else "no protal_sketch.bin.partial written", said[0])
            with open(w.sketch_folder(ps.LOG)) as fh:
                self.assertEqual(fh.read(), "")  # the silent one's
            self.assertEqual(w.build(protal=w.broken)[0], None)
            self.assertIn("cannot read the genome /x.fna.gz", w.said[-1])
            # nothing on disk: none made, protal not run
            for acc in ON_DISK:
                os.remove(w.path(acc))
            n = len(w.ran())
            self.assertEqual(w.build(), (None, None))
            self.assertEqual(len(w.ran()), n)
            self.assertIn("no genome of the package's species on disk: none made", w.said[-1])

    def test_into_the_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            path, about = w.build()
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            collected = dict(genomes={YEAST: dict(GENOMES[YEAST], fetched=YEAST, taxid=559292)},
                             copies={YEAST: {}}, features=[], copy_rows=[], taxids=set())
            out = package.write(os.path.join(tmp, "out", "Eukore_2026.11"), "2026.11", ncbi, collected, [],
                                inputs={"protal_sketch": about}, writers=[lambda part: ps.into_package(path, part)])
            self.assertIsNone(package.check(out, deep=True))
            manifest = package.read(out)
            inside = os.path.join(out, ps.IN_PACKAGE)
            self.assertEqual(manifest["files"][ps.IN_PACKAGE], [os.path.getsize(path), util.md5(path)])
            self.assertEqual(manifest["inputs"]["protal_sketch"], about)
            self.assertTrue(os.path.samefile(inside, path))  # a hard link: no second copy
            self.assertIsNone(release.bundle_of(ps.IN_PACKAGE))  # in no bundle: a download of its own

    def test_a_release_download_of_its_own(self):
        # a few GB for all eukaryotes: not in the core, a download of its own as the sourmash sketches are, a hard link
        import tarfile  # noqa: PLC0415
        import test_release  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            pkg = test_release.toy_package(os.path.join(tmp, "out", "Eukore_2026.11"))
            with open(os.path.join(pkg, ps.IN_PACKAGE), "wb") as fh:
                fh.write(b"PSK not BGZF")
            out, said = os.path.join(tmp, "out", "release"), []
            release.write(pkg, out, log=said.append)
            folder = os.path.join(out, "downloads", "Eukore_2026.11")
            with open(os.path.join(out, "data", "downloads.json")) as fh:
                files = {e["key"]: e for e in json.load(fh)["files"]}
            e = files["protal_sketch"]
            path = os.path.join(folder, "Eukore_2026.11_protal_sketch.bin")
            self.assertEqual((e["file"], e["bytes"], e["md5"], e["files"]),
                             (os.path.basename(path), 12, util.md5(path), None))
            self.assertTrue(os.path.samefile(path, os.path.join(pkg, ps.IN_PACKAGE)))
            self.assertEqual(list(files).index("protal_sketch"), list(files).index("sketches") + 1)
            with tarfile.open(os.path.join(folder, "Eukore_2026.11_core.tar.gz"), "r:gz") as tar:
                self.assertNotIn(f"Eukore_2026.11/{ps.IN_PACKAGE}", tar.getnames())
            with open(os.path.join(folder, "MD5SUMS")) as fh:
                self.assertIn(f"{e['md5']}  {e['file']}\n", fh.read())
            # a rerun keeps it (no BGZF end block asked of it)
            stamp = os.stat(path).st_mtime_ns
            said.clear()
            release.write(pkg, out, log=said.append)
            self.assertEqual((os.stat(path).st_mtime_ns, len(said)), (stamp, 1))

    def test_the_builds_step(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            w = World(tmp)
            collected = dict(genomes=GENOMES, dropped=[POOR])
            opts = types.SimpleNamespace(protal=None, protal_sketch_kmers=7, protal_sketch_mb=100, protal_sketch_exponent=1,
                                         threads=3, genomes=w.folder, work=w.work, report_minutes=10)
            stages = []
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(eukore_build.protal_step(opts, collected, stages), (None, None))
            self.assertIn("protal sketch: none made (no --protal)", out.getvalue())
            self.assertEqual((stages, w.ran()), ([], []))  # not timed, nothing run
            opts.protal = os.path.join(tmp, "no_such_protal")
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(eukore_build.protal_step(opts, collected, stages), (None, None))
            self.assertIn("no such program", out.getvalue())
            opts.protal = w.protal
            with contextlib.redirect_stdout(io.StringIO()) as out:
                writer, about = eukore_build.protal_step(opts, collected, stages)
            self.assertEqual((about["kmers"], about["species"], [s for s, _t in stages]),
                             (7, 2, ["protal genome sketch"] * 2))
            self.assertEqual(w.ran()[-1][-12:], ["--kmers", "7", "--kmers_mb", "100", "--kmers_exponent", "1.0",
                                                "--kmers_min", "7", "--kmers_max", "0", "-t", "3"])
            part = os.path.join(tmp, "part")
            writer(part)
            self.assertTrue(os.path.isfile(os.path.join(part, ps.IN_PACKAGE)))

    def test_the_options(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stderr(io.StringIO()) as err, \
                self.assertRaises(SystemExit):
            eukore_build.main(["--build_dir", tmp, "--out", tmp, "--protal_sketch_kmers", "0"])
        self.assertIn("--protal_sketch_kmers 0, --protal_sketch_mb 4: 1 or more each", err.getvalue())
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stderr(io.StringIO()) as err, \
                self.assertRaises(SystemExit):
            eukore_build.main(["--build_dir", tmp, "--out", tmp, "--protal_sketch_exponent", "0"])
        self.assertIn("--protal_sketch_exponent 0.0: above 0", err.getvalue())


if __name__ == "__main__":
    unittest.main()
