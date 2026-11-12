#!/usr/bin/env python3
"""Tests of a release among releases (benchmarks/eukore_bench/history.py; since 2026-11-12): the package's parts'
sizes, its genomes by the year NCBI released them, the changes from an earlier package.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import importlib.util
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore import benchmarks  # noqa: E402
from eukore_bench import history  # noqa: E402
import fixtures as fx  # noqa: E402


def package(folder, species, genomes, files=None, version="x"):
    """A package's species.tsv ([(taxid, name, group, representative)]), genomes.tsv ([(accession, fetched, species
    taxid)]) and manifest."""
    os.makedirs(os.path.join(folder, "metadata"), exist_ok=True)
    with open(os.path.join(folder, "species.tsv"), "w") as fh:
        fh.write("taxid\tname\tgroup\trepresentative\n" + "".join("\t".join(r) + "\n" for r in species))
    with open(os.path.join(folder, "genomes.tsv"), "w") as fh:
        fh.write("accession\tfetched\tspecies_taxid\trepresentative\n"
                 + "".join(f"{a}\t{f}\t{t}\tt\n" for a, f, t in genomes))
    with open(os.path.join(folder, "manifest.json"), "w") as fh:
        json.dump({"version": version, "built": "2026-11-12T10:00:00", "files": files or {}}, fh)
    return folder


class HistoryTest(unittest.TestCase):
    def test_the_parts_sizes(self):
        rows = history.size_rows({"files": {"species.tsv": [100, "m"], "busco/copies/1.fna.gz": [5000, "m"],
                                            "busco/copies/2.fna.gz": [3000, "m"], "busco/cds/1.fna.gz": [900, "m"],
                                            "busco/copies.tsv": [700, "m"], "rrna/records.fna.gz": [2000, "m"],
                                            "metadata/sketches/a.msh": [50, "m"]}})
        self.assertEqual([(r["part"], r["files"], r["bytes"]) for r in rows],
                         [("busco/copies/", 2, 8000), ("rrna/", 1, 2000), ("busco/cds/", 1, 900),
                          ("busco/ tables", 1, 700), ("species.tsv", 1, 100), ("metadata/", 1, 50),
                          ("all", 7, 11750)])
        self.assertEqual(history.size_rows({}), [])

    def test_one_assembly_two_versions(self):
        self.assertEqual(history.assembly_of("GCA_000146045.2"), history.assembly_of("GCF_000146045.3"))
        self.assertNotEqual(history.assembly_of("GCA_000146045.2"), history.assembly_of("GCA_000146046.2"))
        self.assertEqual(history.assembly_of("TOGA2_A1"), "TOGA2_A1")

    def test_changes_growth_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            old = package(os.path.join(tmp, "old"), [
                ("4932", "Saccharomyces cerevisiae", "Fungi", "GCA_000146045.2"),
                ("27291", "Saccharomyces paradoxus", "Fungi", "GCA_000000013.1"),
                ("28985", "Kluyveromyces lactis", "Fungi", "-"),
                ("5476", "Candida albicans", "Fungi", "GCA_000182965.3"),
                ("6253", "Ascaris lumbricoides", "Metazoa", "GCA_000000001.1")], [], version="2026.05")
            new = package(os.path.join(tmp, "new"), [
                ("4932", "Saccharomyces cerevisiae", "Fungi", "GCA_000146045.3"),     # a newer version
                ("27291", "Saccharomyces paradoxus", "Fungi", "GCA_999999999.1"),     # another assembly
                ("28985", "Kluyveromyces lactis", "Fungi", "GCA_000002515.1"),        # a genome at last
                ("5476", "Candida albicans", "Fungi", "GCA_000182965.3"),             # as it was
                ("9606", "Homo sapiens", "Vertebrata", "GCA_000001405.29")],          # new; Ascaris gone
                [("GCA_000146045.3", "GCF_000146045.3", "4932"), ("GCA_999999999.1", "GCA_999999999.1", "27291"),
                 ("GCA_000002515.1", "GCF_000002515.2", "28985"), ("GCA_000182965.3", "GCA_000182965.3", "5476"),
                 ("GCA_000001405.29", "GCA_000001405.29", "9606")],
                files={"species.tsv": [10, "m"], "rrna/records.tsv": [990, "m"]})
            with gzip.open(os.path.join(new, "metadata", "ncbi.tsv.gz"), "wt") as fh:
                fh.write("x\tassembly_accession\tseq_rel_date\n-\tGCA_000146045.3\t2011/04/18\n"
                         "-\tGCA_999999999.1\t2026/01/02\n-\tGCF_000002515.2\t2004/07/02\n"    # by the one fetched
                         "-\tGCA_000001405.29\t2026/03/01\n-\tGCA_000182965.3\tna\n")           # no date: not counted
            species = history.species_of(new)
            rows, listed = history.change_rows(history.species_of(old), species)
            at = {r["supergroup"]: r for r in rows}
            self.assertEqual(at["all"], dict(supergroup="all", previous=5, species=5, added=1, removed=1, kept=4,
                                             replaced=1, updated=1, gained_genome=1, lost_genome=0))
            self.assertEqual((at["Fungi"]["kept"], at["Metazoa"]["removed"], at["Vertebrata"]["added"]), (4, 1, 1))
            self.assertEqual({r["taxid"]: r["change"] for r in listed},
                             {"4932": "updated", "27291": "replaced", "28985": "gained_genome", "6253": "removed",
                              "9606": "added"})
            self.assertEqual(next(r for r in listed if r["taxid"] == "27291")["previous_genome"], "GCA_000000013.1")
            out = os.path.join(tmp, "results")
            os.makedirs(out)
            benchmarks.write(out, "stages.tsv", ("step", "seconds"), [dict(step="a", seconds="60.0"),
                                                                     dict(step="b", seconds="30.5")])
            summary = {"species": 5, "busco_loci": 3, "left_out": {}, "runtime": {"seconds": 7200.0, "threads": "2", "median_seconds": 60.0}}
            history.run(new, out, species, summary, old, log=lambda t: None)
            self.assertEqual(summary["sizes"], {"bytes": 1000, "files": 2,
                                                "largest": {"rrna/": 990, "species.tsv": 10}})
            self.assertEqual(summary["cost"], {"package_stage_seconds": 90.5, "genome_runs_seconds": 7200.0,
                                               "genome_runs_threads": "2"})
            growth = {(r["supergroup"], r["year"]): r["genomes"] for r in fx.read_tsv(os.path.join(out, "growth.tsv"))}
            self.assertEqual(growth, {("Fungi", "2004"): "1", ("Fungi", "2011"): "1", ("Fungi", "2026"): "1",
                                      ("Vertebrata", "2026"): "1", ("all", "2004"): "1", ("all", "2011"): "1",
                                      ("all", "2026"): "2"})
            self.assertEqual(summary["growth"], {"genomes": 4, "first_year": 2004, "last_year": 2026,
                                                 "released_in_last_two_years": 0.5})
            self.assertEqual((summary["changes"]["added"], summary["changes"]["previous_version"]), (1, "2026.05"))
            md = benchmarks.summary_md(summary, "Eukore_x")
            self.assertIn("| the build's cost | genome runs 2.0 h (on 2 threads each); package stage 2 min |", md)
            self.assertIn("| genomes NCBI released in 2025 and 2026 | 50.0% of 4 |", md)
            self.assertIn("| since 2026.05: species added; removed; genomes replaced; updated | 1; 1; 1; 1 "
                          "(5 to 5 species) |", md)
            with self.assertRaises(OSError):
                history.run(new, out, species, summary, os.path.join(tmp, "nothing"), log=lambda t: None)
            if importlib.util.find_spec("matplotlib") is None:
                return
            from eukore_bench import figures  # noqa: PLC0415
            for fn in (figures.sizes, figures.growth, figures.changes):
                self.assertTrue(os.path.isfile(fn(out, new, out, "png")), fn.__name__)

    def test_other_databases_release_dates(self):
        dates = history.release_dates()
        self.assertEqual(dates["UNITE"], (2025, 2))
        self.assertNotIn("Eukore", dates)


if __name__ == "__main__":
    unittest.main()
