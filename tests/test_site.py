#!/usr/bin/env python3
"""Tests of the benchmarks' results for the website (benchmarks/eukore_bench/site.py; since 2026-11-12): the
release's data/benchmarks.json and its figures, the build's own files there untouched.

    python3 -m unittest discover -s tests -t tests
"""

import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore import benchmarks  # noqa: E402

try:
    from eukore_bench import site  # noqa: E402  (its captions are the figures': matplotlib)
except ImportError:
    site = None


@unittest.skipUnless(site, "no matplotlib (the benchmarks' environment)")
class SiteTest(unittest.TestCase):
    def test_the_results_for_the_website(self):
        with tempfile.TemporaryDirectory() as tmp:
            results, release = os.path.join(tmp, "results"), os.path.join(tmp, "release")
            os.makedirs(os.path.join(results, "figures"))
            os.makedirs(os.path.join(release, "data"))
            with open(os.path.join(release, "data", "stats.json"), "w") as fh:
                json.dump({"package": "Eukore_2026.11"}, fh)
            benchmarks.write(results, "breadth_ranks.tsv", ("supergroup", "rank", "eukore"),
                             [dict(supergroup="Fungi", rank="genus", eukore=3),
                              dict(supergroup="all", rank="genus", eukore=5)])
            benchmarks.write(results, "databases.tsv", ("database", "species"),
                             [dict(database="Eukore", species=9), dict(database="PR2", species=7)])
            benchmarks.write(results, "thresholds.tsv", ("marker", "rank", "threshold", "source"),
                             [dict(marker="BUSCO (all)", rank="genus", threshold="98.5", source="package"),
                              dict(marker="1at2759", rank="genus", threshold="97.0", source="package"),
                              dict(marker="ITS (published)", rank="genus", threshold="94.3", source="published: x")])
            with open(os.path.join(results, "summary.json"), "w") as fh:
                json.dump({"species": 9, "tree": {"genomes": 4}, "upset": {"eukore_alone": 2},
                           "completeness": {"median": 0.9}, "left_out": {"reads": "no minimap2"}}, fh)
            for name in ("tree", "selection"):  # a figure the page shows, and one it does not
                with open(os.path.join(results, "figures", f"{name}.png"), "wb") as fh:
                    fh.write(b"PNG")
            path = site.write(results, release, log=lambda t: None)
            with open(path) as fh:
                got = json.load(fh)
            self.assertEqual(got["package"], "Eukore_2026.11")  # the release's own
            self.assertEqual(got["tables"]["breadth_ranks"], [dict(supergroup="all", rank="genus", eukore="5")])
            self.assertEqual(len(got["tables"]["databases"]), 2)  # a table without groups: whole
            self.assertEqual([r["marker"] for r in got["tables"]["thresholds"]], ["BUSCO (all)", "ITS (published)"])
            self.assertNotIn("panels", got["tables"])  # no such table: left out
            self.assertEqual(got["summary"], {"upset": {"eukore_alone": 2}, "tree": {"genomes": 4}})
            self.assertEqual(got["left_out"], {"reads": "no minimap2"})
            self.assertEqual([(f["name"], f["file"]) for f in got["figures"]], [("tree", "figures/tree.png")])
            self.assertTrue(got["figures"][0]["caption"])
            self.assertEqual(sorted(os.listdir(os.path.join(release, "data"))),
                             ["benchmarks.json", "figures", "stats.json"])
            self.assertEqual(os.listdir(os.path.join(release, "data", "figures")), ["tree.png"])
            self.assertEqual(site.write(results, release, package="/x/Eukore_y", log=lambda t: None), path)
            with open(path) as fh:
                self.assertEqual(json.load(fh)["package"], "Eukore_y")
            with self.assertRaises(SystemExit):  # no release there
                site.write(results, os.path.join(tmp, "nothing"), log=lambda t: None)


if __name__ == "__main__":
    unittest.main()
