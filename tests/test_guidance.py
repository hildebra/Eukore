#!/usr/bin/env python3
"""Tests of the identity that parts taxa per marker and rank (benchmarks/eukore_bench/guidance.py and the build's
resolution_hist.tsv; since 2026-11-12).

    python3 -m unittest discover -s tests -t tests
"""

import importlib.util
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))
from eukore import benchmarks  # noqa: E402
from eukore_bench import guidance  # noqa: E402
import fixtures as fx  # noqa: E402


class GuidanceTest(unittest.TestCase):
    def test_the_builds_histogram(self):
        # identities counted in steps of half a point, per shared rank (benchmarks.DEPTHS: 5 genus, 4 family)
        rows = benchmarks.resolution_hist({"m": [(99.26, 5), ("99.4", 5), (87.0, 4), (None, None), (100.0, 5)]})
        self.assertEqual([(r["shared_rank"], r["identity"], r["copies"]) for r in rows],
                         [("family", "87.0", 1), ("genus", "100.0", 1), ("genus", "99.0", 2), ("unmatched", "-", 1)])

    def test_the_identity_that_parts_two_sets(self):
        # the highest identity that keeps all of within and none of between
        self.assertEqual(guidance.best_threshold({99.5: 3, 100.0: 3}, {99.0: 6, 90.0: 5, None: 3}),
                         (99.5, 1.0, 1.0, 1.0))
        # overlapping sets: the best balance
        t, kept, apart, accuracy = guidance.best_threshold({99.0: 2, 97.0: 2}, {98.0: 2, 90.0: 2})
        self.assertEqual((t, kept, apart, accuracy), (99.0, 0.5, 1.0, 0.75))

    def test_thresholds_per_rank(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, package = os.path.join(tmp, "results"), os.path.join(tmp, "Eukore_x")
            os.makedirs(out)
            os.makedirs(os.path.join(package, "intraspecific"))
            hist = [("genus", "99.0", 6), ("genus", "98.5", 2), ("family", "90.0", 5), ("family", "95.0", 1),
                    ("order", "80.0", 6), ("none", "70.0", 2), ("unmatched", "-", 3)]
            benchmarks.write(out, "resolution_hist.tsv", benchmarks.HIST_COLUMNS,
                             [dict(marker=benchmarks.ALL_BUSCOS, shared_rank=r, identity=i, copies=n)
                              for r, i, n in hist] +
                             [dict(marker="18S", shared_rank="genus", identity="99.5", copies=2)])  # too few
            with open(os.path.join(package, "intraspecific", "barcode_gap.tsv"), "w") as fh:
                fh.write("marker\tspecies_taxid\tintra_min\n" + "".join(
                    f"1at2759\t{i}\t{v}\n" for i, v in enumerate(("99.5", "99.5", "99.5", "100.0", "100.0", "100.0")))
                    + "1at2759\t9\t-\n")
            summary = {"species": 6, "busco_loci": 1, "left_out": {}}
            guidance.run(package, out, summary, log=lambda t: None)
            rows = {(r["marker"], r["rank"]): r for r in fx.read_tsv(os.path.join(out, "thresholds.tsv"))}
            buscos = benchmarks.ALL_BUSCOS
            # a species: the panel's own assemblies (99.5 or more alike) against every copy's nearest other species
            self.assertEqual([rows[(buscos, "species")][k] for k in guidance.COLUMNS[2:]],
                             ["99.5", "6", "1.000", "25", "1.000", "1.000", "package"])
            # a genus: the copies whose nearest is a congener against those whose nearest is not
            self.assertEqual([rows[(buscos, "genus")][k] for k in ("threshold", "within", "between")],
                             ["98.5", "8", "17"])
            self.assertEqual((rows[(buscos, "family")]["threshold"], rows[(buscos, "order")]["threshold"]),
                             ("90.0", "80.0"))
            self.assertNotIn(("18S", "genus"), rows)  # two copies: no threshold said
            self.assertEqual(rows[("ITS (published)", "species")]["threshold"], "99.6")
            self.assertTrue(rows[("ITS (published)", "species")]["source"].startswith("published: "))
            self.assertEqual(summary["thresholds"][buscos]["genus"], {"threshold": 98.5, "accuracy": 1.0})
            self.assertIn("| identity that parts taxa, BUSCOs | species 99.5% (100.0% right); genus 98.5%",
                          benchmarks.summary_md(summary, "Eukore_x"))
            if importlib.util.find_spec("matplotlib") is not None:
                from eukore_bench import figures  # noqa: PLC0415
                self.assertTrue(os.path.isfile(figures.thresholds(out, package, out, "png")))
            # a package without the histogram: said, nothing written
            other = os.path.join(tmp, "none")
            os.makedirs(other)
            left = guidance.run(package, other, {"left_out": {}}, log=lambda t: None)["left_out"]
            self.assertIn("no resolution_hist.tsv", left["thresholds"])
            self.assertFalse(os.listdir(other))


if __name__ == "__main__":
    unittest.main()
