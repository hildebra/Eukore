#!/usr/bin/env python3
"""Tests of the rRNA part screen (eukore/screen.py): a record whose 28S is another family's (a mosaic) is flagged at
that part and rank; records whose parts agree with their lineages, and a matched record whose clade has no other
genome, are not; parts cut from regions.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi as ncbi_mod, screen  # noqa: E402
import fixtures as fx  # noqa: E402


def mutate(seq, rate, rng):
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in seq)


class ScreenTest(unittest.TestCase):
    def test_parts(self):
        seq = "A" * 4000
        got = screen.parts_of(seq, [("18S", 0, 1800), ("ITS1", 1800, 2000), ("5.8S", 2000, 2150),
                                    ("ITS2", 2150, 2400), ("28S", 2400, 4000)])
        self.assertEqual({k: len(v) for k, v in got.items()},
                         {"18S": 1800, "junction": 700, "28S 5'": 800, "28S 3'": 800})
        self.assertEqual(list(screen.parts_of("A" * 1000, None)), ["18S"])
        self.assertEqual(screen.parts_of("A" * 100, None), {})

    def test_a_mosaic_flagged(self):
        rng = random.Random(7)
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            ssu = {"sacch": "".join(rng.choice("ACGT") for _ in range(1000)),
                   "kluy": None, "cand": "".join(rng.choice("ACGT") for _ in range(1000))}
            ssu["kluy"] = mutate(ssu["sacch"], 0.03, rng)  # same family: 97%
            lsu = {"sacch": "".join(rng.choice("ACGT") for _ in range(1600)),
                   "cand": "".join(rng.choice("ACGT") for _ in range(1600))}
            lsu["kluy"] = mutate(lsu["sacch"], 0.03, rng)
            regions = [("18S", 0, 1000), ("28S", 1000, 2600)]
            records = []
            taxa = {"sacch": [4932, 27291], "kluy": [28985], "cand": [5476]}
            for group, taxids in taxa.items():
                for i, taxid in enumerate(taxids * 2):
                    s = mutate(ssu[group], 0.002, rng) + mutate(lsu[group], 0.002, rng)
                    records.append((f"{group}{i}", s, dict(accession=f"G_{group}{i}", taxid=taxid,
                                                           species_taxid=taxid, regions_list=regions)))
            # a S. cerevisiae genome whose 28S is Candida's (another family's), its 18S its own
            records.append(("mosaic", mutate(ssu["sacch"], 0.002, rng) + mutate(lsu["cand"], 0.002, rng),
                            dict(accession="G_mosaic", taxid=4932, species_taxid=4932, regions_list=regions)))
            flags, counts = screen.screen(records, ncbi, fx.vsearch(tmp))
            self.assertEqual(flags.get("mosaic"), "28S 5':family,28S 3':family")
            self.assertEqual(set(flags), {"mosaic"})
            self.assertTrue(counts)


if __name__ == "__main__":
    unittest.main()
