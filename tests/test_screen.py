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


CMSEARCH = '''#!/bin/sh
# cmsearch --tblout's stand-in: s0 scored best by the eukaryote model, s1 by the bacterial one by 40 bits, s2 by the
# bacterial one by 10 bits (under the margin: still eukaryote), s3 by none
while [ "$#" -gt 0 ]; do [ "$1" = "--tblout" ] && out=$2; shift; done
cat > "$out" <<EOT
#target name accession query name accession mdl mdl from mdl to seq from seq to strand trunc pass gc bias score ...
s0 - SSU_rRNA_eukarya RF01960 hmm 1 1800 1 1800 + - 6 0.5 0 1500.0 0 ! -
s0 - SSU_rRNA_bacteria RF00177 hmm 1 1500 1 1500 + - 6 0.5 0 900.0 0 ! -
s1 - SSU_rRNA_eukarya RF01960 hmm 1 1800 1 1500 + - 6 0.5 0 700.0 0 ! -
s1 - SSU_rRNA_bacteria RF00177 hmm 1 1500 1 1500 + - 6 0.5 0 740.0 0 ! -
s2 - SSU_rRNA_eukarya RF01960 hmm 1 1800 1 1500 + - 6 0.5 0 700.0 0 ! -
s2 - SSU_rRNA_bacteria RF00177 hmm 1 1500 1 1500 + - 6 0.5 0 710.0 0 ! -
EOT
'''


class SsuModelTest(unittest.TestCase):
    def test_an_18s_told_from_a_bacterial_ssu(self):
        # since 2026-11-09 (regions.ssu_models): a bacterium's or an organelle's SSU is no nucleus' 18S
        from eukore import regions, rrna  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            cmsearch = fx.tool(tmp, "cmsearch", CMSEARCH)
            seqs = {"euk": "ACGT" * 10, "bac": "CCGT" * 10, "near": "GCGT" * 10, "none": "TCGT" * 10,
                    "twin": "ACGT" * 10}
            got = regions.ssu_models(seqs, os.path.join(tmp, "x.cm"), cmsearch)
            self.assertEqual(got, {"euk": "eukarya", "twin": "eukarya", "bac": "bacteria", "near": "eukarya",
                                   "none": "-"})
        row = dict(kind="18S", own_species="t", screen="-", ambiguous="0", gaps=0, ssu_model="bacteria")
        self.assertEqual(rrna.confidence(row), ("low", "ssu_prokaryote"))
        self.assertEqual(rrna.confidence(dict(row, ssu_model="eukarya")), ("high", "-"))
        self.assertEqual(rrna.confidence(dict(row, ssu_model="-", own_species="-")),
                         ("unverified", "own_species_undecided"))

    def test_an_assembly_one_genome_whatever_its_copy(self):
        # since 2026-11-09: ROD's operons of the assembly a record comes from (another version, RefSeq's copy) are
        # that genome's, no confirmation of it
        self.assertEqual(screen.genome_of("a", {"accession": "GCF_000146045.2"}),
                         screen.genome_of("b", {"accession": "GCA_000146045.1"}))
        self.assertNotEqual(screen.genome_of("a", {"accession": "GCA_000146045.2"}),
                            screen.genome_of("b", {"accession": "GCA_000146046.2"}))
        self.assertEqual(screen.genome_of("EUK_1", {"accession": "-"}), "EUK_1")


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
            # the mosaic genome holding 40 copies of its unit (a long-read assembly's rDNA array): its own copies
            # fill vsearch's 32 first hits unless searched once (screen.representatives); every copy flagged
            unit = records[-1][1]
            array = [(f"mosaic{i}", unit, dict(records[-1][2])) for i in range(40)]  # an array's identical copies
            flags, _counts = screen.screen(records[:-1] + array, ncbi, fx.vsearch(tmp))
            self.assertEqual(set(flags), {name for name, _s, _r in array})
            # the own-species check alike: the array's 18S (S. cerevisiae's) agree with its lineage, not undecided
            from eukore import rrna  # noqa: PLC0415
            rows = [(n, s, dict(r, source="genome", kind="rDNA")) for n, s, r in records[:-1] + array]
            own = rrna.own_species(rows, ncbi, fx.vsearch(tmp))
            self.assertEqual({own[n][0] for n, _s, _r in array}, {"t"})  # (verdict, cross) since 2026-11-09

    def test_a_protist_without_a_phylum_is_tested_deeper(self):
        # NCBI names its class but no phylum: its clade at class is (kingdom, class), not None
        lineage = ["Eukaryota", None, "Someclass", "Someorder", None]
        self.assertEqual(screen.clade(lineage, 2), ("Eukaryota", "Someclass"))
        self.assertIsNone(screen.clade(lineage, 1))
        self.assertEqual(screen.agrees_to(lineage, ["Eukaryota", None, "Someclass", "Other", None]), 2)


if __name__ == "__main__":
    unittest.main()
