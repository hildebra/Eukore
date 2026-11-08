#!/usr/bin/env python3
"""Tests of the in-silico PCR (eukore/primers.py; since 2026-11-07): every pair of data/primer_pairs.tsv amplifies its
reference record (data/primer_references.fna), TAReuk's V4 where Stoeck et al. 2010 place it on yeast's 18S, the
mismatch and 3'-end rules, and which records a pair runs on.

    python3 -m unittest discover -s tests -t tests
"""

import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import primers  # noqa: E402

REFERENCES = os.path.join(os.path.dirname(HERE), "data", "primer_references.fna")


def references():
    out, name = {}, None
    with open(REFERENCES) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                out[name] = line.strip()
    return out


class PrimerTest(unittest.TestCase):
    def test_every_pair_amplifies_its_reference(self):
        found, pairs = primers.load()
        refs = references()
        for pair in pairs:
            for name in (pair["forward"], pair["reverse"]):
                self.assertIn(name, found, pair["pair"])
            got = primers.amplify(refs[pair["reference"]], pair, found)
            self.assertEqual(got["amplifies"], "t", (pair["pair"], got))

    def test_v4_on_yeast(self):
        found, pairs = primers.load()
        v4 = next(p for p in pairs if p["pair"] == "V4_TAReuk")
        got = primers.amplify(references()["Sc_18S"], v4, found)
        # Stoeck et al. 2010: S. cerevisiae 565-584 and 964-981 (1-based)
        self.assertEqual((got["start"], got["end"], got["forward_mismatches"], got["reverse_mismatches"]),
                         (564, 981, 0, 0))

    def test_rules(self):
        pair = dict(pair="p", forward="F", reverse="R", min_length="20", max_length="60")
        found = {"F": dict(sequence="GATCCTAGGA"), "R": dict(sequence=primers.revcomp("GGCCTTAAGG"))}
        record = "T" * 5 + "GATCCTAGGA" + "C" * 20 + "GGCCTTAAGG" + "T" * 5
        got = primers.amplify(record, pair, found)
        self.assertEqual((got["start"], got["length"], got["amplifies"]), (5, 40, "t"))
        three = record.replace("GATCCTAGGA", "GATCCTAGGT")  # the forward primer's 3'-most base wrong
        self.assertEqual(primers.amplify(three, pair, found)["amplifies"], "f")
        four = record.replace("GATCCTAGGA", "CTAGCTAGGA")  # four mismatches away from its 3' end
        got = primers.amplify(four, pair, found)
        self.assertEqual((got["forward_mismatches"], got["amplifies"]), (4, "f"))
        self.assertEqual(primers.amplify(record, dict(pair, max_length="30"), found)["amplifies"], "f")  # too long
        self.assertEqual(primers.amplify("ACGT", pair, found), None)  # shorter than a primer
        n = record.replace("GATCCTAGGA", "GATCNTAGGA")  # a record's N matches nothing
        self.assertEqual(primers.amplify(n, pair, found)["forward_mismatches"], 1)

    def test_records_a_pair_runs_on(self):
        _found, pairs = primers.load()
        by = {p["pair"]: p for p in pairs}
        self.assertTrue(primers.runs_on(by["V4_TAReuk"], {}, "rrna"))
        self.assertFalse(primers.runs_on(by["COI_Folmer"], {}, "rrna"))
        self.assertTrue(primers.runs_on(by["COI_Folmer"], dict(gene="COX1", organelle="mitochondrion"), "organelle"))
        self.assertTrue(primers.runs_on(by["MiFish-U"], dict(gene="SSU", organelle="mitochondrion"), "organelle"))
        self.assertFalse(primers.runs_on(by["MiFish-U"], dict(gene="SSU", organelle="plastid"), "organelle"))

    def test_records_in_a_pool_as_on_one_core(self):
        # run's pool (since 2026-11-08): the rows of one core's run, in the records' order
        from unittest import mock  # noqa: PLC0415
        yeast = references()["Sc_18S"]
        records = [(f"r{i}", yeast[i % 50:], dict(kind="18S")) for i in range(7)]
        alone = primers.run(records, "rrna")
        with mock.patch.object(primers, "CHUNK", 2):
            self.assertEqual(primers.run(records, "rrna", threads=3), alone)
        self.assertEqual(len({r["record"] for r in alone}), 7)


if __name__ == "__main__":
    unittest.main()
