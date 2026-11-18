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

    def test_a_barcode_cut_at_its_primers_holds_no_site(self):
        # since 2026-11-17: a published barcode is the amplicon without its primers (MIDORI2's COX1), so the pair
        # amplifies no such record; the species whose record holds both sites are counted apart
        found, pairs = primers.load()
        folmer = next(p for p in pairs if p["pair"] == "COI_Folmer")
        limits = primers.site_limits(found, pairs)
        self.assertEqual(limits["COI_Folmer"], (len(found["LCO1490"]["sequence"]) // 3,
                                                len(found["HCO2198"]["sequence"]) // 3))
        whole = references()[folmer["reference"]]
        got = primers.amplify(whole, folmer, found)
        self.assertEqual(got["amplifies"], "t")
        self.assertTrue(primers.holds_sites(got, limits))
        inside = whole[got["start"] + len(found["LCO1490"]["sequence"]):got["end"] - len(found["HCO2198"]["sequence"])]
        cut = primers.amplify(inside, folmer, found)  # the barcode alone: each primer's best place a chance one
        self.assertEqual(cut["amplifies"], "f")
        self.assertFalse(primers.holds_sites(cut, limits))
        self.assertFalse(primers.holds_sites(dict(pair="COI_Folmer", forward_mismatches=0), limits))  # no reverse site
        self.assertFalse(primers.holds_sites(dict(pair="unknown", forward_mismatches=0, reverse_mismatches=0), limits))

    def test_the_benchmark_counts_the_species_with_both_sites(self):
        import tempfile  # noqa: PLC0415
        from eukore import benchmarks  # noqa: PLC0415

        def tsv(path, head, rows):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w", newline="\n") as fh:
                fh.write("\t".join(head) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
        with tempfile.TemporaryDirectory() as tmp:
            tsv(os.path.join(tmp, "organelle", "genes.tsv"), ("record", "species_taxid"),
                [("whole", 1), ("off", 2), ("barcode", 3), ("barcode2", 4)])
            head = ("record", "pair", "forward_mismatches", "reverse_mismatches", "amplifies")
            tsv(os.path.join(tmp, "organelle", "amplicons.tsv"), head,
                [("whole", "COI_Folmer", 1, 2, "t"),      # amplified
                 ("off", "COI_Folmer", 5, 4, "f"),        # both sites there, too many mismatches
                 ("barcode", "COI_Folmer", 11, 12, "f"),  # the amplicon alone: no site
                 ("barcode2", "COI_Folmer", 10, "-", "f")])
            species = {str(i): dict(supergroup="Metazoa") for i in (1, 2, 3, 4)}
            rows = benchmarks.primer_rows(tmp, species)
            row = next(r for r in rows if r["supergroup"] == "all")
            self.assertEqual((row["species"], row["amplified"], row["share"], row["with_sites"], row["share_of_sites"]),
                             (4, 1, "0.250", 2, "0.500"))
            target = benchmarks.primers_on_target(rows)["COI_Folmer"]
            self.assertEqual((target["share"], target["with_sites"], target["share_of_sites"]), (0.25, 2, 0.5))
            # a table of before (no with_sites): the share of all species alone
            old = [{k: r[k] for k in ("pair", "supergroup", "species", "amplified", "share")} for r in rows]
            self.assertNotIn("share_of_sites", benchmarks.primers_on_target(old)["COI_Folmer"])
            self.assertIn("COI_Folmer 50.0% of 2", benchmarks.summary_md(
                dict(species=4, busco_loci=1, primers_target=benchmarks.primers_on_target(rows)), "p"))


if __name__ == "__main__":
    unittest.main()
