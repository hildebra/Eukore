#!/usr/bin/env python3
"""Tests of EUKDB's own locus calls (odb.py, loci.py; since 2026-11-05): the dataset read (cutoffs, the HMMs' own
names, proteins by species, the nearest species' queries first); miniprot's GFF with --trans parsed; overlapping
alignments one region; a region confirmed by its locus' HMM at the cutoff; a locus Single, Duplicated (two complete
regions), Fragmented (a short HMM alignment) or Missing; the copy whole on its own strand with its exons and introns, its
soft-masked and low-complexity runs; the windows around it.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukdb import loci, ncbi as ncbi_mod, odb as odb_mod  # noqa: E402
import fixtures as fx  # noqa: E402


class OdbTest(unittest.TestCase):
    def test_the_dataset(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = fx.odb_world(tmp)
            odb = odb_mod.Odb(root, fx.LOCI)
            self.assertEqual(odb.cutoff["2at2759"], 50.0)
            self.assertEqual(odb.locus_of_hmm["2"], "2at2759")  # BUSCO's profile names
            self.assertEqual(odb.length["1at2759"], 120)
            self.assertEqual(odb.species(), [4932, 6253, 9606])
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            lineages = odb_mod.species_lineages(odb, ncbi)
            # a Kluyveromyces genome: the yeast's proteins first, then the nematode's or human's, then the ancestral
            q = odb.queries(list(ncbi.ancestors(28985)), lineages, nearest=1)
            self.assertEqual(sorted(q), sorted([f"{l}~{l}_4932_0:000004" for l in fx.LOCI] +
                                               [f"{l}~{l}_0" for l in fx.LOCI]))
            q = odb.queries(list(ncbi.ancestors(9606)), lineages, nearest=1)
            self.assertIn("1at2759~1at2759_9606_0:000002", q)
            with self.assertRaises(ValueError):
                odb_mod.Odb(root, ["1at2759", "9at2759"])
            path = odb.hmm_db(os.path.join(tmp, "all.hmm"))
            with open(path) as fh:
                self.assertEqual(fh.read().count("NAME"), 3)
            versions = os.path.join(tmp, "v.tsv")
            with open(versions, "w") as fh:
                fh.write("eukaryota_odb12\t2025-01-15\taaa\tEukaryota\tlineages\n"
                         "eukaryota_odb12\t2026-05-22\tbbb\tEukaryota\tlineages\n"
                         "eukaryota_odb12.2\t2026-05-13\tccc\tEukaryota\tlineages\n")
            self.assertEqual(odb_mod.latest(versions), ("2026-05-22", "bbb"))


class LociTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.odb = odb_mod.Odb(fx.odb_world(self.dir), fx.LOCI)
        self.seq, self.gff, self.spans = fx.genome()
        self.path = os.path.join(self.dir, "g.fna.gz")
        with gzip.open(self.path, "wt") as fh:
            fh.write(f">c1 a test contig\n{self.seq}\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_parsing_and_regions(self):
        got = loci.parse_alignments(self.gff)
        self.assertEqual(len(got), 6)
        self.assertEqual((got[0]["locus"], got[0]["strand"], len(got[0]["exons"]), got[0]["protein"]),
                         ("1at2759", "+", 2, fx.PROTEINS["1at2759"]))
        # each alignment its own translation (miniprot writes it before the alignment)
        self.assertTrue(all(a["protein"] == fx.PROTEINS[a["locus"]] for a in got))
        lengths = {a["query"]: 120 for a in got}
        reg = loci.regions(got, lengths)
        self.assertEqual(len(reg["1at2759"]), 2)  # two places: two regions
        self.assertEqual(len(reg["2at2759"]), 1)
        self.assertTrue(all(r["support"] == 2 and r["score"] == 900 for rs in reg.values() for r in rs))
        self.assertTrue(all(len(r["top"]) == 2 for rs in reg.values() for r in rs))  # each HMM-checked
        # under the identity or coverage gate: left out
        low = [dict(a, identity=0.1) for a in got]
        self.assertEqual(loci.regions(low, lengths), {})
        self.assertEqual(loci.regions(got, {q: 1000 for q in lengths}), {})

    def test_a_genome_called(self):
        got = loci.call(self.path, {}, self.odb, fx.miniprot(self.dir, self.gff), fx.hmmsearch(self.dir), 1, self.dir)
        self.assertIsNone(got["failed"])
        one, two, three = (got["loci"][l] for l in fx.LOCI)
        self.assertEqual((one["status"], one["copies"], two["status"], three["status"]),
                         ("Duplicated", 2, "Single", "Missing"))
        # locus 1: the first region (both score alike), whole with its intron; its coding sequence the protein's
        a, b = self.spans["1at2759"]
        self.assertEqual(one["copy"], self.seq[a - 1:b].upper())
        self.assertEqual(one["spliced"], fx.cds_of(fx.PROTEINS["1at2759"]))
        kinds = [(k, s, e) for k, s, e in one["features"] if k in ("exon", "intron")]
        self.assertEqual(kinds, [("exon", 0, 180), ("intron", 180, 480), ("exon", 480, b - a + 1)])
        soft = [(s, e) for k, s, e in one["features"] if k == "soft_masked"]
        self.assertEqual(soft, [(180 + 62, 180 + 142)])
        low = [(s, e) for k, s, e in one["features"] if k == "low_complexity"]
        self.assertTrue(any(s <= 180 + 150 and e >= 180 + 200 for s, e in low))  # the poly-A run (by 64-base windows)
        # locus 2: on the reverse strand, the protein's coding sequence read forward
        a, b = self.spans["2at2759"]
        self.assertEqual(two["copy"], self.seq[a - 1:b].translate(str.maketrans("ACGT", "TGCA"))[::-1])
        self.assertTrue(two["copy"].startswith("ATG"))
        self.assertEqual((two["hmm_score"], two["cutoff"], two["strand"]), (300.0, 50.0, "-"))
        self.assertEqual(set(n.split("|")[0] for n in got["windows"]), {"1at2759", "2at2759"})

    def test_a_fragment_and_a_failure(self):
        short = self.gff.replace(fx.PROTEINS["2at2759"], fx.PROTEINS["2at2759"][:90])  # 3/4 of the profile
        mp = fx.miniprot(self.dir, short)
        got = loci.call(self.path, {}, self.odb, mp, fx.hmmsearch(self.dir), 1, self.dir)
        self.assertEqual((got["loci"]["2at2759"]["status"], round(got["loci"]["2at2759"]["hmm_coverage"], 2)),
                         ("Fragmented", 0.75))
        other = self.gff.replace(fx.PROTEINS["2at2759"], fx.PROTEINS["3at2759"])  # another locus' protein: not 2
        got = loci.call(self.path, {}, self.odb, fx.miniprot(self.dir, other), fx.hmmsearch(self.dir), 1, self.dir)
        self.assertEqual(got["loci"]["2at2759"]["status"], "Missing")
        failing = fx.tool(self.dir, "mp_fail", "#!/bin/sh\necho 'out of memory' >&2\nexit 1\n")
        got = loci.call(self.path, {}, self.odb, failing, fx.hmmsearch(self.dir), 1, self.dir)
        self.assertTrue(got["failed"].startswith("miniprot failed: out of memory"))

    def test_masks(self):
        self.assertEqual(loci.soft_runs("ACgtacgtacgtacgtacgtacgtAC", 20), [(2, 24)])
        self.assertEqual(loci.soft_runs("ACgtAC", 20), [])
        self.assertEqual(loci.on_copy([(0, 10)], 101, 200, "-"), [(90, 100)])
        self.assertEqual(loci.dust("ACGT" * 5), [])
        rng = __import__("random").Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(300)) + "A" * 200 + "".join(rng.choice("ACGT") for _ in range(300))
        runs = loci.dust(seq)
        self.assertTrue(runs and runs[0][0] <= 320 and runs[-1][1] >= 480 and all(e - s <= 400 for s, e in runs))


if __name__ == "__main__":
    unittest.main()
