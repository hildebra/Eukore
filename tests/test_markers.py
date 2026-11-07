#!/usr/bin/env python3
"""Tests of Eukore's own marker loci (eukore/markers.py; since 2026-11-07): the table, KOfam's rule of a candidate's
KO, the members kept, the dataset read beside BUSCO's (its decoys among the HMMs, its species apart), a region whose
best profile is a decoy's not called, and the markers' own folder in the package.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import sys
import tempfile
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import loci, markers, odb as odb_mod  # noqa: E402
import fixtures as fx  # noqa: E402


def marker_dataset(folder):
    """A marker dataset of one locus, MK1 (two target profiles), and one decoy."""
    for sub in ("hmms", "decoys"):
        os.makedirs(os.path.join(folder, sub))
    with open(os.path.join(folder, "hmms", "MK1.hmm"), "w") as fh:
        fh.write(fx.hmm_text("MK1", [{"L": 0.9}] * 3) + fx.hmm_text("MK1.2", [{"S": 0.9}] * 3))
    with open(os.path.join(folder, "decoys", "K99999.hmm"), "w") as fh:
        fh.write(fx.hmm_text("decoy-K99999", [{"A": 0.9}] * 3))
    with open(os.path.join(folder, "scores_cutoff"), "w") as fh:
        fh.write("MK1\t100.0\n")
    with open(os.path.join(folder, "links_to_ODB12.txt"), "w") as fh:
        fh.write("MK1\ta marker\n")
    with open(os.path.join(folder, "hmm_spans"), "w") as fh:
        fh.write("MK1\t0.600\nMK1.2\t0.900\n")
    with gzip.open(os.path.join(folder, "refseq_db.faa.gz"), "wt") as fh:
        fh.write(">MK1_28985_P12345\nMKLLL\n>MK1_6253_Q1\nMKLSL\n")
    return folder


class MarkerTest(unittest.TestCase):
    def test_table(self):
        got = {m["locus"]: m for m in markers.read_table()}
        self.assertTrue({"RPB1", "RPB2", "EF2", "TEF1"} <= set(got))
        self.assertEqual((got["RPB1"]["targets"], got["RPB1"]["kind"]), (["K03006"], "clean"))
        self.assertIn("K02999", got["RPB1"]["decoys"])  # RNA polymerase I's largest subunit
        self.assertTrue(all("_" not in m for m in got))  # refseq_db's names: <locus>_<taxid>_<id>

    def test_kofam_rule_and_members(self):
        kos = {"K1": (100.0, "full", "target"), "K2": (100.0, "domain", "paralog")}
        self.assertEqual(markers.classify({"K1": (300.0, 290.0), "K2": (200.0, 190.0)}, kos, {"K1", "K2"}), "K1")
        self.assertIsNone(markers.classify({"K1": (90.0, 90.0)}, kos, {"K1", "K2"}))  # under its threshold
        self.assertEqual(markers.classify({"K1": (150.0, 140.0), "K2": (400.0, 120.0)}, kos, {"K1", "K2"}), "K2")
        actin = {"KA": (903.0, 0.0, "beta/gamma"), "KB": (831.0, 0.0, "other eukaryote"), "KD": (400.0, 0.0, "ARP2")}
        actin = {k: (t, "full" if k == "KA" else "domain", d) for k, (t, _s, d) in actin.items()}
        plant = {"KA": (880.0, 870.0), "KB": (858.0, 850.0), "KD": (500.0, 490.0)}
        self.assertEqual(markers.classify(plant, actin, {"KA", "KB"}), "KA")  # under KA's threshold, over KB's
        self.assertIsNone(markers.classify(dict(plant, KD=(990.0, 980.0)), actin, {"KA", "KB"}))  # a decoy's best
        gapdh = [dict(accession="P25858", organism_id="3702", sequence="M" * 338),
                 dict(accession="Q9SAJ6", organism_id="3702", sequence="M" * 422),
                 dict(accession="P04406", organism_id="9606", sequence="M" * 335)]
        kept = markers.members(gapdh, {"P25858", "Q9SAJ6", "P04406"}, {"P25858": 712.3, "Q9SAJ6": 714.0,
                                                                        "P04406": 787.0})
        self.assertEqual([r["accession"] for r in kept], ["P25858", "P04406"])  # cytosolic GAPC, not GAPCp
        rows = [dict(accession="A", organism_id="1", sequence="M" * 100), dict(accession="B", organism_id="1",
                                                                              sequence="M" * 100),
                dict(accession="C", organism_id="2", sequence="M" * 10), dict(accession="D", organism_id="3",
                                                                             sequence="M" * 105)]
        kept = markers.members(rows, {"A", "B", "C", "D"}, {"A": 10.0, "B": 20.0, "C": 30.0, "D": 5.0})
        self.assertEqual([r["accession"] for r in kept], ["B", "D"])  # one per species, a fragment left out
        self.assertEqual(markers.renamed("HMMER3/f\nNAME  K03006\nLENG  3\n", "RPB1"), "HMMER3/f\nNAME  RPB1\nLENG  3\n")

    def test_read_beside_buscos(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = fx.odb_world(tmp)
            mk = marker_dataset(os.path.join(tmp, "eukore_markers"))
            odb = odb_mod.Odb(root, fx.LOCI + ["MK1"], [mk])
            self.assertEqual((odb.markers, odb.cutoff["MK1"], odb.names["MK1"]), ({"MK1"}, 100.0, "a marker"))
            self.assertEqual((odb.locus_of_hmm["MK1"], odb.locus_of_hmm["MK1.2"]), ("MK1", "MK1"))
            self.assertEqual(odb.span, {"MK1": 0.6, "MK1.2": 0.9})
            self.assertNotIn(28985, odb.species())  # the markers' UniProt species: no held-out calls, no benchmark
            self.assertIn(28985, odb.all_species())  # but their lineages rank the queries
            db = odb.hmm_db(os.path.join(tmp, "all.hmm"))
            with open(db) as fh:
                text = fh.read()
            self.assertIn("NAME  decoy-K99999", text)
            self.assertEqual(text.count("NAME "), len(fx.LOCI) + 3)
            self.assertIn("MK1~MK1_6253_Q1", odb.queries([1, 2759, 33208, 6253], odb_mod.species_lineages(
                odb, __import__("eukore.ncbi", fromlist=["Ncbi"]).Ncbi(fx.taxdump(os.path.join(tmp, "t.tar.gz"))))))
            self.assertNotEqual(odb.digest(), odb_mod.Odb(root, fx.LOCI).digest())  # a new hmm_db's name

    def test_ncbi_choice_and_spans(self):
        summaries = [("XP_1.1", "9606", 860), ("XP_2.1", "9606", 858), ("XP_3.1", "9606", 400), ("XP_4.1", "9606", 861),
                     ("XP_5.1", "9606", 2000), ("EAA1.1", "5833", 850), ("AAB1.1", "10090", 858)]
        self.assertEqual(markers.chosen(summaries), ["AAB1.1", "EAA1.1", "XP_1.1", "XP_2.1", "XP_4.1"])
        rows = [dict(accession="A", organism_id="559292", sequence="M" * 100),
                dict(accession="B", organism_id="4932", sequence="M" * 100)]
        kept = markers.members(rows, {"A", "B"}, {"A": 10.0, "B": 20.0}, lambda t: "4932")
        self.assertEqual([r["accession"] for r in kept], ["B"])  # a strain's and its species' proteins: one species
        odb = types.SimpleNamespace(locus_of_hmm={"TEF1": "TEF1", "TEF1.2": "TEF1"}, span={"TEF1": 0.6, "TEF1.2": 0.9})
        got = odb_mod.Odb.by_locus(odb, {"TEF1": (400.0, 0.55), "TEF1.2": (300.0, 0.9), "decoy-K1": (200.0, 0.7)})
        self.assertEqual(got["TEF1"][0], 400.0)
        self.assertAlmostEqual(got["TEF1"][1], 0.55 / 0.6)  # over its span share: complete as a BUSCO's would be
        self.assertEqual(got["decoy-K1"], (200.0, 0.7))
        self.assertEqual(odb_mod.Odb.by_locus(odb, {"TEF1": (400.0, 0.7)})["TEF1"][1], 1.0)  # at most 1

    def test_member_checks(self):
        self.assertEqual([markers.source_rank(a) for a in ("P04050", "NP_000928.1", "XP_1.1", "EAA1.1")], [0, 1, 2, 3])
        rows = [dict(accession="XP_9.1", organism_id="9606", sequence="M" * 100),
                dict(accession="NP_1.1", organism_id="9606", sequence="M" * 100),
                dict(accession="EAA1.1", organism_id="5833", sequence="M" * 100),
                dict(accession="EAB1.1", organism_id="5811", sequence="M" * 100)]
        own = {"XP_9.1": 900.0, "NP_1.1": 800.0, "EAA1.1": 700.0, "EAB1.1": 700.0}
        kept = markers.members(rows, set(own), own, complete={"XP_9.1", "NP_1.1", "EAA1.1"})
        self.assertEqual([r["accession"] for r in kept], ["EAA1.1", "NP_1.1"])  # curated first; a partial model out

        class Depth:  # taxon t's phylum: t // 10; the same phylum: depth 3
            def __call__(self, a, b):
                return 3 if a // 10 == b // 10 else 0

            def phylum(self, t):
                return t // 10

            def kingdom(self, t):
                return -1
        import random  # noqa: PLC0415
        rnd = random.Random(1)
        aa = "ACDEFGHIKLMNPQRSTVWY"

        def seq():
            return "".join(rnd.choice(aa) for _ in range(300))

        def mutated(s, n):
            s = list(s)
            for i in rnd.sample(range(len(s)), n):
                s[i] = rnd.choice(aa)
            return "".join(s)
        b = seq()  # phylum 2's protein (taxa 20, 21) also in a taxon of phylum 3 (31): its source's contaminant; 20
        proteins = [b, mutated(b, 6), mutated(b, 1)] + [seq() for _ in range(30)]  # is nearer 31's than 21's, kept
        taxa = [20, 21, 31] + [100 + 10 * i for i in range(30)]  # 30 more phyla, each its own protein
        self.assertEqual(markers.foreign(proteins, taxa, Depth()), {2})  # 20 and 21 have each other
        self.assertEqual(markers.foreign(proteins, taxa, Depth(), trusted=[2]), set())  # a curated record: kept
        h3 = seq()
        same = [mutated(h3, 2) for _ in range(30)]  # H3 of 30 lone phyla of one kingdom (animals'): no evidence

        class OneKingdom(Depth):
            def kingdom(self, t):
                return 1
        self.assertEqual(markers.foreign(same, [10 * i for i in range(30)], OneKingdom()), set())

    def test_scores_coverage(self):
        with tempfile.TemporaryDirectory() as tmp:
            binary, cache = fx.hmmsearch(tmp), os.path.join(tmp, "s.tsv.gz")  # a protein's start: its share covered
            p = fx.PROTEINS["1at2759"]
            got = markers.score({"whole": p, "part": p[:len(p) // 2]}, "db", binary, 1, cache)
            self.assertEqual(got["whole"]["1"], (300.0, 300.0, 1.0))
            self.assertLess(got["part"]["1"][2], 0.55)  # a partial model: not complete
            self.assertEqual(markers.score({}, "db", binary, 1, cache)["whole"], got["whole"])  # read from the cache

    def test_family_duplicated_is_no_flag(self):
        from eukore import package  # noqa: PLC0415
        r = dict(status="Duplicated", own_species="t", length_check="ok", gaps=0, contig_edge="-", rescued="f")
        self.assertIn("duplicated", package.copy_confidence(r)[1])
        self.assertEqual(package.copy_confidence(r, family=True), ("high", "-"))  # actin's copies: its nature
        with tempfile.TemporaryDirectory() as tmp:
            folder = os.path.join(tmp, "m")
            os.makedirs(folder)
            with open(os.path.join(folder, markers.STATE), "w") as fh:
                fh.write('{"markers": {"EF2": {"kind": "clean"}, "ACT": {"kind": "family"}}}')
            self.assertEqual(markers.kinds(folder), {"EF2": "clean", "ACT": "family"})
            state = markers.state_of(folder)  # the manifest's inputs.markers
            self.assertEqual((state["markers"]["ACT"]["kind"], state["members_md5"]), ("family", None))
        self.assertEqual(markers.kinds(None), {})
        self.assertIsNone(markers.state_of(None))
        self.assertEqual(markers.release_id("UniProt Knowledgebase Release 2026_04 consists of:"), "2026_04")

    def test_a_decoy_best_is_no_call(self):
        odb = types.SimpleNamespace(cutoff={"RPB1": 100.0}, locus_of_hmm={"RPB1": "RPB1", "decoy-K02999": "decoy-K02999"},
                                    span={})
        odb.by_locus = lambda hits: odb_mod.Odb.by_locus(odb, hits)
        aln = dict(query="RPB1~q", locus="RPB1", contig="c", start=1, end=300, strand="+", score=500.0, identity=0.9,
                   qstart=1, qend=100, exons=[(1, 300)], protein="M" * 100)
        saved = loci.hmmsearch
        try:
            loci.hmmsearch = lambda proteins, *a, **k: {n: {"RPB1": (900.0, 1.0), "decoy-K02999": (1200.0, 1.0)}
                                                         for n in proteins}
            self.assertEqual(dict(loci.confirm([aln], {"RPB1~q": 100}, odb, "db", "hmmsearch", 1, None)), {})
            loci.hmmsearch = lambda proteins, *a, **k: {n: {"RPB1": (900.0, 1.0), "decoy-K02999": (600.0, 1.0)}
                                                         for n in proteins}
            got = loci.confirm([aln], {"RPB1~q": 100}, odb, "db", "hmmsearch", 1, None)
            self.assertEqual(got["RPB1"][0]["hmm_score"], 900.0)
        finally:
            loci.hmmsearch = saved


if __name__ == "__main__":
    unittest.main()
