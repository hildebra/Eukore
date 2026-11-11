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

    def test_any_failure_of_the_fetch_is_an_oserror(self):
        # an NCBI reply without its count (AttributeError) or of bad XML: the fetch stage goes on without the markers
        # (until 2026-11-08 it stopped, every source after them unfetched)
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(
                markers, "fetch_inputs", side_effect=AttributeError("'NoneType' has no attribute 'group'")):
            with self.assertRaisesRegex(OSError, "marker set's inputs not fetched: AttributeError"):
                markers.fetch(tmp)

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
        # since 2026-11-10: Swiss-Prot's by UniProt's accession pattern; NCBI's records without a version (a PDB
        # chain, PIR's, PRF's) are no curated ones (until then rank 0: trusted, the basis of the median length)
        self.assertEqual([markers.source_rank(a) for a in ("A0A0B4J2F0", "Q9SAJ6", "4CBU_A", "S66116", "2209341A")],
                         [0, 0, 3, 3, 3])
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


RECHECK_HMMSEARCH = '''#!/usr/bin/env python3
"""hmmsearch's stand-in for markers.recheck: an EF-1a-like protein (MGKEK...) scores 500 against TEF1 and 300
against decoy-EFL, an EFL-like one (MSRQV...) 400 and 800; any other profile 10."""
import sys
args = sys.argv[1:]
dom, hmm, faa = args[args.index("--domtblout") + 1], args[-2], args[-1]
names = [line.split()[1] for line in open(hmm) if line.startswith("NAME")]
proteins, name = {}, None
for line in open(faa):
    if line.startswith(">"):
        name = line[1:].split()[0]
        proteins[name] = ""
    else:
        proteins[name] += line.strip()
scores = {"MGKEK": {"TEF1": 500.0, "decoy-EFL": 300.0}, "MSRQV": {"TEF1": 400.0, "decoy-EFL": 800.0}}
with open(dom, "w") as out:
    for p, seq in proteins.items():
        for h in names:
            s = scores.get(seq[:5], {}).get(h, 10.0)
            f = [p, "-", str(len(seq)), h, "-", "100", "1e-10", str(s), "0", "1", "1", "1e-10", "1e-10", str(s), "0",
                 "1", "100", "1", str(len(seq)), "1", str(len(seq)), "0.9", "-"]
            out.write(" ".join(f) + "\\n")
'''


class RecheckTest(unittest.TestCase):
    def test_a_copy_whose_best_profile_is_a_decoy(self):
        # since 2026-11-09: a genome's TEF1 called before the EFL decoy joined the set is EFL: low (`paralog`)
        with tempfile.TemporaryDirectory() as tmp:
            folder = os.path.join(tmp, "eukore_markers")
            for sub, name in (("hmms", "TEF1"), ("decoys", "EFL"), ("decoys", "K03267")):
                os.makedirs(os.path.join(folder, sub), exist_ok=True)
                profile = "decoy-" + name if sub == "decoys" else name
                with open(os.path.join(folder, sub, f"{name}.hmm"), "w") as fh:
                    fh.write(f"HMMER3/f [3.4]\nNAME  {profile}\nLENG  100\n//\n")
            check = markers.recheck(folder, fx.tool(tmp, "hmmsearch", RECHECK_HMMSEARCH))
            cds = {("G1", "TEF1"): fx.cds_of("MSRQV" + "A" * 40), ("G2", "TEF1"): fx.cds_of("MGKEK" + "A" * 40)}
            self.assertEqual(check(cds, {"G1": {"gencode": 1}, "G2": {"gencode": 1}}), {("G1", "TEF1")})
            self.assertEqual(check({}, {}), set())
            self.assertIsNone(markers.recheck(None))

    def test_the_scores_are_kept(self):
        # since 2026-11-10: the same copies against the same profiles are not searched again (until then at every
        # package run); other copies, or other profiles, are, and the calls are the same either way
        with tempfile.TemporaryDirectory() as tmp:
            folder = os.path.join(tmp, "eukore_markers")
            for sub, name in (("hmms", "TEF1"), ("decoys", "EFL")):
                os.makedirs(os.path.join(folder, sub))
                with open(os.path.join(folder, sub, f"{name}.hmm"), "w") as fh:
                    fh.write(f"HMMER3/f [3.4]\nNAME  {'decoy-' * (sub == 'decoys')}{name}\nLENG  100\n//\n")
            calls = os.path.join(tmp, "calls")
            binary = fx.tool(tmp, "hmmsearch", RECHECK_HMMSEARCH + f'open({calls!r}, "a").write("x")\n')
            check = markers.recheck(folder, binary)
            cds = {("G1", "TEF1"): fx.cds_of("MSRQV" + "A" * 40), ("G2", "TEF1"): fx.cds_of("MGKEK" + "A" * 40)}
            rows = {"G1": {"gencode": 1}, "G2": {"gencode": 1}}

            def searches():
                with open(calls) as fh:
                    return len(fh.read())

            def caches():
                return [f for f in os.listdir(folder) if f.startswith("recheck.")]
            self.assertEqual((check(cds, rows), searches(), len(caches())), ({("G1", "TEF1")}, 1, 1))
            self.assertEqual((check(cds, rows), searches()), ({("G1", "TEF1")}, 1))  # read, not searched
            self.assertEqual((markers.recheck(folder, "no-hmmsearch")(cds, rows)), {("G1", "TEF1")})
            more = {**cds, ("G3", "TEF1"): fx.cds_of("MSRQV" + "C" * 40)}
            self.assertEqual((check(more, rows), searches()), ({("G1", "TEF1"), ("G3", "TEF1")}, 2))
            self.assertEqual(len(caches()), 1)  # the last search's alone
            with open(os.path.join(folder, "decoys", "EFL.hmm"), "a") as fh:  # another profile set: searched
                fh.write("HMMER3/f [3.4]\nNAME  decoy-EFL2\nLENG  100\n//\n")
            self.assertEqual((check(more, rows), searches()), ({("G1", "TEF1"), ("G3", "TEF1")}, 3))
            from unittest import mock  # noqa: PLC0415
            with mock.patch.object(markers, "keep_scores", side_effect=PermissionError("read-only")):
                self.assertEqual((check(cds, rows), searches()), ({("G1", "TEF1")}, 4))  # a folder not to be written


def eutils_of(replies, summaries=(("XP_1.1", "9606", 12), ("XP_2.1", "4932", 12))):
    """markers.eutils' stand-in: an esearch of these summaries [(accession, taxid, length)], and efetch's replies in
    their order (the last one again where they run out) -> (the function, the list of efetch's calls' ids)."""
    asked = []

    def eutils(endpoint, **params):
        if endpoint == "esearch.fcgi":
            return f"<Count>{len(summaries)}</Count><WebEnv>E</WebEnv><QueryKey>1</QueryKey>"
        if endpoint == "esummary.fcgi":
            return "<eSummaryResult>" + "".join(
                f'<DocSum><Item Name="AccessionVersion">{a}</Item><Item Name="TaxId">{t}</Item>'
                f'<Item Name="Length">{n}</Item></DocSum>' for a, t, n in summaries) + "</eSummaryResult>"
        asked.append(params["id"])
        return replies[min(len(asked), len(replies)) - 1]
    return eutils, asked


WHOLE = ">XP_1.1 5'->3' exoribonuclease 2 [Homo sapiens]\nMKVLAA\nGGSSTT\n>XP_2.1 a protein [yeast]\nMKVLAAGGSSTT\n"


class NcbiProteinsTest(unittest.TestCase):
    """markers.ncbi_proteins' cache (since 2026-11-10): written only of whole answers."""

    def fetched(self, tmp, replies):
        from unittest import mock  # noqa: PLC0415
        eutils, asked = eutils_of(replies)
        with mock.patch.object(markers, "eutils", eutils):
            return markers.ncbi_proteins("a query", os.path.join(tmp, "q.faa.gz"), log=lambda *a: None), asked

    def test_a_reply_in_part_is_asked_for_again(self):
        with tempfile.TemporaryDirectory() as tmp:
            # the first reply ends after its first record, the second inside its last record's sequence
            got, asked = self.fetched(tmp, [WHOLE.split(">XP_2.1")[0], WHOLE[:-7] + "\n", WHOLE])
            self.assertEqual(asked, ["XP_1.1,XP_2.1"] * 3)
            self.assertEqual([(r["accession"], r["organism_id"], r["sequence"]) for r in got],
                             [("XP_1.1", "9606", "MKVLAAGGSSTT"), ("XP_2.1", "4932", "MKVLAAGGSSTT")])
            # a title's ">" starts no record (until 2026-11-10 XP_1.1 was cached with an empty sequence)
            self.assertEqual(markers.fasta_records(WHOLE)[0], ("XP_1.1 5'->3' exoribonuclease 2 [Homo sapiens]",
                                                               "MKVLAAGGSSTT"))
            self.assertEqual(self.fetched(tmp, ["no reply"])[0], got)  # the cache read, NCBI not asked

    def test_no_cache_of_a_reply_that_stays_in_part(self):
        for reply in (WHOLE.split(">XP_2.1")[0], WHOLE[:-7] + "\n", ""):
            with tempfile.TemporaryDirectory() as tmp:
                with self.assertRaisesRegex(OSError, "answered in part 3 times"):
                    self.fetched(tmp, [reply])
                self.assertEqual(os.listdir(tmp), [])  # no cache a later build would read as whole, no half file


BUILD_HMMSEARCH = RECHECK_HMMSEARCH.replace('"TEF1"', '"K1"').replace('"decoy-EFL"', '"K2"')
SPECIES = (4932, 27291, 28985, 5476, 40410, 6253, 9606)  # fixtures' taxdump's


class BuildTest(unittest.TestCase):
    """markers.build (since 2026-11-10): the set written beside the folder and put in its place when whole, and what
    ready() holds it against."""

    def world(self, tmp):
        """A sources folder whose marker folder holds the build's inputs (nothing is fetched): KOfam's list and the
        profiles K1 (the marker MK1's) and K2 (its decoy), a Swiss-Prot of seven species' MK1 -> (sources, the
        marker table with the decoy, the one without, a decoy table that is not there, hmmsearch's stand-in)."""
        import random  # noqa: PLC0415
        src = os.path.join(tmp, "sources")
        folder = os.path.join(src, "eukore_markers")
        os.makedirs(os.path.join(folder, "kofam"))
        fx.taxdump(os.path.join(src, "taxdump.tar.gz"))
        tables = []
        for name, decoys in (("with.tsv", "K2"), ("without.tsv", "")):
            tables.append(os.path.join(tmp, name))
            with open(tables[-1], "w") as fh:
                fh.write(f"locus\tname\ttargets\tdecoys\tkind\nMK1\ta marker\tK1\t{decoys}\tclean\n")
        with gzip.open(os.path.join(folder, "ko_list.gz"), "wt") as fh:
            fh.write("knum\tthreshold\tscore_type\tdefinition\nK1\t100\tfull\tthe marker\nK2\t100\tfull\tits paralog\n")
        for ko in ("K1", "K2"):
            with open(os.path.join(folder, "kofam", f"{ko}.hmm"), "w") as fh:
                fh.write(f"HMMER3/f [3.4]\nNAME  {ko}\nLENG  100\n//\n")
        rnd = random.Random(3)
        with gzip.open(os.path.join(folder, "uniprot_sprot.fasta.gz"), "wt") as fh:
            for i, taxid in enumerate(SPECIES):
                seq = "MGKEK" + "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(100))
                fh.write(f">sp|P0000{i}|MK1_X A marker OS=a species OX={taxid} PE=1 SV=1\n{seq}\n")
        with open(os.path.join(folder, "profiles.tar.gz"), "w") as fh:  # a copy of KOfam's archive: kept
            fh.write("x")
        binary = fx.tool(tmp, "hmmsearch", BUILD_HMMSEARCH)
        return src, tables[0], tables[1], os.path.join(tmp, "no_decoys.tsv"), binary

    def build(self, src, table, binary, decoys):
        from unittest import mock  # noqa: PLC0415
        with mock.patch.object(markers, "MIN_MEMBERS", 3), mock.patch.object(markers, "release_of",
                                                                             return_value="2026_04"):
            return markers.build(src, table, binary, 1, lambda *a: None, decoys)

    def ready(self, src, table, decoys=None):
        from unittest import mock  # noqa: PLC0415
        with mock.patch.object(markers, "MIN_MEMBERS", 3):
            return markers.ready(src, table, decoys)

    @staticmethod
    def text(*path):
        with open(os.path.join(*path)) as fh:
            return fh.read()

    def test_the_set_is_put_in_place_whole(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, table, _without, decoys, binary = self.world(tmp)
            folder = self.build(src, table, binary, decoys)
            self.assertEqual(folder, os.path.join(src, "eukore_markers"))
            self.assertEqual(sorted(os.listdir(src)), ["eukore_markers", "taxdump.tar.gz"])  # nothing left beside it
            self.assertTrue(set(markers.SET) <= set(os.listdir(folder)))
            # the inputs and caches are the new folder's
            self.assertTrue({"ko_list.gz", "kofam", "kofam.hmm", "uniprot_sprot.fasta.gz", "profiles.tar.gz"}
                            <= set(os.listdir(folder)))
            self.assertTrue(any(f.startswith("scores.") for f in os.listdir(folder)))
            self.assertIn("NAME  MK1\n", self.text(folder, "hmms", "MK1.hmm"))
            self.assertIn("NAME  decoy-K2\n", self.text(folder, "decoys", "K2.hmm"))
            self.assertEqual(self.text(folder, "scores_cutoff"), "MK1\t250.0\n")
            self.assertEqual(markers.kinds(folder), {"MK1": "clean"})
            self.assertEqual(markers.state_of(folder)["markers"]["MK1"]["members"], len(SPECIES))
            self.assertEqual((self.ready(src, table), self.ready(src, table, decoys)), (folder, folder))

    def test_a_build_that_fails_leaves_the_old_set(self):
        # until 2026-11-10 hmms/ and decoys/ were emptied in the folder itself before the new set was made: a build
        # that failed after it left the old markers.json over no profiles, which ready() took
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src, table, without, decoys, binary = self.world(tmp)
            folder = self.build(src, table, binary, decoys)
            was = {name: self.text(folder, *name.split("/")) for name in ("hmms/MK1.hmm", "decoys/K2.hmm",
                                                                         "scores_cutoff", "markers.json")}
            with mock.patch.object(markers, "foreign", side_effect=RuntimeError("killed")):
                with self.assertRaisesRegex(RuntimeError, "killed"):
                    self.build(src, without, binary, decoys)
            self.assertEqual({name: self.text(folder, *name.split("/")) for name in was}, was)
            self.assertEqual(sorted(os.listdir(src)), ["eukore_markers", "taxdump.tar.gz"])  # no half set left
            self.assertEqual((self.ready(src, table), self.ready(src, without)), (folder, None))
            self.assertTrue({"ko_list.gz", "kofam", "uniprot_sprot.fasta.gz"} <= set(os.listdir(folder)))
            self.build(src, without, binary, decoys)  # and one that ends: the new table's set, its decoy gone
            self.assertEqual((os.listdir(os.path.join(folder, "decoys")), self.ready(src, table),
                              self.ready(src, without, decoys)), ([], None, folder))
            self.assertTrue({"ko_list.gz", "kofam", "uniprot_sprot.fasta.gz", "profiles.tar.gz"}
                            <= set(os.listdir(folder)))

    def test_stale_by_the_rules_and_the_inputs(self):
        # until 2026-11-10 the table's (and the decoy table's) MD5 alone: another PER_TAXON or threshold, another
        # Swiss-Prot, did not make the set stale. The other stages (no decoys given) take it as before
        import json  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src, table, _without, decoys, binary = self.world(tmp)
            folder = self.build(src, table, binary, decoys)
            for name, value in (("PER_TAXON", 5), ("NCBI_SCOPE", "refseq[filter]"), ("CUTOFF_SHARE", 0.4),
                                ("RECIPE", markers.RECIPE + 1)):
                with mock.patch.object(markers, name, value):
                    self.assertEqual((self.ready(src, table, decoys), self.ready(src, table)), (None, folder), name)
            self.assertEqual(self.ready(src, table, decoys), folder)
            sprot = os.path.join(folder, "uniprot_sprot.fasta.gz")
            with gzip.open(sprot, "at") as fh:  # another release put there
                fh.write(">sp|P00009|MK1_Y A marker OS=a species OX=952 PE=1 SV=1\nMGKEKAAAA\n")
            self.assertEqual((self.ready(src, table, decoys), self.ready(src, table)), (None, folder))
            self.assertEqual(self.build(src, table, binary, decoys), folder)  # built from it
            self.assertEqual(self.ready(src, table, decoys), folder)
            os.remove(sprot)  # an input deleted: fetched by the next fetch stage, built from by the process stage
            self.assertEqual((self.ready(src, table, decoys), self.ready(src, table)), (None, folder))
            state = json.loads(self.text(folder, "markers.json"))
            with open(os.path.join(folder, "markers.json"), "w") as fh:  # a set of before: neither recipe nor inputs
                json.dump({k: v for k, v in state.items() if k not in ("recipe", "inputs")}, fh)
            self.assertEqual((self.ready(src, table, decoys), self.ready(src, table)), (None, folder))

    def test_what_a_cut_build_left(self):
        from eukore import util  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src, table, _without, decoys, binary = self.world(tmp)
            folder = self.build(src, table, binary, decoys)

            def dead(name, n=1):  # a run's folder on a host that is heard of no more
                return os.path.join(src, f".gone-{n}.{name}.partial")

            def aged(path):
                os.utime(path, (0, 0))
                return path
            # cut between the swap's two renames: no folder, the new set whole beside it, and the old one
            os.makedirs(os.path.join(dead("eukore_markers.old"), "hmms"))
            os.rename(folder, dead("eukore_markers"))
            aged(dead("eukore_markers"))
            aged(dead("eukore_markers.old"))
            markers.recover(src, lambda *a: None)
            self.assertEqual(sorted(os.listdir(src)), ["eukore_markers", "taxdump.tar.gz"])
            self.assertEqual(self.ready(src, table, decoys), folder)
            # cut while the caches moved into the new set: they are taken back, the half swap's folder removed
            os.makedirs(os.path.join(dead("eukore_markers", 2), "hmms"))
            os.rename(os.path.join(folder, "ko_list.gz"), os.path.join(dead("eukore_markers", 2), "ko_list.gz"))
            aged(dead("eukore_markers", 2))
            mine = util.part_of(folder)  # a build going on (this process's): not touched
            os.makedirs(mine)
            markers.recover(src, lambda *a: None)
            self.assertEqual(sorted(os.listdir(src)), sorted(["eukore_markers", "taxdump.tar.gz",
                                                              os.path.basename(mine)]))
            self.assertIn("ko_list.gz", os.listdir(folder))
            self.assertEqual(self.ready(src, table, decoys), folder)

    def test_the_fetch_fetches_and_the_build_computes(self):
        # since 2026-11-11 the fetch stage fetches the set's inputs alone (until then it built the set too: its
        # hmmsearch of ~230,000 candidates ran on the download node), and the build, the process stage's, fetches
        # nothing
        import json  # noqa: PLC0415
        import shutil  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src, table, _without, _none, binary = self.world(tmp)
            folder = os.path.join(src, "eukore_markers")
            decoys = os.path.join(tmp, "decoys.tsv")
            with open(decoys, "w") as fh:
                fh.write("locus\tdecoy\tname\tncbi\nMK1\tEFL\tits look-alike\tefl[title]\n")
            sprot = os.path.join(folder, "uniprot_sprot.fasta.gz")
            shutil.move(sprot, os.path.join(tmp, "sprot.gz"))
            urls, queries = [], []

            def download(url, dest):
                urls.append(url)
                if url == markers.SWISSPROT:
                    shutil.copy(os.path.join(tmp, "sprot.gz"), dest)
                else:
                    with open(dest, "w") as fh:
                        fh.write("UniProt Knowledgebase Release 2026_04 consists of:\n")

            def proteins(term, path, log=print, scope=markers.NCBI_SCOPE):
                queries.append((term, scope))
                with gzip.open(path, "wt") as fh:
                    fh.write(">XP_9.1 9606\nMSRQVAAAAA\n")
            with mock.patch.object(markers, "score", side_effect=AssertionError("scored")), \
                    mock.patch.object(markers, "ncbi_proteins", side_effect=proteins):
                with self.assertRaisesRegex(OSError, r"2 of its inputs not fetched into .*--stage fetch"):
                    markers.build(src, table, binary, 1, lambda *a: None, decoys)  # nothing fetched there
                self.assertEqual(markers.fetch_inputs(src, table, download, lambda *a: None, decoys), folder)
                self.assertEqual((urls[1:], queries), ([markers.SWISSPROT], [("efl[title]", markers.DECOY_SCOPE)]))
                self.assertTrue(urls[0].endswith("/reldate.txt"))
                self.assertEqual((markers.inputs_missing(src, table, decoys), self.ready(src, table)), ([], None))
                urls.clear()
                markers.fetch_inputs(src, table, download, lambda *a: None, decoys)  # all there: nothing fetched
                self.assertEqual((urls, len(queries)), ([], 1))
            self.assertEqual(markers.release_of(os.path.join(folder, markers.RELDATE)), "2026_04")
            self.assertEqual(markers.release_of(os.path.join(folder, "none.txt")), "?")
            with mock.patch.object(markers, "eutils", side_effect=AssertionError("NCBI asked")), \
                    mock.patch.object(markers.urllib.request, "urlopen", side_effect=AssertionError("fetched")), \
                    mock.patch.object(markers.urllib.request, "urlretrieve", side_effect=AssertionError("fetched")):
                self.assertEqual(self.build(src, table, binary, decoys), folder)
            self.assertEqual(self.ready(src, table, decoys), folder)
            self.assertIn(markers.decoy_proteins({"decoy": "EFL", "ncbi": "efl[title]"}),
                          json.loads(self.text(folder, "markers.json"))["inputs"])

    def test_one_run_builds_the_set_the_others_wait(self):
        # the process stage's (since 2026-11-11): a SLURM array's tasks start together, one builds, the others wait
        # for its set and take it
        import threading  # noqa: PLC0415
        from unittest import mock  # noqa: PLC0415
        from eukore import util  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src, table, _without, decoys, binary = self.world(tmp)
            folder = os.path.join(src, "eukore_markers")
            said, got = [], []
            with mock.patch.object(markers, "MIN_MEMBERS", 3), mock.patch.object(markers, "release_of",
                                                                                 return_value="2026_04"):
                with util.single_run(src, "markers"):  # another run's build
                    waiter = threading.Thread(target=lambda: got.append(
                        markers.prepare(src, binary, 1, said.append, table, decoys, wait=0.02)))
                    waiter.start()
                    waiter.join(0.3)
                    self.assertTrue(waiter.is_alive())  # it waits
                    markers.build(src, table, binary, 1, lambda *a: None, decoys)
                waiter.join(10)
                self.assertEqual(got, [folder])
                self.assertEqual(len(said), 1)  # its one line: waiting; it built nothing
                self.assertIn("another run builds the set", said[0])
                self.assertFalse(os.path.exists(os.path.join(src, ".markers.lock")))
                self.assertEqual(markers.prepare(src, binary, 1, said.append, table, decoys), folder)  # built: taken
            # an input the fetch stage never fetched: said, no lock or half set left
            os.remove(os.path.join(folder, "markers.json"))
            os.remove(os.path.join(folder, "uniprot_sprot.fasta.gz"))
            with self.assertRaisesRegex(OSError, r"marker set not built: 1 of its inputs not fetched .*--stage fetch"):
                markers.prepare(src, binary, 1, said.append, table, decoys)
            self.assertEqual(sorted(os.listdir(src)), ["eukore_markers", "taxdump.tar.gz"])


if __name__ == "__main__":
    unittest.main()
