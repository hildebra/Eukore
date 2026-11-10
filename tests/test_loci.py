#!/usr/bin/env python3
"""Tests of Eukore's own locus calls (odb.py, loci.py; since 2026-11-05): the dataset read (cutoffs, the HMMs' own
names, proteins by species, the nearest species' queries first); miniprot's GFF with --trans parsed; overlapping
alignments one region; a region confirmed by its locus' HMM at the cutoff; a locus Single, Duplicated (two complete
regions), Fragmented (a short HMM alignment) or Missing; the copy whole on its own strand with its exons and introns, its
soft-masked and low-complexity runs; the windows around it.

    python3 -m unittest discover -s tests -t tests
"""

import concurrent.futures
import gzip
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import loci, ncbi as ncbi_mod, odb as odb_mod, placement, qc  # noqa: E402
import fixtures as fx  # noqa: E402

INDEXING_MINIPROT = '''#!/usr/bin/env python3
"""A miniprot that knows -d (its saved index: here a line of what it was made of, written to the -d file as miniprot
writes it) and notes each command line in a log; every other call is fixtures.miniprot's, which reads its proteins
alone. The mode file's word, where it is there: "fail" (-d writes half an index and fails, as one killed), "empty"
(-d writes nothing and exits 0, as miniprot does where it cannot open its file)."""
import os
import sys
LOG, MODE, REAL = %r, %r, %r
args = sys.argv[1:]
with open(LOG, "a") as fh:
    fh.write("\\t".join(args) + "\\n")
if "-d" in args:
    mode = open(MODE).read().strip() if os.path.exists(MODE) else ""
    if mode == "empty":
        sys.exit(0)
    with open(args[args.index("-d") + 1], "w") as fh:
        fh.write("half an index, ")
        if mode == "fail":
            sys.stderr.write("reading the genome\\nkilled while indexing\\n")
            sys.exit(1)
        fh.write("whole: " + " ".join(a for a in args if a != "/dev/stdout") + "\\n")
    sys.exit(0)
os.execv(REAL, [REAL] + args)
'''


def indexing_miniprot(folder, gff):
    """A miniprot stand-in that saves an index -> (its path, its log of command lines, its mode file)."""
    log, mode = os.path.join(folder, "mp.log"), os.path.join(folder, "mp.mode")
    return fx.tool(folder, "miniprot_d", INDEXING_MINIPROT % (log, mode, fx.miniprot(folder, gff))), log, mode


def commands(log, clear=True):
    """The command lines a stand-in noted (their arguments), the log emptied."""
    with open(log) as fh:
        out = [line.rstrip("\n").split("\t") for line in fh]
    if clear:
        open(log, "w").close()
    return out


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
        got = loci.call(self.path, fx.QUERIES, self.odb, fx.miniprot(self.dir, self.gff), fx.hmmsearch(self.dir), 1,
                        self.dir)
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
        got = loci.call(self.path, fx.QUERIES, self.odb, mp, fx.hmmsearch(self.dir), 1, self.dir)
        self.assertEqual((got["loci"]["2at2759"]["status"], round(got["loci"]["2at2759"]["hmm_coverage"], 2)),
                         ("Fragmented", 0.75))
        other = self.gff.replace(fx.PROTEINS["2at2759"], fx.PROTEINS["3at2759"])  # another locus' protein: not 2
        got = loci.call(self.path, fx.QUERIES, self.odb, fx.miniprot(self.dir, other), fx.hmmsearch(self.dir), 1,
                        self.dir, rescue_missing=False)
        self.assertEqual(got["loci"]["2at2759"]["status"], "Missing")
        failing = fx.tool(self.dir, "mp_fail", "#!/bin/sh\necho 'out of memory' >&2\nexit 1\n")
        got = loci.call(self.path, fx.QUERIES, self.odb, failing, fx.hmmsearch(self.dir), 1, self.dir)
        self.assertTrue(got["failed"].startswith("miniprot failed: out of memory"))

    def test_a_profile_s_coverage_the_union_of_its_domains(self):
        # three domains of a profile of 120: 40-80, 1-50 (overlapping), 100-120 (apart): 101 positions, each once
        # (since 2026-11-10 by the merged spans; the same value as the set of positions until then)
        rows = "".join(f"p - 100 7 - 120 1e-9 55.0 0 {i} 3 1e-5 1e-5 20.0 0 {a} {b} 1 9 1 9 0.9 -\n"
                       for i, (a, b) in enumerate([(40, 80), (1, 50), (100, 120)], 1))
        tool = fx.tool(self.dir, "hs", "#!/bin/sh\nwhile [ \"$1\" != --domtblout ]; do shift; done\n"
                                       f"printf '%s' '{rows}' > \"$2\"\n")
        got = loci.hmmsearch({"p": "MKV"}, os.path.join(self.dir, "x.hmm"), tool, 1, self.dir)
        self.assertEqual(got, {"p": {"7": (55.0, 101 / 120)}})

    def test_frameshifts_left_out_of_the_coding_sequence(self):
        # miniprot's CIGAR: amino acids (M, D: 3 bases each), a protein insertion (I: none), a phase-0 intron (N: not
        # coding), a phase-1 or 2 one (U, V: its split codon's 3 bases), frameshifts (F, G: their bases, out of frame)
        self.assertEqual(loci.frameshifts("3M1F2M35V4M2G1I1M"), [(9, 1), (31, 2)])
        self.assertEqual(loci.frameshifts("100M500U50M1F49M"), [(453, 1)])
        self.assertEqual(loci.frameshifts("100M500N50M1F49M"), [(450, 1)])
        self.assertEqual(loci.in_frame("AAACCCGGGtTTTAAA", [(9, 1)]), "AAACCCGGGTTTAAA")
        self.assertEqual(loci.in_frame("ACGT", []), "ACGT")
        text = ("##PAF\t1at2759~q\t120\t0\t120\t+\tc1\t9000\t99\t500\t1\t1\t0\tAS:i:1\tcg:Z:3M1F37M\tcs:Z::1\n"
                "##STA\tMAK\nc1\tminiprot\tmRNA\t100\t500\t9\t+\t.\tID=MP1;Identity=0.9;Target=1at2759~q 1 120\n"
                "c1\tminiprot\tCDS\t100\t500\t9\t+\t0\tParent=MP1;Target=1at2759~q 1 120\n")
        self.assertEqual(loci.parse_alignments(text)[0]["shifts"], [(9, 1)])
        # the stop codon miniprot writes apart from the CDS is the last exon's end, on either strand
        self.assertEqual(loci.parse_alignments(text + "c1\tminiprot\tstop_codon\t501\t503\t9\t+\t0\tParent=MP1\n")[0]
                         ["exons"], [(100, 503)])
        rev = text.replace("\t+\t", "\t-\t") + "c1\tminiprot\tstop_codon\t97\t99\t9\t-\t0\tParent=MP1\n"
        self.assertEqual(loci.parse_alignments(rev)[0]["exons"], [(97, 500)])
        # on the copy: a frameshift's two bases across an intron (exons 0-10 and 20-40), the copy's coding sequence
        # (qc.spliced, the package's) the same as the call's
        copy = "".join("ACGT"[i % 4] for i in range(40))
        exons = [(0, 10), (20, 40)]
        feats = loci.frameshift_features(exons, [(9, 2)])
        self.assertEqual(feats, [("frameshift", 9, 10), ("frameshift", 20, 21)])
        joined = copy[0:10] + copy[20:40]
        self.assertEqual(qc.spliced(copy, [("exon", s, e) for s, e in exons] + feats), loci.in_frame(joined, [(9, 2)]))

    def test_the_coding_sequence_checked(self):
        # since 2026-11-09: a copy's coding sequence as it ships, against the protein its HMM confirmed
        protein = "MKLVWHGTEERRASDPLLKKEQWTYV"
        cds = fx.cds_of(protein)  # its stop codon last
        self.assertEqual(loci.cds_state(cds, protein), "clean")
        self.assertEqual(loci.cds_state(cds[:-3], protein), "clean")  # without the stop that ends it
        self.assertEqual(loci.cds_state(cds[:30] + "TAA" + cds[30:], protein), "stops")  # one inside
        self.assertNotEqual(loci.cds_state(cds[:12] + cds[13:], protein), "clean")  # out of frame after a base
        self.assertEqual(loci.cds_state(fx.cds_of("MKLVWHGTCCCCCCCCCCPLLKKEQWTYV"), protein), "mismatch")
        self.assertEqual(loci.cds_state(cds, ""), "clean")  # nothing to set it beside: its stops alone
        # a region's clean alignment before its best HMM score's broken one (a base lost: out of frame); complete
        # where that is
        contigs = {"c1": "A" * 30 + cds + "A" * 30}
        broken = dict(contig="c1", start=31, end=30 + len(cds), strand="+", exons=[(31, 40), (42, 30 + len(cds))],
                      protein=protein, hmm_score=90.0, score=50.0, hmm_coverage=0.9)
        good = dict(broken, exons=[(31, 30 + len(cds))], hmm_score=80.0)
        fragment = dict(good, hmm_score=85.0, hmm_coverage=0.5)
        got = loci.cleanest(dict(broken, alts=[fragment, good]), contigs, 1)
        self.assertEqual((got["cds_check"], got["hmm_score"], got["exons"]), ("recovered", 80.0, good["exons"]))
        self.assertEqual(loci.cleanest(dict(good, alts=[broken]), contigs, 1)["cds_check"], "clean")
        alone = loci.cleanest(dict(broken, alts=[fragment]), contigs, 1)  # no complete clean one: as it was
        self.assertEqual(alone["hmm_score"], 90.0)
        self.assertIn(alone["cds_check"], ("stops", "mismatch"))
        # the locus' copy: a clean one before a broken one of a higher score
        chosen, statuses = loci.choose({"1": [alone, dict(good, cds_check="clean")]}, ["1", "2"])
        self.assertEqual((chosen["1"]["hmm_score"], statuses["1"][0], len(statuses["1"][1])), (80.0, "Duplicated", 1))

    def test_an_annotated_transcript_read_whole(self):
        # miniprot's span 1000-5000 misses NCBI's first exon (500-600): the transcript keeps it; one off the spans not
        from eukore import annotation  # noqa: PLC0415
        gff = os.path.join(self.dir, "t.gff")
        with open(gff, "w") as fh:
            fh.write("c1\tRefSeq\tCDS\t500\t600\t.\t+\t0\tParent=rna-1;gene=ONE\n"
                     "c1\tRefSeq\tCDS\t1200\t1500\t.\t+\t0\tParent=rna-1;gene=ONE\n"
                     "c1\tRefSeq\tCDS\t4800\t6000\t.\t+\t0\tParent=rna-1;gene=ONE\n"
                     "c1\tRefSeq\tCDS\t9000\t9100\t.\t+\t0\tParent=rna-2;gene=TWO\n")
        got = annotation.read_cds(gff, {"c1": [(1000, 5000)]})
        self.assertEqual(list(got), ["rna-1"])
        self.assertEqual(got["rna-1"]["cds"], [(500, 600), (1200, 1500), (4800, 6000)])

    def test_a_genome_in_pieces_filtered_over_them_all(self):
        # piece 1's alignment of q scores 900, piece 2's 450 (a weak paralogue, its piece's best): --outs over both
        hits = [dict(query="q", score=900.0), dict(query="q", score=450.0), dict(query="r", score=10.0)]
        hits += [dict(query="p", score=1000.0 - i) for i in range(loci.OUTN + 3)]
        got = loci.best_of_each(hits)
        self.assertEqual([a["score"] for a in got if a["query"] == "q"], [900.0])
        self.assertEqual([a["score"] for a in got if a["query"] == "r"], [10.0])
        self.assertEqual(len([a for a in got if a["query"] == "p"]), loci.OUTN)

    def test_annotated_models_taken_where_their_hmm_holds(self):
        a, b = self.spans["1at2759"]
        c, d = self.spans["2at2759"]
        gff = os.path.join(self.dir, "ann.gff.gz")
        with gzip.open(gff, "wt") as fh:  # locus 1: one base before the gene, phase 1 (its translation the protein)
            fh.write("##gff-version 3\n"
                     f"c1\tRefSeq\tCDS\t{a - 1}\t{a + 179}\t.\t+\t1\tID=c1;Parent=rna-1;gene=ONE\n"
                     f"c1\tRefSeq\tCDS\t{a + 480}\t{b}\t.\t+\t0\tID=c1;Parent=rna-1;gene=ONE\n"
                     f"c1\tRefSeq\tCDS\t{c}\t{d - 1}\t.\t-\t0\tID=c2;Parent=rna-2;gene=TWO\n")  # locus 2: out of frame
        got = loci.call(self.path, fx.QUERIES, self.odb, fx.miniprot(self.dir, self.gff), fx.hmmsearch(self.dir), 1,
                        self.dir, annotation=(gff, None))
        one, two = got["loci"]["1at2759"], got["loci"]["2at2759"]
        self.assertEqual((one["model"], one["transcript"], one["start"], one["exons"]),
                         ("ncbi", "rna-1", a, [(a, a + 179), (a + 480, b)]))
        self.assertEqual((one["spliced"], one["miniprot"]["start"]), (fx.cds_of(fx.PROTEINS["1at2759"]), a))
        self.assertEqual((two["model"], two["transcript"], two["start"], two["end"]), ("miniprot", "-", c, d))
        plain = loci.call(self.path, fx.QUERIES, self.odb, fx.miniprot(self.dir, self.gff), fx.hmmsearch(self.dir), 1,
                          self.dir)
        self.assertEqual(plain["loci"]["1at2759"]["model"], "miniprot")

    def test_a_row_s_status_its_shipped_model_s(self):
        # since 2026-11-10 (the audit's B8): miniprot's model of locus 2 a fragment (3/4 of the profile), the
        # annotation's the whole gene: taken, and the locus Single (until then Fragmented beside a coverage of 1);
        # locus 1's annotated model its first exon alone (half the profile) where miniprot's is complete: not taken
        # (until then taken: Duplicated beside a coverage of 0.5)
        a, _b = self.spans["1at2759"]
        c, d = self.spans["2at2759"]
        short = self.gff.replace(fx.PROTEINS["2at2759"], fx.PROTEINS["2at2759"][:90])
        gff = os.path.join(self.dir, "ann.gff.gz")
        with gzip.open(gff, "wt") as fh:
            fh.write("##gff-version 3\n"
                     f"c1\tRefSeq\tCDS\t{a}\t{a + 179}\t.\t+\t0\tID=c1;Parent=rna-1;gene=ONE\n"
                     f"c1\tRefSeq\tCDS\t{c}\t{d}\t.\t-\t0\tID=c2;Parent=rna-2;gene=TWO\n")
        args = (self.path, fx.QUERIES, self.odb, fx.miniprot(self.dir, short), fx.hmmsearch(self.dir), 1, self.dir)
        plain = loci.call(*args)["loci"]
        got = loci.call(*args, annotation=(gff, None))["loci"]
        was, two = plain["2at2759"], got["2at2759"]
        self.assertEqual((was["status"], round(was["hmm_coverage"], 2), was["model"]), ("Fragmented", 0.75, "miniprot"))
        self.assertEqual((two["status"], two["copies"], two["hmm_coverage"], two["model"], two["transcript"]),
                         ("Single", 1, 1.0, "ncbi", "rna-2"))
        self.assertEqual(two["spliced"], fx.cds_of(fx.PROTEINS["2at2759"]))
        # how the gene was found stays miniprot's alignment's: its query, identity, coverage of the query, score
        self.assertEqual([two[f] for f in ("query", "identity", "coverage", "score")],
                         [plain["2at2759"][f] for f in ("query", "identity", "coverage", "score")])
        one = got["1at2759"]
        self.assertEqual((one["status"], one["hmm_coverage"], one["model"], one["transcript"]),
                         ("Duplicated", 1.0, "miniprot", "-"))

    def test_a_missed_gene_rescued(self):
        # locus 3's gene after the others: miniprot's first pass (the fixed GFF) does not align it, its HMM finds it in
        # the six-frame ORFs, and the sensitive pass on that place aligns it (the stand-in searches with -M0)
        rng = __import__("random").Random(9)
        pad = "".join(rng.choice("ACGT") for _ in range(400))
        gene = fx.cds_of(fx.PROTEINS["3at2759"])
        path = os.path.join(self.dir, "g3.fna.gz")
        with gzip.open(path, "wt") as fh:
            fh.write(f">c1 a test contig\n{self.seq + pad + gene + pad}\n")
        start = len(self.seq) + len(pad) + 1
        args = (path, fx.QUERIES | {k.replace("1at", "3at"): fx.PROTEINS["3at2759"] for k in fx.QUERIES
                                     if k.startswith("1at")}, self.odb, fx.miniprot(self.dir, self.gff),
                fx.hmmsearch(self.dir), 1, self.dir)
        translate = fx.esl_translate(self.dir)
        got = loci.call(*args, translate=translate)
        three = got["loci"]["3at2759"]
        self.assertEqual((three["status"], three["rescued"], three["start"], three["end"], three["contig"]),
                         ("Single", True, start, start + len(gene) - 1, "c1"))  # its codons, its stop
        self.assertEqual(three["spliced"], gene)
        self.assertEqual((got["loci"]["1at2759"]["rescued"], got["loci"]["2at2759"]["status"]), (False, "Single"))
        self.assertEqual(loci.call(*args, translate=translate, rescue_missing=False)["loci"]["3at2759"]["status"],
                         "Missing")
        self.assertEqual(loci.call(*args, translate="no-such-tool")["loci"]["3at2759"]["status"], "Missing")

    def test_a_gene_rescued_after_an_assembly_gap(self):
        # a gap of 310 kb of N right before locus 3's gene: six ORFs over hmmsearch's 100,000 residues, split at their
        # X runs; the gene's piece found and placed (Takifugu rubripes' build failed on such a gap, 2026-11-06)
        rng = __import__("random").Random(9)
        pad = "".join(rng.choice("ACGT") for _ in range(400))
        gene = fx.cds_of(fx.PROTEINS["3at2759"])
        path = os.path.join(self.dir, "g3.fna.gz")
        with gzip.open(path, "wt") as fh:
            fh.write(f">c1 a test contig\n{self.seq + pad + 'N' * 310_000 + gene + pad}\n")
        start = len(self.seq) + len(pad) + 310_000 + 1
        queries = fx.QUERIES | {k.replace("1at", "3at"): fx.PROTEINS["3at2759"] for k in fx.QUERIES
                                if k.startswith("1at")}
        got = loci.call(path, queries, self.odb, fx.miniprot(self.dir, self.gff), fx.hmmsearch(self.dir), 1, self.dir,
                        translate=fx.esl_translate(self.dir))
        self.assertIsNone(got["failed"])
        three = got["loci"]["3at2759"]
        self.assertEqual((three["status"], three["rescued"], three["start"], three["end"]),
                         ("Single", True, start, start + len(gene) - 1))  # its stop codon in the copy
        self.assertEqual(got["loci"]["1at2759"]["status"], "Duplicated")

    def test_orfs_split_at_gaps(self):
        # esl-translate's own output (a contig of 24 bases, 57 N, 42 bases): both strands' coords kept on the pieces
        text = (">orf2 source=ctg coords=1..120 length=40 frame=1 desc=\n"
                "MKPGFKPGXXXXXXXXXXXXXXXXXXXKPGFKPGFKPGFK\n"
                ">orf6 source=ctg coords=123..1 length=41 frame=4 desc=\n"
                "LFKPGFKPGFKPGFXXXXXXXXXXXXXXXXXXXPGFKPGFH\n")
        self.assertEqual(loci.orf_pieces("MKPGFKPGXXXXXXXXXXXXXXXXXXXKPGFKPGFKPGFK", 1, 120), [])  # pieces under 20
        old = loci.ORF_MIN
        loci.ORF_MIN = 5
        try:
            out = __import__("io").StringIO()
            loci.split_orfs(text.splitlines(True), out)
        finally:
            loci.ORF_MIN = old
        self.assertEqual(out.getvalue(), ">orf2.0 source=ctg coords=1..24 length=8\nMKPGFKPG\n"
                                         ">orf2.1 source=ctg coords=82..120 length=13\nKPGFKPGFKPGFK\n"
                                         ">orf6.0 source=ctg coords=123..82 length=14\nLFKPGFKPGFKPGF\n"
                                         ">orf6.1 source=ctg coords=24..1 length=8\nPGFKPGFH\n")
        seq = "M" * 30 + "X" * 9 + "M" * 30  # a short run of X (an ambiguous base or two) kept
        self.assertEqual(loci.orf_pieces(seq, 1, 3 * len(seq)), [(seq, 1, 3 * len(seq))])
        self.assertEqual(loci.orf_pieces("M" * (loci.HMMSEARCH_MAX + 1), 1, 3 * (loci.HMMSEARCH_MAX + 1)), [])

    def test_located(self):
        def hit(start, score, hf, ht, protein, contig="c1", strand="+"):
            return dict(contig=contig, strand=strand, start=start, end=start + 3 * len(protein) - 1, score=score,
                        hmm_from=hf, hmm_to=ht, protein=protein)
        hits = {"1at2759": [hit(5000, 80.0, 40, 90, "BBBB"), hit(4000, 30.0, 1, 35, "AAAA"),
                            hit(4500, 20.0, 30, 60, "XXXX"),                 # overlaps the best on the HMM
                            hit(90000, 25.0, 95, 120, "FAR"),                # too far from the best
                            hit(5200, 50.0, 95, 120, "OTHER", strand="-"),   # the other strand
                            hit(5300, 40.0, 95, 120, "CCCC")]}
        self.assertEqual(loci.located(hits), {"1at2759": ("c1", 4000, 5311, "AAAABBBBCCCC")})

    def test_masks(self):
        self.assertEqual(loci.soft_runs("ACgtacgtacgtacgtacgtacgtAC", 20), [(2, 24)])
        self.assertEqual(loci.soft_runs("ACgtAC", 20), [])
        self.assertEqual(loci.on_copy([(0, 10)], 101, 200, "-"), [(90, 100)])
        self.assertEqual(loci.dust("ACGT" * 5), [])
        rng = __import__("random").Random(1)
        seq = "".join(rng.choice("ACGT") for _ in range(300)) + "A" * 200 + "".join(rng.choice("ACGT") for _ in range(300))
        runs = loci.dust(seq)
        self.assertTrue(runs and runs[0][0] <= 320 and runs[-1][1] >= 480 and all(e - s <= 400 for s, e in runs))

    # miniprot's index saved (since 2026-11-10; the audit's S1): a genome called more than once indexed once

    def aligned(self, target, table=1, options=("-I",)):
        """The alignment's command line as loci.align issues it on a target (without the proteins' file)."""
        return ["-t", "1", *options, "-T", str(table), "--gff", "--trans", f"--outs={loci.OUTS}",
                f"--outn={loci.OUTN}", target]

    def test_the_index_made_once_for_a_genome_s_calls(self):
        mp, log, _mode = indexing_miniprot(self.dir, self.gff)
        idx = os.path.join(self.dir, "idx")
        args = (self.odb, mp, fx.hmmsearch(self.dir), 1, self.dir)
        plain = loci.call(self.path, fx.QUERIES, *args, rescue_missing=False)
        # without index_dir: the command as it was, on the FASTA, no index made or kept
        (old,) = commands(log)
        self.assertEqual((old[:-1], os.path.basename(old[-1])), (self.aligned(self.path), "q.faa"))
        self.assertFalse(os.path.exists(idx))
        # the main call and two hold-outs (other proteins, the same genome): one -d, each alignment on its index
        held = {q: p for q, p in fx.QUERIES.items() if "_4932_" not in q}
        got = [loci.call(self.path, q, *args, rescue_missing=False, index_dir=idx) for q in (fx.QUERIES, held, held)]
        self.assertEqual(got[0], plain)
        self.assertEqual(got[1]["loci"]["2at2759"]["status"], "Single")
        self.assertEqual(got[1], got[2])
        runs = commands(log)
        (index,) = [f for f in os.listdir(idx) if f.endswith(".mpi")]
        index = os.path.join(idx, index)
        self.assertEqual([c for c in runs if "-d" in c], [["-t", "1", "-T1", "-d", "/dev/stdout", self.path]])
        self.assertEqual(runs[0][3], "-d")  # before the first alignment
        self.assertEqual([c[:-1] for c in runs[1:]], [self.aligned(index)] * 3)
        with open(index) as fh:  # whole, and nothing else left beside it but the pieces' number
            self.assertEqual(fh.read(), f"half an index, whole: -t 1 -T1 -d {self.path}\n")
        self.assertEqual(sorted(f.rsplit(".", 1)[1] for f in os.listdir(idx)), ["mpi", "pieces"])
        # a genome's file replaced under its name: not the old one's index
        with gzip.open(self.path, "wt") as fh:
            fh.write(f">c1 another genome\n{self.seq}AC\n")
        loci.call(self.path, held, *args, rescue_missing=False, index_dir=idx)
        self.assertEqual(len([c for c in commands(log) if "-d" in c]), 1)

    def test_the_index_of_each_genetic_code_and_seeding(self):
        mp, log, _mode = indexing_miniprot(self.dir, self.gff)
        idx = os.path.join(self.dir, "idx")
        self.assertEqual(loci.split_options(["-G", "5000", "-k", "5", "-M0", "-n1", "-l4", "-T", "6", "-k7"]),
                         (["-k7", "-M0", "-T6"], ["-G", "5000", "-n1", "-l4"]))
        self.assertEqual(loci.split_options(["-I"]), ([], ["-I"]))
        for table in (1, 4, 1, 4):
            self.assertIsNone(loci.align(self.path, fx.QUERIES, mp, 1, table, self.dir, index_dir=idx)[1])
        runs = commands(log)
        built = [c for c in runs if "-d" in c]
        self.assertEqual(built, [["-t", "1", f"-T{t}", "-d", "/dev/stdout", self.path] for t in (1, 4)])
        indexes = {}  # the table of each index, by its own text
        for name in os.listdir(idx):
            if name.endswith(".mpi"):
                with open(os.path.join(idx, name)) as fh:
                    indexes[fh.read().split()[6]] = os.path.join(idx, name)
        self.assertEqual(sorted(indexes), ["-T1", "-T4"])
        self.assertEqual([c[:-1] for c in runs if "-d" not in c],
                         [self.aligned(indexes[f"-T{t}"], t) for t in (1, 4, 1, 4)])
        # the seeding options shape the index (theirs another one), the intron's the alignment alone
        for intron in ("5000", "9000"):
            loci.align(self.path, fx.QUERIES, mp, 1, 1, self.dir, ["-G", intron, "-k5", "-M0", "-n1"], index_dir=idx)
        runs = commands(log)
        self.assertEqual([c for c in runs if "-d" in c], [["-t", "1", "-k5", "-M0", "-T1", "-d", "/dev/stdout",
                                                           self.path]])
        seeded = [c[-2] for c in runs if "-d" not in c]
        self.assertEqual([c[:-1] for c in runs if "-d" not in c],
                         [self.aligned(seeded[0], 1, ("-G", g, "-n1")) for g in ("5000", "9000")])
        self.assertNotIn(seeded[0], indexes.values())

    def test_a_failed_index_is_no_index(self):
        mp, log, mode = indexing_miniprot(self.dir, self.gff)
        idx = os.path.join(self.dir, "idx")
        args = (self.odb, mp, fx.hmmsearch(self.dir), 1, self.dir)
        for word, why in (("fail", "miniprot failed: killed while indexing"),
                          ("empty", "miniprot failed: its index was not written (empty)")):
            with open(mode, "w") as fh:
                fh.write(word)
            got = loci.call(self.path, fx.QUERIES, *args, rescue_missing=False, index_dir=idx)
            self.assertEqual(got, dict(loci={}, failed=why, windows={}))  # as a failed miniprot run is reported
            self.assertEqual(os.listdir(idx), [])  # the half index gone, none a later call would read
            self.assertEqual([c[3] for c in commands(log)], ["-d"])  # and nothing aligned
        os.remove(mode)
        got = loci.call(self.path, fx.QUERIES, *args, rescue_missing=False, index_dir=idx)
        self.assertEqual(got, loci.call(self.path, fx.QUERIES, *args, rescue_missing=False))
        self.assertEqual(len([c for c in commands(log) if "-d" in c]), 1)

    def test_an_index_of_each_piece_made_once_by_threads(self):
        mp, log, _mode = indexing_miniprot(self.dir, self.gff)
        idx = os.path.join(self.dir, "idx")
        path = os.path.join(self.dir, "three.fna")
        with open(path, "w") as fh:
            fh.writelines(f">c{i}\n{'ACGT' * 250}\n" for i in range(3))

        def run(k):
            work = os.path.join(self.dir, f"work{k}")
            os.makedirs(work)
            return loci.align(path, fx.QUERIES, mp, 1, 1, work, index_dir=idx)

        with mock.patch.object(placement, "CHUNK_BASES", 1500):  # a contig a piece
            with concurrent.futures.ThreadPoolExecutor(4) as pool:
                got = list(pool.map(run, range(4)))
            runs = commands(log)
            built = [c for c in runs if "-d" in c]
            self.assertEqual([os.path.basename(c[-1]) for c in built], [f"piece{i}.fna" for i in range(3)])
            self.assertEqual(len({os.path.dirname(c[-1]) for c in built}), 1)  # one thread's pieces
            pieces = sorted(f for f in os.listdir(idx) if f.endswith(".mpi"))
            self.assertEqual([f.split(".")[2] for f in pieces], ["piece0", "piece1", "piece2"])
            by_call = {}
            for c in runs:
                if "-d" not in c:
                    by_call.setdefault(os.path.dirname(c[-1]), []).append(os.path.basename(c[-2]))
            self.assertEqual(list(by_call.values()), [pieces] * 4)  # each call: every piece's index, in order
            self.assertTrue(all(g == got[0] and g[1] is None for g in got))
            self.assertFalse([f for k in range(4) for f in os.listdir(os.path.join(self.dir, f"work{k}"))
                              if f.startswith("piece")])  # the pieces' FASTA deleted, their indexes read
            # a later call: the genome not cut again, nothing built
            with mock.patch.object(placement, "chunks", side_effect=AssertionError("cut again")):
                self.assertEqual(run(4), got[0])
            self.assertEqual([os.path.basename(c[-2]) for c in commands(log)], pieces)
        # the same genome whole (not over CHUNK_BASES): another index, of the file itself
        self.assertIsNone(run(5)[1])
        self.assertEqual([c[-1] for c in commands(log) if "-d" in c], [path])

    # a contig over CHUNK_BASES in windows (since 2026-11-10; the audit's C8): until then alone in its piece, whole

    def long_contig(self):
        """A genome of a contig of 8,000 bases and one of 1,000 (pieces of 3,000, windows 800 or more into each
        other, cut on multiples of 256 from the contig's ends: 0-2880, 2048-4928, 4096-6976, 6144-8000) with four
        genes of one exon -> (its path, {gene: (contig, start, end)} 1-based with the stop, the contigs): locus 1's
        inside the first window alone, locus 2's across the first window's end (whole in the second alone), locus
        3's inside the second and third windows' overlap (whole in both), locus 1's again on the second contig (the
        stand-in finds a protein once a record)."""
        rng = __import__("random").Random(11)
        c1 = [rng.choice("ACGT") for _ in range(8000)]
        c2 = [rng.choice("ACGT") for _ in range(1000)]
        places = {}
        for gene, locus, seq, at in (("a", "1at2759", c1, 500), ("b", "2at2759", c1, 2700), ("c", "3at2759", c1, 4500),
                                     ("d", "1at2759", c2, 300)):
            cds = fx.cds_of(fx.PROTEINS[locus])
            seq[at:at + len(cds)] = cds
            places[gene] = ("c1" if seq is c1 else "c2", at + 1, at + len(cds))
        path = os.path.join(self.dir, "long.fna")
        with open(path, "w") as fh:  # lines of 70 bases: a window's ends inside lines
            for name, seq in (("c1", c1), ("c2", c2)):
                fh.write(f">{name} a contig\n")
                fh.writelines("".join(seq[i:i + 70]) + "\n" for i in range(0, len(seq), 70))
        return path, places, "".join(c1), "".join(c2)

    def test_a_long_contig_in_windows(self):
        self.assertEqual(placement.windows(8000, 3000, 800), [(0, 2880), (2048, 4928), (4096, 6976), (6144, 8000)])
        self.assertEqual(placement.windows(3000, 3000, 800), [(0, 3000)])   # not over a piece: whole
        self.assertEqual(placement.windows(3001, 3000, 800), [(0, 2745), (1792, 3001)])
        big = placement.windows(10_000_123, 4_000_000, 500_000)  # on miniprot's blocks of 256, from either end
        self.assertEqual((len(big), big[0][0], big[-1][1]), (3, 0, 10_000_123))
        self.assertTrue(all(a % 256 == 0 and (10_000_123 - b) % 256 == 0 and b - a <= 4_000_000 for a, b in big))
        self.assertTrue(all(500_000 <= b - c < 500_256 for (_a, b), (c, _d) in zip(big, big[1:])))
        self.assertEqual(placement.windows(8000, 3000, 0), [(0, 8000)])     # no overlap asked: never cut
        self.assertEqual(placement.window_of("c1" + placement.WINDOW_MARK + "2200"), ("c1", 2200))
        self.assertEqual(placement.window_of("c1|x"), ("c1|x", 0))
        path, _places, c1, c2 = self.long_contig()
        work = os.path.join(self.dir, "pieces")
        os.makedirs(work)
        # as until 2026-11-10 (no overlap given: the six-frame ORFs' pieces): the long contig whole, alone in its piece
        old = placement.chunks(path, work, 3000)
        self.assertEqual([[(n, len(s)) for n, s in odb_mod.read_fasta(p)] for p in old],
                         [[("c1", 8000)], [("c2", 1000)]])
        parts = placement.chunks(path, work, 3000, 800)
        got = [list(odb_mod.read_fasta(p)) for p in parts]
        mark = placement.WINDOW_MARK
        self.assertEqual([[n for n, _s in piece] for piece in got],
                         [[f"c1{mark}0"], [f"c1{mark}2048"], [f"c1{mark}4096"], [f"c1{mark}6144", "c2"]])
        self.assertTrue(all(sum(len(s) for _n, s in piece) <= 3000 for piece in got))  # what the build reserved
        for piece in got:
            for name, seq in piece:
                contig, start = placement.window_of(name)
                self.assertEqual(seq, {"c1": c1, "c2": c2}[contig][start:start + len(seq)])
        self.assertEqual(sum(len(s) for n, s in sum(got, []) if n != "c2"), 8000 + 3 * 832)
        # the genome's bases where the caller knows them (the audit's L5): under a piece, not cut, whatever the file's
        # size says; over it, cut
        self.assertEqual(placement.chunks(path, work, 3000, 800, bases=2500), [path])
        self.assertEqual(len(placement.chunks(path, work, 20_000, 800, bases=30_000)), 1)  # (one piece holds it all)
        self.assertEqual(placement.chunks(path, work, 20_000, 800), [path])

    def test_a_long_contig_s_alignments_those_of_the_whole(self):
        # the stand-in searches each piece (-M0): a gene is found in every window that holds it whole
        path, places, _c1, _c2 = self.long_contig()
        mp = fx.miniprot(self.dir, self.gff)
        queries = fx.QUERIES | {k.replace("1at", "3at"): fx.PROTEINS["3at2759"] for k in fx.QUERIES
                                if k.startswith("1at")}
        folders = [os.path.join(self.dir, f"w{k}") for k in range(2)]
        for folder in folders:
            os.makedirs(folder)
        whole, why = loci.align(path, queries, mp, 1, 1, folders[0], ["-G", "5000", "-M0"])
        self.assertIsNone(why)
        with mock.patch.object(placement, "CHUNK_BASES", 3000), mock.patch.object(placement, "CHUNK_OVERLAP", 800):
            cut, why = loci.align(path, queries, mp, 1, 1, folders[1], ["-G", "5000", "-M0"])
        self.assertIsNone(why)
        self.assertEqual(sorted(f for f in os.listdir(folders[1]) if f.startswith("piece")),
                         [f"piece{k}.fna" for k in range(4)])

        def key(a):
            return a["query"], a["contig"], a["start"]
        self.assertEqual(sorted(cut, key=key), sorted(whole, key=key))  # places, exons, proteins: all of them
        self.assertEqual(sorted({(a["contig"], a["start"], a["end"]) for a in cut}), sorted(places.values()))
        self.assertEqual(len(cut), 2 * len(places))  # two proteins a locus: the overlap's gene once for each
        # of two windows' alignments of a protein that overlap, the best score: one cut short at a window's end is
        # not a second alignment; one beside it, or of the same window, stays
        mark = placement.WINDOW_MARK

        def hit(window, start, end, score, query="q"):
            return dict(query=query, contig=f"c1{mark}{window}", start=start, end=end, strand="+", score=score,
                        exons=[(start, end)])
        got = loci.whole_contigs([hit(0, 2500, 3000, 400.0), hit(2200, 300, 1100, 900.0),    # cut short; whole
                                  hit(2200, 1500, 1900, 300.0), hit(2200, 1800, 2000, 290.0),  # one window's two
                                  hit(4400, 100, 500, 50.0, "r"), dict(hit(0, 5, 9, 1.0), contig="c2")])
        self.assertEqual(sorted((a["query"], a["contig"], a["start"], a["end"], a["exons"], a["score"]) for a in got),
                         [("q", "c1", 2500, 3300, [(2500, 3300)], 900.0),
                          ("q", "c1", 3700, 4100, [(3700, 4100)], 300.0),
                          ("q", "c1", 4000, 4200, [(4000, 4200)], 290.0), ("q", "c2", 5, 9, [(5, 9)], 1.0),
                          ("r", "c1", 4500, 4900, [(4500, 4900)], 50.0)])

    def test_a_long_contig_s_windows_each_their_index(self):
        path, _places, _c1, _c2 = self.long_contig()
        mp, log, _mode = indexing_miniprot(self.dir, self.gff)
        idx = os.path.join(self.dir, "idx")
        with mock.patch.object(placement, "CHUNK_BASES", 3000), mock.patch.object(placement, "CHUNK_OVERLAP", 800):
            for k in range(2):
                work = os.path.join(self.dir, f"work{k}")
                os.makedirs(work)
                self.assertIsNone(loci.align(path, fx.QUERIES, mp, 1, 1, work, index_dir=idx)[1])
                self.assertFalse([f for f in os.listdir(work) if f.startswith("piece")])
        runs = commands(log)
        built = [c for c in runs if "-d" in c]
        self.assertEqual([os.path.basename(c[-1]) for c in built], [f"piece{k}.fna" for k in range(4)])  # once
        pieces = sorted(f for f in os.listdir(idx) if f.endswith(".mpi"))
        self.assertEqual([f.split(".")[2] for f in pieces], [f"piece{k}" for k in range(4)])
        self.assertEqual([os.path.basename(c[-2]) for c in runs if "-d" not in c], pieces * 2)

    def test_the_genome_s_size_where_the_caller_knows_it(self):
        # since 2026-11-10 (the audit's L5): the rescue's gate and the intron bound of the rescue by the genome's true
        # bases where the caller gives them; else by the file's size, as until then
        self.assertEqual([loci.intron_bound(n) for n in (1_000_000, 12_071_326, 64_000_000, 10 ** 11)],
                         [10_000, 12_507, 28_800, 300_000])  # miniprot -I's: 3.6 x the root, 10 kb to 300 kb
        rng = __import__("random").Random(9)
        pad = "".join(rng.choice("ACGT") for _ in range(400))
        path = os.path.join(self.dir, "g3.fna.gz")
        with gzip.open(path, "wt") as fh:
            fh.write(f">c1 a test contig\n{self.seq + pad + fx.cds_of(fx.PROTEINS['3at2759']) + pad}\n")
        mp, log, _mode = indexing_miniprot(self.dir, self.gff)
        queries = fx.QUERIES | {k.replace("1at", "3at"): fx.PROTEINS["3at2759"] for k in fx.QUERIES
                                if k.startswith("1at")}
        args = (path, queries, self.odb, mp, fx.hmmsearch(self.dir), 1, self.dir)

        def rescued(**more):
            got = loci.call(*args, translate=fx.esl_translate(self.dir), **more)
            bound = [c[c.index("-G") + 1] for c in commands(log) if "-G" in c]
            return got["loci"]["3at2759"]["status"], bound
        with mock.patch.object(loci, "RESCUE_MAX_BASES", 100_000_000):
            self.assertEqual(rescued(), ("Single", ["10000"]))  # the file's size / 0.3: a few kb
            self.assertEqual(rescued(bases=64_000_000), ("Single", ["28800"]))
            self.assertEqual(rescued(bases=100_000_001), ("Missing", []))  # over the gate: not rescued
        with mock.patch.object(loci, "RESCUE_MAX_BASES", 1000):  # by the file's size over the gate, by its bases not
            self.assertEqual(rescued(), ("Missing", []))
            self.assertEqual(rescued(bases=900), ("Single", ["10000"]))


if __name__ == "__main__":
    unittest.main()
