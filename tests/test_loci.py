#!/usr/bin/env python3
"""Tests of Eukore's own locus calls (odb.py, loci.py; since 2026-11-05): the dataset read (cutoffs, the HMMs' own
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
from eukore import loci, ncbi as ncbi_mod, odb as odb_mod, qc  # noqa: E402
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


if __name__ == "__main__":
    unittest.main()
