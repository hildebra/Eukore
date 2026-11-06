#!/usr/bin/env python3
"""Tests of the data the paper's figures need (since 2026-11-06): a call laid on NCBI's annotation (exons, splice
sites, bases); what became of each NCBI species (the selection funnel); the read benchmark's assignment (hits within
NEAR of the best, their lowest common ancestor, the read's own species held out) and its per-rank verdicts; the test
genomes chosen; the breadth table and the distance to BUSCO's dataset.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import annotation, benchmarks, genomes, labels, ncbi as ncbi_mod, reads, sources  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST, PARADOXUS, LACTIS = "GCA_000146045.2", "GCA_000000013.1", "GCA_000002515.1"


class AnnotationTest(unittest.TestCase):
    def test_a_call_on_ncbi_gene_models(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = os.path.join(tmp, "done")
            os.makedirs(folder)
            # a copy on - at 1001-2000 (contig c1), three exons on the copy (0-based, end excluded, its own strand)
            with open(os.path.join(folder, "loci.tsv"), "w") as fh:
                fh.write("locus\tstatus\tcontig\tcontig_length\tstart\tend\tstrand\n"
                         "1at2759\tSingle\tc1\t5000\t1001\t2000\t-\n2at2759\tMissing\t-\t-\t-\t-\t-\n"
                         "3at2759\tSingle\tc2\t5000\t101\t400\t+\n")
            with open(os.path.join(folder, "features.tsv"), "w") as fh:
                fh.write("locus\tkind\tstart\tend\n1at2759\texon\t0\t200\n1at2759\tintron\t200\t300\n"
                         "1at2759\texon\t300\t700\n1at2759\tintron\t700\t800\n1at2759\texon\t800\t1000\n"
                         "3at2759\texon\t0\t300\n")
            calls = annotation.calls_of(folder)
            self.assertEqual(calls["1at2759"]["exons"], [(1001, 1200), (1301, 1700), (1801, 2000)])
            self.assertNotIn("2at2759", calls)
            # NCBI's transcript: the same introns, its last exon 3 bases longer (the stop); another transcript off
            gff = os.path.join(tmp, "a.gff.gz")
            with gzip.open(gff, "wt") as fh:
                fh.write("##gff-version 3\n")
                for a, b in ((998, 1200), (1301, 1700), (1801, 2000)):
                    fh.write(f"c1\tGnomon\tCDS\t{a}\t{b}\t.\t-\t0\tID=cds1;Parent=rna-XM_1;gene=ABC1\n")
                fh.write("c1\tGnomon\tCDS\t4000\t4100\t.\t-\t0\tID=cds2;Parent=rna-XM_2;gene=FAR\n")
                fh.write("c1\tGnomon\tCDS\t1100\t1150\t.\t+\t0\tID=cds3;Parent=rna-XM_3;gene=OTHER%2C1\n")
            self.assertEqual(annotation.annotate(folder, gff), 2)
            rows = {r["locus"]: r for r in fx.read_tsv(os.path.join(folder, annotation.TABLE))}
            r = rows["1at2759"]
            self.assertEqual((r["transcript"], r["gene"], r["exons"], r["ann_exons"], r["exact_exons"], r["introns"],
                              r["exact_introns"], r["bases"], r["ann_bases"], r["shared_bases"]),
                             ("rna-XM_1", "ABC1", "3", "3", "2", "2", "2", "800", "803", "800"))
            self.assertEqual(rows["3at2759"]["transcript"], "-")  # no annotated gene there
            # the GFF names the RefSeq copy's sequence (NC_1), the calls its GenBank name (c1): the report joins them
            renamed, report = os.path.join(tmp, "b.gff.gz"), os.path.join(tmp, "report.txt")
            with gzip.open(gff, "rt") as src, gzip.open(renamed, "wt") as fh:
                fh.write(src.read().replace("c1\t", "NC_1\t"))
            with open(report, "w") as fh:
                fh.write("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
                         "GenBank-Accn\tRelationship\tRefSeq-Accn\n1\tassembled-molecule\t1\tChromosome\tc1\t=\tNC_1\n")
            self.assertEqual(annotation.annotate(folder, renamed), 2)
            self.assertEqual(fx.read_tsv(os.path.join(folder, annotation.TABLE))[0]["transcript"], "-")
            annotation.annotate(folder, renamed, report)
            self.assertEqual(fx.read_tsv(os.path.join(folder, annotation.TABLE))[0]["shared_bases"], "800")
            self.assertEqual(annotation.gff_url(dict(fetch="GCA_1", url="https://x/GCA_1_a_genomic.fna.gz")),
                             "https://x/GCA_1_a_genomic.gff.gz")  # not known: asked for by default
            self.assertEqual(annotation.gff_url(dict(fetch="GCA_1", url="https://x/GCA_1_a_genomic.fna.gz",
                                                     annotated=True)), "https://x/GCA_1_a_genomic.gff.gz")
            self.assertIsNone(annotation.gff_url(dict(fetch="GCA_1", url=None)))  # records without an assembly
            # a GFF without coding sequences (an unannotated assembly's regions alone): no annotation, no rows
            bare = os.path.join(tmp, "bare.gff.gz")
            with gzip.open(bare, "wt") as fh:
                fh.write("##gff-version 3\nc1\tGenbank\tregion\t1\t5000\t.\t+\t.\tID=c1\n")
            self.assertEqual(annotation.annotate(folder, bare), 0)
            self.assertIsNone(annotation.gff_url(dict(fetch="GCF_1", url="https://x/GCF_1_a_genomic.fna.gz",
                                                      annotated=False)))
            self.assertEqual(annotation.gff_url(dict(fetch="GCF_1", url="https://x/GCF_1_a_genomic.fna.gz")),
                             "https://x/GCF_1_a_genomic.gff.gz")
            # no annotation at NCBI: a table without rows, so that it is not asked for again
            self.assertEqual(annotation.annotate(folder, None), 0)
            self.assertEqual(fx.read_tsv(os.path.join(folder, annotation.TABLE)), [])


class FatesTest(unittest.TestCase):
    def test_what_became_of_each_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK
                                                            + [fx.assembly("GCA_000009999.1", 5475, 5475, "Candida")]
                                                            + [fx.assembly("GCA_000009998.1", 6253, 6253, "Ascaris",
                                                                           flags="contaminated")])])
            fates = genomes.fates(ncbi, gb, genomes.jobs(ncbi, gb))
        self.assertEqual((fates[4932], fates[953], fates[1080349], fates[5475], fates[6253]),
                         ("chosen", "placeholder", "hybrid", "not a species", "every assembly excluded"))
        self.assertIn(fates[27291], ("chosen", "chosen, flagged"))
        self.assertTrue(set(fates.values()) <= set(genomes.FATES))


class ReadsTest(unittest.TestCase):
    def test_assignment_and_its_verdicts(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            taxa = {"cer": 4932, "par": 27291, "lac": 28985, "asc": 6253}
            hits = [(0.99, "cer"), (0.985, "par"), (0.95, "lac")]
            self.assertEqual(reads.assign(hits, taxa, ncbi), 4930)  # S. cerevisiae and S. paradoxus: their genus
            self.assertEqual(reads.assign(hits[:1] + hits[2:], taxa, ncbi), 4932)
            self.assertEqual(reads.assign(hits, taxa, ncbi, leave_out=4932), 27291)  # its own species held out
            self.assertIsNone(reads.assign([(0.85, "asc")], taxa, ncbi))  # under MIN_ID
            self.assertIsNone(reads.assign([(0.99, "unknown")], taxa, ncbi))
            cache = {}
            v = reads.score(4932, 27291, ncbi, cache)  # a read of S. cerevisiae given to S. paradoxus
            self.assertEqual((v["species"], v["genus"], v["family"], v["phylum"]),
                             ("wrong", "correct", "correct", "correct"))
            v = reads.score(4932, 4930, ncbi, cache)  # given to its genus: right, and above the species
            self.assertEqual((v["species"], v["genus"]), ("above", "correct"))
            v = reads.score(4932, 6253, ncbi, cache)  # given to a nematode
            self.assertEqual((v["genus"], v["phylum"]), ("wrong", "wrong"))
            # minimap2's PAF, grouped per read; a hit covering under MIN_COVER of the read left out
            paf = os.path.join(tmp, "h.paf")
            with open(paf, "w") as fh:
                fh.write("r1\t150\t0\t150\t+\tcer\t900\t0\t150\t148\t150\t60\n"
                         "r1\t150\t0\t60\t+\tlac\t900\t0\t60\t60\t60\t0\n"
                         "r2\t150\t0\t150\t-\tpar\t900\t0\t150\t140\t150\t0\n")
            got = dict(reads.paf_groups(paf))
            self.assertEqual(list(got), ["r1", "r2"])
            self.assertEqual([t for _i, t in got["r1"]], ["cer"])
            # a read: READ_LEN bases, substitutions at ERROR, either strand
            rng = random.Random(3)
            contig = "".join(rng.choice("ACGT") for _ in range(1000))
            read = reads.read_of(contig, 100, rng)
            self.assertEqual(len(read), reads.READ_LEN)
            self.assertTrue(read in contig or reads.revcomp(read) in contig
                            or sum(a != b for a, b in zip(read, contig[100:250])) <= 5
                            or sum(a != b for a, b in zip(reads.revcomp(read), contig[100:250])) <= 5)
            self.assertIsNone(reads.read_of(contig, 900, rng))  # off the end

    def test_test_genomes_chosen(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            reps = {YEAST: (4932, "y"), PARADOXUS: (27291, "p"), LACTIS: (28985, "k"), "GCA_9": (6253, "a")}
            got = reads.choose(ncbi, reps, 1, exists=lambda p: p != "a")
        # one fungus: a Saccharomyces (two of its genus; the first by accession), the nematode's genome gone
        self.assertEqual(got, [PARADOXUS])


class TablesTest(unittest.TestCase):
    def test_breadth_and_distance_to_the_dataset(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            species = {"4932": dict(representative=YEAST, supergroup="Fungi"),
                       "28985": dict(representative=LACTIS, supergroup="Fungi"),
                       "6253": dict(representative="-", supergroup="Metazoa")}
            fates = {4932: "chosen", 28985: "chosen, flagged", 27291: "every assembly excluded", 953: "placeholder"}
            rows = {r["supergroup"]: r for r in benchmarks.breadth_rows(species, fates, {4932, 27291}, ncbi)}
            self.assertEqual((rows["Fungi"]["ncbi_species"], rows["Fungi"]["eukore"], rows["Fungi"]["compared"],
                              rows["Fungi"]["both"], rows["Fungi"]["eukore_only"], rows["Fungi"]["compared_only"]),
                             (3, 2, 2, 1, 1, 1))
            self.assertEqual((rows["Metazoa"]["eukore_rrna_only"], rows["all"]["ncbi_species"]), (1, 3))

            class Odb:  # BUSCO's dataset of S. cerevisiae (a strain's taxid) and the nematode
                def species(self):
                    return [559292, 6253]
            got = benchmarks.odb_distance(Odb(), ncbi, {4932, 27291, 28985, 6253, 9606})
        self.assertEqual(got[27291], ("genus", "f"))      # S. paradoxus: S. cerevisiae is of its genus
        self.assertEqual(got[28985], ("family", "f"))     # K. lactis: of its family
        self.assertEqual(got[4932], ("none", "t"))        # S. cerevisiae itself: no other fungus in the dataset
        self.assertEqual(got[9606], ("none", "f"))        # human: the nematode shares no phylum (Metazoa a kingdom)

    def test_other_databases_and_rrna_copies(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            # PR2 names two NCBI species (one by its own name, one as written in PR2), a placeholder and an unknown one
            pr2 = labels.Reference("PR2", "x", labels.PR2_RANKS)
            for sp in ("Saccharomyces_cerevisiae", "Saccharomyces_paradoxus", "Saccharomyces_sp.", "Novus_ignotus"):
                pr2.add(["Eukaryota", "Obazoa", "Opisthokonta", "Fungi", "x", "x", "x", sp.split("_")[0], sp])
            got = benchmarks.reference_species([pr2], ncbi)
            self.assertEqual(got, {"PR2": {4932, 27291}})
            self.assertEqual(benchmarks.reference_species([pr2], ncbi, inside=lambda t: t == 4932), {"PR2": {4932}})
            species = {"4932": dict(supergroup="Fungi"), "28985": dict(supergroup="Fungi"),
                       "6253": dict(supergroup="Metazoa")}
            rows = {(r["database"], r["supergroup"]): r for r in
                    benchmarks.overlap_rows(species, {"PR2": {4932, 27291}, "ROD": {6253}}, ncbi)}
            self.assertEqual({k: rows[("PR2", "Fungi")][k] for k in ("species", "both", "eukore_only",
                                                                     "database_only")},
                             dict(species=2, both=1, eukore_only=1, database_only=1))
            self.assertEqual((rows[("PR2", "all")]["eukore_only"], rows[("ROD", "Metazoa")]["both"]), (2, 1))
        genomes_ = {YEAST: dict(species_taxid="4932", level="Chromosome", contigs="17", rrna_18s="2", rdna_units="1",
                                rdna_complete="1", source="model"),
                    "EUK_6253": dict(species_taxid="6253", source="no_assembly")}
        rows = benchmarks.copy_rows(genomes_, species)
        self.assertEqual([(r["accession"], r["supergroup"], r["rrna_18s"], r["rdna_units"]) for r in rows],
                         [(YEAST, "Fungi", "2", "1")])

    def test_organelle_gene_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "rrna"))
            os.makedirs(os.path.join(tmp, "organelle"))
            with open(os.path.join(tmp, "rrna", "records.tsv"), "w") as fh:
                fh.write("record\tspecies_taxid\tsource\tregions\nr1\t4932\tgenome\t18S:0-1800,ITS1:1800-2100\n")
            with open(os.path.join(tmp, "organelle", "genes.tsv"), "w") as fh:
                fh.write("record\tgene\torganelle\tspecies_taxid\tlength\tpartial\n"
                         "m1_SSU\tSSU\tmitochondrion\t4932\t1650\tf\nm1_COX1\tCOX1\tmitochondrion\t4932\t12880\tf\n"
                         "m2_SSU\tSSU\tmitochondrion\t4932\t300\tt\n")
            rows = {(r["source"], r["region"]): r for r in
                    benchmarks.region_rows(tmp, {"4932": dict(supergroup="Fungi")})}
        self.assertEqual(rows[("genome", "18S")]["median"], 1800)
        self.assertEqual((rows[("mitochondrion", "SSU")]["records"], rows[("mitochondrion", "SSU")]["median"]),
                         (1, 1650))  # the partial one left out
        self.assertEqual(rows[("mitochondrion", "COX1")]["median"], 12880)


if __name__ == "__main__":
    unittest.main()
