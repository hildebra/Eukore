#!/usr/bin/env python3
"""Tests of the data the paper's figures need (since 2026-11-06): a call laid on NCBI's annotation (exons, splice
sites, bases); what became of each NCBI species (the selection funnel); the read benchmark's assignment (hits within
NEAR of the best, their lowest common ancestor, the read's own species held out) and its per-rank verdicts; the test
genomes chosen; the breadth table and the distance to BUSCO's dataset.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import random
import shutil
import struct
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "benchmarks"))  # the benchmarks (since 2026-11-12)
from eukore import annotation, benchmarks, genomes, labels, ncbi as ncbi_mod, sources, twobit  # noqa: E402
from eukore_bench import reads  # noqa: E402
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

    def test_ensembl_and_toga2(self):
        with tempfile.TemporaryDirectory() as tmp:
            # Ensembl's index: the newest build, Ensembl's own on a tie; Prokka's (bacteria) and RefSeq's left out
            builds = {"braker": {"2023_07": "B"}, "ensembl": {"2023_07": "E", "2021_01": "old"},
                      "refseq": {"2026_01": "R"}}
            species = {"species": {"Sp": {"assemblies": {
                "GCA_1.1": {"genebuild_providers": {p: {d: {"paths": {"genebuild": {"files": {"annotations": {
                    "genes.gff3.gz": f"GCA/{x}/genes.gff3.gz"}}}}} for d, x in ds.items()}
                    for p, ds in builds.items()}},
                "GCA_2.1": {"genebuild_providers": {"prokka": {"2024_01": {"paths": {"genebuild": {"files": {
                    "annotations": {"genes.gff3.gz": "GCA/P/genes.gff3.gz"}}}}}}}}}}}}
            path = os.path.join(tmp, "species.json")
            with open(path, "w") as fh:
                json.dump(species, fh)
            self.assertEqual(annotation.ensembl_index(path),
                             {"GCA_1.1": ("ensembl", "2023_07", f"{annotation.ENSEMBL_FTP}/GCA/E/genes.gff3.gz")})
            # TOGA2's: the first reference of TOGA2_REFERENCES listing the assembly; none without an accession
            table = os.path.join(tmp, "toga2.tsv")
            with open(table, "w") as fh:
                fh.write("Directory Name\tSpecies\tSpecies Taxonomy ID\tAssembly name\tNCBI accession\t"
                         "Source of assembly\tscaffold N50 (bp)\n"
                         "Bird__b__GCA_3.1\tAves b\t9\tHLb1\tGCA_3.1\tNCBI\t500\n"
                         "Zoo__z__HL1\tMus z\t7\tHL1\t\tDNA Zoo Consortium\t900\n")
            listings = {"reference_zebrafinch_X": ["Bird__b__GCA_3.1"], "reference_chicken_Y": ["Bird__b__GCA_3.1"],
                        "reference_human_hg38": ["Zoo__z__HL1"]}
            rows = {r["accession"]: r for r in annotation.toga2_rows(table, listings, {"Mammalia": ["Zoo__z__HL1"]})}
            self.assertEqual((rows["TOGA2_HLb1"]["ncbi_accession"], rows["TOGA2_HLb1"]["reference"],
                              rows["TOGA2_HLb1"]["annotation"]), ("GCA_3.1", "reference_chicken_Y",
                             f"{annotation.TOGA2_URL}/TOGA2/reference_chicken_Y/Bird__b__GCA_3.1/query_annotation.gtf.gz"))
            self.assertEqual((rows["TOGA2_HL1"]["ncbi_accession"], rows["TOGA2_HL1"]["genome"],
                              rows["TOGA2_HLb1"]["genome"]),
                             ("-", f"{annotation.TOGA2_URL}/Genomes/Mammalia/Zoo__z__HL1/HL1.2bit", "-"))
            # the index given to the jobs: by GenBank accession, or the RefSeq copy fetched
            os.makedirs(os.path.join(tmp, annotation.INDEX_DIR))
            with open(annotation.index_file(tmp), "w") as fh:
                fh.write("accession\tsource\tbuild\turl\nGCA_3.1\ttoga2\tref\thttps://t\nGCA_3.1\tensembl\te\thttps://e\n")
            jobs = {"GCA_3.1": dict(accession="GCA_3.1", url="https://x/GCA_3.1_a_genomic.fna.gz", annotated=False),
                    "GCA_4.1": dict(accession="GCA_4.1", url="https://x/GCA_4.1_a_genomic.fna.gz")}
            self.assertEqual(annotation.assign(jobs, annotation.read_index(tmp)), 1)
            self.assertEqual(annotation.candidates(jobs["GCA_3.1"]), [("ensembl", "https://e"), ("toga2", "https://t")])
            self.assertEqual(annotation.candidates(jobs["GCA_4.1"]), [("ncbi", "https://x/GCA_4.1_a_genomic.gff.gz")])
            # names: Ensembl's chromosome `2`, TOGA2's accession without its version, a RefSeq copy's sequences
            report = os.path.join(tmp, "report.txt")
            with open(report, "w") as fh:
                fh.write("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
                         "GenBank-Accn\tRelationship\tRefSeq-Accn\n"
                         "chr2\tassembled-molecule\t2\tChromosome\tCM2.1\t=\tNC_2.1\n"
                         "u1\tunplaced-scaffold\tna\tna\tJA1.1\t<>\tna\n")
            alias = annotation.aliases(report)
            self.assertEqual((alias["2"], alias["CM2"], alias["NC_2.1"], alias["CM2.1"], alias["JA1"]),
                             ("CM2.1", "CM2.1", "CM2.1", "NC_2.1", "JA1.1"))
            gtf = os.path.join(tmp, "a.gtf.gz")
            with gzip.open(gtf, "wt") as fh:
                fh.write('CM2\tbed2gtf\tCDS\t101\t250\t.\t-\t2\tgene_id "G_a"; transcript_id "T#ABC#1";\n'
                         '2\tensembl\tCDS\t101\t250\t.\t-\t0\tID=CDS:P1;Parent=transcript:T2;protein_id=P1\n')
            on_genbank = annotation.read_cds(gtf, {"CM2.1": [(1, 1000)]}, alias)
            self.assertEqual(sorted(on_genbank), ["T#ABC#1", "transcript:T2"])
            self.assertEqual((on_genbank["T#ABC#1"]["gene"], on_genbank["T#ABC#1"]["phase"][(101, 250)]), ("G_a", 2))
            on_refseq = annotation.read_cds(gtf, {"NC_2.1": [(1, 1000)]}, alias)  # two steps: 2 -> CM2.1 -> NC_2.1
            self.assertEqual({t["contig"] for t in on_refseq.values()}, {"NC_2.1"})
            # fetched in order: NCBI's 404, then a GFF without coding sequences, then one with them
            genomes_dir = os.path.join(tmp, "genomes")
            os.makedirs(genomes_dir)
            served = {"https://n/GCA_5.1_a_genomic.gff.gz": None, "https://e5": "c1\tx\tregion\t1\t9\t.\t+\t.\tID=c1\n",
                      "https://t5": 'c1\tx\tCDS\t1\t9\t.\t+\t0\ttranscript_id "t";\n'}

            def fetch_file(url, dest, **_how):
                if served.get(url) is None:
                    return "not found"
                with gzip.open(dest, "wt") as fh:
                    fh.write(served[url])
                return ""
            job = dict(accession="GCA_5.1", url="https://n/GCA_5.1_a_genomic.fna.gz",
                       annotations=[("ensembl", "https://e5"), ("toga2", "https://t5")])
            real, annotation.util.fetch_file = annotation.util.fetch_file, fetch_file
            try:
                path, source = annotation.fetch(job, genomes_dir)
            finally:
                annotation.util.fetch_file = real
            self.assertEqual((path, source), (annotation.gff_file(genomes_dir, "GCA_5.1", "toga2"), "toga2"))
            self.assertFalse(os.path.isfile(annotation.gff_file(genomes_dir, "GCA_5.1", "ensembl")))
            self.assertEqual(annotation.found(genomes_dir, "GCA_5.1"), (path, "toga2"))

    def test_veupathdb(self):
        with tempfile.TemporaryDirectory() as tmp:
            # its organism table: the GFF's URL moved to the file service; a genome without an accession left out
            table = os.path.join(tmp, "organisms.json")
            with open(table, "w") as fh:
                json.dump({"records": [
                    {"attributes": {"genome_version": "GCA_6.1", "annotation_source": "VEuPathDB",
                                    "annotation_version": "2024",
                                    "URLgff": "http://PlasmoDB.org/common/downloads/Current_Release/P/gff/data/P.gff"}},
                    {"attributes": {"genome_version": "2014-05-13", "URLgff": "http://x/common/downloads/y.gff"}}]}, fh)
            self.assertEqual(annotation.veupathdb_index(table), {"GCA_6.1": (
                "VEuPathDB 2024", "https://PlasmoDB.org/a/service/raw-files/Current_Release/P/gff/data/P.gff")})
            # without an API key its annotations are not given to the jobs
            v = "https://plasmodb.org/a/service/raw-files/P.gff"
            index = {"GCA_6.1": [("veupathdb", v), ("ensembl", "https://e")]}
            jobs = {"GCA_6.1": dict(accession="GCA_6.1", url="https://n/GCA_6.1_a_genomic.fna.gz", annotated=False)}
            annotation.assign(jobs, index)
            self.assertEqual(annotation.candidates(jobs["GCA_6.1"]), [("ensembl", "https://e")])
            annotation.assign(jobs, index, veupathdb=True)
            self.assertEqual(annotation.candidates(jobs["GCA_6.1"])[0], ("veupathdb", v))
            # with it: asked for with the key as a bearer token, kept plain; its own sequence names matched by length
            asked = {}

            def fetch_file(url, dest, **how):
                asked[url] = how
                if url != v:
                    return "not found"
                with open(dest, "w") as fh:
                    fh.write("##gff-version 3\n##sequence-region Pf_01 1 5000\n##sequence-region Pf_02 1 7000\n"
                             "Pf_01\tVEuPathDB\tCDS\t101\t400\t.\t+\t0\tID=c1;Parent=PF01_t1\n")
                return ""
            genomes_dir = os.path.join(tmp, "genomes")
            os.makedirs(genomes_dir)
            real, annotation.util.fetch_file = annotation.util.fetch_file, fetch_file
            try:
                path, source = annotation.fetch(jobs["GCA_6.1"], genomes_dir, "KEY")
            finally:
                annotation.util.fetch_file = real
            self.assertEqual((path, source), (annotation.gff_file(genomes_dir, "GCA_6.1", "veupathdb"), "veupathdb"))
            self.assertTrue(path.endswith(".veupathdb.gff"))
            self.assertEqual((asked[v]["headers"], asked[v]["gzip_check"]), ({"Authorization": "Bearer KEY"}, False))
            self.assertNotIn("headers", asked["https://n/GCA_6.1_a_assembly_report.txt"])  # the key to VEuPathDB alone
            # a URL its table names on another host, or over http: never asked with the key (since 2026-11-08)
            for url in ("https://cdn.example.com/P.gff", "http://plasmodb.org/P.gff"):
                asked.clear()
                elsewhere = dict(jobs["GCA_6.1"], accession="GCA_7.1", annotations=[("veupathdb", url)])
                real, annotation.util.fetch_file = annotation.util.fetch_file, fetch_file
                try:
                    annotation.fetch(elsewhere, genomes_dir, "KEY")
                finally:
                    annotation.util.fetch_file = real
                self.assertNotIn(url, asked)
            report = os.path.join(tmp, "report.txt")
            with open(report, "w") as fh:
                fh.write("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
                         "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\n"
                         "1\tassembled-molecule\t1\tChromosome\tLR1.1\t<>\tna\tPrimary Assembly\t5000\n"
                         "2\tassembled-molecule\t2\tChromosome\tLR2.1\t<>\tna\tPrimary Assembly\t7000\n"
                         "3\tassembled-molecule\t3\tChromosome\tLR3.1\t<>\tna\tPrimary Assembly\t7000\n")
            alias = annotation.aliases(report, path)
            self.assertEqual(alias["Pf_01"], "LR1.1")
            self.assertNotIn("Pf_02", alias)  # two sequences of its length: not told apart
            self.assertEqual(sorted(annotation.read_cds(path, {"LR1.1": [(1, 5000)]}, alias)), ["PF01_t1"])


class AnnotationSwapTest(unittest.TestCase):
    def test_a_model_against_an_independent_truth(self):
        # tools/annotation_swap.py (since 2026-11-09): miniprot's model and the copy shipped, each laid on the truth
        sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tools"))
        import annotation_swap  # noqa: PLC0415
        truth = {"t1": dict(contig="c1", strand="+", gene="g1", cds=[(100, 200), (300, 400)], phase={})}
        exact = annotation_swap.measure(dict(contig="c1", strand="+", status="Single",
                                             exons=[(100, 200), (300, 400)]), truth)
        off = annotation_swap.measure(dict(contig="c1", strand="+", status="Single", exons=[(100, 200), (310, 400)]),
                                      truth)
        self.assertEqual(exact[:4], ("t1", 1, 1, True))
        self.assertEqual(off[:4], ("t1", 1, 0, False))
        self.assertLess(off[4], exact[4])
        self.assertEqual(annotation_swap.verdict(off, exact), "better")
        self.assertEqual(annotation_swap.verdict(exact, off), "worse")
        self.assertEqual(annotation_swap.verdict(exact, exact), "same")
        self.assertEqual(annotation_swap.measure(dict(contig="c2", strand="+", status="Single", exons=[(1, 9)]),
                                                 truth)[0], "-")


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
        # since 2026-11-12: why a taxon's assemblies were left out, each taxon with NCBI's assemblies of it
        self.assertEqual(fates.reasons[6253], (1, {"contaminated": 1}))
        self.assertEqual(fates.reasons[4932][1], {})
        rows = {r["taxid"]: r for r in benchmarks.taxa_rows(fates, ncbi)}
        self.assertEqual((rows[6253]["fate"], rows[6253]["assemblies"], rows[6253]["reasons"]),
                         ("every assembly excluded", 1, "contaminated"))
        self.assertEqual(rows[4932]["reasons"], "-")

    def test_toga2_genomes_of_species_ncbi_has_none_of(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            job_map = genomes.jobs(ncbi, gb)
            row = dict(ncbi_accession="-", organism="Ascaris suum", scaffold_n50="10", reference="r",
                       annotation="https://t/a.gtf.gz", genome="https://t/A1.2bit", assembly_source="DNA Zoo")
            rows = [dict(row, accession="TOGA2_A1", species_taxid="6253"),
                    dict(row, accession="TOGA2_A2", species_taxid="6253", scaffold_n50="99", genome="https://t/A2.2bit"),
                    dict(row, accession="TOGA2_A3", species_taxid="6253", scaffold_n50="999", genome="-"),
                    dict(row, accession="TOGA2_Y", species_taxid="4932")]  # NCBI has yeast's genome
            got = genomes.toga2_jobs(ncbi, rows, job_map)
            self.assertEqual(list(got), ["TOGA2_A2"])  # the longest scaffold N50 with a genome to fetch
            job = got["TOGA2_A2"]
            self.assertEqual((job["why"], job["url"], job["annotations"], job["annotated"]),
                             (["toga2"], "https://t/A2.2bit", [("toga2", "https://t/a.gtf.gz")], False))
            self.assertEqual(genomes.toga2_jobs(ncbi, rows, job_map, groups=("fungi",)), {})
            self.assertEqual(genomes.fates(ncbi, gb, dict(job_map, **got))[6253], "TOGA2's genome")
            # in place of NCBI's genome where its contig N50 is twice NCBI's or more and its length 90% of it or more
            par = next(a for a, j in job_map.items() if j["species_taxid"] == 27291)
            better = [dict(row, accession="TOGA2_P", species_taxid="27291", contig_n50="2000")]
            ncbi_stats = dict(contig_n50=1000, total_length=10_000)
            swaps = genomes.prefer_toga2(ncbi, better, job_map, lambda job: ncbi_stats, lambda r: 9_000)
            self.assertEqual(list(swaps), [par])
            self.assertEqual((swaps[par]["accession"], swaps[par]["why"], swaps[par]["others"][0]),
                             ("TOGA2_P", ["toga2"], par))
            for theirs, length in ((dict(ncbi_stats, contig_n50=1001), 9_000),  # under twice the contig N50
                                   (ncbi_stats, 8_999),  # contiguous but short: a part of the genome
                                   (None, 9_000), (ncbi_stats, None)):  # either not known: NCBI's kept
                self.assertEqual(genomes.prefer_toga2(ncbi, better, job_map, lambda job: theirs, lambda r: length),
                                 {})
            stats = os.path.join(tmp, "stats.txt")
            with open(stats, "w") as fh:
                fh.write("# Assembly name: x\nall\tall\tall\tall\ttotal-length\t9000\n"
                         "all\tall\tall\tall\tscaffold-N50\t900\nall\tall\tall\tall\tcontig-N50\t450\n")
            self.assertEqual(genomes.assembly_stats(stats), dict(contig_n50=450, total_length=9000))
            self.assertEqual(genomes.stats_url(dict(url="https://x/GCA_1.1_a_genomic.fna.gz")),
                             "https://x/GCA_1.1_a_assembly_stats.txt")
            sizes = os.path.join(tmp, "chrom.sizes")
            with open(sizes, "w") as fh:
                fh.write("c1\t700\nc2\t300\n")
            self.assertEqual(genomes.total_length(sizes), 1000)
            self.assertEqual(genomes.chrom_sizes_url(dict(genome="https://t/G/A/x.2bit")), "https://t/G/A/chrom.sizes")

    def test_two_bit_genomes(self):
        # a .2bit of two sequences (version 0, little-endian): N blocks and soft-masked blocks
        seqs = [("chr1", "ACGTTGCAAC", [(2, 3)], [(6, 4)]), ("u2", "GGGTA", [], [])]
        code = {"T": 0, "C": 1, "A": 2, "G": 3}
        records = []
        for name, s, ns, masks in seqs:
            packed = bytearray()
            padded = s + "T" * (-len(s) % 4)
            for i in range(0, len(padded), 4):
                b = 0
                for c in padded[i:i + 4]:
                    b = b << 2 | code[c]
                packed.append(b)
            rec = struct.pack("<II", len(s), len(ns)) + b"".join(struct.pack("<I", x[0]) for x in ns)
            rec += b"".join(struct.pack("<I", x[1]) for x in ns) + struct.pack("<I", len(masks))
            rec += b"".join(struct.pack("<I", x[0]) for x in masks) + b"".join(struct.pack("<I", x[1]) for x in masks)
            records.append(rec + struct.pack("<I", 0) + bytes(packed))
        head = struct.pack("<IIII", twobit.SIGNATURE, 0, len(seqs), 0)
        index_size = sum(1 + len(n) + 4 for n, *_x in seqs)
        offset, index = len(head) + index_size, b""
        for (name, *_x), rec in zip(seqs, records):
            index += bytes([len(name)]) + name.encode() + struct.pack("<I", offset)
            offset += len(rec)
        with tempfile.TemporaryDirectory() as tmp:
            path, fasta = os.path.join(tmp, "g.2bit"), os.path.join(tmp, "g.fna.gz")
            with open(path, "wb") as fh:
                fh.write(head + index + b"".join(records))
            self.assertEqual(list(twobit.records(path)), [("chr1", b"ACNNNGcaac"), ("u2", b"GGGTA")])
            self.assertEqual(twobit.to_fasta(path, fasta), 15)
            with gzip.open(fasta, "rt") as fh:
                self.assertEqual(fh.read(), ">chr1\nACNNNGcaac\n>u2\nGGGTA\n")
            # TOGA2's genome fetched: its .2bit in --local_tmp (since 2026-11-08), none left anywhere after
            from unittest import mock  # noqa: PLC0415
            out, scratch = os.path.join(tmp, "genomes"), os.path.join(tmp, "scratch")
            os.makedirs(out)
            os.makedirs(scratch)
            got = []

            def fetch(url, dest, **_how):
                got.append(dest)
                shutil.copyfile(path, dest)
                return ""
            with mock.patch.object(genomes.util, "fetch_file", fetch):
                job = dict(url="https://t/A1.2bit")
                self.assertEqual(genomes.fetch_genome(job, os.path.join(out, "TOGA2_A1.fna.gz"), scratch=scratch), "")
            # since 2026-11-10 under this process' own name (C9), still a *.2bit for the fetch stage's leftover glob
            self.assertEqual(got, [genomes.util.fetch_part(os.path.join(scratch, "TOGA2_A1.fna.gz"))[:-5] + ".2bit"])
            self.assertEqual((os.listdir(out), os.listdir(scratch)), (["TOGA2_A1.fna.gz"], []))


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


class CompareReferencesTest(unittest.TestCase):
    def test_metaphlan_and_pr2_labelled_by_ncbi(self):
        from eukore_bench import references as cr  # noqa: PLC0415
        info = ["27987__gene-1\t{'clade': 't__EUK27987', 'ext': [], 'len': 259, 'taxon': 'k__Eukaryota|t__EUK27987'}",
                "UniRef90_X|1__2|SGB1\t{'clade': 't__SGB1', 'ext': [], 'len': 900, 'taxon': 'k__Bacteria|t__SGB1'}"]
        self.assertEqual(cr.euk_markers(info), {"27987__gene-1": 27987})
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            index = cr.rrna.name_index(ncbi)
            # PR2 5.1.1's taxo_long titles: record|gene|organelle||its 9 ranks
            lineage = "Eukaryota|Obazoa|Opisthokonta|Fungi|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|"
            got = [cr.pr2_taxid(f"R{i}|18S_rRNA|{where}||{lineage}{x}", index, ncbi) for i, (where, x) in enumerate((
                ("nucleus", "Saccharomyces|Saccharomyces_cerevisiae"),        # its species
                ("nucleus", "Saccharomyces|Saccharomyces_novus"),             # NCBI has no such species: its genus
                ("nucleus", "Saccharomyces|Saccharomyces_sp."),               # no species named: its genus
                ("nucleus", "Novus|Novus_ignotus"),                           # neither
                ("plastid", "Saccharomyces|Saccharomyces_cerevisiae")))]      # an organelle's
        self.assertEqual(got, [("R0", 4932), ("R1", 4930), ("R2", 4930), None, None])

    def test_detection_f1_and_misses(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        cache = {}
        # a reference of S. cerevisiae (a strain's record) and K. lactis
        held = reads.holdings({"r1": 559292, "r2": 28985}, ncbi, cache)
        self.assertTrue(reads.in_reference(4932, "in", held, ncbi, cache))
        self.assertFalse(reads.in_reference(27291, "in", held, ncbi, cache))      # S. paradoxus: not in it
        self.assertTrue(reads.in_reference(27291, "held_out", held, ncbi, cache))  # its genus is
        self.assertFalse(reads.in_reference(4932, "held_out", held, ncbi, cache))  # no other Saccharomyces
        # a sample of S. cerevisiae: 4 reads to it, 3 to K. lactis (a wrong species, detected), 2 of a bacterium's to
        # S. paradoxus (under MIN_READS); depth 10: the test's reads 0..9, 5 of the one background genome's
        placed = {"T": {0: 4932, 1: 559292, 2: 4932, 3: 4932, 4: 28985, 5: 28985, 6: 28985, 12: 4932},
                  "B": {0: 27291, 1: 27291, 7: 28985}}
        self.assertEqual(reads.detect(placed, "T", 4932, "species", 10, ["B"], ncbi, cache), (4, True, 1))
        self.assertEqual(reads.detect(placed, "T", 4932, "species", 3, ["B"], ncbi, cache), (3, True, 0))
        # a deep sample: 1% of the reads placed is the floor (3 reads of 1,000 placed are under it)
        deep = {"T": {i: 4932 for i in range(997)} | {997: 28985, 998: 28985, 999: 28985}}
        self.assertEqual(reads.detect(deep, "T", 4932, "species", 1000, [], ncbi, cache), (997, True, 0))
        deep["T"].update({i: 28985 for i in range(1000, 1012)})  # 15 K. lactis reads of 1,012: detected
        self.assertEqual(reads.detect(deep, "T", 4932, "species", 2000, [], ncbi, cache), (997, True, 1))
        rows = [dict(reference="R", mode="in", depth=10, congeners=1, detected="t", in_reference="t", false_taxa=1),
                dict(reference="R", mode="in", depth=10, congeners=0, detected="f", in_reference="t", false_taxa=0),
                dict(reference="R", mode="in", depth=10, congeners=1, detected="f", in_reference="f", false_taxa=0),
                dict(reference="R", mode="held_out", depth=10, congeners=0, detected="f", in_reference="f",
                     false_taxa=0)]  # no congeners: left out
        got = reads.detection_summary(rows)["R"]
        self.assertEqual({k: got["in_10"][k] for k in ("tp", "fn_missed", "fn_absent", "fp", "f1")},
                         dict(tp=1, fn_missed=1, fn_absent=1, fp=1, f1=0.4))
        self.assertNotIn("held_out_10", got)

    def test_a_harder_benchmark(self):
        # since 2026-11-09: a test's whole order out of the reference (novel), and a eukaryote background (the host)
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            table = os.path.join(tmp, "background.tsv")
            with open(table, "w") as fh:
                fh.write("# a comment\naccession\tassembly\torganism\ttaxid\nGCF_000005845.2\tASM584v2\tE. coli\t\n"
                         "NC_000021.9\tGRCh38.p14\tHomo sapiens chromosome 21\t9606\n")
            self.assertEqual(reads.background_rows(table), [
                ("GCF_000005845.2", "ASM584v2", "E. coli", None),
                ("NC_000021.9", "GRCh38.p14", "Homo sapiens chromosome 21", 9606)])
            self.assertTrue(reads.background_url("NC_000021.9", "x").endswith("id=NC_000021.9&rettype=fasta&retmode=text"))
            open(os.path.join(tmp, "NC_000021.9.fna.gz"), "w").close()
            self.assertEqual(reads.read_background(table, tmp), [("Homo_sapiens_chromosome_21",
                                                                   os.path.join(tmp, "NC_000021.9.fna.gz"), 9606)])
        cache = {}
        taxa = {"y": 4932, "p": 27291, "k": 28985, "c": 5476, "h": 9606}
        hits = [(0.99, "y"), (0.985, "p"), (0.95, "c")]
        self.assertEqual(reads.assign(hits, taxa, ncbi), 4930)  # yeasts: their genus
        # Saccharomycetales out: the read goes nowhere it is right below its class; Candida, its own order's, is out too
        self.assertIsNone(reads.assign(hits, taxa, ncbi, clade=("order", 4892), cache=cache))
        self.assertEqual(reads.assign(hits + [(0.97, "h")], taxa, ncbi, clade=("order", 4892), cache=cache), 9606)
        held = reads.holdings({"r1": 559292, "r2": 28985, "r3": 5476}, ncbi, cache)
        self.assertFalse(reads.in_reference(4932, "novel", held, ncbi, cache))  # one order in Saccharomycetes alone
        # the host's reads on its own lineage: no false detection; on another species: one
        placed = {"T": {i: 4932 for i in range(6)}, "H": {0: 9606, 1: 9606, 2: 9606, 3: 9606}}
        self.assertEqual(reads.detect(placed, "T", 4932, "species", 8, ["H"], ncbi, cache, {9606}), (6, True, 0))
        self.assertEqual(reads.detect(placed, "T", 4932, "species", 8, ["H"], ncbi, cache), (6, True, 1))
        self.assertEqual(reads.MODES, ("in", "held_out", "novel"))
        self.assertEqual(reads.score(9606, 9606, ncbi, cache)["species"], "correct")  # the host's read on itself

    def test_several_compared_references(self):
        self.assertEqual(reads.compared(None), [])
        self.assertEqual(reads.compared(("A", "a.fa", "a.tsv")), [("A", "a.fa", "a.tsv")])
        self.assertEqual(reads.compared([("A", "a.fa", "a.tsv"), ("B", "b.fa", "b.tsv")]),
                         [("A", "a.fa", "a.tsv"), ("B", "b.fa", "b.tsv")])


if __name__ == "__main__":
    unittest.main()
