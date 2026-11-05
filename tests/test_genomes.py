#!/usr/bin/env python3
"""Tests of EUKDB's genome set and package (since 2026-11-05: NCBI's genomes, no EukDetect): a species' genome by
NCBI's flags (contaminated never, a metagenome's only without better), RefSeq's category, level, contigs and date; its
RefSeq copy fetched where it has one; NCBI's placeholder species left out; a model species' reference whatever the
ranking; a species without an assembly from its records; --species. One genome's run with stand-in tools into its done
folder, then the package (format 2): genomes, copies whole, copies.tsv with the checks across species, features,
species, manifest; the prune.

    python3 -m unittest discover -s EUKDBbuild/tests -t EUKDBbuild/tests
"""

import functools
import gzip
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukdb import genomes, ncbi as ncbi_mod, odb as odb_mod, package, qc, rrna, sources  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST = "GCA_000146045.2"


class SetTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        self.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(d, "taxdump.tar.gz")))
        self.gb = sources.read_assemblies([fx.assembly_table(os.path.join(d, "genbank_fungi.txt"), fx.GENBANK)])
        self.rs = sources.read_assemblies([fx.assembly_table(os.path.join(d, "refseq_fungi.txt"), fx.REFSEQ)])

    def tearDown(self):
        self.tmp.cleanup()

    def test_one_genome_a_species(self):
        reps = genomes.representatives(self.gb, self.ncbi)
        best = {sp: row["accession"] for sp, (row, _o) in reps.items()}
        self.assertEqual(best[4932], YEAST)  # the reference genome
        self.assertEqual(best[27291], "GCA_000000013.1")  # contaminated never; a metagenome's only without better
        self.assertEqual(best[28985], "GCA_000002516.1")  # alike but newer
        self.assertNotIn(953, best)  # NCBI's placeholder
        self.assertNotIn(1080349, best)  # NCBI's hybrid
        self.assertEqual({r["accession"] for r in reps[27291][1]}, {"GCA_000000012.1"})
        self.assertEqual(self.gb[0]["group"], "fungi")
        self.assertNotIn("GCA_000000560.1", {r["accession"] for r in self.gb})  # a gene's record, no genome

    def test_jobs(self):
        model = [dict(species="Saccharomyces cerevisiae", ncbi_taxid="4932", domain="Eukaryota", genbank=YEAST,
                      refseq="GCF_000146045.2", assembly_name="R64", level="Complete Genome", genome_bases="12100000",
                      rdna_unit="J01355.1", ftp="-")]
        records = os.path.join(self.tmp.name, "na.fna.gz")
        with gzip.open(records, "wt") as fh:
            fh.write(">6253|Ascaris_lumbricoides|worm-Ascaris-1at2759-S1\nACGT\n")
        no_assembly = genomes.no_assembly_genomes(records, os.path.join(self.tmp.name, "na"), self.ncbi)
        self.assertEqual(list(no_assembly), ["EUK_6253"])
        jobs = genomes.jobs(self.ncbi, self.gb, self.rs, model, None, no_assembly)
        self.assertEqual(jobs[YEAST]["why"], ["model"])
        self.assertEqual((jobs[YEAST]["fetch"], jobs[YEAST]["rdna_unit"]), ("GCF_000146045.2", "J01355.1"))
        self.assertIn("GCA_000000003.1", jobs[YEAST]["others"])
        self.assertEqual(jobs["GCA_000000013.1"]["why"], ["ncbi"])
        self.assertEqual(jobs["GCA_000000013.1"]["itsx"], "F")
        self.assertEqual(jobs["EUK_6253"]["why"], ["no_assembly"])
        self.assertTrue(os.path.isfile(jobs["EUK_6253"]["path"]))
        # a GenBank genome with a RefSeq copy: the copy fetched
        jobs = genomes.jobs(self.ncbi, self.gb, self.rs)
        self.assertEqual(jobs[YEAST]["fetch"], "GCF_000146045.2")
        self.assertTrue(jobs[YEAST]["url"].endswith("GCF_000146045.2/GCF_000146045.2_genomic.fna.gz"))
        # --species
        inside = lambda t: 4930 in set(self.ncbi.ancestors(int(t)))  # noqa: E731
        self.assertEqual(sorted(genomes.jobs(self.ncbi, self.gb, self.rs, inside=inside)),
                         ["GCA_000000013.1", YEAST])


class RunTest(unittest.TestCase):
    def test_a_genome_into_the_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            root = fx.odb_world(tmp)
            odb = odb_mod.Odb(root, fx.LOCI)
            seq, gff, spans = fx.genome()
            genomes_dir = os.path.join(tmp, "out", "genomes")
            os.makedirs(genomes_dir)
            with gzip.open(genomes.genome_file(genomes_dir, YEAST), "wt") as fh:
                fh.write(f">c1\n{seq}\n")
            work = os.path.join(tmp, "work")
            settings = dict(work=work, genomes=genomes_dir, miniprot=fx.miniprot(tmp, gff),
                            hmmsearch=fx.hmmsearch(tmp), threads=1, staged_only=True, odb=root, loci=fx.LOCI,
                            hmm_db=odb.hmm_db(os.path.join(tmp, "loci.hmm")),
                            odb_lineages=odb_mod.species_lineages(odb, ncbi), tiara=None, rrna=None)
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            jobs = genomes.jobs(ncbi, gb)
            acc, summary, why = genomes.run((jobs[YEAST], settings))
            self.assertEqual((acc, why, summary["loci"], summary["status"]),
                             (YEAST, None, 2, {"Duplicated": 1, "Single": 1, "Missing": 1}))
            self.assertTrue(genomes.is_done(work, YEAST))
            self.assertEqual(genomes.run((jobs["GCA_000000013.1"], settings))[2], "not fetched yet")
            # the same copy in a second species (S. paradoxus' genome = the yeast's): the checks find them alike
            other = "GCA_000000013.1"
            with gzip.open(genomes.genome_file(genomes_dir, other), "wt") as fh:
                fh.write(f">c1\n{seq}\n")
            self.assertIsNone(genomes.run((jobs[other], settings))[2])
            collected = package.collect(ncbi, jobs, gb, work, fx.vsearch(tmp), 1)
            rows = {(r["accession"], r["locus"]): r for r in collected["copy_rows"]}
            self.assertEqual(rows[(YEAST, "1at2759")]["status"], "Duplicated")
            self.assertEqual(rows[(YEAST, "1at2759")]["same_as"], "4930")  # identical in two species: their genus
            self.assertEqual(rows[(YEAST, "1at2759")]["own_species"], "t")  # nearest agree to the genus
            loci_table = os.path.join(tmp, "loci.tsv")
            with open(loci_table, "w") as fh:
                fh.write("locus\tfungi_share\tprotist_share\tmedian_length\n")
                fh.writelines(f"{locus}\t1\t1\t400\n" for locus in fx.LOCI)
            out = package.write(os.path.join(tmp, "out", "EUKDB_2026.11"), "2026.11", ncbi, collected,
                                loci_table, None, genomes_dir, {"test": True})
            self.assertIsNone(package.check(out, deep=True))
            g = {r["accession"]: r for r in package.read_table(out, "genomes.tsv")}
            self.assertEqual((g[YEAST]["source"], g[YEAST]["representative"], g[YEAST]["loci"], g[YEAST]["local"]),
                             ("ncbi", "t", "2", f"genomes/{YEAST}.fna.gz"))
            self.assertEqual((g["GCA_000000003.1"]["source"], g["GCA_000000003.1"]["loci"]), ("assembly", "0"))
            a, b = spans["1at2759"]
            self.assertEqual(package.read_copies(out, "1at2759")[YEAST], seq[a - 1:b].upper())
            feats = [(f["kind"], f["start"], f["end"]) for f in package.read_table(out, "euk30/features.tsv")
                     if f["accession"] == YEAST and f["locus"] == "1at2759"]
            self.assertIn(("intron", "180", "480"), feats)
            self.assertIn("soft_masked", {k for k, _s, _e in feats})
            copies = package.read_table(out, "euk30/copies.tsv")
            self.assertEqual({c["locus"] for c in copies if c["accession"] == YEAST}, {"1at2759", "2at2759"})
            self.assertEqual(set(copies[0]), set(package.COPY_COLUMNS))
            with open(os.path.join(out, "manifest.json")) as fh:
                manifest = json.load(fh)
            self.assertEqual((manifest["format"], manifest["counts"]["species"], manifest["counts"]["copies"]),
                             (2, 2, 4))
            sp = {r["taxid"]: r for r in package.read_table(out, "species.tsv")}
            self.assertEqual((sp["4932"]["representative"], sp["4932"]["source"]), (YEAST, "ncbi"))
            # the prune: the model species' genome kept first
            size = os.path.getsize(genomes.genome_file(genomes_dir, YEAST))
            kept = genomes.prune(genomes_dir, work, size, set(), {YEAST})
            self.assertEqual(sorted(kept), sorted([(YEAST, size, "kept"), (other, size, "deleted")]))


class QcTest(unittest.TestCase):
    def test_checks_across_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            rng = __import__("random").Random(3)
            base = "".join(rng.choice("ACGT") for _ in range(300))
            near = base[:150] + "".join("A" if c != "A" else "C" for c in base[150:155]) + base[155:]  # 98.3%
            foreign = base[:297] + "TTT"  # 99%: Ascaris' copy as the yeast's (a contaminant)
            cds = {("Y", "1"): base, ("P", "1"): near, ("A", "1"): foreign, ("K", "1"): "".join(
                rng.choice("ACGT") for _ in range(300))}
            species = {"Y": 4932, "P": 27291, "A": 6253, "K": 28985}
            got = qc.check(cds, species, ncbi, fx.vsearch(tmp))
            self.assertEqual(got[("Y", "1")]["own_species"], "t")  # its nearest: S. paradoxus and Ascaris; genus
            self.assertEqual(got[("A", "1")]["own_species"], "f")  # Ascaris' nearest: yeasts, another kingdom
            self.assertEqual(got[("A", "1")]["same_as"], "2759")  # 99% to the yeast: their common ancestor
            self.assertEqual(got[("K", "1")]["own_species"], "-")  # no hit
            lengths = {(str(i), "1"): 100 + i for i in range(40)}
            flags = qc.check({k: "A" * 3 for k in lengths}, {str(i): 4932 for i in range(40)}, ncbi, None, 1, lengths)
            self.assertEqual((flags[("0", "1")]["length"], flags[("39", "1")]["length"], flags[("20", "1")]["length"]),
                             ("short", "long", "ok"))
            self.assertEqual(qc.lca(ncbi, [4932, 27291]), 4930)
            self.assertEqual(qc.depth(ncbi, 4932, 28985), 4)  # to the family


if __name__ == "__main__":
    unittest.main()
