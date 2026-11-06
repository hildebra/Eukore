#!/usr/bin/env python3
"""Tests of EUKDB's genome set and package (since 2026-11-05: NCBI's genomes, no EukDetect): a species' genome by
NCBI's flags (contaminated never, a metagenome's only without better), RefSeq's category, level, contigs and date; its
RefSeq copy fetched where it has one; NCBI's placeholder species left out; a model species' reference whatever the
ranking; a species without an assembly from its records; --species. One genome's run with stand-in tools into its done
folder, then the package (format 4): genomes, copies whole, copies.tsv with the checks across species and the
confidence, features, busco/loci.tsv's call rates, species with their PR2 lineage, manifest; the prune.

    python3 -m unittest discover -s tests -t tests
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
from eukdb import genomes, labels, ncbi as ncbi_mod, odb as odb_mod, package, qc, rrna, sources  # noqa: E402
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

    def test_small_for_its_clade(self):
        class Tax:  # genus 10 (species 1-4, ~1 Gb), genus 20 in family 100 with genus 10 (species 5: one 4 Mb genome)
            parent = {1: 10, 2: 10, 3: 10, 4: 10, 5: 20, 10: 100, 20: 100, 100: 1000}
            rank = {**{t: "species" for t in (1, 2, 3, 4, 5)}, 10: "genus", 20: "genus", 100: "family", 1000: "order"}

            def ancestors(self, t):
                while True:
                    yield t
                    if t not in self.parent:
                        return
                    t = self.parent[t]
        rows = {sp: [dict(accession=f"A{sp}", size=size)] for sp, size in ((1, 900e6), (2, 1.1e9), (3, 1e9))}
        rows[4] = [dict(accession="A4a", size=40e6), dict(accession="A4b", size=950e6)]
        rows[5] = [dict(accession="A5", size=4e6, flags=["partial genome"])]
        out = genomes.flag_small(rows, Tax())
        self.assertEqual(out[4][0]["flags"], [genomes.SMALL])  # a 40 Mb part of a 1 Gb genus
        self.assertNotIn("flags", out[4][1])
        self.assertEqual(out[5][0]["flags"], ["partial genome", genomes.SMALL])  # its genus alone: by its family
        self.assertNotIn("flags", out[1][0])

    def test_partial_genome_flag(self):
        table = os.path.join(self.tmp.name, "genbank_vertebrate_other.txt")
        with open(table, "w") as fh:
            fh.write("# assembly_accession\ttaxid\tspecies_taxid\torganism_name\tgenome_rep\tgenome_size\n"
                     "GCA_1.1\t4932\t4932\tSaccharomyces cerevisiae\tPartial\t5000000\n")
        self.assertEqual(sources.read_assemblies([table])[0]["flags"], ["partial genome"])

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
            self.assertTrue(genomes.is_done(work, YEAST, fx.LOCI))
            self.assertFalse(genomes.is_done(work, YEAST, fx.LOCI + ["4at2759"]))  # a BUSCO it was not searched for
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
            loci = [dict(locus=l, name=odb.names.get(l, "-"), cutoff=odb.cutoff[l], hmm_length=odb.length[l])
                    for l in odb.loci]
            pr2 = os.path.join(tmp, "pr2.tax")
            with open(pr2, "w") as fh:
                fh.write("A1\tEukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;Saccharomycetales;"
                         "Saccharomycetaceae;Saccharomyces;Saccharomyces_cerevisiae;\n")
            out = package.write(os.path.join(tmp, "out", "EUKDB_2026.11"), "2026.11", ncbi, collected,
                                loci, None, genomes_dir, {"test": True}, None, [labels.read_pr2(pr2)])
            self.assertIsNone(package.check(out, deep=True))
            g = {r["accession"]: r for r in package.read_table(out, "genomes.tsv")}
            self.assertEqual((g[YEAST]["source"], g[YEAST]["representative"], g[YEAST]["loci"], g[YEAST]["local"]),
                             ("ncbi", "t", "2", f"genomes/{YEAST}.fna.gz"))
            self.assertEqual((g["GCA_000000003.1"]["source"], g["GCA_000000003.1"]["loci"]), ("assembly", "0"))
            a, b = spans["1at2759"]
            self.assertEqual(package.read_copies(out, "1at2759")[YEAST], seq[a - 1:b].upper())
            feats = [(f["kind"], f["start"], f["end"]) for f in package.read_table(out, "busco/features.tsv")
                     if f["accession"] == YEAST and f["locus"] == "1at2759"]
            self.assertIn(("intron", "180", "480"), feats)
            self.assertIn("soft_masked", {k for k, _s, _e in feats})
            copies = package.read_table(out, "busco/copies.tsv")
            self.assertEqual({c["locus"] for c in copies if c["accession"] == YEAST}, {"1at2759", "2at2759"})
            self.assertEqual(set(copies[0]), set(package.COPY_COLUMNS))
            one = next(c for c in copies if c["accession"] == YEAST and c["locus"] == "1at2759")
            self.assertEqual((one["contig_edge"], one["gaps"]), (str(spans["1at2759"][0] - 1), "0"))
            self.assertEqual((one["confidence"], one["flags"]), ("medium", "duplicated"))  # the yeast's two copies
            with open(os.path.join(out, "manifest.json")) as fh:
                manifest = json.load(fh)
            self.assertEqual((manifest["format"], manifest["counts"]["species"], manifest["counts"]["copies"],
                              manifest["counts"]["busco_loci"]), (4, 2, 4, 3))
            self.assertEqual(manifest["inputs"]["references"]["PR2"]["species"], 1)
            # each BUSCO's calls over the package's species: two fungi, 1at2759 and 2at2759 called in both
            got = {r["locus"]: r for r in package.read_table(out, "busco/loci.tsv")}
            self.assertEqual([got[l]["called"] for l in fx.LOCI], ["2", "2", "0"])
            self.assertEqual((got["1at2759"]["duplicated"], got["1at2759"]["call_rate"], got["1at2759"]["call_Fungi"],
                              got["3at2759"]["call_rate"], got["1at2759"]["cutoff"], got["1at2759"]["hmm_length"]),
                             ("2", "1.000", "1.000", "0.000", "50.0", "120"))
            sp = {r["taxid"]: r for r in package.read_table(out, "species.tsv")}
            self.assertEqual((sp["4932"]["representative"], sp["4932"]["source"]), (YEAST, "ncbi"))
            self.assertEqual((sp["4932"]["supergroup"], sp["4932"]["pr2_match"], sp["4932"]["ksgp"]),
                             ("Fungi", "species", "-"))
            self.assertTrue(sp["4932"]["pr2"].endswith(";Saccharomyces;Saccharomyces_cerevisiae"))
            self.assertEqual(sp["27291"]["pr2_match"], "genus")  # S. paradoxus: PR2 names its genus alone
            self.assertEqual(sp["4932"]["ncbi_lineage"], "Eukaryota;Fungi;Ascomycota;Saccharomycetes;"
                                                         "Saccharomycetales;Saccharomycetaceae;Saccharomyces;"
                                                         "Saccharomyces cerevisiae")
            # the prune: the model species' genome kept first
            size = os.path.getsize(genomes.genome_file(genomes_dir, YEAST))
            kept = genomes.prune(genomes_dir, work, size, set(), {YEAST})
            self.assertEqual(sorted(kept), sorted([(YEAST, size, "kept"), (other, size, "deleted")]))


class FcsTest(unittest.TestCase):
    def test_fcs_gx_on_the_windows(self):
        # run_gx's stand-in: a TRIM of 1,000-3,000 on the first window (its report beside the FASTA it was given)
        with tempfile.TemporaryDirectory() as tmp:
            work = os.path.join(tmp, "work")
            folder = genomes.done_dir(work, YEAST)
            os.makedirs(folder)
            with gzip.open(os.path.join(folder, "windows.fna.gz"), "wt") as fh:
                fh.write(">1at2759|c1|1001-61000\n" + "A" * 60000 + "\n>2at2759|c1|9001-70000\n" + "C" * 61000 + "\n")
            run_gx = fx.tool(tmp, "run_gx", "#!/bin/sh\nwhile [ $# -gt 0 ]; do case $1 in --out-dir) out=$2;; "
                                            "--tax-id) t=$2;; esac; shift; done\n"
                                            "printf '#seq_id\\tstart_pos\\tend_pos\\tseq_len\\taction\\tdiv\\t"
                                            "agg_cont_cov\\ttop_tax_name\\nw0\\t1000\\t3000\\t60000\\tTRIM\\t"
                                            "prok:CFB group bacteria\\t3\\tBacteroides\\n' "
                                            "> $out/w.$t.fcs_gx_report.txt\n")
            ranges = genomes.fcs_screen(work, YEAST, 559292, "/db/gxdb/all", run_gx)
            self.assertEqual([(r["window"], r["action"]) for r in ranges], [("1at2759|c1|1001-61000", "TRIM")])
            self.assertTrue(os.path.isfile(os.path.join(folder, "fcs.tsv")))
            # the range is 2,000-4,000 on the contig: a copy over it flagged, one elsewhere or another locus not
            self.assertEqual(genomes.fcs_of_copy(ranges, "1at2759", "c1", 3500, 5000),
                             "TRIM:prok:CFB group bacteria:Bacteroides")
            self.assertEqual(genomes.fcs_of_copy(ranges, "1at2759", "c1", 5000, 6000), "-")
            self.assertEqual(genomes.fcs_of_copy(ranges, "2at2759", "c1", 3500, 5000), "-")
            failing = fx.tool(tmp, "run_gx_fail", "#!/bin/sh\nexit 3\n")
            self.assertIsNone(genomes.fcs_screen(work, YEAST, 559292, "/db", failing))


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
