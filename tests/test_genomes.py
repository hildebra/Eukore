#!/usr/bin/env python3
"""Tests of Eukore's genome set and package (since 2026-11-05: NCBI's genomes, no EukDetect): a species' genome by
NCBI's flags (contaminated never, a metagenome's only without better), RefSeq's category, level, contigs and date; its
RefSeq copy fetched where it has one; NCBI's placeholder species left out; a model species' reference whatever the
ranking; a species without an assembly from its records; --species. One genome's run with stand-in tools into its done
folder, then the package (format 4): genomes, copies whole, copies.tsv with the checks across species and the
confidence, features, busco/loci.tsv's call rates, species with their PR2 lineage, manifest; the prune.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import annotation, genomes, labels, metadata, ncbi as ncbi_mod, odb as odb_mod, package, qc  # noqa: E402
from eukore import sources  # noqa: E402
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

    def test_named_hybrid_species_stay(self):
        # NCBI's "Genus x epithet" is a species (the strawberry); a cross of two named species is left out
        for name in ("Fragaria x ananassa", "Mentha x piperita", "Pyrus × bretschneideri"):
            self.assertIsNone(genomes.HYBRID.search(name), name)
        for name in ("Saccharomyces cerevisiae x Saccharomyces paradoxus", "Saccharomyces cerevisiae x S. eubayanus"):
            self.assertIsNotNone(genomes.HYBRID.search(name), name)

    def test_the_best_genome_that_can_be_fetched(self):
        # Kluyveromyces lactis' best-ranked assembly without an FTP path: its next one, not the species lost
        rows = [dict(r, ftp="na") if r["accession"] == "GCA_000002516.1" else r for r in self.gb]
        jobs = genomes.jobs(self.ncbi, rows, self.rs)
        self.assertIn("GCA_000002515.1", jobs)
        self.assertEqual(jobs["GCA_000002515.1"]["others"], ["GCA_000002516.1"])
        # a species the taxdump does not hold (a newer assembly table): its own fate
        fates = genomes.fates(self.ncbi, rows + [dict(rows[0], accession="GCA_9.1", species_taxid=999999)], jobs)
        self.assertEqual(fates[999999], "not in the taxdump")

    def test_a_source_that_cannot_be_fetched(self):
        # its site moved: the source is not ready, the fetch goes on (until 2026-11-06 the stage stopped)
        def gone(url, dest, md5=None):
            raise OSError(f"{url}: not found")
        state = sources.fetch(self.tmp.name, "ROD", gone, log=lambda *_a: None)
        self.assertFalse(state["ready"])
        self.assertIn("not found", state["problems"][0])
        self.assertIsNone(sources.ready(self.tmp.name, "ROD"))
        from eukore import quality  # noqa: PLC0415
        self.assertEqual([quality.number({"x": v}, "x") for v in ("0.25", "", "na", "-")], [0.25, 0.0, 0.0, 0.0])

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
    def test_copy_numbers(self):
        genomes_ = {"G1": dict(representative="t", species_taxid=1, source="ncbi", _searched={"A", "B"},
                               _rrna_copies={"rDNA_unit": (3, 1), "18S": (4, 0)}, _organelle_molecules={"plastid": 1}),
                    "T1": dict(representative="t", species_taxid=2, source="tsa", _searched={"A"}),
                    "E1": dict(representative="t", species_taxid=3, source="no_assembly", _searched={"A"}),
                    "G2": dict(representative="f", species_taxid=1, source="assembly")}
        copy_rows = [dict(accession="G1", locus="A", status="Duplicated", copies="3"),
                     dict(accession="G1", locus="B", status="Fragmented", copies="0"),
                     dict(accession="T1", locus="A", status="Duplicated", copies="2")]
        organelle = [("x", "ACGT", dict(source="assembly", accession="G1", organelle="plastid", gene="rbcL", copies=1,
                                        partial=False)),
                     ("y", "ACGT", dict(source="refseq", accession="-", organelle="plastid", gene="psbA", copies=2,
                                        partial=False))]
        loci = [dict(locus="A"), dict(locus="B"), dict(locus="M1")]  # M1: a marker none was searched for
        got = {(r["accession"], r["set"], r["element"]): (r["copies"], r["partial"], r["basis"])
               for r in package.copy_numbers(genomes_, copy_rows, loci, {"A": "busco", "B": "busco", "M1": "markers"},
                                             organelle)}
        self.assertEqual(got, {("G1", "busco", "A"): (3, 0, "assembly"), ("G1", "busco", "B"): (0, 1, "assembly"),
                               ("G1", "rrna", "rDNA_unit"): (3, 1, "assembly"), ("G1", "rrna", "18S"): (4, 0, "assembly"),
                               ("G1", "plastid", "molecules"): (1, 0, "assembly"),
                               ("G1", "plastid", "rbcL"): (1, 0, "assembly"),  # RefSeq's psbA: not the genome's
                               ("T1", "busco", "A"): (2, 0, "transcripts")})  # isoforms; EukDetect2's records: none
        from eukore import rrna as rrna_mod  # noqa: PLC0415
        c = dict(contig="c1", strand="+")
        rows = [dict(c, kind="18S", start=1, end=1800), dict(c, kind="18S", start=20001, end=21800),
                dict(c, kind="18S", start=50001, end=51800),  # an 18S alone
                dict(c, kind="rDNA", lsu_complete=True, start=1, end=6120, regions=[
                    ("18S", 0, 1800), ("5.8S", 2000, 2160), ("28S", 2300, 5600), ("5S", 6000, 6120)]),
                dict(c, kind="rDNA", lsu_complete=False, start=20001, end=24000, regions=[
                    ("18S", 0, 1800), ("5.8S", 2000, 2160), ("28S", 2300, 4000)]),
                dict(c, kind="rDNA", lsu_complete=True, start=90001, end=95000, regions=[  # its 18S under MIN_LEN
                    ("18S", 0, 700), ("5.8S", 900, 1060), ("28S", 1200, 4500)])]
        self.assertEqual(rrna_mod.element_counts(rows), {"rDNA_unit": (3, 1), "18S": (4, 1), "5.8S": (3, 0),
                                                         "28S": (3, 1), "5S": (1, 0)})

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
                            odb_lineages=odb_mod.species_lineages(odb, ncbi), tiara=None, rrna=None,
                            vsearch=fx.vsearch(tmp))
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            jobs = genomes.jobs(ncbi, gb)
            acc, summary, why = genomes.run((jobs[YEAST], settings))
            self.assertEqual((acc, why, summary["loci"], summary["status"]),
                             (YEAST, None, 2, {"Duplicated": 1, "Single": 1, "Missing": 1}))
            self.assertTrue(genomes.is_done(work, YEAST))
            # the duplicated locus' two copies alike (loci.dup_identity: the genome's two identical copies)
            done = {r["locus"]: r for r in package.done_table(work, YEAST, "loci.tsv")}
            self.assertEqual((done["1at2759"]["dup_identity"], done["2at2759"]["dup_identity"],
                              done["1at2759"]["copies"], done["2at2759"]["copies"]), ("100.0", "-", "2", "1"))
            # its assembly's metadata measured by the run (metadata.py; no sourmash: no sketch, none missing)
            folder = genomes.done_dir(work, YEAST)
            self.assertEqual(metadata.missing(folder), set())
            self.assertEqual(metadata.missing(folder, "sourmash"), {metadata.SKETCH})
            with open(os.path.join(folder, metadata.STATS)) as fh:
                stats = json.load(fh)
            self.assertEqual((stats["sequences"], stats["length"], stats["fetched"]), (1, len(seq), YEAST))
            self.assertTrue(genomes.is_done(work, YEAST, fx.LOCI))
            self.assertFalse(genomes.is_done(work, YEAST, fx.LOCI + ["4at2759"]))  # a BUSCO it was not searched for
            self.assertEqual(genomes.run((jobs["GCA_000000013.1"], settings))[2], "not fetched yet")
            # the same copy in a second species (S. paradoxus' genome = the yeast's): the checks find them alike
            other = "GCA_000000013.1"
            with gzip.open(genomes.genome_file(genomes_dir, other), "wt") as fh:
                fh.write(f">c1\n{seq}\n")
            self.assertIsNone(genomes.run((jobs[other], settings))[2])
            collected = package.collect(ncbi, jobs, gb, work, fx.vsearch(tmp), 1)
            package.genetics(collected, ncbi, odb.hmm, fx.hmmalign(tmp), work, log=lambda *_: None)
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
            out = package.write(os.path.join(tmp, "out", "Eukore_2026.11"), "2026.11", ncbi, collected,
                                loci, None, genomes_dir, {"test": True}, None, [labels.read_pr2(pr2)])
            self.assertIsNone(package.check(out, deep=True))
            g = {r["accession"]: r for r in package.read_table(out, "genomes.tsv")}
            self.assertEqual((g[YEAST]["source"], g[YEAST]["representative"], g[YEAST]["loci"], g[YEAST]["local"]),
                             ("ncbi", "t", "2", f"genomes/{YEAST}.fna.gz"))
            self.assertEqual((g["GCA_000000003.1"]["source"], g["GCA_000000003.1"]["loci"]), ("assembly", "0"))
            self.assertEqual((g[YEAST]["type_material"], g["GCA_000000003.1"]["type_material"]),
                             ("assembly from type material", "-"))  # NCBI's relation_to_type_material
            a, b = spans["1at2759"]
            self.assertEqual(package.read_copies(out, "1at2759")[YEAST], seq[a - 1:b].upper())
            # its coding sequence (the intron out) and protein (the stop out), format 5
            self.assertEqual(package.read_copies(out, "1at2759", "cds")[YEAST], fx.cds_of(fx.PROTEINS["1at2759"]))
            self.assertEqual(package.read_copies(out, "1at2759", "proteins")[YEAST], fx.PROTEINS["1at2759"])
            self.assertEqual(package.read_copies(out, "2at2759", "proteins")[YEAST], fx.PROTEINS["2at2759"])
            self.assertEqual((g[YEAST]["gencode"], g[YEAST]["duplication"]), ("1", "-"))  # 2 BUSCOs: too few
            # its intron: between codons 60 and 61 (180 coding bases before it), on the HMM's state 61 (hmmalign's
            # stand-in: every residue on its own position); both genomes hold it there
            (intron,) = [r for r in package.read_table(out, "busco/introns.tsv") if r["accession"] == YEAST]
            self.assertEqual((intron["locus"], intron["start"], intron["end"], intron["phase"], intron["class"],
                              intron["u12"], intron["aa_position"], intron["hmm_position"], intron["site"]),
                             ("1at2759", "180", "480", "0", "GT-AG", "f", "60", "61", "61.0"))
            (site,) = package.read_table(out, "busco/intron_sites.tsv")
            self.assertEqual((site["site"], site["genomes"], site["share"], site["share_Fungi"]),
                             ("61.0", "2", "1.000", "1.000"))
            self.assertEqual((g[YEAST]["introns"], g[YEAST]["median_intron_length"], g[YEAST]["gencode_check"]),
                             ("1", "300", "undecided"))  # the stand-in HMMs have no emissions: no codon seen
            self.assertTrue(os.path.isfile(os.path.join(out, "metadata", "codons.tsv.gz")))
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
            self.assertEqual((one["internal_stops"], one["dup_identity"]), ("0", "100.0"))
            with open(os.path.join(out, "manifest.json")) as fh:
                manifest = json.load(fh)
            self.assertEqual((manifest["format"], manifest["counts"]["species"], manifest["counts"]["copies"],
                              manifest["counts"]["busco_loci"]), (package.FORMAT, 2, 4, 3))
            self.assertEqual(manifest["inputs"]["references"]["PR2"]["species"], 1)
            # every genome's copy numbers in one table (since 2026-11-08): its BUSCOs, Missing ones 0; its organelle
            # molecules (none in the stand-in genome); no rRNA searched: no rRNA rows
            cn = {(r["accession"], r["set"], r["element"]): r
                  for r in package.read_gz(os.path.join(out, package.COPY_NUMBERS))}
            self.assertEqual([cn[(YEAST, "busco", l)]["copies"] for l in fx.LOCI], ["2", "1", "0"])
            self.assertEqual((cn[(YEAST, "mitochondrion", "molecules")]["copies"], cn[(YEAST, "busco", "1at2759")]["basis"]),
                             ("0", "assembly"))
            self.assertNotIn((YEAST, "rrna", "18S"), cn)
            self.assertEqual(manifest["counts"]["copy_numbers"], {"busco": 6, "mitochondrion": 2, "plastid": 2})
            # each BUSCO's calls over the package's species: two fungi, 1at2759 and 2at2759 called in both
            got = {r["locus"]: r for r in package.read_table(out, "busco/loci.tsv")}
            self.assertEqual([got[l]["called"] for l in fx.LOCI], ["2", "2", "0"])
            self.assertEqual((got["1at2759"]["duplicated"], got["1at2759"]["call_rate"], got["1at2759"]["call_Fungi"],
                              got["3at2759"]["call_rate"], got["1at2759"]["cutoff"], got["1at2759"]["hmm_length"]),
                             ("2", "1.000", "1.000", "0.000", "50.0", "120"))
            sp = {r["taxid"]: r for r in package.read_table(out, "species.tsv")}
            self.assertEqual((sp["4932"]["representative"], sp["4932"]["source"], sp["4932"]["type_material"]),
                             (YEAST, "ncbi", "assembly from type material"))
            self.assertEqual((sp["4932"]["supergroup"], sp["4932"]["group"], sp["4932"]["pr2_match"], sp["4932"]["ksgp"]),
                             ("Obazoa", "Fungi", "species", "-"))
            self.assertTrue(sp["4932"]["pr2"].endswith(";Saccharomyces;Saccharomyces_cerevisiae"))
            self.assertEqual(sp["27291"]["pr2_match"], "genus")  # S. paradoxus: PR2 names its genus alone
            self.assertEqual(sp["4932"]["ncbi_lineage"], "Eukaryota;Fungi;Ascomycota;Saccharomycetes;"
                                                         "Saccharomycetales;Saccharomycetaceae;Saccharomyces;"
                                                         "Saccharomyces cerevisiae")
            # a locus of Eukore's own markers (since 2026-11-07) in its own folder, with the same files
            marked = [dict(i, set="markers") if i["locus"] == "2at2759" else i for i in loci]
            out2 = package.write(os.path.join(tmp, "out", "Eukore_2026.12"), "2026.12", ncbi, collected, marked,
                                 None, genomes_dir, {"test": True}, None, [labels.read_pr2(pr2)])
            self.assertIsNone(package.check(out2, deep=True))
            self.assertEqual(set(package.read_copies(out2, "2at2759", "proteins", "markers")), {YEAST, other})
            self.assertFalse(os.path.isfile(os.path.join(out2, "busco", "copies", "2at2759.fna.gz")))
            self.assertEqual([r["locus"] for r in package.read_table(out2, "markers/loci.tsv")], ["2at2759"])
            self.assertNotIn("2at2759", {r["locus"] for r in package.read_table(out2, "busco/copies.tsv")})
            self.assertEqual({r["locus"] for r in package.read_table(out2, "markers/copies.tsv")}, {"2at2759"})
            with open(os.path.join(out2, "manifest.json")) as fh:
                counts = json.load(fh)["counts"]
            self.assertEqual((counts["busco_loci"], counts["marker_loci"]), (2, 1))
            # the prune: the model species' genome kept first
            size = os.path.getsize(genomes.genome_file(genomes_dir, YEAST))
            kept = genomes.prune(genomes_dir, work, size, set(), {YEAST})
            self.assertEqual(sorted(kept), sorted([(YEAST, size, "kept"), (other, size, "deleted")]))


    def test_a_genome_with_its_ncbi_annotation(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            root = fx.odb_world(tmp)
            odb = odb_mod.Odb(root, fx.LOCI)
            seq, gff, spans = fx.genome()
            genomes_dir = os.path.join(tmp, "out", "genomes")
            os.makedirs(genomes_dir)
            with gzip.open(genomes.genome_file(genomes_dir, YEAST), "wt") as fh:
                fh.write(f">c1\n{seq}\n")
            a, b = spans["1at2759"]
            ann = annotation.gff_file(genomes_dir, YEAST)  # the fetch stage's, beside the genome
            with gzip.open(ann, "wt") as fh:
                fh.write(f"c1\tRefSeq\tCDS\t{a}\t{a + 179}\t.\t+\t0\tParent=rna-1\n"
                         f"c1\tRefSeq\tCDS\t{a + 480}\t{b}\t.\t+\t0\tParent=rna-1\n")
            work = os.path.join(tmp, "work")
            settings = dict(work=work, genomes=genomes_dir, miniprot=fx.miniprot(tmp, gff),
                            hmmsearch=fx.hmmsearch(tmp), threads=1, staged_only=True, odb=root, loci=fx.LOCI,
                            hmm_db=odb.hmm_db(os.path.join(tmp, "loci.hmm")),
                            odb_lineages=odb_mod.species_lineages(odb, ncbi), tiara=None, rrna=None)
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            job = dict(genomes.jobs(ncbi, gb)[YEAST], annotated=True)  # the summary: protein-coding genes
            self.assertTrue(annotation.gff_url(job).endswith("_genomic.gff.gz"))
            acc, summary, why = genomes.run((job, settings))
            self.assertEqual((why, summary["annotated"], summary["annotation"], summary["annotated_models"]),
                             (None, True, "ncbi", 1))
            done = genomes.done_dir(work, YEAST)
            models = {r["locus"]: r["model"] for r in fx.read_tsv(os.path.join(done, "loci.tsv"))}
            self.assertEqual(models, {"1at2759": "ncbi", "2at2759": "miniprot", "3at2759": "-"})
            compared = {r["locus"]: r for r in fx.read_tsv(os.path.join(done, annotation.TABLE))}
            self.assertEqual((compared["1at2759"]["exact_introns"], compared["2at2759"]["transcript"],
                              compared["1at2759"]["source"]), ("1", "-", "ncbi"))
            self.assertFalse(os.path.isfile(ann))  # taken and compared: deleted
            self.assertEqual(summary["annotation_none"], [])
            self.assertIsNone(genomes.better_annotation(work, job))  # NCBI's taken: none better
            self.assertIsNone(annotation.gff_url(dict(job, annotated=False)))
            # no NCBI annotation, Ensembl's (a GTF here: TOGA2's form) beside the genome: its models, `ensembl`
            other = dict(job, annotated=False, annotations=[("ensembl", "https://x/genes.gff3.gz")])
            ann = annotation.gff_file(genomes_dir, YEAST, "ensembl")
            with gzip.open(ann, "wt") as fh:
                fh.write(f'c1\tensembl\tCDS\t{a}\t{a + 179}\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n'
                         f'c1\tensembl\tCDS\t{a + 480}\t{b}\t.\t+\t0\tgene_id "g1"; transcript_id "t1";\n')
            # as if its run had found no NCBI models: Ensembl's are better, here, so it is run again (staged)
            genomes.note_annotation_none(work, other, ["ncbi"])
            with open(os.path.join(done, "done.json")) as fh:
                record = json.load(fh)
            with open(os.path.join(done, "done.json"), "w") as fh:
                json.dump(dict(record, annotation="-"), fh)
            self.assertEqual(genomes.better_annotation(work, other), "ensembl")
            self.assertTrue(genomes.rerun_for_annotation(work, other, genomes_dir, False))
            with open(os.path.join(done, metadata.STATS)) as fh:
                measured = fh.read()
            acc, summary, why = genomes.run((other, settings))
            with open(os.path.join(done, metadata.STATS)) as fh:
                self.assertEqual(fh.read(), measured)  # its metadata kept, not measured again
            self.assertIsNone(genomes.better_annotation(work, other))
            self.assertEqual((why, summary["annotation"], summary["annotated_models"]), (None, "ensembl", 1))
            models = {r["locus"]: r["model"] for r in fx.read_tsv(os.path.join(done, "loci.tsv"))}
            self.assertEqual(models["1at2759"], "ensembl")
            compared = {r["locus"]: r for r in fx.read_tsv(os.path.join(done, annotation.TABLE))}
            self.assertEqual((compared["1at2759"]["source"], compared["1at2759"]["transcript"]), ("ensembl", "t1"))
            self.assertFalse(os.path.isfile(ann))


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

    def test_duplication(self):
        def rows(single, dups):
            return [dict(status="Single")] * single + [dict(status="Duplicated", dup_identity=d) for d in dups]
        got = package.duplication_of(rows(60, ["99.9"] * 30 + ["98.0"] * 10))  # a diploid's two haplotypes
        self.assertEqual((got["duplicated_share"], got["dup_identity_median"], got["duplication"]),
                         ("0.400", "99.9", "haplotypes"))
        self.assertEqual(package.duplication_of(rows(60, ["85.0"] * 20 + ["<50"] * 20))["duplication"],
                         "wgd_or_polyploid")
        self.assertEqual(package.duplication_of(rows(95, ["99.9"] * 5))["duplication"], "low")
        self.assertEqual(package.duplication_of(rows(5, ["99.9"] * 5))["duplication"], "-")  # too few called
        self.assertEqual(package.duplication_of(rows(60, ["-"] * 30))["duplication"], "-")  # run before 0.2


class RecheckTest(unittest.TestCase):
    """A rerun on genomes already done (since 2026-11-07): a done genome run again where better gene models than its
    run took are there; a species' genome another assembly took the place of deleted once that one is done."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.work, self.genomes = os.path.join(self.tmp.name, "work"), os.path.join(self.tmp.name, "genomes")
        os.makedirs(self.genomes)

    def tearDown(self):
        self.tmp.cleanup()

    def done(self, acc, summary, table=None):
        folder = genomes.done_dir(self.work, acc)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "done.json"), "w") as fh:
            json.dump(dict(summary, accession=acc), fh)
        if table is not None:
            with open(os.path.join(folder, annotation.TABLE), "w") as fh:
                fh.write("\t".join(annotation.ANNOTATION_COLUMNS) + "\n" + "".join(table))

    @staticmethod
    def job(acc, others=(("ensembl", "https://ensembl.example/x.gff3.gz"),)):
        return dict(accession=acc, fetch=acc, url=f"https://ftp.example/{acc}_x_genomic.fna.gz", annotated=None,
                    annotations=list(others))

    def test_better_gene_models(self):
        a, j = "GCA_000000001.1", self.job("GCA_000000001.1")
        self.assertIsNone(genomes.better_annotation(self.work, j))  # not done
        # a run before 2026-11-07 without NCBI's models: NCBI never asked (no annotation.tsv), asked first
        self.done(a, {"annotated": False})
        self.assertEqual(genomes.annotation_record(self.work, j), ("-", set()))
        self.assertEqual(genomes.better_annotation(self.work, j), "ncbi")
        # NCBI's 404 (an empty annotation.tsv): Ensembl's, which an index names since
        self.done(a, {"annotated": False}, table=[])
        self.assertEqual(genomes.annotation_record(self.work, j), ("-", {"ncbi"}))
        self.assertEqual(genomes.better_annotation(self.work, j), "ensembl")
        # one with NCBI's models: none better
        self.done(a, {"annotated": True, "ncbi_models": 3}, table=["x\n"])
        self.assertEqual(genomes.annotation_record(self.work, j), ("ncbi", set()))
        self.assertIsNone(genomes.better_annotation(self.work, j))
        # a run with Ensembl's; VEuPathDB's (a key given since) comes first
        self.done(a, {"annotation": "ensembl", "annotation_none": ["ncbi"]})
        self.assertIsNone(genomes.better_annotation(self.work, j))
        veu = self.job(a, (("veupathdb", "https://veupathdb.example/x.gff"), ("ensembl", "https://e.example/x")))
        self.assertEqual(genomes.better_annotation(self.work, veu), "veupathdb")
        # found to have none (the fetch stage's): not asked again; an old record carried over
        genomes.note_annotation_none(self.work, veu, ["veupathdb"])
        self.assertIsNone(genomes.better_annotation(self.work, veu))
        self.done(a, {"annotated": False}, table=[])
        genomes.note_annotation_none(self.work, j, ["ensembl"])
        self.assertEqual(genomes.annotation_record(self.work, j), ("-", {"ncbi", "ensembl"}))
        self.assertIsNone(genomes.better_annotation(self.work, j))
        self.assertEqual(annotation.none_before(veu, "ensembl"), ["ncbi", "veupathdb"])
        self.assertEqual(annotation.none_before(veu, None), ["ncbi", "veupathdb", "ensembl"])

    def test_run_again_where_its_genome_and_annotation_are_here(self):
        a, j = "GCA_000000002.1", self.job("GCA_000000002.1")
        self.done(a, {"annotation": "-", "annotation_none": ["ncbi"]})
        self.assertFalse(genomes.rerun_for_annotation(self.work, j, self.genomes, True))  # no genome here
        with gzip.open(genomes.genome_file(self.genomes, a), "wt") as fh:
            fh.write(">c1\nACGT\n")
        self.assertTrue(genomes.rerun_for_annotation(self.work, j, self.genomes, True))  # the run fetches it
        self.assertFalse(genomes.rerun_for_annotation(self.work, j, self.genomes, False))  # staged: not here yet
        with gzip.open(annotation.gff_file(self.genomes, a, "ensembl"), "wt") as fh:
            fh.write("c1\tx\tCDS\t1\t3\t.\t+\t0\tParent=t1\n")
        self.assertTrue(genomes.rerun_for_annotation(self.work, j, self.genomes, False))

    def test_a_superseded_genome_deleted_once_its_replacement_is_done(self):
        old, new, other = "GCA_000000003.1", "TOGA2_HLtesSpe1", "GCA_000000004.1"
        self.assertEqual(genomes.supersede(self.work, {old: dict(species_taxid=1), other: dict(species_taxid=2)}), [])
        for acc in (old, other):
            with gzip.open(genomes.genome_file(self.genomes, acc), "wt") as fh:
                fh.write(">c1\nACGT\n")
        with gzip.open(annotation.gff_file(self.genomes, old), "wt") as fh:
            fh.write("#\n")
        # TOGA2's in place of NCBI's; a --species plan without species 2 keeps its row
        plan = {new: dict(species_taxid=1)}
        self.assertEqual(genomes.supersede(self.work, plan), [(old, new)])
        self.assertEqual(genomes.supersede(self.work, plan), [])  # said once
        self.assertEqual(genomes.read_rows(os.path.join(self.work, genomes.GENOME_SET)), [["1", new], ["2", other]])
        self.assertEqual(genomes.drop_superseded(self.genomes, self.work, plan), [])  # the new one not done yet
        self.assertTrue(os.path.isfile(genomes.genome_file(self.genomes, old)))
        self.done(new, {"annotation": "toga2"})
        self.assertEqual(genomes.drop_superseded(self.genomes, self.work, plan), [old])
        self.assertFalse(os.path.isfile(genomes.genome_file(self.genomes, old)))
        self.assertFalse(os.path.isfile(annotation.gff_file(self.genomes, old)))
        self.assertTrue(os.path.isfile(genomes.genome_file(self.genomes, other)))
        self.assertEqual(genomes.drop_superseded(self.genomes, self.work, plan), [])
        # chosen again (a swap undone): its file kept
        with gzip.open(genomes.genome_file(self.genomes, old), "wt") as fh:
            fh.write(">c1\nACGT\n")
        self.assertEqual(genomes.drop_superseded(self.genomes, self.work, {old: dict(species_taxid=1)}), [])


if __name__ == "__main__":
    unittest.main()
