#!/usr/bin/env python3
"""Tests of the intraspecific panel and its barcode gap (since 2026-11-07): the further assemblies chosen (other
samples, unflagged first), a marker's identity within and between species, and the gap's verdict.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, intraspecific, ncbi as ncbi_mod, sources  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST = "GCA_000146045.2"


def mutate(seq, n, seed):
    rng = random.Random(seed)
    s = list(seq)
    for i in rng.sample(range(len(s)), n):
        s[i] = "A" if s[i] != "A" else "C"
    return "".join(s)


class PanelTest(unittest.TestCase):
    def test_further_assemblies_chosen(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            job_map = genomes.jobs(ncbi, gb)
            extras = genomes.intraspecific_jobs(ncbi, gb, [], job_map, {4932: None, 6253: 2}, extras=3)
            self.assertTrue(extras)
            self.assertNotIn(YEAST, extras)  # the species' genome is not a further one
            self.assertTrue(all(j["why"] == ["intraspecific"] and j["species_taxid"] == 4932 for j in extras.values()))
            self.assertEqual(genomes.intraspecific_jobs(ncbi, gb, [], job_map, {4932: None}, extras=0), {})
            self.assertTrue(genomes.same_sample(dict(biosample="SAMN1"), dict(biosample="samn1")))
            self.assertFalse(genomes.same_sample(dict(strain="na"), dict(strain="na")))
            self.assertIn(4932, genomes.read_panel())

    def test_barcode_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            rng = random.Random(1)
            base = "".join(rng.choice("ACGT") for _ in range(400))
            # yeast: its genome and two strains 99.5% alike; its nearest other species 95% (resolved). A marker the
            # strains differ in more than the species do: not resolved
            rep = {"Y": {"m1": base, "m2": base}}
            extra = {"Y2": {"m1": mutate(base, 2, 2), "m2": mutate(base, 40, 3)},
                     "Y3": {"m1": mutate(base, 2, 4), "m2": mutate(base, 40, 5)}}
            species_of = {"Y": 4932, "Y2": 4932, "Y3": 4932}
            others = {"m1": {"27291|P": mutate(base, 20, 6)}, "m2": {"27291|P": mutate(base, 20, 7)}}
            rows = intraspecific.barcode_gap({4932: "Y"}, rep, extra, species_of, {}, others, ncbi, fx.vsearch(tmp))
            got = {r["marker"]: r for r in rows}
            self.assertEqual((got["m1"]["assemblies"], got["m1"]["inter_max"], got["m1"]["resolved"]), (3, "95.0", "t"))
            self.assertEqual(got["m2"]["resolved"], "f")
            self.assertGreater(float(got["m1"]["gap"]), 0)
            # a BUSCO's nearest other species given (qc.check's): not measured again
            rows = intraspecific.barcode_gap({4932: "Y"}, rep, extra, species_of, {(4932, "m1"): 99.9}, {}, ncbi,
                                             fx.vsearch(tmp))
            self.assertEqual({r["marker"]: r["resolved"] for r in rows}, {"m1": "f", "m2": "-"})

    def test_markers_of_a_genome(self):
        unit = dict(kind="rDNA", sequence="A" * 10 + "C" * 20 + "G" * 5 + "T" * 15,
                    regions=[("18S", 0, 10), ("ITS1", 10, 20), ("5.8S", 20, 30), ("ITS2", 30, 35)], lsu_complete=True)
        got = intraspecific.first_rrna([dict(kind="18S", sequence="TTTT", regions=None), unit])
        self.assertEqual(got, {"18S": "A" * 10, "ITS": "C" * 20 + "G" * 5})
        org = intraspecific.organelle_seqs([dict(gene="COX1", organelle="mitochondrion", sequence="ATGAAA", exons="1",
                                                 start="1", end="6", parts="-"),
                                            dict(gene="SSU", organelle="mitochondrion", sequence="ACGT")])
        self.assertEqual(org, {"COX1": "ATGAAA", "mito_SSU": "ACGT"})


class TsaTest(unittest.TestCase):
    """Transcriptomes for the species without a genome (genomes.tsa_jobs; since 2026-11-07)."""

    def test_transcriptomes_chosen(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            job_map = genomes.jobs(ncbi, gb)
            job_map["EUK_6253"] = dict(accession="EUK_6253", species_taxid=6253, why=["no_assembly"])
            path = os.path.join(tmp, "tsa_sets.tsv")
            with open(path, "w") as fh:
                fh.write("\t".join(genomes.TSA_FIELDS) + "\n")
                for row in (("GAAA01", "6253", "Ascaris lumbricoides", "5000000", "2015-01-01", "", "false",
                             "ftp.ebi.ac.uk/pub/databases/ena/tsa/public/gaa/GAAA01.fasta.gz"),
                            ("GAAB01", "6253", "Ascaris lumbricoides", "9000000", "2014-01-01", "", "false",
                             "ftp.ebi.ac.uk/pub/databases/ena/tsa/public/gaa/GAAB01.fasta.gz"),
                            ("GAAC01", "6253", "Ascaris lumbricoides", "99000000", "2016-01-01", "", "true",
                             "ftp.ebi.ac.uk/x/GAAC01.fasta.gz"),  # an environmental sample: not taken
                            ("GAAD01", "4932", "Saccharomyces cerevisiae", "9000000", "2016-01-01", "", "false",
                             "ftp.ebi.ac.uk/x/GAAD01.fasta.gz")):  # a species with a genome: not taken
                    fh.write("\t".join(row) + "\n")
            got = genomes.tsa_jobs(ncbi, genomes.read_tsa(path), job_map)
            self.assertEqual(sorted(got), ["GAAB01"])  # the most bases
            job = got["GAAB01"]
            self.assertEqual((job["why"], job["kind"], job["species_taxid"], job["url"]),
                             (["tsa"], "transcriptome", 6253,
                              "https://ftp.ebi.ac.uk/pub/databases/ena/tsa/public/gaa/GAAB01.fasta.gz"))
            self.assertNotIn("EUK_6253", job_map)  # the transcriptome stands for EukDetect2's records
            os.makedirs(os.path.join(tmp, "tsa"))
            os.replace(path, genomes.tsa_path(tmp))
            state = genomes.tsa_state(tmp)  # the manifest's inputs.tsa
            self.assertEqual((state["file"], state["sets"], len(state["md5"])), ("tsa_sets.tsv", 3, 32))
            self.assertIsNone(genomes.tsa_state(os.path.join(tmp, "none")))

    def test_a_transcripts_confidence(self):
        from eukore import metadata, package
        r = dict(status="Duplicated", own_species="t", transcript="t", gaps=0, contig_edge="5000")
        self.assertEqual(package.copy_confidence(r), ("medium", "transcript"))  # isoforms: not "duplicated"
        self.assertEqual(package.copy_confidence(dict(r, own_species="f"))[0], "low")
        self.assertEqual(metadata.genome_source({}, dict(why=["tsa"])), "ENA TSA")


if __name__ == "__main__":
    unittest.main()
