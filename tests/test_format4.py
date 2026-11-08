#!/usr/bin/env python3
"""Tests of format 4's additions (since 2026-11-05): each species' PR2, EUKARYOME and KSGP lineage beside its NCBI taxid
(by its name, a synonym, its genus; a homonym's lineage by NCBI's; organelle and placeholder records left out); an
18S's introns on its records; a genome's rRNA saved and read back (and as saved before format 4); the records' gaps,
contig ends and confidence; the copies' confidence.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import labels, ncbi as ncbi_mod, package, rrna  # noqa: E402
import fixtures as fx  # noqa: E402

PR2 = [
    "A1\tEukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces;"
    "Saccharomyces_cerevisiae;",
    "A2\tEukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces;"
    "Saccharomyces_cerevisiae;",
    # K. lactis by its anamorph's name, as PR2 may name it
    "A3\tEukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Kluyveromyces;"
    "Candida_sphaerica;",
    # a Candida not named to its species: its genus alone
    "A4\tEukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;Saccharomycetales;Debaryomycetaceae;Candida;Candida_sp.;",
    # a homonym: Mysterium filed twice, the more records under another kingdom
    "A5\tEukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Insecta;Insecta_X;Mysterium;Mysterium_sp.;",
    "A6\tEukaryota;Obazoa;Opisthokonta;Metazoa;Arthropoda;Insecta;Insecta_X;Mysterium;Mysterium_sp.;",
    "A7\tEukaryota;Excavata;Mysteriozoa;Mysteriozoa_X;Mysteriozoa_XX;Mysteriozoa_XXX;Mysteriozoa_XXXX;Mysterium;"
    "Mysterium_sp.;",
    # a plastid's record and a placeholder: no species
    "A8\tEukaryota;Archaeplastida;Chlorophyta:plas;Chlorophyta_X:plas;x;x;x;Homo;Homo_sapiens;",
    "A9\tEukaryota;Obazoa;Opisthokonta;Metazoa;Chordata;Mammalia;Hominidae;Homo_X;Homo_X_sp.;",
]
KSGP = [
    "K1\tk__Obazoa;p__Fungi;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;"
    "s__Saccharomyces_cerevisiae",
    "K2\tk__Bacteria;p__Pseudomonadota;c__G;o__E;f__P;g__Kluyveromyces;s__Kluyveromyces lactis",  # a prokaryote's
    "K3\tk__Obazoa;p__Fungi;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Candida;"
    "s__GCA_000182965",  # a genome's 18S: names its genus, no species
    "K4\tk__Obazoa;p__Metazoa;c__Mammalia;o__Primates;f__Hominidae;g__Homo;s__Homo_sapiens:mito",
]
EUKARYOME = [("superkingdom", "kingdom", "subreg", "phylum", "subphy", "class", "order", "family", "gen",
              "kgd-specific"),
             ("Obazoa", "Fungi", "Dikarya", "Ascomycota", "Saccharomycotina", "Saccharomycetes", "Saccharomycetales",
              "Saccharomycetaceae", "Kluyveromyces", ""),
             ("Obazoa", "Metazoa", "", "Nematoda", "", "Chromadorea", "Ascaridida", "Ascarididae", "Ascaris", "")]
LONG_READS = [("EUK1", "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;"
                       "g__Kluyveromyces;s__lactis"),  # EUKARYOME 2.1's s__: the epithet alone
              ("EUK2", "k__Fungi;p__Ascomycota;c__S;o__S;f__S;g__Saccharomyces;s__unclassified")]


def lines(path, rows):
    with open(path, "w") as fh:
        fh.writelines(r + "\n" for r in rows)
    return path


class LabelTest(unittest.TestCase):
    def test_species_named_in_the_three_references(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            pr2 = labels.read_pr2(lines(os.path.join(tmp, "pr2.tax"), PR2))
            ksgp = labels.read_ksgp(lines(os.path.join(tmp, "KSGP_v4.4_plus2.tax"), KSGP))
            euk = labels.read_eukaryome(EUKARYOME, LONG_READS)
            got = labels.label_species(ncbi, [4932, 28985, 5476, 952, 9606, 6253], [pr2, euk, ksgp])
            self.assertEqual(got[4932]["PR2"], ("Eukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;Saccharomycetales;"
                                                "Saccharomycetaceae;Saccharomyces;Saccharomyces_cerevisiae", "species"))
            self.assertEqual(got[28985]["PR2"][1], "synonym")  # by NCBI's synonym Candida sphaerica
            self.assertEqual(got[5476]["PR2"], ("Eukaryota;Obazoa;Opisthokonta;Fungi;Saccharomycetes;"
                                                "Saccharomycetales;Debaryomycetaceae;Candida", "genus"))
            # the homonym: the lineage that shares NCBI's phylum Mysteriozoa, not the one of more records
            self.assertEqual(got[952]["PR2"][0].split(";")[1:3], ["Excavata", "Mysteriozoa"])
            self.assertEqual(got[9606]["PR2"], ("-", "-"))  # its plastid record and its Homo_X placeholder name none
            # KSGP: its own prefixes; a prokaryote's record and an accession-named one name no species
            self.assertEqual(got[4932]["KSGP"], ("k__Obazoa;p__Fungi;c__Saccharomycetes;o__Saccharomycetales;"
                                                 "f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae",
                                                 "species"))
            self.assertEqual(got[28985]["KSGP"], ("-", "-"))
            self.assertEqual(got[5476]["KSGP"][1], "genus")
            self.assertEqual(got[9606]["KSGP"], ("-", "-"))  # a mitochondrion's record
            self.assertEqual(ksgp.version, "4.4")
            # EUKARYOME: a long read's species with its genus's lineage of the classification; a genus of it alone
            self.assertEqual(got[28985]["EUKARYOME"], ("Obazoa;Fungi;Dikarya;Ascomycota;Saccharomycotina;"
                                                       "Saccharomycetes;Saccharomycetales;Saccharomycetaceae;"
                                                       "Kluyveromyces;Kluyveromyces_lactis", "species"))
            self.assertEqual(got[6253]["EUKARYOME"][1], "genus")
            self.assertEqual(got[4932]["EUKARYOME"], ("-", "-"))  # s__unclassified names no species
            # format 5: two levels (data/supergroups.tsv), the most specific clade's row
            self.assertEqual((labels.supergroup(ncbi, 4932), labels.group(ncbi, 4932)), ("Obazoa", "Fungi"))
            self.assertEqual((labels.supergroup(ncbi, 6253), labels.group(ncbi, 6253)), ("Obazoa", "Metazoa"))
            self.assertEqual((labels.supergroup(ncbi, 952), labels.group(ncbi, 952)),
                             ("incertae sedis", "other eukaryotes"))
            # NCBI's empty ranks filled: an unranked clade where one lies between, else the name above with _X
            self.assertEqual(labels.ncbi_lineage_filled(ncbi, 6253),
                             "Eukaryota;Metazoa;Nematoda;Nematoda_X;Nematoda_XX;Nematoda_XXX;Ascaris;Ascaris lumbricoides")
            self.assertEqual(labels.ncbi_lineage_filled(ncbi, 4932), labels.ncbi_lineage(ncbi, 4932))  # none empty
            from eukore import groups
            self.assertTrue(groups.check(ncbi, lambda _t: None))  # the stand-in taxonomy lacks most clades

    def test_names(self):
        self.assertTrue(labels.is_species_name("Saccharomyces_cerevisiae"))
        for name in ("Saccharomyces_sp.", "Homo_X_sp.", "GCA_000182965", "Neurospora_?", "Fungi_X_XX",
                     "Blastocystis_sp._subtype_1", "unclassified", "Candida cf. albicans"):
            self.assertFalse(labels.is_species_name(name), name)
        self.assertEqual(labels.key("Lankesterella(Viridiplantae)"), "lankesterella")
        self.assertFalse(labels.is_genus_name("Boletales.fam.incertae.sedis"))


class RrnaTest(unittest.TestCase):
    def test_an_18s_cut_by_an_intron(self):
        rng = random.Random(3)
        seq = "".join(rng.choice("ACGT") for _ in range(20000))
        feats = {"18S_rRNA": [(1001, 1700, "+", False), (1901, 2800, "+", False)],
                 "5_8S_rRNA": [(3201, 3350, "+", False)], "28S_rRNA": [(3601, 9000, "+", False)]}
        item = rrna.units_of(seq, feats)[0]
        self.assertEqual(item["introns"], [("18S", 700, 900)])  # on the span: its bases 701-900
        self.assertEqual(item["unit_introns"], [("18S", 700, 900)])  # the unit starts at the 18S
        self.assertEqual(item["exons"], seq[1000:1700] + seq[1900:2800])
        n = len(seq)
        flip = {k: [(n - b + 1, n - a + 1, "-", p) for a, b, _s, p in v] for k, v in feats.items()}
        back = rrna.units_of(rrna.revcomp(seq), flip)[0]
        self.assertEqual((back["introns"], back["unit_introns"]), ([("18S", 700, 900)], [("18S", 700, 900)]))

    def test_saved_records_and_their_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = [dict(kind="rDNA", contig="c1", start=501, end=6500, strand="+", regions=[("18S", 0, 1800)],
                         lsu_complete=True, sequence="ACGT" * 10 + "N" * 12 + "ACGT", contig_length=7000,
                         introns=[("18S", 700, 900)]),
                    dict(kind="18S", contig="c2", start=1, end=1800, strand="-", regions=None, lsu_complete=False,
                         sequence="ACGT" * 450, contig_length=90000, introns=[])]
            self.assertEqual(rrna.save(tmp, rows), {"18S": 1, "rDNA": 1, "rDNA_complete": 1, "rDNA_contigs": 1})
            back = rrna.load(tmp)
            self.assertEqual((back[0]["contig_length"], back[0]["introns"], back[1]["introns"]),
                             (7000, [("18S", 700, 900)], None))
            # one saved before format 4: no contig length, no introns
            with open(os.path.join(tmp, "rrna.tsv"), "w") as fh:
                fh.write("kind\tcontig\tstart\tend\tstrand\tregions\tlsu_complete\tsequence\n"
                         "18S\tc1\t1\t1800\t+\t-\tf\tACGT\n")
            old = rrna.load(tmp)[0]
            self.assertEqual((old["contig_length"], old["introns"], old["sequence"]), (None, None, "ACGT"))
            records = rrna.genome_records("GCA_1", back, 4932, 4932)
            done = {n: r for n, _s, r in rrna.finish(records, None)}
            unit = done["GCA_1_rDNA_1"]
            self.assertEqual((unit["contig"], unit["start"], unit["contig_edge"], unit["introns"], unit["gaps"]),
                             ("c1", 501, 500, "18S:700-900", 1))
            self.assertEqual(unit["confidence"], "low")  # 12 Ns of 56 bases: ambiguous
            self.assertIn("contig_edge", unit["flags"].split(","))
            self.assertEqual(rrna.confidence(dict(kind="rDNA", own_species="t", screen="-", ambiguous="0", gaps=0,
                                                  homopolymer=5, contig_edge=5000, lsu_complete="t")),
                             ("high", "-"))
            self.assertEqual(rrna.confidence(dict(kind="18S", own_species="t", screen="18S:phylum", ambiguous="0",
                                                  gaps=0, homopolymer=5, contig_edge="-"))[0], "low")


class CopyTest(unittest.TestCase):
    def test_a_copys_confidence(self):
        base = dict(status="Single", own_species="t", fcs="-", window_class="eukarya", length_check="ok", gaps=0,
                    contig_edge=5000)
        self.assertEqual(package.copy_confidence(base), ("high", "-"))
        self.assertEqual(package.copy_confidence(dict(base, fcs="EXCLUDE:prok:Escherichia")), ("low", "fcs"))
        self.assertEqual(package.copy_confidence(dict(base, status="Fragmented", contig_edge=10)),
                         ("medium", "fragmented,contig_edge"))
        self.assertEqual(package.copy_confidence(dict(base, rescued="t")), ("medium", "rescued"))
        self.assertEqual(package.copy_confidence(dict(base, rescued="f")), ("high", "-"))
        self.assertEqual(package.gap_runs("ACGT" + "N" * 10 + "AC" + "N" * 9 + "A"), [(4, 14)])

    def test_call_rates_of_the_genomes_searched(self):
        # a marker added after one genome's run: that genome counts for it nowhere (until 2026-11-08 as not calling
        # it, its call rate halved)
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        genomes = {"G1": dict(representative="t", species_taxid=4932, _searched={"busco1", "marker1"}),
                   "G2": dict(representative="t", species_taxid=28985, _searched={"busco1"})}
        rows = [dict(accession=a, locus=lo, status="Single", length="100", cds_length="90", protein_length="30")
                for a, lo in (("G1", "busco1"), ("G2", "busco1"), ("G1", "marker1"))]
        got, columns = package.loci_rows([dict(locus="busco1"), dict(locus="marker1")],
                                         dict(genomes=genomes, copy_rows=rows), ncbi)
        self.assertEqual([(r["species"], r["called"], r["call_rate"]) for r in got], [(2, 2, "1.000"), (1, 1, "1.000")])
        self.assertEqual({r[columns[0]] for r in got}, {"1.000"})


if __name__ == "__main__":
    unittest.main()
