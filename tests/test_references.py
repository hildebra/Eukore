#!/usr/bin/env python3
"""Tests of the references added in format 5 (eukore/labels.py; since 2026-11-07): SILVA's paths of any depth, UNITE's
k__ lineages and SH codes, MIDORI2's own NCBI taxids, and crossrefs.tsv (Open Tree's ids by NCBI taxid, UNITE's SH
codes and BOLD's BINs by name).

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import io
import os
import sys
import tarfile
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import labels, ncbi as ncbi_mod  # noqa: E402
import fixtures as fx  # noqa: E402

UNITE = [">Saccharomyces_cerevisiae|AB000001|SH0001.10FU|reps|k__Fungi;p__Ascomycota;c__Saccharomycetes;"
         "o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae\n",
         ">Saccharomyces_cerevisiae|AB000002|SH0002.10FU|refs|k__Fungi;p__Ascomycota;c__Saccharomycetes;"
         "o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae\n", "ACGT\n"]


def tgz(path, name, text):
    with tarfile.open(path, "w:gz") as tar:
        data = text.encode()
        info = tarfile.TarInfo(name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    return path


class ReferenceTest(unittest.TestCase):
    def test_silva_unite_midori2(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            silva = os.path.join(tmp, "silva.fasta.gz")
            with gzip.open(silva, "wt") as fh:
                fh.write(">AB1.1.1800 Eukaryota;Amorphea;Obazoa;Opisthokonta;Nucletmycea;Fungi;Dikarya;Ascomycota;"
                         "Saccharomycotina;Saccharomycetes;Saccharomyces;Saccharomyces cerevisiae (baker's yeast)\nACGU\n"
                         ">AB2.1.1800 Eukaryota;Amorphea;Obazoa;Opisthokonta;Holozoa;Metazoa;Nematoda;Ascaris;"
                         "Ascaris sp. KX1\nACGU\n"
                         ">CP1.1.1500 Bacteria;Cyanobacteria;Chloroplast;Zea mays\nACGU\n")
            ref = labels.read_silva(silva)
            lineage, how = ref.label(ncbi, 4932)
            self.assertEqual((how, ref.text(lineage).split(";")[-2:]), ("species", ["Saccharomyces",
                                                                                   "Saccharomyces_cerevisiae"]))
            self.assertEqual(ref.label(ncbi, 6253)[1], "genus")  # "Ascaris sp." names its genus alone
            unite = labels.read_unite(UNITE)
            lineage, how = unite.label(ncbi, 4932)
            self.assertEqual((how, unite.text(lineage)[:12]), ("species", "k__Fungi;p__"))
            self.assertEqual(unite.ids_of(ncbi, 4932), ["SH0001.10FU", "SH0002.10FU"])
            midori = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB264_CO1_RAW.fasta")
            with open(midori, "w") as fh:
                fh.write(">MK1.1.1.658###root_1;Eukaryota_2759;Nematoda_6231;Chromadorea_119089;Ascaridida_6249;"
                         "Ascarididae_6250;Ascaris_6251;Ascaris lumbricoides_6253\nACGT\n")
            ref = labels.read_midori2(midori)
            self.assertEqual(ref.version, "GB264")
            lineage, how = ref.label(ncbi, 6253)
            self.assertEqual((how, lineage[-1]), ("taxid", "Ascaris_lumbricoides"))

    def test_crossrefs(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            ott = tgz(os.path.join(tmp, "ott3.7.3.tgz"), "ott3.7.3/taxonomy.tsv",
                      "uid\t|\tparent_uid\t|\tname\t|\trank\t|\tsourceinfo\t|\tuniqname\t|\tflags\t|\t\n"
                      "356221\t|\t908546\t|\tSaccharomyces cerevisiae\t|\tspecies\t|\t"
                      "ncbi:4932,worms:147593,gbif:2599166,irmng:10384091\t|\t\t|\t\t|\t\n"
                      "1\t|\t\t|\tHomo sapiens\t|\tspecies\t|\tncbi:9606,gbif:2436436\t|\t\t|\t\t|\t\n")
            bold = os.path.join(tmp, "bold.tsv")
            with open(bold, "w") as fh:
                fh.write("processid\tbin_uri\tspecies\nP1\tBOLD:AAA0001\tAscaris lumbricoides\n"
                         "P2\tNone\tAscaris lumbricoides\n")
            rows = labels.crossrefs(ncbi, {4932, 6253}, [labels.read_unite(UNITE)], ott, bold, log=lambda *_: None)
            got = {(r["taxid"], r["reference"]): r["id"] for r in rows}
            self.assertEqual(got[(4932, "OTT")], "356221")
            self.assertEqual((got[(4932, "GBIF")], got[(4932, "WoRMS")], got[(4932, "IRMNG")]),
                             ("2599166", "147593", "10384091"))
            self.assertEqual(got[(6253, "BOLD")], "BOLD:AAA0001")
            self.assertEqual(sorted(r["id"] for r in rows if r["reference"] == "UNITE"), ["SH0001.10FU", "SH0002.10FU"])
            self.assertNotIn((9606, "OTT"), got)  # not a species of the package


if __name__ == "__main__":
    unittest.main()
