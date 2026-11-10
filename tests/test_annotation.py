#!/usr/bin/env python3
"""Tests of the 2026-11-10 audit's fixes in eukore/annotation.py and eukore/twobit.py: the assembly report's
submitter's and UCSC names as aliases (C20), a GTF's stop_codon lines made coding (C21), a .2bit record of length 0
(C23). The modules' older tests are in test_paper.py and test_loci.py.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import struct
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import annotation, twobit  # noqa: E402

REPORT_HEAD = ("# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\t"
               "Relationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n")
REPORT = (
    # a scaffold its submitter named as chromosome 4 is named (its Assigned-Molecule): `4` stays the chromosome's
    "4\tunplaced-scaffold\tna\tna\tJA17.1\t<>\tna\tPrimary Assembly\t400\tchrUn_JA17v1\n"
    "chrIV\tassembled-molecule\t4\tChromosome\tCM4.1\t=\tNC_4.1\tPrimary Assembly\t5000\tchr4\n"
    "scaffold_12\tunplaced-scaffold\tna\tna\tJA12.1\t=\tNW_12.1\tPrimary Assembly\t900\tchrUn_JA12v1\n"
    "scaffold_13\tunplaced-scaffold\tna\tna\tJA13.1\t<>\tna\tPrimary Assembly\t800\tna\n"
    # one name, two sequences (two units' scaffolds named alike)
    "dup\tunplaced-scaffold\tna\tna\tJA14.1\t<>\tna\tPrimary Assembly\t700\tna\n"
    "dup\talt-scaffold\tna\tna\tJA15.1\t<>\tna\tALT_1\t600\tna\n"
    "na\tunplaced-scaffold\tna\tna\tJA16.1\t<>\tna\tPrimary Assembly\t500\tna\n"
    # named as another sequence's accession without its version
    "JA13\tunplaced-scaffold\tna\tna\tJA18.1\t<>\tna\tPrimary Assembly\t450\tna\n"
    "MT\tassembled-molecule\tMT\tMitochondrion\tna\t<>\tNC_9.1\tnon-nuclear\t300\tchrM\n")


def gtf_line(contig, kind, start, end, strand, phase, attrs):
    return f"{contig}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"


def twobit_file(path, seqs):
    """A .2bit file (version 0, little-endian) of [(name, sequence of ACGT)]: no N block, no masked block."""
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    bodies = []
    for _name, seq in seqs:
        padded = seq + "T" * (-len(seq) % 4)
        packed = bytes(sum(code[c] << (6 - 2 * i) for i, c in enumerate(padded[j:j + 4]))
                       for j in range(0, len(padded), 4))
        bodies.append(struct.pack("<IIII", len(seq), 0, 0, 0) + packed)  # length, N blocks, masked blocks, reserved
    head = struct.pack("<IIII", twobit.SIGNATURE, 0, len(seqs), 0)
    offset, index = len(head) + sum(1 + len(n) + 4 for n, _s in seqs), b""
    for (name, _s), body in zip(seqs, bodies):
        index += bytes([len(name)]) + name.encode() + struct.pack("<I", offset)
        offset += len(body)
    with open(path, "wb") as fh:
        fh.write(head + index + b"".join(bodies))


class AliasTest(unittest.TestCase):
    def test_the_submitters_and_ucsc_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = os.path.join(tmp, "report.txt")
            with open(report, "w") as fh:
                fh.write(REPORT_HEAD + REPORT)
            alias = annotation.aliases(report)
            self.assertEqual({n: alias.get(n) for n in ("chrIV", "chr4", "scaffold_12", "chrUn_JA12v1", "scaffold_13",
                                                        "chrUn_JA17v1", "MT", "chrM")},
                             {"chrIV": "CM4.1", "chr4": "CM4.1", "scaffold_12": "JA12.1", "chrUn_JA12v1": "JA12.1",
                              "scaffold_13": "JA13.1", "chrUn_JA17v1": "JA17.1", "MT": "NC_9.1", "chrM": "NC_9.1"})
            # a name of two sequences, and `na`, lead nowhere; a name that led elsewhere before still does
            self.assertNotIn("dup", alias)
            self.assertNotIn("na", alias)
            self.assertEqual((alias["4"], alias["JA13"]), ("CM4.1", "JA13.1"))
            # as before: the accessions' pairs, the unversioned accession
            self.assertEqual((alias["NC_4.1"], alias["CM4.1"], alias["JA12"]), ("CM4.1", "NC_4.1", "JA12.1"))
            # a gene set naming its scaffolds the submitter's way, on GenBank's sequences and on the RefSeq copy's
            gff = os.path.join(tmp, "a.gff")
            with open(gff, "w") as fh:
                fh.write(gtf_line("scaffold_12", "CDS", 101, 250, "+", 0, "ID=c1;Parent=t1")
                         + gtf_line("chr4", "CDS", 101, 250, "+", 0, "ID=c2;Parent=t2")
                         + gtf_line("dup", "CDS", 101, 250, "+", 0, "ID=c3;Parent=t3")
                         + gtf_line("na", "CDS", 101, 250, "+", 0, "ID=c4;Parent=t4"))
            spans = {n: [(1, 1000)] for n in ("JA12.1", "CM4.1", "JA14.1", "JA15.1", "JA16.1")}
            self.assertEqual({n: t["contig"] for n, t in annotation.read_cds(gff, spans, alias).items()},
                             {"t1": "JA12.1", "t2": "CM4.1"})
            self.assertEqual({n: t["contig"] for n, t in annotation.read_cds(
                gff, {"NW_12.1": [(1, 1000)], "NC_4.1": [(1, 1000)]}, alias).items()},
                {"t1": "NW_12.1", "t2": "NC_4.1"})  # two steps: scaffold_12 -> JA12.1 -> NW_12.1
            self.assertEqual(annotation.read_cds(gff, spans), {})  # without the report's names: none

    def test_a_report_without_those_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            report = os.path.join(tmp, "report.txt")
            with open(report, "w") as fh:
                fh.write("chr2\tassembled-molecule\t2\tChromosome\tCM2.1\t=\tNC_2.1\n")  # seven columns: no UCSC name
            alias = annotation.aliases(report)
            self.assertEqual((alias["2"], alias["chr2"], alias["CM2"]), ("CM2.1", "CM2.1", "CM2.1"))


class StopCodonTest(unittest.TestCase):
    def test_a_gtfs_stop_codon_is_coding(self):
        with tempfile.TemporaryDirectory() as tmp:
            gtf = os.path.join(tmp, "a.gtf.gz")
            with gzip.open(gtf, "wt") as fh:
                fh.write(
                    # the plus strand: the stop after the last CDS
                    gtf_line("c1", "CDS", 101, 200, "+", 0, 'gene_id "G1"; transcript_id "T1";')
                    + gtf_line("c1", "CDS", 301, 397, "+", 2, 'gene_id "G1"; transcript_id "T1";')
                    + gtf_line("c1", "stop_codon", 398, 400, "+", 0, 'gene_id "G1"; transcript_id "T1";')
                    # the minus strand: the stop before its first CDS on the sequence, its line before the CDS'
                    + gtf_line("c1", "stop_codon", 1001, 1003, "-", 0, 'gene_id "G2"; transcript_id "T2";')
                    + gtf_line("c1", "CDS", 1004, 1100, "-", 1, 'gene_id "G2"; transcript_id "T2";')
                    + gtf_line("c1", "CDS", 1201, 1300, "-", 0, 'gene_id "G2"; transcript_id "T2";')
                    # a stop an intron splits: one base ends the exon, two begin the next
                    + gtf_line("c1", "CDS", 2001, 2100, "+", 0, 'gene_id "G3"; transcript_id "T3";')
                    + gtf_line("c1", "stop_codon", 2101, 2101, "+", 0, 'gene_id "G3"; transcript_id "T3";')
                    + gtf_line("c1", "stop_codon", 2201, 2202, "+", 2, 'gene_id "G3"; transcript_id "T3";')
                    # a GFF3's CDS holds its stop: as it is
                    + gtf_line("c1", "CDS", 3001, 3100, "+", 0, "ID=cds-4;Parent=rna-4")
                    + gtf_line("c1", "stop_codon", 3098, 3100, "+", 0, "ID=stop-4;Parent=rna-4")
                    # no stop line (a partial model): as it is; a stop of no transcript, or off the spans: not read
                    + gtf_line("c1", "CDS", 4001, 4100, "+", 0, 'gene_id "G5"; transcript_id "T5";')
                    + gtf_line("c1", "stop_codon", 4101, 4103, "+", 0, 'gene_id "G5";')
                    + gtf_line("c9", "stop_codon", 4101, 4103, "+", 0, 'gene_id "G5"; transcript_id "T5";')
                    + gtf_line("c1", "stop_codon", 9001, 9003, "+", 0, 'gene_id "G6"; transcript_id "T6";'))
            got = annotation.read_cds(gtf, {"c1": [(1, 5000)]})
            self.assertEqual({n: t["cds"] for n, t in got.items()},
                             {"T1": [(101, 200), (301, 400)], "T2": [(1001, 1100), (1201, 1300)],
                              "T3": [(2001, 2101), (2201, 2202)], "rna-4": [(3001, 3100)], "T5": [(4001, 4100)]})
            # the interval a stop ends keeps its phase; a stop's own interval has its line's
            self.assertEqual(got["T1"]["phase"], {(101, 200): 0, (301, 400): 2})
            self.assertEqual(got["T2"]["phase"], {(1001, 1100): 1, (1201, 1300): 0})
            self.assertEqual(got["T3"]["phase"], {(2001, 2101): 0, (2201, 2202): 2})


class TwoBitTest(unittest.TestCase):
    def test_a_record_of_length_0(self):
        seqs = [("empty0", ""), ("a", "ACGTACGTAC"), ("empty", ""), ("b", "GGA")]
        with tempfile.TemporaryDirectory() as tmp:
            path, dest = os.path.join(tmp, "g.2bit"), os.path.join(tmp, "g.fna.gz")
            twobit_file(path, seqs)
            self.assertEqual([(n, bytes(s).decode()) for n, s in twobit.records(path)], seqs)
            self.assertEqual(twobit.to_fasta(path, dest), 13)
            with gzip.open(dest, "rt") as fh:
                self.assertEqual(fh.read(), ">empty0\n>a\nACGTACGTAC\n>empty\n>b\nGGA\n")


if __name__ == "__main__":
    unittest.main()
