#!/usr/bin/env python3
"""Tests of barrnap's input (rrna.barrnap_fasta; since 2026-11-12): nhmmer guesses its target's alphabet from the start
of the FASTA's first sequence and stops where a base is missing there. Every genome the backfill of 2026-11-11 tried
failed so: Eupithecia's first window began at a chromosome's end, 4,000 bases of the telomere's (TTAGG)n with no G.
The window whose start holds every base now goes first (else a guard sequence), and the guard's hits are dropped.

    python3 -m unittest discover -s tests -t tests
"""

import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import organelle, rrna  # noqa: E402
import fixtures as fx  # noqa: E402

# barrnap's stand-in: fails as nhmmer does where the first sequence's first 4,000 bases lack one of A, C, G and T;
# else a GFF of one 18S on every sequence named in HITS (the guard's too, which barrnap_gff drops)
BARRNAP = r'''#!/usr/bin/env python3
import sys
path = sys.argv[-1]
names, first = [], None
with open(path) as fh:
    for line in fh:
        if line.startswith(">"):
            names.append(line[1:].split()[0])
        elif first is None:
            first = line.strip()
if first is not None and any(b not in first[:4000].upper() for b in "ACGT"):
    sys.stderr.write("[barrnap] ERROR: nhmmer failed to run - Error: Unable to guess alphabet for target sequence "
                     f"database file {path}\n")
    sys.exit(2)
print("##gff-version 3")
for n in names:
    print(f"{n}\tbarrnap:0.9\trRNA\t1\t100\t1e-50\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA (partial);"
          "note=aligned only 5 percent of the 18S ribosomal RNA")
'''

TELOMERE = "CCTAA"  # the insects' (TTAGG)n on the reverse strand
G_FREE_SEED = "CAATTCCTTTAA"  # a reverse seed of the 18S with no G (rrna.SEEDS)


def dna(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


class AlphabetFirstTest(unittest.TestCase):
    def test_the_window_with_every_base_first(self):
        telomeric = ("w0", TELOMERE * 1000)
        records = [telomeric, ("w1", dna(5000, 1)), ("w2", dna(5000, 2))]
        self.assertEqual([n for n, _s in rrna.alphabet_first(records)], ["w1", "w0", "w2"])
        # a good first sequence stays where it is; a single G in 4,000 bases is not enough
        self.assertEqual([n for n, _s in rrna.alphabet_first(records[1:] + [telomeric])], ["w1", "w2", "w0"])
        some_g = ("w3", TELOMERE * 700 + "G" + TELOMERE * 300)
        self.assertEqual([n for n, _s in rrna.alphabet_first([some_g, ("w1", dna(5000, 1))])], ["w1", "w3"])
        # none holds every base: the guard first, a sequence of all four
        got = rrna.alphabet_first([telomeric, some_g, ("w4", "A" * 9000)])
        self.assertEqual([n for n, _s in got], [rrna.ALPHABET_GUARD, "w0", "w3", "w4"])
        self.assertTrue(rrna.guessable(got[0][1]))
        self.assertEqual(rrna.alphabet_first([]), [])

    def test_barrnap_gets_a_guessable_file_and_the_guards_hits_are_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            barrnap = fx.tool(tmp, "barrnap", BARRNAP)
            naive = os.path.join(tmp, "naive.fna")
            records = [("w0", TELOMERE * 1000), ("w1", dna(5000, 1))]
            with open(naive, "w") as fh:  # as written until 2026-11-12: the contigs' order
                fh.writelines(f">{n}\n{s}\n" for n, s in records)
            with self.assertRaisesRegex(rrna.BarrnapError, "Unable to guess alphabet"):
                rrna.barrnap_gff(barrnap, naive)
            fasta = os.path.join(tmp, "w.fna")
            rrna.barrnap_fasta(fasta, records)
            self.assertEqual(sorted(rrna.parse_gff(rrna.barrnap_gff(barrnap, fasta))), ["w0", "w1"])
            rrna.barrnap_fasta(fasta, [("w0", TELOMERE * 1000)])  # with the guard
            gff = rrna.barrnap_gff(barrnap, fasta)
            self.assertNotIn(rrna.ALPHABET_GUARD, gff)
            self.assertEqual(sorted(rrna.parse_gff(gff)), ["w0"])

    def test_a_genome_whose_first_window_is_a_telomere(self):
        # chr1's 18S seed sits in its telomere's repeat: its window, the genome's first, begins with 4,000 bases of no G
        chr1 = TELOMERE * 600 + G_FREE_SEED + TELOMERE * 600 + dna(3000, 3)
        chr2 = dna(4000, 4) + "AAACTTAAAGGA" + dna(20000, 5)
        with tempfile.TemporaryDirectory() as tmp:
            barrnap = fx.tool(tmp, "barrnap", BARRNAP)
            genome = os.path.join(tmp, "g.fna")
            with open(genome, "w") as fh:
                fh.write(f">chr1 telomeric\n{chr1}\n>chr2\n{chr2}\n")
            self.assertEqual(rrna.windows(chr1)[0][0], 0)  # the window from the chromosome's first base
            # until 2026-11-12: BarrnapError, Unable to guess alphabet; the stand-in's 100-base hits are under MIN_LEN
            self.assertEqual(rrna.extract_genome(genome, barrnap), [])

    def test_the_organelle_molecules_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            barrnap = fx.tool(tmp, "barrnap", BARRNAP)
            at_rich = ("mt1", "mitochondrion", "AT" * 3000)
            got = organelle.barrnap_genes([at_rich, ("mt2", "mitochondrion", dna(6000, 6))], barrnap, "mito")
            self.assertIsInstance(got, dict)


if __name__ == "__main__":
    unittest.main()
