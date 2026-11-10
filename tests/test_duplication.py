#!/usr/bin/env python3
"""Tests of the 2026-11-10 audit's B3-B5 (eukore/package.py, eukore/release.py): a genome's duplication class where
its Duplicated BUSCOs' copies lie on short contigs (`haplotigs`: an unpurged diploid, not a polyploid; the
simulation's core cases, ~/eukore-sim/duplication); a Duplicated copy's tier in a genome duplicated genome-wide; the
package's `unverified` copies of lone lineages counted (the manifest's counts.copy_unverified, stats.json's
unverified_alone); each BUSCO's call rate over the genomes alone, the transcriptomes' apart.

    python3 -m unittest discover -s tests -t tests
"""

import json
import os
import random
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi as ncbi_mod, package, release  # noqa: E402
import fixtures as fx  # noqa: E402
from test_release import toy_package  # noqa: E402

LONG, SHORT = 1_000_000, 60_000


def genome(single, dups, identity, chosen=LONG, other=LONG):
    """A genome's BUSCO rows and its alternates: `single` Single copies on contigs around LONG; `dups` Duplicated
    ones `identity` alike, the copy chosen on a contig of `chosen` bases (a function of the locus' number, or a
    length), its other copy on one of `other`."""
    def of(length, i):
        return str(length(i) if callable(length) else length)
    rows = [dict(locus=f"s{i}", status="Single", contig_length=str(LONG + 10_000 * (i % 7))) for i in range(single)]
    rows += [dict(locus=f"d{i}", status="Duplicated", dup_identity=identity, contig_length=of(chosen, i))
             for i in range(dups)]
    return rows, [dict(locus=f"d{i}", contig_length=of(other, i)) for i in range(dups)]


class DuplicationTest(unittest.TestCase):
    def test_an_unpurged_diploid_is_not_a_polyploid(self):
        # 2% heterozygosity: the alleles 98% alike, under HAPLOTYPE_ID (until 2026-11-10 wgd_or_polyploid); the copy
        # chosen is the haplotig's in half of the loci
        rows, alts = genome(60, 40, "98.0", chosen=lambda i: SHORT if i % 2 else LONG,
                            other=lambda i: LONG if i % 2 else SHORT)
        got = package.duplication_of(rows)
        self.assertEqual((got["duplicated_share"], got["dup_identity_median"], got["duplication"],
                          got["dup_short_contigs"]), ("0.400", "98.0", "haplotigs", "0.50/0.00"))
        # with the other copies' contigs (alternates.tsv's rows) every locus shows its haplotig
        got = package.duplication_of(rows, alts)
        self.assertEqual((got["duplication"], got["dup_short_contigs"]), ("haplotigs", "1.00/0.00"))
        # the chosen copies all on the long contigs: the rows alone do not tell, the alternates do
        rows, alts = genome(60, 40, "98.0", chosen=LONG, other=SHORT)
        self.assertEqual(package.duplication_of(rows)["duplication"], "wgd_or_polyploid")
        self.assertEqual(package.duplication_of(rows, alts)["duplication"], "haplotigs")

    def test_a_polyploid_stays_one(self):
        # a recent polyploid's copies are as alike as alleles, but lie on chromosomes like any gene
        rows, alts = genome(60, 40, "98.0")
        for got in (package.duplication_of(rows), package.duplication_of(rows, alts)):
            self.assertEqual((got["duplication"], got["dup_short_contigs"]), ("wgd_or_polyploid", "0.00/0.00"))
        # two haplotypes assembled equally long are not told from it: the class says so (docs/FORMAT.md)
        # a contaminant's copies lie on short contigs too, but are no alleles (under HAPLOTIG_ID alike)
        rows, alts = genome(60, 40, "82.0", chosen=LONG, other=SHORT)
        self.assertEqual(package.duplication_of(rows, alts)["duplication"], "wgd_or_polyploid")
        # 99.5% alike or more: haplotypes, wherever they lie (as before)
        rows, alts = genome(60, 40, "99.8", chosen=LONG, other=SHORT)
        self.assertEqual(package.duplication_of(rows, alts)["duplication"], "haplotypes")
        self.assertEqual(package.duplication_of(*genome(60, 40, "99.8"))["duplication"], "haplotypes")

    def test_too_little_to_tell(self):
        # under SHORT_SINGLE_MIN Single copies: what a short contig is is not known, the class by identity alone
        rows, alts = genome(5, 40, "98.0", chosen=LONG, other=SHORT)
        got = package.duplication_of(rows, alts)
        self.assertEqual((got["duplication"], got["dup_short_contigs"]), ("wgd_or_polyploid", "-"))
        # rows without contig lengths (a test's, an old run's): as until 2026-11-10
        rows = [dict(status="Single")] * 60 + [dict(status="Duplicated", dup_identity="98.0")] * 40
        self.assertEqual(package.duplication_of(rows)["duplication"], "wgd_or_polyploid")
        self.assertEqual(package.duplication_of(rows[:10])["dup_short_contigs"], "-")  # too few called
        # a fragmented assembly, every contig short: short is by the Single copies' own contigs, none stands out
        rows, alts = genome(60, 40, "98.0")
        for r in rows + alts:
            r["contig_length"] = str(int(r["contig_length"]) // 50)
        self.assertEqual(package.duplication_of(rows, alts)["duplication"], "wgd_or_polyploid")

    def test_the_simulations_core_cases(self):
        # ~/eukore-sim/duplication/sim_dup.py: 120 BUSCOs, a locus' contig a lognormal draw around the assembly's
        # typical contig. An unpurged diploid (60% assembled twice, the second copy on haplotigs a tenth as long, the
        # copy chosen the haplotig's in half) against a recent polyploid (60% retained twice, 1-3% apart)
        def simulated(rng, second, identity):
            rows, alts = [], []
            for i in range(120):
                first = int(1e6 * rng.lognormvariate(0, 1.0))
                if rng.random() >= 0.6:
                    rows.append(dict(locus=str(i), status="Single", contig_length=str(first)))
                    continue
                other = int(second * rng.lognormvariate(0, 1.0))
                if rng.random() < 0.5:
                    first, other = other, first
                rows.append(dict(locus=str(i), status="Duplicated", contig_length=str(first),
                                 dup_identity=f"{identity(rng):.1f}"))
                alts.append(dict(locus=str(i), contig_length=str(other)))
            return rows, alts
        rng = random.Random(3)
        diploid = [simulated(rng, 1e5, lambda r: r.uniform(97.5, 98.5)) for _ in range(40)]
        polyploid = [simulated(rng, 1e6, lambda r: r.uniform(97.0, 99.0)) for _ in range(40)]
        with_alts = [package.duplication_of(rows, alts)["duplication"] for rows, alts in diploid]
        alone = [package.duplication_of(rows)["duplication"] for rows, _alts in diploid]
        self.assertGreaterEqual(with_alts.count("haplotigs"), 38)
        self.assertGreaterEqual(alone.count("haplotigs"), 32)  # the rows alone: half of the haplotigs seen
        self.assertEqual(set(with_alts) | set(alone), {"haplotigs", "wgd_or_polyploid"})
        calls = [package.duplication_of(rows, alts)["duplication"] for rows, alts in polyploid] + \
            [package.duplication_of(rows)["duplication"] for rows, _alts in polyploid]
        self.assertGreaterEqual(calls.count("wgd_or_polyploid"), 78)

    def test_a_duplicated_copy_in_a_duplicated_genome(self):
        r = dict(status="Duplicated", own_species="t", fcs="-", window_class="eukarya", length_check="ok", gaps=0,
                 contig_edge=5000)
        self.assertEqual(package.copy_confidence(r), ("medium", "duplicated"))
        # a polyploid's, an unpurged diploid's: two copies are the genome's truth (the caller says so)
        self.assertEqual(package.copy_confidence(r, duplicated_genome=True), ("high", "-"))
        self.assertEqual(package.copy_confidence(dict(r, gaps=1), duplicated_genome=True), ("medium", "gap"))
        self.assertEqual(package.copy_confidence(dict(r, status="Fragmented"), duplicated_genome=True),
                         ("medium", "fragmented"))


class CountsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(self.tmp.name, "taxdump.tar.gz")))

    def tearDown(self):
        self.tmp.cleanup()

    def test_call_rates_over_the_genomes_alone(self):
        # two yeasts' genomes, a Candida's transcriptome (source tsa), EukDetect2's records of a nematode
        genomes = {"G1": dict(representative="t", species_taxid=4932, source="ncbi"),
                   "G2": dict(representative="t", species_taxid=28985, source="model"),
                   "T1": dict(representative="t", species_taxid=5476, source="tsa"),
                   "E1": dict(representative="t", species_taxid=6253, source="no_assembly"),
                   "A1": dict(representative="f", species_taxid=4932, source="assembly")}
        for g in genomes.values():
            g["_searched"] = {"b1", "b2"}

        def row(acc, locus, status="Single"):
            return dict(accession=acc, locus=locus, status=status, length="100", cds_length="90", protein_length="30")
        rows = [row("G1", "b1"), row("G2", "b1"), row("T1", "b1", "Fragmented"), row("E1", "b1"),
                row("G1", "b2"), row("T1", "b2")]
        got, columns = package.loci_rows([dict(locus="b1"), dict(locus="b2")], dict(genomes=genomes, copy_rows=rows),
                                         self.ncbi)
        b1, b2 = got
        # the columns until 2026-11-10 as they were: every record counts
        self.assertEqual((b1["species"], b1["called"], b1["fragmented"], b1["call_rate"], b1["call_Fungi"],
                          b1["call_Metazoa"]), (4, 3, 1, "0.750", "0.667", "1.000"))
        self.assertEqual((b2["called"], b2["call_rate"], b2["call_Fungi"]), (2, "0.500", "0.667"))
        # the genomes alone: the transcriptome's fragment is no gene lost
        self.assertEqual((b1["genomes"], b1["genome_call_rate"], b1["genome_call_Fungi"], b1["transcriptomes"],
                          b1["transcriptome_call_rate"]), (2, "1.000", "1.000", 1, "0.000"))
        self.assertEqual((b2["genome_call_rate"], b2["genome_call_Fungi"], b2["transcriptome_call_rate"]),
                         ("0.500", "0.500", "1.000"))
        self.assertEqual(columns, ["call_Fungi", "call_Metazoa", "genomes", "genome_call_rate", "transcriptomes",
                                   "transcriptome_call_rate", "genome_call_Fungi"])  # no genome among the Metazoa
        self.assertFalse([c for c in columns[2:] if c.startswith("call_")])  # a reader's groups: call_<group> alone
        # a package of genomes alone: its rates the old ones, the transcriptomes' not rated
        only = {a: g for a, g in genomes.items() if a in ("G1", "G2")}
        collected = dict(genomes=only, copy_rows=[r for r in rows if r["accession"] in only])
        (b1, _b2), _columns = package.loci_rows([dict(locus="b1"), dict(locus="b2")], collected, self.ncbi)
        self.assertEqual((b1["call_rate"], b1["genome_call_rate"], b1["transcriptomes"], b1["transcriptome_call_rate"]),
                         ("1.000", "1.000", 0, "-"))

    def copies(self):
        """Two Saccharomycetaceae (a family), a Candida (alone in its family), two Mysterium (NCBI names the genus
        no family): their copies of two BUSCOs and a marker."""
        genomes = {a: dict(representative="t", species_taxid=sp) for a, sp in (
            ("G1", 4932), ("G2", 28985), ("G3", 5476), ("G4", 952), ("G5", 953))}
        rows = [dict(accession=a, locus=lo, status="Single", confidence=c) for a, lo, c in (
            ("G1", "b1", "unverified"), ("G2", "b1", "high"), ("G3", "b1", "unverified"), ("G4", "b1", "unverified"),
            ("G5", "b1", "unverified"), ("G1", "b2", "unverified"), ("G3", "b2", "medium"), ("G3", "m1", "unverified"))]
        return genomes, rows

    def test_the_unverified_of_lone_lineages(self):
        genomes, rows = self.copies()
        got = package.unverified_alone(rows, genomes, self.ncbi, {"b1": "busco", "b2": "busco", "m1": "markers"})
        # b1: the Candida is alone in its family; the yeast has its family's Kluyveromyces, the two Mysterium each
        # other (their genus). b2: the yeast is its family's only holder of it. m1: the Candida alone
        self.assertEqual(got, {"busco": dict(copies=5, alone_in_family=2),
                               "markers": dict(copies=1, alone_in_family=1)})
        self.assertEqual(package.unverified_alone([r for r in rows if r["confidence"] != "unverified"], genomes,
                                                  self.ncbi), {})

    def test_the_count_in_the_manifest_and_the_stats(self):
        genomes, rows = self.copies()
        for acc, g in genomes.items():
            g.update(accession=acc, organism="x", source="ncbi", loci=2, _searched={"b1", "b2", "m1"})
        for r in rows:
            r.update(copies="1", length="100", cds_length="90", protein_length="30")
        loci = [dict(locus="b1"), dict(locus="b2"), dict(locus="m1", set="markers")]
        collected = dict(genomes=genomes, copy_rows=rows, taxids=set(), cds={}, proteins={}, work=None, done=[],
                         copy_count=len(rows))
        out = package.write(os.path.join(self.tmp.name, "Eukore_x"), "x", self.ncbi, collected, loci)
        with open(os.path.join(out, package.MANIFEST)) as fh:
            counts = json.load(fh)["counts"]
        self.assertEqual(counts["copy_confidence"], {"unverified": 6, "high": 1, "medium": 1})
        self.assertEqual(counts["copy_unverified"], {"busco": {"copies": 5, "alone_in_family": 2},
                                                     "markers": {"copies": 1, "alone_in_family": 1}})
        # busco/loci.tsv with the genomes' call rates after the groups' (B5)
        head = list(package.read_table(out, "busco/loci.tsv")[0])
        self.assertEqual(head[:len(package.LOCI_COLUMNS)], list(package.LOCI_COLUMNS))
        self.assertIn("genome_call_rate", head)
        # stats.json's locus categories carry the count; a package written before has none to carry
        toy = toy_package(os.path.join(self.tmp.name, "toy"))
        stats, _tree = release.stats(toy)
        self.assertEqual({c["key"]: c.get("unverified_alone", "no key") for c in stats["markers"]},
                         {"busco": None, "markers": None, "rrna": "no key", "mitochondrion": "no key",
                          "plastid": "no key"})
        with open(os.path.join(toy, package.MANIFEST)) as fh:
            manifest = json.load(fh)
        manifest["counts"]["copy_unverified"] = {"busco": {"copies": 5, "alone_in_family": 2}}
        with open(os.path.join(toy, package.MANIFEST), "w") as fh:
            json.dump(manifest, fh)
        stats, _tree = release.stats(toy)
        self.assertEqual({c["key"]: c["unverified_alone"] for c in stats["markers"] if "unverified_alone" in c},
                         {"busco": 2, "markers": 0})


if __name__ == "__main__":
    unittest.main()
