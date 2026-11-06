#!/usr/bin/env python3
"""Tests of the benchmarks (eukore/benchmarks.py, since 2026-11-06): a copy translated by its genetic code and its
5-mers against OrthoDB's protein; the resolution table's ranks and shares; each species' best 18S; and a package
written with its benchmarks (the build's way: before its manifest), then measured again from its files alone.

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
from eukore import benchmarks, genomes, ncbi as ncbi_mod, odb as odb_mod, package, sources  # noqa: E402
import fixtures as fx  # noqa: E402

YEAST = "GCA_000146045.2"
PARADOXUS = "GCA_000000013.1"


class HelpersTest(unittest.TestCase):
    def test_translation_and_kmers(self):
        p = fx.PROTEINS["1at2759"]
        self.assertEqual(benchmarks.translate(fx.cds_of(p)), p + "*")
        self.assertEqual(benchmarks.translate("TAATAGTGA", 6), "QQ*")  # the ciliates' code
        self.assertEqual(benchmarks.translate("ATGNNNTGG"), "MXW")
        self.assertEqual(benchmarks.kmers("MKVLA*"), {"MKVLA"})
        self.assertEqual(benchmarks.kmers("MKVLAXPPQRS"), {"MKVLA", "PPQRS"})
        a, b = benchmarks.kmers(p), benchmarks.kmers(p[:60])
        self.assertEqual((benchmarks.shared(b, a), round(benchmarks.shared(a, b), 2)), (1.0, 0.48))  # 56 of 116

    def test_resolution_rows(self):
        rows = benchmarks.resolution_rows({"m": [(100.0, 5), (98.0, 5), (90.0, 2), (None, None)]})
        got = {r["shared_rank"]: r for r in rows}
        self.assertEqual(set(got), {"genus", "class", "all", "unmatched"})
        self.assertEqual((got["genus"]["copies"], got["genus"]["under_99"], got["genus"]["under_97"]), (2, 1, 0))
        self.assertEqual((got["all"]["copies"], got["all"]["under_99"], got["all"]["under_97"]), (4, 3, 2))
        self.assertEqual(got["unmatched"]["median"], "-")

    def test_best_ssu(self):
        rows = [dict(record="a", kind="rDNA", species_taxid="1", regions="18S:0-1800,ITS1:1800-2000", length="6000",
                     confidence="medium", source="genome"),
                dict(record="b", kind="18S", species_taxid="1", regions="-", length="1750", confidence="high",
                     source="EUKARYOME"),
                dict(record="c", kind="18S", species_taxid="2", regions="-", length="600", confidence="high",
                     source="genome"),  # too short
                dict(record="d", kind="rDNA", species_taxid="3", regions="18S:100-1900", length="9000",
                     confidence="high", source="genome")]
        self.assertEqual(benchmarks.best_ssu(rows), {"b": (1, 0, 1750), "d": (3, 100, 1900)})

    def test_queries_leave_out_own_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            odb = odb_mod.Odb(fx.odb_world(tmp), fx.LOCI)
            own = odb_mod.own_taxids(odb, ncbi)
            self.assertEqual(own[4932], [4932])
            lineage = list(ncbi.ancestors(4932))
            every = odb.queries(lineage, odb_mod.species_lineages(odb, ncbi))
            held = odb.queries(lineage, odb_mod.species_lineages(odb, ncbi), leave_out=own[4932])
            self.assertTrue(any("_4932_" in q for q in every))
            self.assertFalse(any("_4932_" in q for q in held))
            self.assertEqual(len(every) - len(held), len(fx.LOCI))  # one protein of each locus left out


class PackageTest(unittest.TestCase):
    def test_a_package_measured(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            root = fx.odb_world(tmp)
            odb = odb_mod.Odb(root, fx.LOCI)
            seq, gff, _spans = fx.genome()
            genomes_dir = os.path.join(tmp, "out", "genomes")
            os.makedirs(genomes_dir)
            for acc in (YEAST, PARADOXUS):  # S. paradoxus' genome = the yeast's: their copies identical
                with gzip.open(genomes.genome_file(genomes_dir, acc), "wt") as fh:
                    fh.write(f">c1\n{seq}\n")
            work = os.path.join(tmp, "work")
            vsearch = fx.vsearch(tmp)
            settings = dict(work=work, genomes=genomes_dir, miniprot=fx.miniprot(tmp, gff),
                            hmmsearch=fx.hmmsearch(tmp), threads=1, staged_only=True, odb=root, loci=fx.LOCI,
                            hmm_db=odb.hmm_db(os.path.join(tmp, "loci.hmm")),
                            odb_lineages=odb_mod.species_lineages(odb, ncbi), odb_own=odb_mod.own_taxids(odb, ncbi),
                            tiara=None, rrna=None)
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            jobs = {a: j for a, j in genomes.jobs(ncbi, gb).items() if a in (YEAST, PARADOXUS)}
            for acc in jobs:
                self.assertIsNone(genomes.run((jobs[acc], settings))[2])
            collected = package.collect(ncbi, jobs, gb, work, vsearch, 1)
            nearest = {(r["accession"], r["locus"]): (r["nearest_id"], r["nearest_depth"])
                       for r in collected["copy_rows"]}
            loci = [dict(locus=l, name="-", cutoff=odb.cutoff[l], hmm_length=odb.length[l]) for l in odb.loci]

            def measure(part):
                return benchmarks.run(part, ncbi, root, vsearch, 1, work, nearest, [("a step", 2.0)], log=lambda t: 0)
            out = package.write(os.path.join(tmp, "out", "Eukore_2026.11"), "2026.11", ncbi, collected, loci,
                                before_manifest=measure)
            self.assertIsNone(package.check(out, deep=True))
            with open(os.path.join(out, "manifest.json")) as fh:
                manifest = json.load(fh)
            self.assertIn("benchmarks/completeness.tsv", manifest["files"])
            self.assertEqual(manifest["benchmarks"]["completeness"]["median"], 0.667)
            folder = os.path.join(out, "benchmarks")

            def table(name):
                return fx.read_tsv(os.path.join(folder, name))
            # each genome: 1at2759 Duplicated, 2at2759 Single, 3at2759 Missing
            c = {r["accession"]: r for r in table("completeness.tsv")}
            self.assertEqual((c[YEAST]["single"], c[YEAST]["duplicated"], c[YEAST]["missing"], c[YEAST]["complete"],
                              c[YEAST]["supergroup"]), ("1", "1", "1", "0.667", "Fungi"))
            # the yeast is in BUSCO's dataset (its proteins the fixtures' genes), S. paradoxus is not
            o = {r["locus"]: r for r in table("orthodb.tsv")}
            self.assertEqual({r["accession"] for r in table("orthodb.tsv")}, {YEAST})
            self.assertEqual((o["1at2759"]["verdict"], o["1at2759"]["recall"], o["1at2759"]["guided"]),
                             ("concordant", "1.000", "t"))
            self.assertEqual((o["3at2759"]["status"], o["3at2759"]["verdict"]), ("Missing", "missed"))
            # its held-out calls, queried without the yeast's own proteins: the next species' (6253) find the genes
            self.assertTrue(os.path.isfile(os.path.join(genomes.done_dir(work, YEAST), genomes.HELD_OUT)))
            self.assertFalse(os.path.isfile(os.path.join(genomes.done_dir(work, PARADOXUS), genomes.HELD_OUT)))
            self.assertEqual((o["1at2759"]["held_out_verdict"], o["1at2759"]["held_out_query"].split("_")[1],
                              o["3at2759"]["held_out_verdict"]), ("concordant", "6253", "missed"))
            self.assertEqual(manifest["benchmarks"]["orthodb"]["held_out"]["verdicts"],
                             {"concordant": 2, "missed": 1})
            # identical copies in two species of one genus: nothing told apart at 99%
            r = {(x["marker"], x["shared_rank"]): x for x in table("resolution.tsv")}
            self.assertEqual((r[(benchmarks.ALL_BUSCOS, "genus")]["copies"], r[(benchmarks.ALL_BUSCOS, "all")]
                              ["under_99"], r[("1at2759", "genus")]["median"]), ("4", "0", "100.00"))
            self.assertEqual(len(table("runtime.tsv")), 2)
            self.assertEqual([s["step"] for s in table("stages.tsv")], ["a step", "benchmarks"])
            self.assertEqual({(x["supergroup"], x["reference"]) for x in table("labels.tsv")} >= {("all", "PR2")}, True)
            # measured again from the package's files alone (its taxonomy/, vsearch run here): the same figures
            again = os.path.join(tmp, "again")
            summary = benchmarks.run(out, ncbi_mod.Ncbi(os.path.join(out, "taxonomy")), root, vsearch, 1, work,
                                     out=again, log=lambda t: 0)
            self.assertEqual(summary["resolution"], manifest["benchmarks"]["resolution"])
            self.assertEqual(summary["orthodb"], manifest["benchmarks"]["orthodb"])
            self.assertFalse(os.path.isfile(os.path.join(again, "stages.tsv")))


if __name__ == "__main__":
    unittest.main()
