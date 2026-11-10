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
        # a EUKARYOME read labelled by its genus alone is no species' (its genus' taxid its species_taxid)
        genus = dict(rows[1], record="g", species_taxid="4930", label_from="eukaryome_genus")
        self.assertNotIn("g", benchmarks.best_ssu(rows + [genus]))

    def test_reads_simulated_and_their_hits(self):
        # genome reads from anywhere, each within its contig (bisect over the contigs' ends), in random order: the
        # first `depth` a sample of the whole genome, not of its first contigs
        import random  # noqa: PLC0415
        from eukore import reads  # noqa: PLC0415
        rng = random.Random(3)
        with tempfile.TemporaryDirectory() as tmp:
            genome = os.path.join(tmp, "g.fna")
            contigs = {f"c{i}": "".join(rng.choice("ACGT") for _ in range(2000)) for i in range(5)}
            with open(genome, "w") as fh:
                fh.writelines(f">{n}\n{s[:1000]}\n{s[1000:]}\n" for n, s in contigs.items())
            from unittest import mock  # noqa: PLC0415
            with mock.patch.object(reads, "ERROR", 0):  # exact reads: each found in its contig
                _m, anywhere = reads.simulate(genome, [], random.Random(1), 400)
            self.assertGreater(len(anywhere), 350)
            of = [[n for n, s in contigs.items() if r in s or reads.revcomp(r) in s] for r in anywhere]
            self.assertTrue(all(len(n) == 1 for n in of))  # each read within one contig
            self.assertGreaterEqual(len({n[0] for n in of[:40]}), 4)  # the first reads from across the genome
            # a read's hits together (minimap2's --split-prefix merges an index's parts): one group a read
            paf = os.path.join(tmp, "h.paf")
            with open(paf, "w") as fh:
                for read, target in (("r1", "t1"), ("r1", "t2"), ("r2", "t1")):
                    fh.write(f"{read}\t150\t0\t150\t+\t{target}\t1000\t0\t150\t148\t150\t60\n")
            self.assertEqual([(r, len(h)) for r, h in reads.paf_groups(paf)], [("r1", 2), ("r2", 1)])

    def test_place_of(self):
        call = dict(contig="c1", start="101", end="200")
        self.assertEqual(benchmarks.place_of(call, dict(contig="c1", start="101", end="200")), ("same", "1.000"))
        self.assertEqual(benchmarks.place_of(call, dict(contig="c1", start="151", end="250")), ("same", "0.333"))
        self.assertEqual(benchmarks.place_of(call, dict(contig="c2", start="101", end="200")), ("elsewhere", "0.000"))
        self.assertEqual(benchmarks.place_of(call, dict(contig="c1", start="201", end="300")), ("elsewhere", "0.000"))
        self.assertEqual(benchmarks.place_of(call, dict(contig="c1", start="-", end="-")), ("-", "-"))  # older file

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
            # on the package call's gene, with its span
            self.assertEqual((o["1at2759"]["held_out_place"], o["1at2759"]["held_out_span"],
                              o["1at2759"]["held_out_rescued"], o["3at2759"]["held_out_place"]),
                             ("same", "1.000", "f", "-"))
            self.assertEqual(manifest["benchmarks"]["orthodb"]["held_out"]["place"],
                             {"both_called": 2, "same_gene": 1.0, "median_span": 1.0, "span_at_least_90": 1.0,
                              "rescued": 0})
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


class CladeHoldOutTest(unittest.TestCase):
    def test_the_clade_left_out_rank_by_rank(self):
        # since 2026-11-09: a dataset species' genome called again without its genus', family's, order's species
        ranks = {4932: dict(genus=4930, family=4893, order=4892, **{"class": 4891}, phylum=4890),
                 559292: dict(genus=4930, family=4893, order=4892, **{"class": 4891}, phylum=4890),
                 27291: dict(genus=4930, family=4893, order=4892, **{"class": 4891}, phylum=4890),
                 28985: dict(genus=4910, family=4893, order=4892, **{"class": 4891}, phylum=4890),
                 5476: dict(genus=5475, family=1535326, order=4892, **{"class": 4891}, phylum=4890),
                 6253: dict(genus=6251, phylum=6231)}
        job = dict(accession="Y", species_taxid=4932, ranks=ranks[4932])
        settings = dict(odb_own={4932: [559292]}, odb_ranks={t: r for t, r in ranks.items() if t != 4932},
                        held_out_ranks=list(genomes.HELD_OUT_RANKS))
        left = {r: genomes.clade_left_out(job, settings, r) for r in genomes.HELD_OUT_RANKS}
        self.assertEqual(left, {"genus": {559292, 27291}, "family": {559292, 27291, 28985},
                                "order": {559292, 27291, 28985, 5476}})
        self.assertEqual([genomes.nearest_kept(job, settings, s) for s in ({559292}, *left.values())], [5, 4, 3, 0])
        # an order that leaves out no more than the family: not run again; a lineage without an order: none
        settings["odb_ranks"].pop(5476)
        got, why = genomes.clade_held_out("-", job, settings, None, "-", "-", "order")
        self.assertEqual((got, why), (dict(left_out=3, nearest=0, same_as="family"), None))
        self.assertEqual(genomes.clade_held_out("-", dict(job, ranks={}), settings, None, "-", "-", "order")[0],
                         dict(left_out=0, no_taxon=True))
        self.assertEqual(genomes.held_out_file("genus"), "held_out_genus.tsv")
        self.assertEqual(genomes.held_out_file(), genomes.HELD_OUT)

    def test_the_levels_measured(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = genomes.done_dir(tmp, "Y")
            os.makedirs(folder)
            with open(os.path.join(folder, "done.json"), "w") as fh:
                json.dump({"held_out": {"left_out": [559292], "nearest": 5, "ranks": {
                    "genus": {"left_out": 2, "nearest": 4}, "family": {"left_out": 2, "nearest": 4, "same_as": "genus"},
                    "order": {"left_out": 0, "no_taxon": True}}}}, fh)
            for name in ("held_out.tsv", "held_out_genus.tsv"):
                with open(os.path.join(folder, name), "w") as fh:
                    fh.write("\t".join(genomes.HELD_OUT_COLUMNS) + "\n1at2759\tSingle\tq\tc1\t1\t90\tf\tATG\n")
            levels = benchmarks.clade_levels(tmp, "Y")
            self.assertEqual(sorted(levels), ["family", "genus", "species"])  # the family's: the genus' calls
            self.assertEqual(levels["family"][0]["1at2759"]["status"], "Single")
        tallies = {("Y", 4932, "species"): dict(loci=10, missed=1, concordant=8, similar=1, same_gene=9, span_90=8),
                   ("Y", 4932, "genus"): dict(loci=10, missed=4, concordant=3, similar=2, discordant=1, same_gene=5,
                                              span_90=3)}
        import collections  # noqa: PLC0415
        tallies = {k: collections.Counter(v) for k, v in tallies.items()}
        rows = benchmarks.holdout_rows(tallies, {("Y", "species"): dict(left_out=1, nearest=5),
                                                 ("Y", "genus"): dict(left_out=2, nearest=4)}, {})
        self.assertEqual([(r["level"], r["nearest"], r["called"]) for r in rows],
                         [("species", "genus", 9), ("genus", "family", 6)])
        summary = benchmarks.holdout_summary(rows)
        self.assertEqual((summary["genus"]["called"], summary["genus"]["concordant_of_called"]), (0.6, 0.5))


CENSOR_MINIMAP2 = '''#!/usr/bin/env python3
"""minimap2's stand-in for the read benchmark's censored reads: read 0's primary alignment and 50 secondary ones are
all to S. cerevisiae's records (y0 .. y50); read 1's too, but for one to S. paradoxus' (p); read 2 has one alignment."""
for read, targets in (("Y|marker|0", [f"y{i}" for i in range(51)]),
                      ("Y|marker|1", [f"y{i}" for i in range(50)] + ["p"]), ("Y|marker|2", ["y0"])):
    for i, target in enumerate(targets):
        print(read, 150, 0, 150, "+", target, 1000, 0, 150, 148, 150, 60, "tp:A:" + "SP"[i == 0], sep="\\t")
'''


class CensoredReadsTest(unittest.TestCase):
    """reads.py's `censored` (since 2026-11-10): a read held_out or novel does not assign because the records left out
    took minimap2's whole list of secondary alignments; counted apart, the classes as they were."""

    def test_a_full_list_of_records_left_out(self):
        from eukore import reads  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            taxa = dict({f"y{i}": 4932 for i in range(60)}, p=27291, h=9606)
            yeast, full = {f"y{i}" for i in range(51)}, reads.MAX_SECONDARY
            order = ("order", 4892)  # the yeasts': S. paradoxus' too
            self.assertTrue(reads.censored(yeast, full, taxa, ncbi, leave_out=4932))
            self.assertTrue(reads.censored(yeast, full, taxa, ncbi, clade=order))
            self.assertFalse(reads.censored(yeast, full, taxa, ncbi))  # nothing left out (`in`)
            self.assertFalse(reads.censored(yeast, full - 1, taxa, ncbi, leave_out=4932))  # the list not full
            self.assertFalse(reads.censored(yeast | {"p"}, full, taxa, ncbi, leave_out=4932))  # a record kept
            self.assertTrue(reads.censored(yeast | {"p"}, full, taxa, ncbi, clade=order))
            self.assertFalse(reads.censored(yeast | {"h"}, full, taxa, ncbi, clade=order))
            self.assertFalse(reads.censored(yeast | {"unknown"}, full, taxa, ncbi, leave_out=4932))
            paf = os.path.join(tmp, "h.paf")
            with open(paf, "w") as fh:
                for i in range(12):  # a read's alignments: the primary and 11 secondary, 2 of them over half the read
                    fh.write(f"r1\t150\t0\t{150 if i < 10 else 75}\t+\ty{i}\t1000\t0\t150\t148\t150\t60\t"
                             f"NM:i:2\ttp:A:{'S' if i else 'P'}\n")
            self.assertEqual([(r, len(h), len(t), n) for r, h, t, n in reads.paf_reads(paf)], [("r1", 10, 12, 11)])
            self.assertEqual([(r, len(h)) for r, h in reads.paf_groups(paf)], [("r1", 10)])  # as before

    def test_counted_apart_in_the_table(self):
        from unittest import mock  # noqa: PLC0415
        from eukore import reads  # noqa: PLC0415
        self.assertEqual((reads.COLUMNS[9], reads.COLUMNS[10], reads.COLUMNS[-1], len(reads.COLUMNS)),
                         ("assigned", "species_correct", "censored", 10 + 2 * len(reads.RANKS) + 1))
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            taxa = dict({f"y{i}": 4932 for i in range(60)}, p=27291)
            ref = os.path.join(tmp, "ref.fna")
            with open(ref, "w") as fh:
                fh.write(">y0\nACGT\n")
            tables = {"genomes.tsv": [dict(accession="Y", representative="t", species_taxid="4932")],
                      "species.tsv": [dict(taxid="4932", group="Fungi")]}
            with mock.patch.object(reads.package, "read_table", lambda folder, name: tables[name]), \
                    mock.patch.object(reads, "references", return_value={"Eukore": (ref, taxa)}), \
                    mock.patch.object(reads, "marker_spans", return_value={}), \
                    mock.patch.object(reads, "simulate", return_value=(["ACGT" * 40] * 3, [])):
                rows = reads.run(tmp, ncbi, tmp, ["Y"], fx.tool(tmp, "minimap2", CENSOR_MINIMAP2), 1,
                                 os.path.join(tmp, "work"), log=lambda *a: None)
            got = {r["mode"]: (r["n"], r["assigned"], r["censored"], r["species_correct"], r["species_wrong"])
                   for r in rows if r["reads"] == "marker"}
            # in: all three assigned (read 1 to the genus, its congener's record as good as its own). held out: read 1
            # its congener's (wrong at the species), read 0 censored, read 2 not assigned and not censored (its list
            # was not full). novel: reads 0 and 1 censored
            self.assertEqual(got, {"in": (3, 3, 0, 2, 0), "held_out": (3, 1, 1, 0, 1), "novel": (3, 0, 2, 0, 0)})
            self.assertTrue(all(r["censored"] == 0 for r in rows if r["reads"] == "genome"))


if __name__ == "__main__":
    unittest.main()
