#!/usr/bin/env python3
"""Tests of the same-species merge (eukore/same_species.py, eukore_build.py's species_merge_step; since 2026-11-20):
sourmash's signature files read as sourmash writes them, containment and ANI as sourmash gives them; within a genus an
isolate-named species merged into the named species it is a strain of, a named pair and an isolate pair only
reviewed, a hybrid (by its genomes: two species wholly in it, it in neither; by its name) kept and reviewed with its
parents, an isolate whose genome is only in part a named species' reviewed; genomes without a sketch not compared;
the genera compared by a pool alike; the step off without sketches and with --no_species_merge. Then the merge in the
package: the merged genome under the named species from the copies' checks on (its copies kept under its accession,
not checked against their own species), a merge into a species whose genome is not done undone, the rRNA records
relabelled, the protal sketch's table taking the merged genome as the named species' further genome, genomes.tsv,
species.tsv, copy numbers, metadata and metadata/species_merges.tsv.

    python3 -m unittest discover -s tests -t tests
"""

import contextlib
import functools
import gzip
import io
import json
import os
import random
import sys
import tempfile
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, metadata, ncbi as ncbi_mod, odb as odb_mod, package, protal_sketch as ps  # noqa: E402
from eukore import rrna, same_species as ss, sources, store  # noqa: E402
import fixtures as fx  # noqa: E402

# the world's further species: isolates, a hybrid of two Saccharomyces, one named a hybrid, a congener of K. lactis
EXTRA_NODES = [
    (900001, 4930, "species", "Saccharomyces sp. M14"),
    (900002, 4930, "species", "Saccharomyces eubayanus"),
    (900003, 4930, "species", "Saccharomyces pastorianus"),
    (900005, 4930, "species", "Saccharomyces x bayanus"),
    (900008, 4930, "species", "Saccharomyces sp. P1"),
    (900009, 4930, "species", "Saccharomyces sp. N0"),
    (900010, 4930, "species", "Saccharomyces sp. Q"),
    (900011, 4930, "species", "Saccharomyces sp. tiny"),
    (900004, 4910, "species", "Kluyveromyces marxianus"),
    (900006, 5475, "species", "Candida sp. A"),
    (900007, 5475, "species", "Candida sp. B"),
]
SCALED = 10_000
MAX_HASH = round(2 ** 64 / SCALED)


def taxdump(path, extra=()):
    """The fixtures' taxonomy with EXTRA_NODES (and extra)."""
    nodes = fx.NODES
    fx.NODES = nodes + EXTRA_NODES + list(extra)
    try:
        return fx.taxdump(path)
    finally:
        fx.NODES = nodes


def signature(name, mins31, mins21=None, max_hash=MAX_HASH, gz=True, one_record=False):
    """A sourmash signature file's bytes as sourmash 4 writes `sketch dna -p k=21,k=31,scaled=10000,noabund`: a
    signature a ksize (one_record: sourmash 2's two sketches under one signature)."""
    def sketch(k, mins):
        return dict(num=0, ksize=k, seed=42, max_hash=max_hash, mins=sorted(set(mins)), md5sum="0" * 32,
                    molecule="DNA")
    sketches = [sketch(21, mins21 if mins21 is not None else [h // 3 for h in mins31]), sketch(31, mins31)]

    def record(parts):
        return {"class": "sourmash_signature", "email": "", "hash_function": "0.murmur64", "filename": f"{name}.fna.gz",
                "name": name, "license": "CC0", "signatures": parts, "version": 0.4}
    text = json.dumps([record(sketches)] if one_record else [record([s]) for s in sketches]).encode()
    return gzip.compress(text) if gz else text


class Hashes:
    """Random hashes under MAX_HASH, none drawn twice."""

    def __init__(self, seed=7):
        self.rng, self.seen = random.Random(seed), set()

    def __call__(self, n):
        out = []
        while len(out) < n:
            h = self.rng.randrange(1, MAX_HASH)
            if h not in self.seen:
                self.seen.add(h)
                out.append(h)
        return out


def accession(sp):
    return f"GCA_{sp:09d}.1"


def job(sp, why="ncbi", path=None):
    return dict(accession=accession(sp), species_taxid=sp, taxid=sp, why=[why], path=path, others=[])


def put_sketch(work, acc, data, pack=False):
    folder = genomes.done_dir(work, acc)
    store.put(folder, {metadata.SKETCH: data})
    if pack:
        store.pack(folder)


def sketches_world():
    """{species taxid: its genome's k-31 hashes} (None: no sketch) of three genera (module docstring's cases)."""
    h = Hashes()
    cer, eub, par = h(2000), h(2000), h(2000)
    lactis, cand_a = h(2000), h(2000)
    return {
        4932: cer,                                # S. cerevisiae
        900001: cer[:1600] + h(400),              # Saccharomyces sp. M14: a strain of S. cerevisiae
        900002: eub,                              # S. eubayanus
        900003: cer[:1600] + eub[:1600],          # S. pastorianus: S. cerevisiae x S. eubayanus
        27291: par,                               # S. paradoxus
        900005: par[:1800] + h(200),              # "Saccharomyces x bayanus": named a hybrid
        900008: cer[:1600] + h(4000),             # Saccharomyces sp. P1: S. cerevisiae's genome and much else
        900009: None,                             # Saccharomyces sp. N0: no sketch
        900011: h(50),                            # Saccharomyces sp. tiny: too few hashes
        28985: lactis,                            # K. lactis
        900004: lactis[:1200] + h(800),           # K. marxianus, 60% of K. lactis' hashes: a named pair
        5476: h(2000),                            # C. albicans
        900006: cand_a,                           # Candida sp. A
        900007: cand_a[:1700] + h(300),           # Candida sp. B, Candida sp. A's genome: an isolate pair
    }


class SketchTest(unittest.TestCase):
    def test_a_sourmash_signature_read(self):
        h = Hashes()
        mins = h(300)
        for data in (signature("x", mins), signature("x", mins, gz=False), signature("x", mins, one_record=True)):
            max_hash, got = ss.read_sketch(data)
            self.assertEqual((max_hash, got.tolist()), (MAX_HASH, sorted(mins)))
            self.assertEqual(ss.read_sketch(data, k=21)[1].tolist(), sorted({m // 3 for m in mins}))
        self.assertIsNone(ss.read_sketch(signature("x", mins), k=51))  # no sketch of that k
        for bad in (None, b"", b"\x1f\x8b not gzip", b"[not json", json.dumps({"signatures": 3}).encode()):
            self.assertIsNone(ss.read_sketch(bad))

    def test_containment_and_ani_as_sourmash_gives_them(self):
        h = Hashes()
        a = h(1000)
        b = a[:600] + h(1400)
        sa, sb = ss.read_sketch(signature("a", a)), ss.read_sketch(signature("b", b))
        self.assertEqual(ss.containments(sa, sb), (0.6, 0.3))
        self.assertAlmostEqual(ss.ani_of(0.6), 0.6 ** (1 / 31))
        self.assertEqual(ss.ani_of(0.0), 0.0)
        # another scale: both cut to the coarser one (the smaller max_hash) first
        coarse = MAX_HASH // 2
        sc = ss.read_sketch(signature("c", [x for x in a if x <= coarse], max_hash=coarse))
        below = [x for x in a if x <= coarse]
        self.assertEqual(ss.containments(sa, sc), (1.0, 1.0))
        self.assertEqual(ss.containments(sb, sc)[0], sum(x <= coarse for x in a[:600]) / sum(x <= coarse for x in b))
        self.assertEqual(ss.containments(sc, sb)[0], sum(x <= coarse for x in a[:600]) / len(below))

    def test_names(self):
        for name in ("Saccharomyces sp. M14", "Saccharomyces sp. KACC 48362", "Blastocystis sp. subtype 4",
                     "Blastocystis sp. NMH", "Fusarium sp.", "Candidatus Nitrosopumilus sp. X"):
            self.assertEqual(ss.kind_of(name), "isolate", name)
        for name in ("Saccharomyces x bayanus", "Fragaria x ananassa", "Mentha × piperita",
                     "× Triticosecale rimpaui"):
            self.assertEqual(ss.kind_of(name), "hybrid", name)
        for name in ("Saccharomyces cerevisiae", "Saccharomyces pastorianus", "Spathaspora xylofermentans",
                     "Pichia spartinae", "Blastocystis hominis", "Aspergillus cf. niger"):
            self.assertEqual(ss.kind_of(name), "named", name)


class RulesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        d = self.tmp.name
        self.ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(d, "taxdump.tar.gz")))
        self.work = os.path.join(d, "work")
        self.world = sketches_world()
        self.jobs = {accession(sp): job(sp) for sp in self.world}
        # a transcriptome and EukProt's records of the genus: not compared
        self.jobs["GAAA01000000"] = dict(job(900010, "tsa"), accession="GAAA01000000")
        self.jobs["EUK_1"] = dict(job(900010, "no_assembly", path="x.fna.gz"), accession="EUK_1")
        for i, (sp, mins) in enumerate(sorted(self.world.items())):
            if mins is not None:  # one genome in two packed, as WORK holds them since 2026-11-10
                put_sketch(self.work, accession(sp), signature(accession(sp), mins), pack=i % 2 == 0)

    def tearDown(self):
        store.close()
        self.tmp.cleanup()

    def found(self, threads=1):
        said = []
        merges = ss.find(self.jobs, self.ncbi, self.work, threads, log=said.append)
        return merges, {(r["from_taxid"], r["rule"]): r for r in merges.rows}, said

    def test_the_rules(self):
        merges, rows, said = self.found()
        # the isolate of S. cerevisiae's genome merged into it; nothing else merged
        self.assertEqual(merges.species, {900001: 4932})
        self.assertTrue(merges)
        self.assertEqual(merges.accessions(), {accession(900001): 4932})
        m = rows[(900001, ss.MERGED)]
        self.assertEqual((m["from_name"], m["from_accession"], m["to_taxid"], m["to_name"], m["to_accession"]),
                         ("Saccharomyces sp. M14", accession(900001), 4932, "Saccharomyces cerevisiae",
                          accession(4932)))
        self.assertEqual((m["ani"], m["containment"], m["share"], m["note"]),
                         (f"{0.8 ** (1 / 31):.4f}", "0.8000", "1.000", "-"))
        # the hybrid by its genomes: both parents', the named first, the nearer (here alike: the lower taxid) as `to`;
        # never a merge target (M14 is wholly in it too)
        p = rows[(900003, ss.HYBRID)]
        self.assertEqual((p["to_taxid"], p["other_taxid"], p["other_name"], p["containment"], p["share"]),
                         (4932, 900002, "Saccharomyces eubayanus", "0.5000", "0.625"))
        self.assertIn("Saccharomyces eubayanus: ANI", p["note"])
        # the species named a hybrid: kept, reviewed with its nearest
        x = rows[(900005, ss.HYBRID)]
        self.assertEqual((x["to_taxid"], x["other_taxid"], x["note"]), (27291, "-", "named a hybrid"))
        # an isolate whose genome is S. cerevisiae's and much else: not merged, reviewed
        q = rows[(900008, ss.PARTIAL)]
        self.assertEqual((q["to_taxid"], q["share"]), (4932, f"{1600 / 5600 / 0.8:.3f}"))
        # two named species of one genome: reviewed, not merged
        n = rows[(28985, ss.NAMED_PAIR)]
        self.assertEqual((n["to_taxid"], n["containment"], n["share"]), (900004, "0.6000", "1.000"))
        # two isolates of one genome, no named species to merge into: reviewed
        i = rows[(900006, ss.ISOLATE_PAIR)]
        self.assertEqual((i["to_taxid"], i["containment"], i["note"]), (900007, "0.8500", "-"))
        # the merged isolate and P1, which holds its genome: reviewed, P1 not merged through it
        i = rows[(900001, ss.ISOLATE_PAIR)]
        self.assertEqual((i["to_taxid"], i["note"]),
                         (900008, "Saccharomyces sp. M14 merged into Saccharomyces cerevisiae (4932)"))
        # no other row: S. paradoxus, S. eubayanus, C. albicans are no one's; N0 (no sketch) and tiny not compared
        self.assertEqual(sorted(rows), sorted([(900001, ss.MERGED), (900003, ss.HYBRID), (900005, ss.HYBRID),
                                               (900008, ss.PARTIAL), (28985, ss.NAMED_PAIR),
                                               (900006, ss.ISOLATE_PAIR), (900001, ss.ISOLATE_PAIR)]))
        self.assertEqual(merges.about, dict(state="made", table=ss.TABLE, k=31, ani_min=0.95, share_min=0.7,
                                            hashes_min=100, genera=3, genomes=12, without_sketch=2,
                                            pairs=21 + 1 + 3,
                                            rules={ss.MERGED: 1, ss.NAMED_PAIR: 1, ss.HYBRID: 2,
                                                   ss.ISOLATE_PAIR: 2, ss.PARTIAL: 1}))
        self.assertEqual(len(said), 1)
        self.assertIn("1 isolate-named species merged into a named species (e.g. Saccharomyces sp. M14 -> "
                      "Saccharomyces cerevisiae), to review: 1 named pair, 2 hybrid, 2 isolate pair, 1 partial",
                      said[0])
        self.assertIn("2 without a sketch, not compared", said[0])
        # the genera compared by a pool: the same rows
        self.assertEqual(self.found(threads=3)[0].rows, merges.rows)

    def test_merges_never_chain_nor_cross_a_genus(self):
        # Kluyveromyces sp. Z's genome is Saccharomyces cerevisiae's: another genus, never compared with it.
        # Saccharomyces sp. M14b: M14's own hashes and 80 of those it shares with S. cerevisiae, wholly in M14 and not
        # near S. cerevisiae (0.17 of it there): never merged through M14, reviewed beside it
        self.ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(self.tmp.name, "taxdump2.tar.gz"), [
            (900012, 4910, "species", "Kluyveromyces sp. Z"), (900013, 4930, "species", "Saccharomyces sp. M14b")]))
        self.jobs[accession(900012)] = job(900012)
        self.jobs[accession(900013)] = job(900013)
        m14 = self.world[900001]
        put_sketch(self.work, accession(900012), signature("z", self.world[4932]))
        put_sketch(self.work, accession(900013), signature("b", m14[1600:] + m14[:80]))
        merges, rows, _said = self.found()
        self.assertEqual(merges.species, {900001: 4932})
        self.assertFalse([k for k in rows if k[0] == 900012])  # Z: no row (its genus' species are not its genome)
        b = rows[(900013, ss.ISOLATE_PAIR)]
        self.assertEqual((b["to_taxid"], b["containment"], b["share"], b["note"]),
                         (900001, "1.0000", "1.000",
                          "Saccharomyces sp. M14 merged into Saccharomyces cerevisiae (4932)"))
        self.assertNotIn((900013, ss.MERGED), rows)
        self.assertNotIn((900013, ss.PARTIAL), rows)

    def test_merges_relabel_undo(self):
        merges, _rows, _said = self.found()
        records = [("a", "ACGT", dict(species_taxid=900001)), ("b", "ACGT", dict(species_taxid="900001")),
                   ("c", "ACGT", dict(species_taxid=27291)), ("d", "ACGT", dict(species_taxid="-"))]
        self.assertIs(merges.relabel(records), records)
        self.assertEqual([r["species_taxid"] for _n, _s, r in records], [4932, "4932", 27291, "-"])
        jobs = merges.jobs({"x": dict(species_taxid=900001, a=1), "y": dict(species_taxid=27291)})
        self.assertEqual((jobs["x"], jobs["y"]), (dict(species_taxid=4932, a=1), dict(species_taxid=27291)))
        self.assertEqual((merges.of(900001), merges.of("27291")), (4932, 27291))
        merges.undo(900001, "the named species' genome is not done")
        self.assertEqual((merges.species, bool(merges)), ({}, False))
        row = next(r for r in merges.rows if r["from_taxid"] == 900001)
        self.assertEqual((row["rule"], row["note"]), (ss.UNDONE, "the named species' genome is not done"))
        self.assertEqual(merges.about["rules"].get(ss.UNDONE), 1)
        self.assertNotIn(ss.MERGED, merges.about["rules"])
        records = [("a", "ACGT", dict(species_taxid=900001))]
        self.assertEqual(merges.relabel(records)[0][2]["species_taxid"], 900001)  # undone: as it was


class StepTest(unittest.TestCase):
    def test_off_without_sketches_or_by_its_flag(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            work = os.path.join(tmp, "work")
            world = sketches_world()
            jobs = {accession(sp): job(sp) for sp in world}
            for sp, mins in world.items():
                if mins is not None:
                    put_sketch(work, accession(sp), signature(accession(sp), mins))
            sourmash = fx.tool(tmp, "sourmash", "#!/bin/sh\nexit 0\n")

            def step(**given):
                opts = types.SimpleNamespace(**dict(dict(no_species_merge=False, no_sketches=False, sourmash=sourmash,
                                                         work=work, threads=1, report_minutes=10), **given))
                stages, printed = [], io.StringIO()
                with contextlib.redirect_stdout(printed):
                    merges = eukore_build.species_merge_step(opts, jobs, ncbi, stages)
                return merges, stages, printed.getvalue()
            merges, stages, said = step(no_species_merge=True)
            self.assertEqual((merges.species, merges.rows, merges.about, stages),
                             ({}, [], {"state": "off (--no_species_merge)"}, []))
            self.assertIn("same-species merge: none (--no_species_merge)", said)
            merges, stages, said = step(no_sketches=True)  # the sketches there are not taken
            self.assertEqual((merges.species, merges.rows, merges.about["state"], stages),
                             ({}, [], "no sketches (--no_sketches, or no sourmash)", []))
            self.assertIn("no sourmash sketches of the genomes to compare them by", said)
            merges, stages, said = step(sourmash=os.path.join(tmp, "no_such_sourmash"))
            self.assertEqual((merges.species, merges.about["state"]),
                             ({}, "no sketches (--no_sketches, or no sourmash)"))
            merges, stages, said = step()
            self.assertEqual((merges.species, merges.about["state"], [s for s, _t in stages]),
                             ({900001: 4932}, "made", ["same-species merge"]))
            self.assertIn("package: same-species merge ..", said)
            store.close()

    def test_the_flag_and_the_banner(self):
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()) as printed:
            # a run that stops at its sources (none in WORK) after the banner
            with self.assertRaises(SystemExit):
                eukore_build.main(["--build_dir", tmp, "--out", tmp, "--stage", "package", "--no_species_merge"])
        self.assertRegex(printed.getvalue(), r"off\s+.*same-species merge")


class PackageTest(unittest.TestCase):
    def test_the_merge_in_the_package(self):
        ISO, Q, PARADOXUS, YEAST = "GCA_900000001.1", "GCA_900000010.1", "GCA_000000013.1", "GCA_000146045.2"
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            root = fx.odb_world(tmp)
            odb = odb_mod.Odb(root, fx.LOCI)
            seq, gff, _spans = fx.genome()
            genomes_dir = os.path.join(tmp, "out", "genomes")
            os.makedirs(genomes_dir)
            for acc in (YEAST, ISO):  # the isolate's genome is S. cerevisiae's
                with gzip.open(genomes.genome_file(genomes_dir, acc), "wt") as fh:
                    fh.write(f">c1\n{seq}\n")
            work = os.path.join(tmp, "work")
            settings = dict(work=work, genomes=genomes_dir, miniprot=fx.miniprot(tmp, gff),
                            hmmsearch=fx.hmmsearch(tmp), threads=1, staged_only=True, odb=root, loci=fx.LOCI,
                            hmm_db=odb.hmm_db(os.path.join(tmp, "loci.hmm")),
                            odb_lineages=odb_mod.species_lineages(odb, ncbi), tiara=None, rrna=None,
                            vsearch=fx.vsearch(tmp))
            rows = fx.GENBANK + [fx.assembly(ISO, 900001, 900001, "Saccharomyces sp. M14", "Contig", 10),
                                 fx.assembly(Q, 900010, 900010, "Saccharomyces sp. Q", "Contig", 10)]
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), rows)])
            jobs = genomes.jobs(ncbi, gb)
            for acc in (YEAST, ISO):
                self.assertIsNone(genomes.run((jobs[acc], settings))[2])
            # the sketches: the isolate S. cerevisiae's strain; Q S. paradoxus', whose genome is not done
            h = Hashes()
            cer, par = h(2000), h(2000)
            for acc, mins in ((YEAST, cer), (ISO, cer[:1600] + h(400)), (PARADOXUS, par), (Q, par[:1800] + h(200))):
                put_sketch(work, acc, signature(acc, mins))
            merges = ss.find(jobs, ncbi, work, log=lambda *_: None)
            self.assertEqual(merges.species, {900001: 4932, 900010: 27291})
            # without the merge: the isolate another species of the genus, S. cerevisiae's copies identical to its
            plain = package.collect(ncbi, jobs, gb, work, fx.vsearch(tmp), 1)
            first = {(r["accession"], r["locus"]): r for r in plain["copy_rows"]}
            self.assertEqual((first[(YEAST, "1at2759")]["same_as"], plain["genomes"][ISO]["species_taxid"],
                              plain["genomes"][ISO]["representative"]), ("4930", 900001, "t"))
            collected = package.collect(ncbi, jobs, gb, work, fx.vsearch(tmp), 1, merges=merges)
            # the merge into S. paradoxus undone: its genome is not done
            self.assertEqual(merges.species, {900001: 4932})
            undone = next(r for r in merges.rows if r["from_taxid"] == 900010)
            self.assertEqual((undone["rule"], undone["note"]), (ss.UNDONE, "the named species' genome is not done"))
            # the isolate's genome S. cerevisiae's further genome, its NCBI taxid and name kept; its copies kept under
            # its accession and checked as S. cerevisiae's own: no other species to be the same as
            g = collected["genomes"][ISO]
            self.assertEqual((g["species_taxid"], g["representative"], g["source"], g["taxid"], g["organism"]),
                             (4932, "f", ss.SOURCE, 900001, "Saccharomyces sp. M14"))
            self.assertEqual(ss.SOURCE, package.MERGED)
            self.assertNotEqual(g.get("tier", "-"), "-")  # measured as a species' genome is
            copies = {(r["accession"], r["locus"]): r for r in collected["copy_rows"]}
            self.assertEqual(copies[(YEAST, "1at2759")]["same_as"], "-")
            self.assertIn((ISO, "1at2759"), copies)
            self.assertEqual(collected["genomes"][YEAST]["representative"], "t")
            # the protal sketch's table: the merged genome S. cerevisiae's further genome under its key
            table, missing = ps.table_rows(collected["genomes"], functools.partial(genomes.genome_file, genomes_dir),
                                           collected["dropped"])
            self.assertEqual([(k, a) for k, a, _p in table], [(YEAST, YEAST), (YEAST, ISO)])
            # the isolate's rRNA records S. cerevisiae's (the build relabels them before their checks)
            records = merges.relabel(rrna.genome_records(ISO, [dict(kind="18S", sequence="ACGT" * 450, regions=[],
                                                                    lsu_complete=False, contig="c1", start=1,
                                                                    end=1800, strand="+")], 900001, 900001))
            self.assertEqual([(r["species_taxid"], r["taxid"]) for _n, _s, r in records], [(4932, 900001)])
            # the package: genomes.tsv, species.tsv, copies, copy numbers, metadata, the merges' table
            loci = [dict(locus=l, name=odb.names.get(l, "-"), cutoff=odb.cutoff[l], hmm_length=odb.length[l])
                    for l in odb.loci]
            meta = functools.partial(metadata.write, work=work, tables=[], jobs=jobs)
            out = package.write(os.path.join(tmp, "out", "Eukore_2026.11"), "2026.11", ncbi, collected, loci,
                                records, genomes_dir, {"species_merge": merges.about}, None, [], metadata=meta,
                                tables=[(ss.TABLE, ss.COLUMNS, merges.rows)])
            self.assertIsNone(package.check(out, deep=True))
            gt = {r["accession"]: r for r in package.read_table(out, "genomes.tsv")}
            self.assertEqual((gt[ISO]["species_taxid"], gt[ISO]["representative"], gt[ISO]["source"],
                              gt[ISO]["taxid"]), ("4932", "f", "merged", "900001"))
            species = {r["taxid"]: r for r in package.read_table(out, "species.tsv")}
            self.assertNotIn("900001", species)
            self.assertEqual((species["4932"]["representative"], species["4932"]["rrna_records"]), (YEAST, "1"))
            self.assertIn(ISO, {r["accession"] for r in package.read_table(out, "busco/copies.tsv")})
            self.assertIn(ISO, package.read_copies(out, "1at2759"))
            numbers = package.read_gz(os.path.join(out, package.COPY_NUMBERS))
            self.assertEqual({r["species_taxid"] for r in numbers if r["accession"] == ISO}, {"4932"})
            assembly = {r["accession"]: r for r in package.read_gz(os.path.join(out, "metadata", "assembly.tsv.gz"))}
            self.assertEqual(sorted(assembly), [YEAST, ISO])
            table = package.read_table(out, ss.TABLE)
            self.assertEqual([(r["from_taxid"], r["to_taxid"], r["rule"]) for r in table],
                             [("900001", "4932", ss.MERGED), ("900010", "27291", ss.UNDONE)])
            manifest = package.read(out)
            self.assertEqual(manifest["inputs"]["species_merge"]["rules"], {ss.MERGED: 1, ss.UNDONE: 1})
            self.assertIn(ss.TABLE, manifest["files"])
            store.close()


if __name__ == "__main__":
    unittest.main()
