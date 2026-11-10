#!/usr/bin/env python3
"""Tests of BGZF (util.bgzf_open: every file Eukore compresses, since 2026-11-09) and of a package's release
(release.py: the website's data and the downloads).

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import io
import json
import os
import random
import struct
import sys
import tarfile
import tempfile
import unittest
import zipfile
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import package, release, util  # noqa: E402


def blocks(path):
    """A BGZF file's blocks -> [(compressed size, uncompressed size)], each checked as htslib reads it (its BC
    subfield's size, its CRC32), raising AssertionError where one is not BGZF."""
    out = []
    with open(path, "rb") as fh:
        data = fh.read()
    at = 0
    while at < len(data):
        id1, id2, cm, flg, _mtime, _xfl, _os, xlen = struct.unpack("<4BIBBH", data[at:at + 12])
        assert (id1, id2, cm, flg & 4) == (31, 139, 8, 4), f"no gzip header with FEXTRA at {at}"
        si1, si2, slen, bsize = struct.unpack("<2BHH", data[at + 12:at + 18])
        assert (si1, si2, slen, xlen) == (66, 67, 2, 6), f"no BC subfield at {at}"
        body = data[at + 12 + xlen:at + bsize + 1 - 8]
        crc, isize = struct.unpack("<II", data[at + bsize + 1 - 8:at + bsize + 1])
        plain = zlib.decompress(body, -15)
        assert len(plain) == isize and zlib.crc32(plain) == crc, f"block at {at}: size or CRC wrong"
        assert bsize + 1 <= 65536 and isize <= util.BGZF_BLOCK
        out.append((bsize + 1, isize))
        at += bsize + 1
    assert out and out[-1] == (28, 0), "no BGZF EOF block"
    return out


def dna(n, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def content(path):
    with open(path, "rb") as fh:
        return fh.read()


class Bgzf(unittest.TestCase):
    def tearDown(self):
        util.compression_threads(1)

    def test_a_text_file_read_back_by_gzip(self):
        text = "".join(f">seq{i}\n{dna(500, i)}\n" for i in range(400))  # ~200 kB: four blocks
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x.fna.gz")
            with util.bgzf_open(path) as fh:
                fh.write(text)
            got = blocks(path)
            self.assertEqual(len(got), 5)
            self.assertTrue(all(isize == util.BGZF_BLOCK for _c, isize in got[:-2]))
            self.assertTrue(util.is_bgzf(path))
            with gzip.open(path, "rt") as fh:
                self.assertEqual(fh.read(), text)
            with util.open_text(path) as fh:
                self.assertEqual(fh.readline(), ">seq0\n")

    def test_threads_write_the_same_file(self):
        data = dna(1_000_000, 7).encode()
        with tempfile.TemporaryDirectory() as tmp:
            one, many = os.path.join(tmp, "one.gz"), os.path.join(tmp, "many.gz")
            with util.bgzf_open(one, "wb") as fh:
                for i in range(0, len(data), 7777):  # writes of any size
                    fh.write(data[i:i + 7777])
            util.compression_threads(4)
            with util.bgzf_open(many, "wb") as fh:
                fh.write(data)
            self.assertEqual(content(one), content(many))
            self.assertEqual(gzip.decompress(content(many)), data)

    def test_an_empty_file_and_the_levels(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = os.path.join(tmp, "empty.gz")
            with util.bgzf_open(empty, "wb"):
                pass
            self.assertEqual(content(empty), util.BGZF_EOF)
            self.assertEqual(gzip.decompress(util.BGZF_EOF), b"")
            text = ("ACGTTGCA" * 20 + "\n") * 5000
            sizes = {}
            for level in (util.WORK_LEVEL, util.BGZF_LEVEL):
                path = os.path.join(tmp, f"{level}.gz")
                with util.bgzf_open(path, level=level) as fh:
                    fh.write(text)
                sizes[level] = os.path.getsize(path)
                blocks(path)
            self.assertLessEqual(sizes[util.BGZF_LEVEL], sizes[util.WORK_LEVEL])
            with self.assertRaises(ValueError):
                util.bgzf_open(os.path.join(tmp, "x.gz"), "rt")

    def test_a_copy_of_plain_gzip_and_bgzf(self):
        text = "a\tb\n" * 10000
        with tempfile.TemporaryDirectory() as tmp:
            plain, gz = os.path.join(tmp, "t.tsv"), os.path.join(tmp, "t.tsv.gz")
            with open(plain, "w", newline="\n") as fh:
                fh.write(text)
            with gzip.open(gz, "wt") as fh:
                fh.write(text)
            self.assertFalse(util.is_bgzf(gz))
            for src in (plain, gz):
                dest = os.path.join(tmp, "out.gz")
                util.bgzf_copy(src, dest)
                blocks(dest)
                self.assertEqual(gzip.decompress(content(dest)).decode(), text)
                self.assertFalse(os.path.exists(util.part_of(dest)))

    def test_the_package_tables_are_bgzf(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.tsv.gz")
            package.write_gz(path, ("a", "b"), [dict(a=1, b=None), dict(a="x", b="y")])
            blocks(path)
            self.assertEqual(package.read_gz(path), [dict(a="1", b="-"), dict(a="x", b="y")])


SPECIES = ("taxid", "name", "representative", "source", "loci", "rrna_records", "organelle_genes", "supergroup", "group",
           "ncbi_lineage", "ncbi_lineage_filled")


def tsv(path, columns, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        fh.write("\t".join(columns) + "\n")
        fh.writelines("\t".join(str(r.get(c, "-")) for c in columns) + "\n" for r in rows)


def toy_package(folder, sketches=True):
    """A small package: three species with a genome (RefSeq, GenBank, a transcriptome), one with rRNA records alone,
    two BUSCOs and a marker locus, rDNA records, organelle genes, metadata tables and two sketches."""
    def lineage(*names):
        return ";".join(names)
    species = [
        dict(taxid=4932, name="Saccharomyces cerevisiae", representative="GCA_000146045.2", source="model",
             supergroup="Obazoa", group="Fungi", ncbi_lineage_filled=lineage(
                 "Eukaryota", "Fungi", "Ascomycota", "Saccharomycetes", "Saccharomycetales", "Saccharomycetaceae",
                 "Saccharomyces", "Saccharomyces cerevisiae")),
        dict(taxid=4930, name="Saccharomyces paradoxus", representative="GCA_002079055.1", source="ncbi",
             supergroup="Obazoa", group="Fungi", ncbi_lineage_filled=lineage(
                 "Eukaryota", "Fungi", "Ascomycota", "Saccharomycetes", "Saccharomycetales", "Saccharomycetaceae",
                 "Saccharomyces", "Saccharomyces paradoxus")),
        dict(taxid=2850, name="Phaeodactylum tricornutum", representative="TSA_GAAB01", source="tsa",
             supergroup="TSAR", group="Stramenopiles", ncbi_lineage_filled=lineage(
                 "Eukaryota", "Sar", "Stramenopiles", "Bacillariophyceae", "Naviculales", "Phaeodactylaceae",
                 "Phaeodactylum", "Phaeodactylum tricornutum")),
        dict(taxid=5911, name="Tetrahymena thermophila", representative="-", source="rrna", supergroup="TSAR",
             group="Alveolata", ncbi_lineage_filled=lineage(
                 "Eukaryota", "Sar", "Ciliophora", "Oligohymenophorea", "Hymenostomatida", "",
                 "Tetrahymena", "Tetrahymena thermophila")),
    ]
    tsv(os.path.join(folder, "species.tsv"), SPECIES, species)
    genomes = [dict(accession="GCA_000146045.2", fetched="GCF_000146045.2", source="model", species_taxid=4932,
                    representative="t"),
               dict(accession="GCA_002079055.1", fetched="GCA_002079055.1", source="ncbi", species_taxid=4930,
                    representative="t"),
               dict(accession="GCA_000000001.1", fetched="GCA_000000001.1", source="assembly", species_taxid=4930,
                    representative="f"),
               dict(accession="TSA_GAAB01", fetched="TSA_GAAB01", source="tsa", species_taxid=2850,
                    representative="t")]
    tsv(os.path.join(folder, "genomes.tsv"), ("accession", "fetched", "source", "species_taxid", "representative"),
        genomes)
    tsv(os.path.join(folder, "busco", "loci.tsv"), ("locus", "name", "call_rate"),
        [dict(locus="1at2759", name="RPS3", call_rate=1.0), dict(locus="2at2759", name="EF", call_rate=0.333)])
    tsv(os.path.join(folder, "busco", "copies.tsv"), ("accession", "locus", "status", "confidence"),
        [dict(accession="GCA_000146045.2", locus="1at2759", status="Single", confidence="high"),
         dict(accession="GCA_002079055.1", locus="1at2759", status="Single", confidence="medium"),
         dict(accession="TSA_GAAB01", locus="1at2759", status="Fragmented", confidence="medium"),
         dict(accession="GCA_000146045.2", locus="2at2759", status="Duplicated", confidence="low"),
         dict(accession="GCA_002079055.1", locus="2at2759", status="Missing", confidence="-")])
    tsv(os.path.join(folder, "markers", "loci.tsv"), ("locus", "name", "call_rate"),
        [dict(locus="RPB1", name="RPB1", call_rate=0.66)])
    tsv(os.path.join(folder, "markers", "copies.tsv"), ("accession", "locus", "status", "confidence"),
        [dict(accession="GCA_000146045.2", locus="RPB1", status="Single", confidence="high")])
    for locus_set, locus in (("busco", "1at2759"), ("busco", "2at2759"), ("markers", "RPB1")):
        for sub, ext in (("copies", "fna"), ("cds", "fna"), ("proteins", "faa")):
            path = os.path.join(folder, locus_set, sub, f"{locus}.{ext}.gz")
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with util.bgzf_open(path) as fh:
                fh.write(f">GCA_000146045.2\n{dna(300) if ext == 'fna' else 'MKV'}\n")
    tsv(os.path.join(folder, "busco", "features.tsv"), ("accession", "locus", "kind", "start", "end"), [])
    tsv(os.path.join(folder, "rrna", "records.tsv"), ("record", "kind", "species_taxid", "source", "regions",
                                                       "confidence"),
        [dict(record="r1", kind="rDNA", species_taxid=4932, source="genome",
              regions="18S:1-1800,ITS1:1800-2000,5.8S:2000-2158,ITS2:2158-2400,28S:2400-5800", confidence="high"),
         dict(record="r2", kind="18S", species_taxid=5911, source="EUKARYOME", regions="-",  # an 18S span's
              confidence="medium")])
    with util.bgzf_open(os.path.join(folder, "rrna", "records.fna.gz")) as fh:
        fh.write(f">r1\n{dna(5800)}\n>r2\n{dna(1750)}\n")
    tsv(os.path.join(folder, "organelle", "genes.tsv"), ("record", "gene", "organelle", "species_taxid",
                                                         "confidence"),
        [dict(record="m1", gene="SSU", organelle="mitochondrion", species_taxid=4932, confidence="high"),
         dict(record="m2", gene="COX1", organelle="mitochondrion", species_taxid=4932, confidence="high"),
         dict(record="m3", gene="COX1", organelle="mitochondrion", species_taxid=4930, confidence="low"),
         dict(record="p1", gene="rbcL", organelle="plastid", species_taxid=2850, confidence="medium")])
    for name in ("genes.fna.gz", "molecules.fna.gz"):
        with util.bgzf_open(os.path.join(folder, "organelle", name)) as fh:
            fh.write(f">m1\n{dna(900)}\n")
    os.makedirs(os.path.join(folder, "metadata"))
    with gzip.open(os.path.join(folder, "metadata", "assembly.tsv.gz"), "wt") as fh:  # an older package's gzip
        fh.write("accession\tgenome_source\nGCA_000146045.2\tNCBI RefSeq\n")
    if sketches:
        with zipfile.ZipFile(os.path.join(folder, "metadata", "sketches.zip"), "w", zipfile.ZIP_STORED) as z:
            for acc in ("GCA_000146045.2", "GCA_002079055.1"):
                sig = [{"class": "sourmash_signature", "name": acc, "signatures": [{"ksize": 31, "mins": [1, 2]}]}]
                z.writestr(f"signatures/{acc}.sig.gz", gzip.compress(json.dumps(sig).encode()))
    os.makedirs(os.path.join(folder, "benchmarks"))
    with open(os.path.join(folder, "benchmarks", "summary.json"), "w") as fh:
        json.dump({"completeness": {"median": 0.8, "at_least_90": 0.1, "by_supergroup": {"Fungi": 0.9}}}, fh)
    with open(os.path.join(folder, package.MANIFEST), "w") as fh:
        json.dump({"format": 5, "version": "2026.11", "eukore": "0.2.2", "built": "2026-11-09T10:00:00",
                   "inputs": {}, "counts": {"genomes": 4}, "files": {}}, fh)
    return folder


def walk(node, path=()):
    """A taxonomy tree's nodes -> {path of names: node}."""
    out = {path + (node[0],): node}
    if isinstance(node[4], list):
        for child in node[4]:
            out.update(walk(child, path + (node[0],)))
    return out


class Release(unittest.TestCase):
    def test_a_species_path_collapses_and_fills(self):
        row = dict(supergroup="TSAR", group="Stramenopiles", name="X y",
                   ncbi_lineage_filled="Eukaryota;Sar;Stramenopiles;Bigyra;;Opalinidae;X;X y")
        self.assertEqual(release._path(row), [("supergroup", "TSAR"), ("group", "Stramenopiles"), ("class", "Bigyra"),
                                              ("order", "Bigyra_X"), ("family", "Opalinidae"), ("genus", "X"),
                                              ("species", "X y")])
        old = dict(supergroup="Obazoa", group="Fungi", name="A b", ncbi_lineage="Eukaryota;Fungi;Ascomycota")
        self.assertEqual(release._path(old)[2:4], [("phylum", "Ascomycota"), ("class", "Ascomycota_X")])

    def test_the_website_data_and_the_downloads(self):
        with tempfile.TemporaryDirectory() as tmp:
            pkg = toy_package(os.path.join(tmp, "Eukore_2026.11"))
            out = os.path.join(tmp, "release")
            said = []
            util.compression_threads(2)
            try:
                release.write(pkg, out, log=said.append)
            finally:
                util.compression_threads(1)
            with open(os.path.join(out, "data", "stats.json")) as fh:
                stats = json.load(fh)
            with open(os.path.join(out, "data", "taxonomy.json")) as fh:
                taxonomy = json.load(fh)
            with open(os.path.join(out, "data", "downloads.json")) as fh:
                listed = json.load(fh)
            # the species per genome source and supergroup
            self.assertEqual(stats["totals"]["species"], 4)
            self.assertEqual(stats["totals"]["species_with_genome"], 3)
            self.assertEqual(stats["sources"]["columns"], ["NCBI RefSeq", "NCBI GenBank", "ENA TSA",
                                                           "rRNA records only"])
            by_name = {s["name"]: s for s in stats["sources"]["supergroups"]}
            self.assertEqual(by_name["Obazoa"]["counts"], {"NCBI RefSeq": 1, "NCBI GenBank": 1, "ENA TSA": 0,
                                                            "rRNA records only": 0})
            self.assertEqual([g["name"] for g in by_name["TSAR"]["groups"]], ["Alveolata", "Stramenopiles"])
            self.assertEqual(by_name["TSAR"]["counts"]["rRNA records only"], 1)
            # the tree: counts, a name its parent's once, an empty rank filled
            nodes = walk(taxonomy["tree"])
            root = taxonomy["tree"]
            self.assertEqual(root[:4], ["Eukaryota", -1, 4, 3])
            fungi = nodes[("Eukaryota", "Obazoa", "Fungi")]
            self.assertEqual(fungi[1:4], [1, 2, 2])
            genus = nodes[("Eukaryota", "Obazoa", "Fungi", "Ascomycota", "Saccharomycetes", "Saccharomycetales",
                           "Saccharomycetaceae", "Saccharomyces")]
            self.assertEqual([c[0] for c in genus[4]], ["Saccharomyces cerevisiae", "Saccharomyces paradoxus"])
            self.assertEqual(genus[4][0][4], "GCA_000146045.2")
            self.assertIn(("Eukaryota", "TSAR", "Stramenopiles", "Bacillariophyceae"), nodes)
            self.assertIn(("Eukaryota", "TSAR", "Alveolata", "Ciliophora", "Oligohymenophorea", "Hymenostomatida",
                           "Hymenostomatida_X", "Tetrahymena", "Tetrahymena thermophila"), nodes)
            self.assertEqual(nodes[("Eukaryota", "TSAR", "Alveolata")][2:4], [1, 0])
            # the marker categories
            cats = {c["key"]: c for c in stats["markers"]}
            self.assertEqual(list(cats), ["busco", "markers", "rrna", "mitochondrion", "plastid"])
            self.assertEqual((cats["busco"]["loci"], cats["busco"]["sequences"], cats["busco"]["species"]), (2, 4, 3))
            self.assertEqual(cats["busco"]["confidence"], {"high": 1, "unverified": 0, "medium": 2, "low": 1})
            self.assertEqual(sum(i["sequences"] for i in cats["busco"]["items"]), 2)  # the loci by call rate
            self.assertEqual(cats["markers"]["items"], [dict(name="RPB1", locus="RPB1", sequences=1, species=1)])
            rrna_items = {i["name"]: i for i in cats["rrna"]["items"]}
            self.assertEqual((rrna_items["18S"]["sequences"], rrna_items["28S"]["sequences"]), (2, 1))
            self.assertEqual(rrna_items["whole units (18S to 28S)"]["sequences"], 1)
            self.assertEqual([i["name"] for i in cats["mitochondrion"]["items"]], ["SSU rRNA (12S)", "COX1"])
            self.assertEqual(cats["mitochondrion"]["items"][1], dict(name="COX1", sequences=2, species=2))
            self.assertEqual(stats["totals"]["marker_categories"], 5)
            self.assertEqual(stats["benchmarks"]["completeness"]["median"], 0.8)
            # the downloads, each BGZF, listed with their sizes and MD5s
            folder = os.path.join(out, "downloads", "Eukore_2026.11")
            self.assertEqual(listed["folder"], "downloads/Eukore_2026.11/")
            files = {e["key"]: e for e in listed["files"]}
            self.assertEqual(list(files), ["complete", "marker_genes", "marker_cds", "rrna_organelle", "sketches",
                                           "species", "genomes", "assembly"])
            for e in listed["files"]:
                path = os.path.join(folder, e["file"])
                blocks(path)
                self.assertEqual((e["bytes"], e["md5"]), (os.path.getsize(path), util.md5(path)))
            with open(os.path.join(folder, "MD5SUMS")) as fh:
                self.assertEqual(len(fh.read().splitlines()), len(files))

            def members(key):
                with tarfile.open(os.path.join(folder, files[key]["file"]), "r:gz") as tar:
                    return tar.getnames()
            complete = members("complete")
            self.assertEqual(complete[0], "Eukore_2026.11/manifest.json")
            self.assertEqual(len(complete), files["complete"]["files"])
            genes = members("marker_genes")
            self.assertIn("Eukore_2026.11/busco/copies/1at2759.fna.gz", genes)
            self.assertIn("Eukore_2026.11/markers/copies/RPB1.fna.gz", genes)
            self.assertIn("Eukore_2026.11/species.tsv", genes)
            self.assertFalse([m for m in genes if "/cds/" in m or "/proteins/" in m or "rrna" in m])
            cds = members("marker_cds")
            self.assertIn("Eukore_2026.11/busco/proteins/2at2759.faa.gz", cds)
            self.assertFalse([m for m in cds if "/copies/" in m or "features" in m])
            rrna = members("rrna_organelle")
            self.assertIn("Eukore_2026.11/rrna/records.fna.gz", rrna)
            self.assertIn("Eukore_2026.11/organelle/genes.tsv", rrna)
            self.assertNotIn("Eukore_2026.11/organelle/molecules.fna.gz", rrna)
            # unpacked side by side, the bundles are a package's folders as its reader reads them
            unpacked = os.path.join(tmp, "unpacked")
            with tarfile.open(os.path.join(folder, files["marker_cds"]["file"]), "r:gz") as tar:
                tar.extractall(unpacked, filter="data")
            self.assertEqual(package.read_copies(os.path.join(unpacked, "Eukore_2026.11"), "2at2759", "proteins"),
                             {"GCA_000146045.2": "MKV"})
            with gzip.open(os.path.join(folder, files["sketches"]["file"]), "rt") as fh:
                self.assertEqual([s["name"] for s in json.load(fh)], ["GCA_000146045.2", "GCA_002079055.1"])
            self.assertEqual(files["sketches"]["files"], 2)
            with gzip.open(os.path.join(folder, "Eukore_2026.11_assembly.tsv.gz"), "rt") as fh:  # gzip made BGZF
                self.assertEqual(fh.readline(), "accession\tgenome_source\n")
            with gzip.open(os.path.join(folder, "Eukore_2026.11_species.tsv.gz"), "rt") as fh:
                self.assertEqual(fh.read(), content(os.path.join(pkg, "species.tsv")).decode())
            self.assertTrue(any("Eukore_2026.11.tar.gz" in s for s in said))
            # a rerun keeps the downloads of the same package; another package's are made again
            stamp = os.path.getmtime(os.path.join(folder, "Eukore_2026.11.tar.gz"))
            said.clear()
            release.write(pkg, out, log=said.append)
            self.assertEqual(os.path.getmtime(os.path.join(folder, "Eukore_2026.11.tar.gz")), stamp)
            self.assertEqual(len(said), 1)  # the release's own line alone
            with open(os.path.join(pkg, package.MANIFEST), "a") as fh:
                fh.write(" ")
            os.remove(os.path.join(pkg, "metadata", "sketches.zip"))
            release.write(pkg, out, log=said.append)
            with open(os.path.join(out, "data", "downloads.json")) as fh:
                self.assertNotIn("sketches", [e["key"] for e in json.load(fh)["files"]])
            self.assertFalse(os.path.exists(os.path.join(folder, "Eukore_2026.11_sketches.sig.gz")))

    def test_no_release_of_an_unfinished_package(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit):
                release.write(os.path.join(tmp, "Eukore_2026.11"), os.path.join(tmp, "release"))

    def test_the_release_stage(self):
        import contextlib  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            toy_package(os.path.join(tmp, "out", "Eukore_2026.11"), sketches=False)
            printed = io.StringIO()
            with contextlib.redirect_stdout(printed):
                self.assertEqual(eukore_build.main(["--build_dir", os.path.join(tmp, "work"), "--out",
                                                    os.path.join(tmp, "out"), "--stage", "release",
                                                    "--version", "2026.11", "-t", "2"]), 0)
            self.assertIn("OUT/release/ of OUT/Eukore_2026.11/", printed.getvalue())  # the banner
            self.assertTrue(os.path.isfile(os.path.join(tmp, "out", "release", "downloads", "Eukore_2026.11",
                                                        "Eukore_2026.11.tar.gz")))
            self.assertFalse(os.path.exists(os.path.join(tmp, "work")))  # nothing of the build's own made


def listed(out):
    with open(os.path.join(out, "data", "downloads.json")) as fh:
        return json.load(fh)["files"]


class Downloads(unittest.TestCase):
    """The release's downloads since 2026-11-10: the same bytes from the same package whenever and on however many
    threads they are made, whole or made again, on the disk before their rename, their MD5s not taken twice."""

    def tearDown(self):
        util.compression_threads(1)

    def test_the_same_bytes_whatever_the_files_times_and_the_threads(self):
        import datetime  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            pkg = toy_package(os.path.join(tmp, "Eukore_2026.11"))
            built = int(datetime.datetime(2026, 11, 9, 10, tzinfo=datetime.timezone.utc).timestamp())
            self.assertEqual(release.built_time(pkg), built)  # the manifest's `built`, read as UTC on every machine
            one, many = os.path.join(tmp, "one"), os.path.join(tmp, "many")
            release.write(pkg, one, log=lambda line: None)
            # the package's files copied or written again (other mtimes), the downloads made four at once, their
            # blocks compressed on three threads
            for folder, _dirs, names in os.walk(pkg):
                for name in names:
                    os.utime(os.path.join(folder, name), (1_000_000_000, 1_000_000_000))
            util.compression_threads(3)
            said = []
            release.write(pkg, many, log=said.append, threads=4)
            self.assertEqual([e["file"] for e in listed(many)], [e["file"] for e in listed(one)])  # the plan's order
            self.assertEqual([(e["md5"], e["bytes"], e["files"]) for e in listed(many)],
                             [(e["md5"], e["bytes"], e["files"]) for e in listed(one)])
            self.assertEqual(len(said), len(listed(many)) + 1)
            folders = [os.path.join(out, "downloads", "Eukore_2026.11") for out in (one, many)]
            self.assertEqual(sorted(os.listdir(folders[0])), sorted(os.listdir(folders[1])))
            for name in os.listdir(folders[0]):
                if name != release.STATE:  # the state: the files' own mtimes
                    self.assertEqual(content(os.path.join(folders[0], name)), content(os.path.join(folders[1], name)),
                                     name)
            for e in listed(many):
                self.assertEqual(e["md5"], util.md5(os.path.join(folders[1], e["file"])))
            with tarfile.open(os.path.join(folders[1], "Eukore_2026.11.tar.gz"), "r:gz") as tar:
                self.assertEqual({(m.mtime, m.uid, m.gid, m.uname, m.gname) for m in tar.getmembers()},
                                 {(built, 0, 0, "", "")})
            # no build time in the manifest: one fixed time
            with open(os.path.join(pkg, package.MANIFEST), "w") as fh:
                json.dump({"format": 5, "version": "2026.11", "files": {}}, fh)
            self.assertEqual(release.built_time(pkg), release.TAR_EPOCH)
            self.assertEqual(release.built_time(os.path.join(tmp, "none")), release.TAR_EPOCH)

    def test_a_cut_download_is_made_again_and_a_kept_one_not_read_again(self):
        from unittest import mock  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            pkg = toy_package(os.path.join(tmp, "Eukore_2026.11"))
            out = os.path.join(tmp, "release")
            with mock.patch.object(util, "publish", wraps=util.publish) as published:
                release.write(pkg, out, log=lambda line: None)
            first = listed(out)
            folder = os.path.join(out, "downloads", "Eukore_2026.11")
            # every file of the release went into place through util.publish: on the disk before its rename
            self.assertEqual(sorted(os.path.basename(c.args[1]) for c in published.call_args_list),
                             sorted([e["file"] for e in first] + ["MD5SUMS", "stats.json", "taxonomy.json",
                                                                  "downloads.json"] + 2 * [release.STATE]))
            # a node that went down: a bundle cut in a block, a table without its end block
            cut = {"Eukore_2026.11.tar.gz": lambda data: data[:len(data) // 2],
                   "Eukore_2026.11_species.tsv.gz": lambda data: data[:-len(util.BGZF_EOF)]}
            for name, short in cut.items():
                data = content(os.path.join(folder, name))
                with open(os.path.join(folder, name), "wb") as fh:
                    fh.write(short(data))
            kept = os.path.getmtime(os.path.join(folder, "Eukore_2026.11_marker_cds.tar.gz"))
            said = []
            release.write(pkg, out, log=said.append)
            self.assertEqual(sorted(s for s in said if "cut short" in s),
                             [f"release: {name} of an earlier run is cut short (no BGZF end block): made again"
                              for name in sorted(cut)])
            self.assertEqual(len(said), 5)  # each said and made, and the release's own line
            self.assertEqual(listed(out), first)  # the same bytes again
            for e in first:
                blocks(os.path.join(folder, e["file"]))
                self.assertEqual(e["md5"], util.md5(os.path.join(folder, e["file"])))
            self.assertEqual(os.path.getmtime(os.path.join(folder, "Eukore_2026.11_marker_cds.tar.gz")), kept)
            # a rerun: no download read again for its MD5, the sketches not loaded to count them
            sketches = os.path.join(folder, "Eukore_2026.11_sketches.sig.gz")
            with mock.patch.object(util, "md5", wraps=util.md5) as md5, \
                    mock.patch.object(gzip, "open", wraps=gzip.open) as opened:
                release.write(pkg, out, log=said.append)
            self.assertEqual([c.args[0] for c in md5.call_args_list], [os.path.join(pkg, package.MANIFEST)])
            self.assertNotIn(sketches, [c.args[0] for c in opened.call_args_list])
            self.assertEqual(listed(out), first)
            # but a file changed since (its mtime) is: its state is not its own any more
            os.utime(sketches, (1_000_000_000, 1_000_000_000))
            with mock.patch.object(util, "md5", wraps=util.md5) as md5, \
                    mock.patch.object(gzip, "open", wraps=gzip.open) as opened:
                release.write(pkg, out, log=said.append)
            self.assertEqual([c.args[0] for c in md5.call_args_list], [os.path.join(pkg, package.MANIFEST), sketches])
            self.assertIn(sketches, [c.args[0] for c in opened.call_args_list])
            self.assertEqual(listed(out), first)
            # an earlier Eukore's state (the manifest's MD5 alone): the downloads kept, measured once more
            with open(os.path.join(folder, release.STATE)) as fh:
                state = json.load(fh)
            with open(os.path.join(folder, release.STATE), "w") as fh:
                json.dump({"manifest_md5": state["manifest_md5"], "package": state["package"]}, fh)
            said.clear()
            release.write(pkg, out, log=said.append)
            self.assertEqual((len(said), listed(out)), (1, first))

    def test_a_cut_package_table_is_no_download(self):
        with tempfile.TemporaryDirectory() as tmp:
            pkg = toy_package(os.path.join(tmp, "Eukore_2026.11"))
            table = os.path.join(pkg, "metadata", "ncbi.tsv.gz")
            with util.bgzf_open(table) as fh:
                fh.write("genome\tfetched\n" + "GCA_000146045.2\tt\n" * 20000)
            data = content(table)
            with open(table, "wb") as fh:
                fh.write(data[:-len(util.BGZF_EOF)])  # the copy of it broke off at a block's end: gzip reads it
            out = os.path.join(tmp, "release")
            with self.assertRaises(OSError) as failed:
                release.write(pkg, out, log=lambda line: None)
            self.assertIn("ncbi.tsv.gz is cut short (no BGZF end block)", str(failed.exception))
            folder = os.path.join(out, "downloads", "Eukore_2026.11")
            self.assertNotIn("Eukore_2026.11_ncbi.tsv.gz", os.listdir(folder))
            # the website's data is written after the downloads: none of a release whose downloads failed
            self.assertEqual(os.listdir(os.path.join(out, "data")), [])
            with open(table, "wb") as fh:
                fh.write(data)
            release.write(pkg, out, log=lambda line: None)
            self.assertEqual(content(os.path.join(folder, "Eukore_2026.11_ncbi.tsv.gz")), data)
            self.assertEqual(sorted(os.listdir(os.path.join(out, "data"))),
                             ["downloads.json", "stats.json", "taxonomy.json"])

    def test_the_partials_of_killed_runs_are_cleared(self):
        import socket  # noqa: PLC0415
        import subprocess  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            pkg = toy_package(os.path.join(tmp, "Eukore_2026.11"))
            out = os.path.join(tmp, "release")
            release.write(pkg, out, log=lambda line: None)
            folder, data = os.path.join(out, "downloads", "Eukore_2026.11"), os.path.join(out, "data")
            before = sorted(os.listdir(folder)), sorted(os.listdir(data))
            gone = subprocess.Popen(["true"])
            gone.wait()
            dead, live = f".{socket.gethostname()}-{gone.pid}", ".elsewhere-7"  # a killed run's; another host's, now
            for where, name in ((folder, "Eukore_2026.11.tar.gz"), (folder, "MD5SUMS"), (data, "stats.json")):
                for whose in (dead, live):
                    with open(os.path.join(where, f"{whose}.{name}.partial"), "w") as fh:
                        fh.write("half")
            release.write(pkg, out, log=lambda line: None)
            self.assertEqual(sorted(os.listdir(folder)), sorted(
                before[0] + [f"{live}.Eukore_2026.11.tar.gz.partial", f"{live}.MD5SUMS.partial"]))
            self.assertEqual(sorted(os.listdir(data)), sorted(before[1] + [f"{live}.stats.json.partial"]))


if __name__ == "__main__":
    unittest.main()
