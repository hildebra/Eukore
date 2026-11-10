#!/usr/bin/env python3
"""Tests of the organelle genes' models (eukore/organelle.py; since 2026-11-10, the audit's B8 and B2): a
kinetoplastid maxicircle's 12S as its LSU and its 9S as its SSU (rrna_gene, rrna_naming, RefSeq's parsed genes);
miniprot's whole model in place of nhmmer's single hit that is not the gene's coding sequence (with_miniprot: one
exon of a gene with introns, an end reached past the gene's stop), by rows as the searches give them and through an
assembly's search (stand-in nhmmer, miniprot and barrnap).

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, ncbi as ncbi_mod, organelle  # noqa: E402
import fixtures as fx  # noqa: E402
from test_organelle import feature, gbff_record, random_seq  # noqa: E402

YEAST = "GCA_000146045.2"
RELATIVE = "NC_000028.1"  # Kluyveromyces lactis' mitochondrion: RefSeq's, the yeast's queries

NHMMER = r'''#!/usr/bin/env python3
"""nhmmer's stand-in: of each query one hit, its longest stretch from its first base (in steps of 3 bases, 60 at
least) that lies on the target's forward strand: a gene's first exon alone where an intron follows, the gene as far
as the query's where it ends before the query does."""
import sys
args = sys.argv[1:]
tbl = args[args.index("--tblout") + 1]


def read(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        else:
            out[name] += line.strip()
    return out


queries, targets = read(args[-2]), read(args[-1])
with open(tbl, "w") as fh:
    fh.write("# target name\n")
    for q, s in queries.items():
        for t, seq in targets.items():
            k = len(s)
            while k >= 60 and s[:k] not in seq:
                k -= 3
            if k >= 60:
                at = seq.find(s[:k]) + 1
                fh.write(f"{t} - {q} - 1 {k} {at} {at + k - 1} {at} {at + k - 1} {len(seq)} + 1e-50 {k}.0 0 -\n")
'''
BARRNAP = "#!/bin/sh\necho '##gff-version 3'\n"  # no rRNA: the protein-coding genes alone


def coding(seed):
    """A gene's coding sequence of 363 bases, ATG to TAA, no stop inside by any code (fixtures' codons)."""
    return fx.cds_of(fx.protein(seed))


def alignment(n, query, contig, exons, first, last, score=900):
    """miniprot's GFF lines of one alignment on the forward strand: exons [(start, end)], the last with its stop
    codon; the protein's residues first..last aligned (of 120)."""
    a, b = exons[0][0], exons[-1][1]
    lines = [f"##PAF\t{query}\t120", "##STA\tM",
             f"{contig}\tminiprot\tmRNA\t{a}\t{b}\t{score}\t+\t.\tID=MP{n};Rank=1;Identity=0.9500;"
             f"Target={query} {first} {last}"]
    for x, y in exons:
        lines.append(f"{contig}\tminiprot\tCDS\t{x}\t{y - 3 if y == b else y}\t{score}\t+\t0\tParent=MP{n};Rank=1;"
                     f"Identity=0.95;Target={query} {first} {last}")
    lines.append(f"{contig}\tminiprot\tstop_codon\t{b - 2}\t{b}\t{score}\t+\t0\tParent=MP{n};Rank=1")
    return lines

# Leishmania tarentolae under Kinetoplastea, beside the fixtures' taxa (fx.NODES holds no kinetoplastid)
KINETOPLASTIDS = [(5653, 2759, "class", "Kinetoplastea"), (5658, 5653, "genus", "Leishmania"),
                  (5689, 5658, "species", "Leishmania tarentolae")]
MAXICIRCLE = random_seq(3000, 31)
# RefSeq's NC_000894 (Leishmania tarentolae's maxicircle): its rRNA features as it names them, the 12S first
MAXICIRCLE_RRNA = (("438..1610", "12S"), ("1639..2249", "9S"))


def taxonomy(tmp):
    with mock.patch.object(fx, "NODES", fx.NODES + KINETOPLASTIDS):
        return ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))


class KinetoplastidTest(unittest.TestCase):
    """The audit's B8: "12S" was read as the small subunit's rRNA in every lineage; a kinetoplastid's is the LSU."""

    def test_rrna_named_by_lineage(self):
        k = organelle.KINETOPLASTID
        self.assertEqual([organelle.rrna_gene(t, k) for t in ("12s rrna 12s ribosomal rna",
                                                              "9s rrna 9s ribosomal rna")], ["LSU", "SSU"])
        # a feature that names its subunit is read by it, whatever its number; 5S and smaller ones as elsewhere
        self.assertEqual([organelle.rrna_gene(t, k) for t in ("12s small subunit ribosomal rna", "rnl 9s", "rrn5",
                                                              "4.5s ribosomal rna")], ["SSU", "LSU", "5S", None])
        # elsewhere as before: 12S and 9S a small subunit's, 16S an animal mitochondrion's LSU alone
        self.assertEqual([organelle.rrna_gene("12s ribosomal rna", how) for how in (True, False)], ["SSU", "SSU"])
        self.assertEqual([organelle.rrna_gene("16s ribosomal rna", how) for how in (True, False, k)],
                         ["LSU", "SSU", "SSU"])
        self.assertEqual(organelle.rrna_gene("9s ribosomal rna"), "SSU")
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = taxonomy(tmp)
        self.assertEqual([organelle.rrna_naming(ncbi, t) for t in (5689, 5653, 9606, 4932, None)],
                         [k, k, True, False, False])

    def test_a_maxicircles_record(self):
        features = [("rRNA", [location], [("gene", f"{name} rRNA"), ("product", f"{name} ribosomal RNA")])
                    for location, name in MAXICIRCLE_RRNA]
        rec = dict(seqid="NC_7.1", taxid=5689, length=len(MAXICIRCLE), circular=True, sequence=MAXICIRCLE,
                   features=features)
        # record_genes hands its `animal` on to rrna_gene for a mitochondrion's features: rrna_naming's value
        got = organelle.record_genes(rec, "mitochondrion", organelle.KINETOPLASTID)
        self.assertEqual({g: (r["start"], r["end"], r["copies"]) for g, r in got.items()},
                         {"LSU": (438, 1610, 1), "SSU": (1639, 2249, 1)})
        self.assertEqual(got["SSU"]["sequence"], MAXICIRCLE[1638:2249])
        # until 2026-11-10 (the record read as any protist's): the 12S its "SSU", the longer of two, and no LSU
        old = organelle.record_genes(rec, "mitochondrion", False)
        self.assertEqual({g: (r["start"], r["end"], r["copies"]) for g, r in old.items()}, {"SSU": (438, 1610, 2)})

    def test_refseq_parsed_by_lineage(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = taxonomy(tmp)
            src = os.path.join(tmp, "sources")
            folder = organelle.folder_of(src)
            os.makedirs(folder)
            rrnas = [feature("rRNA", location, f'/gene="{name} rRNA"', f'/product="{name} ribosomal RNA"')
                     for location, name in MAXICIRCLE_RRNA]
            with gzip.open(os.path.join(folder, "mitochondrion.1.genomic.gbff.gz"), "wt") as fh:
                fh.write(gbff_record("NC_000894.1", MAXICIRCLE, 5689, "mitochondrion:kinetoplast", rrnas))
                # an animal's beside it: its 12S the SSU, its 16S the LSU, as before
                fh.write(gbff_record("NC_000001.1", MAXICIRCLE, 9606, "mitochondrion", [
                    feature("rRNA", "101..500", '/product="12S ribosomal RNA"'),
                    feature("rRNA", "601..1500", '/product="16S ribosomal RNA"')]))
            with open(os.path.join(folder, "release.txt"), "w") as fh:
                fh.write("mitochondrion.1.genomic.gbff.gz\n")
            rows = organelle.refseq_genes(src, os.path.join(tmp, "work"), ncbi, lambda _t: 0)
            self.assertEqual({(r["seqid"], r["gene"]): (r["start"], r["end"]) for r in rows},
                             {("NC_000894.1", "LSU"): ("438", "1610"), ("NC_000894.1", "SSU"): ("1639", "2249"),
                              ("NC_000001.1", "SSU"): ("101", "500"), ("NC_000001.1", "LSU"): ("601", "1500")})
            # the species' genes, and the query of a relative's assembly: its SSU is the 9S
            near = organelle.Nearest(organelle.choose(rows, ncbi), ncbi)
            queries = dict(near.queries(5689, "mitochondrion"))
            self.assertEqual(queries["SSU|NC_000894.1"], MAXICIRCLE[1638:2249])
            self.assertEqual(queries["LSU|NC_000894.1"], MAXICIRCLE[437:1610])


def hit(gene, molecule, start, end, cover):
    """A gene's row as nhmmer_genes gives its single hit: the hit's span reached to the query's ends (gene_span)."""
    return dict(gene=gene, organelle="mitochondrion", seqid="m1", start=start, end=end, strand="+", exons=1, copies=1,
                partial=cover < organelle.MIN_COVER, sequence=molecule[start - 1:end], parts=f"{start}-{end}:+",
                model="nhmmer", cover=cover)


def model(gene, molecule, exons, partial=False):
    """A gene's row as miniprot_genes gives it: its exons [(start, end)] on the forward strand."""
    start, end = exons[0][0], exons[-1][1]
    return dict(gene=gene, organelle="mitochondrion", seqid="m1", start=start, end=end, strand="+", exons=len(exons),
                copies=1, partial=partial, sequence=molecule[start - 1:end],
                parts=organelle.parts_text([(a, b, "+") for a, b in exons]), model="miniprot")


def taken(nhmmer, miniprot, table=1):
    return organelle.with_miniprot({nhmmer["gene"]: nhmmer}, {miniprot["gene"]: miniprot}, table)[nhmmer["gene"]]


class ModelsTest(unittest.TestCase):
    """The audit's B2: nhmmer's single hit was kept whatever it was - one exon of a gene with introns, reached to the
    query's length into the intron and shipped as the coding sequence."""

    def test_a_single_hit_that_is_not_the_coding_sequence(self):
        cds, flank = coding(41), random_seq(200, 43)
        intron = ("C" if cds[180] != "C" else "G") + random_seq(299, 42)
        split = flank + cds[:180] + intron + cds[180:] + flank  # exons 201-380 and 681-863
        # nhmmer hits the first exon alone (the second below its threshold): half the query, the span reached by the
        # query's other half - 183 bases of the intron. miniprot aligned the whole protein: its model
        one_exon, whole = hit("COX1", split, 201, 563, 180 / 363), model("COX1", split, [(201, 380), (681, 863)])
        self.assertTrue(one_exon["partial"])
        got = taken(one_exon, whole)
        self.assertEqual((got["model"], got["exons"], organelle.coding_of(got, got["sequence"])), ("miniprot", 2, cds))
        # miniprot's model of a part of the protein never replaces nhmmer's hit, as before
        self.assertEqual(taken(one_exon, dict(whole, partial=True))["model"], "nhmmer")
        # by the share of the query the hit covers alone (no stop in either): under NHMMER_WHOLE miniprot's model
        gene = flank + cds + flank
        shorter = model("CYTB", gene, [(210, 563)])  # miniprot's: three codons short of the gene's start
        for cover, by in ((organelle.NHMMER_WHOLE - 0.01, "miniprot"), (organelle.NHMMER_WHOLE, "nhmmer"),
                          (1.0, "nhmmer")):
            self.assertEqual(taken(hit("CYTB", gene, 201, 563, cover), shorter)["model"], by, cover)

    def test_a_single_hit_with_more_stops_than_miniprots_model(self):
        cds, flank = coding(41), random_seq(200, 43)
        # the gene ends 9 bases before the query does: the hit's end is reached past the gene's stop, which is then
        # a stop inside the coding sequence; miniprot's model ends at it
        short = flank + cds[:351] + "TAA" + random_seq(200, 44)
        past, ends = hit("nad1", short, 201, 563, 351 / 363), model("nad1", short, [(201, 554)])
        self.assertEqual((organelle.stops_of(past, 1) > 0, organelle.stops_of(ends, 1)), (True, 0))
        self.assertEqual(taken(past, ends)["model"], "miniprot")
        # as many stops in miniprot's model (a stop that RNA editing removes: in the gene, whoever finds it): kept
        edited = flank + cds[:150] + "TAA" + cds[153:] + flank
        self.assertEqual(taken(hit("nad5", edited, 201, 563, 1.0), model("nad5", edited, [(201, 563)]))["model"],
                         "nhmmer")
        # by the organelle's code: TGA is tryptophan in a yeast's mitochondrion (table 3), a stop by the standard one
        tga = flank + cds[:150] + "TGA" + cds[153:] + flank
        whole, around = hit("COX3", tga, 201, 563, 1.0), model("COX3", tga, [(201, 350), (354, 563)])
        self.assertEqual([organelle.stops_of(whole, t) for t in (3, 1)], [0, 1])
        self.assertEqual([taken(whole, around, t)["model"] for t in (3, 1)], ["nhmmer", "miniprot"])
        # no coding sequence to read (a gene across the origin: start after end; a row without its sequence, an
        # rRNA): kept, as before
        gene = flank + cds + flank
        across = dict(hit("atp6", gene, 500, 100, 1.0))
        self.assertIsNone(organelle.stops_of(across, 1))
        self.assertEqual(taken(across, model("atp6", gene, [(201, 563)]))["model"], "nhmmer")
        self.assertIsNone(organelle.stops_of(dict(gene="CYTB", exons=1, model="nhmmer"), 1))
        self.assertIsNone(organelle.stops_of(dict(hit("SSU", gene, 201, 563, 1.0)), 1))

    def test_through_an_assemblys_search(self):
        def laid(pieces):  # [(name or None, sequence)] -> the molecule, {name: (start, end)} 1-based
            seq, at = "", {}
            for name, s in pieces:
                if name:
                    at[name] = (len(seq) + 1, len(seq) + len(s))
                seq += s
            return seq, at
        cox1, cytb, nad1 = coding(51), coding(52), coding(53)
        gaps = [random_seq(150, s) for s in (61, 62, 63, 64)]
        # RefSeq's relative: the three genes whole, one exon each. The yeast's mitochondrion: its COX1 in two exons,
        # its CYTB as the relative's, its nad1 ending 9 bases before the relative's
        relative, r = laid([(None, gaps[0]), ("COX1", cox1), (None, gaps[1]), ("CYTB", cytb), (None, gaps[2]),
                            ("nad1", nad1), (None, gaps[3])])
        intron = ("C" if cox1[180] != "C" else "G") + random_seq(299, 66)
        yeast, y = laid([(None, gaps[3]), ("COX1a", cox1[:180]), (None, intron), ("COX1b", cox1[180:]),
                         (None, gaps[2]), ("CYTB", cytb), (None, gaps[1]), ("nad1", nad1[:351] + "TAA"),
                         (None, gaps[0])])
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            src, work = os.path.join(tmp, "sources"), os.path.join(tmp, "work", "organelle")
            folder = organelle.folder_of(src)
            os.makedirs(folder)
            with gzip.open(os.path.join(folder, "mitochondrion.1.genomic.gbff.gz"), "wt") as fh:
                fh.write(gbff_record(RELATIVE, relative, 28985, "mitochondrion",
                                     [feature("CDS", f"{a}..{b}", f'/gene="{g}"') for g, (a, b) in r.items()]))
            with open(os.path.join(folder, "release.txt"), "w") as fh:
                fh.write("mitochondrion.1.genomic.gbff.gz\n")
            chosen = organelle.choose(organelle.refseq_genes(src, work, ncbi, lambda _t: 0), ncbi)
            genomes_dir = os.path.join(tmp, "genomes")
            os.makedirs(genomes_dir)
            with gzip.open(genomes.genome_file(genomes_dir, YEAST), "wt") as fh:
                fh.write(f">c1 Saccharomyces cerevisiae chromosome I\n{random_seq(500, 67)}\n"
                         f">m1 Saccharomyces cerevisiae mitochondrion, complete genome\n{yeast}\n")
            jobs = {YEAST: dict(accession=YEAST, species_taxid=4932, taxid=559292)}
            barrnap, nhmmer = fx.tool(tmp, "barrnap", BARRNAP), fx.tool(tmp, "nhmmer", NHMMER)
            # miniprot's alignments: COX1's two exons, CYTB three codons short of its start, nad1 to its own stop
            gff = ["##gff-version 3"]
            gff += alignment(1, f"COX1~{RELATIVE}", "m1", [y["COX1a"], y["COX1b"]], 1, 120)
            gff += alignment(2, f"CYTB~{RELATIVE}", "m1", [(y["CYTB"][0] + 9, y["CYTB"][1])], 4, 120)
            gff += alignment(3, f"nad1~{RELATIVE}", "m1", [y["nad1"]], 1, 117)
            miniprot = fx.miniprot(tmp, "\n".join(gff) + "\n")
            # nhmmer's rows: COX1's first exon reached to the query's length (into the intron), nad1's hit past its
            # stop; until 2026-11-10 both were the genes' records and coding sequences
            queries = organelle.Nearest(chosen, ncbi).queries(4932, "mitochondrion")
            alone = organelle.nhmmer_genes([("m1", "mitochondrion", yeast)], queries, nhmmer)
            a1, n1 = y["COX1a"][0], y["nad1"][0]
            self.assertEqual({g: (h["start"], h["end"], h["exons"], h["partial"], round(h["cover"], 3))
                              for g, h in alone.items()},
                             {"COX1": (a1, a1 + 362, 1, True, 0.496), "CYTB": y["CYTB"] + (1, False, 1.0),
                              "nad1": (n1, n1 + 362, 1, False, 0.967)})

            def run():
                return organelle.assembly_pass(work, jobs, lambda a: genomes.genome_file(genomes_dir, a), chosen,
                                               ncbi, barrnap, nhmmer, 1, 1, lambda _t: 0, miniprot=miniprot)
            self.assertEqual((run(), run()), ((1, 0), (0, 0)))
            rows = {r["gene"]: r for r in organelle.assembly_rows(organelle.assembly_file(work, YEAST))}
            self.assertEqual({g: (r["model"], int(r["start"]), int(r["end"]), r["exons"], r["partial"])
                              for g, r in rows.items()},
                             {"COX1": ("miniprot", a1, y["COX1b"][1], "2", "f"),
                              "CYTB": ("nhmmer",) + y["CYTB"] + ("1", "f"),  # whole and without a stop: its ends kept
                              "nad1": ("miniprot",) + y["nad1"] + ("1", "f")})
            self.assertEqual(rows["COX1"]["parts"], "{}-{}:+,{}-{}:+".format(*y["COX1a"], *y["COX1b"]))
            # into the package: the coding sequences are the genes', read by the yeast's mitochondrial code
            records = organelle.check(organelle.collect(ncbi, jobs, {4932}, chosen, work), ncbi)
            self.assertEqual({r["gene"]: (r["cds"], r["internal_stops"]) for _n, _s, r in records},
                             {"COX1": (cox1, 0), "CYTB": (cytb, 0), "nad1": (nad1[:351] + "TAA", 0)})
            # the threshold is one of the search's inputs: another one, searched again
            with mock.patch.object(organelle, "NHMMER_WHOLE", organelle.NHMMER_WHOLE + 0.01):
                self.assertEqual(run(), (1, 0))


if __name__ == "__main__":
    unittest.main()
