#!/usr/bin/env python3
"""Tests of format 7 (since 2026-11-14): EUKARYOME's taxonomy as one tree and every taxon placed in it
(eukore/eukaryome.py), labels reaching above the genus (rrna.label_taxon), EUKARYOME's long reads and SSU set as the
only published 18S, their audit (eukore/audit.py: ends, flanks; eukore/chimeras.py: KSGP 4.5's window screen, 4.2's
UCHIME), the links to PR2's and KSGP's records (eukore/links.py), the telomere check, and a package of them: records
and organelle genes with their rank and two lineages, species of organelle genes alone, taxonomy/eukaryome.tsv.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import json
import os
import random
import shutil
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import audit, chimeras, eukaryome, labels, links, ncbi as ncbi_mod, package, rrna, screen  # noqa: E402
import fixtures as fx  # noqa: E402

HEADER = ("superkingdom", "kingdom", "subreg", "phylum", "subphy", "class", "order", "family", "gen", "kgd-specific")
DIKARYA = ("Obazoa", "Fungi", "Dikarya", "Ascomycota", "Saccharomycotina", "Saccharomycetes", "Saccharomycetales")
CLASSIFICATION = [
    HEADER,
    DIKARYA + ("Saccharomycetaceae", "Saccharomyces", ""),
    DIKARYA + ("Saccharomycetaceae", "Kazachstania", ""),
    DIKARYA + ("Debaryomycetaceae", "Candida", ""),
    # Ascomycota under another subkingdom in one row: the parent of most rows taken
    ("Obazoa", "Fungi", "Fungi.subreg.incertae.sedis", "Ascomycota", "Saccharomycotina", "Saccharomycetes",
     "Saccharomycetales", "Saccharomycetaceae", "Lachancea", ""),
    ("Obazoa", "Metazoa", "", "Chordata", "", "Mammalia", "Primates", "Hominidae", "Homo", ""),
    # a genus of two kingdoms: the plant's written with its kingdom (EUKARYOME's kingdom-specific column)
    ("Archaeplastida", "Viridiplantae", "", "Tracheophyta", "Tracheophytina", "Angiospermae", "Rosales", "Moraceae",
     "Morus", "Morus(Viridiplantae)"),
    ("Obazoa", "Metazoa", "", "Chordata", "", "Aves", "Suliformes", "Sulidae", "Morus", ""),
    HEADER,  # the header again inside the table
    # a genus named as its family: no genus
    ("Archaeplastida", "Viridiplantae", "", "Tracheophyta", "Tracheophytina", "Angiospermae", "Solanales",
     "Solanaceae", "Solanaceae", ""),
    ("Obazoa", "Metazoa", "", "Nematoda", "", "Chromadorea", "Rhabditida", "Toxocaridae", "Toxocara", ""),
    ("Obazoa", "Metazoa", "", "Nematoda", "", "Chromadorea", "Rhabditida", "Anisakidae", "Anisakis", ""),
    # a row a column off (its phylum in the subkingdom's column)
    ("Obazoa", "Metazoa", "Nematoda", "Chromadorea", "Rhabditida", "Ascarididae", "Ascaris", "", "", ""),
]
MORE_NODES = [(33090, 2759, "kingdom", "Viridiplantae"), (3497, 33090, "genus", "Morus"),
              (3498, 3497, "species", "Morus alba"), (2, 1, "superkingdom", "Bacteria"),
              (265, 2, "genus", "Paracoccus"), (9900001, 4930, "species", "Saccharomyces novus")]


def ncbi_world(tmp):
    with mock.patch.object(fx, "NODES", fx.NODES + MORE_NODES):
        return ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))


def tree_of():
    tree = eukaryome.Tree.from_classification(CLASSIFICATION)
    for label in ("k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;"
                  "g__Saccharomyces;s__cerevisiae",
                  "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Candida;"
                  "s__sphaerica"):
        tree.place_label(label)
    return tree


class TreeTest(unittest.TestCase):
    def test_the_classification_one_tree(self):
        tree = eukaryome.Tree.from_classification(CLASSIFICATION)
        self.assertEqual(tree.text((8, "Saccharomyces")),
                         ";".join(DIKARYA) + ";Saccharomycetaceae;Saccharomyces;")
        self.assertIn(";Dikarya;", tree.text((8, "Lachancea")))  # its Ascomycota the one of most rows
        # an empty rank filled as EUKARYOME fills its own
        self.assertEqual(tree.text((8, "Homo")), "Obazoa;Metazoa;Metazoa.subreg.incertae.sedis;Chordata;"
                                                 "Chordata.subphy.incertae.sedis;Mammalia;Primates;Hominidae;Homo;")
        # the homonym genus of two kingdoms: two taxa, the plant's by its kingdom-specific name
        self.assertIn("Viridiplantae", tree.text((8, "Morus(Viridiplantae)")))
        self.assertIn(";Aves;", tree.text((8, "Morus")))
        self.assertNotIn((8, "Solanaceae"), tree.parent)  # a genus named as its family: none
        self.assertIn((7, "Solanaceae"), tree.parent)
        self.assertFalse([n for n in tree.parent if n[1] in HEADER])  # the header again: no taxa
        # the row a column off put back: Ascaris a genus of a family of Rhabditida, Nematoda a phylum alone
        self.assertEqual(tree.text((8, "Ascaris")), "Obazoa;Metazoa;Metazoa.subreg.incertae.sedis;Nematoda;"
                                                    "Nematoda.subphy.incertae.sedis;Chromadorea;Rhabditida;"
                                                    "Ascarididae;Ascaris;")
        self.assertEqual({n[0] for n in tree.by_name["Nematoda"]}, {3})
        # every name under one lineage
        for name, nodes in tree.by_name.items():
            for a in nodes:
                for b in nodes:
                    self.assertTrue(a in tree.lineage(b) or b in tree.lineage(a), name)

    def test_labels_placed(self):
        tree = eukaryome.Tree.from_classification(CLASSIFICATION)
        base = "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;"
        self.assertEqual(tree.place_label(base + "f__Saccharomycetaceae;g__Saccharomyces;s__cerevisiae"),
                         (9, "Saccharomyces_cerevisiae"))
        # a family the table lacks (EUKARYOME's own Saccharomycetales.fam08): added under its order
        node = tree.place_label(base + "f__Saccharomycetales.fam08;g__unclassified;s__unclassified")
        self.assertEqual((node, tree.parent[node]), ((7, "Saccharomycetales.fam08"), (6, "Saccharomycetales")))
        # a label of a homonym genus takes its own kingdom's
        node = tree.place_label("k__Viridiplantae;p__Tracheophyta;c__Angiospermae;o__Rosales;f__Moraceae;g__Morus;"
                                "s__alba")
        self.assertEqual(tree.parent[node], (8, "Morus(Viridiplantae)"))
        self.assertEqual(tree.place_label("k__Metazoa;p__Nematoda;c__unclassified;o__unclassified;f__unclassified;"
                                          "g__unclassified;s__unclassified"), (3, "Nematoda"))
        self.assertIsNone(tree.place_label("k__unclassified;p__unclassified"))

    def test_ncbi_taxa_placed(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_world(tmp)
        tree = tree_of()
        self.assertEqual(tree.place(ncbi, 4932), ((9, "Saccharomyces_cerevisiae"), "species"))
        self.assertEqual(tree.place(ncbi, 559292)[0], (9, "Saccharomyces_cerevisiae"))  # a strain: its species
        self.assertEqual(tree.place(ncbi, 28985), ((9, "Candida_sphaerica"), "synonym"))  # NCBI's anamorph
        node, how = tree.place(ncbi, 27291)  # S. paradoxus: its genus EUKARYOME's, the species NCBI's
        self.assertEqual((node, how, tree.parent[node], tree.source[node]),
                         ((9, "Saccharomyces_paradoxus"), "genus", (8, "Saccharomyces"), "ncbi"))
        self.assertEqual(tree.place(ncbi, 6253)[1], "genus")  # Ascaris, of the row a column off
        # Cryptococcus: none of its family, order, class or phylum EUKARYOME's; under Fungi, NCBI's names below
        node, how = tree.place(ncbi, 40410)
        self.assertEqual((how, tree.text(node)), (
            "upward", "Obazoa;Fungi;Fungi.subreg.incertae.sedis;Basidiomycota;Basidiomycota.subphy.incertae.sedis;"
                      "Tremellomycetes;Tremellales;Cryptococcaceae;Cryptococcus;Cryptococcus_deneoformans"))
        node, how = tree.place(ncbi, 952)  # Mysterium: none of its clades
        self.assertEqual((how, tree.text(node).split(";")[0], node), ("none", eukaryome.TOP,
                                                                     (9, "Mysterium_obscurum")))
        # the plant Morus: the plant's genus, never the bird's
        self.assertIn("Viridiplantae", tree.text(tree.place(ncbi, 3498)[0]))
        # a taxon above species by its own name, to its rank
        self.assertEqual(tree.place(ncbi, 4893), ((7, "Saccharomycetaceae"), "name"))
        self.assertEqual(tree.place(ncbi, 5206)[0], (8, "Cryptococcus"))
        # a record of a genus whose label names a species EUKARYOME has not: its label's, deeper
        label = ("k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;"
                 "g__Saccharomyces;s__novum")
        self.assertEqual(eukaryome.record_node(tree, ncbi, 4930, label), ((9, "Saccharomyces_novum"), "label"))
        self.assertEqual(eukaryome.record_node(tree, ncbi, 4932, label)[1], "species")  # a species: its own
        rows = tree.rows([tree.place(ncbi, 40410)[0]])
        self.assertEqual(rows[0]["name"], "Eukaryota")
        added = {r["name"]: r for r in rows}
        self.assertEqual((added["Cryptococcus"]["source"], added["Cryptococcus"]["ncbi_taxid"],
                          added["Cryptococcus"]["parent"]), ("ncbi", 5206, "Cryptococcaceae"))
        self.assertEqual(tree.stamp()["placed"]["genus"] >= 1, True)


HARMONY_ROWS = [
    HEADER,
    ("Obazoa", "Metazoa", "", "Cnidaria", "", "Anthozoa", "Alcyonacea", "Alcyoniidae", "Alcyonium", ""),
    ("Obazoa", "Metazoa", "", "Cnidaria", "", "Anthozoa", "Alcyonacea", "Nephtheidae", "Nephthea", ""),
    ("Obazoa", "Metazoa", "", "Cnidaria", "", "Anthozoa", "Alcyonacea", "Alcyoniidae", "Zebra", ""),
    # EUKARYOME 2.1's golden moles: Mammalia, but the corals' order
    ("Obazoa", "Metazoa", "Bilateria", "Chordata", "", "Mammalia", "Alcyonacea", "Chrysochloridae", "Amblysomus", ""),
    ("Obazoa", "Metazoa", "Bilateria", "Chordata", "", "Mammalia", "Primates", "Hominidae", "Homo", ""),
    ("Obazoa", "Metazoa", "Bilateria", "Chordata", "", "Actinopteri", "Characiformes", "Chilodontidae", "Caenotropus",
     ""),
    ("Obazoa", "Metazoa", "Bilateria", "Mollusca", "", "Gastropoda", "Seguenziida", "Chilodontidae.Seguenziida",
     "Turcica", ""),
    ("Obazoa", "Metazoa", "Bilateria", "Arthropoda", "Crustacea", "Hexanauplia", "Calanoida", "Calanidae", "Calanus",
     ""),
    ("Obazoa", "Metazoa", "Bilateria", "Arthropoda", "Crustacea", "Ichthyostraca", "Siphonostomatoida", "Caligidae",
     "Caligus", ""),
    ("Obazoa", "Metazoa", "Bilateria", "Arthropoda", "Crustacea", "Hexanauplia", "Porocephalida", "Linguatulidae",
     "Linguatula", ""),
    ("Archaeplastida", "Rhodoplantae", "Rhodophyta", "Cyanidiophyta", "Cyanidiophytina", "Cyanidiophyceae",
     "Cyanidiales", "Cyanidiaceae", "Galdieria", ""),
    ("Archaeplastida", "Rhodoplantae", "Rhodophyta", "Rhodophyta", "Rhodophytina", "Bangiophyceae", "Bangiales",
     "Bangiaceae", "Bangia", ""),
    ("Obazoa", "Fungi", "Rozellomyceta", "Rozellomycota", "Microsporidia", "Microsporidea", "Amblyosporida",
     "Amblyosporidae", "Amblyospora", ""),
]
HARMONY_NODES = [
    (6073, 33208, "phylum", "Cnidaria"), (6101, 6073, "class", "Anthozoa"), (40674, 7711, "class", "Mammalia"),
    (9347, 40674, "order", "Afrosoricida"), (9389, 9347, "family", "Chrysochloridae"),
    (9390, 9389, "genus", "Amblysomus"), (9391, 9390, "species", "Amblysomus hottentotus"),
    (9400, 9389, "genus", "Zebra"), (9401, 9400, "species", "Zebra striata"),
    (6447, 33208, "phylum", "Mollusca"), (6448, 6447, "class", "Gastropoda"), (6449, 6448, "family", "Chilodontidae"),
    (6450, 6449, "genus", "Clypeostoma"), (6451, 6450, "species", "Clypeostoma nortoni"),
    (6656, 33208, "phylum", "Arthropoda"), (6657, 6656, "class", "Copepoda"), (6658, 6657, "order", "Calanoida"),
    (6659, 6658, "family", "Calanidae"), (6660, 6659, "genus", "Calanus"),
    (6661, 6660, "species", "Calanus finmarchicus"), (6662, 6657, "order", "Mormonilloida"),
    (6663, 6662, "family", "Mormonillidae"), (6664, 6663, "genus", "Mormonilla"),
    (6665, 6664, "species", "Mormonilla phasma"), (6666, 6657, "order", "Siphonostomatoida"),
    (6667, 6666, "family", "Caligidae"), (6668, 6667, "genus", "Caligus"), (6669, 6668, "species", "Caligus elongatus"),
    (6670, 6656, "class", "Ichthyostraca"), (6671, 6670, "order", "Cephalobaenida"),
    (6672, 6671, "family", "Linguatulida"), (6673, 6672, "genus", "Linguatula"),
    (6674, 6673, "species", "Linguatula serrata"), (2763, 2759, "phylum", "Rhodophyta"),
    (2764, 2763, "class", "Galdieriophyceae"), (2765, 2764, "order", "Galdieriales"),
    (2766, 2765, "family", "Galdieriaceae"), (2767, 2766, "genus", "Galdieria"),
    (2768, 2767, "species", "Galdieria sulphuraria"), (6029, 4751, "phylum", "Microsporidia"),
    (6030, 6029, "genus", "Amblyospora"), (6031, 6030, "species", "Amblyospora connecticus"),
    (6032, 6029, "genus", "Nosema"), (6033, 6032, "species", "Nosema ceranae"),
]


class HarmonyTest(unittest.TestCase):
    """One taxonomic framework (since 2026-11-15): a species' higher taxa those of its relatives, and of its NCBI
    lineage where EUKARYOME holds a name in another lineage alone."""

    def setUp(self):
        tmp = self.enterContext(tempfile.TemporaryDirectory())
        with mock.patch.object(fx, "NODES", fx.NODES + MORE_NODES + HARMONY_NODES):
            self.ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))

    def lineage(self, tree, taxid):
        node, how = tree.place(self.ncbi, taxid)
        return how, tree.text(node)

    def test_a_name_misfiled_under_another_phylum_is_its_own_taxon(self):
        tree = eukaryome.Tree.from_classification(HARMONY_ROWS)
        self.assertIn((6, "Alcyonacea(Mammalia)"), tree.parent)  # the corals' order keeps its name
        how, text = self.lineage(tree, 9391)
        self.assertEqual(how, "genus")
        self.assertIn(";Chordata;", text)
        self.assertNotIn("Anthozoa", text)

    def test_a_genus_name_of_another_phylum_is_not_the_species(self):
        # EUKARYOME holds Zebra as a coral's genus alone: the mammal goes under its NCBI family, Chrysochloridae, its
        # genus written with its parent beside the coral's (the tree's homonyms)
        how, text = self.lineage(eukaryome.Tree.from_classification(HARMONY_ROWS), 9401)
        self.assertEqual(how, "upward")
        self.assertTrue(text.endswith(";Mammalia;Alcyonacea(Mammalia);Chrysochloridae;Zebra(Chrysochloridae);"
                                      "Zebra_striata"), text)

    def test_a_genus_classified_otherwise_within_ncbis_taxon_is_eukaryomes(self):
        # EUKARYOME's classification is Eukore's where NCBI's taxon holds it: Galdieria in Cyanidiophyta, under its
        # subkingdom Rhodophyta (NCBI's phylum); Linguatula in Hexanauplia, a class of the subphylum of NCBI's
        # Ichthyostraca, though no name below their class is shared
        tree = eukaryome.Tree.from_classification(HARMONY_ROWS)
        how, text = self.lineage(tree, 2768)
        self.assertEqual(how, "genus")
        self.assertIn(";Cyanidiophyta;", text)
        how, text = self.lineage(tree, 6674)
        self.assertEqual(how, "genus")
        self.assertIn(";Hexanauplia;Porocephalida;", text)

    def test_a_label_of_a_homonym_takes_its_own_lineage(self):
        # EUKARYOME files Triadinium as a ciliate's and a dinoflagellate's genus: two taxa, and its reads' label of the
        # dinoflagellate (g__Triadinium(Dinoflagellata)) goes to the dinoflagellates' by the label's higher names
        tree = eukaryome.Tree.from_classification([
            HEADER,
            ("Harosa", "Alveolata", "", "Ciliophora", "Intramacronucleata", "Litostomatea", "Trichostomatia",
             "Ditoxidae", "Triadinium", ""),
            ("Harosa", "Alveolata", "", "Dinoflagellata", "Dinophytina", "Dinophyceae", "Gonyaulacales",
             "Goniodomataceae", "Triadinium", "")])
        node = tree.place_label("k__Alveolata;p__Dinoflagellata;c__Dinophyceae;o__Gonyaulacales;f__Goniodomataceae;"
                                "g__Triadinium(Dinoflagellata);s__polyedricum")
        self.assertIn(";Dinoflagellata;", tree.text(node))
        self.assertTrue(tree.text(node).endswith(";Triadinium_polyedricum"), tree.text(node))
        node = tree.place_label("k__Alveolata;p__Ciliophora;c__Litostomatea;o__Trichostomatia;f__Ditoxidae;"
                                "g__Triadinium;s__caudatum")
        self.assertIn(";Ciliophora;", tree.text(node))

    def test_a_dotted_homonym_found_by_its_name(self):
        # NCBI's gastropod Chilodontidae: EUKARYOME's Chilodontidae.Seguenziida, not the fish family
        how, text = self.lineage(eukaryome.Tree.from_classification(HARMONY_ROWS), 6451)
        self.assertEqual(how, "upward")
        self.assertIn(";Mollusca;", text)
        self.assertIn("Chilodontidae.Seguenziida;Clypeostoma", text)

    def test_a_clade_named_otherwise_is_its_relatives_taxon(self):
        # NCBI's Copepoda is EUKARYOME's Hexanauplia (its Calanus): a copepod EUKARYOME lacks goes there, not under a
        # second class named Copepoda
        plain = eukaryome.Tree.from_classification(HARMONY_ROWS)
        self.assertIn(";Copepoda;", self.lineage(plain, 6665)[1])
        tree = eukaryome.Tree.from_classification(HARMONY_ROWS)
        learned = tree.learn(self.ncbi, [6661, 6665])
        self.assertEqual(learned[6657], (5, "Hexanauplia"))
        how, text = self.lineage(tree, 6665)
        self.assertEqual(how, "upward")
        self.assertIn(";Hexanauplia;Mormonilloida;Mormonillidae;Mormonilla;Mormonilla_phasma", text)
        self.assertEqual(self.lineage(tree, 6661)[0], "genus")  # placed by its own name as before
        # since 2026-11-15 (the audit's T1): the same where the copepod was placed before learning (the audit places
        # the published reads before the package learns); until then its NCBI-named class voted, nothing was learned
        early = eukaryome.Tree.from_classification(HARMONY_ROWS)
        self.assertIn(";Copepoda;", self.lineage(early, 6665)[1])
        self.assertEqual(early.learn(self.ncbi, [6661, 6665])[6657], (5, "Hexanauplia"))
        self.assertEqual(self.lineage(early, 6665), (how, text))
        self.assertEqual(early.learn(self.ncbi, [6661, 6665])[6657], (5, "Hexanauplia"))  # learned twice: the same

    def test_a_learned_clade_yields_to_a_deeper_taxon_of_its_name(self):
        # NCBI's phylum Microsporidia is EUKARYOME's phylum Rozellomycota by its species, and EUKARYOME's subphylum
        # Microsporidia below it: a microsporidian EUKARYOME lacks goes there, not under an empty subphylum
        tree = eukaryome.Tree.from_classification(HARMONY_ROWS)
        learned = tree.learn(self.ncbi, [6031, 6033])
        self.assertEqual(learned[6029], (3, "Rozellomycota"))
        how, text = self.lineage(tree, 6033)
        self.assertEqual(how, "upward")
        self.assertIn(";Rozellomycota;Microsporidia;", text)

    def test_a_clade_split_across_taxa_is_not_learned(self):
        # with Caligus EUKARYOME puts NCBI's copepods into two classes: a copepod it lacks goes neither under the class
        # most of them stand in nor under their subphylum with its ranks empty, but by its names as before
        tree = eukaryome.Tree.from_classification(HARMONY_ROWS)
        learned = tree.learn(self.ncbi, [6661, 6665, 6669])
        self.assertNotIn(6657, learned)
        self.assertEqual(learned[6668], (8, "Caligus"))
        how, text = self.lineage(tree, 6665)
        self.assertEqual(how, "upward")
        self.assertIn(";Arthropoda;Arthropoda.subphy.incertae.sedis;Copepoda;Mormonilloida;", text)


class LoadTest(unittest.TestCase):
    def test_the_source_folders_tree_holds_its_reads_species(self):
        # the build's tree: the classification table's xlsx (a formula's saved value read, not its text) and the long
        # reads' labels, so that an NCBI species is found by EUKARYOME's species (until 2026-11-14's fix: by its
        # genus alone, the reads' species never placed)
        try:
            import openpyxl  # noqa: PLC0415
        except ImportError:
            self.skipTest("no openpyxl")
        with tempfile.TemporaryDirectory() as tmp:
            book = openpyxl.Workbook()
            for row in CLASSIFICATION:
                book.active.append(list(row))
            book.save(os.path.join(tmp, "EUKARYOME2_1_classification-2.xlsx"))
            with gzip.open(os.path.join(tmp, "General_EUK_longread_v2.1.fasta.gz"), "wt") as fh:
                fh.write(">EUK1;k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;"
                         "f__Saccharomycetaceae;g__Saccharomyces;s__cerevisiae\nACGT\n")
            tree = eukaryome.load(tmp, log=lambda *_: None)
            ncbi = ncbi_world(tmp)
        self.assertEqual(tree.place(ncbi, 4932), ((9, "Saccharomyces_cerevisiae"), "species"))
        self.assertEqual(tree.source[(9, "Saccharomyces_cerevisiae")], "eukaryome_label")
        self.assertEqual(tree.md5 is not None, True)


class LabelTest(unittest.TestCase):
    def test_labels_above_the_genus(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_world(tmp)
        index = rrna.name_index(ncbi)
        base = "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;"
        self.assertEqual(rrna.label_taxon(base + "f__Saccharomycetaceae;g__Saccharomyces;s__cerevisiae", index, ncbi),
                         (4932, "eukaryome_name", "species"))
        self.assertEqual(rrna.label_taxon(base + "f__unclassified;g__unclassified;s__unclassified", index, ncbi),
                         (4892, "eukaryome_upward", "order"))
        # a genus NCBI lacks: its family's, the label's rank the genus
        self.assertEqual(rrna.label_taxon(base + "f__Saccharomycetaceae;g__Nogenus;s__unclassified", index, ncbi),
                         (4893, "eukaryome_upward", "genus"))
        # EUKARYOME's own lineages NCBI names none of: Eukaryota
        self.assertEqual(rrna.label_taxon("k__Eukaryota.reg112b;p__GS90;c__unclassified", index, ncbi),
                         (2759, "eukaryome_upward", "phylum"))
        # an organelle's set, and a label naming nothing: none
        self.assertEqual(rrna.label_taxon("k___mitochondrion;p__unclassified;g__Drosophila(Metazoa);s__unclassified",
                                          index, ncbi)[:2], (None, None))
        self.assertEqual(rrna.label_taxon("k__unclassified;p__unclassified", index, ncbi), (None, None, None))
        # a genus a bacterium holds is not taken: the label's kingdom
        self.assertEqual(rrna.label_taxon("k__Metazoa;p__Arthropoda;g__Paracoccus;s__unclassified", index, ncbi),
                         (33208, "eukaryome_upward", "genus"))
        self.assertEqual(rrna.label_taxon(base + "f__unclassified;g__unclassified;s__unclassified", index, ncbi,
                                          "ksgp")[1], "ksgp_upward")

    def test_ksgp_labels(self):
        self.assertEqual(rrna.ksgp_label("k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;"
                                         "f__Cryptococcaceae;g__Cryptococcus(Fungi);s__Cryptococcus(Fungi)_neoformans"),
                         "k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;"
                         "g__Cryptococcus(Fungi);s__neoformans")
        self.assertTrue(rrna.ksgp_label("k__Fungi;p__Ascomycota;c__Pichiomycetes;o__Serinales;f__Debaryomycetaceae;"
                                        "g__(Candida)(Debaryomycetaceae);s__(Candida)(Debaryomycetaceae)_auris")
                        .endswith("g__Candida(Debaryomycetaceae);s__auris"))
        self.assertEqual(rrna.ksgp_label("k__Eukaryota:plas;p__Chlorophyta:plas;c__Chlorophyceae:plas;"
                                         "o__Oedogoniales:plas;f__Oedogoniales_X:plas;g__Oedogonium:plas;"
                                         "s__Oedogonium_cardiacum:plas"),
                         "k__Eukaryota;p__Chlorophyta;c__Chlorophyceae;o__Oedogoniales;f__unclassified;"
                         "g__Oedogonium;s__cardiacum")
        self.assertEqual(rrna.ksgp_label("k__Metazoa;p__Nematoda;c__?;o__?;f__?;g__Ascaris;s__Ascaris_?"),
                         "k__Metazoa;p__Nematoda;c__unclassified;o__unclassified;f__unclassified;g__Ascaris;"
                         "s__unclassified")
        self.assertEqual(rrna.label_rank(rrna.ksgp_label("k__Metazoa;p__Nematoda;c__?;o__?;f__?;g__Ascaris;"
                                                         "s__Ascaris_?")), "genus")

    def test_ncbi_lineage_to_its_rank(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_world(tmp)
        self.assertEqual(labels.ncbi_lineage_to(ncbi, 4930), "Eukaryota;Fungi;Ascomycota;Saccharomycetes;"
                                                             "Saccharomycetales;Saccharomycetaceae;Saccharomyces;")
        self.assertEqual(labels.ncbi_lineage_to(ncbi, 559292).split(";")[-1], "Saccharomyces cerevisiae")
        self.assertEqual(labels.ncbi_lineage_to(ncbi, 2759), "Eukaryota;;;;;;;")
        self.assertEqual(labels.rank_of(ncbi, 952), "species")


class CheckTest(unittest.TestCase):
    def test_telomere_repeats(self):
        flank = "ACGTTGCA" * 20
        self.assertEqual(rrna.telomere(flank + "TTAGGG" * 6 + flank), 36)
        self.assertEqual(rrna.telomere(flank + "CCCTAA" * 7 + flank), 42)  # the other strand
        self.assertEqual(rrna.telomere(flank + "AGGGTT" * 6 + flank), 36)  # another phase
        self.assertEqual(rrna.telomere(flank + "TTAGGG" * 5 + flank), 0)
        self.assertEqual(rrna.telomere(flank + "TTTAGGG" * 6 + flank), 42)  # a plant's

    def test_tiers_of_format_7(self):
        row = dict(kind="18S", own_species="t", ambiguous="0", homopolymer=3, gaps=0)
        self.assertEqual(rrna.confidence(dict(row, label_from="eukaryome_upward")), ("medium", "higher_label"))
        self.assertEqual(rrna.confidence(dict(row, label_from="eukaryome_genus")), ("medium", "genus_label"))
        self.assertEqual(rrna.confidence(dict(row, telomere=42)), ("low", "telomere"))
        self.assertEqual(rrna.confidence(dict(row, flags_extra=["chimera"])), ("low", "chimera"))
        self.assertEqual(rrna.confidence(dict(row, flags_extra=["chimera_cut"])), ("medium", "chimera_cut"))
        self.assertEqual(rrna.confidence(dict(row, flags_extra=["flank", "vector"])), ("low", "vector,flank"))
        self.assertEqual(rrna.confidence(dict(row, flags_extra=["chimera_weak"])), ("medium", "chimera_weak"))
        self.assertEqual(rrna.confidence(dict(row, label_from="eukaryome_name")), ("high", "-"))


def seq_of(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng, s, share):
    return "".join(c if rng.random() >= share else rng.choice("ACGT".replace(c, "")) for c in s)


class AuditTest(unittest.TestCase):
    def test_vectors_but_low_complexity(self):
        rng = random.Random(3)
        vector = seq_of(rng, 120)
        vectors = audit.Vectors([("NGB1", vector), ("oligo_dT", "AAGCAGTG" + "T" * 40)])
        reads = [seq_of(rng, 300) + vector[:60] + seq_of(rng, 300),
                 seq_of(rng, 200) + "A" * 40 + seq_of(rng, 200),  # a poly-A run: no oligo-dT match
                 seq_of(rng, 500)]
        hits = vectors.hits(reads)
        self.assertEqual(hits[0], [(300, 360, "NGB1")])
        self.assertNotIn(1, hits)
        self.assertNotIn(2, hits)

    def test_end_cuts(self):
        rng = random.Random(4)
        rrna_part, vec = seq_of(rng, 1500), seq_of(rng, 80)
        # a vector before the rRNA: cut through it; inside the rRNA, a long match flagged, a short one not
        a, b, why, inside = audit.end_cuts("C" + vec + rrna_part, (81, 1581), [(1, 81, "pUC")], {})
        self.assertEqual((a, b, why, inside), (81, 0, ["5':vector pUC"], None))
        _a, _b, _w, inside = audit.end_cuts(rrna_part, (0, 1500), [(700, 760, "pUC")], {})
        self.assertEqual(inside, "pUC")
        _a, _b, _w, inside = audit.end_cuts(rrna_part, (0, 1500), [(700, 725, "pUC")], {})
        self.assertIsNone(inside)
        # a barcode before a primer at the start, no rRNA: cut, the primer kept; a primer at the end: kept
        primers = audit.rdna_primers()
        nsl = next(p for n, p, e in primers if n == "NS1short")
        read = "ACGTACGTAC" + nsl + rrna_part + "ACGT"
        ends = audit.primer_ends([read], primers, {0: (10, 4)})
        self.assertEqual(ends[0][5][:2], ("NS1short", 10))
        a, b, why, _i = audit.end_cuts(read, (10, len(read) - 4), [], ends[0])
        self.assertEqual((a, b, why), (10, 0, ["5':leader before primer NS1short"]))
        # a telomere beside the rRNA, Ns at an end, a poly-A tail
        read = "TTAGGG" * 8 + "ACGTTA" * 5 + rrna_part + "A" * 20 + "NNNN"
        a, b, why, _i = audit.end_cuts(read, (78, 78 + 1500), [], {})
        self.assertEqual((a, b, why), (48, 24, ["5':telomere", "3':Ns", "3':homopolymer"]))

    def test_flank_cut(self):
        lengths = {"18S_rRNA": 1800, "5_8S_rRNA": 150}
        # an 18S from its model's first position, 200 bases before it: beyond the gene, cut (10 kept)
        hits = [("18S_rRNA", 1, 1800, 200, 2000, "+")]
        self.assertEqual(audit.flank_cut(2100, (200, 2000), hits, lengths, 5), (190, "flank_beyond"))
        self.assertEqual(audit.flank_cut(2100, (200, 2000), hits, lengths, 3), (90, "flank_beyond"))
        # nhmmer's 18S reaching 30 bases further than barrnap's: the cut counted from nhmmer's
        hits = [("18S_rRNA", 1, 1800, 170, 2000, "+")]
        self.assertEqual(audit.flank_cut(2100, (200, 2000), hits, lengths, 5), (160, "flank_beyond"))
        hits = [("18S_rRNA", 1, 1800, 30, 2000, "+")]  # within FLANK of the read's start: nothing beside it
        self.assertEqual(audit.flank_cut(2100, (200, 2000), hits, lengths, 5), (0, None))
        # the 18S from its model's 300th position: the flank in place of the gene, kept and flagged
        hits = [("18S_rRNA", 300, 1800, 200, 1700, "+")]
        self.assertEqual(audit.flank_cut(1800, (200, 1700), hits, lengths, 5), (0, "flank"))
        # a 5.8S after an 18S record's 18S, a reversed piece: another gene or two parts joined, cut
        hits = [("18S_rRNA", 1, 1790, 0, 1790, "+"), ("5_8S_rRNA", 1, 150, 2000, 2150, "+")]
        self.assertEqual(audit.flank_cut(2300, (0, 1790), hits, lengths, 3), (500, "flank_gene"))
        hits = [("18S_rRNA", 1, 1790, 0, 1790, "+"), ("18S_rRNA", 1500, 1800, 1900, 2200, "-")]
        self.assertEqual(audit.flank_cut(2300, (0, 1790), hits, lengths, 3), (500, "flank_joined"))
        self.assertEqual(audit.flank_cut(1830, (0, 1790), hits[:1], lengths, 3), (0, None))

    def test_a_read_audited_and_applied(self):
        rng = random.Random(5)
        vectors = audit.Vectors([("NGB1", "GACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTAC")])
        read = "GACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC" + seq_of(rng, 1500)
        got = audit.audit_reads([("EUK_X", read, [("18S", 36, 1536)])], vectors)
        self.assertEqual(got["EUK_X"], (36, 1536, ["5':vector NGB1"], []))
        self.assertEqual(audit.parse_audit(audit.audit_text(got["EUK_X"])), (36, 1536, ["5':vector NGB1"], []))
        self.assertIsNone(audit.parse_audit("0|1|2||"))
        seq, row, trim = audit.apply("EUK_X", read, dict(source="EUKARYOME", regions_list=[("18S", 36, 1536)]),
                                     got["EUK_X"])
        self.assertEqual((len(seq), row["regions_list"], row["trimmed"], trim["start"]),
                         (1500, [("18S", 0, 1500)], "37-1536 of 1536: 5':vector NGB1", 37))
        self.assertEqual(audit.audit_reads([("EUK_Y", read, None)])["EUK_Y"][3], ["no_rrna"])
        # a cut that would leave too little: none, flagged
        got = audit.audit_reads([("EUK_Z", read[:300], [("18S", 36, 300)])], vectors)
        self.assertEqual(got["EUK_Z"], (0, 300, [], ["trim_short"]))

    @unittest.skipUnless(shutil.which("nhmmer") and audit.barrnap_hmm(shutil.which("barrnap")),
                         "no nhmmer or barrnap's models")
    def test_flanks_by_barrnaps_models(self):
        # C. elegans' rDNA unit: its 18S with 300 bases of its spacer before it and 150 of its ITS1 after it, both
        # beside a whole 18S: cut to the 18S and 10 bases
        unit = next(s for n, s in fx_units() if n.startswith("X03680"))
        hmm, nhmmer = audit.barrnap_hmm(shutil.which("barrnap")), shutil.which("nhmmer")
        hits = audit.nhmmer_hits([("unit", unit)], nhmmer, hmm)["unit"]
        s18 = max((h for h in hits if h[0] == "18S_rRNA" and h[5] == "+"), key=lambda h: h[4] - h[3])
        read = unit[s18[3] - 300:s18[4] + 150]
        got = audit.audit_reads([("r", read, [("18S", 300, 300 + s18[4] - s18[3])])], None, None, nhmmer, hmm)["r"]
        self.assertEqual([w.split(":")[1] for w in got[2]], ["flank_beyond", "flank_beyond"])
        self.assertLessEqual(got[0], 300 - audit.KEEP)
        self.assertGreaterEqual(got[0], 300 - audit.KEEP - audit.MODEL_END)


def fx_units():
    with gzip.open(os.path.join(os.path.dirname(HERE), "data", "model_rdna_units.fna.gz"), "rt") as fh:
        name, parts = None, []
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(parts).upper()
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
        if name:
            yield name, "".join(parts).upper()


class ChimeraTest(unittest.TestCase):
    @unittest.skipUnless(shutil.which("minimap2"), "no minimap2")
    def test_window_screen(self):
        rng = random.Random(7)
        a, b = seq_of(rng, 1800), seq_of(rng, 1800)
        refs = {f"g{i}": (mutate(rng, a if i < 3 else b, 0.004), "Ascomycota" if i < 3 else "Chordata", f"G{i}")
                for i in range(1, 5)}
        queries = {"chimera": (mutate(rng, a[:1000] + b[1000:], 0.002), "Ascomycota"),
                   "mislabel": (mutate(rng, b, 0.002), "Ascomycota"),
                   "clean": (mutate(rng, a, 0.002), "Ascomycota"),
                   "placeholder": (mutate(rng, b, 0.002), None)}
        with tempfile.TemporaryDirectory() as tmp:
            calls = chimeras.window_screen(queries, refs, shutil.which("minimap2"), 2, tmp)
        self.assertEqual({q: c[0] for q, c in calls.items()},
                         {"chimera": "chimera", "mislabel": "mislabel", "placeholder": "unplaced"})
        keep = calls["chimera"][1]  # its own part: the Ascomycota stretch
        self.assertTrue(keep[0] == 0 and 900 <= keep[1] <= 1050, keep)

    @unittest.skipUnless(shutil.which("vsearch"), "no vsearch")
    def test_uchime(self):
        rng = random.Random(8)
        a = seq_of(rng, 1500)
        b = mutate(rng, a, 0.15)
        others = {"A": a, "B": b, "C": mutate(rng, a, 0.2), "D": mutate(rng, b, 0.05)}
        screened = {"EUK_chim": mutate(rng, a[:700], 0.003) + mutate(rng, b[700:], 0.003),
                    "EUK_fine": mutate(rng, a, 0.004)}
        genera = {"A": ("Alpha", "Fam"), "B": ("Beta", "Fam"), "C": ("Alpha", "Fam"), "D": ("Beta", "Fam"),
                  "EUK_chim": ("Alpha", "Fam"), "EUK_fine": ("Alpha", "Fam")}
        with tempfile.TemporaryDirectory() as tmp:
            searcher = screen.Searcher(shutil.which("vsearch"), 2, tmp)
            calls = chimeras.uchime_screen(screened, others, genera, searcher, shutil.which("vsearch"), 2, tmp)
        self.assertEqual(set(calls), {"EUK_chim"})
        call, keep, evidence = calls["EUK_chim"]
        self.assertEqual(call, "chimera")
        self.assertTrue(keep[0] == 0 and 600 <= keep[1] <= 800, (keep, evidence))  # its own genus' part

    @unittest.skipUnless(shutil.which("vsearch"), "no vsearch")
    def test_uchime_parents_among_the_screened(self):
        # since 2026-11-15 (the audit's R2): a chimera whose parents are published records screened beside it (one
        # part: every record of the set screened at once) is found; until then a part's records were taken out of its
        # database and could be no parent of each other
        rng = random.Random(9)
        a = seq_of(rng, 1500)
        b = mutate(rng, a, 0.15)
        screened = {"EUK_a": a, "EUK_b": b, "EUK_chim": mutate(rng, a[:700], 0.003) + mutate(rng, b[700:], 0.003)}
        others = {"G": seq_of(rng, 1500)}
        genera = {"EUK_a": ("Alpha", "Fam"), "EUK_b": ("Beta", "Fam"), "EUK_chim": ("Alpha", "Fam"), "G": None}
        with tempfile.TemporaryDirectory() as tmp:
            searcher = screen.Searcher(shutil.which("vsearch"), 2, tmp)
            calls = chimeras.uchime_screen(screened, others, genera, searcher, shutil.which("vsearch"), 2, tmp)
        self.assertIn("EUK_chim", calls)
        self.assertNotIn("EUK_a", calls)

    def test_references_one_assembly_few_copies(self):
        # since 2026-11-15 (the audit's R1, R3): a genome's copies one genome whatever their names (its own records,
        # ROD's copy, EUKARYOME's EUK_GCA_ and _dup records), its distinct 18S, at most audit.REFS_PER_GENOME of them
        rows = {"g1": dict(source="genome", accession="GCA_000340665.1"),
                "g2": dict(source="genome", accession="GCA_000340665.1"),
                "r1": dict(source="ROD", accession="GCF_000340665.2"),
                "EUK_GCA_000340665_dup1": dict(source="EUKARYOME"),
                "EUK_GCA_000340665": dict(source="EUKARYOME")}
        self.assertEqual({audit.genome_id(n, r) for n, r in rows.items()}, {"GC_000340665"})
        self.assertEqual(screen.genome_of("x", {"accession": "GCA_000340665_dup1"}), "GC_000340665")
        seqs = {"g1": "A" * 1800, "g2": "A" * 1800, "r1": "C" * 1790, "EUK_GCA_000340665_dup1": "G" * 1700,
                "EUK_GCA_000340665": "T" * 1600}
        records = [(n, seqs[n], r) for n, r in rows.items()]
        lin = {n: [""] * (audit.CLADE + 1) for n in rows}
        for n in rows:
            lin[n][audit.CLADE] = "Fungi"
        refs = audit.references(records, seqs, lin)
        self.assertEqual(sorted(refs), ["EUK_GCA_000340665_dup1", "g1", "r1"])  # g2 a copy of g1; the 3 longest
        self.assertEqual({g for _s, _c, g in refs.values()}, {"GC_000340665"})

    def test_crossovers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "alns.txt")
            with open(path, "w") as fh:
                fh.write("------\nQuery   (   12 nt) q1\nParentA (   12 nt) A\nParentB (   12 nt) B\n\n"
                         "A     1 ACGTAC-GTACG 11\nQ     1 ACGTACAGTACG 12\nB     1 TTGTACAGTTTT 12\n"
                         "Diffs   \nVotes   \nModel   AAAAAAxxBBBB\n")
            self.assertEqual(chimeras.crossovers(path), {"q1": (7, "A")})


class LinksTest(unittest.TestCase):
    def test_records_linked_to_pr2_and_ksgp(self):
        with tempfile.TemporaryDirectory() as tmp:
            pr2, ksgp_tax = os.path.join(tmp, "pr2.tax.gz"), os.path.join(tmp, "KSGP_v5.2_plus2.tax")
            with gzip.open(pr2, "wt") as fh:
                fh.write("AB000001.1.1500_U\tEukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Saccharomycetes;"
                         "Saccharomycetales;Saccharomyces;Saccharomyces_paradoxus;\n")
            with open(ksgp_tax, "w") as fh:
                fh.write("EUK1\tk__Fungi;p__Ascomycota;c__?;o__?;f__?;g__Saccharomyces;s__Saccharomyces_?\n"
                         "AB000001_dup1\tk__Fungi;p__Ascomycota;c__?;o__?;f__?;g__Saccharomyces;s__S_paradoxus\n"
                         "RS_GCF_1~NZ_1\tk__Bacteria;p__Bacillota;c__B;o__B;f__B;g__B;s__B b\n")
            indexes = {"PR2": links.read_index(pr2, "PR2"), "KSGP": links.read_index(ksgp_tax, "KSGP")}
        records = [("EUK_EUK1", "ACGT", dict(source="EUKARYOME")), ("EUK_AB000001", "ACGT", dict(source="EUKARYOME")),
                   ("GCA_1|c1", "ACGT", dict(source="genome"))]
        rows, counts = links.rows(records, indexes)
        self.assertEqual(sorted((r["record"], r["reference"], r["id"]) for r in rows),
                         [("EUK_AB000001", "KSGP", "AB000001_dup1"), ("EUK_AB000001", "PR2", "AB000001.1.1500_U"),
                          ("EUK_EUK1", "KSGP", "EUK1")])
        self.assertEqual(counts, {"KSGP": 2, "PR2": 1})
        self.assertEqual(links.accession_key("GB_GCA_7.1~LC000007.1.1009_U"), "LC000007")


class BuildTest(unittest.TestCase):
    def test_published_reads_of_eukaryome(self):
        # the build's published reads: EUKARYOME's long reads and SSU set (an 18S of a long read taken whole, the
        # compromised left out, labelled above the genus too); its mitochondrial SSU organelle genes
        import argparse  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import sources  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_world(tmp)
            src = os.path.join(tmp, "sources")
            for name in ("EUKARYOME", "EUKARYOME_SSU"):
                folder = sources.folder_of(src, name)
                os.makedirs(folder)
                with open(os.path.join(folder, sources.STATE), "w") as fh:
                    json.dump({"ready": True, "version": sources.SOURCES[name]["version"]}, fh)
            euk, ssu = sources.folder_of(src, "EUKARYOME"), sources.folder_of(src, "EUKARYOME_SSU")
            with open(os.path.join(euk, "EUK2_1_compromised_quality.tsv"), "w") as fh:
                fh.write("accession\nEUK6\n")
            order = "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;"
            yeast = order + "f__Saccharomycetaceae;g__Saccharomyces;s__cerevisiae"
            with gzip.open(os.path.join(euk, "General_EUK_longread_v2.1.fasta.gz"), "wt") as fh:
                fh.write(f">EUK1;{order}f__unclassified;g__unclassified;s__unclassified\n{'ACGT' * 10}\n"
                         f">EUK6;{yeast}\n{'ACGT' * 10}\n"
                         # a mitochondrial long read: no nucleus' read; its SSU set's SSU an organelle gene
                         ">EUK8;k___mitochondrion;p__unclassified;c__unclassified;o__unclassified;"
                         f"f__unclassified;g__Homo;s__sapiens\n{'ACGA' * 30}\n"
                         # a bacterium's SSU among the long reads: left out
                         ">EUK9;k___Bacteria;p__unclassified;c__unclassified;o__unclassified;"
                         f"f__unclassified;g__unclassified;s__unclassified\n{'ACGA' * 30}\n")
            with gzip.open(os.path.join(ssu, "General_EUK_SSU_v2.1.fasta.gz"), "wt") as fh:
                fh.write(f">EUK1;{yeast}\n{'GGGG' * 10}\n"            # a long read's 18S: the read taken
                         f">EUK6;{yeast}\n{'GGGG' * 10}\n"            # compromised
                         f">EUK7;{yeast}\n{'CCCC' * 10}\n"
                         f">AB000001;{order}f__Saccharomycetaceae;g__Saccharomyces;s__paradoxus\n{'TTTT' * 10}\n"
                         ">EUK8;k___mitochondrion;p__unclassified;c__unclassified;o__unclassified;"
                         f"f__unclassified;g__Homo;s__sapiens\n{'ACGA' * 10}\n"
                         ">EUK9;k___Bacteria;p__unclassified;c__unclassified;o__unclassified;"
                         f"f__unclassified;g__unclassified;s__unclassified\n{'ACGA' * 10}\n"
                         ">EUK10;k___nucleomorph;p__unclassified;c__unclassified;o__unclassified;"
                         f"f__unclassified;g__Saccharomyces;s__cerevisiae\n{'AGGA' * 10}\n")
            opts = argparse.Namespace(inside=lambda t: True)
            with mock.patch.object(eukore_build, "say", lambda *_a, **_k: None):
                reads = {r[0]: r for r in eukore_build.published_reads(opts, src, ncbi)}
                organelle = eukore_build.organelle_published(opts, src, ncbi, {})
        self.assertEqual(set(reads), {"EUK_EUK1", "EUK_EUK7", "EUK_AB000001"})
        self.assertEqual((reads["EUK_EUK1"][1], reads["EUK_EUK1"][4], reads["EUK_EUK1"][6]["via"]),
                         ("ACGT" * 10, "eukaryome_upward", "EUKARYOME_long"))
        self.assertEqual(reads["EUK_AB000001"][2:6], (27291, "EUKARYOME", "eukaryome_name", "AB000001"))
        self.assertEqual(reads["EUK_EUK7"][6]["via"], "EUKARYOME_SSU")
        genes = {n: row for n, _s, row in organelle}
        # its organelle SSU organelle genes (a nucleomorph's 18S too); a bacterium's SSU nowhere
        self.assertEqual(set(genes), {"EUK_EUK8", "EUK_EUK10"})
        self.assertEqual((genes["EUK_EUK8"]["species_taxid"], genes["EUK_EUK8"]["name"],
                          genes["EUK_EUK8"]["source"]), (9606, "12S", "EUKARYOME"))
        self.assertEqual((genes["EUK_EUK10"]["organelle"], genes["EUK_EUK10"]["name"]), ("nucleomorph", "18S"))
        self.assertEqual([eukore_build.eukaryome_kind(x) for x in ("k__Fungi;p__X", "k___plastid;p__X",
                                                                   "k___Archaea;p__X", "k___pseudogene.Fungi;p__X")],
                         ["nuclear", "plastid", None, None])


class Midori2Test(unittest.TestCase):
    def test_midori2_ssu_one_record_a_species(self):
        # MIDORI2's mitochondrial SSU as organelle genes (since 2026-11-14): the screen's drops left out, its cuts
        # taken as cut, a record for review only where its species has no other (flagged), one record a species (its
        # longest), none of a species whose RefSeq mitochondrial genome gives its genes
        import argparse  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import organelle  # noqa: PLC0415

        def head(acc, taxid, name):
            return f">{acc}###root_1;Eukaryota_2759;X_1;Y_2;Z_3;W_4;V_5;{name}_{taxid}\n"
        tmp = self.enterContext(tempfile.TemporaryDirectory())
        ncbi = ncbi_world(tmp)
        path = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_srRNA_RAW.fasta.gz")
        with gzip.open(path, "wt") as fh:
            fh.write(head("A1.1", 4932, "Saccharomyces cerevisiae") + "ACGT" * 100 + "\n"
                     + head("A2.1", 4932, "Saccharomyces cerevisiae") + "ACGT" * 150 + "\n"   # longest, dropped
                     + head("A3.1", 4932, "Saccharomyces cerevisiae") + "ACGA" * 120 + "\n"   # cut to 101-400
                     + head("B1.1", 27291, "Saccharomyces paradoxus") + "ACGT" * 90 + "\n"    # review alone
                     + head("C1.1", 9606, "Homo sapiens") + "ACGT" * 90 + "\n"                # RefSeq's genome
                     + head("D1.1", 2, "Bacteria") + "ACGT" * 90 + "\n")                      # no eukaryote
        calls = {("SSU", "A2.1"): {"action": "drop"},
                 ("SSU", "A3.1"): {"action": "trim", "keep_start": "101", "keep_end": "400"},
                 ("SSU", "B1.1"): {"action": "review"}}
        opts = argparse.Namespace(midori2=[path], inside=lambda t: True)
        with mock.patch.object(eukore_build, "say", lambda *_a, **_k: None):
            got = eukore_build.midori2_published(opts, ncbi, {(9606, "mitochondrion"): []}, calls)
        rows = {n: (s, r) for n, s, r in got}
        self.assertEqual(set(rows), {"MIDORI2_SSU_A1.1", "MIDORI2_SSU_B1.1"})  # A1 (400) beats A3's part (300)
        self.assertEqual((rows["MIDORI2_SSU_A1.1"][1]["source"], rows["MIDORI2_SSU_A1.1"][1]["via"],
                          rows["MIDORI2_SSU_A1.1"][1]["label_from"], rows["MIDORI2_SSU_A1.1"][1]["organelle"]),
                         ("MIDORI2", "MIDORI2_GB273", "midori2_taxid", "mitochondrion"))
        self.assertEqual(rows["MIDORI2_SSU_B1.1"][1]["flags_extra"], ["screen_weak"])
        row = dict(rows["MIDORI2_SSU_B1.1"][1], own_species="t", length_check="ok")
        self.assertEqual(organelle.confidence(row), ("medium", "screen_weak"))
        calls[("SSU", "A1.1")] = {"action": "drop"}
        with mock.patch.object(eukore_build, "say", lambda *_a, **_k: None):
            got = {n: s for n, s, _r in eukore_build.midori2_published(opts, ncbi, {}, calls)}
        self.assertEqual(got["MIDORI2_SSU_A3.1"], ("ACGA" * 120)[100:400])  # the screen's cut, its kept part
        self.assertIn("MIDORI2_SSU_C1.1", got)  # no RefSeq genome given now

    def test_midori2_cox1_and_cytb(self):
        # MIDORI2's COX1 and CYTB (since 2026-11-15): a record's reading frame found (its first codon not marked), a
        # record with a stop in every frame of both strands (a NUMT, a pseudogene) not taken, nor one of the wrong
        # length; the coding sequence the frame's, translated by the organelle check without stops
        import argparse  # noqa: PLC0415
        import random  # noqa: PLC0415
        import eukore_build  # noqa: PLC0415
        from eukore import gencode, organelle  # noqa: PLC0415
        rng = random.Random(15)
        sense = [c for c in (a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT")
                 if c not in ("TAA", "TAG", "TGA")]

        def orf(n):
            return "".join(rng.choice(sense) for _ in range(n))

        def head(acc, taxid, name):
            return f">{acc}###root_1;Eukaryota_2759;X_1;Y_2;Z_3;W_4;V_5;{name}_{taxid}\n"
        tmp = self.enterContext(tempfile.TemporaryDirectory())
        ncbi = ncbi_world(tmp)
        cox1 = "A" + orf(200) + "GC"                                    # its frame one base in
        numt = orf(60) + "TAAATAAAATAA" + orf(60) + "TTATTTATTTTA" + orf(60)  # stops in every frame, both strands
        cytb = orf(130)
        co1 = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_CO1_RAW.fasta.gz")
        with gzip.open(co1, "wt") as fh:
            fh.write(head("E1.1", 4932, "Saccharomyces cerevisiae") + cox1 + "\n"
                     + head("E2.1", 27291, "Saccharomyces paradoxus") + numt + "\n"
                     + head("E3.1", 27291, "Saccharomyces paradoxus") + orf(50) + "\n")  # 150 bases: too short
        cyt = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_Cytb_RAW.fasta.gz")
        with gzip.open(cyt, "wt") as fh:
            fh.write(head("F1.1", 4932, "Saccharomyces cerevisiae") + cytb + "\n")
        opts = argparse.Namespace(midori2=[co1, cyt], inside=lambda t: True)
        with mock.patch.object(eukore_build, "say", lambda *_a, **_k: None):
            got = eukore_build.midori2_published(opts, ncbi, {}, None)
        rows = {n: (s, r) for n, s, r in got}
        self.assertEqual(set(rows), {"MIDORI2_COX1_E1.1", "MIDORI2_CYTB_F1.1"})
        seq, row = rows["MIDORI2_COX1_E1.1"]
        self.assertEqual((row["gene"], row["organelle"], seq, row["cds"]), ("COX1", "mitochondrion", cox1, cox1[1:601]))
        self.assertEqual(gencode.protein_of(row["cds"], 1)[1], 0)
        self.assertEqual(organelle.reading_frame(numt)[1] > 0, True)
        _out, counts = organelle.midori2_records([co1, cyt], ncbi)
        self.assertEqual((counts["COX1: no reading frame (a NUMT, a pseudogene)"],
                          counts["COX1: not the gene's length"]), (1, 1))
        checked = organelle.check(got, ncbi)
        self.assertEqual({r["record"]: r["internal_stops"] for _n, _s, r in checked},
                         {"MIDORI2_COX1_E1.1": 0, "MIDORI2_CYTB_F1.1": 0})

    def test_midori2_records_on_their_genes_strand(self):
        # KSGP 5.1 turned 721 of MIDORI's mitochondrial SSU stored reversed: a coding record by its frame, an rRNA by
        # the strand model of the records taken (most of them forward)
        import random  # noqa: PLC0415
        from eukore import organelle  # noqa: PLC0415
        rng = random.Random(16)
        sense = [c for c in (a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT")
                 if c not in ("TAA", "TAG", "TGA")]

        def head(acc, taxid, name):
            return f">{acc}###root_1;Eukaryota_2759;X_1;Y_2;Z_3;W_4;V_5;{name}_{taxid}\n"

        def mutated(s, share):
            return "".join(c if rng.random() > share else rng.choice("ACGT") for c in s)
        tmp = self.enterContext(tempfile.TemporaryDirectory())
        ncbi = ncbi_world(tmp)
        cox1 = "".join(rng.choice(sense) for _ in range(200))
        co1 = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_CO1_RAW.fasta.gz")
        with gzip.open(co1, "wt") as fh:
            fh.write(head("R1.1", 4932, "Saccharomyces cerevisiae") + rrna.revcomp(cox1) + "\n")
        ssu = "".join(rng.choice("ACGT") for _ in range(900))
        sr = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_srRNA_RAW.fasta.gz")
        with gzip.open(sr, "wt") as fh:
            for acc, taxid, name in (("S1.1", 4932, "Saccharomyces cerevisiae"), ("S2.1", 27291, "Saccharomyces "
                                     "paradoxus"), ("S3.1", 9900001, "Saccharomyces novus")):
                fh.write(head(acc, taxid, name) + mutated(ssu, 0.03) + "\n")
            fh.write(head("S4.1", 9606, "Homo sapiens") + rrna.revcomp(mutated(ssu, 0.03)) + "\n")
        got, counts = organelle.midori2_records([co1, sr], ncbi)
        rows = {n: s for n, s, _r in got}
        self.assertEqual(rows["MIDORI2_COX1_R1.1"], cox1)  # turned by its frame
        self.assertEqual(counts["SSU: turned to its gene's strand"], 1)
        self.assertEqual(counts["COX1: turned to its gene's strand"], 1)
        self.assertGreater(organelle.strand_score(rows["MIDORI2_SSU_S4.1"], organelle.strand_model([ssu])), 0)

    def test_midori2_lsu(self):
        # MIDORI2's mitochondrial LSU (since 2026-11-16: protal's fourth gene, which protal takes from the package
        # alone): one record a species, of the gene's length; a stretch MIDORI2 files under the srRNA and the lrRNA of
        # one species taken once, as its SSU
        from eukore import organelle  # noqa: PLC0415
        rng = random.Random(17)

        def head(acc, taxid, name):
            return f">{acc}###root_1;Eukaryota_2759;X_1;Y_2;Z_3;W_4;V_5;{name}_{taxid}\n"
        tmp = self.enterContext(tempfile.TemporaryDirectory())
        ncbi = ncbi_world(tmp)
        lsu = "".join(rng.choice("ACGT") for _ in range(1200))
        both = "".join(rng.choice("ACGT") for _ in range(1500))
        lr = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_lrRNA_RAW.fasta.gz")
        with gzip.open(lr, "wt") as fh:
            fh.write(head("L1.1", 4932, "Saccharomyces cerevisiae") + lsu + "\n"
                     + head("L2.1", 9606, "Homo sapiens") + lsu[:250] + "\n"             # too short
                     + head("T1.1", 27291, "Saccharomyces paradoxus") + both + "\n")    # filed under both genes
        sr = os.path.join(tmp, "MIDORI2_UNIQ_NUC_GB273_srRNA_RAW.fasta.gz")
        with gzip.open(sr, "wt") as fh:
            fh.write(head("T1.1", 27291, "Saccharomyces paradoxus") + both + "\n")
        self.assertIn("LSU", organelle.MIDORI2_GENES)
        got, counts = organelle.midori2_records([lr, sr], ncbi)
        rows = {n: (s, r) for n, s, r in got}
        self.assertEqual(set(rows), {"MIDORI2_LSU_L1.1", "MIDORI2_SSU_T1.1"})
        self.assertEqual(rows["MIDORI2_LSU_L1.1"][1]["gene"], "LSU")
        self.assertIn(rows["MIDORI2_LSU_L1.1"][0], (lsu, rrna.revcomp(lsu)))  # the strand model of one record
        self.assertEqual(counts["LSU: the species' record of another gene too (taken as that one)"], 1)
        self.assertEqual(counts["LSU: not the gene's length"], 1)


class PackageTest(unittest.TestCase):
    def test_a_package_of_format_7(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_world(tmp)
            tree = tree_of()
            label = ("k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetales.fam08;"
                     "g__unclassified;s__unclassified")
            base = dict(kind="18S", accession="-", regions="-", length=40, lsu_complete="f", confidence="high",
                        flags="-")
            records = [("EUK_EUK1", "ACGT" * 10, dict(base, record="EUK_EUK1", source="EUKARYOME", taxid=4932,
                                                       species_taxid=4932, label_from="eukaryome_name",
                                                       source_label="-")),
                       ("EUK_EUK2", "ACGT" * 10, dict(base, record="EUK_EUK2", source="EUKARYOME", taxid=4892,
                                                       species_taxid=4892, label_from="eukaryome_upward",
                                                       source_label=label, confidence="medium",
                                                       flags="higher_label")),
                       ("EUK_AB000001", "ACGT" * 10, dict(base, record="EUK_AB000001", source="EUKARYOME",
                                                           taxid=27291, species_taxid=27291,
                                                           label_from="eukaryome_name", via="EUKARYOME_SSU",
                                                           trimmed="11-50 of 50: 5':flank_beyond"))]
            organelle = [("NC_1_SSU", "ACGT" * 10, dict(record="NC_1_SSU", gene="SSU", organelle="mitochondrion",
                                                         taxid=9606, species_taxid=9606, source="refseq",
                                                         confidence="high", flags="-")),
                         ("EUK_X", "ACGT" * 10, dict(record="EUK_X", gene="SSU", organelle="mitochondrion",
                                                      taxid=6253, species_taxid=6253, source="EUKARYOME",
                                                      confidence="unverified", flags="-"))]
            collected = {"genomes": {}, "copy_rows": [], "taxids": set(), "work": tmp, "done": []}
            out = package.write(os.path.join(tmp, "Eukore_2026.11"), "2026.11", ncbi, collected, [], records,
                                organelle=organelle, eukaryome=tree)
            self.assertIsNone(package.check(out, deep=True))
            species = {r["taxid"]: r for r in package.read_table(out, "species.tsv")}
            # species alone: the order a read is labelled by is no row; a species of organelle genes alone is
            self.assertEqual({t: r["source"] for t, r in species.items()},
                             {"4932": "rrna", "27291": "rrna", "9606": "organelle", "6253": "organelle"})
            self.assertEqual((species["27291"]["eukaryome_match"], species["9606"]["eukaryome_match"]),
                             ("genus", "genus"))
            self.assertTrue(species["27291"]["eukaryome"].endswith(";Saccharomyces;Saccharomyces_paradoxus"))
            got = {r["record"]: r for r in package.read_table(out, "rrna/records.tsv")}
            self.assertEqual((got["EUK_EUK2"]["rank"], got["EUK_EUK2"]["eukaryome_match"]), ("family", "label"))
            self.assertTrue(got["EUK_EUK2"]["eukaryome_lineage"].endswith(";Saccharomycetales;"
                                                                          "Saccharomycetales.fam08;;"))
            self.assertEqual(got["EUK_EUK2"]["ncbi_lineage"], "Eukaryota;Fungi;Ascomycota;Saccharomycetes;"
                                                              "Saccharomycetales;;;")
            self.assertEqual((got["EUK_AB000001"]["rank"], got["EUK_AB000001"]["via"],
                              got["EUK_AB000001"]["trimmed"]), ("species", "EUKARYOME_SSU",
                                                                 "11-50 of 50: 5':flank_beyond"))
            genes = {r["record"]: r for r in package.read_table(out, "organelle/genes.tsv")}
            self.assertTrue(genes["EUK_X"]["eukaryome_lineage"].endswith(";Ascaris;Ascaris_lumbricoides"))
            # the EUKARYOME taxonomy, cut to the package's lineages: every name of every lineage a taxon of it
            nodes = {(r["rank"], r["name"]): r for r in package.read_table(out, eukaryome.TABLE)}
            for lineage in [r["eukaryome_lineage"] for r in got.values()] + [r["eukaryome"] for r in
                                                                             species.values()]:
                for rank, name in zip(eukaryome.RANKS, lineage.split(";")):
                    if name:
                        self.assertIn((rank, name), nodes)
            with open(os.path.join(out, "manifest.json")) as fh:
                manifest = json.load(fh)
            counts = manifest["counts"]
            self.assertEqual((manifest["format"], counts["species_rrna_only"], counts["species_organelle_only"],
                              counts["species_organelle_only_by_organelle"]), (7, 2, 2, {"mitochondrion": 2}))
            self.assertEqual(counts["rrna_sources"], {"EUKARYOME": {"high": 2, "medium": 1}})
            self.assertEqual(counts["rrna_ranks"]["family"], {"medium": 1})
            self.assertEqual(manifest["inputs"]["eukaryome_taxonomy"]["version"], "2.1")


if __name__ == "__main__":
    unittest.main()
