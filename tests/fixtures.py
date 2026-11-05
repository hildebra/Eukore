"""A small world for EUKDB's tests (copied from protal's scripts/mini_db/test_eukaryotes.py): six NCBI species of fungi,
one of an unknown phylum and a nematode, their EukDetect2 genomes and loci, EukDetect1's records, GenBank assemblies."""

import io
import os
import random
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukdb import sources as eukaryotes  # noqa: E402


# NCBI's taxonomy: (taxid, parent, rank, name)
NODES = [
    (1, 1, "no rank", "root"),
    (2759, 1, "superkingdom", "Eukaryota"),
    (4751, 2759, "kingdom", "Fungi"),
    (4890, 4751, "phylum", "Ascomycota"),
    (4891, 4890, "class", "Saccharomycetes"),
    (4892, 4891, "order", "Saccharomycetales"),
    (4893, 4892, "family", "Saccharomycetaceae"),
    (4930, 4893, "genus", "Saccharomyces"),
    (4932, 4930, "species", "Saccharomyces cerevisiae"),
    (559292, 4932, "strain", "Saccharomyces cerevisiae S288C"),
    (4910, 4893, "genus", "Kluyveromyces"),
    (28985, 4910, "species", "Kluyveromyces lactis"),
    (900, 4893, "genus", "Newgenus"),
    (901, 900, "species", "Newgenus novus"),
    (1535326, 4892, "family", "Debaryomycetaceae"),
    (5475, 1535326, "genus", "Candida"),
    (5476, 5475, "species", "Candida albicans"),
    (5204, 4751, "phylum", "Basidiomycota"),
    (155616, 5204, "class", "Tremellomycetes"),
    (5234, 155616, "order", "Tremellales"),
    (1884633, 5234, "family", "Cryptococcaceae"),
    (5206, 1884633, "genus", "Cryptococcus"),
    (40410, 5206, "species", "Cryptococcus deneoformans"),
    (960, 4892, "species", "[Candida] bracketa"),  # no true Candida: under the order, NCBI gives it no genus
    (950, 2759, "phylum", "Mysteriozoa"),
    (951, 950, "genus", "Mysterium"),
    (952, 951, "species", "Mysterium obscurum"),
    # EukDetect1's alone: a nematode, and a protist it has one mapped locus of
    (33208, 2759, "kingdom", "Metazoa"),
    (6231, 33208, "phylum", "Nematoda"),
    (6251, 6231, "genus", "Ascaris"),
    (6253, 6251, "species", "Ascaris lumbricoides"),
    (970, 950, "species", "Mysterium rarum"),
]


SYNONYMS = [(28985, "Candida sphaerica", "anamorph")]


KSGP = [
    "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae",
    "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Saccharomyces;s__Saccharomyces_cerevisiae",
    "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Debaryomycetaceae;g__Lodderomyces;s__Lodderomyces_albicans",
    "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Candida;s__Candida_sphaerica",
    "k__Fungi;p__Basidiomycota;c__Tremellomycetes;o__Tremellales;f__Cryptococcaceae;g__Cryptococcus(Fungi);s__Cryptococcus(Fungi)_neoformans",
    "k__Eukaryota:mito;p__Ascomycota;c__x;o__x;f__x;g__Mitoella;s__Mitoella_mito",  # an organelle: left out
]


# EukDetect2's genomes: name -> (taxid, group, loci)
GENOMES = {
    "Saccharomyces_cerevisiae_S288C": (559292, "Fungi", ["1at2759", "2at2759", "3at2759"]),
    "Saccharomyces_cerevisiae_YJM789": (4932, "Fungi", ["1at2759", "2at2759"]),
    "Kluyveromyces_lactis": (28985, "Fungi", ["1at2759", "2at2759", "3at2759"]),
    "Candida_albicans_SC5314": (5476, "Fungi", ["1at2759", "2at2759"]),
    "Newgenus_novus": (901, "Fungi", ["1at2759", "3at2759"]),
    "Cryptococcus_deneoformans": (40410, "Fungi", ["1at2759"]),
    "Mysterium_obscurum": (952, "Protist", ["2at2759", "4at2759"]),
}


ASSEMBLY_COLUMNS = ("assembly_accession", "taxid", "species_taxid", "organism_name", "infraspecific_name", "isolate",
                    "version_status", "assembly_level", "asm_name", "genome_size", "contig_count", "excluded_from_refseq",
                    "ftp_path")


ASSEMBLIES = [
    ("GCA_000146045.2", 559292, 4932, "Saccharomyces cerevisiae S288C", "strain=S288C", "", "latest", "Complete Genome", "R64", 12100000, 17, ""),
    ("GCA_000000003.1", 4932, 4932, "Saccharomyces cerevisiae", "strain=YJM789", "", "latest", "Scaffold", "YJM", 11900000, 80, ""),
    ("GCA_000000004.1", 4932, 4932, "Saccharomyces cerevisiae", "strain=X", "", "latest", "Contig", "X", 12000000, 300, ""),
    ("GCA_000000005.1", 4932, 4932, "Saccharomyces cerevisiae", "strain=old", "", "replaced", "Contig", "old", 12000000, 300, ""),
    # no genomes: a gene's record as a "complete genome" (GenBank has such of Aspergillus flavus), and a partial assembly
    ("GCA_000000560.1", 4932, 4932, "Saccharomyces cerevisiae S288C", "strain=S288C", "", "latest", "Complete Genome", "J", 560, 1, ""),
    ("GCA_000000006.1", 4932, 4932, "Saccharomyces cerevisiae", "strain=P", "", "latest", "Contig", "P", 2500000, 900, ""),
    ("GCA_000002515.1", 28985, 28985, "Kluyveromyces lactis", "strain=CBS 2359", "", "latest", "Chromosome", "ASM251v1", 10700000, 6, ""),
    ("GCA_000182965.3", 5476, 5476, "Candida albicans SC5314", "strain=SC5314", "", "latest", "Chromosome", "ASM18296v3", 14300000, 8, ""),
    ("GCA_000000901.1", 901, 901, "Newgenus novus", "", "", "latest", "Scaffold", "NN1", 30000000, 40, ""),
    ("GCA_000000952.1", 952, 952, "Mysterium obscurum", "", "", "latest", "Contig", "MO1", 200000000, 5000, ""),
]


def taxdump(path):
    with tarfile.open(path, "w:gz") as tar:
        for name, lines in (("nodes.dmp", [f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\n" for t, p, r, _ in NODES]),
                            ("names.dmp", [f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, _, _, n in NODES]
                             + [f"{t}\t|\t{n}\t|\t\t|\t{c}\t|\n" for t, n, c in SYNONYMS])):
            data = "".join(lines).encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def world(folder, seed=3):
    """The EukDetect2 folder, KSGP's tax file, an assembly table and a GTDB release folder in `folder`."""
    rng = random.Random(seed)
    eukdetect = os.path.join(folder, "eukdetect2")
    os.makedirs(eukdetect)
    taxdump(os.path.join(eukdetect, "taxdump.tar.gz"))
    with open(os.path.join(eukdetect, "busco_taxid_genome_link.txt"), "w") as link, \
            open(os.path.join(eukdetect, "eukdb.fasta"), "w") as fasta:
        for genome, (taxid, group, loci) in GENOMES.items():
            for locus in loci:
                record = f"{group}-{genome}-{locus}-S1"
                link.write(f"{record}\t{taxid}\t{genome}\n")
                seq = "".join(rng.choice("ACGT") for _ in range(400))
                if genome == "Kluyveromyces_lactis" and locus == "3at2759":
                    seq = "A" * 17000  # longer than MAX_LOCUS: left out
                if genome == "Newgenus_novus" and locus == "1at2759":
                    seq = "NNNN" + seq + "NN"  # Ns at its ends: cut
                fasta.write(f">{record}\n{seq}\n")
                if genome == "Candida_albicans_SC5314" and locus == "1at2759":  # a second, shorter piece
                    link.write(f"{group}-{genome}-{locus}-S2\t{taxid}\t{genome}\n")
                    fasta.write(f">{group}-{genome}-{locus}-S2\n{seq[:100]}\n")
    ksgp = os.path.join(folder, "ksgp.tax")
    with open(ksgp, "w") as fh:
        fh.writelines(f"R{i}\t{label}\n" for i, label in enumerate(KSGP))
    assemblies = os.path.join(folder, "assembly_summary_fungi.txt")
    with open(assemblies, "w") as fh:
        fh.write("##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n")
        fh.write("#" + "\t".join(ASSEMBLY_COLUMNS) + "\n")
        for row in ASSEMBLIES:
            fh.write("\t".join(map(str, row)) + f"\thttps://ftp.example/{row[0]}\n")
    gtdb = os.path.join(folder, "gtdb_226")
    os.makedirs(gtdb)
    with open(os.path.join(gtdb, "bac120_taxonomy_r226.tsv"), "w") as fh:
        fh.write("RS_GCF_000000001.1\td__Bacteria;p__P;c__C;o__O;f__F;g__Newgenus;s__Newgenus bacterialis\n")
    return eukdetect, ksgp, assemblies, gtdb


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def read_fasta_file(path):
    """[(header, sequence)] of a FASTA of one sequence line per record, in its order."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 2)]


MAP1 = {"1at2759": "901at2759", "2at2759": "902at2759", "3at2759": "903at2759"}  # euk30 locus -> EukDetect1 locus


def eukdetect1_world(folder, eukdetect):
    """EukDetect1's two files in folder/eukdetect1: S. cerevisiae S288C's copies as EukDetect2 has them (under
    EukDetect1's loci, for map_loci); Ascaris lumbricoides with two mapped loci single, one duplicated and one unmapped;
    Mysterium rarum with one; a collapsed record."""
    genomes, records = eukaryotes.read_link(os.path.join(eukdetect, "busco_taxid_genome_link.txt"))
    copies = eukaryotes.read_loci(os.path.join(eukdetect, "eukdb.fasta"), records)
    rng = random.Random(9)
    out = os.path.join(folder, "eukdetect1")
    os.makedirs(out)
    entries = [(f"fungi-Saccharomyces_cerevisiae_S288C-{MAP1[l]}-S1", 559292, s)
               for l, s in copies["Saccharomyces_cerevisiae_S288C"].items()]
    for locus, copy in (("901at2759", "S1"), ("902at2759", "S1"), ("903at2759", "D1"), ("903at2759", "D2"),
                        ("999at2759", "S1")):
        entries.append((f"worm-Ascaris_lumbricoides-{locus}-{copy}", 6253, "".join(rng.choice("ACGT") for _ in range(300))))
    entries.append(("protist-Mysterium_rarum-901at2759-S1", 970, "ACGT" * 75))
    entries.append(("other_SGCollapse_Mysterium-902at2759-S1", 951, "ACGT" * 75))
    with open(os.path.join(out, "busco_taxid_link.txt"), "w") as link, \
            open(os.path.join(out, "ncbi_eukprot_met_arch_markers.fna"), "w") as fasta:
        for record, taxid, seq in entries:
            link.write(f"{record}\t{taxid}\n")
            fasta.write(f">{record}\n{seq}\n")
    return out


