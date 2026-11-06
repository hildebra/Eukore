"""A small world for EUKDB's tests (since 2026-11-05: NCBI's genomes, no EukDetect): NCBI's taxonomy of a few fungi,
an unknown phylum, a nematode and human; NCBI's assembly tables; a stand-in BUSCO dataset (odb.py's layout) of three
loci; a genome holding two genes; stand-in miniprot, hmmsearch and vsearch scripts (Python, as the tools write)."""

import gzip
import io
import os
import random
import stat
import tarfile

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
    (27291, 4930, "species", "Saccharomyces paradoxus"),
    (1080349, 4930, "species", "Saccharomyces cerevisiae x Saccharomyces paradoxus"),  # NCBI's hybrid: no species
    (4910, 4893, "genus", "Kluyveromyces"),
    (28985, 4910, "species", "Kluyveromyces lactis"),
    (1535326, 4892, "family", "Debaryomycetaceae"),
    (5475, 1535326, "genus", "Candida"),
    (5476, 5475, "species", "Candida albicans"),
    (5204, 4751, "phylum", "Basidiomycota"),
    (155616, 5204, "class", "Tremellomycetes"),
    (5234, 155616, "order", "Tremellales"),
    (1884633, 5234, "family", "Cryptococcaceae"),
    (5206, 1884633, "genus", "Cryptococcus"),
    (40410, 5206, "species", "Cryptococcus deneoformans"),
    (950, 2759, "phylum", "Mysteriozoa"),
    (951, 950, "genus", "Mysterium"),
    (952, 951, "species", "Mysterium obscurum"),
    (953, 951, "species", "uncultured Mysterium"),
    (33208, 2759, "kingdom", "Metazoa"),
    (6231, 33208, "phylum", "Nematoda"),
    (6251, 6231, "genus", "Ascaris"),
    (6253, 6251, "species", "Ascaris lumbricoides"),
    (7711, 33208, "phylum", "Chordata"),
    (9605, 7711, "genus", "Homo"),
    (9606, 9605, "species", "Homo sapiens"),
]
SYNONYMS = [(28985, "Candida sphaerica", "anamorph")]


CODES = {28985: 6}  # a stand-in: Kluyveromyces lactis given the ciliates' code, to test that it is read


def taxdump(path):
    with tarfile.open(path, "w:gz") as tar:
        # nodes.dmp: taxid, parent, rank, EMBL code, division, inherited, genetic code (a ciliate's 6, else 1), ...
        for name, lines in (("nodes.dmp", [f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\t0\t|\t1\t|\t{CODES.get(t, 1)}\t|\t1\t|\n"
                                           for t, p, r, _ in NODES]),
                            ("names.dmp", [f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, _, _, n in NODES]
                             + [f"{t}\t|\t{n}\t|\t\t|\t{c}\t|\n" for t, n, c in SYNONYMS])):
            data = "".join(lines).encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


ASSEMBLY_COLUMNS = ("assembly_accession", "refseq_category", "taxid", "species_taxid", "organism_name",
                    "infraspecific_name", "isolate", "version_status", "assembly_level", "seq_rel_date", "asm_name",
                    "gbrs_paired_asm", "ftp_path", "excluded_from_refseq", "genome_size", "contig_count")


def assembly(acc, taxid, species, organism, level="Scaffold", contigs=50, size=12_000_000, category="na",
             paired="na", flags="na", date="2020/01/01", status="latest"):
    return dict(assembly_accession=acc, refseq_category=category, taxid=taxid, species_taxid=species,
                organism_name=organism, infraspecific_name="", isolate="", version_status=status,
                assembly_level=level, seq_rel_date=date, asm_name=acc.replace(".", "_"), gbrs_paired_asm=paired,
                ftp_path=f"https://ftp.example/{acc}", excluded_from_refseq=flags, genome_size=size,
                contig_count=contigs)


GENBANK = [
    assembly("GCA_000146045.2", 559292, 4932, "Saccharomyces cerevisiae S288C", "Complete Genome", 17,
             category="reference genome", paired="GCF_000146045.2"),
    assembly("GCA_000000003.1", 4932, 4932, "Saccharomyces cerevisiae", "Contig", 300),
    assembly("GCA_000000011.1", 27291, 27291, "Saccharomyces paradoxus", "Contig", 900, flags="contaminated"),
    assembly("GCA_000000012.1", 27291, 27291, "Saccharomyces paradoxus", "Scaffold", 400,
             flags="derived from metagenome"),
    assembly("GCA_000000013.1", 27291, 27291, "Saccharomyces paradoxus", "Contig", 2000),
    assembly("GCA_000002515.1", 28985, 28985, "Kluyveromyces lactis", "Chromosome", 6),
    assembly("GCA_000002516.1", 28985, 28985, "Kluyveromyces lactis", "Chromosome", 6, date="2024/05/01"),
    assembly("GCA_000182965.3", 5476, 5476, "Candida albicans SC5314", "Chromosome", 8),
    assembly("GCA_000000952.1", 952, 952, "Mysterium obscurum", "Contig", 5000, size=200_000_000),
    assembly("GCA_000000953.1", 953, 953, "uncultured Mysterium", "Contig", 100),  # a placeholder: no species
    assembly("GCA_000000954.1", 1080349, 1080349, "Saccharomyces cerevisiae x Saccharomyces paradoxus", "Contig", 100),
    assembly("GCA_000000560.1", 4932, 4932, "Saccharomyces cerevisiae", "Complete Genome", 1, size=560),
]
REFSEQ = [assembly("GCF_000146045.2", 559292, 4932, "Saccharomyces cerevisiae S288C", "Complete Genome", 17,
                   category="reference genome", paired="GCA_000146045.2")]


def assembly_table(path, rows):
    with open(path, "w") as fh:
        fh.write("##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt\n")
        fh.write("#" + "\t".join(ASSEMBLY_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in ASSEMBLY_COLUMNS) + "\n")
    return path


# The stand-in BUSCO dataset: three loci, each a protein a test genome's gene encodes; species 4932 and 6253 hold them
LOCI = ["1at2759", "2at2759", "3at2759"]
CODE = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA",
        "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
        "W": "TGG", "Y": "TAT"}


def protein(seed, n=120):
    rng = random.Random(seed)
    return "M" + "".join(rng.choice("ACDEFGHIKLNPQRSTVWY") for _ in range(n - 1))


PROTEINS = {locus: protein(i + 1) for i, locus in enumerate(LOCI)}


def odb_world(folder):
    """The dataset folder (odb.py's layout) -> its path."""
    root = os.path.join(folder, "eukaryota_odb12")
    os.makedirs(os.path.join(root, "hmms"))
    with open(os.path.join(root, "scores_cutoff"), "w") as fh:
        fh.writelines(f"{locus}\t50.0\n" for locus in LOCI)
    for locus in LOCI:
        with open(os.path.join(root, "hmms", f"{locus}.hmm"), "w") as fh:
            fh.write(f"HMMER3/f [3.4]\nNAME  {locus.split('at')[0]}\nLENG  120\nHMM          A        C\n//\n")
    with gzip.open(os.path.join(root, "refseq_db.faa.gz"), "wt") as fh:
        for locus, p in PROTEINS.items():
            for taxid in (4932, 6253, 9606):
                fh.write(f">{locus}_{taxid}_0:00000{taxid % 7}\n{p}\n")
    with open(os.path.join(root, "ancestral_variants"), "w") as fh:
        fh.writelines(f">{locus}_0\n{p}\n" for locus, p in PROTEINS.items())
    return root


def cds_of(p):
    return "".join(CODE[a] for a in p) + "TAA"


def genome(seed=5):
    """A contig 'c1' with locus 1's gene on + (two exons, an intron of 300 with a soft-masked run and a poly-A), locus
    2's gene on - (one exon), and locus 1 again (a second copy) further on; and the GFF miniprot writes for it.
    -> (contig sequence, GFF text, {locus: (start, end)} 1-based)."""
    rng = random.Random(seed)

    def bases(n):
        return "".join(rng.choice("ACGT") for _ in range(n))
    c1, c2 = cds_of(PROTEINS["1at2759"]), cds_of(PROTEINS["2at2759"])
    intron = "GT" + bases(60) + "a" * 80 + "A" * 70 + bases(86) + "AG"  # 300: a soft-masked run, a low-complexity run
    seq = bases(1000)
    g1 = (len(seq) + 1, len(seq) + len(c1) + len(intron))
    seq += c1[:180] + intron + c1[180:]
    seq += bases(2000)
    rc = c2.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    g2 = (len(seq) + 1, len(seq) + len(rc))
    seq += rc + bases(3000)
    g3 = (len(seq) + 1, len(seq) + len(c1))
    seq += c1 + bases(1000)
    lines = ["##gff-version 3"]
    n = 0
    for locus, (a, b), strand, exons in (
            ("1at2759", g1, "+", [(g1[0], g1[0] + 179), (g1[0] + 180 + len(intron), g1[1])]),
            ("2at2759", g2, "-", [g2]),
            ("1at2759", g3, "+", [g3])):
        for taxid in (4932, 6253):
            n += 1
            q = f"{locus}~{locus}_{taxid}_0:00000{taxid % 7}"
            score = 900 if taxid == 4932 else 700
            lines.append(f"##PAF\t{q}\t120")
            lines.append(f"##STA\t{PROTEINS[locus]}")  # miniprot's order: PAF, STA, then the alignment
            lines.append(f"c1\tminiprot\tmRNA\t{a}\t{b}\t{score}\t{strand}\t.\tID=MP{n};Rank=1;Identity=0.9500;"
                         f"Positive=0.97;Target={q} 1 120")
            for x, y in exons:
                lines.append(f"c1\tminiprot\tCDS\t{x}\t{y}\t{score}\t{strand}\t0\tParent=MP{n};Rank=1;"
                             f"Identity=0.95;Target={q} 1 120")
    return seq, "\n".join(lines) + "\n", {"1at2759": g1, "2at2759": g2, "1at2759b": g3}


def tool(folder, name, text):
    path = os.path.join(folder, name)
    with open(path, "w", newline="\n") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
    return path


def miniprot(folder, gff):
    """A miniprot that prints the GFF (its arguments ignored)."""
    path = os.path.join(folder, "mp.gff")
    with open(path, "w") as fh:
        fh.write(gff)
    return tool(folder, "miniprot", f"#!/bin/sh\ncat '{path}'\n")


HMMSEARCH = '''#!/usr/bin/env python3
"""hmmsearch's stand-in: a protein that is a locus' (fixtures.PROTEINS), or its first 30 residues or more, scores 300
against it, 20 against the others; its domain covers the profile in proportion to its share of the locus' protein."""
import sys
args = sys.argv[1:]
dom = args[args.index("--domtblout") + 1]
faa = args[-1]
proteins = {}
name = None
for line in open(faa):
    if line.startswith(">"):
        name = line[1:].split()[0]
        proteins[name] = ""
    else:
        proteins[name] += line.strip()
LOCI = %r
with open(dom, "w") as out:
    out.write("# target\\n")
    for name, seq in proteins.items():
        for locus, p in LOCI.items():
            score = 300.0 if len(seq) >= 30 and p.startswith(seq) else 20.0
            to = max(1, int(120 * min(len(seq), len(p)) / len(p)))
            f = [name, "-", str(len(seq)), locus.split("at")[0], "-", "120", "1e-50", str(score), "0", "1", "1",
                 "0", "0", str(score), "0", "1", str(to), "1", str(len(seq)), "1", str(len(seq)), "0.99", "-"]
            out.write(" ".join(f) + "\\n")
'''


def hmmsearch(folder):
    return tool(folder, "hmmsearch", HMMSEARCH % PROTEINS)


VSEARCH = '''#!/usr/bin/env python3
"""vsearch --usearch_global's stand-in: the identity of two sequences is their share of equal positions (100 - 2 per
difference for those of one length), written to --userout (query, target, id), self hits left out."""
import sys
args = sys.argv[1:]
q, db, out = (args[args.index(o) + 1] for o in ("--usearch_global", "--db", "--userout"))


def read(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs


seqs, targets = read(q), read(db)
with open(out, "w") as fh:
    for a, s in seqs.items():
        for b, t in targets.items():
            if a == b:
                continue
            n = max(len(s), len(t))
            same = sum(x == y for x, y in zip(s, t))
            ident = 100.0 * same / n
            if ident >= 70:
                fh.write(f"{a}\\t{b}\\t{ident:.1f}\\n")
'''


def vsearch(folder):
    return tool(folder, "vsearch", VSEARCH)


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
