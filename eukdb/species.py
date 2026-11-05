#!/usr/bin/env python3
"""model_species.py - the model and host species protal makes an extra effort for (since 2026-11-05; docs/databases.md).

EukDetect2 and EukDetect1, euk30's sources, hold no vertebrate and few plants: a host's reads (human, mouse, a farm
animal's) were unknown microbes in a profile, and the species euk30 has from EukDetect1 (Drosophila melanogaster,
Caenorhabditis elegans) are on 11-16 of its 30 loci. model_species.tsv (beside this script) names them, with their
reference assembly and a curated GenBank rDNA unit (model_rdna_units.fna.gz). For each eukaryote of it,
eukdetect1_genomic.py --models places the euk30 loci its euk30 copies lack in its assembly, whatever the assembly's size
(the download's caps, --euk_max_genome_mb and --rdna_max_genome_mb, are for the ordinary eukaryotes):

- the proteins: a locus's copies of the species' nearest relatives (NCBI's lineage: the deepest shared ancestor first;
  euk30_exons.Relatives), 4 of them: EukDetect1's coding sequences (the 16 loci euk30_loci.tsv maps to EukDetect1's),
  EukDetect2's genomic copies that are one open frame (orf_cds; 36% of them) and the placed model species' coding
  sequences (their exons joined, so that a mammal placed after human and mouse takes theirs), translated;
- placing: miniprot, as for EukDetect1's genomes (eukdetect1_genomic.place); of a locus's proteins the alignment of the
  highest identity times coverage, taken at MIN_IDENTITY and MIN_COVERAGE (a relative's protein, not the species' own:
  BUSCO's single-copy genes of a vertebrate are 50-75% identical to a lancelet's, 35-60% to a fungus's);
- the copy: the genome from the alignment's first codon to its last, introns and all; a gene longer than
  eukaryotes.MAX_LOCUS bases (many a vertebrate's) fitted into it by cutting its introns from their middles, every exon
  kept (euk30_copies.fit; window: the run of its exons with the most coding bases, where they alone exceed it).

eukaryotes.build_release adds a species euk30 lacks on its placed loci (MIN_EUKDETECT1_LOCI at least), its reference
assembly its representative (assembly_row: the table's accession, size and level, no other GenBank assembly of the
species), and gives a species euk30 has the loci its representative lacks. The assembly's rRNA goes into the download's
store (rdna_store.py) under its accession, which the freeze reads as a euk30 genome's; where it gives no operon (GRCm39
holds few rDNA units) the curated unit's rRNA is kept with it (unit_fasta). The prokaryotes of the table are checked only
(their GTDB species, 16S and operons).
"""

import collections
import gzip
import os
import re

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TABLE = os.path.join(DATA, "model_species.tsv")
UNITS = os.path.join(DATA, "model_rdna_units.fna.gz")
GROUP = "model"            # the placements' group (eukdetect1_genomic.py's table), the records' prefix
MODEL_FASTA = "model_species.genomic.fna"  # eukdetect1_genomic.py --models' copies and table in its --out folder
MODEL_TABLE = "model_species_placed.tsv"
SOURCE = "model species"   # euk30's eukdetect_source and euk_match of a species added from the table
MIN_COVERAGE = 0.7         # of the protein, by its best alignment
MIN_IDENTITY = 0.35        # miniprot's Identity of it
MIN_WINDOW_SHARE = 0.3     # a gene cut to eukaryotes.MAX_LOCUS bases: its coding bases the window holds, at least
MIN_ORF_SHARE = 0.9       # an EukDetect2 copy as one open frame: the frame's codons this share of the copy at least
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
CODE = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_table(path=TABLE):
    """model_species.tsv -> [row dict], the comment lines left out."""
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    return [dict(zip(lines[0], f)) for f in lines[1:]]


def eukaryotes(rows):
    return [r for r in rows if r["domain"] == "Eukaryota"]


def genome_name(row):
    """The placements' genome name of a model species: model_<Genus_species>, apart from EukDetect1's genome names."""
    return "model_" + re.sub(r"[^A-Za-z0-9_.-]", "_", row["species"])


def record(row, locus):
    """A placed copy's record: model-<Genus_species>-<locus>-S1, EukDetect1's form (eukaryotes.RECORD_RE)."""
    return f"{GROUP}-{re.sub(r'[^A-Za-z0-9_.]', '_', row['species'])}-{locus}-S1"


def accessions(row):
    """The accessions under which a download may hold a model species' assembly: its GenBank and RefSeq ones."""
    return [a for a in (row["genbank"], row.get("refseq", "-")) if a and a != "-"]


def ftp_path(accession, name):
    """An assembly's folder on NCBI's FTP server."""
    digits = re.match(r"(GC[AF])_(\d{3})(\d{3})(\d{3})", accession)
    return f"{FTP}/{digits.group(1)}/{digits.group(2)}/{digits.group(3)}/{digits.group(4)}/{accession}_{name}"


def assembly_row(row):
    """A model species' reference assembly as eukaryotes.read_assemblies writes a row; "model" marks it (no size cap)."""
    return dict(accession=row["genbank"], taxid=int(row["ncbi_taxid"]), species_taxid=int(row["ncbi_taxid"]),
                organism=row["species"], strain="", isolate="", level=row["level"], asm_name=row["assembly_name"],
                size=int(row["genome_bases"]) if row["genome_bases"].isdigit() else None, contigs=None,
                category="isolate", model=True,
                ftp=row["ftp"] if row.get("ftp", "-") not in ("", "-", "na") else ftp_path(row["genbank"],
                                                                                         row["assembly_name"]),
                rdna_unit=row.get("rdna_unit", "-"))


def translate(cds):
    return "".join(CODE.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def orf_cds(seq):
    """A genomic copy that is one open frame (no stop but a last one, on either strand) -> its coding bases from the
    frame's first codon (MIN_ORF_SHARE of the copy at least), else None: a copy with introns has none."""
    seq = seq.upper()
    found = []
    for strand in (seq, seq.translate(COMPLEMENT)[::-1]):
        for frame in range(3):
            cds = strand[frame:]
            cds = cds[:len(cds) - len(cds) % 3]
            protein = translate(cds)
            body = protein[:-1] if protein.endswith("*") else protein
            if "*" not in body and "N" not in cds and 3 * len(body) >= MIN_ORF_SHARE * len(seq):
                found.append((not protein.endswith("*"), cds))  # a frame that ends at a stop codon first
    return min(found, key=lambda f: f[0])[1] if found else None


def locus_of(query):
    """A query record's copy record (record()), its locus being the fourth part of it."""
    return query.split("~", 1)[0]


def taken(identity, coverage, share):
    """Whether a model species' placed copy is taken: its protein's alignment at MIN_IDENTITY and MIN_COVERAGE, and the
    copy (a gene's window where it is longer than eukaryotes.MAX_LOCUS) holding MIN_WINDOW_SHARE of its coding bases."""
    return coverage >= MIN_COVERAGE and identity >= MIN_IDENTITY and share >= MIN_WINDOW_SHARE


def best(results):
    """eukdetect1_genomic.place's results of a model species' queries ({query: (copy, fields[, features])}) -> {copy
    record: (copy, fields + [the relative], features)}: per locus the placed query of the highest identity times
    coverage, its status by taken(); a locus none of whose queries was placed "not placed", or the failure of its run
    (miniprot failed: placed again by a later run), without the relative."""
    by_locus = collections.defaultdict(list)
    for query, value in results.items():
        copy, fields, features = (tuple(value) + (None,))[:3]
        by_locus[locus_of(query)].append((query, copy, fields, features))
    out = {}
    for rec, hits in sorted(by_locus.items()):
        placed = [(float(f[8]) * float(f[9]), q, c, f, x) for q, c, f, x in hits if c is not None]
        if not placed:
            out[rec] = (None, next((f for _q, _c, f, _x in hits if f[0] != "not placed"), ["not placed"]), None)
            continue
        _score, query, copy, fields, features = max(placed, key=lambda p: (p[0], p[1]))
        fields = list(fields)
        fields[0] = "genomic" if taken(float(fields[8]), float(fields[9]), float(fields[10])) else "poor placement"
        out[rec] = (copy, fields + [query.split("~", 2)[2]], features)
    return out


def window(exons, max_span):
    """Of a gene's exons [(start, end)] in genome order, the run with the most coding bases whose span (the first's start
    to the last's end) is max_span bases or fewer -> that run (all of them where the gene fits)."""
    best, best_bases = exons[:1], 0
    for i in range(len(exons)):
        bases = 0
        for j in range(i, len(exons)):
            if exons[j][1] - exons[i][0] + 1 > max_span:
                break
            bases += exons[j][1] - exons[j][0] + 1
            if bases > best_bases:
                best, best_bases = exons[i:j + 1], bases
    return best


def unit_fasta(accession, out, units=UNITS):
    """A curated rDNA unit of model_rdna_units.fna.gz written as FASTA into `out` -> out, or None where it has none."""
    name = None
    with gzip.open(units, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif name == accession:
                with open(out, "w", newline="\n") as fa:
                    fa.write(f">{accession}\n{line.strip()}\n")
                return out
    return None


# The RefSeq species euk30 lacks or has on few loci (since 2026-11-05: "the RefSeq tier"; download_ref_data.py): their
# RefSeq reference or representative genome joins the download's genome set, for their euk30 loci (placed as the model
# species' are) and their rRNA. RefSeq's eukaryote groups and the euk_group each is in.
REFSEQ_GROUPS = {"vertebrate_mammalian": "Metazoa", "vertebrate_other": "Metazoa", "invertebrate": "Metazoa",
                 "plant": "Archaeplastida", "fungi": "Fungi", "protozoa": "Protist"}
GOOD_LOCI = 20  # a species of euk30 with fewer of the 30 loci has its RefSeq genome placed for the others
TABLE_COLUMNS = ("species", "ncbi_taxid", "domain", "group", "genbank", "refseq", "assembly_name", "level",
                 "genome_bases", "rdna_unit", "note", "ftp")
CATEGORY_RANK = {"reference genome": 0, "representative genome": 1}
LEVEL_RANK = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def refseq_genomes(path, group):
    """A RefSeq assembly_summary.txt -> {species taxid: row (model species table columns)} of each species' reference
    genome, else its representative one (then the better assembly level, fewer contigs), the latest version."""
    best = {}
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                continue
            f = dict(zip(header or (), line.rstrip("\n").split("\t")))
            if f.get("version_status") != "latest" or f.get("refseq_category") not in CATEGORY_RANK:
                continue
            key = (CATEGORY_RANK[f["refseq_category"]], LEVEL_RANK.get(f.get("assembly_level"), 4),
                   int(f["contig_count"]) if f.get("contig_count", "").isdigit() else 1 << 30, f["assembly_accession"])
            taxid = int(f["species_taxid"])
            if taxid in best and best[taxid][0] <= key:
                continue
            paired = f.get("gbrs_paired_asm", "na")
            best[taxid] = (key, {"species": re.sub(r"\s+", " ", f.get("organism_name", "")).strip(),
                                 "ncbi_taxid": str(taxid), "domain": "Eukaryota", "group": REFSEQ_GROUPS[group],
                                 "genbank": paired if paired.startswith("GCA_") else f["assembly_accession"],
                                 "refseq": f["assembly_accession"], "assembly_name": f.get("asm_name", ""),
                                 "level": f.get("assembly_level", ""), "genome_bases": f.get("genome_size", "-"),
                                 "rdna_unit": "-", "note": f"RefSeq {f['refseq_category']} ({group})",
                                 "ftp": f.get("ftp_path", "").replace("ftp://", "https://")})
    return {t: row for t, (_key, row) in best.items()}


def gap_rows(refseq, have, ncbi, exclude=(), good=GOOD_LOCI):
    """The RefSeq species ({species taxid: row}) euk30's sources lack or have on fewer than `good` loci (have: {NCBI
    species taxid: loci}), those of the model species table (exclude: taxids) left out -> rows, by species name. The
    species name is NCBI's (the organism name without its strain)."""
    out = []
    excluded = {ncbi.species_of(int(t)) for t in exclude}
    for taxid, row in refseq.items():
        sp = ncbi.species_of(taxid)
        if sp in excluded or len(have.get(sp, ())) >= good or ncbi.rank.get(sp) != "species":
            continue
        out.append({**row, "species": ncbi.name.get(sp, row["species"]), "ncbi_taxid": str(sp)})
    return sorted(out, key=lambda r: (r["species"], r["ncbi_taxid"]))


def write_table(path, rows):
    """Rows in the model species table's format (TABLE_COLUMNS), under the writer's own name first."""
    from .util import part_of  # noqa: PLC0415
    with open(part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(TABLE_COLUMNS) + "\n")
        fh.writelines("\t".join(str(r.get(c, "-") or "-") for c in TABLE_COLUMNS) + "\n" for r in rows)
    os.replace(part_of(path), path)


def host_lineage(host_path, names_path, table=TABLE):
    """The euk30 lineage of a host genome's species (since 2026-11-05: a host scenario's samples then hold it, in their
    truth): the table's species whose GenBank or RefSeq accession names the host's FASTA (OUT/host/<accession>.fna.gz)
    and that euk30_names_r<R>.tsv lists by its NCBI species taxid -> its lineage, or None (another host, a release
    without it)."""
    if not host_path or not names_path or not os.path.isfile(names_path):
        return None
    base = os.path.basename(host_path)
    row = next((r for r in read_table(table) if any(base.startswith(a + ".") or base.startswith(a + "_")
                                                     for a in accessions(r))), None)
    if row is None:
        return None
    with open(names_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(header, line.rstrip("\n").split("\t")))
            if f.get("ncbi_species_taxid") == row["ncbi_taxid"]:
                return f.get("lineage")
    return None


def copies_of(paths):
    """The model species' placed copies in eukdetect1_genomic's FASTA files (plain or .zst; those that exist) ->
    {genome name: {locus: copy}}, a later file's in place of an earlier's."""
    from . import sources as eukaryotes  # noqa: PLC0415
    out = collections.defaultdict(dict)
    for path in paths:
        if not path or not os.path.isfile(path):
            continue
        with eukaryotes.open_text(path) as fh:
            for name, seq in zip(*[iter(line.strip() for line in fh)] * 2):
                m = eukaryotes.RECORD_RE.match(name[1:])
                if m and m.group("group") == GROUP:
                    out[f"model_{m.group('genome')}"][m.group("locus")] = seq
    return out
