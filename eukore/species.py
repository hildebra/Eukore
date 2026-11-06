"""The model and host species Eukore makes an extra effort for (data/model_species.tsv, since 2026-11-05; from protal's
model_species.py): human, mouse, rat, zebrafish, cow, pig, chicken, Drosophila, C. elegans, yeast and Arabidopsis.
Each takes its reference assembly whatever genomes.rank would choose, and its curated rDNA unit
(data/model_rdna_units.fna.gz) where its assembly gives no operon. The table's prokaryotes are checked by protal, not
held by Eukore.
"""

import gzip
import os
import re

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
TABLE = os.path.join(DATA, "model_species.tsv")
UNITS = os.path.join(DATA, "model_rdna_units.fna.gz")
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def read_table(path=TABLE):
    """model_species.tsv -> [row dict], the comment lines left out."""
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    return [dict(zip(lines[0], f)) for f in lines[1:]]


def eukaryotes(rows):
    return [r for r in rows if r["domain"] == "Eukaryota"]


def accessions(row):
    """The accessions under which a user may hold a model species' assembly: its GenBank and RefSeq ones."""
    return [a for a in (row["genbank"], row.get("refseq", "-")) if a and a != "-"]


def ftp_path(accession, name):
    """An assembly's folder on NCBI's FTP server."""
    digits = re.match(r"(GC[AF])_(\d{3})(\d{3})(\d{3})", accession)
    return f"{FTP}/{digits.group(1)}/{digits.group(2)}/{digits.group(3)}/{digits.group(4)}/{accession}_{name}"


def assembly_row(row):
    """A model species' reference assembly as sources.read_assemblies writes a row."""
    return dict(accession=row["genbank"], taxid=int(row["ncbi_taxid"]), species_taxid=int(row["ncbi_taxid"]),
                organism=row["species"], strain="", isolate="", level=row["level"], asm_name=row["assembly_name"],
                size=int(row["genome_bases"]) if row["genome_bases"].isdigit() else None, contigs=None,
                category="isolate", model=True, flags=[], refseq_category="reference genome",
                ftp=row["ftp"] if row.get("ftp", "-") not in ("", "-", "na") else ftp_path(row["genbank"],
                                                                                         row["assembly_name"]),
                rdna_unit=row.get("rdna_unit", "-"))


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
