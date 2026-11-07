#!/usr/bin/env python3
"""data/primer_references.fna (since 2026-11-07): the records each primer pair of data/primer_pairs.tsv amplifies in
tests/test_primers.py, and the S. cerevisiae 18S whose coordinates data/rrna_regions.tsv gives. Made from NCBI's
records (fetched by the user with efetch, as below), the organelle genes cut by eukore/organelle.py's own parser:

    Sc_18S        Z75578.1 (S. cerevisiae NCYC 505 18S; Hadziavdic et al. 2014's coordinates)
    Sc_rDNA       NR_132207.1 (S. cerevisiae S288C 35S pre-rRNA, data/model_rdna_units.fna.gz)
    Bm_COX1       NC_002355.1's COX1 (Bombyx mori mitochondrion: Folmer's primers fit it)
    Dr_COX1       NC_002333.2's COX1 (Danio rerio mitochondrion)
    Dr_12S, Dr_16S  NC_002333.2's SSU and LSU
    At_rbcL, At_matK  NC_000932.1's (Arabidopsis thaliana plastid)
    Nt_trnH-psbA  NC_001879.2's (Nicotiana tabacum plastid: Arabidopsis' trnH holds an indel under trnHf_05)

    efetch: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<accession>&rettype=gbwithparts
    python3 tools/primer_references.py --z75578 Z75578.fna --bm NC_002355.gb --dr NC_002333.gb --at NC_000932.gb --nt NC_001879.gb
"""

import argparse
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import organelle  # noqa: E402

DATA = os.path.join(os.path.dirname(HERE), "data")


def fasta(path):
    out, name = {}, None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = []
            elif name:
                out[name].append(line.strip())
    return {n: "".join(p).upper() for n, p in out.items()}


def genes(path, kind, animal=False):
    (rec,) = list(organelle.gbff_records(path))
    return organelle.record_genes(rec, kind, animal)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("z75578", "bm", "dr", "at", "nt"):
        p.add_argument(f"--{name}", required=True)
    opts = p.parse_args(argv)
    bm, dr = genes(opts.bm, "mitochondrion", True), genes(opts.dr, "mitochondrion", True)
    at, nt = genes(opts.at, "plastid"), genes(opts.nt, "plastid")
    records = {"Sc_18S": next(iter(fasta(opts.z75578).values())),
               "Sc_rDNA": fasta(os.path.join(DATA, "model_rdna_units.fna.gz"))["NR_132207.1"],
               "Bm_COX1": bm["COX1"]["sequence"], "Dr_COX1": dr["COX1"]["sequence"],
               "Dr_12S": dr["SSU"]["sequence"], "Dr_16S": dr["LSU"]["sequence"],
               "At_rbcL": at["rbcL"]["sequence"], "At_matK": at["matK"]["sequence"],
               "Nt_trnH-psbA": nt["trnH-psbA"]["sequence"]}
    path = os.path.join(DATA, "primer_references.fna")
    with open(path, "w", newline="\n") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")
    print({n: len(s) for n, s in records.items()}, "regions of trnH-psbA:", nt["trnH-psbA"]["regions"])


if __name__ == "__main__":
    main()
