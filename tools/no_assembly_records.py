#!/usr/bin/env python3
"""data/no_assembly.fna.gz: the euk30 copies of the species Eukore finds no NCBI assembly for, from EukDetect2 (Shih,
Zhao, Pollard & Lind 2026, bioRxiv 10.64898/2026.06.24.734308; Zenodo 19056625, CC BY 4.0), the one remnant of EukDetect
in Eukore (since 2026-11-05). EukDetect2 took them from EukProt's genomes and transcriptomes, and re-assembled
Dientamoeba fragilis' transcriptome; NCBI has no assembly of them. The build calls the loci in these records as in a
genome (loci.py), which gives them exon models and the same checks as every other copy.

    python3 tools/no_assembly_records.py --eukdetect2 DIR --assemblies WORK/sources/ncbi_assemblies \\
        --taxdump WORK/sources/taxdump.tar.gz [--out data/no_assembly.fna.gz]

Records: >species taxid|organism|EukDetect2 record (its locus in it); one per EukDetect2 record of the species, of
every BUSCO EukDetect2 holds (the caller keeps the loci it calls). Records naming no species (EukDetect2's copies
collapsed across species) are left out.
"""

import argparse
import collections
import glob
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import genomes, ncbi as ncbi_mod, odb as odb_mod, sources  # noqa: E402

LOCUS_RE = re.compile(r"-(\d+at\d+)-S\d+$")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eukdetect2", required=True, help="EukDetect2's eukdb.fasta and busco_taxid_genome_link.txt")
    p.add_argument("--assemblies", required=True, help="NCBI's assembly tables (genbank_*.txt, refseq_*.txt)")
    p.add_argument("--taxdump", required=True)
    p.add_argument("--out", default=os.path.join(sources.DATA, "no_assembly.fna.gz"))
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(opts.taxdump)
    # the species the build takes an assembly of (genomes.representatives: an excluded one is none; until 2026-11-06
    # any assembly counted, and a species whose only one is flagged got neither a genome nor these records)
    have = set(genomes.representatives(
        sources.read_assemblies(sorted(glob.glob(os.path.join(opts.assemblies, "genbank_*.txt")))), ncbi))
    owner = {}
    with open(os.path.join(opts.eukdetect2, "busco_taxid_genome_link.txt")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and LOCUS_RE.search(f[0]) and f[1].isdigit():
                owner[f[0]] = ncbi.species_of(int(f[1]))
    wanted = {r: sp for r, sp in owner.items() if sp not in have and ncbi.rank.get(sp) == "species"}
    n, species = 0, collections.Counter()
    with gzip.open(opts.out + ".part", "wt", compresslevel=9) as out:
        # every line of a record (until 2026-11-06 its first alone: a wrapped FASTA's records cut)
        for name, seq in odb_mod.read_fasta(os.path.join(opts.eukdetect2, "eukdb.fasta")):
            if name not in wanted:
                continue
            sp = wanted[name]
            organism = re.sub(r"\s+", "_", ncbi.name.get(sp, str(sp)))
            out.write(f">{sp}|{organism}|{name}\n{seq.upper()}\n")
            n += 1
            species[sp] += 1
    os.replace(opts.out + ".part", opts.out)
    print(f"{opts.out}: {n} records of {len(species)} species: "
          f"{', '.join(sorted(ncbi.name.get(s, str(s)) for s in species))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
