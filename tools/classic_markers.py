#!/usr/bin/env python3
"""Which of BUSCO's eukaryota loci are the classic phylogenetic and barcode markers (since 2026-11-07): RPB1, RPB2,
EF-1 alpha (TEF1), beta-tubulin (TUB2), actin, HSP90, histone H3, GAPDH, EF-2, calmodulin. Their proteins of yeast and
human (UniProt, MARKERS) are searched with every BUSCO HMM (hmmsearch); a marker is the BUSCO whose HMM scores its
proteins at BUSCO's cutoff or more, best. Writes data/classic_markers.faa (the proteins, fetched once) and
data/classic_markers.tsv (marker, busco, its OrthoDB name, the proteins' best score, the cutoff; `-` for a marker no
BUSCO is), which busco/loci.tsv's `classic_name` reads.

    python3 tools/classic_markers.py --odb WORK/sources/eukaryota_odb12
"""

import argparse
import os
import subprocess
import sys
import tempfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import odb as odb_mod  # noqa: E402

DATA = os.path.join(os.path.dirname(HERE), "data")
FAA, TABLE = os.path.join(DATA, "classic_markers.faa"), os.path.join(DATA, "classic_markers.tsv")
# marker: UniProt accessions (Saccharomyces cerevisiae, Homo sapiens)
MARKERS = {"RPB1": ("P04050", "P24928"), "RPB2": ("P08518", "P30876"), "TEF1": ("P02994", "P68104"),
           "TUB2": ("P02557", "P07437"), "ACT1": ("P60010", "P60709"), "HSP90": ("P02829", "P07900"),
           "H3": ("P61830", "P68431"), "GAPDH": ("P00359", "P04406"), "EF2": ("P32324", "P13639"),
           "CMD1": ("P06787", "P0DP23")}
UNIPROT = "https://rest.uniprot.org/uniprotkb/{}.fasta"


def fetch(path):
    """The markers' proteins from UniProt into path (>marker|accession), kept when there."""
    if os.path.isfile(path):
        return
    out = []
    for marker, accessions in MARKERS.items():
        for acc in accessions:
            with urllib.request.urlopen(UNIPROT.format(acc), timeout=60) as r:
                lines = r.read().decode().splitlines()
            out.append(f">{marker}|{acc}\n" + "".join(lines[1:]) + "\n")
    with open(path, "w", newline="\n") as fh:
        fh.writelines(out)


def best_loci(faa, odb, hmmsearch="hmmsearch"):
    """{marker: (busco, score)} of the BUSCO whose HMM scores the marker's proteins best, at its cutoff or more."""
    with tempfile.TemporaryDirectory() as tmp:
        db, tbl = odb.hmm_db(os.path.join(tmp, "all.hmm")), os.path.join(tmp, "hits.tbl")
        subprocess.run([hmmsearch, "--noali", "--tblout", tbl, "-o", os.devnull, db, faa], check=True)
        best = {}
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                marker, locus, score = f[0].split("|")[0], odb.locus_of_hmm.get(f[2], f[2]), float(f[5])
                if locus in odb.cutoff and score >= odb.cutoff[locus] and score > best.get(marker, ("", 0.0))[1]:
                    best[marker] = (locus, score)
    return best


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--odb", required=True, help="BUSCO's eukaryota_odb12 folder (WORK/sources/eukaryota_odb12)")
    p.add_argument("--hmmsearch", default="hmmsearch")
    opts = p.parse_args(argv)
    fetch(FAA)
    odb = odb_mod.Odb(opts.odb, odb_mod.all_loci(opts.odb))
    best = best_loci(FAA, odb, opts.hmmsearch)
    with open(TABLE, "w", newline="\n") as fh:
        fh.write("# BUSCO's eukaryota loci that are classic markers: tools/classic_markers.py "
                 f"({os.path.basename(opts.odb.rstrip('/'))})\n")
        fh.write("marker\tbusco\tname\tscore\tcutoff\n")
        for marker in MARKERS:
            locus, score = best.get(marker, ("-", None))
            fh.write(f"{marker}\t{locus}\t{odb.names.get(locus, '-')}\t"
                     f"{'-' if score is None else f'{score:.1f}'}\t{odb.cutoff.get(locus, '-')}\n")
    print(open(TABLE).read())


if __name__ == "__main__":
    main()
