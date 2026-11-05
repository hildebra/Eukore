"""Eukaryote supergroups by NCBI taxid (since 2026-11-05): how the package's busco/loci.tsv and tools/rank_loci.py,
tools/eukaryote_panel.py split the species when they measure how widely a BUSCO is found. A species belongs to the
first group of GROUPS whose taxid is among its ancestors, else to OTHER."""

GROUPS = [  # (name, NCBI taxid)
    ("Vertebrata", 7742), ("Metazoa", 33208), ("Fungi", 4751), ("Viridiplantae", 33090), ("Rhodophyta", 2763),
    ("SAR", 2698737), ("Discoba", 2611352), ("Metamonada", 2611341), ("Amoebozoa", 554915), ("Haptista", 2608109),
    ("Cryptophyceae", 3027)]
OTHER = "other eukaryotes"


def group_of(ncbi, taxid):
    taxid = int(taxid)
    anc = set(ncbi.ancestors(ncbi.merged.get(taxid, taxid)))
    return next((name for name, t in GROUPS if t in anc), OTHER)
