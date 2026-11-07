"""Eukaryote supergroups and groups by NCBI taxid (since 2026-11-05; since 2026-11-07, format 5, two levels from
data/supergroups.tsv): a species takes the row of its most specific NCBI clade in the table. supergroup: the deep
lineages (Obazoa, Amoebozoa, Archaeplastida, TSAR, Haptista, Cryptista, Discoba, Metamonada, CRuMs, ...);
group: the finer level that busco/loci.tsv's call rates, the introns' site shares, tools/rank_loci.py,
tools/eukaryote_panel.py and the benchmarks split the species by (Vertebrata, Metazoa, Fungi, Viridiplantae,
Stramenopiles, Alveolata, Rhizaria, ...). Until 2026-11-07 one level, GROUPS' eleven, SAR one of them."""

import functools
import os

TABLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "supergroups.tsv")
OTHER = "other eukaryotes"
INCERTAE = "incertae sedis"


@functools.lru_cache(maxsize=None)
def rows(path=TABLE):
    """The table -> {taxid: (NCBI name, supergroup, group)}."""
    out = {}
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    for f in lines[1:]:
        out[int(f[0])] = (f[1], f[2], f[3])
    return out


def classify(ncbi, taxid):
    """A taxon's (supergroup, group): its most specific clade's in the table; (incertae sedis, other eukaryotes) for
    a eukaryote of none of them (or not in the taxonomy)."""
    taxid = int(taxid)
    table = rows()
    for t in ncbi.ancestors(ncbi.merged.get(taxid, taxid)):
        if t in table:
            return table[t][1], table[t][2]
    return INCERTAE, OTHER


def group_of(ncbi, taxid):
    return classify(ncbi, taxid)[1]


def supergroup_of(ncbi, taxid):
    return classify(ncbi, taxid)[0]


def check(ncbi, log=print):
    """The table's taxids against the build's taxonomy: one missing or renamed is said (NCBI moves clades)."""
    bad = [f"{t} ({name}: {ncbi.name.get(t, 'not in the taxdump')})" for t, (name, _s, _g) in rows().items()
           if ncbi.name.get(t) != name]
    if bad:
        log(f"WARNING: data/supergroups.tsv names {len(bad)} clades otherwise than the taxdump: {', '.join(bad[:5])}")
    return bad
