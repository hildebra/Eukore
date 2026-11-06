"""The BUSCO copies' checks across species (since 2026-11-05: EukDetect's per-gene QC, done on Eukore's own copies).
Each a column of busco/copies.tsv; a copy failing one is kept and flagged, so that a user chooses.

- own_species: a copy's coding sequence against the same locus' copies of every other species (vsearch, 70% at
  least, the first 32 hits): its nearest (within NEAR points of the best) agreeing with its NCBI lineage to the
  family or deeper: t; sharing no more than the phylum at FOREIGN_ID or more: f (another clade's gene in the assembly:
  a contaminant, a misidentified genome; EukDetect's alien index did this against nr); else -.
- same_as: the copy is SAME_ID identical to copies of other species: the lowest common ancestor of them all
  (EukDetect's collapse of copies 99% identical across species: a reader may give their reads to that taxon), else -.
- length: a copy whose protein is in the bottom or top LENGTH_TAIL of its locus' copies: short or long (EukDetect
  dropped them), else ok.
- nearest_id, nearest_depth: the identity of the nearest copy of another species and how deep the nearest species'
  lineages agree with the copy's (depth: 0 none .. 5 genus), None where none is 70% identical: not a column, the
  benchmarks' (benchmarks.py's resolution.tsv).
"""

import collections
import os
import subprocess
import tempfile

NEAR = 1.0
FOREIGN_ID = 90.0
SAME_ID = 99.0
MIN_ID = 0.7
LENGTH_TAIL = 0.05
RANKS = ("phylum", "class", "order", "family", "genus")


def spliced(copy, features):
    """A copy's coding sequence: its exons joined (features [(kind, start, end)] on the copy), the bases of its
    frameshifts (loci.frameshift_features, since 2026-11-06) left out."""
    shifted = {i for k, s, e in features if k == "frameshift" for i in range(int(s), int(e))}
    exons = sorted((int(s), int(e)) for k, s, e in features if k == "exon")
    if not shifted:
        return "".join(copy[s:e] for s, e in exons)
    return "".join(copy[i] for s, e in exons for i in range(s, e) if i not in shifted)


def lca(ncbi, taxids):
    """The lowest common ancestor of taxids."""
    paths = [list(ncbi.ancestors(ncbi.merged.get(t, t))) for t in taxids]
    if not paths:
        return None
    common = set(paths[0]).intersection(*map(set, paths[1:]))
    return next((t for t in paths[0] if t in common), None)


def depth(ncbi, a, b):
    """How deep two species' NCBI lineages agree: 0 none, 1 phylum ... 5 genus."""
    ra, rb = ncbi.ranks(a), ncbi.ranks(b)
    return max([0] + [i + 1 for i, r in enumerate(RANKS) if ra.get(r) and ra.get(r) == rb.get(r)])


def check(copies, species_of, ncbi, vsearch=None, threads=1, protein_length=None):
    """copies: {(accession, locus): coding sequence}; species_of: {accession: species taxid}; protein_length:
    {(accession, locus): length} -> {(accession, locus): dict(own_species, same_as, length)}."""
    out = {k: dict(own_species="-", same_as="-", length="ok", nearest_id=None, nearest_depth=None) for k in copies}
    by_locus = collections.defaultdict(dict)
    for (acc, locus), cds in copies.items():
        if cds:
            by_locus[locus][acc] = cds
    lengths = protein_length or {k: len(v) // 3 for k, v in copies.items()}
    for locus in by_locus:
        values = sorted(lengths[(a, locus)] for a in by_locus[locus])
        if len(values) >= 20:
            lo, hi = values[int(LENGTH_TAIL * len(values))], values[int((1 - LENGTH_TAIL) * len(values)) - 1]
            for a in by_locus[locus]:
                n = lengths[(a, locus)]
                out[(a, locus)]["length"] = "short" if n < lo else "long" if n > hi else "ok"
    if not vsearch:
        return out
    for locus, seqs in by_locus.items():
        if len({species_of[a] for a in seqs}) < 2:
            continue
        hits = collections.defaultdict(list)
        with tempfile.TemporaryDirectory() as tmp:
            fasta, table = os.path.join(tmp, "c.fna"), os.path.join(tmp, "hits.tsv")
            with open(fasta, "w") as fh:
                fh.writelines(f">{a}\n{s}\n" for a, s in seqs.items())
            subprocess.run([vsearch, "--usearch_global", fasta, "--db", fasta, "--id", str(MIN_ID), "--iddef", "2",
                            "--maxaccepts", "32", "--maxrejects", "256", "--self", "--userout", table, "--userfields",
                            "query+target+id", "--threads", str(max(1, threads)), "--quiet"], check=True)
            with open(table) as fh:
                for line in fh:
                    q, t, ident = line.rstrip("\n").split("\t")
                    if species_of[q] != species_of[t]:
                        hits[q].append((float(ident), t))
        for acc in seqs:
            found = hits.get(acc, [])
            if not found:
                continue
            top = max(i for i, _t in found)
            near = [t for i, t in found if i >= top - NEAR]
            agree = max(depth(ncbi, species_of[acc], species_of[t]) for t in near)
            row = out[(acc, locus)]
            row["nearest_id"], row["nearest_depth"] = top, agree  # benchmarks.py's resolution
            row["own_species"] = "t" if agree >= 4 else "f" if agree <= 1 and top >= FOREIGN_ID else "-"
            same = {species_of[t] for i, t in found if i >= SAME_ID}
            if same:
                row["same_as"] = str(lca(ncbi, same | {species_of[acc]}))
    return out
