"""The BUSCO copies' checks across species (since 2026-11-05: EukDetect's per-gene QC, done on Eukore's own copies).
Each a column of busco/copies.tsv; a copy failing one is kept and flagged, so that a user chooses.

- own_species: a copy's coding sequence against the same locus' copies of every other species (vsearch, 70% at
  least, the first 32 hits): its nearest (within NEAR points of the best) agreeing with its NCBI lineage to the
  family or deeper: t; at FOREIGN_ID or more and known to part from its lineage at the class or above (apart: a phylum
  or class both name and differ in, else two groups; since 2026-11-09 two lineages NCBI names no rank of to compare,
  sister protists, are not): f (another clade's gene in the assembly: a contaminant, a misidentified genome;
  EukDetect's alien index did this against nr); else -.
  Cross-contamination (since 2026-11-09; until then both genomes holding one sequence were called f, its true owner
  too). A copy near another clade's (its `cross` partner) is told by the two clades' other members (verdict): the
  sequence is the clade's whose other species hold it nearer; where neither holds it near, the clade that holds the
  locus elsewhere and not this sequence is not its owner; where neither tells, the copy is undecided and keeps its
  partner (`cross`: the package's cross_shared flag). A genome with CROSS_MIN copies or more found another clade's of
  one partner phylum (a contaminated genome) has its undecided copies shared with that phylum called f too.
- same_as: the copy is SAME_ID identical to copies of other species: the lowest common ancestor of them all
  (EukDetect's collapse of copies 99% identical across species: a reader may give their reads to that taxon), else -.
- length: a copy's protein against its relatives' (its vsearch hits of other species, LENGTH_RELATIVES or more: their
  median, RELATIVE_RANGE), else its locus' (the median of LOCUS_MIN copies or more, LOCUS_RANGE): short or long, else
  ok. Until 2026-11-09 the bottom and top 5% of its locus (EukDetect dropped them), a tenth of every locus flagged
  whatever its lengths, a clade of longer proteins (a non-animal mitochondrion's rRNA among animals') all of them.
- nearest_id, nearest_depth: the identity of the nearest copy of another species and how deep the nearest species'
  lineages agree with the copy's (depth: 0 none .. 5 genus), None where none is 70% identical: not a column, the
  benchmarks' (benchmarks.py's resolution.tsv).
"""

import collections

from . import util

NEAR = 1.0
FOREIGN_ID = 90.0
SAME_ID = 99.0
MIN_ID = 0.7
RANKS = ("phylum", "class", "order", "family", "genus")
LENGTH_RELATIVES = 3
RELATIVE_RANGE = (0.75, 1.33)
LOCUS_MIN = 5
LOCUS_RANGE = (0.5, 2.0)
CROSS_MIN = 3


def spliced(copy, features):
    """A copy's coding sequence: its exons joined (features [(kind, start, end)] on the copy), the bases of its
    frameshifts (loci.frameshift_features, since 2026-11-06) left out."""
    shifted = {i for k, s, e in features if k == "frameshift" for i in range(int(s), int(e))}
    exons = sorted((int(s), int(e)) for k, s, e in features if k == "exon")
    if not shifted:
        return "".join(copy[s:e] for s, e in exons)
    return "".join(copy[i] for s, e in exons for i in range(s, e) if i not in shifted)


def lca(ncbi, taxids):
    """The lowest common ancestor of taxids, those not in the taxonomy left out (a reference's taxid NCBI deleted since:
    until 2026-11-06 it made every LCA it was in None)."""
    paths = [p for p in (list(ncbi.ancestors(ncbi.merged.get(t, t))) for t in taxids) if p]
    if not paths:
        return None
    common = set(paths[0]).intersection(*map(set, paths[1:]))
    return next((t for t in paths[0] if t in common), None)


def _ranks(ncbi, t, cache):
    return cache[t] if t in cache else cache.setdefault(t, ncbi.rank_ids(t))


def depth(ncbi, a, b, cache=None):
    """How deep two species' NCBI lineages agree: 0 none .. 5 genus; by taxid, to the first rank both name and differ
    at (until 2026-11-06 by name and the deepest match: two homonym genera of two kingdoms agreed to the genus). 0 says
    no more than that no named rank is shared (apart() whether they are known to part). cache: a dict for
    ncbi.rank_ids of many pairs."""
    cache = {} if cache is None else cache
    ra, rb = _ranks(ncbi, a, cache), _ranks(ncbi, b, cache)
    out = 0
    for i, r in enumerate(RANKS):
        x, y = ra.get(r), rb.get(r)
        if x is None or y is None:
            continue
        if x != y:
            break
        out = i + 1
    return out


def apart(ncbi, a, b, cache=None):
    """Whether two species' lineages are known to part at the class or above (since 2026-11-09): a phylum or class
    both name and differ in; naming none of them to compare, two groups (groups.group_of: Fungi and Metazoa,
    Stramenopiles and Alveolata), unless they share a named order, family or genus. Two sister protists whose genera
    alone NCBI ranks are not (depth gives them 0: until 2026-11-09 another clade's)."""
    cache = {} if cache is None else cache
    ra, rb = _ranks(ncbi, a, cache), _ranks(ncbi, b, cache)
    compared = False
    for r in ("phylum", "class"):
        x, y = ra.get(r), rb.get(r)
        if x is not None and y is not None:
            if x != y:
                return True
            compared = True
    if compared or any(ra.get(r) is not None and ra.get(r) == rb.get(r) for r in ("order", "family", "genus")):
        return False
    from .groups import group_of  # noqa: PLC0415
    key = ("group", a), ("group", b)
    ga = cache[key[0]] if key[0] in cache else cache.setdefault(key[0], group_of(ncbi, a))
    gb = cache[key[1]] if key[1] in cache else cache.setdefault(key[1], group_of(ncbi, b))
    return ga != gb


def clade_key(ncbi, taxid, cache=None):
    """A species' phylum taxid, else its group (cross-contamination's partners counted by it)."""
    cache = {} if cache is None else cache
    ph = _ranks(ncbi, taxid, cache).get("phylum")
    if ph is not None:
        return ph
    from .groups import group_of  # noqa: PLC0415
    return group_of(ncbi, taxid)


def families(ncbi, species, cache=None):
    """{family taxid: {species}} of species (those NCBI names a family of): who else of a species' family holds a
    locus (verdict's own and partner clades)."""
    cache = {} if cache is None else cache
    out = collections.defaultdict(set)
    for sp in species:
        fam = _ranks(ncbi, sp, cache).get("family")
        if fam is not None:
            out[fam].add(sp)
    return out


def verdict(sp, found, ncbi, cache, members, near=NEAR, foreign_id=FOREIGN_ID, related=0.0):
    """A copy's own-species verdict (module docstring), sp: its species; found: its hits of other species [(identity
    %, species)]; members: families() of the species holding the locus -> (own_species t, f or -, cross: the partner
    species where the copy is undecided between its clade and another's, or another's, else None, the top identity
    or None, the agreement depth of the nearest or None). related: the identity under which a clade's other member
    does not hold the sequence when the two clades are weighed (since 2026-11-10, the 18S' rrna.RELATED: at 80% an
    18S is any eukaryote's, and a copy 99% another family's was its own clade's where that family had no other genome
    and its own clade's 18S was 80% like it)."""
    if not found:
        return "-", None, None, None
    top = max(i for i, _s in found)
    near_ = {s for i, s in found if i >= top - near}
    agree = max(depth(ncbi, sp, s, cache) for s in near_)
    if agree >= 4:
        return "t", None, top, agree
    if top < foreign_id or not all(apart(ncbi, sp, s, cache) for s in near_):
        return "-", None, top, agree
    partner = max(found)[1]  # the nearest other clade's species

    def held(x, ranks):
        fam = ranks.get("family")
        return fam is not None and bool(members.get(fam, set()) - {x, sp, partner})
    own = max((i for i, s in found if i >= related and depth(ncbi, sp, s, cache) >= 4), default=0.0)
    theirs = max((i for i, s in found if i >= related and s not in (sp, partner)
                  and depth(ncbi, partner, s, cache) >= 4), default=0.0)
    if theirs > own:  # the partner's clade holds the sequence nearer than the copy's own: another organism's
        return "f", partner, top, agree
    if own > theirs:  # its own clade holds it nearer: the copy is the sequence's owner, the partner's a contaminant
        return "t", None, top, agree
    mine, its = held(sp, _ranks(ncbi, sp, cache)), held(partner, _ranks(ncbi, partner, cache))
    if mine and not its:  # its own clade holds the locus, not near this sequence; the partner's cannot tell
        return "f", partner, top, agree
    if its and not mine:  # the partner's clade holds the locus, not near this sequence: not theirs, the copy's
        return "-", None, top, agree
    return "-", partner, top, agree  # neither tells: undecided, shared with another clade


def check(copies, species_of, ncbi, vsearch=None, threads=1, protein_length=None, extra=None):
    """copies: {(accession, locus): coding sequence}; species_of: {accession: species taxid}; protein_length:
    {(accession, locus) or extra's key: length}; extra: {(accession, locus, k): coding sequence} checked against the
    copies but not among them (since 2026-11-09: a Duplicated locus' other copies, package.collect's choice by
    lineage) -> {key: dict(own_species, cross, same_as, length, nearest_id, nearest_depth)} of the copies and extra."""
    extra = extra or {}
    keys = list(copies) + list(extra)
    out = {k: dict(own_species="-", cross=None, same_as="-", length="ok", nearest_id=None, nearest_depth=None)
           for k in keys}
    by_locus, extra_of = collections.defaultdict(dict), collections.defaultdict(dict)
    for (acc, locus), cds in copies.items():
        if cds:
            by_locus[locus][acc] = cds
    for (acc, locus, k), cds in extra.items():
        if cds:
            extra_of[locus][f"{acc}#{k}"] = cds
    lengths = dict(protein_length or {})
    for k in keys:
        lengths.setdefault(k, len(copies.get(k) or extra.get(k) or "") // 3)

    def key_of(name, locus):
        acc, _hash, k = name.partition("#")
        return (acc, locus, int(k)) if k else (acc, locus)

    def species(name):
        return species_of[name.partition("#")[0]]
    hits_of = collections.defaultdict(dict)  # locus -> {query name: [(identity, target)]}, other species
    if vsearch:
        for locus, seqs in by_locus.items():
            if len({species_of[a] for a in seqs}) < 2:
                continue
            queries = {**seqs, **extra_of.get(locus, {})}
            for q, t, ident in util.vsearch_global(vsearch, queries, seqs if extra_of.get(locus) else queries,
                                                   MIN_ID, threads, 32, 256, ["--self"]):
                if species(q) != species(t):
                    hits_of[locus].setdefault(q, []).append((ident, t))
    for locus, seqs in by_locus.items():  # lengths: against the copy's relatives', else the locus'
        values = sorted(lengths[(a, locus)] for a in seqs)
        middle = values[len(values) // 2] if len(values) >= LOCUS_MIN else None
        for name in list(seqs) + list(extra_of.get(locus, ())):
            key, n = key_of(name, locus), lengths[key_of(name, locus)]
            relatives = sorted(lengths[(t, locus)] for _i, t in hits_of[locus].get(name, ()))
            if len(relatives) >= LENGTH_RELATIVES:
                ref, (lo, hi) = relatives[len(relatives) // 2], RELATIVE_RANGE
            elif middle:
                ref, (lo, hi) = middle, LOCUS_RANGE
            else:
                continue
            out[key]["length"] = "short" if n < lo * ref else "long" if n > hi * ref else "ok"
    if not vsearch:
        return out
    cache = {}
    for locus, seqs in by_locus.items():
        members = families(ncbi, {species_of[a] for a in seqs}, cache)
        for name, found in hits_of[locus].items():
            sp = species(name)
            got = verdict(sp, [(i, species_of[t]) for i, t in found], ncbi, cache, members)
            row = out[key_of(name, locus)]
            row["own_species"], row["cross"], row["nearest_id"], row["nearest_depth"] = got
            same = {species_of[t] for i, t in found if i >= SAME_ID}
            if same:
                row["same_as"] = str(lca(ncbi, same | {sp}))
    # a contaminated genome (CROSS_MIN copies or more another clade's, of one partner phylum): its undecided copies
    # shared with that phylum are its contaminant's too
    foreign = collections.Counter((k[0], clade_key(ncbi, r["cross"], cache)) for k, r in out.items()
                                  if len(k) == 2 and r["own_species"] == "f" and r["cross"])
    for k, r in out.items():
        if r["own_species"] == "-" and r["cross"] and foreign[(k[0], clade_key(ncbi, r["cross"], cache))] >= CROSS_MIN:
            r["own_species"] = "f"
    return out


def cross_summary(rows, ncbi):
    """A genome's copies' rows (own_species, cross) -> (copies called another clade's, the LCA of their partners or
    -): genomes.tsv's foreign_copies and foreign_source (since 2026-11-09)."""
    partners = [int(r["cross"]) for r in rows if r.get("own_species") == "f" and r.get("cross") not in (None, "-", "")]
    n = sum(1 for r in rows if r.get("own_species") == "f")
    return n, (lca(ncbi, set(partners)) or "-") if partners else "-"
