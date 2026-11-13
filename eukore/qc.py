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

import argparse
import collections
import glob
import gzip
import hashlib
import os
import re
import shutil
import sys
import tempfile
import time

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
SEARCHED = "loci searched (vsearch)"
SEARCH_LIMITS = (32, 256)       # vsearch's --maxaccepts and --maxrejects of a locus' search
SEARCH_OPTIONS = ("--self",)
SEARCH_RECORD = 1               # what makes a kept search (searched_loci): raised when its options or files change
SEARCH_POLL = 30                # seconds between looks at the loci other nodes search


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
    and its own clade's 18S was 80% like it). The 18S' check has a rule of its own since that day (rrna.verdict: a
    lone lineage's 18S is not another clade's, of two species at the top identity the partner is not the higher
    taxid): this one is the coding sequences', as it was."""
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


def check(copies, species_of, ncbi, vsearch=None, threads=1, protein_length=None, extra=None, progress=None,
          searches=None):
    """copies: {(accession, locus): coding sequence}; species_of: {accession: species taxid}; protein_length:
    {(accession, locus) or extra's key: length}; extra: {(accession, locus, k): coding sequence} checked against the
    copies but not among them (since 2026-11-09: a Duplicated locus' other copies, package.collect's choice by
    lineage) -> {key: dict(own_species, cross, same_as, length, nearest_id, nearest_depth)} of the copies and extra.
    progress: a function (what, weights) -> util.Progress, the loci's vsearch searches reported by it (since
    2026-11-11: hours of them said nothing), each locus weighing its copies squared. searches: dict(folder, helpers)
    (searched_loci; since 2026-11-11): each locus' search kept in the folder, shared with other nodes."""
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
        todo = [(locus, {**seqs, **extra_of.get(locus, {})}, seqs if extra_of.get(locus) else None)
                for locus, seqs in by_locus.items() if len({species_of[a] for a in seqs}) >= 2]
        if searches:  # each locus' search kept, and shared with other nodes (since 2026-11-11)
            found = searched_loci(todo, vsearch, threads, searches["folder"], searches.get("helpers"), progress,
                                  scratch=searches.get("scratch"))
        else:
            found = {}
            shown = progress(SEARCHED, [weight(q) for _l, q, _d in todo]) if progress and todo else None
            for locus, queries, db in todo:
                found[locus] = util.vsearch_global(vsearch, queries, queries if db is None else db, MIN_ID, threads,
                                                   *SEARCH_LIMITS, list(SEARCH_OPTIONS))
                if shown:
                    shown.step(weight(queries))
        for locus, hits in found.items():
            for q, t, ident in hits:
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


# --- the loci's searches, kept and shared by several nodes (since 2026-11-11) --------------------------------------
# A locus' search of all its copies (some 30,000) took the package stage's node for hours, one locus after another.
# With searches (check's), each search is kept in a folder of the WORK (WORK/qc/<kind>/), named by the locus and a
# digest of its sequences and the search's options: a rerun with the same copies reads it (the package stage again,
# after a failure later on), a locus whose copies changed is searched again. helpers, given, are other nodes (a SLURM
# array, tools/eukore.sh's qc): the loci's input files are written into the folder first, the helpers started, and
# every node then claims a locus at a time (genomes.claimed_file: a killed node's locus taken over), the largest first.
# The node that wrote them waits for the others' last loci, then reads them all. Files of a folder:
#   <locus>.<digest>.q.fna, .db.fna   a locus' queries and (where not the queries) its database, gone once searched
#   <locus>.<digest>.hits.gz          its hits: query, target, identity %
#   .<locus>.<digest>.claim           the node searching it
#   db.<digest>.fna                   a database several searches share (Shared; since 2026-11-13), gone at the end
#   todo.tsv                          the loci to search, the largest first (what the helpers read): a line a search,
#                                     its name and weight, and (since 2026-11-13) its identity, --maxaccepts,
#                                     --maxrejects, further options (- for none) and shared database file (- for none)
# Since 2026-11-13 the rRNA checks' searches go the same way (screen.Searcher): each its own settings (Params), a large
# one split by its queries into searches against one shared database (a query's hits are the same whoever else is
# searched with it).


class Params(tuple):
    """A search's settings: (identity, --maxaccepts, --maxrejects, further options); DEFAULT the loci's."""
    __slots__ = ()

    def __new__(cls, ident=MIN_ID, accepts=SEARCH_LIMITS[0], rejects=SEARCH_LIMITS[1], options=SEARCH_OPTIONS):
        return super().__new__(cls, (ident, accepts, rejects, tuple(options)))

    def text(self):
        """What search_key hashes: for DEFAULT the text of before 2026-11-13 (the loci's kept searches stay)."""
        ident, accepts, rejects, options = self
        return f"{SEARCH_RECORD}\t{ident}\t{(accepts, rejects)}\t{options}\n"

    def columns(self):
        ident, accepts, rejects, options = self
        return [str(ident), str(accepts), str(rejects), " ".join(options) or "-"]

    @classmethod
    def of_columns(cls, fields):
        """todo.tsv's columns after a search's weight (none: DEFAULT)."""
        if len(fields) < 4:
            return DEFAULT
        return cls(float(fields[0]), int(fields[1]), int(fields[2]), () if fields[3] == "-" else fields[3].split())


DEFAULT = Params()
CM_RECORD = 1  # what makes a kept cmsearch part (CmParams): raised when util.CMSEARCH_OPTIONS or its table change


class CmParams(tuple):
    """A cmsearch part's settings (since 2026-11-13: the SSU models' search, regions.ssu_models): the models' file (on
    the file system the helpers share) and its digest, the options util.CMSEARCH_OPTIONS'. Its hits: (sequence, model
    accession, bits) (util.cmsearch_hits)."""
    __slots__ = ()

    def __new__(cls, cm, digest):
        return super().__new__(cls, (cm, digest))

    def text(self):
        return f"cmsearch\t{CM_RECORD}\t{self[1]}\t{util.CMSEARCH_OPTIONS}\n"

    def columns(self):
        return ["cmsearch", self[0], self[1], "-"]


class Shared:
    """A database several searches of a folder share ({name: sequence}), written once as db.<digest>.fna."""
    def __init__(self, seqs):
        self.seqs = seqs
        h = hashlib.sha1()
        for name in sorted(seqs):
            h.update(f"{name}\t{seqs[name]}\n".encode())
        self.digest = h.hexdigest()[:20]
        self.file = f"db.{self.digest}.fna"


def weight(queries):
    """A locus' share of the search's work: its sequences squared."""
    return len(queries) ** 2


def search_key(queries, db=None, params=DEFAULT):
    """The digest of a locus' search: SEARCH_RECORD, the options, its queries and database (None: the queries; a
    Shared one by its digest)."""
    h = hashlib.sha1(params.text().encode())
    for part in (queries,) + (() if db is None else (db,)):
        h.update(b">\n")
        if isinstance(part, Shared):
            h.update(f"shared {part.digest}\n".encode())
            continue
        for name in sorted(part):
            h.update(f"{name}\t{part[name]}\n".encode())
    return h.hexdigest()[:20]


def search_base(locus, queries, db=None, params=DEFAULT):
    """A locus' search's name in its folder: the locus (its characters a file name takes) and search_key."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", str(locus)) + "." + search_key(queries, db, params)


def search_files(folder, base):
    """(queries, database, hits, claim) of a locus' search named base."""
    p = os.path.join(folder, base)
    return p + ".q.fna", p + ".db.fna", p + ".hits.gz", os.path.join(folder, f".{base}.claim")


def _write_fasta(path, seqs):
    if os.path.isfile(path):
        return
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in seqs.items() if s)
    os.replace(util.part_of(path), path)


def write_inputs(folder, base, queries, db=None):
    """A locus' input files, each written whole then renamed (a helper reads none half written), the database
    first (a Shared one once for all its searches): a helper takes a locus by its queries' file."""
    q, d, _hits, _claim = search_files(folder, base)
    if isinstance(db, Shared):
        _write_fasta(os.path.join(folder, db.file), db.seqs)
    elif db is not None:
        _write_fasta(d, db)
    _write_fasta(q, queries)


def search_one(folder, base, vsearch, threads, inputs=None, params=DEFAULT, shared=None, inputs_dir=None,
               timing=None):
    """A locus' search, where no other node holds it (its claim; a dead node's taken over): its input files (written
    from inputs (queries, db) where not there; a helper, given none, leaves a locus without them), vsearch with
    params (Params), the hits written whole and the inputs removed (a shared database, its file name `shared`, kept
    for the folder's other searches) -> whether its hits are there now (False: another node searches it, or no
    inputs). inputs_dir: a folder for the inputs and vsearch's output other than the kept folder (a node-local one,
    where no helper reads them). timing: a dict given vsearch's times, "setup" (the database read and indexed, on
    one thread; None where vsearch did not say) and "all" (util.vsearch_files' timed; since 2026-11-13). With
    CmParams, vsearch is cmsearch, the queries the sequences searched with its models (no database)."""
    from . import genomes  # noqa: PLC0415
    q, d, hits, claim = search_files(folder, base)
    if os.path.isfile(hits):
        return True
    with genomes.claimed_file(claim) as held:
        if held is None or os.path.isfile(hits):  # another node's, or searched by one before this claim
            return os.path.isfile(hits)
        where = inputs_dir or folder
        if inputs_dir:  # the inputs and vsearch's own output on the node's disk (since 2026-11-13), the hits kept
            q, d = search_files(inputs_dir, base)[:2]
        if inputs is not None:
            write_inputs(where, base, *inputs)
        db = os.path.join(where, shared) if shared else d if os.path.isfile(d) else q
        if not os.path.isfile(q) or not os.path.isfile(db):
            return False
        out = util.part_of(os.path.join(where, os.path.basename(hits)[:-len(".gz")]))
        try:
            if isinstance(params, CmParams):
                began = time.time()
                util.cmsearch_files(vsearch, q, params[0], out, threads)
                took = (None, time.time() - began)
            else:
                ident, accepts, rejects, options = params
                took = util.vsearch_files(vsearch, q, db, out, ident, threads, accepts, rejects, options,
                                          timed=timing is not None)
            if timing is not None:
                timing["setup"], timing["all"] = took
            with open(out, "rb") as src, gzip.open(util.part_of(hits), "wb", compresslevel=3) as dst:
                shutil.copyfileobj(src, dst)
            os.replace(util.part_of(hits), hits)
        finally:
            for left in (out, util.part_of(hits)):
                if os.path.exists(left):
                    os.remove(left)
        for path in (q, d):
            if os.path.exists(path):
                os.remove(path)
    return True


def read_hits(folder, base, params=DEFAULT):
    with gzip.open(search_files(folder, base)[2], "rt") as fh:
        return util.cmsearch_hits(fh) if isinstance(params, CmParams) else util.vsearch_hits(fh)


def todo_line(base, size, params=DEFAULT, shared=None):
    """A search's line of todo.tsv: its name and weight, and where not the loci's DEFAULT without a shared database
    (since 2026-11-13) its settings and shared database file (-: none)."""
    fields = [base, str(size)]
    if params != DEFAULT or shared:
        fields += params.columns() + [shared or "-"]
    return "\t".join(fields) + "\n"


def todo_entry(line):
    """A line of todo.tsv -> (name, params, shared database file or None); a line of before 2026-11-13: DEFAULT."""
    f = line.rstrip("\n").split("\t")
    if len(f) >= 7 and f[2] == "cmsearch":  # a cmsearch part (since 2026-11-13): the models' file and digest
        return f[0], CmParams(f[3], f[4]), None
    if len(f) >= 7:
        return f[0], Params.of_columns(f[2:6]), None if f[6] == "-" else f[6]
    return f[0], DEFAULT, None


def searched_loci(todo, vsearch, threads, folder, helpers=None, progress=None, poll=SEARCH_POLL, params=DEFAULT,
                  what=None, scratch=None, log=None):
    """The loci's searches todo [(locus, queries, db or None or a Shared one)], kept in folder: those kept read, the
    others searched here and, with helpers (a function of the number of loci to search -> a function that stops the
    helpers, or None), by other nodes too -> {locus: [(query, target, identity %)]}. params: the searches' settings
    (Params; since 2026-11-13, the loci's DEFAULT), written into todo.tsv for the helpers; what: the progress' name
    (default SEARCHED); scratch: a node-local folder for the inputs where no helper searches with this node (since
    2026-11-13); log: each locus searched here said with its times (searched_line; since 2026-11-13). The folder's
    files of other searches (an older set of copies) are removed at the end."""
    os.makedirs(folder, exist_ok=True)
    base = {locus: search_base(locus, q, d, params) for locus, q, d in todo}
    inputs = {locus: (q, d) for locus, q, d in todo}
    shared = {locus: d.file if isinstance(d, Shared) else None for locus, _q, d in todo}
    order = sorted(todo, key=lambda t: -weight(t[1]))
    shown = progress(what or SEARCHED, [weight(q) for _l, q, _d in order]) if progress and order else None
    left = []
    for locus, q, _d in order:
        if os.path.isfile(search_files(folder, base[locus])[2]):
            if shown:
                shown.step(weight(q), skipped=True)  # kept from an earlier run: out of the time left
        else:
            left.append(locus)
    stop = None
    if left and helpers:
        for locus in left:
            write_inputs(folder, base[locus], *inputs[locus])
        listing = os.path.join(folder, "todo.tsv")
        with open(util.part_of(listing), "w", newline="\n") as fh:
            fh.writelines(todo_line(base[locus], weight(inputs[locus][0]), params, shared[locus]) for locus in left)
        os.replace(util.part_of(listing), listing)
        stop = helpers(len(left))
    # searched here alone: the inputs on the node's disk (scratch: --local_tmp; since 2026-11-13), the WORK's file
    # system given the hits alone
    local = tempfile.mkdtemp(dir=scratch, prefix="search.") if left and not stop and scratch else None
    try:
        while left:
            for locus in list(left):
                timing = {} if log else None
                if search_one(folder, base[locus], vsearch, threads, inputs[locus], params, shared[locus], local,
                              timing):
                    left.remove(locus)
                    if timing:
                        log(f"  {searched_line(base[locus], timing)}")
                    if shown:
                        shown.step(weight(inputs[locus][0]))
            if left:  # the others' loci: their hits awaited (a dead node's claim taken over by search_one)
                if shown:
                    shown.tick(f"; {len(left):,} searched by other nodes")
                time.sleep(poll)
    finally:
        if stop:
            stop()
        if local:
            shutil.rmtree(local, ignore_errors=True)
    found = {locus: read_hits(folder, base[locus], params) for locus, _q, _d in todo}
    keep = set(base.values())
    for path in glob.glob(os.path.join(folder, "*")) + glob.glob(os.path.join(folder, ".*.claim")):
        name = os.path.basename(path).lstrip(".")
        if name.endswith(".partial") or os.path.isdir(path):  # a helper's file being written: its own
            continue
        if name == "todo.tsv" or not any(name.startswith(b + ".") for b in keep):
            os.remove(path)
    return found


def searched_line(base, timing):
    """'screen_18S-3.<key>: searched in 2m 10s (database read and indexed in 25s)': a part's times (search_one's
    timing)."""
    setup = timing.get("setup")
    return (f"{base}: searched in {util.duration(timing['all'])}"
            + ("" if setup is None else f" (database read and indexed in {util.duration(setup)})"))


def help_search(root, vsearch, threads, log=print, cmsearch=None):
    """A helper node's part (tools/eukore.sh qc): the loci of every folder under root (its todo.tsv) whose inputs are
    there and hits not yet, each claimed and searched, until none is left that no other node holds -> loci searched.
    A cmsearch part (CmParams; since 2026-11-13) by cmsearch, left to the other nodes where it is None."""
    done, idle = 0, False
    while not idle:
        idle = True
        for listing in sorted(glob.glob(os.path.join(root, "*", "todo.tsv"))):
            folder = os.path.dirname(listing)
            try:
                with open(listing) as fh:
                    entries = [todo_entry(line) for line in fh if line.strip()]
            except OSError:
                continue
            for base, params, shared in entries:
                q, _d, hits, _claim = search_files(folder, base)
                tool = cmsearch if isinstance(params, CmParams) else vsearch
                if os.path.isfile(hits) or not os.path.isfile(q) or not tool:
                    continue
                timing = {}
                if search_one(folder, base, tool, threads, params=params, shared=shared, timing=timing):
                    done += 1
                    idle = False
                    if timing:
                        log(searched_line(base, timing))
    log(f"{done:,} loci searched here; none left")
    return done


def main(argv=None):
    p = argparse.ArgumentParser(description="A helper node's searches of the package stage's loci (help_search).")
    p.add_argument("command", choices=["helper"])
    p.add_argument("root", help="WORK/qc: the folders of the loci's searches")
    p.add_argument("--vsearch", default="vsearch")
    p.add_argument("--cmsearch", default="cmsearch", help="the SSU models' parts' (none found: left to other nodes)")
    p.add_argument("-t", "--threads", type=int, default=1)
    a = p.parse_args(argv)
    vsearch = shutil.which(a.vsearch)
    if not vsearch:
        sys.exit(f"no {a.vsearch}")
    help_search(a.root, vsearch, a.threads,
                log=lambda text: print(f"[{time.strftime('%H:%M:%S')}] {text}", flush=True),
                cmsearch=shutil.which(a.cmsearch))
    return 0


if __name__ == "__main__":
    sys.exit(main())
