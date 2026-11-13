"""The rRNA part screen (since 2026-11-05; protal's rdna_screen.py on NCBI's lineages): a record holding a part of
another organism (a contaminant contig's operon, a misassembled mosaic, a chimeric read) is flagged.

Each record is cut into its parts: its 18S, its junction (ITS1-5.8S-ITS2 with FLANK bases of the 18S and the 28S) and
the two halves of its 28S; a record without regions is its 18S. Each part is placed among the same part of the other
records (vsearch, SCREEN_ID or more, the first SCREEN_ACCEPTS). A part is another organism's where
- its nearest records (within NEAR points of the best; records of its own genome count on neither side) name no taxon
  of its own clade at the rank tested: the deepest of kingdom, phylum, class, order and family at which its NCBI
  lineage names a taxon and the records hold that part of another genome of it;
- it is OTHER_ID or more identical to the best of them, and MARGIN points or more less identical to every record of
  its own clade (searched again among its clade alone; SISTER_MARGIN where they share its order: rRNA often cannot tell
  a family from its sister);
- and the matched record lies nearer its own clade than the part does (a contaminant among the records would
  otherwise make its source's own records look foreign); a matched record whose clade has no other genome vouches for
  nothing.
The record's `screen` column names the parts and ranks (`28S 3':family`); `-` for none. Eukore flags, the user drops.
"""

import collections
import re
import tempfile
import time

from . import util

FLANK = 50
MIN_PART = 300
SCREEN_ID, SCREEN_ACCEPTS, SCREEN_REJECTS = 0.75, 32, 64
NEAR = 1.0
RANKS = ("kingdom", "phylum", "class", "order", "family")
OTHER_ID = 99.0
MARGIN, SISTER_MARGIN = 5.0, 10.0
PARTS = ("18S", "junction", "28S 5'", "28S 3'")


def parts_of(seq, regions):
    """A record's parts {name: sequence} (those of MIN_PART bases or more); regions [(name, start, end)] or None."""
    if not regions:
        return {"18S": seq} if len(seq) >= MIN_PART else {}
    ssu = next((r for r in regions if r[0] == "18S"), None)
    lsu = next((r for r in regions if r[0] == "28S"), None)
    out = {}
    if ssu:
        out["18S"] = seq[ssu[1]:ssu[2]]
    if ssu and lsu:
        out["junction"] = seq[max(ssu[1], ssu[2] - FLANK):min(lsu[2], lsu[1] + FLANK)]
    if lsu:
        mid = (lsu[1] + lsu[2]) // 2
        out["28S 5'"], out["28S 3'"] = seq[lsu[1]:mid], seq[mid:lsu[2]]
    return {k: v for k, v in out.items() if len(v) >= MIN_PART}


def lineage_of(ncbi, taxid):
    """A taxid's names at RANKS (None where NCBI names none there)."""
    at = ncbi.ranks(int(taxid))
    kingdom = at.get("kingdom") or "Eukaryota"  # protists: no NCBI kingdom
    return [kingdom] + [at.get(r) for r in RANKS[1:]]


def clade(lineage, rank):
    """A lineage's clade at rank (0 the kingdom .. 4 the family): the names it holds down to there, or None where it
    names none at that rank (until 2026-11-06 also where a rank above was empty: a protist without an NCBI phylum was
    tested at its kingdom alone)."""
    if not lineage[rank]:
        return None
    return tuple(x for x in lineage[:rank + 1] if x)


def agrees_to(a, b):
    """The deepest rank at which two lineages name the same clade; -1 where none (ranks a names none at skipped)."""
    depth = -1
    for rank in range(len(RANKS)):
        k = clade(a, rank)
        if k is None:
            continue
        if k != clade(b, rank):
            break
        depth = rank
    return depth


def search(queries, db_seqs, vsearch, threads, tmp, ident=SCREEN_ID, accepts=SCREEN_ACCEPTS, rejects=SCREEN_REJECTS):
    """vsearch --usearch_global (util.vsearch_global; the same set as queries and db written once) -> {query:
    [(identity, target)]}, self hits left out."""
    if not queries or not db_seqs:
        return {}
    hits = collections.defaultdict(list)
    for a, b, found in util.vsearch_global(vsearch, queries, db_seqs, ident, threads, accepts, rejects, tmp=tmp):
        if a != b:
            hits[a].append((found, b))
    return hits


PER_NODE = 3  # a kept search's parts for each node searching it (the final job's and its helpers'; since 2026-11-13)
SEARCH_MIN_PART = 5000  # queries a part has at least: each part reads and indexes the whole database again (1 thread)
HELPED = 4  # parts a search has at least for the helper nodes to be started for it


def search_parts(names, nodes=1):
    """names (sorted) -> the parts of a kept search, [[name]]: PER_NODE for each of `nodes` nodes searching, fewer
    where a part would have less than SEARCH_MIN_PART, interleaved (part i every n-th name from the i-th) so that each
    holds some of every clade (since 2026-11-13: until then consecutive runs of 10,000, at most 32 parts: a part of one
    dense clade's sequences, many hits each, took 5 times as long as others, the search waiting for it)."""
    n = min(max(1, PER_NODE * max(1, nodes)), max(1, len(names) // SEARCH_MIN_PART))
    return [names[i::n] for i in range(n)] if names else []


class Searcher:
    """The rRNA checks' vsearch searches (own_species', screen's; since 2026-11-13; until then each was run here, one
    after another, and again in full by every package stage: one node for hours).

    one(): a large search. With a folder (WORK/qc/rrna) it is kept in <folder>_<phase> (qc.searched_loci), its
    queries in parts (search_parts: PER_NODE for each of `nodes`, this one and its helpers) against one shared database,
    each part's hits those of the whole search (a query's hits do not depend on the other queries), searched by the
    helper nodes too where it has HELPED parts or more: a rerun on the same records (and nodes) reads them. many(): the
    small searches that follow up on it (a clade's records each), several at once here, each with a share of the
    threads (until then one after another, each with them all, the parent writing and reading their files on one
    core). log: each search's time, and each part's searched here."""

    def __init__(self, vsearch, threads=1, tmp=None, folder=None, helpers=None, progress=None, log=None, nodes=1):
        self.vsearch, self.threads, self.tmp, self.folder = vsearch, max(1, threads), tmp, folder
        self.helpers, self.progress, self.log, self.nodes = helpers, progress, log, max(1, nodes)

    def _say(self, text):
        if self.log:
            self.log(text)

    def one(self, phase, queries, db, ident, accepts, rejects):
        """queries {name: sequence} against db (queries itself: the same object) -> {query: [(identity, target)]},
        self hits left out (search's)."""
        if not queries or not db:
            return {}
        from . import qc  # noqa: PLC0415
        began = time.time()
        phase = re.sub(r"[^A-Za-z0-9_-]", "_", phase)
        if not self.folder:
            got = search(queries, db, self.vsearch, self.threads, self.tmp, ident, accepts, rejects)
        else:
            names = sorted(n for n, s in queries.items() if s)
            shared = qc.Shared({n: s for n, s in db.items() if s})
            todo = [(f"{phase}-{i + 1}", {n: queries[n] for n in part}, shared)
                    for i, part in enumerate(search_parts(names, self.nodes))]
            found = qc.searched_loci(todo, self.vsearch, self.threads, f"{self.folder}_{phase}",
                                     self.helpers if len(todo) >= HELPED else None, self.progress,
                                     params=qc.Params(ident, accepts, rejects, ()), what=f"{phase} searched (vsearch)",
                                     scratch=self.tmp, log=self.log)
            got = collections.defaultdict(list)
            for hits in found.values():
                for a, b, i in hits:
                    if a != b:
                        got[a].append((i, b))
        self._say(f"  rRNA {phase}: {len(queries):,} against {len(db):,} (vsearch at {ident:g}): took "
                  f"{util.duration(time.time() - began)}")
        return got

    def many(self, what, items, ident, accepts, rejects):
        """[(queries, db)] -> [search's {query: [(identity, target)]}], in their order: up to `threads` at once,
        each with its share of the threads."""
        if not items:
            return []
        import concurrent.futures  # noqa: PLC0415
        began = time.time()
        workers = min(len(items), self.threads)
        each = max(1, self.threads // workers)

        def one(item):
            return search(item[0], item[1], self.vsearch, each, self.tmp, ident, accepts, rejects)
        with concurrent.futures.ThreadPoolExecutor(workers) as pool:
            out = list(pool.map(one, items))
        self._say(f"  rRNA {what}: {len(items):,} searches, {workers} at once: took "
                  f"{util.duration(time.time() - began)}")
        return out


def verdict_of(depth, top, own, theirs):
    """A suspect part's verdict: "flagged" or "kept: why" (protal's rdna_screen.verdict_of)."""
    margin = SISTER_MARGIN if depth >= 3 else MARGIN
    if theirs == 0.0:
        return "kept: the matched record's clade has no other genome's record that close"
    if top < OTHER_ID:
        return f"kept: under {OTHER_ID:g}% to the other clade"
    if top - max(own, 100 * SCREEN_ID) < margin:
        return f"kept: under {margin:g} points nearer the other clade"
    if own >= theirs:
        return "kept: the matched record lies farther from its own clade"
    return "flagged"


ASSEMBLY = re.compile(r"^GC[AF]_(\d+)(?:\.\d+)?$")


def genome_of(name, row):
    """A record's genome: its assembly, one name for its GenBank and RefSeq copies and every version of it (since
    2026-11-09: ROD's operons are NCBI assemblies' own, and its copy of the genome a record is extracted from, another
    version's or RefSeq's, confirmed the record as another genome's would, a contaminant's operon too), else its
    accession, else its name."""
    acc = row.get("accession") or "-"
    if acc == "-":
        return name
    m = ASSEMBLY.match(acc)
    return f"GC_{m.group(1)}" if m else acc


def representatives(seqs, owner, threads=1, log=None):
    """Each owner's sequences ({name: sequence}; owner: {name: its genome}) in 99% types (rrna.types_of) -> {the type's
    representative: [the names of its type]}, the representative the name of its longest sequence, of several the
    first (since 2026-11-10; until then the type's first name whatever its length: a copy cut at a contig's end, of
    one type with the whole copies as rrna.identity measures over the shorter, was searched in their place and its
    verdict theirs). A genome's near-identical copies (an rDNA array: ~400 in a T2T human, its 18S once as an 18S
    record and once in its rDNA record) are searched once: until 2026-11-06 they filled vsearch's first hits, which
    were then dropped as the genome's own, and the screen and the own-species check found nothing. The owners
    typed by a pool of `threads` processes (rrna.species_types; since 2026-11-13: one after another in this process,
    over an hour for the 4.6 million records of all eukaryotes, again for each part of the screen); log: an owner whose
    typing took a minute or more."""
    from .rrna import species_types  # noqa: PLC0415 (rrna imports this module's users)
    by_owner = collections.defaultdict(list)
    for name in sorted(seqs):
        by_owner[owner[name]].append(name)
    typed = species_types({o: [seqs[k] for k in names] for o, names in by_owner.items()}, threads, log, "genome")
    out = {}
    for o, names in by_owner.items():
        types = typed[o]
        rep = {}  # type -> its representative (names are sorted: the first of the longest stays)
        for k in names:
            t = types[seqs[k]]
            if t not in rep or len(seqs[k]) > len(seqs[rep[t]]):
                rep[t] = k
        for k in names:  # the types in their first names' order, as before
            out.setdefault(rep[types[seqs[k]]], []).append(k)
    return out


def screen(records, ncbi, vsearch, threads=1, searcher=None, log=None):
    """records [(name, sequence, row)] (row: regions_list, accession, taxid or species_taxid) -> ({name: screen text},
    counts). A genome's copies of a part are screened as one per 99% type (representatives), its verdict theirs.
    searcher: a Searcher (rrna.finish's: kept and shared searches; since 2026-11-13), else one of this function's
    own; log: each part's steps and times."""
    if not vsearch or not records:
        return {}, collections.Counter()
    lineage = {n: lineage_of(ncbi, r.get("taxid") or r["species_taxid"]) for n, _s, r in records}
    genome = {n: genome_of(n, r) for n, _s, r in records}
    parts = {n: parts_of(s, r.get("regions_list")) for n, s, r in records}
    flags = collections.defaultdict(list)
    n = collections.Counter()
    with tempfile.TemporaryDirectory() as tmp:
        searcher = searcher or Searcher(vsearch, threads, tmp, log=log)
        for part in PARTS:
            began = time.time()
            every = {k: p[part] for k, p in parts.items() if part in p}
            reps = representatives(every, genome, searcher.threads, log)
            seqs = {k: every[k] for k in reps}
            if log:
                log(f"  rRNA screen {part}: {len(every):,} parts, {len(seqs):,} 99% types: took "
                    f"{util.duration(time.time() - began)}")
            if len(seqs) < 2:
                continue
            members = [collections.defaultdict(set) for _ in RANKS]
            for k in seqs:
                for rank in range(len(RANKS)):
                    key = clade(lineage[k], rank)
                    if key is not None:
                        members[rank][key].add(k)
            genomes_of = [{key: {genome[k] for k in names} for key, names in m.items()} for m in members]
            hits = searcher.one(f"screen {part}", seqs, seqs, SCREEN_ID, SCREEN_ACCEPTS, SCREEN_REJECTS)

            def tested_rank(k):
                for rank in range(len(RANKS) - 1, -1, -1):
                    key = clade(lineage[k], rank)
                    if key is not None and genomes_of[rank].get(key, set()) - {genome[k]}:
                        return rank, key
                return None, None

            def own_best(k, rank, key, found, skip=()):
                return max((i for i, t in found if clade(lineage[t], rank) == key
                            and genome[t] not in (genome[k], *skip)), default=0.0)

            suspects = {}
            for k in seqs:
                found = [(i, t) for i, t in hits.get(k, []) if genome[t] != genome[k]]
                if not found:
                    continue
                top = max(i for i, _t in found)
                near = [t for i, t in found if i >= top - NEAR]
                depth = max(agrees_to(lineage[k], lineage[t]) for t in near)
                rank, key = tested_rank(k)
                if rank is None or depth >= rank:
                    continue
                matched = max((it for it in found if it[1] in near), key=lambda it: (it[0], it[1]))[1]
                suspects[k] = (rank, key, top, matched, depth, own_best(k, rank, key, found))
            again = collections.defaultdict(set)
            for k, (rank, key, _top, matched, _d, _own) in suspects.items():
                again[(rank, key)].add(k)
                mkey = clade(lineage[matched], rank)
                if mkey is not None:
                    again[(rank, mkey)].add(matched)
            more = collections.defaultdict(list)
            asked = sorted(again.items())
            found = searcher.many(f"screen {part}: suspects among their clades", [
                ({t: seqs[t] for t in sorted(names)}, {t: seqs[t] for t in sorted(members[rank][key])})
                for (rank, key), names in asked], SCREEN_ID, SCREEN_ACCEPTS, SCREEN_REJECTS)
            for ((rank, _key), _names), got in zip(asked, found):
                for k, v in got.items():
                    more[(rank, k)] += v
            for k, (rank, key, top, matched, depth, own) in sorted(suspects.items()):
                own = max(own, own_best(k, rank, key, more[(rank, k)]))
                mkey = clade(lineage[matched], rank)
                theirs = own_best(matched, rank, mkey, hits.get(matched, []) + more[(rank, matched)],
                                  (genome[k],)) if mkey is not None else 0.0
                verdict = verdict_of(depth, top, own, theirs)
                n[f"{part}: {verdict.split(':')[0]}"] += len(reps[k])
                if verdict == "flagged":
                    for member in reps[k]:
                        flags[member].append(f"{part}:{RANKS[rank]}")
    return {k: ",".join(v) for k, v in flags.items()}, n
