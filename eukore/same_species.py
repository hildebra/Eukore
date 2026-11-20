"""Eukore's same-species merge (since 2026-11-20): an NCBI species named for an isolate ("Saccharomyces sp. M14",
"Blastocystis sp. NMH") whose genome is a named species' genome of its genus is that species' further genome, not a
species of its own.

Why: Eukore keeps every NCBI species with an assembly but NCBI's placeholders and crosses named "A b x C d"
(genomes.representatives). An isolate NCBI never named stays a species of its own, so a strain of Saccharomyces
cerevisiae filed as "Saccharomyces sp. M14" or "Saccharomyces sp. KACC 48362" was a second Saccharomyces species
beside S. cerevisiae: a rare-taxa audit of protal (2026-11-20) found such genomes sharing 58-89% of a named species'
protal-sketch entries (two species' genomes share under 1%). In protal their marker, rRNA and sketch evidence split
between two labels, and the package's own protal sketch (protal_sketch.py) left out the windows the two species'
representatives share, so that both lost entries.

What is compared (find): within each genus (NCBI's genus of the species; never across genera), the representative
genomes of the genus' species (the plan's jobs: NCBI's, a model species' and TOGA2's genomes, not a transcriptome nor
EukProt's records), by their sourmash sketches (metadata.SKETCH, made at each genome's run or by the metadata pass:
FracMinHash, k 21 and 31, scaled 10,000), at k = K (31), the genera a pool of processes at once. Of two genomes A and B:
- containment of A in B, c(A,B): the share of A's hashes B holds (sourmash's containment; at the coarser scale of the
  two where their scales differ);
- ANI: sourmash's max-containment ANI, max(c(A,B), c(B,A)) ** (1 / K): the identity of what the two share, read from
  the genome more wholly held in the other (sourmash's point estimate, 1 - its distance);
- share of A in B: c(A,B) / max(c(A,B), c(B,A)), the share of A's genome B holds at that identity (an aligned
  fraction): 1 where A is wholly in B. A k-mer containment mixes the two: at 95% ANI a whole genome's c is 0.95 ** 31
  = 0.20, and half a genome at 99% (a hybrid's parent) gives 0.5 x 0.73 = 0.37; the share tells them apart.
The sketches are read as sourmash writes them (gzipped JSON, read_sketch), no sourmash process a genus: sourmash's
containment of two scaled sketches is their shared hashes over one's own, which this module counts the same way.

The rules (the user's decisions of 2026-11-20; decide):
- an isolate-named species (ISOLATE: "<Genus> sp. <something>", or a bare "<Genus> sp." under its genus) whose
  genome is in a named species' at ANI_MIN (95%) or more with a share of SHARE_MIN (0.7) or more is merged into it:
  its genome becomes a further genome of the named species (genomes.tsv's species_taxid the named one's,
  representative f, source SOURCE), its copies, rRNA records and organelle genes the named species', and its NCBI
  taxid stays in genomes.tsv's taxid. The named species keeps its own representative. Of several named species it
  fits, the nearest (ANI, then containment) takes it, the others named in the row's note. Several isolates may merge
  into one species; a merge never chains (an isolate merges into a named species alone, and a named species never
  merges) and never crosses a genus;
- two named species at ANI_MIN or more: no merge, a row to review (NAMED_PAIR);
- a hybrid, kept as its own species, a row to review with both parents (HYBRID): a species named as one (HYBRID_NAME,
  "Saccharomyces x bayanus", "Mentha x piperita": a name genomes.HYBRID's "A b x C d" leaves in), or a genome that
  wholly holds two species of its genus that are not one species, without being wholly in either (Saccharomyces
  pastorianus: S. cerevisiae's and S. eubayanus' genomes each wholly in it, it half in each). The user has not
  decided what becomes of hybrids: they stay species. A hybrid is no merge target, and an isolate that is one is not
  merged;
- not decided either, and so rows to review, never merges: two isolate-named species at ANI_MIN or more, not both
  merged (ISOLATE_PAIR: Blastocystis's subtype isolates, of no named species; or an isolate near one that merged but
  not near enough the named species itself, never merged through it: the note says which merged), and an isolate at
  ANI_MIN or more to a named species but under SHARE_MIN of it held there (PARTIAL: an assembly holding another genome
  too, or a named species' assembly that is a part of its genome).

Thresholds: ANI_MIN 0.95 is the species boundary the user set (the prokaryotes' 95% ANI). SHARE_MIN 0.7: an assembly
of the same species holds 0.8-1.0 of another's (completeness, a second haplotype, a contaminant's contigs); a hybrid
holds half of each parent's genome, at most two thirds of one of two parents of 2:1 genome sizes. A raw containment
cut (0.5, say) would have meant ANI 97.8% or more at k 31 and let a recent hybrid's 0.4-0.5 near it. A sketch of
fewer than HASHES_MIN hashes (a genome under ~1 Mb at scaled 10,000) is too few to compare: its genome is left out of
the comparisons.

What the merge changes: package.collect labels the merged genome with the named species before the copies' checks
across species (qc.check: a copy is not confirmed by its own species' second genome, nor same_as its genus for it);
a merge into a named species whose representative the package does not hold (not done) or whose copies it leaves out
(tier poor, partial or contaminated) is undone there (UNDONE: the isolate stays a species, so that its copies are not
filed under a species of none). The build relabels the rRNA records and organelle genes of a merged species before
their checks (Merges.relabel), the intraspecific panel's jobs too, and writes TABLE (metadata/species_merges.tsv,
COLUMNS: a row per merge and per review) and the manifest's inputs.species_merge (Merges.about). The protal sketch
takes a merged genome as a further assembly of its species (protal_sketch.FURTHER). Without sketches (--no_sketches,
or no sourmash) or with --no_species_merge, nothing is compared and nothing merged: the package as before.
"""

import collections
import gzip
import json
import re
import zlib

import numpy as np

from . import store, util

K = 31                 # the sketches' k compared (metadata.SKETCH_PARAMS: k 21 and 31)
ANI_MIN = 0.95         # the species boundary: max-containment ANI
SHARE_MIN = 0.7        # an isolate's share of its genome held by the named species' (an aligned fraction)
HASHES_MIN = 100       # a sketch of fewer hashes is not compared (a genome under ~1 Mb at scaled 10,000)
CLOSE = ANI_MIN ** K   # the containment (of either genome in the other) of a pair at ANI_MIN: 0.204 at k 31
TABLE = "metadata/species_merges.tsv"
COLUMNS = ("from_taxid", "from_name", "from_accession", "to_taxid", "to_name", "to_accession", "ani", "containment",
           "share", "rule", "other_taxid", "other_name", "note")
SOURCE = "merged"      # genomes.tsv's source of a merged genome (package.MERGED, package.holds_copies)
FATE = "merged into a named species"  # benchmarks' selection funnel (genomes.fates)
MERGED = "merged"
NAMED_PAIR = "review: named pair"
HYBRID = "review: hybrid"
ISOLATE_PAIR = "review: isolate pair"
PARTIAL = "review: partial"
UNDONE = "review: target left out"
RULES = (MERGED, NAMED_PAIR, HYBRID, ISOLATE_PAIR, PARTIAL, UNDONE)
NOT_GENOMES = ("tsa", "no_assembly")  # a transcriptome, EukProt's records: not compared
# NCBI's name of a species it knows by an isolate alone: "<Genus> sp. <isolate>", or a genus' bare "<Genus> sp." (a
# bare "sp." above the genus is a placeholder, never a species of the package: ncbi.BARE_SP)
ISOLATE = re.compile(r"^(?:Candidatus\s+)?\S+\s+sp\.(?:\s|$)")
# a hybrid's name: an " x " (or "×") between two names ("Saccharomyces x bayanus", "Fragaria x ananassa"), or a
# nothogenus' "× Genus"
HYBRID_NAME = re.compile(r"(?:^|\s)[x×]\s")


def isolate_named(name):
    """Whether NCBI's name of a species is an isolate's (ISOLATE)."""
    return bool(ISOLATE.match(name or ""))


def hybrid_named(name):
    """Whether NCBI's name of a species names a hybrid (HYBRID_NAME)."""
    return bool(HYBRID_NAME.search(name or ""))


def kind_of(name):
    """A species' kind by its name: hybrid, isolate or named."""
    return "hybrid" if hybrid_named(name) else "isolate" if isolate_named(name) else "named"


# --- the sketches ---------------------------------------------------------------------------------------------------

def read_sketch(data, k=K):
    """A sourmash signature file's bytes (metadata.SKETCH: gzipped JSON, or plain: a list of signatures, each with its
    sketches under `signatures`, one a ksize, as sourmash 4 writes them; sourmash 2's several sketches under one
    signature too) -> (max_hash, its hashes: a sorted numpy uint64 array) of its scaled DNA sketch of k, or None (no
    such sketch, an unreadable file, None)."""
    if not data:
        return None
    try:
        got = json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)
        for sig in got if isinstance(got, list) else [got]:
            sketches = sig.get("signatures") if isinstance(sig, dict) else None
            for sketch in sketches if isinstance(sketches, list) else ():
                if isinstance(sketch, dict) and int(sketch.get("ksize") or 0) == k \
                        and str(sketch.get("molecule", "DNA")).upper() == "DNA" \
                        and int(sketch.get("max_hash") or 0) > 0:
                    return int(sketch["max_hash"]), np.unique(np.asarray(sketch.get("mins") or [], dtype=np.uint64))
    except (OSError, EOFError, ValueError, TypeError, OverflowError, zlib.error):  # not a signature file sourmash wrote
        return None
    return None


def containments(a, b):
    """Two sketches (read_sketch's) -> (the containment of a in b, of b in a): their shared hashes over each one's own,
    both cut to the coarser scale (the smaller max_hash) where their scales differ; 0 for an empty one."""
    top = min(a[0], b[0])
    x = a[1][:np.searchsorted(a[1], np.uint64(top), "right")]
    y = b[1][:np.searchsorted(b[1], np.uint64(top), "right")]
    if not len(x) or not len(y):
        return 0.0, 0.0
    shared = np.intersect1d(x, y, assume_unique=True).size
    return shared / len(x), shared / len(y)


def ani_of(c, k=K):
    """sourmash's ANI point estimate of a containment: c ** (1 / k) (0 for none)."""
    return c ** (1.0 / k) if c > 0 else 0.0


def genus_part(args):
    """A pool's task: (WORK, [(species taxid, its genome's accession)] of a genus) -> ({(a, b): containment of a's
    genome in b's} of the pairs whose max-containment ANI is ANI_MIN or more, both orders; [the species without a sketch
    of HASHES_MIN hashes or more]; the pairs compared). Each genome's sketch read from its done genome (store.members:
    a folder's or a packed genome's)."""
    from .genomes import done_dir  # noqa: PLC0415
    from .metadata import SKETCH  # noqa: PLC0415
    work, members = args
    sketches, without = {}, []
    for sp, acc in members:
        got = store.members(done_dir(work, acc), (SKETCH,))
        sketch = read_sketch(got.get(SKETCH)) if got else None
        if sketch is None or len(sketch[1]) < HASHES_MIN:
            without.append(sp)
            continue
        sketches[sp] = sketch
    held, n = {}, 0
    order = sorted(sketches)
    for i, a in enumerate(order):
        for b in order[i + 1:]:
            n += 1
            c_ab, c_ba = containments(sketches[a], sketches[b])
            if max(c_ab, c_ba) >= CLOSE:
                held[(a, b)], held[(b, a)] = c_ab, c_ba
    return held, without, n


# --- the rules --------------------------------------------------------------------------------------------------------

def decide(members, names, held):
    """A genus' rows (COLUMNS dicts) by the rules (module docstring). members: {species taxid: its genome's accession}
    of the genus' compared species; names: {species taxid: NCBI's name}; held: genus_part's {(a, b): containment of a in
    b} of its pairs at ANI_MIN or more (both orders). A species is never merged twice, nor merged and a target."""
    def close(a, b):
        return (a, b) in held

    def ani(a, b):
        return ani_of(max(held[(a, b)], held[(b, a)]))

    def share(a, b):
        top = max(held[(a, b)], held[(b, a)])
        return held[(a, b)] / top if top else 0.0

    def within(a, b):  # a's genome wholly (SHARE_MIN) in b's at ANI_MIN or more
        return close(a, b) and share(a, b) >= SHARE_MIN

    def row(a, b, rule, other=None, note="-"):
        out = dict(from_taxid=a, from_name=names.get(a, "-"), from_accession=members[a], rule=rule, note=note,
                   other_taxid="-", other_name="-", to_taxid="-", to_name="-", to_accession="-", ani="-",
                   containment="-", share="-")
        if b is not None:
            out.update(to_taxid=b, to_name=names.get(b, "-"), to_accession=members[b], ani=f"{ani(a, b):.4f}",
                       containment=f"{held[(a, b)]:.4f}", share=f"{share(a, b):.3f}")
        if other is not None:
            out.update(other_taxid=other, other_name=names.get(other, "-"))
            if b is None or other != b:
                out["note"] = (f"{names.get(other, other)}: ANI {ani(a, other):.4f}, containment "
                               f"{held[(a, other)]:.4f}, share {share(a, other):.3f}"
                               + ("" if note == "-" else f"; {note}"))
        return out

    def said(xs, a):
        return "; ".join(f"also {names.get(x, x)} ({x}): ANI {ani(a, x):.4f}, share {share(a, x):.3f}" for x in xs) \
            or "-"
    kind = {s: kind_of(names.get(s, "")) for s in members}
    nearest = {s: sorted((p for p in members if p != s and close(p, s)), key=lambda p: (-ani(s, p), p))
               for s in members}
    # the hybrids: two species of the genus wholly in a genome that is wholly in neither, the two not one species
    # (close to each other: two genomes of one species are each in the other, never a cross); the named parents first
    parents = {}
    for h in sorted(members):
        inside = sorted(p for p in members if p != h and within(p, h) and not within(h, p))
        best = None
        for i, p1 in enumerate(inside):
            for p2 in inside[i + 1:]:
                if close(p1, p2):
                    continue
                key = ((kind[p1] == "named") + (kind[p2] == "named"), min(ani(p1, h), ani(p2, h)), -p1, -p2)
                if best is None or key > best[0]:
                    best = (key, p1, p2)
        if best is not None:  # the nearer parent first
            parents[h] = tuple(sorted(best[1:], key=lambda p: (-ani(p, h), p)))
        elif kind[h] == "hybrid" and nearest[h]:  # named a hybrid: its nearest one or two species
            parents[h] = (nearest[h] + [None])[:2]
    hybrids = set(parents) | {s for s in members if kind[s] == "hybrid"}
    rows, merged = [], {}
    for q in sorted(members):
        if kind[q] != "isolate" or q in hybrids:
            continue
        named = [r for r in members if r != q and kind[r] == "named" and r not in hybrids]
        targets = sorted((r for r in named if within(q, r)), key=lambda r: (-ani(q, r), -held[(q, r)], r))
        if targets:
            merged[q] = targets[0]
            rows.append(row(q, targets[0], MERGED, note=said(targets[1:], q)))
            continue
        near = sorted((r for r in named if close(q, r)), key=lambda r: (-ani(q, r), -held[(q, r)], r))
        if near:  # its genome only in part in the named species': an assembly of another genome too, or theirs a part
            rows.append(row(q, near[0], PARTIAL, note=said(near[1:], q)))
    # two isolates of one genome, not both merged (one merged: the other, not near enough the named species itself,
    # is never merged through it: no chain), the more wholly held one first
    isolates = sorted(s for s in members if kind[s] == "isolate" and s not in hybrids)
    for i, a in enumerate(isolates):
        for b in isolates[i + 1:]:
            if close(a, b) and not (a in merged and b in merged):
                x, y = (a, b) if (share(a, b), -a) >= (share(b, a), -b) else (b, a)
                rows.append(row(x, y, ISOLATE_PAIR, note="; ".join(
                    f"{names.get(z, z)} merged into {names.get(merged[z], merged[z])} ({merged[z]})"
                    for z in (x, y) if z in merged) or "-"))
    named = sorted(s for s in members if kind[s] == "named" and s not in hybrids)
    for i, a in enumerate(named):
        for b in named[i + 1:]:
            if close(a, b):
                x, y = (a, b) if (share(a, b), -a) >= (share(b, a), -b) else (b, a)
                rows.append(row(x, y, NAMED_PAIR))
    for h in sorted(hybrids):
        p1, p2 = parents.get(h, (None, None))
        if p1 is not None:
            rows.append(row(h, p1, HYBRID, other=p2,
                            note="named a hybrid" if kind[h] == "hybrid" else "-"))
    return rows


# --- the build's step -------------------------------------------------------------------------------------------------

class Merges:
    """The step's outcome: rows (COLUMNS dicts, TABLE's), species {isolate species taxid: the named species' taxid} of
    the merges in force, about (the manifest's inputs.species_merge: state, thresholds, counts). False where nothing is
    merged."""

    def __init__(self, rows=(), about=None):
        self.rows = list(rows)
        self.species = {int(r["from_taxid"]): int(r["to_taxid"]) for r in self.rows if r["rule"] == MERGED}
        self.about = dict(about or {})
        self._count()

    def _count(self):
        if self.rows or "rules" in self.about:
            self.about["rules"] = {rule: n for rule in RULES
                                   if (n := sum(r["rule"] == rule for r in self.rows))}

    def __bool__(self):
        return bool(self.species)

    def of(self, taxid):
        """A species taxid as the package files it: the named species' where it was merged, else itself."""
        return self.species.get(int(taxid), int(taxid))

    def accessions(self):
        """{the merged genomes' accessions: the named species' taxid}."""
        return {r["from_accession"]: int(r["to_taxid"]) for r in self.rows if r["rule"] == MERGED}

    def undo(self, taxid, why):
        """A merge undone (package.collect: its named species' genome is not the package's, or its copies are left
        out): the isolate stays a species of its own, its row a review (UNDONE) saying why."""
        taxid = int(taxid)
        if taxid not in self.species:
            return
        del self.species[taxid]
        for r in self.rows:
            if r["rule"] == MERGED and int(r["from_taxid"]) == taxid:
                r["rule"], r["note"] = UNDONE, why + ("" if r["note"] in ("", "-") else f"; {r['note']}")
        self._count()

    def relabel(self, records):
        """Records [(name, sequence, row)] (rRNA records, organelle genes) with a merged species' species_taxid set to
        the named species' (its `taxid`, NCBI's label of the record, kept) -> the same list, changed in place."""
        if self.species:
            for _name, _seq, row in records:
                sp = row.get("species_taxid")
                if sp not in (None, "", "-") and int(sp) in self.species:
                    row["species_taxid"] = self.species[int(sp)] if isinstance(sp, int) else str(self.species[int(sp)])
        return records

    def jobs(self, jobs):
        """{accession: job} -> the same, a merged species' jobs copies with the named species' species_taxid (the
        intraspecific panel's further assemblies)."""
        if not self.species:
            return jobs
        return {a: dict(j, species_taxid=self.of(j["species_taxid"])) if int(j["species_taxid"]) in self.species
                else j for a, j in jobs.items()}


def off(state):
    """The step not run: no merge; the manifest's inputs.species_merge says why."""
    return Merges(about=dict(state=state))


def candidates(job_map, ncbi):
    """The genera the step compares -> {genus taxid: [(species taxid, its genome's accession)]}: those holding two
    species' genomes or more (the plan's jobs, not a transcriptome nor EukProt's records, nor a species of no genus)."""
    by_genus = collections.defaultdict(dict)
    for acc, job in sorted(job_map.items()):
        if job["why"][0] in NOT_GENOMES or job.get("path"):
            continue
        sp = int(job["species_taxid"])
        genus = ncbi.rank_ids(sp).get("genus")
        if genus is None or ncbi.rank.get(sp) != "species":
            continue
        by_genus[genus].setdefault(sp, acc)
    return {g: sorted(members.items()) for g, members in sorted(by_genus.items()) if len(members) >= 2}


def find(job_map, ncbi, work, threads=1, log=print):
    """The same-species merge's comparisons and rules over the plan's genomes (module docstring) -> Merges. threads:
    the genera compared that many at once by a pool of processes (genus_part), the result the same."""
    groups = candidates(job_map, ncbi)
    items = [(work, members) for members in groups.values()]
    rows, without, compared, sketched = [], 0, 0, 0
    parts = util.ordered(genus_part, items, max(1, min(threads, len(items)))) if items else ()
    for members, (held, lacking, n) in zip(groups.values(), parts):
        compared += n
        without += len(lacking)
        here = {sp: acc for sp, acc in members if sp not in lacking}
        sketched += len(here)
        if len(here) >= 2:
            rows += decide(here, {sp: ncbi.name.get(sp, str(sp)) for sp in here}, held)
    about = dict(state="made", table=TABLE, k=K, ani_min=ANI_MIN, share_min=SHARE_MIN, hashes_min=HASHES_MIN,
                 genera=len(groups), genomes=sketched, without_sketch=without, pairs=compared, rules={})
    merges = Merges(rows, about)
    counts = merges.about["rules"]
    examples = [f"{r['from_name']} -> {r['to_name']}" for r in rows if r["rule"] == MERGED][:3]
    log(f"same-species merge: {len(groups):,} genera of two species' genomes or more, {sketched:,} genomes "
        f"compared ({compared:,} pairs at k {K}" + (f"; {without:,} without a sketch, not compared" if without else "")
        + f"): {counts.get(MERGED, 0):,} isolate-named species merged into a named species"
        + (f" (e.g. {'; '.join(examples)})" if examples else "")
        + (", to review: " + ", ".join(f"{n:,} {rule[len('review: '):]}" for rule, n in counts.items()
                                       if rule != MERGED) if any(rule != MERGED for rule in counts) else "")
        + f" ({TABLE})")
    return merges
