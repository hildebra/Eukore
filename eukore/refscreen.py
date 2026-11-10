"""The reference screen (since 2026-11-10; KSGP 4.5's and 5.2's window screens, on NCBI's lineages and for any gene):
the records of a reference database Eukore reads (MIDORI2's genes) that are another organism's gene (a contaminant, a
misidentified specimen: `misnamed`), that hold a piece of another organism's gene (`chimera`) or that are another gene
(`other_gene`) are found, so that such a record names no species in species.tsv, and a user of the reference can leave
it out or cut it (tools/screen_reference.py writes KSGP's lists and the screened FASTA).

A record is told by its windows (WINDOW bases every STEP, the last one ending at the record's end; a shorter record is
one window), each set beside the same gene of other species:
- The references: Eukore's own curated records of the gene (`trusted`: RefSeq's organelle genomes' genes, by NCBI
  taxid) and the screened database's own records, one of each species (its longest: a species' other records are not
  independent of the record screened, and a database of every haplotype is searched as fast as one of a record a
  species). Records of other genes are references too (RefSeq's other organelle genes, the other gene files screened
  with it): a window nearer one of them is another gene's.
- Pass 1 (minimap2 -x asm20, its SECONDARY best hits): a record is flagged where a window is TRIGGER or more identical
  (over COVER of it) to a record of another gene, or of a species NCBI's lineage parts from its own at the order or
  above (RANKS; the first rank both name and differ at).
- Pass 2, the check, on the flagged records' windows every VERIFY_STEP: minimap2 gives each window's copies (records
  SAME or more identical), vsearch its relatives (FLOOR or more), and each window that looks another clade's is
  searched again among its own clade's records alone (a large clade's relatives fill a search's hits: KSGP's check of
  the label's phylum). A window is `claimed` by another clade where a record of it is a copy and GAP points nearer
  than the nearest record of the window's own clade (its own species' trusted record holding a copy settles it: its
  own). A window whose own clade's records are about as near tells the two apart no more than a conserved stretch does
  (plants' mitochondrial genes are 97-100% alike across orders, and their cox1's intron spreads by horizontal
  transfer, converting the exon beside it).
- Whose a claimed sequence is. The copy alone does not say: it may be the misnamed record (KSGP's second-genome rule),
  and so may several copies be (a batch of contaminated specimens are copies of their source under wrong names). A
  clade holds a sequence where it holds relatives of it, graded as a lineage's are: the other clade's relatives are
  its HOLD-th best species under SAME (a trusted record counts as two), the record's own clade's its nearest record.
  The consecutive windows one clade claims are a stretch, decided by all its windows (single windows' relatives at
  80-90% differ by chance; each value no lower than FLOOR, the search's limit): `foreign` where the other clade's
  relatives are REL_MARGIN points nearer on average than its own clade's nearest records, and no third clade holds
  half its windows about as well (within CONSERVED_MARGIN, or copies of two species); `own` where its own clade's are
  nearer by as much (the copy is the other clade's misnamed record); else `undecided`.
  The rank tested is the deepest of RANKS at which the references hold another species of the record's clade (a lone
  order's record is tested at its class: nothing of its order could be near).
- The call. A foreign stretch is certain where it is STRONG_RUN windows or more, or its copy STRONG or more identical.
  A record with a certain stretch is a `chimera` where STRONG_RUN other windows are its own clade's (cut to its longest
  part outside the stretches where MIN_KEEP bases or more stay, else dropped) or the stretches are two clades';
  `misnamed` where they are not and the stretches are MISNAMED_SHARE of its windows or more; `other_gene` likewise
  (whole, or a piece). `weak` (for review, kept as it is): an undecided stretch (one of two records is misnamed,
  which is not told), a foreign one that is not certain, or a piece in a record nothing else places. Everything else
  the first pass flagged is `cleared`.

What the screen does not see: clades NCBI does not part at the order or above (a misidentified congener, a chimera of
two species of an order); a contaminant of a clade the references hold no relatives of; records of one sequence under
names of the record's own clade vouch for it (a batch of contaminated specimens of one order: the window's own clade
"holds a copy").
"""

import collections
import os
import shutil
import subprocess
import tempfile
import threading

from . import util

WINDOW, STEP, VERIFY_STEP = 300, 150, 50
MIN_LENGTH = 150      # shorter records are not screened
COVER = 0.75          # of a window a hit covers at least
TRIGGER = 95.0        # pass 1: a window this identical to another clade's or gene's record flags its record
SECONDARY = 8         # pass 1: minimap2's secondary hits a window
SAME = 97.0           # a record this identical to a window is a copy of its sequence
FLOOR = 80.0          # the check's hits (vsearch --id): a relative is this identical at least
STRONG = 99.0         # a foreign stretch of one window whose copy is this identical is certain
STRONG_RUN = 2        # ... or one of this many windows
GAP = 5.0             # a window's copy in another clade this much nearer than its own clade's nearest record
HOLD = 2              # another clade's relatives: its best species under SAME, this many (a trusted one: two)
REL_MARGIN = 3.0      # ... nearer on average than the record's own clade's nearest records by this many points
CONSERVED_MARGIN = 2.0  # a third clade holding it this near: a conserved stretch
GENE_MARGIN = 5.0     # another gene's record this much nearer than any record of the window's gene
OWN_ID = 85.0         # a window this identical to a record of its own clade, no other clade's much nearer: its own
MIN_KEEP = 300        # a chimera cut to its own part keeps this many bases at least
MISNAMED_SHARE = 0.5  # foreign windows this share of a record's windows, none its own: misnamed
COPIES = 200         # the check: minimap2's secondary hits a window (its copies)
ACCEPTS, REJECTS = 100, 400  # the check: vsearch's hits a window (its relatives), and the rejects it tries
RANKS = ("phylum", "class", "order")
VERSION = "1"         # raised where the screen's rules change: a build's kept screen of another version is run again
CALLS = ("misnamed", "chimera", "other_gene", "weak", "cleared")
CALL_COLUMNS = ("record", "gene", "taxid", "species_taxid", "length", "call", "action", "keep_start", "keep_end",
                "clade", "label_clade", "identity", "relatives", "own_clade", "windows", "layout", "targets")
TABLE_COLUMNS = ("reference",) + CALL_COLUMNS + ("reason",)  # the package's reference_screen.tsv

Record = collections.namedtuple("Record", "name taxid gene sequence")
Verdict = collections.namedtuple("Verdict", "kind clade identity target held own third",
                                 defaults=(None, 0.0, None, 0.0, 0.0, False))


def window_starts(length, size=WINDOW, step=STEP):
    """A record's window starts: every `step`, the last window ending at its end; one window for a short record."""
    if length <= size:
        return [0]
    starts = list(range(0, length - size + 1, step))
    if starts[-1] + size < length:
        starts.append(length - size)
    return starts


def parting(a, b):
    """The first of RANKS (its index) at which two lineages both name a taxon and differ, else None."""
    for i, (x, y) in enumerate(zip(a, b)):
        if x is not None and y is not None and x != y:
            return i
    return None


class Taxa:
    """NCBI's species and lineage at RANKS of taxids (cached); a taxid the taxonomy lacks has neither."""

    def __init__(self, ncbi):
        self.ncbi, self._species, self._lineage = ncbi, {}, {}

    def species(self, taxid):
        if taxid not in self._species:
            t = self.ncbi.merged.get(taxid, taxid)
            self._species[taxid] = self.ncbi.species_of(t) if t in self.ncbi.parent else None
        return self._species[taxid]

    def lineage(self, species):
        if species not in self._lineage:
            at = self.ncbi.rank_ids(species)
            self._lineage[species] = tuple(at.get(r) for r in RANKS)
        return self._lineage[species]

    def clade_name(self, clade):
        return "-" if clade is None else f"{RANKS[clade[0]]} {self.ncbi.name.get(clade[1], clade[1])}"


class References:
    """The records windows are set beside: their species, gene, lineage and whether Eukore trusts them, by index; the
    species of each (gene, rank, taxon) among them (members: who else of a clade holds the gene)."""

    def __init__(self):
        self.name, self.species, self.gene, self.lineage, self.trusted, self.twice = [], [], [], [], [], []
        self.members = collections.defaultdict(set)

    def add(self, name, species, gene, lineage, trusted, twice=False):
        """twice: a database's record filed under two genes (one accession and span in two gene files): no witness
        of either gene."""
        self.name.append(name)
        self.species.append(species)
        self.gene.append(gene)
        self.lineage.append(lineage)
        self.trusted.append(trusted)
        self.twice.append(twice)
        for rank, taxon in enumerate(lineage):
            if taxon is not None:
                self.members[(gene, rank, taxon)].add(species)
        return len(self.name) - 1

    def tested(self, species, gene, lineage):
        """The deepest of RANKS (its index) at which the references hold the gene of another species of the
        lineage's clade, else None: the rank a record of that species is tested at."""
        for rank in range(len(RANKS) - 1, -1, -1):
            if lineage[rank] is not None and self.members.get((gene, rank, lineage[rank]), set()) - {species}:
                return rank
        return None


def kth_species(hits, ref, k=HOLD):
    """The identity of the k-th best species among hits [(identity, target)] (a species its best hit; one of a
    trusted record counts as two), 0.0 where they hold fewer."""
    best = {}
    for identity, t in hits:
        weight = 2 if ref.trusted[t] else 1
        sp = ref.species[t]
        if sp not in best or (identity, weight) > best[sp]:
            best[sp] = (identity, weight)
    n = 0
    for identity, weight in sorted(best.values(), reverse=True):
        n += weight
        if n >= k:
            return identity
    return 0.0


def judge(species, gene, lineage, tested, hits, ref):
    """A window's verdict (module docstring) from its hits [(identity %, reference index)]: Verdict(kind: `own` (its
    own clade's records the nearest), `claim` (a clade apart from its own holds a copy of it, GAP nearer than its own
    clade's nearest record: clade (rank index, taxon), identity and target its best copy's, held the clade's
    relatives, own its own clade's nearest record, third whether a third clade holds it about as well), `gene`
    (another gene's: clade the gene) or `none`)."""
    hits = [(i, t) for i, t in hits if ref.trusted[t] or ref.species[t] != species or ref.gene[t] != gene]
    mine = [(i, t) for i, t in hits if ref.gene[t] == gene]
    other = [(i, t) for i, t in hits if ref.gene[t] != gene]
    if other:
        # the trusted records tell the gene where they are near: a database's records can be misfiled alike (MIDORI2
        # files one record under lrRNA and srRNA: 9,189 of GB273's); else its records filed once, of two species
        sure_other = max(((i, t) for i, t in other if ref.trusted[t]), default=(0.0, None))
        sure_mine = max((i for i, t in mine if ref.trusted[t]), default=0.0)
        if sure_other[0] >= sure_mine + GENE_MARGIN:
            return Verdict("gene", ref.gene[sure_other[1]], sure_other[0], sure_other[1])
        once = [(i, t) for i, t in other if not ref.twice[t]]
        if not sure_mine and once:
            top, target = max(once)
            of_g = {ref.species[t] for _i, t in once if ref.gene[t] == ref.gene[target]}
            if len(of_g) >= 2 and top >= max((i for i, t in mine if not ref.twice[t]), default=0.0) + GENE_MARGIN:
                return Verdict("gene", ref.gene[target], top, target)
    if not mine:
        return Verdict("none")

    def clade_of(t):
        p = parting(lineage, ref.lineage[t])
        return None if p is None or tested is None or p > tested else (p, ref.lineage[t][p])
    own, foreign = [], collections.defaultdict(list)
    for i, t in mine:
        k = clade_of(t)
        (own if k is None else foreign[k]).append((i, t))
    own_top = max((i for i, _t in own), default=0.0)
    other_top = max((i for hs in foreign.values() for i, _t in hs), default=0.0)
    claims = {k: hs for k, hs in foreign.items() if max(hs)[0] >= SAME}
    trusted_own = any(i >= SAME and ref.trusted[t] and ref.species[t] == species for i, t in own)
    if not claims or trusted_own or max(max(hs)[0] for hs in claims.values()) - own_top < GAP:
        # no other clade's copy, a trusted record of its own species holding it, or its own clade's records about as
        # near as the copy (the window does not tell the two apart)
        return Verdict("own", None, own_top) if own_top >= OWN_ID and own_top >= other_top - GAP else Verdict("none")

    def relatives(hs):
        return kth_species([(i, t) for i, t in hs if i < SAME], ref)

    def copies(hs):
        return len({ref.species[t] for i, t in hs if i >= SAME})
    held = {k: relatives(hs) for k, hs in claims.items()}
    best = max(claims, key=lambda k: (held[k], copies(claims[k]), max(claims[k])[0], k))
    top, target = max(claims[best])
    third = any(k != best and ((held[k] > 0 and held[k] >= held[best] - CONSERVED_MARGIN) or copies(claims[k]) >= 2)
                for k in claims)
    return Verdict("claim", best, top, target, held[best], own_top, third)


def kept_part(removed, length):
    """The longest stretch (start, end) of a record outside the removed spans [(start, end)]."""
    pieces, at = [], 0
    for start, end in sorted(removed):
        if start > at:
            pieces.append((at, start))
        at = max(at, end)
    if at < length:
        pieces.append((at, length))
    return max(pieces, key=lambda p: p[1] - p[0], default=(0, 0))


def stretches(verdicts, step=VERIFY_STEP):
    """A record's windows' verdicts {start: Verdict} as stretches of consecutive windows of one kind (and one clade or
    gene) -> [dict(kind, clade, first, last: window starts, n, identity, target, held, own, third)]. A claim's
    stretch is decided by all its windows (each one's relatives at 80-90% differ by chance): `foreign` where the other
    clade's relatives are REL_MARGIN nearer on average than its own clade's nearest records (each no lower than
    FLOOR, the search's limit) and a third clade holds under half of its windows; `own` where its own clade's are
    (the copy is the other clade's misnamed record); else `undecided`."""
    out = []
    for s in sorted(verdicts):
        v = verdicts[s]
        key = (v.kind, v.clade if v.kind in ("claim", "gene") else None)
        if out and (out[-1]["kind"], out[-1]["clade"]) == key and s - out[-1]["last"] <= step:
            out[-1]["last"] = s
            out[-1]["windows"].append(v)
        else:
            out.append(dict(kind=key[0], clade=key[1], first=s, last=s, windows=[v]))
    for r in out:
        vs = r.pop("windows")
        best = max(vs, key=lambda v: v.identity)
        r.update(n=len(vs), identity=best.identity, target=best.target, held=0.0, own=0.0, third=False)
        if r["kind"] == "claim":
            r["held"] = sum(max(v.held, FLOOR) for v in vs) / len(vs)
            r["own"] = sum(max(v.own, FLOOR) for v in vs) / len(vs)
            r["third"] = 2 * sum(v.third for v in vs) >= len(vs)
            r["kind"] = ("foreign" if r["held"] - r["own"] >= REL_MARGIN and not r["third"] else
                         "own" if r["own"] - r["held"] >= REL_MARGIN else "undecided")
    return out


def call_of(length, verdicts, step=VERIFY_STEP, window=WINDOW):
    """A record's call from its windows' verdicts {start: Verdict} -> dict(call, action: drop, trim, review or -,
    keep (start, end) or None, clade, identity, relatives, own_clade: its main stretch's, layout, targets, windows:
    the stretches' kinds' windows)."""
    runs = stretches(verdicts, step)
    n_windows = len(verdicts)
    counts = collections.Counter()
    for r in runs:
        counts[r["kind"]] += r["n"]
    layout = " ".join(f"{r['kind']}:{r['first']}-{min(r['last'] + window, length)}" for r in runs)
    out = dict(call="cleared", action="-", keep=None, clade=None, identity=0.0, relatives=0.0, own_clade=0.0,
               layout=layout, targets=[], windows=dict(counts))

    def sure(r):
        return r["n"] >= STRONG_RUN or n_windows == 1 or (r["kind"] == "foreign" and r["identity"] >= STRONG)
    certain = [r for r in runs if r["kind"] in ("foreign", "gene") and sure(r)]
    doubt = [r for r in runs if r["kind"] in ("foreign", "gene", "undecided") and r not in certain]
    main = max(certain or doubt, key=lambda r: (r["n"], r["identity"]), default=None)
    if main is None:
        return out
    out.update(clade=main["clade"], identity=main["identity"], relatives=main["held"], own_clade=main["own"],
               targets=sorted({r["target"] for r in (certain or [main]) if r["target"] is not None}))
    if not certain:
        out.update(call="weak", action="review")
        return out
    genes = [r for r in certain if r["kind"] == "gene"]
    kind = "other_gene" if len(genes) == len(certain) else "chimera"
    removed = [(r["first"], min(r["last"] + window, length)) for r in certain]
    if counts["own"] >= STRONG_RUN:  # a part of its own clade as certain as a piece: one window alone is not
        start, end = kept_part(removed, length)
        out.update(call=kind, action="trim", keep=(start, end)) if end - start >= MIN_KEEP else \
            out.update(call=kind, action="drop")
    elif len({r["clade"] for r in certain if r["kind"] == "foreign"}) >= 2:
        out.update(call="chimera", action="drop")
    elif sum(r["n"] for r in certain) >= MISNAMED_SHARE * n_windows:
        out.update(call="other_gene" if genes else "misnamed", action="drop")
    else:  # a piece of another clade's in a record nothing else places
        out.update(call="weak", action="review")
    return out


# --- the searches ---------------------------------------------------------------------------------------------------

def minimap2_hits(minimap2, reference, windows, threads=1, log=None, secondary=SECONDARY):
    """minimap2 -x asm20 of windows ((name, sequence), fed through a pipe) on a reference FASTA -> (name, {target:
    (identity %, share of the window covered)}) of each window with hits (a window's pieces on one target summed;
    minimap2 writes a window's alignments together)."""
    cmd = [minimap2, "-x", "asm20", "-c", "-t", str(max(1, threads)), "--secondary=yes", "-N", str(secondary),
           "-p", "0.5", reference, "-"]
    with open(log or os.devnull, "a") as err:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err, text=True,
                                bufsize=1 << 20)
        failed = []

        def feed():
            try:
                for name, seq in windows:
                    proc.stdin.write(f">{name}\n{seq}\n")
            except (BrokenPipeError, OSError) as e:  # minimap2 died: its exit code says so below
                failed.append(e)
            finally:
                try:
                    proc.stdin.close()
                except OSError:
                    pass
        feeder = threading.Thread(target=feed, daemon=True)
        feeder.start()
        name, sums = None, {}

        def result(sums):
            return {t: (100.0 * m / b, min(1.0, c)) for t, (m, b, c) in sums.items() if b}
        for line in proc.stdout:
            f = line.split("\t", 12)
            if len(f) < 11:
                continue
            if f[0] != name:
                if name is not None:
                    yield name, result(sums)
                name, sums = f[0], {}
            m, b, c = sums.get(f[5], (0, 0, 0.0))
            sums[f[5]] = (m + int(f[9]), b + int(f[10]), c + (int(f[3]) - int(f[2])) / int(f[1]))
        if name is not None:
            yield name, result(sums)
        feeder.join()
        proc.stdout.close()
        if proc.wait():
            raise OSError(f"minimap2 failed (exit {proc.returncode})" + (f", see {log}" if log else ""))


def vsearch_hits(vsearch, queries, reference, ident, threads=1, accepts=ACCEPTS, rejects=REJECTS, maxid=None,
                 tmp=None):
    """vsearch --usearch_global of a query FASTA on a reference FASTA (hits covering COVER of the query, of `ident` %
    or more, under `maxid` % where given) -> {query: [(identity %, target)]}."""
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        out = os.path.join(work, "hits.tsv")
        cmd = [vsearch, "--usearch_global", queries, "--db", reference, "--id", f"{ident / 100:.4f}", "--iddef", "2",
               "--maxaccepts", str(accepts), "--maxrejects", str(rejects), "--query_cov", str(COVER),
               "--strand", "plus", "--userout", out, "--userfields", "query+target+id",
               "--threads", str(max(1, threads)), "--quiet"]
        if maxid is not None:
            cmd += ["--maxid", f"{maxid / 100:.4f}"]
        run = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if run.returncode != 0 or not os.path.isfile(out):
            raise OSError(f"vsearch failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        hits = collections.defaultdict(list)
        with open(out) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3:
                    hits[f[0]].append((float(f[2]), f[1]))
        return hits


# --- the screen -----------------------------------------------------------------------------------------------------

def screen(sources, trusted, ncbi, minimap2, vsearch, threads=1, tmp=None, log=None, keep=None):
    """sources: functions, each yielding a database's Records (name, taxid, gene, sequence; read three times: they
    are not held in memory); trusted: Eukore's curated Records of the genes -> (calls {(gene, record name): dict of
    CALL_COLUMNS} of every record the first pass flagged, counts). A record is one gene's (a name two gene files hold
    is two records). keep: a folder for the searches' files (else a temporary one, deleted)."""
    log = log or (lambda _text: None)
    taxa, ref = Taxa(ncbi), References()
    counts = collections.Counter()
    with tempfile.TemporaryDirectory(dir=tmp) as temporary:
        work = keep or temporary
        os.makedirs(work, exist_ok=True)
        reference = os.path.join(work, "reference.fna")
        # the references: the trusted records, and one record (the longest) of each species and gene of the sources
        # (gene, name) -> (species, gene, length, taxid); (species, gene) -> (length, name, sequence)
        meta, longest = {}, {}
        for source in sources:
            for r in source():
                counts["records"] += 1
                sp = taxa.species(r.taxid) if r.taxid is not None else None
                if sp is None:
                    counts["no_taxid"] += 1
                    continue
                if len(r.sequence) < MIN_LENGTH:
                    counts["short"] += 1
                    continue
                if (r.gene, r.name) in meta:
                    counts["duplicate_names"] += 1
                    continue
                meta[(r.gene, r.name)] = (sp, r.gene, len(r.sequence), r.taxid)
                held = longest.get((sp, r.gene))
                if held is None or len(r.sequence) > held[0]:
                    longest[(sp, r.gene)] = (len(r.sequence), r.name, r.sequence.upper())
        with open(reference, "w") as fh:
            for r in trusted:
                sp = taxa.species(r.taxid) if r.taxid is not None else None
                if sp is None or len(r.sequence) < MIN_LENGTH:
                    continue
                i = ref.add(r.name, sp, r.gene, taxa.lineage(sp), True)
                fh.write(f">{i}\n{r.sequence.upper()}\n")
            counts["trusted"] = len(ref.name)
            filed = collections.Counter(name for _gene, name in meta)  # a name under two genes: filed twice
            counts["filed_twice"] = sum(1 for n in filed.values() if n > 1)
            for (sp, gene), (_n, name, seq) in sorted(longest.items()):
                i = ref.add(name, sp, gene, taxa.lineage(sp), False, filed[name] > 1)
                fh.write(f">{i}\n{seq}\n")
        del longest
        counts["references"] = len(ref.name)
        counts["screened"] = len(meta)
        log(f"reference screen: {counts['screened']:,} records of {counts['records']:,} to screen ("
            f"{counts['no_taxid']:,} of no taxid NCBI has, {counts['short']:,} under {MIN_LENGTH} bases), against "
            f"{counts['trusted']:,} trusted records and {counts['references'] - counts['trusted']:,} species' own")

        # pass 1: every window's best hits; a record with a window near another gene's or another clade's is flagged
        names = sorted(meta)
        index = {name: i for i, name in enumerate(names)}

        def windows():
            seen = set()
            for source in sources:
                for r in source():
                    key = (r.gene, r.name)
                    if key in seen or key not in meta or meta[key][2] != len(r.sequence):
                        continue
                    seen.add(key)
                    seq = r.sequence.upper()
                    for s in window_starts(len(seq)):
                        counts["windows"] += 1
                        yield f"{index[key]}:{s}", seq[s:s + WINDOW]
        flagged = set()
        for name, hits in minimap2_hits(minimap2, reference, windows(), threads, os.path.join(work, "minimap2.log")):
            record = names[int(name.split(":", 1)[0])]
            if record in flagged:
                continue
            sp, gene, _n, _t = meta[record]
            lineage = taxa.lineage(sp)
            for target, (identity, cover) in hits.items():
                t = int(target)
                if identity < TRIGGER or cover < COVER or (ref.species[t] == sp and not ref.trusted[t]
                                                           and ref.gene[t] == gene):
                    continue
                if ref.gene[t] != gene or parting(lineage, ref.lineage[t]) is not None:
                    flagged.add(record)
                    break
        counts["flagged"] = len(flagged)
        log(f"reference screen: {counts['windows']:,} windows; the first pass flagged {len(flagged):,} records")
        calls = check(sources, flagged, meta, ref, taxa, reference, minimap2, vsearch, threads, work, counts)
    for c in calls.values():
        counts[c["call"]] += 1
        counts[f"action_{c['action']}"] += 1
    return calls, counts


def check(sources, flagged, meta, ref, taxa, reference, minimap2, vsearch, threads, work, counts):
    """Pass 2: the flagged records' windows every VERIFY_STEP against the references (minimap2: their copies;
    vsearch: the relatives of the windows of records still in doubt), each judged, and each record called ->
    {(gene, record name): dict of CALL_COLUMNS}."""
    if not flagged:
        return {}
    queries = os.path.join(work, "check.fna")
    order, lengths = [], {}
    with open(queries, "w") as fh:
        for source in sources:
            for r in source():
                key = (r.gene, r.name)
                if key not in flagged or meta[key][2] != len(r.sequence) or key in lengths:
                    continue
                seq = r.sequence.upper()
                lengths[key] = len(seq)
                for s in window_starts(len(seq), WINDOW, VERIFY_STEP):
                    order.append((key, s))
                    fh.write(f">{len(order) - 1}\n{seq[s:s + WINDOW]}\n")
    counts["check_windows"] = len(order)
    # each window's copies and near records (minimap2: every hit of TRIGGER or more); the records with a window that a
    # clade apart from theirs holds a copy of, and no other species of their own: their windows' relatives (vsearch)
    hits = collections.defaultdict(list)
    with open(queries) as fh:
        pairs = ((line[1:].strip(), next(fh).strip()) for line in fh)
        for q, found in minimap2_hits(minimap2, reference, pairs, threads, os.path.join(work, "minimap2.log"),
                                      secondary=COPIES):
            hits[int(q)] = [(identity, int(t)) for t, (identity, cover) in found.items()
                            if cover >= COVER and identity >= TRIGGER]

    def judged(i, found):
        name, _s = order[i]
        sp, gene, _n, _t = meta[name]
        lineage = taxa.lineage(sp)
        return judge(sp, gene, lineage, ref.tested(sp, gene, lineage), found, ref)
    doubted = {order[i][0] for i in range(len(order)) if judged(i, hits.get(i, ())).kind in ("claim", "gene")}
    counts["checked"] = len(doubted)
    if doubted:
        again = os.path.join(work, "relatives.fna")
        with open(queries) as fh, open(again, "w") as out:
            for line in fh:
                seq = next(fh)
                if order[int(line[1:])][0] in doubted:
                    out.write(line + seq)
                    counts["relatives_windows"] += 1
        for q, found in vsearch_hits(vsearch, again, reference, FLOOR, threads, ACCEPTS, REJECTS, maxid=TRIGGER - 0.01,
                                     tmp=work).items():
            hits[int(q)] += [(identity, int(t)) for identity, t in found]
    verdicts = collections.defaultdict(dict)
    for i, (name, s) in enumerate(order):
        verdicts[name][s] = judged(i, hits.get(i, ()))
    # a window's hits can all be other clades' (a large clade's relatives fill them): each window that is not its
    # own clade's is searched again among its own clade's records alone (KSGP's check of the label's phylum)
    groups = collections.defaultdict(list)  # (gene, rank, taxon) -> windows
    for i, (name, s) in enumerate(order):
        if verdicts[name][s].kind == "claim":
            sp, gene, _n, _t = meta[name]
            lineage = taxa.lineage(sp)
            rank = ref.tested(sp, gene, lineage)
            groups[(gene, rank, lineage[rank])].append(i)
    counts["own_clade_windows"] = sum(len(v) for v in groups.values())
    if groups:
        wanted = {i for v in groups.values() for i in v}
        windows = {}
        with open(queries) as fh:
            for line in fh:
                seq = next(fh)
                if int(line[1:]) in wanted:
                    windows[int(line[1:])] = seq
        folder = os.path.join(work, "own_clades")
        os.makedirs(folder, exist_ok=True)
        files = {key: open(os.path.join(folder, f"{n}.fna"), "w") for n, key in enumerate(sorted(groups))}
        with open(reference) as fh:
            for line in fh:
                seq, t = next(fh), int(line[1:])
                for rank, taxon in enumerate(ref.lineage[t]):
                    fh_out = files.get((ref.gene[t], rank, taxon))
                    if fh_out is not None:
                        fh_out.write(line + seq)
        for fh in files.values():
            fh.close()
        for n, key in enumerate(sorted(groups)):
            own_queries = os.path.join(folder, f"{n}.q.fna")
            with open(own_queries, "w") as fh:
                fh.writelines(f">{i}\n{windows[i]}" for i in groups[key])
            for q, found in vsearch_hits(vsearch, own_queries, os.path.join(folder, f"{n}.fna"), FLOOR, threads,
                                         ACCEPTS, REJECTS, tmp=work).items():
                name, s = order[int(q)]
                hits[int(q)] += [(identity, int(t)) for identity, t in found]
                verdicts[name][s] = judged(int(q), hits[int(q)])
        shutil.rmtree(folder, ignore_errors=True)
    calls = {}
    for name in sorted(verdicts):
        sp, gene, length, taxid = meta[name]
        c = call_of(length, verdicts[name])
        lineage = taxa.lineage(sp)
        rank = c["clade"][0] if isinstance(c["clade"], tuple) else None
        keep = c["keep"] or ("-", "-")
        calls[name] = dict(
            record=name[1], gene=gene, taxid=taxid, species_taxid=sp, length=length, call=c["call"], action=c["action"],
            keep_start=keep[0] + 1 if c["keep"] else "-", keep_end=keep[1],
            clade=taxa.clade_name(c["clade"]) if isinstance(c["clade"], tuple) else (c["clade"] or "-"),
            label_clade=taxa.clade_name((rank, lineage[rank])) if rank is not None and lineage[rank] else "-",
            identity=f"{c['identity']:.1f}" if c["identity"] else "-",
            relatives=f"{c['relatives']:.1f}" if c["relatives"] else "-",
            own_clade=f"{c['own_clade']:.1f}" if c["own_clade"] else "-",
            windows=",".join(f"{k}:{n}" for k, n in sorted(c["windows"].items())), layout=c["layout"],
            targets=" ".join(ref.name[t] for t in c["targets"] if t is not None) or "-")
    return calls


def trusted_records(rows, organelle="mitochondrion"):
    """RefSeq's organelle genes (organelle.refseq_genes' rows) of an organelle as trusted Records: a gene's exons
    joined where it has introns (a database's records of a coding gene are its coding sequence)."""
    for r in rows:
        if r["organelle"] == organelle and str(r["taxid"]).isdigit():
            yield Record(f"{r['seqid']}_{r['gene']}", int(r["taxid"]), r["gene"], r.get("spliced") or r["sequence"])


def midori2_source(path):
    """A MIDORI2 RAW FASTA as a source of screen(): a function yielding its Records, its gene by the file's name."""
    from . import labels  # noqa: PLC0415
    gene = labels.midori2_gene(path)

    def read():
        for accession, taxid, seq in labels.midori2_records(path):
            yield Record(accession, taxid, gene, seq)
    return read


def read_calls(path):
    """A screen's calls.tsv -> {(gene, record name): dict of CALL_COLUMNS} (write_calls')."""
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    return {(r["gene"], r["record"]): r for r in rows}


def table_rows(reference, calls):
    """The package's reference_screen.tsv rows (TABLE_COLUMNS) of a reference's calls: those not cleared, each with
    its reason in words."""
    return [dict(c, reference=reference, reason=reason(c))
            for c in sorted(calls.values(), key=lambda c: (CALLS.index(c["call"]), c["gene"], c["record"]))
            if c["call"] != "cleared"]


def stamp(counts):
    """A screen's counts for the manifest (a reference's `screen`)."""
    keys = ("records", "screened", "flagged", "misnamed", "chimera", "other_gene", "weak")
    out = {k: int(counts.get(k, 0)) for k in keys}
    out.update(version=VERSION, dropped=int(counts.get("action_drop", 0)), cut=int(counts.get("action_trim", 0)),
               review=int(counts.get("action_review", 0)))
    return out


def reason(c):
    """A call in words (the lists' reason)."""
    windows = f"windows of {WINDOW} bases every {VERIFY_STEP}: {c['layout']}"
    near = f"{c['identity']}% to {c['targets']}"
    if c.get("relatives", "-") != "-":
        near += f"; that clade's relatives {c['relatives']}%, its own clade's nearest {c['own_clade']}%"
    if c["call"] == "misnamed":
        return f"the {c['clade']}'s gene under a name of the {c['label_clade']} ({near}); {windows}"
    if c["call"] == "chimera":
        cut = f", cut to {c['keep_start']}-{c['keep_end']}" if c["action"] == "trim" else ""
        return f"a piece of the {c['clade']}'s gene in a record of the {c['label_clade']} ({near}){cut}; {windows}"
    if c["call"] == "other_gene":
        cut = f", cut to {c['keep_start']}-{c['keep_end']}" if c["action"] == "trim" else ""
        return f"{c['clade']}, not {c['gene']} ({near}){cut}; {windows}"
    if c["call"] == "weak":
        return f"{c['clade']} ({near}), not certain; {windows}"
    return windows


def dropped(calls):
    """The records a screen leaves out: {(gene, record name)}."""
    return {key for key, c in calls.items() if c["action"] == "drop"}


def write_calls(path, calls):
    """calls -> a table of CALL_COLUMNS, the flagged records by call and name."""
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(CALL_COLUMNS) + "\n")
        for c in sorted(calls.values(), key=lambda c: (CALLS.index(c["call"]), c["gene"], c["record"])):
            fh.write("\t".join(str(c[k]) for k in CALL_COLUMNS) + "\n")
    os.replace(util.part_of(path), path)


def write_lists(folder, calls):
    """calls.tsv and KSGP's lists of a screen: drop.tsv, trim.tsv (the part kept, 1-based and inclusive) and
    review.tsv."""
    os.makedirs(folder, exist_ok=True)
    write_calls(os.path.join(folder, "calls.tsv"), calls)
    lists = {"drop": ("record\tgene\tcall\treason", []),
             "trim": ("record\tgene\tstart\tend\tlength\treason", []),
             "review": ("record\tgene\tcall\treason", [])}
    for key in sorted(calls):
        c = calls[key]
        name = f"{c['record']}\t{c['gene']}"
        if c["action"] == "trim":
            lists["trim"][1].append(f"{name}\t{c['keep_start']}\t{c['keep_end']}\t{c['length']}\t{reason(c)}")
        elif c["action"] in lists:
            lists[c["action"]][1].append(f"{name}\t{c['call']}\t{reason(c)}")
    for kind, (head, lines) in lists.items():
        path = os.path.join(folder, f"{kind}.tsv")
        with open(util.part_of(path), "w", newline="\n") as fh:
            fh.write(head + "\n" + "".join(line + "\n" for line in lines))
        os.replace(util.part_of(path), path)


def summary(calls, counts, files=()):
    """A screen's counts in words (summary.txt, the build's log)."""
    by = collections.Counter((c["gene"], c["call"], c["action"]) for c in calls.values())
    genes = sorted({g for g, _c, _a in by})
    lines = [f"reference screen: {', '.join(os.path.basename(f) for f in files)}"] if files else []
    lines.append(f"records {counts['records']:,}: screened {counts['screened']:,} ({counts['no_taxid']:,} of no taxid "
                 f"NCBI has, {counts['short']:,} under {MIN_LENGTH} bases, {counts['duplicate_names']:,} named twice "
                 f"in a gene's file; {counts['filed_twice']:,} names under two genes); "
                 f"references: {counts['trusted']:,} trusted, {counts['references'] - counts['trusted']:,} species' "
                 "own")
    lines.append(f"windows {counts['windows']:,}; the first pass flagged {counts['flagged']:,} records, the check "
                 f"({counts['check_windows']:,} windows) cleared {counts['cleared']:,}")
    for gene in genes:
        lines.append(f"  {gene}: " + ", ".join(
            f"{call} {sum(n for (g, c, _a), n in by.items() if g == gene and c == call):,}"
            + (" (" + ", ".join(f"{a} {by[(gene, call, a)]:,}" for a in ("drop", "trim", "review")
                                if by[(gene, call, a)]) + ")" if call != "cleared" else "")
            for call in CALLS if any(g == gene and c == call for g, c, _a in by)))
    lines.append(f"dropped {counts['action_drop']:,}, cut {counts['action_trim']:,}, for review "
                 f"{counts['action_review']:,}")
    return "\n".join(lines) + "\n"


def write_screened(path, out, calls, gene):
    """A gene's FASTA (its records' names their headers' first word, to a tab or ###) without the records its screen
    drops, those it cuts cut (" screen:start-end" after their header) -> the records written (BGZF)."""
    n, keep, span = 0, False, None
    parts, head = [], None

    def flush(fh):
        seq = "".join(parts)
        if span:
            seq = seq[span[0] - 1:span[1]]
        fh.write(head + (f" screen:{span[0]}-{span[1]}" if span else "") + "\n")
        fh.writelines(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))
    with util.open_text(path) as src, util.bgzf_open(util.part_of(out)) as fh:
        for line in src:
            if line.startswith(">"):
                if keep:
                    flush(fh)
                    n += 1
                head, parts = line.rstrip("\n"), []
                name = head[1:].split("###")[0].split("\t")[0].split(" ")[0]
                c = calls.get((gene, name))
                keep = not c or c["action"] != "drop"
                span = (int(c["keep_start"]), int(c["keep_end"])) if c and c["action"] == "trim" else None
            elif keep:
                parts.append(line.strip())
        if keep and head is not None:
            flush(fh)
            n += 1
    os.replace(util.part_of(out), out)
    return n
