"""The chimera screens of the published 18S (since 2026-11-14, format 7): KSGP 4.5's window screen and KSGP 4.2's
UCHIME check with 4.3's parts, rebuilt for Eukore's records (KSGP5.0/scripts/chimera_screen.py; KSGP 5.2's README).
The clades are EUKARYOME's (the build's tree: eukaryome.py), every record's placed in it.

Window screen (chimeras of two clades and mislabels; KSGP 4.5 found that UCHIME, which needs both parents among the
references and looks for one crossover, missed the chimeras behind most of ribotal's wrong-phylum placements). Its
"phylum" is the clade the caller gives (audit.CLADE: EUKARYOME's kingdom, the 18S's match of KSGP's 16S phylum).
  The references are the 18S of genomes (Eukore's own, ROD's, EUKARYOME's GCA_ records), each with its phylum and
  genome. Each published 18S is cut into WINDOW-base windows every STEP bases (the last ending at its end) and each is
  mapped against them (minimap2 -x asm20 --secondary=no -c). A window counts for its best hit's phylum where that hit
  covers COVER of it at TRIGGER or more; a record is flagged where its windows name two phyla, or one other than its
  label's. The check: a flagged record's windows every VERIFY_STEP bases mapped against each phylum it involves (its
  label's and its windows'), one phylum at a time, secondary hits included. A window stands for a phylum where that
  phylum's identity is IDENTITY or more and MARGIN above every other phylum's: the label's phylum with its best hit,
  another with the best hit of a second genome (one reference, a chimera or a contaminant itself, does not make a
  window another phylum's). Another phylum's piece is certain from STRONG_RUN consecutive windows or one window at
  STRONG_IDENTITY; one window below that counts where LABEL_WINDOWS windows stand for the label's phylum. Calls:
  chimera (windows stand for two phyla, or for another while one leans to its label's), mislabel (another phylum
  only, its windows over MISLABEL_SHARE of the 18S: KSGP 4.5's relabel rule), unplaced (the label's phylum a
  placeholder or no reference's), weak (one uncertain window, or another phylum's windows over less of a lineage
  that has no genome relative: a conserved stretch), cleared.
  A chimera is cut to its longest stretch outside the windows of another phylum where MIN_PART bases or more stay
  (an 18S record; flagged `chimera_cut`, medium), unless RECURRENT or more calls share that stretch's start (within
  100 bases) and reference (a trait of the lineage: `chimera_recurrent`, medium); else flagged `chimera` (low), as
  an operon read always is (its 18S another organism's in part). A mislabel: `mislabel` (low); weak: `chimera_weak`
  (medium).
UCHIME (chimeras of two genera; KSGP 4.2's class Q, 4.3's parts)
  Every distinct 18S is screened against all the others (vsearch --uchime_ref; the queries of a part left out of its
  database, a part's queries every n-th distinct 18S: KSGP's ten folds, as many as the nodes take). A call is
  confirmed where each part (the query cut at UCHIME's crossover) is PART_ID or more like a record and their best
  records name no genus in common, and the whole query is WHOLE_MARGIN points or more below its weaker part (KSGP
  4.2's rule); where the parts' genera meet or are unnamed, `chimera_weak` (medium). A confirmed chimera is cut to
  the part of its own genus (else family) where that holds MIN_PART bases or more (`chimera_cut`), else flagged
  `chimera` (low).
Distributed: the window screen's first pass and the UCHIME parts are searches of qc's protocol (qc.MapParams,
qc.UchimeParams): kept in WORK/qc/<phase>/, taken by the helper nodes (--qc_helpers) where there are HELPED parts or
more; the checks that follow (a few thousand records) run on the final node.
"""

import collections
import os
import subprocess
import tempfile
import threading
import zlib

from . import util

WINDOW, STEP, VERIFY_STEP = 300, 150, 50
TRIGGER, IDENTITY, COVER, MARGIN = 0.95, 0.97, 0.75, 0.02
STRONG_RUN, STRONG_IDENTITY, LABEL_WINDOWS = 2, 0.99, 2
SECONDARY = 30      # the check's secondary alignments a window (past a genome's copies to a second genome)
MIN_LENGTH = 300    # a shorter 18S is not screened
MIN_PART = 500      # a chimera's own part kept from this long (KSGP 4.3's rule)
RECURRENT = 5
MISLABEL_SHARE = 0.5  # a mislabel: another clade's windows over this share of the 18S (KSGP 4.5's relabel rule)
PART_ID = 95.0
WHOLE_MARGIN = 2.0
PRESETS = {"asm20": ["-x", "asm20"], "sr": ["-x", "sr"]}
AUDIT_COLUMNS = ("record", "check", "call", "action", "evidence")


# --- mapping ---------------------------------------------------------------------------------------------------------

def window_starts(length, size=WINDOW, step=STEP):
    if length <= size:
        return [0]
    starts = list(range(0, length - size + 1, step))
    if starts[-1] + size < length:
        starts.append(length - size)
    return starts


def windows_of(seqs, step=STEP, size=WINDOW):
    """{name: sequence} -> (window name `name:start`, its bases)."""
    for name, seq in seqs.items():
        for s in window_starts(len(seq), size, step):
            yield f"{name}:{s}", seq[s:s + size]


def minimap2_hits(minimap2, threads, reference, windows, secondary=0, preset="asm20", log=None):
    """Windows ((name, sequence), fed through a pipe) mapped against reference (a FASTA file) -> (name, {target:
    (matches, block, share of the window covered)}) of each window's hits: its primary alignment and up to `secondary`
    secondary ones, the pieces of a window aligned to one target summed (KSGP's run_minimap2)."""
    cmd = [minimap2] + PRESETS[preset] + ["-c", "-t", str(max(1, threads))]
    cmd += ["--secondary=yes", "-N", str(secondary), "-p", "0.8"] if secondary else ["--secondary=no"]
    cmd += [reference, "-"]
    with open(log or os.devnull, "a") as err:
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err, text=True,
                                bufsize=1 << 20)

        def feed():
            try:
                for name, seq in windows:
                    proc.stdin.write(f">{name}\n{seq}\n")
            finally:
                proc.stdin.close()
        feeder = threading.Thread(target=feed, daemon=True)
        feeder.start()
        name, hits = None, {}
        for line in proc.stdout:
            f = line.split("\t", 12)
            if f[0] != name:
                if name is not None:
                    yield name, hits
                name, hits = f[0], {}
            matches, block, covered = hits.get(f[5], (0, 0, 0.0))
            hits[f[5]] = (matches + int(f[9]), block + int(f[10]),
                          min(1.0, covered + (int(f[3]) - int(f[2])) / max(1, int(f[1]))))
        if name is not None:
            yield name, hits
        feeder.join()
        proc.stdout.close()
        if proc.wait():
            raise RuntimeError(f"minimap2 failed (exit {proc.returncode}): {' '.join(cmd)}")


def read_fasta(path, wanted=None):
    out, name, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    out[name] = "".join(parts)
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        out[name] = "".join(parts)
    return out


def map_file(minimap2, queries, reference, out, threads=1, preset="asm20", step=STEP):
    """A window-screen part (qc.MapParams): the queries' (a FASTA file) windows every `step` mapped against
    reference -> out: a line a window with a hit, `window <tab> target <tab> matches <tab> block <tab> covered` of its
    best hit."""
    seqs = read_fasta(queries)
    with open(out, "w", newline="\n") as fh:
        for name, hits in minimap2_hits(minimap2, threads, reference, windows_of(seqs, step), 0, preset):
            target, (m, b, c) = max(hits.items(), key=lambda kv: (kv[1][0], kv[0]))
            fh.write(f"{name}\t{target}\t{m}\t{b}\t{c:.3f}\n")


def map_hits(fh):
    """map_file's lines -> [(window, target, matches, block, covered)]."""
    out = []
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) == 5:
            out.append((f[0], f[1], int(f[2]), int(f[3]), float(f[4])))
    return out


# --- the window screen -----------------------------------------------------------------------------------------------

def _first_pass(queries, refs, minimap2, threads, tmp, searcher=None, log=None):
    """Every query's windows against every reference -> {query: {phyla its windows count for}}."""
    standing = collections.defaultdict(set)
    phylum = {r: p for r, (_s, p, _g) in refs.items()}
    rows = None
    if searcher is not None and getattr(searcher, "folder", None):
        rows = searcher.mapped("window_screen", {q: s for q, (s, _p) in queries.items()},
                               {r: s for r, (s, _p, _g) in refs.items()}, minimap2)
    if rows is None:
        with tempfile.TemporaryDirectory(dir=tmp) as d:
            ref = os.path.join(d, "reference.fna")
            with open(ref, "w") as fh:
                fh.writelines(f">{r}\n{s}\n" for r, (s, _p, _g) in refs.items())
            rows = []
            for name, hits in minimap2_hits(minimap2, threads, ref, windows_of({q: s for q, (s, _p) in
                                                                                 queries.items()})):
                target, (m, b, c) = max(hits.items(), key=lambda kv: (kv[1][0], kv[0]))
                rows.append((name, target, m, b, c))
    for window, target, m, b, c in rows:
        q = window.rsplit(":", 1)[0]
        if b and m / b >= TRIGGER and c >= COVER and q in queries:
            standing[q].add(phylum.get(target))
    if log:
        log(f"  window screen, first pass: {len(queries):,} 18S, {len(rows):,} windows with a hit")
    return standing


def _scores(queries, flagged, refs, minimap2, threads, tmp):
    """The check: each flagged query's windows every VERIFY_STEP against each phylum it involves -> {(query, start):
    {phylum: (identity, target)}}: the label's phylum its best hit, another phylum a second genome's best."""
    involved = collections.defaultdict(set)
    ref_phyla = {p for _s, p, _g in refs.values()}
    for q, phyla in flagged.items():
        for p in phyla | {queries[q][1]}:
            if p in ref_phyla:
                involved[p].add(q)
    by_phylum = collections.defaultdict(dict)
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        ref = os.path.join(d, "phylum.fna")
        for p, qs in sorted(involved.items(), key=lambda kv: str(kv[0])):
            with open(ref, "w") as fh:
                fh.writelines(f">{r}\n{s}\n" for r, (s, rp, _g) in refs.items() if rp == p)
            wins = windows_of({q: queries[q][0] for q in sorted(qs)}, VERIFY_STEP)
            for name, hits in minimap2_hits(minimap2, threads, ref, wins, SECONDARY):
                q, start = name.rsplit(":", 1)
                per_genome = {}
                for target, (m, b, c) in hits.items():
                    ident = m / b if b and c >= COVER else 0.0
                    g = refs[target][2]
                    if g not in per_genome or (ident, target) > per_genome[g]:
                        per_genome[g] = (ident, target)
                ranked = sorted(per_genome.values(), reverse=True)
                if p != queries[q][1]:
                    ranked = ranked[1:2]  # a second genome's best
                by_phylum[(q, int(start))][p] = ranked[0] if ranked else (0.0, "")
    return by_phylum, ref_phyla


def _call(q, seq, own, phyla, by_phylum, ref_phyla, refs):
    """One flagged query's call -> (call, foreign intervals [(start, end)], standing counts, evidence)."""
    stands, leans = collections.Counter(), 0
    run, longest, best_foreign, foreign = 0, 0, 0.0, []
    targets = set()
    for s in window_starts(len(seq), WINDOW, VERIFY_STEP):
        scores = dict(by_phylum.get((q, s), {}))
        for p in phyla | {own}:
            if p in ref_phyla:
                scores.setdefault(p, (0.0, ""))
        ranked = sorted(scores.items(), key=lambda kv: (-kv[1][0], str(kv[0])))
        if not ranked:
            continue
        (top, (ident, target)), second = ranked[0], ranked[1][1][0] if len(ranked) > 1 else 0.0
        stand = top if ident >= IDENTITY and ident - second >= MARGIN else None
        lean = top == own and ident - second >= MARGIN
        if stand:
            stands[stand] += 1
        if stand and stand != own:
            run += 1
            longest, best_foreign = max(longest, run), max(best_foreign, ident)
            foreign.append((s, s + WINDOW))
            targets.add(target)
        else:
            run = 0
        leans += lean
    others = {p for p in stands if p != own}
    strong = longest >= STRONG_RUN or best_foreign >= STRONG_IDENTITY
    if others and not strong and stands[own] < LABEL_WINDOWS:
        call = "weak"
    elif own not in ref_phyla:
        call = "chimera" if len(stands) >= 2 else "unplaced" if stands else "cleared"
    elif len(stands) >= 2 or (others and leans):
        call = "chimera"
    elif others:
        # another clade's windows over half of it (KSGP 4.5's relabel rule), else a conserved stretch of a lineage
        # with no genome relative of its own (EUKARYOME 2.1: Sarcomonadea reads with two windows at 98% to
        # stramenopile genomes; 2026-11-14)
        covered = sum(b - a for a, b in merge(foreign))
        call = "mislabel" if covered >= MISLABEL_SHARE * len(seq) else "weak"
    else:
        call = "cleared"
    evidence = (",".join(f"{p}:{n}" for p, n in stands.most_common()) + f"; run {longest}, best foreign "
                f"{best_foreign:.3f}; targets " + ",".join(sorted(targets)[:3]))
    return call, merge(foreign), sorted(targets), evidence


def merge(intervals):
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


def own_stretch(length, foreign):
    """The longest stretch of a query outside its foreign intervals -> (start, end)."""
    edges, at, best = [], 0, (0, 0)
    for a, b in foreign:
        edges.append((at, a))
        at = max(at, b)
    edges.append((at, length))
    for a, b in edges:
        if b - a > best[1] - best[0]:
            best = (a, b)
    return best


def window_screen(queries, refs, minimap2, threads=1, tmp=None, searcher=None, log=None):
    """queries {name: (18S, its label's phylum or None)}, refs {name: (18S, phylum, genome)} -> {name: (call, keep
    (start, end) on its 18S or None, evidence)} of the calls but cleared."""
    queries = {q: v for q, v in queries.items() if len(v[0]) >= MIN_LENGTH}
    refs = {r: v for r, v in refs.items() if v[1]}
    if not queries or not refs or not minimap2:
        return {}
    standing = _first_pass(queries, refs, minimap2, threads, tmp, searcher, log)
    flagged = {q: {p for p in phyla if p} for q, phyla in standing.items()
               if len(phyla) >= 2 or phyla != {queries[q][1]}}
    if log:
        log(f"  window screen: {len(flagged):,} flagged by the first pass")
    by_phylum, ref_phyla = _scores(queries, flagged, refs, minimap2, threads, tmp)
    calls = {}
    for q, phyla in flagged.items():
        seq, own = queries[q]
        call, foreign, targets, evidence = _call(q, seq, own, phyla, by_phylum, ref_phyla, refs)
        if call == "cleared":
            continue
        keep = None
        if call == "chimera":
            a, b = own_stretch(len(seq), foreign)
            keep = (a, b, targets[0] if targets else "") if b - a >= MIN_PART else None
        calls[q] = (call, keep, evidence)
    # recurrent: RECURRENT or more chimeras kept to a stretch of the same start (within 100 bases) and reference
    shared = collections.Counter((k[0] // 100, k[2]) for c, k, _e in calls.values() if c == "chimera" and k)
    for q, (c, k, e) in list(calls.items()):
        if c == "chimera" and k and shared[(k[0] // 100, k[2])] >= RECURRENT:
            calls[q] = ("recurrent", None, e)
        elif k:
            calls[q] = (c, k[:2], e)
    if log:
        log("  window screen: " + ", ".join(f"{k} {v:,}" for k, v in sorted(collections.Counter(
            c for c, _k, _e in calls.values()).items())))
    return calls


# --- UCHIME ----------------------------------------------------------------------------------------------------------

def fold_of(seq, folds):
    return zlib.crc32(seq.encode()) % max(1, folds)


def uchime_file(vsearch, queries, database, out, threads=1, tmp=None):
    """A UCHIME part (qc.UchimeParams): the queries (a FASTA file) against the database less themselves (written to
    tmp) -> out: a line a chimera, `query <tab> parent A <tab> parent B <tab> score <tab> crossover (bases of the
    query before it) <tab> which parent the start holds (A or B)`."""
    names = set(read_fasta(queries))
    with tempfile.TemporaryDirectory(dir=tmp) as d:
        db = os.path.join(d, "db.fna")
        with open(database) as src, open(db, "w") as dst:
            keep = True
            for line in src:
                if line.startswith(">"):
                    keep = line[1:].strip() not in names
                if keep:
                    dst.write(line)
        table, alns = os.path.join(d, "uchime.tsv"), os.path.join(d, "alns.txt")
        subprocess.run([vsearch, "--uchime_ref", queries, "--db", db, "--uchimeout", table, "--uchimealns", alns,
                        "--threads", str(max(1, threads)), "--quiet"], check=True)
        crossings = crossovers(alns)
        with open(table) as fh, open(out, "w", newline="\n") as dst:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 18 and f[17] == "Y":
                    at, first = crossings.get(f[1], (None, None))
                    dst.write(f"{f[1]}\t{f[2]}\t{f[3]}\t{f[0]}\t{at if at is not None else '-'}\t{first or '-'}\n")


def crossovers(path):
    """vsearch's --uchimealns -> {query: (its bases before the crossover, the parent its start holds)}: where the
    Model line turns from one parent to the other (the middle of the `x` stretch between them)."""
    out = {}
    query, q_line, model = None, [], []

    def done():
        if query is None or not model:
            return
        m, q = "".join(model), "".join(q_line)
        letters = [(i, c) for i, c in enumerate(m) if c in "AB"]
        if not letters:
            return
        first = letters[0][1]
        turn = next((i for i, c in letters if c != first), None)
        if turn is None:
            return
        last_first = max(i for i, c in letters if c == first and i < turn)
        col = (last_first + turn) // 2 + 1
        out[query] = (sum(1 for c in q[:col] if c != "-"), first)
    with open(path) as fh:
        for line in fh:
            if line.startswith("Query "):
                done()
                query, q_line, model = line.split(")", 1)[1].strip(), [], []
            elif line.startswith("Q ") and query is not None:
                f = line.split()
                if len(f) >= 4:
                    q_line.append(f[2])
            elif line.startswith("Model") and query is not None:
                model.append(line[8:].rstrip("\n"))
    done()
    return out


def uchime_hits(fh):
    out = []
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) == 6:
            out.append((f[0], f[1], f[2], float(f[3]), None if f[4] == "-" else int(f[4]), f[5]))
    return out


def uchime_screen(screened, others, genera, searcher, vsearch, threads=1, tmp=None, log=None):
    """screened {name: 18S} (the published records'), others {name: 18S} (the genomes': references only), genera
    {name: (genus, family) of its EUKARYOME lineage, None where unnamed} -> {name: (call, keep (start, end) or None,
    evidence)}: confirmed chimeras (`chimera`) and calls whose parts' genera meet or are unnamed (`weak`)."""
    if not vsearch or not screened:
        return {}
    distinct, members = {}, collections.defaultdict(list)
    for pool in (screened, others):
        for name, seq in pool.items():
            if len(seq) < MIN_LENGTH:
                continue
            rep = distinct.setdefault(seq, f"u{len(distinct)}")
            members[rep].append(name)
    seq_of = {rep: s for s, rep in distinct.items()}
    to_screen = sorted({distinct[s] for n, s in screened.items() if s in distinct}, key=lambda r: int(r[1:]))
    rows = searcher.uchimed("uchime", {r: seq_of[r] for r in to_screen}, seq_of, vsearch) if searcher is not None \
        and getattr(searcher, "folder", None) else None
    if rows is None:
        with tempfile.TemporaryDirectory(dir=tmp) as d:
            q, db, out = os.path.join(d, "q.fna"), os.path.join(d, "db.fna"), os.path.join(d, "out.tsv")
            with open(q, "w") as fh:
                fh.writelines(f">{r}\n{seq_of[r]}\n" for r in to_screen)
            with open(db, "w") as fh:
                fh.writelines(f">{r}\n{s}\n" for r, s in seq_of.items())
            uchime_file(vsearch, q, db, out, threads, d)
            with open(out) as fh:
                rows = uchime_hits(fh)
    if log:
        log(f"  UCHIME: {len(to_screen):,} distinct 18S screened against {len(seq_of):,}, {len(rows):,} called")
    # the check: each part and the whole searched against every distinct 18S
    parts, called = {}, {}
    for rep, a, b, score, at, first in rows:
        s = seq_of[rep]
        if at is None or at <= 0 or at >= len(s):
            continue
        called[rep] = (a, b, score, at, first)
        parts[f"{rep}|L"], parts[f"{rep}|R"], parts[f"{rep}|W"] = s[:at], s[at:], s
    found = searcher.one("uchime_parts", parts, seq_of, 0.8, 32, 64) if parts and searcher is not None else {}

    def best(key):
        own = key.split("|")[0]  # the query's own sequence among the distinct ones: not a hit
        hits = sorted(((i, t) for i, t in found.get(key, ()) if t != own), reverse=True)
        if not hits:
            return 0.0, set(), set()
        top = hits[0][0]
        near = [t for i, t in hits if i >= top - 0.5]
        gen = {genera.get(m, (None, None))[0] for t in near for m in members.get(t, ())}
        fam = {genera.get(m, (None, None))[1] for t in near for m in members.get(t, ())}
        return top, gen, fam
    calls = {}
    for rep, (a, b, score, at, first) in called.items():
        (li, lg, lf), (ri, rg, rf), wi = best(f"{rep}|L"), best(f"{rep}|R"), best(f"{rep}|W")[0]
        if min(li, ri) < PART_ID or wi > min(li, ri) - WHOLE_MARGIN:
            continue  # not confirmed: a part unlike every record, or the whole as like one
        evidence = (f"UCHIME {score:.2f} ({a}|{b}); crossover {at}; parts {li:.1f}/{ri:.1f}, whole {wi:.1f}; genera "
                    + ",".join(sorted(g for g in lg if g)) + " | " + ",".join(sorted(g for g in rg if g)))
        # the named genera of each part's best records (an unnamed one, EUKARYOME's environmental reads', tells
        # nothing; until 2026-11-14 one made the call uncertain: 682 of 1,006 KSGP 5.2 chimeras sampled)
        lg, rg = lg - {None}, rg - {None}
        for m in members[rep]:
            if m not in screened:
                continue
            if not lg or not rg or lg & rg:
                calls[m] = ("weak", None, evidence)
                continue
            genus, family = genera.get(m, (None, None))
            own = [(0, at)] if genus and genus in lg else [(at, len(seq_of[rep]))] if genus and genus in rg else \
                [(0, at)] if family and family in lf and family not in rf else \
                [(at, len(seq_of[rep]))] if family and family in rf and family not in lf else []
            keep = own[0] if own and own[0][1] - own[0][0] >= MIN_PART else None
            calls[m] = ("chimera", keep, evidence)
    if log:
        log("  UCHIME, confirmed: " + ", ".join(f"{k} {v:,}" for k, v in sorted(collections.Counter(
            ("chimera, own part kept" if k else c) for c, k, _e in calls.values()).items())))
    return calls


def write_table(path, rows):
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(AUDIT_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "-")) for c in AUDIT_COLUMNS) + "\n")
    os.replace(util.part_of(path), path)
