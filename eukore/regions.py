"""The 18S variable regions of each rRNA record (since 2026-11-07; rrna/regions.tsv): V1-V5 and V7-V9 as Hadziavdic et
al. 2014 place them on Saccharomyces cerevisiae's 18S (data/rrna_regions.tsv), carried to every record through
Rfam's eukaryote SSU model (RF01960, SSU_rRNA_eukarya): cmalign aligns the yeast 18S once, which gives each region's
first and last model positions, then every record's 18S; a region is the record's bases from the first at its first
position or after to the last at its last position or before (inserts between included). It is complete where the
record holds bases within SLACK positions of both its ends. The alignments are cached by sequence under WORK/cmalign/.
"""

import collections
import hashlib
import os
import shutil
import subprocess
import tempfile
import time

from . import align as align_mod
from . import util

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REGIONS = os.path.join(HERE, "data", "rrna_regions.tsv")
REFERENCE = os.path.join(HERE, "data", "primer_references.fna")
REFERENCE_NAME = "Sc_18S"
SSU_CM = ("RF01960", "https://rfam.org/family/RF01960/cm")  # Rfam's eukaryote SSU model
SLACK = 3
BATCH = 500  # sequences per cmalign run (an alignment's width grows with its inserts)
REGION_COLUMNS = ("record", "region", "start", "end", "length", "complete")


# Rfam's SSU models an 18S is scored by (since 2026-11-09, ssu_models): barrnap's eukaryote models alone call a
# bacterium's or a mitochondrion's SSU an 18S, and an assembly's bacterial contig, a NUMT, an endosymbiont's insertion
# then gave the species a "18S" (and the universal motifs rrna.py seeds by find them). The model of the best score names
# it; a prokaryote model's score must pass the eukaryote models' best by SSU_MARGIN bits.
SSU_MODELS = (("RF01960", "eukarya"), ("RF02542", "eukarya"),  # SSU_rRNA_eukarya, SSU_rRNA_microsporidia
              ("RF00177", "bacteria"), ("RF01959", "archaea"))  # SSU_rRNA_bacteria (a mitochondrion's, a plastid's)
SSU_MARGIN = 20.0
SSU_MODELS_FILE = "ssu_models.cm"
# the SSU models' search in kept parts (since 2026-11-13; until then one cmsearch of every new 18S, an hour or more on
# 60 cores for a release's, kept only once it had ended): SSU_PER_NODE parts for each node searching (the final job and
# its helpers), more where a part would have over SSU_PART_MAX sequences (each part kept as it ends), none under
# SSU_PART_MIN; the helpers started for SSU_HELPED parts or more
SSU_PER_NODE = 3
SSU_PART_MAX = 20000
SSU_PART_MIN = 500
SSU_HELPED = 4


def ssu_parts(names, nodes=1):
    """names (sorted) -> the SSU models' search's parts [[name]], interleaved (part i every n-th name from the i-th)."""
    n = max(SSU_PER_NODE * max(1, nodes), -(-len(names) // SSU_PART_MAX))
    n = max(1, min(n, len(names) // SSU_PART_MIN))
    return [names[i::n] for i in range(n)] if names else []


def cm_path(src_dir):
    return os.path.join(src_dir, "rfam", f"{SSU_CM[0]}.cm")


def models_path(src_dir):
    """The SSU models' file (SSU_MODELS, one Infernal file; ssu_models)."""
    return os.path.join(src_dir, "rfam", SSU_MODELS_FILE)


def fetch(src_dir, download, log=print):
    """Rfam's SSU model into src_dir/rfam (kept when there and an Infernal model: a page saved in its place, e.g. a
    proxy's, is fetched again; since 2026-11-08 it was kept, and every package stage then stopped), and (since
    2026-11-09) SSU_MODELS' joined into models_path's file: written where it is not those models' files one after the
    other (since 2026-11-10; until then kept where it was an Infernal model, so that a model added to SSU_MODELS, or
    one fetched again, never reached cmsearch)."""
    path = cm_path(src_dir)
    folder = os.path.dirname(path)
    os.makedirs(folder, exist_ok=True)
    for acc, _kind in SSU_MODELS:
        one = os.path.join(folder, f"{acc}.cm")
        if not is_cm(one):
            log(f"Rfam: fetching {acc}")
            url = f"https://rfam.org/family/{acc}/cm"
            download(url, one)
            if not is_cm(one):
                os.remove(one)
                raise OSError(f"{url}: not an Infernal model")
    joined = models_path(src_dir)
    text = ""
    for acc, _kind in SSU_MODELS:
        with open(os.path.join(folder, f"{acc}.cm"), newline="") as fh:
            text += fh.read()
    try:
        with open(joined, newline="") as fh:
            same = fh.read() == text
    except OSError:
        same = False
    if not same:
        with open(util.part_of(joined), "w", newline="") as out:
            out.write(text)
        os.replace(util.part_of(joined), joined)
    return path


def _digest(seq):
    return hashlib.md5(seq.encode()).hexdigest()


def ssu_models(seqs, cm, cmsearch="cmsearch", threads=1, tmp=None, work=None, searches=None, progress=None,
               log=None):
    """Each sequence's SSU model of the best score ({name: its 18S} -> {name: eukarya, bacteria, archaea, or - where
    none scores it}; since 2026-11-09): cmsearch of SSU_MODELS (cm: models_path's file) by their profile HMMs alone
    (--hmmonly: an HMM search's speed), each distinct sequence once; a prokaryote model's call needs its score over
    the eukaryote models' best by SSU_MARGIN bits.
    Since 2026-11-10 a sequence's model is kept under work/cmalign/models.<digest of the models' file, SSU_MODELS and
    SSU_MARGIN> (work: the build's WORK; tmp where none is given, which the build gives its WORK as), a line a
    sequence as run's cache: until then every package stage searched every 18S again. -Z 1 makes a sequence's hits
    its own whatever is searched with it, so that the kept call is the one a search would give.
    searches (the build's qc_searches: folder, helpers, nodes, scratch; since 2026-11-13): the sequences not in the
    cache searched in parts (ssu_parts) kept in its folder as each ends (a killed run's read again), shared with the
    helpers as the vsearch checks' are (qc.searched_loci, qc.CmParams); the folder removed once the cache holds their
    calls. progress: the parts' (util.Progress); log: each part's time and the search's."""
    kind_of = dict(SSU_MODELS)
    distinct = {}
    for name, seq in seqs.items():
        if seq:
            distinct.setdefault(seq.upper(), []).append(name)
    if not distinct:
        return {}
    folder = work or tmp
    stem = None
    digest = util.md5(cm) if os.path.isfile(cm) else None
    if folder and digest:
        key = hashlib.md5("\n".join([digest, repr(SSU_MODELS), str(SSU_MARGIN)]).encode())
        stem = os.path.join(folder, "cmalign", f"models.{key.hexdigest()[:10]}")
    cache = util.cache_read(stem) if stem else {}
    order = {_digest(seq): seq for seq in distinct if _digest(seq) not in cache}
    if order:
        began = time.time()
        kept = (searches or {}).get("folder") if digest else None
        if kept:
            from . import qc  # noqa: PLC0415
            todo = [(f"ssu-{i + 1}", {d: order[d] for d in part}, None)
                    for i, part in enumerate(ssu_parts(sorted(order), (searches or {}).get("nodes", 1)))]
            found = qc.searched_loci(todo, cmsearch, threads, kept,
                                     searches.get("helpers") if len(todo) >= SSU_HELPED else None, progress,
                                     params=qc.CmParams(os.path.abspath(cm), digest),
                                     what="SSU models searched (cmsearch)", scratch=searches.get("scratch") or tmp,
                                     log=log)
            hits = [h for part in found.values() for h in part]
        else:
            with tempfile.TemporaryDirectory(dir=tmp) as scratch:
                fasta, table = os.path.join(scratch, "s.fna"), os.path.join(scratch, "hits.tbl")
                with open(fasta, "w") as fh:
                    fh.writelines(f">{d}\n{seq}\n" for d, seq in order.items())
                util.cmsearch_files(cmsearch, fasta, cm, table, threads)  # failed: nothing kept
                with open(table) as fh:
                    hits = util.cmsearch_hits(fh)
        best = collections.defaultdict(dict)  # sequence's digest -> {kind: best score}
        for d, accession, score in hits:
            kind = kind_of.get(accession)
            if kind and score > best[d].get(kind, float("-inf")):
                best[d][kind] = score
        new = {}
        for d in order:
            scores = best.get(d)
            if not scores:
                kind = "-"
            else:
                kind = max(scores, key=scores.get)
                if kind != "eukarya" and scores[kind] < scores.get("eukarya", float("-inf")) + SSU_MARGIN:
                    kind = "eukarya"
            new[d] = kind
        cache.update(new)
        if stem:
            with util.cache_open(stem) as out:
                out.writelines(f"{k}\t{v}\n" for k, v in sorted(new.items()))
            if kept:  # the parts' calls in the cache now
                shutil.rmtree(kept, ignore_errors=True)
        if log:
            log(f"  SSU models: {len(order):,} 18S searched (cmsearch"
                f"{f', {len(todo)} parts' if kept else ''}): took {util.duration(time.time() - began)}")
    return {name: cache[_digest(seq)] for seq, names in distinct.items() for name in names}


def is_cm(path):
    """Whether a file is an Infernal model (its first line INFERNAL1/a ...)."""
    try:
        with open(path, "rb") as fh:
            return fh.read(9) == b"INFERNAL1"
    except OSError:
        return False


def read_regions(path=REGIONS):
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    return [(r[0], int(r[1]), int(r[2])) for r in lines[1:]]


def reference(path=REFERENCE, name=REFERENCE_NAME):
    with open(path) as fh:
        lines = fh.read().split(">")
    return next(chunk.split("\n", 1)[1].replace("\n", "") for chunk in lines if chunk.startswith(name + "\n"))


MAX_SSU = 4000   # an 18S longer than this (a misjoined unit, an intron-laden one) is not aligned
MXSIZE = 2048    # cmalign's matrix limit (Mb, per thread); a sequence needing more is left out


def _program(tool):
    """A program's identity in a cache's key: its file, its size and its time (its name where it is not found)."""
    path = shutil.which(tool) or tool
    try:
        found = os.stat(path)
    except OSError:
        return tool
    return f"{os.path.realpath(path)}\t{found.st_size}\t{found.st_mtime_ns}"


def _cmalign(cm, seqs, names, cmalign_bin, threads, tmp):
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        fasta, got = os.path.join(work, "s.fna"), os.path.join(work, "a.pfam")
        with open(fasta, "w") as fh:
            fh.writelines(f">{n}\n{seqs[n]}\n" for n in names)
        run = subprocess.run([cmalign_bin, "--cpu", str(max(1, threads)), "--outformat", "Pfam", "--noprob",
                              "--mxsize", str(MXSIZE), "-o", got, cm, fasta], capture_output=True, text=True,
                             check=False)
        if run.returncode != 0:  # its own exit (1: a sequence over its matrix limit), or killed (a signal's, a shell's)
            return None, run.stderr.strip()[-300:], 0 < run.returncode < 128
        rows, rf = align_mod.read_pfam(got)
        return {n: align_mod.states_of(row, rf) for n, row in rows.items()}, None, True


def cmalign(cm, seqs, cmalign_bin="cmalign", threads=1, tmp=None, log=None, each=None, fail=None):
    """seqs {name: sequence} aligned to the CM -> {name: [the model position of each base, or None in an insert]}. A
    sequence over MAX_SSU is left out; a batch cmalign fails on (one sequence over its matrix limit) is aligned in
    halves, a half that fails in its halves, down to the sequences that fail alone: left out (said). Since 2026-11-10;
    until then one sequence at a time, BATCH runs for one sequence that fails, where the halves take ~2 log2(BATCH).
    each: a function given each run's alignments as they come (run's cache); fail: one given {name: cmalign's words}
    of the sequences that failed alone by cmalign's own exit (one killed says nothing of its sequence)."""
    out, failed, own = {}, [], {}
    names = sorted(n for n in seqs if len(seqs[n]) <= MAX_SSU)
    for i in range(0, len(names), BATCH):
        parts = [names[i:i + BATCH]]
        while parts:
            part = parts.pop()
            got, why, ended = _cmalign(cm, seqs, part, cmalign_bin, threads, tmp)
            if got is not None:
                out.update(got)
                if each:
                    each(got)
            elif len(part) > 1:
                half = len(part) // 2
                parts += [part[half:], part[:half]]  # the first half next
            else:
                failed.append((part[0], why))
                if ended:
                    own[part[0]] = why
    if fail and own:
        fail(own)
    long_ones = len(seqs) - len(names)
    if log and (failed or long_ones):
        log(f"  cmalign: {long_ones} 18S over {MAX_SSU} bases and {len(failed)} it could not align left out"
            + (f" (e.g. {failed[0][1][-120:]})" if failed else ""))
    return out


def columns(states, regions):
    """The reference's states and the regions on it (1-based, end excluded) -> {region: (first, last model position)}."""
    out = {}
    for name, start, end in regions:
        inside = [s for s in states[start - 1:end - 1] if s is not None]
        if inside:
            out[name] = (min(inside), max(inside))
    return out


def locate(states, cols):
    """A record's 18S states and the regions' model positions -> {region: (start, end, complete)} on the 18S (0-based,
    end excluded)."""
    out = {}
    placed = [(i, s) for i, s in enumerate(states) if s is not None]
    if not placed:
        return out
    low, high = min(s for _i, s in placed), max(s for _i, s in placed)
    for name, (first, last) in cols.items():
        inside = [i for i, s in placed if first <= s <= last]
        if inside:
            out[name] = (min(inside), max(inside) + 1, low <= first + SLACK and high >= last - SLACK)
    return out


def ssu_of(seq, row):
    """A record's 18S and where it starts on the record, or None (no 18S found in it)."""
    if row.get("kind") == "18S":
        return seq, 0
    for n, a, b in row.get("regions_list") or ():
        if n == "18S":
            return seq[a:b], a
    return None


def run(records, cm, cmalign_bin="cmalign", threads=1, work=None, log=print, tmp=None):
    """Every record's 18S variable regions (records [(name, sequence, row)]) -> rows (REGION_COLUMNS). Each distinct 18S
    is aligned once (cached under WORK/cmalign/; one cmalign fails on is not aligned again by that cmalign). tmp:
    cmalign's files (since 2026-11-13; until then in work, the WORK on a shared file system)."""
    tmp = tmp or work
    regions = read_regions()
    ref = cmalign(cm, {"reference": reference()}, cmalign_bin, threads, tmp).get("reference")
    if ref is None:  # a damaged model, a cmalign that fails: no regions, the package goes on (until 2026-11-08 stopped)
        log(f"WARNING: cmalign could not align the yeast 18S reference to {cm}: no 18S variable regions")
        return []
    cols = columns(ref, regions)
    # the cache: a line per distinct 18S, added after each cmalign batch (since 2026-11-08: written once all were
    # aligned, so that a killed run kept none; and an 18S cmalign failed on was cached as one without regions, never
    # aligned again: now left out, tried again by the next run). Its name is the digest of all its regions depend on:
    # the model, the regions table, the reference, SLACK, MAX_SSU (since 2026-11-08: the model's alone, an edited
    # regions table's regions never placed); one file per process (since 2026-11-08: two package stages on one WORK
    # appended to one, their lines interleaved), all read
    key = hashlib.md5("\n".join([util.md5(cm), util.md5(REGIONS), reference(), str(SLACK), str(MAX_SSU)]).encode())
    stem = os.path.join(work, "cmalign", f"ssu.{key.hexdigest()[:10]}") if work else None
    cache = util.cache_read(stem) if stem else {}  # a line a killed run cut: aligned again
    # the 18S this cmalign failed on (since 2026-11-10, a cache of their own: every run aligned them again, each
    # failing batch in its halves, to fail again): kept under the digest of the model, the matrix limit and the
    # program (its file, size and time: another cmalign tries them again), and only where cmalign said so itself of
    # the sequence alone and aligned the reference again after it (a full disk, a node out of memory fail it too)
    failed_stem = None
    if work:
        failed_key = hashlib.md5("\n".join([util.md5(cm), str(MXSIZE), _program(cmalign_bin)]).encode())
        failed_stem = os.path.join(work, "cmalign", f"unaligned.{failed_key.hexdigest()[:10]}")
    unaligned = util.cache_read(failed_stem) if failed_stem else {}
    ssu = {}
    for name, seq, row in records:
        got = ssu_of(seq, row)
        if got:
            ssu[name] = (got[0].upper(), got[1], _digest(got[0].upper()))
    todo = {key: s for s, _o, key in ssu.values() if key not in cache}
    known = sorted(key for key in todo if key in unaligned)
    if known:
        log(f"18S variable regions: {len(known)} distinct 18S cmalign failed on before left out, e.g. "
            f"{unaligned[known[0]][-120:]} ({failed_stem}.*.tsv: removed, they are tried again)")
        todo = {key: s for key, s in todo.items() if key not in unaligned}
    if todo:
        log(f"18S variable regions: {len(todo)} distinct 18S to align (cmalign)")
        out = util.cache_open(stem) if stem else None
        failures = {}

        def keep(aligned):
            got = {key: ",".join(f"{n}:{a}-{b}:{'t' if c else 'f'}" for n, (a, b, c) in
                                 locate(states, cols).items()) or "-" for key, states in aligned.items()}
            cache.update(got)
            if out:
                out.writelines(f"{k}\t{v}\n" for k, v in sorted(got.items()))
                out.flush()
        try:
            keep({key: [] for key, s in todo.items() if len(s) > MAX_SSU})  # never aligned: none, for good
            cmalign(cm, todo, cmalign_bin, threads, tmp, log, each=keep, fail=failures.update)
        finally:
            if out:
                out.close()
        if failures and failed_stem and cmalign(cm, {"reference": reference()}, cmalign_bin, threads, tmp):
            with util.cache_open(failed_stem) as fh:
                fh.writelines(f"{k}\t{' '.join(why.split()) or '-'}\n" for k, why in sorted(failures.items()))
    out = []
    for name, (_s, offset, key) in sorted(ssu.items()):
        for part in cache.get(key, "-").split(","):
            if part == "-":
                continue
            region, span, complete = part.split(":")
            a, b = map(int, span.split("-"))
            out.append(dict(record=name, region=region, start=a + offset, end=b + offset, length=b - a,
                            complete=complete))
    return out
