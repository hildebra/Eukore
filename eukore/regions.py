"""The 18S variable regions of each rRNA record (since 2026-11-07; rrna/regions.tsv): V1-V5 and V7-V9 as Hadziavdic et
al. 2014 place them on Saccharomyces cerevisiae's 18S (data/rrna_regions.tsv), carried to every record through
Rfam's eukaryote SSU model (RF01960, SSU_rRNA_eukarya): cmalign aligns the yeast 18S once, which gives each region's
first and last model positions, then every record's 18S; a region is the record's bases from the first at its first
position or after to the last at its last position or before (inserts between included). It is complete where the
record holds bases within SLACK positions of both its ends. The alignments are cached by sequence under WORK/cmalign/.
"""

import collections
import glob
import hashlib
import os
import socket
import subprocess
import tempfile

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


def cm_path(src_dir):
    return os.path.join(src_dir, "rfam", f"{SSU_CM[0]}.cm")


def models_path(src_dir):
    """The SSU models' file (SSU_MODELS, one Infernal file; ssu_models)."""
    return os.path.join(src_dir, "rfam", SSU_MODELS_FILE)


def fetch(src_dir, download, log=print):
    """Rfam's SSU model into src_dir/rfam (kept when there and an Infernal model: a page saved in its place, e.g. a
    proxy's, is fetched again; since 2026-11-08 it was kept, and every package stage then stopped), and (since
    2026-11-09) SSU_MODELS' joined into models_path's file."""
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
    if not is_cm(joined):
        with open(util.part_of(joined), "w", newline="\n") as out:
            for acc, _kind in SSU_MODELS:
                with open(os.path.join(folder, f"{acc}.cm")) as fh:
                    out.write(fh.read())
        os.replace(util.part_of(joined), joined)
    return path


def ssu_models(seqs, cm, cmsearch="cmsearch", threads=1, tmp=None):
    """Each sequence's SSU model of the best score ({name: its 18S} -> {name: eukarya, bacteria, archaea, or - where
    none scores it}; since 2026-11-09): cmsearch of SSU_MODELS (cm: models_path's file) by their profile HMMs alone
    (--hmmonly: an HMM search's speed), each distinct sequence once; a prokaryote model's call needs its score over
    the eukaryote models' best by SSU_MARGIN bits."""
    kind_of = dict(SSU_MODELS)
    distinct = {}
    for name, seq in seqs.items():
        if seq:
            distinct.setdefault(seq.upper(), []).append(name)
    if not distinct:
        return {}
    best = collections.defaultdict(dict)  # sequence index -> {kind: best score}
    order = list(distinct)
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        fasta, table = os.path.join(work, "s.fna"), os.path.join(work, "hits.tbl")
        with open(fasta, "w") as fh:
            fh.writelines(f">s{i}\n{seq}\n" for i, seq in enumerate(order))
        run = subprocess.run([cmsearch, "--hmmonly", "--noali", "-Z", "1", "--cpu", str(max(1, threads)), "--tblout",
                              table, "-o", os.devnull, cm, fasta], capture_output=True, text=True, check=False)
        if run.returncode != 0:
            raise RuntimeError(f"cmsearch failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        with open(table) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                kind, score = kind_of.get(f[3]), float(f[14])
                if kind and score > best[int(f[0][1:])].get(kind, float("-inf")):
                    best[int(f[0][1:])][kind] = score
    out = {}
    for i, seq in enumerate(order):
        scores = best.get(i)
        if not scores:
            kind = "-"
        else:
            kind = max(scores, key=scores.get)
            if kind != "eukarya" and scores[kind] < scores.get("eukarya", float("-inf")) + SSU_MARGIN:
                kind = "eukarya"
        for name in distinct[seq]:
            out[name] = kind
    return out


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


def _cmalign(cm, seqs, names, cmalign_bin, threads, tmp):
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        fasta, got = os.path.join(work, "s.fna"), os.path.join(work, "a.pfam")
        with open(fasta, "w") as fh:
            fh.writelines(f">{n}\n{seqs[n]}\n" for n in names)
        run = subprocess.run([cmalign_bin, "--cpu", str(max(1, threads)), "--outformat", "Pfam", "--noprob",
                              "--mxsize", str(MXSIZE), "-o", got, cm, fasta], capture_output=True, text=True,
                             check=False)
        if run.returncode != 0:
            return None, run.stderr.strip()[-300:]
        rows, rf = align_mod.read_pfam(got)
        return {n: align_mod.states_of(row, rf) for n, row in rows.items()}, None


def cmalign(cm, seqs, cmalign_bin="cmalign", threads=1, tmp=None, log=None, each=None):
    """seqs {name: sequence} aligned to the CM -> {name: [the model position of each base, or None in an insert]}. A
    sequence over MAX_SSU is left out; a batch cmalign fails on (one sequence over its matrix limit) is aligned one
    sequence at a time, those that fail left out (said). each: a function given each batch's alignments as they come
    (run's cache)."""
    out, failed = {}, []
    names = sorted(n for n in seqs if len(seqs[n]) <= MAX_SSU)
    for i in range(0, len(names), BATCH):
        got, why = _cmalign(cm, seqs, names[i:i + BATCH], cmalign_bin, threads, tmp)
        if got is None:
            got = {}
            for n in names[i:i + BATCH]:
                one, why = _cmalign(cm, seqs, [n], cmalign_bin, threads, tmp)
                if one is None:
                    failed.append((n, why))
                else:
                    got.update(one)
        out.update(got)
        if each:
            each(got)
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


def run(records, cm, cmalign_bin="cmalign", threads=1, work=None, log=print):
    """Every record's 18S variable regions (records [(name, sequence, row)]) -> rows (REGION_COLUMNS). Each distinct 18S
    is aligned once (cached under WORK/cmalign/)."""
    regions = read_regions()
    ref = cmalign(cm, {"reference": reference()}, cmalign_bin, threads, work).get("reference")
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
    cache = {}
    for path in sorted(glob.glob(glob.escape(stem) + "*.tsv")) if stem else ():
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if line.endswith("\n") and len(f) == 2:  # a line a killed run cut: aligned again
                    cache[f[0]] = f[1]
    cache_path = f"{stem}.{socket.gethostname()}-{os.getpid()}.tsv" if stem else None
    ssu = {}
    for name, seq, row in records:
        got = ssu_of(seq, row)
        if got:
            ssu[name] = (got[0].upper(), got[1], hashlib.md5(got[0].upper().encode()).hexdigest())
    todo = {key: s for s, _o, key in ssu.values() if key not in cache}
    if todo:
        log(f"18S variable regions: {len(todo)} distinct 18S to align (cmalign)")
        out = None
        if cache_path:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            cut = os.path.isfile(cache_path) and os.path.getsize(cache_path) and not util.ends_in_newline(cache_path)
            out = open(cache_path, "a", newline="\n")
            if cut:
                out.write("\n")

        def keep(aligned):
            got = {key: ",".join(f"{n}:{a}-{b}:{'t' if c else 'f'}" for n, (a, b, c) in
                                 locate(states, cols).items()) or "-" for key, states in aligned.items()}
            cache.update(got)
            if out:
                out.writelines(f"{k}\t{v}\n" for k, v in sorted(got.items()))
                out.flush()
        try:
            keep({key: [] for key, s in todo.items() if len(s) > MAX_SSU})  # never aligned: none, for good
            cmalign(cm, todo, cmalign_bin, threads, work, log, each=keep)
        finally:
            if out:
                out.close()
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
