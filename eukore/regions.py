"""The 18S variable regions of each rRNA record (since 2026-11-07; rrna/regions.tsv): V1-V5 and V7-V9 as Hadziavdic et
al. 2014 place them on Saccharomyces cerevisiae's 18S (data/rrna_regions.tsv), carried to every record through
Rfam's eukaryote SSU model (RF01960, SSU_rRNA_eukarya): cmalign aligns the yeast 18S once, which gives each region's
first and last model positions, then every record's 18S; a region is the record's bases from the first at its first
position or after to the last at its last position or before (inserts between included). It is complete where the
record holds bases within SLACK positions of both its ends. The alignments are cached by sequence under WORK/cmalign/.
"""

import gzip
import hashlib
import os
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


def cm_path(src_dir):
    return os.path.join(src_dir, "rfam", f"{SSU_CM[0]}.cm")


def fetch(src_dir, download, log=print):
    """Rfam's SSU model into src_dir/rfam (kept when there)."""
    path = cm_path(src_dir)
    if not os.path.isfile(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        log(f"Rfam: fetching {SSU_CM[0]}")
        download(SSU_CM[1], path)
    return path


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


def cmalign(cm, seqs, cmalign_bin="cmalign", threads=1, tmp=None, log=None):
    """seqs {name: sequence} aligned to the CM -> {name: [the model position of each base, or None in an insert]}. A
    sequence over MAX_SSU is left out; a batch cmalign fails on (one sequence over its matrix limit) is aligned one
    sequence at a time, those that fail left out (said)."""
    out, failed = {}, []
    names = sorted(n for n in seqs if len(seqs[n]) <= MAX_SSU)
    for i in range(0, len(names), BATCH):
        got, why = _cmalign(cm, seqs, names[i:i + BATCH], cmalign_bin, threads, tmp)
        if got is None:
            for n in names[i:i + BATCH]:
                one, why = _cmalign(cm, seqs, [n], cmalign_bin, threads, tmp)
                if one is None:
                    failed.append((n, why))
                else:
                    out.update(one)
        else:
            out.update(got)
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
    cache_path = os.path.join(work, "cmalign", f"ssu.{util.md5(cm)[:10]}.tsv.gz") if work else None
    cache = {}
    if cache_path and os.path.isfile(cache_path):
        with gzip.open(cache_path, "rt") as fh:
            for line in fh:
                key, text = line.rstrip("\n").split("\t")
                cache[key] = text
    ref = cmalign(cm, {"reference": reference()}, cmalign_bin, threads, work)["reference"]
    cols = columns(ref, regions)
    ssu = {}
    for name, seq, row in records:
        got = ssu_of(seq, row)
        if got:
            ssu[name] = (got[0].upper(), got[1], hashlib.md5(got[0].upper().encode()).hexdigest())
    todo = {key: s for s, _o, key in ssu.values() if key not in cache}
    if todo:
        log(f"18S variable regions: {len(todo)} distinct 18S to align (cmalign)")
        aligned = cmalign(cm, todo, cmalign_bin, threads, work, log)
        for key in todo:
            found = locate(aligned.get(key, []), cols)
            cache[key] = ",".join(f"{n}:{a}-{b}:{'t' if c else 'f'}" for n, (a, b, c) in found.items()) or "-"
        if cache_path:
            os.makedirs(os.path.dirname(cache_path), exist_ok=True)
            with gzip.open(util.part_of(cache_path), "wt") as fh:
                fh.writelines(f"{k}\t{v}\n" for k, v in sorted(cache.items()))
            os.replace(util.part_of(cache_path), cache_path)
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
