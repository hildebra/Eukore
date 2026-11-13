"""The BUSCO copies' proteins placed on their locus' HMM (since 2026-11-07): hmmalign of each locus' proteins against
BUSCO's profile, read as the HMM match state of each residue (None for a residue in an insert). The match states are
positions shared across species: an intron's place in the gene (introns.py) and the amino acid a codon most likely
codes (codecheck.py) are read through them. Each distinct protein's row on the profile is cached under WORK/align/
(since 2026-11-10; until then each locus' alignment whole); the package does not hold them.
"""

import hashlib
import math
import os
import subprocess
import tempfile

from . import util

AMINO = "ACDEFGHIKLMNPQRSTVWY"  # HMMER's order of the match emissions


def states_of(aligned, rf):
    """A sequence's row of a Pfam-format alignment and the alignment's RF line (x: a match column) -> [the match
    state (1-based) of each residue, or None for one in an insert]."""
    out, state = [], 0
    for c, r in zip(aligned, rf):
        match = r not in ".-~"
        state += match
        if c not in "-.~":  # ~: cmalign's missing data (a truncated sequence's ends)
            out.append(state if match else None)
    return out


def read_pfam(path):
    """hmmalign's Pfam-format output -> ({name: aligned row}, RF line)."""
    rows, rf = {}, ""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#=GC RF"):
                rf += line.split()[-1]
            elif line.strip() and not line.startswith(("#", "//")):
                name, row = line.split()
                rows[name] = rows.get(name, "") + row
    return rows, rf


def a2m_of(aligned, rf):
    """A sequence's row of a Pfam-format alignment and its RF line -> the row on the profile alone (A2M): a residue
    on a match state upper case, a match state it skips -, a residue in an insert lower case; nothing of the other
    sequences' inserts."""
    out = []
    for c, r in zip(aligned, rf):
        if r not in ".-~":
            out.append("-" if c in "-.~" else c.upper())
        elif c not in "-.~":
            out.append(c.lower())
    return "".join(out)


def states_of_a2m(a2m):
    """An A2M row (a2m_of) -> [the match state (1-based) of each residue, or None for one in an insert]."""
    out, state = [], 0
    for c in a2m:
        if c == "-":
            state += 1
        elif c.isupper():
            state += 1
            out.append(state)
        else:
            out.append(None)
    return out


def align_locus(hmm, proteins, hmmalign="hmmalign", cache=None, tmp=None):
    """One locus' proteins ({name: protein}) on its HMM -> {name: [match state or None per residue]}; {} where hmmalign
    could not run. An internal stop (*) is aligned as X. cache: a folder each distinct protein's row on the profile
    is kept in (since 2026-11-10: util.cache_read's files <locus>.<the HMM's MD5>, a line a protein, keyed by its MD5;
    until then the locus' whole alignment, by the MD5 of all its proteins: one protein more, a genome added, and every
    protein of the locus was aligned again). hmmalign places each sequence on the profile alone, whatever else it
    aligns (3,796 of 3,796 proteins of 129 loci the same alone as together, HMMER 3.4), so a kept row is the one a run
    with the others gives; and each distinct protein is aligned once."""
    if not proteins or not hmm or not os.path.isfile(hmm):
        return {}
    key_of = {n: hashlib.md5(p.replace("*", "X").encode()).hexdigest() for n, p in proteins.items() if p}
    stem = os.path.join(cache, f"{os.path.basename(hmm)[:-4]}.{util.md5(hmm)[:16]}") if cache else None
    rows = util.cache_read(stem) if stem else {}
    todo = {k: proteins[n].replace("*", "X") for n, k in key_of.items() if k not in rows}
    if todo:
        with tempfile.TemporaryDirectory(dir=tmp) as work:
            fasta, got = os.path.join(work, "p.faa"), os.path.join(work, "a.pfam")
            with open(fasta, "w") as fh:
                fh.writelines(f">{k}\n{p}\n" for k, p in sorted(todo.items()))
            try:
                run = subprocess.run([hmmalign, "--outformat", "Pfam", "-o", got, hmm, fasta], capture_output=True,
                                     text=True, check=False)
            except OSError:  # no such program (since 2026-11-10: raised, the docstring's {} not given)
                return {}
            if run.returncode != 0 or not os.path.isfile(got):
                return {}
            aligned, rf = read_pfam(got)
        new = {k: a2m_of(row, rf) for k, row in aligned.items() if k in todo}
        rows.update(new)
        if stem:
            with util.cache_open(stem) as fh:
                fh.writelines(f"{k}\t{a}\n" for k, a in sorted(new.items()))
    return {n: states_of_a2m(rows[k]) for n, k in key_of.items() if k in rows}


def align(collected, hmm_of, hmmalign="hmmalign", work=None, log=print, threads=1, tmp=None):
    """Every locus' proteins of the package (collected["proteins"] {(accession, locus): protein}) on their HMMs
    (hmm_of {locus: path}), `threads` loci at once (hmmalign is single-threaded: since 2026-11-08; until then one
    locus at a time, ~7 h for all eukaryotes) -> {(accession, locus): [match state or None per residue]}. tmp:
    hmmalign's files (since 2026-11-13; until then in work, the WORK on a shared file system), the cache in work."""
    import concurrent.futures  # noqa: PLC0415
    by_locus = {}
    for (acc, locus), p in collected.get("proteins", {}).items():
        by_locus.setdefault(locus, {})[acc] = p
    cache = os.path.join(work, "align") if work else None
    out, failed = {}, []
    with concurrent.futures.ThreadPoolExecutor(max(1, threads)) as pool:
        jobs = {locus: pool.submit(align_locus, hmm_of.get(locus), proteins, hmmalign, cache, tmp or work)
                for locus, proteins in sorted(by_locus.items())}
        for locus, job in jobs.items():
            got = job.result()
            if not got:
                failed.append(locus)
            out.update({(acc, locus): s for acc, s in got.items()})
    if failed:
        log(f"WARNING: hmmalign placed no proteins of {len(failed)} BUSCOs (e.g. {failed[:3]}): no intron places, "
            "no genetic code check on them")
    return out


def read_emissions(hmm):
    """An HMMER3 profile's match emissions -> [{amino acid: probability} per match state, 1-based: index 0 is None]."""
    out = [None]
    started = False
    with open(hmm) as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "HMM":
                started = True
            elif f[0] == "//":
                break
            elif started and f[0].isdigit() and len(f) >= 21:  # a node's line: its number, then 20 -ln(p)
                out.append({a: math.exp(-float(v)) if v != "*" else 0.0 for a, v in zip(AMINO, f[1:21])})
    return out
