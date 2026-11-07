"""The BUSCO copies' proteins placed on their locus' HMM (since 2026-11-07): hmmalign of each locus' proteins against
BUSCO's profile, read as the HMM match state of each residue (None for a residue in an insert). The match states are
positions shared across species: an intron's place in the gene (introns.py) and the amino acid a codon most likely
codes (codecheck.py) are read through them. The alignments are cached under WORK/align/ by their input's MD5; the
package does not hold them.
"""

import hashlib
import math
import os
import shutil
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


def align_locus(hmm, proteins, hmmalign="hmmalign", cache=None, tmp=None):
    """One locus' proteins ({name: protein}) on its HMM -> {name: [match state or None per residue]}; {} where hmmalign
    could not run. An internal stop (*) is aligned as X; cache: a folder the alignment is kept in, by its input's
    MD5."""
    if not proteins or not hmm or not os.path.isfile(hmm):
        return {}
    text = "".join(f">{n}\n{p.replace('*', 'X')}\n" for n, p in sorted(proteins.items()) if p)
    digest = hashlib.md5((text + util.md5(hmm)).encode()).hexdigest()[:16]
    out = os.path.join(cache, f"{os.path.basename(hmm)[:-4]}.{digest}.pfam") if cache else None
    if not out or not os.path.isfile(out):
        with tempfile.TemporaryDirectory(dir=tmp) as work:
            fasta, got = os.path.join(work, "p.faa"), os.path.join(work, "a.pfam")
            with open(fasta, "w") as fh:
                fh.write(text)
            run = subprocess.run([hmmalign, "--outformat", "Pfam", "-o", got, hmm, fasta], capture_output=True,
                                 text=True, check=False)
            if run.returncode != 0 or not os.path.isfile(got):
                return {}
            if not out:
                return _states(got)
            os.makedirs(cache, exist_ok=True)
            shutil.move(got, util.part_of(out))  # the temporary folder may be on another file system
            os.replace(util.part_of(out), out)
    return _states(out)


def _states(path):
    rows, rf = read_pfam(path)
    return {n: states_of(row, rf) for n, row in rows.items()}


def align(collected, hmm_of, hmmalign="hmmalign", work=None, log=print):
    """Every locus' proteins of the package (collected["proteins"] {(accession, locus): protein}) on their HMMs
    (hmm_of {locus: path}) -> {(accession, locus): [match state or None per residue]}."""
    by_locus = {}
    for (acc, locus), p in collected.get("proteins", {}).items():
        by_locus.setdefault(locus, {})[acc] = p
    cache = os.path.join(work, "align") if work else None
    out, failed = {}, []
    for locus, proteins in sorted(by_locus.items()):
        got = align_locus(hmm_of.get(locus), proteins, hmmalign, cache, work)
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
