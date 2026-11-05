"""The BUSCO genes called in a genome (since 2026-11-05; compleasm's method on BUSCO's eukaryota_odb12, odb.py): EUKDB's
own copies, in place of EukDetect1's and EukDetect2's.

1. miniprot aligns the genome's queries (odb.Odb.queries: the nearest species' proteins and the ancestral variants of
   each locus), keeping alignments of OUTS of a protein's best score (several per protein: a second copy);
2. per locus, the alignments overlapping on a contig's strand are one candidate region; a region's best alignment (by
   miniprot's score) is its gene, with MIN_IDENTITY and MIN_COVERAGE of its protein (compleasm's 0.4 and 0.6);
3. hmmsearch scores each region's translated gene against every locus' HMM: a region is the locus where its best HMM is
   the locus' own at the locus' BUSCO score cutoff (scores_cutoff) or more;
4. a confirmed region is complete where its HMM alignment covers COMPLETE of the profile, else a fragment. The locus is
   Single (one complete region), Duplicated (more: the best by score is the copy, the others counted), Fragmented
   (fragments alone: the best), or Missing.

A copy is the genome from the gene's first codon to its last on its own strand, introns and all, uncut and in upper
case (copies.fit, no length cap), with its features on it: exons and introns, the assembly's soft-masked runs of
SOFT_MIN bases or more, and low-complexity runs (DUST's triplet score over DUST_SCORE in DUST_WINDOW bases), which the
user masks as it sees fit.
"""

import collections
import os
import subprocess
import tempfile

from . import copies as copies_mod
from . import placement

OUTS = 0.95          # miniprot --outs: a protein's alignments of this share of its best score
OUTN = 10
MIN_IDENTITY = 0.4   # compleasm's defaults
MIN_COVERAGE = 0.6
COMPLETE = 0.8       # of the HMM's length, by the union of its domains' alignments
SOFT_MIN = 20
DUST_WINDOW = 64
DUST_SCORE = 20.0
WINDOW = 50_000      # a copy's flanks held for a contamination screen (windows.fna.gz)
STATUS = ("Single", "Duplicated", "Fragmented", "Missing")


def parse_alignments(text):
    """miniprot's GFF with --trans -> [dict(query, locus, contig, start, end, strand, score, identity, qstart, qend,
    exons [(start, end)] in genome order, protein)] (positions 1-based inclusive)."""
    out, cur = [], None
    for line in text.splitlines():
        if line.startswith("##STA"):
            if cur is not None:
                cur["protein"] = line.split("\t", 1)[1].strip() if "\t" in line else line[5:].strip()
            continue
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        attrs = dict(a.split("=", 1) for a in f[8].split(";") if "=" in a)
        if f[2] == "mRNA":
            target = attrs.get("Target", "").split()
            cur = dict(query=target[0], locus=target[0].split("~", 1)[0], contig=f[0], start=int(f[3]),
                       end=int(f[4]), strand=f[6], score=float(f[5]) if f[5] != "." else 0.0,
                       identity=float(attrs.get("Identity", 0)), qstart=int(target[1]), qend=int(target[2]),
                       exons=[], protein="")
            out.append(cur)
        elif f[2] == "CDS" and cur is not None:
            cur["exons"].append((int(f[3]), int(f[4])))
    for a in out:
        a["exons"].sort()
    return out


def regions(alignments, lengths):
    """Per locus, the alignments' regions: overlapping alignments on a contig's strand are one -> {locus: [best
    alignment of each region, with its region's span]}; alignments under MIN_IDENTITY or MIN_COVERAGE of their protein
    (lengths: {query: protein length}) left out."""
    by_locus = collections.defaultdict(list)
    for a in alignments:
        cov = (a["qend"] - a["qstart"] + 1) / max(lengths.get(a["query"], 0), 1)
        if a["identity"] >= MIN_IDENTITY and cov >= MIN_COVERAGE:
            by_locus[a["locus"]].append(dict(a, coverage=cov))
    out = {}
    for locus, hits in by_locus.items():
        groups = []
        for a in sorted(hits, key=lambda a: (a["contig"], a["strand"], a["start"])):
            g = groups[-1] if groups else None
            if g and g["contig"] == a["contig"] and g["strand"] == a["strand"] and a["start"] <= g["end"]:
                g["end"] = max(g["end"], a["end"])
                g["members"].append(a)
            else:
                groups.append(dict(contig=a["contig"], strand=a["strand"], start=a["start"], end=a["end"],
                                   members=[a]))
        out[locus] = [dict(max(g["members"], key=lambda a: (a["score"], a["identity"], a["query"])),
                           region=(g["start"], g["end"]), support=len(g["members"])) for g in groups]
    return out


def hmmsearch(proteins, hmm_db, binary="hmmsearch", threads=1, tmp=None):
    """Each protein's scores against the HMMs -> {protein: {hmm: (bit score, HMM coverage)}}; HMM coverage: the union
    of its domains' hmm_from..hmm_to over the HMM's length."""
    if not proteins:
        return {}
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        faa, dom = os.path.join(work, "c.faa"), os.path.join(work, "dom.tsv")
        with open(faa, "w") as fh:
            fh.writelines(f">{n}\n{s.replace('*', 'X')}\n" for n, s in proteins.items() if s)
        run = subprocess.run([binary, "--cpu", str(max(1, threads)), "--noali", "-o", os.devnull, "--domtblout", dom,
                              hmm_db, faa], capture_output=True, text=True, check=False)
        if run.returncode != 0:
            raise RuntimeError(f"hmmsearch failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        score, spans, qlen = {}, collections.defaultdict(list), {}
        with open(dom) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                target, query, length, full = f[0], f[3], int(f[5]), float(f[7])
                score[(target, query)] = full
                qlen[query] = length
                spans[(target, query)].append((int(f[15]), int(f[16])))
    out = collections.defaultdict(dict)
    for key, s in score.items():
        covered = set()
        for a, b in spans[key]:
            covered.update(range(a, b + 1))
        out[key[0]][key[1]] = (s, len(covered) / max(qlen[key[1]], 1))
    return out


def dust(seq, window=DUST_WINDOW, threshold=DUST_SCORE):
    """Low-complexity runs of a sequence (DUST: in a window, the sum over its triplets of c(c-1)/2 over the window's
    triplets less one) -> [(start, end)] 0-based, end excluded, merged."""
    seq = seq.upper()
    n = len(seq)
    if n < window:
        return []
    marked = []
    step = window // 2
    for s in range(0, n - window + 1, step):
        counts = collections.Counter(seq[i:i + 3] for i in range(s, s + window - 2))
        score = sum(c * (c - 1) / 2 for c in counts.values()) / (window - 3)
        if score > threshold:
            marked.append((s, s + window))
    if n % step and n > window:  # the tail
        s = n - window
        counts = collections.Counter(seq[i:i + 3] for i in range(s, n - 2))
        if sum(c * (c - 1) / 2 for c in counts.values()) / (window - 3) > threshold:
            marked.append((s, n))
    merged = []
    for a, b in sorted(marked):
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def soft_runs(seq, min_len=SOFT_MIN):
    """Runs of lower-case bases (an assembly's soft mask) of min_len or more -> [(start, end)] 0-based."""
    out, start = [], None
    for i, c in enumerate(seq + "."):
        if c.islower():
            start = i if start is None else start
        else:
            if start is not None and i - start >= min_len:
                out.append((start, i))
            start = None
    return out


def on_copy(intervals, start, end, strand):
    """Genome intervals of a copy's span (0-based on the span, its genome order) on the copy's own strand."""
    length = end - start + 1
    return sorted((a, b) if strand == "+" else (length - b, length - a) for a, b in intervals)


def call(genome_path, queries, odb, miniprot="miniprot", hmmsearch_bin="hmmsearch", threads=1, tmp=None,
         hmm_db=None):
    """The BUSCO genes of a genome -> dict(loci {locus: dict(status, copy, spliced, features [(kind, start, end)],
    contig, start, end, strand, identity, coverage, score, hmm_score, hmm_coverage, cutoff, copies, query, protein)},
    failed (why the genome could not be called, or None), windows {name: sequence} of the copies' flanks)."""
    lengths = {q: len(p) for q, p in queries.items()}
    alignments = []
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        faa = os.path.join(work, "q.faa")
        with open(faa, "w") as fh:
            fh.writelines(f">{q}\n{p}\n" for q, p in queries.items())
        for part in placement.chunks(genome_path, work):
            run = subprocess.run([miniprot, "-t", str(threads), "-I", "--gff", "--trans", f"--outs={OUTS}",
                                  f"--outn={OUTN}", part, faa], capture_output=True, text=True, check=False)
            if run.returncode != 0:
                return dict(loci={}, failed="miniprot failed: " + (run.stderr.strip().splitlines() or ["?"])[-1][:200],
                            windows={})
            alignments += parse_alignments(run.stdout)
        found = regions(alignments, lengths)
        candidates = {f"c{i}": r for i, r in enumerate(r for rs in found.values() for r in rs)}
        scores = hmmsearch({n: r["protein"] for n, r in candidates.items()}, hmm_db or odb.hmm_db(
            os.path.join(work, "loci.hmm")), hmmsearch_bin, threads, work)
    confirmed = collections.defaultdict(list)
    for name, r in candidates.items():
        got = {odb.locus_of_hmm.get(h, h): v for h, v in scores.get(name, {}).items()}
        if not got:
            continue
        best_hmm = max(got, key=lambda h: got[h][0])
        s, cov = got.get(r["locus"], (0.0, 0.0))
        if best_hmm == r["locus"] and s >= odb.cutoff[r["locus"]]:
            confirmed[r["locus"]].append(dict(r, hmm_score=s, hmm_coverage=cov))
    contigs = placement.read_contigs(genome_path, {r["contig"] for rs in confirmed.values() for r in rs})
    loci, windows = {}, {}
    for locus in odb.loci:
        rs = confirmed.get(locus, [])
        complete = [r for r in rs if r["hmm_coverage"] >= COMPLETE]
        if not rs:
            loci[locus] = dict(status="Missing")
            continue
        pool = complete or rs
        best = max(pool, key=lambda r: (r["hmm_score"], r["score"]))
        status = "Single" if len(complete) == 1 else "Duplicated" if complete else "Fragmented"
        seq = contigs[best["contig"]]
        fitted = copies_mod.fit(seq, best["start"], best["end"], best["strand"], best["exons"], 10 ** 12,
                                mask_min=10 ** 12)
        raw = seq[best["start"] - 1:best["end"]]
        feats = [(k, s, e) for k, s, e, _r, _g1, _g2 in fitted["features"] if k in ("exon", "intron")]
        feats += [("soft_masked", a, b) for a, b in on_copy(soft_runs(raw), best["start"], best["end"], best["strand"])]
        feats += [("low_complexity", a, b) for a, b in dust(fitted["copy"])]
        loci[locus] = dict(status=status, copy=fitted["copy"], spliced=fitted["spliced"], features=sorted(
            feats, key=lambda f: (f[1], f[0])), contig=best["contig"], contig_length=len(seq), start=best["start"],
            end=best["end"], strand=best["strand"], identity=best["identity"], coverage=best["coverage"],
            score=best["score"], hmm_score=best["hmm_score"], hmm_coverage=best["hmm_coverage"],
            cutoff=odb.cutoff[locus], copies=len(complete), query=best["query"].split("~", 1)[1],
            protein=best["protein"], exons=best["exons"])
        w1, w2 = max(1, best["start"] - WINDOW), min(len(seq), best["end"] + WINDOW)
        windows[f"{locus}|{best['contig']}|{w1}-{w2}"] = seq[w1 - 1:w2].upper()
    return dict(loci=loci, failed=None, windows=windows)
