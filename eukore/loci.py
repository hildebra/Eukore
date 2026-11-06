"""The BUSCO genes called in a genome (since 2026-11-05; compleasm's method on BUSCO's eukaryota_odb12, odb.py): Eukore's
own copies, in place of EukDetect1's and EukDetect2's.

1. miniprot aligns the genome's queries (odb.Odb.queries: the nearest species' proteins and the ancestral variants of
   each locus) with the genome's genetic code (NCBI's, `table`: a ciliate's or a diplomonad's TAA and TAG are
   glutamine), keeping alignments of OUTS of a protein's best score (several per protein: a second copy);
2. per locus, the alignments overlapping on a contig's strand are one candidate region, of the alignments with
   MIN_IDENTITY and MIN_COVERAGE of their protein;
3. hmmsearch scores the translations of a region's TOP_K best alignments (by miniprot's score) against every locus'
   HMM: a region is the locus where its best HMM is the locus' own at the locus' BUSCO score cutoff (scores_cutoff) or
   more, and its gene the confirmed alignment of the best HMM score;
4. a confirmed region is complete where its HMM alignment covers COMPLETE of the profile, else a fragment. The locus is
   Single (one complete region), Duplicated (more: the best by score is the copy, the others counted), Fragmented
   (fragments alone: the best), or Missing;
5. the rescue pass (rescue(), since 2026-11-06) of a Missing locus: its HMM searched in the genome's six-frame ORFs
   (esl-translate) places it, and miniprot aligns its queries and the HMM's hits joined (the "self" protein) on that
   place again, with sensitive seeding; the HMMs confirm as in 3. Its copy is `rescued` (a medium confidence: its span
   may be short of the gene's ends, and no own-protein run confirmed it);
6. the gene model (annotated(), since 2026-11-06): in a genome with NCBI gene models (annotation.py), a copy takes the
   annotated transcript sharing the most coding bases with miniprot's model where the locus' HMM confirms its
   translation at ANNOTATION_SHARE of miniprot's score or more (`model` ncbi), else keeps miniprot's.

A copy's coding sequence (`spliced`) leaves out the bases of miniprot's frameshifts (its CIGAR's F and G), which its
GFF keeps inside the CDS; they are `frameshift` features of the copy.

A copy is the genome from the gene's first codon to its last on its own strand, introns and all, uncut and in upper
case (copies.fit, no length cap), with its features on it: exons and introns, the assembly's soft-masked runs of
SOFT_MIN bases or more, and low-complexity runs (DUST's triplet score over DUST_SCORE in DUST_WINDOW bases), which the
user masks as it sees fit.
"""

import collections
import os
import re
import shutil
import subprocess
import tempfile

from . import copies as copies_mod
from . import placement

OUTS = 0.95          # miniprot --outs: a protein's alignments of this share of its best score
OUTN = 10
# compleasm's 0.4 and 0.6 lost up to half the BUSCOs of divergent genomes (2026-11-06, of 129, with the translations
# parsed right: Cyanidioschyzon 33 -> 62 with these gates and the 5 best checked, Bodo saltans 62 -> 77, Galdieria
# 105 -> 112, Plasmodium falciparum 86 -> 88, Dictyostelium 123 -> 126; duplicates 0-3 a genome). No gates at all,
# 30 more queries spread over the eukaryote groups, or all ~100,000 dataset proteins added little or nothing. The
# HMM (the locus' own best, at BUSCO's cutoff) stays the test of a call.
MIN_IDENTITY = 0.2
MIN_COVERAGE = 0.4
TOP_K = 5
COMPLETE = 0.8       # of the HMM's length, by the union of its domains' alignments
SOFT_MIN = 20
DUST_WINDOW = 64
DUST_SCORE = 20.0
WINDOW = 50_000      # a copy's flanks held for a contamination screen (windows.fna.gz)
# The rescue pass (since 2026-11-06). In a genome far from every dataset species, miniprot finds no seed for many genes
# with the distant proteins (Blastocystis hominis without its own: 27 of 98 BUSCOs lost; its genes compact, introns
# ~30 bp), but the HMM finds them in the six-frame ORFs (26 of the 27 on the right gene). Aligned again on that place
# with sensitive seeding, 32 missed went to 9 (all on the right gene, spans shorter: median overlap 0.62).
ORF_MIN = 20                 # esl-translate's shortest ORF (residues)
# An assembly gap's Ns are X in esl-translate's ORFs, not stops: a gap of ~300 kb made six ORFs over hmmsearch's limit
# of HMMSEARCH_MAX residues, which aborts it (Takifugu rubripes, 2026-11-06). The ORFs are split at runs of GAP_X X.
GAP_X = 10
HMMSEARCH_MAX = 100_000
RESCUE_DOM_E = 1e-3          # a domain's independent E-value in the ORFs
RESCUE_JOIN = 20_000         # the hits joined to the best one: on its contig and strand, this far from it
RESCUE_FLANK = 5_000         # bases each side of the hits' place aligned again
RESCUE_OPTIONS = ["-k5", "-M0", "-n1", "-l4"]  # miniprot's seeding: shorter k-mers, all of them, a chain of one
RESCUE_MAX_BASES = 1_000_000_000  # a genome over this is not rescued (its six-frame ORFs' search too long)
# Step 6 (since 2026-11-06): the genome's own NCBI annotation's gene model in place of miniprot's, where the locus' HMM
# scores its translation at this share of miniprot's or more. Against NCBI's curated models miniprot's are exact (every
# intron) for 81% of Drosophila's and Arabidopsis' BUSCOs and 38% of C. elegans'; all their annotated models pass the
# cutoff. Blastocystis hominis' 2011 annotation is weaker: miniprot's model outscores it by 5% or more in 23 of 94.
ANNOTATION_SHARE = 0.95
STATUS = ("Single", "Duplicated", "Fragmented", "Missing")


def frameshifts(cigar):
    """miniprot's CIGAR (cg:Z) -> [(index, n)]: n bases at an index of the coding sequence (the exons joined, on the
    gene's strand) that a frameshift (F, G) leaves out of the reading frame; M and D read 3 bases an amino acid, I none,
    introns (N, U, V) are not in the coding sequence."""
    out, at = [], 0
    for n, op in re.findall(r"(\d+)([MIDFGNUV])", cigar):
        if op in "MD":
            at += 3 * int(n)
        elif op in "FG":
            out.append((at, int(n)))
            at += int(n)
    return out


def in_frame(cds, shifts):
    """A coding sequence without its frameshifts' bases (frameshifts()). miniprot's GFF keeps them in the CDS, so that
    the exons joined read out of frame after the first (until 2026-11-06: 17 of 69 held-out Blastocystis hominis calls
    with an internal stop; left out, 3, and 21 -> 24 of its BUSCOs concordant with OrthoDB's). The exons stay as they
    are: a copy's coding sequence can be a base or two shorter than its exons."""
    if not shifts:
        return cds
    parts, at = [], 0
    for i, n in shifts:
        parts.append(cds[at:i])
        at = i + n
    return "".join(parts) + cds[at:]


def parse_alignments(text):
    """miniprot's GFF with --trans -> [dict(query, locus, contig, start, end, strand, score, identity, qstart, qend,
    exons [(start, end)] in genome order, protein, shifts (frameshifts())] (positions 1-based inclusive). miniprot
    writes an alignment's ##PAF (its CIGAR) and ##STA (its translation) before its mRNA line; until 2026-11-06 the
    translation was taken for the alignment before, so that the HMM scored a neighbour's protein (a copy at one place
    confirmed by another's: false duplicates)."""
    out, cur, pending, cigar = [], None, "", ""
    for line in text.splitlines():
        if line.startswith("##PAF"):
            cigar = next((t[5:] for t in line.split("\t") if t.startswith("cg:Z:")), "")
            continue
        if line.startswith("##STA"):
            pending = line.split("\t", 1)[1].strip() if "\t" in line else line[5:].strip()
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
                       exons=[], protein=pending, shifts=frameshifts(cigar))
            out.append(cur)
            pending, cigar = "", ""
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
        out[locus] = []
        for g in groups:
            ranked = sorted(g["members"], key=lambda a: (a["score"], a["identity"], a["query"]), reverse=True)
            out[locus].append(dict(ranked[0], region=(g["start"], g["end"]), support=len(g["members"]),
                                   top=ranked[:TOP_K]))
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


def genome_bases(path):
    """A genome's bases, estimated from its file's size (placement's ratio for a gzipped one)."""
    return os.path.getsize(path) / (placement.GZIP_BYTES_PER_BASE if path.endswith(".gz") else 1)


def align(target, queries, miniprot, threads, table, work, options=None):
    """miniprot's alignments of queries {name: protein} on a genome or a windows file (in pieces where large) ->
    (alignments, why it failed or None). options: its intron and seeding options (default -I: the max intron 3.6
    times the square root of the target's length)."""
    faa = os.path.join(work, "q.faa")
    with open(faa, "w") as fh:
        fh.writelines(f">{q}\n{p}\n" for q, p in queries.items())
    alignments = []
    for part in placement.chunks(target, work):
        run = subprocess.run([miniprot, "-t", str(threads)] + list(options or ["-I"]) + [
            "-T", str(table or 1), "--gff", "--trans", f"--outs={OUTS}", f"--outn={OUTN}", part, faa],
            capture_output=True, text=True, check=False)
        if run.returncode != 0:
            return [], "miniprot failed: " + (run.stderr.strip().splitlines() or ["?"])[-1][:200]
        alignments += parse_alignments(run.stdout)
    return alignments, None


def confirm(alignments, lengths, odb, hmm_db, hmmsearch_bin, threads, work):
    """The alignments' regions the HMMs confirm (steps 2-3) -> {locus: [each region's confirmed alignment of the best
    HMM score, with hmm_score and hmm_coverage]}; lengths: {query: its protein's length}."""
    found = regions(alignments, lengths)
    candidates = {f"c{i}": (j, a) for i, (j, a) in enumerate(
        ((locus, j), a) for locus, rs in found.items() for j, r in enumerate(rs) for a in r["top"])}
    scores = hmmsearch({n: a["protein"] for n, (_j, a) in candidates.items()}, hmm_db, hmmsearch_bin, threads, work)
    by_region = {}  # (locus, region) -> its confirmed alignment of the best HMM score
    for name, (key, r) in candidates.items():
        got = {odb.locus_of_hmm.get(h, h): v for h, v in scores.get(name, {}).items()}
        if not got:
            continue
        best_hmm = max(got, key=lambda h: got[h][0])
        s, cov = got.get(r["locus"], (0.0, 0.0))
        if best_hmm == r["locus"] and s >= odb.cutoff[r["locus"]] and (
                key not in by_region or (s, r["score"]) > (by_region[key]["hmm_score"], by_region[key]["score"])):
            by_region[key] = dict(r, hmm_score=s, hmm_coverage=cov)
    confirmed = collections.defaultdict(list)
    for (locus, _j), r in sorted(by_region.items()):
        confirmed[locus].append(r)
    return confirmed


def orf_pieces(seq, a, b):
    """An ORF (esl-translate's coords a..b, a > b on the reverse strand) split at its runs of GAP_X or more X (an
    assembly gap) -> [(piece, a, b)], the pieces of ORF_MIN to HMMSEARCH_MAX residues (hmmsearch's limit)."""
    step, out, at = (1 if a <= b else -1), [], 0
    for m in [*re.finditer(f"X{{{GAP_X},}}", seq), None]:
        end = m.start() if m else len(seq)
        if ORF_MIN <= end - at <= HMMSEARCH_MAX:
            out.append((seq[at:end], a + step * 3 * at, a + step * (3 * end - 1)))
        at = m.end() if m else at
    return out


def split_orfs(lines, out):
    """esl-translate's FASTA (lines) written to out as orf_pieces' pieces, each on one line: >name.k source=..
    coords=.. length=.."""
    def flush(header, seq):
        name, _, desc = header[1:].partition(" ")
        kv = dict(x.split("=", 1) for x in desc.split() if "=" in x)
        a, b = (int(v) for v in kv["coords"].split(".."))
        out.writelines(f">{name}.{k} source={kv['source']} coords={x}..{y} length={len(p)}\n{p}\n"
                       for k, (p, x, y) in enumerate(orf_pieces("".join(seq), a, b)))
    header, seq = None, []
    for line in lines:
        if line.startswith(">"):
            if header:
                flush(header, seq)
            header, seq = line.strip(), []
        else:
            seq.append(line.strip())
    if header:
        flush(header, seq)


def orf_hits(target, loci_, odb, table, hmmsearch_bin, threads, work, translate="esl-translate"):
    """The loci's HMMs in a genome's six-frame ORFs (esl-translate, ORF_MIN residues or more, split at assembly gaps:
    split_orfs) -> {locus: [dict(contig, strand, start, end (the ORF's span on the contig), score, hmm_from, hmm_to,
    protein (the domain's envelope on the ORF))]}, the domains of RESCUE_DOM_E or less."""
    hmm, wanted = os.path.join(work, "rescue.hmm"), set(loci_)
    with open(hmm, "w") as out:
        for locus in loci_:
            with open(odb.hmm[locus]) as fh:
                out.write(fh.read())
    hits = collections.defaultdict(list)
    for i, part in enumerate(placement.chunks(target, work)):
        orfs, dom = os.path.join(work, f"orfs{i}.faa"), os.path.join(work, f"orfs{i}.dom")
        err = os.path.join(work, "translate.err")
        with open(orfs, "w") as fh, open(err, "w") as eh:
            run = subprocess.Popen([translate, "-l", str(ORF_MIN), "-c", str(table or 1), part],
                                   stdout=subprocess.PIPE, stderr=eh, text=True)
            split_orfs(run.stdout, fh)
            run.stdout.close()
        if run.wait() != 0:
            with open(err) as eh:
                raise RuntimeError(f"esl-translate failed (exit {run.returncode}): {eh.read().strip()[-200:]}")
        run = subprocess.run([hmmsearch_bin, "--cpu", str(max(1, threads)), "--noali", "-o", os.devnull, "--domE",
                              str(RESCUE_DOM_E), "--domtblout", dom, hmm, orfs], capture_output=True, text=True,
                             check=False)
        if run.returncode != 0:
            raise RuntimeError(f"hmmsearch failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        found = collections.defaultdict(list)  # ORF -> [(locus, score, hmm_from, hmm_to, env_from, env_to)]
        with open(dom) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split()
                locus = odb.locus_of_hmm.get(f[3], f[3])
                if locus in wanted and float(f[12]) <= RESCUE_DOM_E:
                    found[f[0]].append((locus, float(f[13]), int(f[15]), int(f[16]), int(f[19]), int(f[20])))
        with open(orfs) as fh:  # the hit ORFs' places and sequences alone (a large genome's ORFs are many)
            lines = iter(fh)
            for line in lines:
                if not line.startswith(">"):
                    continue
                name, _, desc = line[1:].partition(" ")
                if name not in found:
                    continue
                seq = next(lines).strip()  # split_orfs writes an ORF on one line
                kv = dict(x.split("=", 1) for x in desc.split() if "=" in x)
                a, b = (int(v) for v in kv["coords"].split(".."))  # a > b on the reverse strand
                for locus, score, hf, ht, ef, et in found[name]:
                    hits[locus].append(dict(contig=kv["source"], strand="+" if a <= b else "-", start=min(a, b),
                                            end=max(a, b), score=score, hmm_from=hf, hmm_to=ht,
                                            protein=seq[ef - 1:et]))
        os.remove(orfs)
    return hits


def located(hits, join=RESCUE_JOIN):
    """Each locus' place from its ORF hits (orf_hits) -> {locus: (contig, start, end, self protein)}: its best hit and
    the hits on that contig and strand within `join` bases of it, those of the best scores that do not overlap on the
    HMM; the self protein their envelopes joined in the HMM's order (nearly the genome's own protein)."""
    out = {}
    for locus, hs in hits.items():
        best = max(hs, key=lambda h: h["score"])
        near = [h for h in hs if (h["contig"], h["strand"]) == (best["contig"], best["strand"])
                and h["start"] <= best["end"] + join and h["end"] >= best["start"] - join]
        chosen = []
        for h in sorted(near, key=lambda h: -h["score"]):
            if all(h["hmm_to"] < c["hmm_from"] or h["hmm_from"] > c["hmm_to"] for c in chosen):
                chosen.append(h)
        chosen.sort(key=lambda h: h["hmm_from"])
        out[locus] = (best["contig"], min(h["start"] for h in chosen), max(h["end"] for h in chosen),
                      "".join(h["protein"] for h in chosen))
    return out


def rescue(genome_path, missing, queries, odb, miniprot, hmmsearch_bin, threads, table, hmm_db, work,
           translate="esl-translate"):
    """The rescue pass (step 5) of the loci the first pass left Missing -> {locus: [confirmed alignments, rescued]}:
    the locus' HMM in the genome's six-frame ORFs places it (located); on that place and RESCUE_FLANK bases each side,
    miniprot with RESCUE_OPTIONS aligns the locus' queries and its self protein; the HMMs confirm as in the first pass.
    Introns as long as the first pass allows (-G: -I's 3.6 times the square root of the genome's length)."""
    places = located(orf_hits(genome_path, missing, odb, table, hmmsearch_bin, threads, work, translate))
    if not places:
        return {}
    by_contig = collections.defaultdict(list)
    for locus, (contig, start, end, _p) in places.items():
        by_contig[contig].append((locus, start, end))
    contigs = placement.read_contigs(genome_path, set(by_contig))
    target = os.path.join(work, "rescue.fna")
    starts = {}  # window -> (locus, contig, its first base on the contig)
    with open(target, "w") as fh:
        for contig, wanted in sorted(by_contig.items()):
            seq = contigs[contig]
            for locus, start, end in wanted:
                lo, hi = max(1, start - RESCUE_FLANK), min(len(seq), end + RESCUE_FLANK)
                name = f"w{len(starts)}"
                starts[name] = (locus, contig, lo)
                fh.write(f">{name}\n{seq[lo - 1:hi]}\n")
    q = {n: p for n, p in queries.items() if n.split("~", 1)[0] in places}
    q.update({f"{locus}~{locus}_self": p[3] for locus, p in places.items() if p[3]})
    intron = int(3.6 * genome_bases(genome_path) ** 0.5)
    alignments, failed = align(target, q, miniprot, threads, table, work, ["-G", str(intron)] + RESCUE_OPTIONS)
    if failed:
        raise RuntimeError(f"rescue: {failed}")
    kept = []
    for a in alignments:  # back on the genome; another locus' alignment in a window is not this rescue's
        locus, contig, lo = starts[a["contig"]]
        if a["locus"] != locus:
            continue
        shift = lo - 1
        kept.append(dict(a, contig=contig, start=a["start"] + shift, end=a["end"] + shift,
                         exons=[(x + shift, y + shift) for x, y in a["exons"]]))
    confirmed = confirm(kept, {n: len(p) for n, p in q.items()}, odb, hmm_db, hmmsearch_bin, threads, work)
    return {locus: [dict(r, rescued=True) for r in rs] for locus, rs in confirmed.items()
            if locus in places and locus in missing}


def frameshift_features(exons, shifts):
    """A copy's frameshifts (shifts: frameshifts() on its coding sequence) as features: [("frameshift", start, end)]
    on the copy (0-based, end excluded), exons: the copy's [(start, end)] in its order."""
    out = []
    for i, n in shifts or ():
        at = 0
        for s, e in sorted(exons):
            for k in range(n):  # a frameshift's bases, each where the exons put it (rarely across two exons)
                if at <= i + k < at + e - s:
                    p = s + i + k - at
                    if out and out[-1][2] == p:
                        out[-1] = ("frameshift", out[-1][1], p + 1)
                    else:
                        out.append(("frameshift", p, p + 1))
            at += e - s
    return out


def annotated(chosen, contigs, annotation, odb, hmmsearch_bin, threads, work, hmm_db, table):
    """Step 6: the genome's own annotated gene model (annotation: (NCBI GFF, assembly report or None)) in place of a
    call's, where the GFF holds a transcript on the call's contig and strand sharing coding bases with it (the most:
    annotation_mod.best_transcript) whose translation the locus' own HMM confirms (its best HMM, at the cutoff) at
    ANNOTATION_SHARE of the call's HMM score or more. chosen: {locus: alignment} -> the same, an annotated model with
    model="ncbi", the transcript, and miniprot's model kept as miniprot=dict(start, end, exons) (annotation.py's
    benchmark compares that). The annotation's exons are its CDS: its 5' piece's phase left out (a partial gene)."""
    from . import annotation as annotation_mod  # noqa: PLC0415
    from .benchmarks import translate  # noqa: PLC0415 (benchmarks imports this module's users)
    gff, report = annotation
    spans = collections.defaultdict(list)
    for a in chosen.values():
        spans[a["contig"]].append((a["start"], a["end"]))
    transcripts = annotation_mod.read_cds(gff, spans, annotation_mod.aliases(report))
    models, proteins = {}, {}
    for locus, a in chosen.items():
        name = annotation_mod.best_transcript(a, transcripts)
        if not name:
            continue
        t = transcripts[name]
        exons = [list(x) for x in t["cds"]]
        five = exons[0] if t["strand"] == "+" else exons[-1]
        phase = t.get("phase", {}).get(tuple(five), 0)
        if phase:
            five[0 if t["strand"] == "+" else 1] += phase if t["strand"] == "+" else -phase
        exons = [tuple(x) for x in exons if x[0] <= x[1]]
        seq = contigs[a["contig"]]
        cds = "".join(seq[x - 1:y] for x, y in exons).upper()
        cds = cds.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1] if t["strand"] == "-" else cds
        models[locus] = dict(exons=exons, start=exons[0][0], end=exons[-1][1], transcript=name)
        proteins[locus] = translate(cds, table).rstrip("*")
    scores = hmmsearch(proteins, hmm_db, hmmsearch_bin, threads, work)
    out = dict(chosen)
    for locus, m in models.items():
        got = {odb.locus_of_hmm.get(h, h): v for h, v in scores.get(locus, {}).items()}
        if not got or max(got, key=lambda h: got[h][0]) != locus:
            continue
        s, cov = got[locus]
        a = chosen[locus]
        if s >= odb.cutoff[locus] and s >= ANNOTATION_SHARE * a["hmm_score"]:
            out[locus] = dict(a, exons=m["exons"], start=m["start"], end=m["end"], shifts=[], hmm_score=s,
                              hmm_coverage=cov, protein=proteins[locus], model="ncbi", transcript=m["transcript"],
                              miniprot=dict(start=a["start"], end=a["end"], exons=a["exons"]))
    return out


def call(genome_path, queries, odb, miniprot="miniprot", hmmsearch_bin="hmmsearch", threads=1, tmp=None,
         hmm_db=None, table=1, rescue_missing=True, translate="esl-translate", annotation=None):
    """The BUSCO genes of a genome (table: its NCBI genetic code) -> dict(loci {locus: dict(status, copy, spliced,
    features [(kind, start, end)],
    contig, start, end, strand, identity, coverage, score, hmm_score, hmm_coverage, cutoff, copies, query, protein,
    rescued, model (miniprot or ncbi), miniprot (its model: start, end, exons))}, failed (why the genome could not be
    called, or None), windows {name: sequence} of the copies' flanks).
    rescue_missing: the rescue pass (step 5) of the loci left Missing, in a genome of RESCUE_MAX_BASES or less;
    annotation: (the genome's NCBI GFF, its assembly report or None) for step 6 (annotated()), or None."""
    lengths = {q: len(p) for q, p in queries.items()}
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        hmm_db = hmm_db or odb.hmm_db(os.path.join(work, "loci.hmm"))
        alignments, failed = align(genome_path, queries, miniprot, threads, table, work)
        if failed:
            return dict(loci={}, failed=failed, windows={})
        confirmed = confirm(alignments, lengths, odb, hmm_db, hmmsearch_bin, threads, work)
        missing = [locus for locus in odb.loci if not confirmed.get(locus)]
        if rescue_missing and missing and genome_bases(genome_path) <= RESCUE_MAX_BASES and shutil.which(translate):
            confirmed.update(rescue(genome_path, missing, queries, odb, miniprot, hmmsearch_bin, threads, table,
                                    hmm_db, work, translate))
        contigs = placement.read_contigs(genome_path, {r["contig"] for rs in confirmed.values() for r in rs})
        chosen, statuses = {}, {}
        for locus in odb.loci:
            rs = confirmed.get(locus, [])
            complete = [r for r in rs if r["hmm_coverage"] >= COMPLETE]
            if rs:
                chosen[locus] = max(complete or rs, key=lambda r: (r["hmm_score"], r["score"]))
                statuses[locus] = ("Single" if len(complete) == 1 else "Duplicated" if complete else "Fragmented",
                                   complete)
        if annotation and annotation[0] and chosen:
            chosen = annotated(chosen, contigs, annotation, odb, hmmsearch_bin, threads, work, hmm_db, table)
    loci, windows = {}, {}
    for locus in odb.loci:
        if locus not in chosen:
            loci[locus] = dict(status="Missing")
            continue
        best = chosen[locus]
        status, complete = statuses[locus]
        seq = contigs[best["contig"]]
        fitted = copies_mod.fit(seq, best["start"], best["end"], best["strand"], best["exons"], 10 ** 12,
                                mask_min=10 ** 12)
        raw = seq[best["start"] - 1:best["end"]]
        feats = [(k, s, e) for k, s, e, _r, _g1, _g2 in fitted["features"] if k in ("exon", "intron")]
        feats += frameshift_features([(s, e) for k, s, e in feats if k == "exon"], best.get("shifts"))
        feats += [("soft_masked", a, b) for a, b in on_copy(soft_runs(raw), best["start"], best["end"], best["strand"])]
        feats += [("low_complexity", a, b) for a, b in dust(fitted["copy"])]
        loci[locus] = dict(status=status, copy=fitted["copy"], spliced=in_frame(fitted["spliced"], best.get("shifts")),
                           features=sorted(
            feats, key=lambda f: (f[1], f[0])), contig=best["contig"], contig_length=len(seq), start=best["start"],
            end=best["end"], strand=best["strand"], identity=best["identity"], coverage=best["coverage"],
            score=best["score"], hmm_score=best["hmm_score"], hmm_coverage=best["hmm_coverage"],
            cutoff=odb.cutoff[locus], copies=len(complete), others=[(r["contig"], r["start"], r["end"], r["hmm_score"])
            for r in complete if r is not best], query=best["query"].split("~", 1)[1],
            protein=best["protein"], exons=best["exons"], rescued=bool(best.get("rescued")),
            model=best.get("model", "miniprot"), transcript=best.get("transcript", "-"),
            miniprot=best.get("miniprot") or dict(start=best["start"], end=best["end"], exons=best["exons"]))
        w1, w2 = max(1, best["start"] - WINDOW), min(len(seq), best["end"] + WINDOW)
        windows[f"{locus}|{best['contig']}|{w1}-{w2}"] = seq[w1 - 1:w2].upper()
    return dict(loci=loci, failed=None, windows=windows)
