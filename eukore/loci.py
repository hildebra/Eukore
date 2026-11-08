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
6. the gene model (annotated(), since 2026-11-06): in a genome with its own gene models (annotation.py: NCBI's, else
   Ensembl's, else TOGA2's), a copy takes the annotated transcript sharing the most coding bases with miniprot's model
   where the locus' HMM confirms its translation at ANNOTATION_SHARE of miniprot's score or more (`model` the source:
   ncbi, ensembl, toga2), else keeps miniprot's.

A copy's coding sequence (`spliced`) leaves out the bases of miniprot's frameshifts (its CIGAR's F and G), which its
GFF keeps inside the CDS; they are `frameshift` features of the copy.

A copy is the genome from the gene's first codon to its last on its own strand, introns and all, uncut and in upper
case (copies.whole, no length cap), with its features on it: exons and introns, the assembly's soft-masked runs of
SOFT_MIN bases or more, and low-complexity runs (DUST's triplet score over DUST_SCORE in DUST_WINDOW bases), which the
user masks as it sees fit.
"""

import collections
import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
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
# Step 6 (since 2026-11-06): the genome's own annotation's gene model in place of miniprot's, where the locus' HMM
# scores its translation at this share of miniprot's or more. Against NCBI's curated models miniprot's are exact (every
# intron) for 81% of Drosophila's and Arabidopsis' BUSCOs and 38% of C. elegans'; all their annotated models pass the
# cutoff. Blastocystis hominis' 2011 annotation is weaker: miniprot's model outscores it by 5% or more in 23 of 94.
ANNOTATION_SHARE = 0.95
STATUS = ("Single", "Duplicated", "Fragmented", "Missing")


def frameshifts(cigar):
    """miniprot's CIGAR (cg:Z) -> [(index, n)]: n bases at an index of the coding sequence (the exons joined, on the
    gene's strand) that a frameshift (F, G) leaves out of the reading frame; M and D read 3 bases an amino acid, I none;
    a phase-0 intron (N) is not in the coding sequence, a phase-1 or phase-2 one (U, V: "n nucleotides and 1 amino
    acid", miniprot's manual) holds the 3 bases of the codon it splits (until 2026-11-06 left out: each such intron
    put the frameshifts after it 3 bases early)."""
    out, at = [], 0
    for n, op in re.findall(r"(\d+)([MIDFGNUV])", cigar):
        if op in "MD":
            at += 3 * int(n)
        elif op in "UV":
            at += 3
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
    confirmed by another's: false duplicates). Its stop_codon (GENCODE's rule: in the mRNA, not the CDS) is coding, its
    last exon's (until 2026-11-06 left out: every copy ending in a stop held a 3-base "intron" at its end)."""
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
        elif f[2] in ("CDS", "stop_codon") and cur is not None:
            cur["exons"].append((int(f[3]), int(f[4])))
    for a in out:
        merged = []
        for x, y in sorted(a["exons"]):  # a stop codon touching its exon is that exon's end
            if merged and x <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], y))
            else:
                merged.append((x, y))
        a["exons"] = merged
    return out


def best_of_each(alignments):
    """miniprot's --outs and --outn over alignments of several runs (a genome in pieces: placement.chunks) -> each
    query's of OUTS of its best score or more, its OUTN best (until 2026-11-06 applied in each piece alone: a weak
    paralogue, its protein's best in its piece, was kept)."""
    by_query = collections.defaultdict(list)
    for a in alignments:
        by_query[a["query"]].append(a)
    out = []
    for hits in by_query.values():
        top = max(a["score"] for a in hits)
        out += sorted((a for a in hits if a["score"] >= OUTS * top), key=lambda a: -a["score"])[:OUTN]
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


def hmm_shards(hmm_path, n):
    """hmm_path's models in n files (every n-th model each; fewer where it holds fewer models), made once beside it
    (each whole or not there: several runs may make them at once) -> their paths; [hmm_path] for n of 1."""
    if n <= 1:
        return [hmm_path]
    models, cur = [], []
    with open(hmm_path) as fh:
        for line in fh:
            cur.append(line)
            if line.startswith("//"):  # a model's last line
                models.append("".join(cur).rstrip("\n") + "\n")
                cur = []
    n = min(n, len(models))
    if n <= 1:
        return [hmm_path]
    paths = [f"{hmm_path}.{k + 1}of{n}" for k in range(n)]
    for k, path in enumerate(paths):
        if os.path.isfile(path) and os.path.getmtime(path) >= os.path.getmtime(hmm_path):
            continue
        part = f"{path}.{os.getpid()}.partial"
        with open(part, "w") as fh:
            fh.writelines(models[k::n])
        os.replace(part, path)
    return paths


def hmmsearch_runs(hmm_path, target, out_prefix, threads=1, binary="hmmsearch", options=()):
    """hmmsearch of target against hmm_path's models, in `threads` runs of a share of the models each at once (since
    2026-11-08: one run on every thread kept one core busy, its hmmsearch's reading of the models and the targets
    the work; 4 runs of a quarter: 8.0 s -> 2.6 s for Drosophila's 664 candidates). A model's scores and E-values
    do not depend on the other models (the targets are the same: Z, domZ), so that the hits are those of one run
    -> the runs' --domtblout files (out_prefix.<k>)."""
    shards = hmm_shards(hmm_path, threads)
    cpu = max(1, threads // len(shards))

    def one(k):
        dom = f"{out_prefix}.{k}"
        run = subprocess.run([binary, "--cpu", str(cpu), "--noali", "-o", os.devnull, *options, "--domtblout", dom,
                              shards[k], target], capture_output=True, text=True, check=False)
        if run.returncode != 0:
            raise RuntimeError(f"hmmsearch failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        return dom
    if len(shards) == 1:
        return [one(0)]
    with concurrent.futures.ThreadPoolExecutor(len(shards)) as pool:
        doms = list(pool.map(one, range(len(shards))))
    # in one run's order (by model, as in hmm_path; a model's lines as its run gave them), so that ties go as before
    order = model_order(hmm_path)
    lines = sorted(domtbl_lines(doms), key=lambda line: order.get(line.split(None, 4)[3], len(order)))
    with open(out_prefix, "w") as fh:
        fh.writelines(lines)
    return [out_prefix]


_ORDER = {}


def model_order(hmm_path):
    """{model name: its place in hmm_path} (read once per file and its time)."""
    key = (hmm_path, os.path.getmtime(hmm_path))
    if key not in _ORDER:
        with open(hmm_path) as fh:
            names = [line.split()[1] for line in fh if line.startswith("NAME ")]
        if len(_ORDER) > 8:
            _ORDER.clear()
        _ORDER[key] = {name: i for i, name in enumerate(names)}
    return _ORDER[key]


def domtbl_lines(paths):
    """The hit lines of --domtblout files."""
    for path in paths:
        with open(path) as fh:
            yield from (line for line in fh if not line.startswith("#"))


def hmmsearch(proteins, hmm_db, binary="hmmsearch", threads=1, tmp=None):
    """Each protein's scores against the HMMs -> {protein: {hmm: (bit score, HMM coverage)}}; HMM coverage: the union
    of its domains' hmm_from..hmm_to over the HMM's length."""
    if not proteins:
        return {}
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        faa = os.path.join(work, "c.faa")
        with open(faa, "w") as fh:
            fh.writelines(f">{n}\n{s.replace('*', 'X')}\n" for n, s in proteins.items() if s)
        score, spans, qlen = {}, collections.defaultdict(list), {}
        for line in domtbl_lines(hmmsearch_runs(hmm_db, faa, os.path.join(work, "dom"), threads, binary)):
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


def align(target, queries, miniprot, threads, table, work, options=None):
    """miniprot's alignments of queries {name: protein} on a genome or a windows file (in pieces where large) ->
    (alignments, why it failed or None). options: its intron and seeding options (default -I: the max intron 3.6
    times the square root of the target's length)."""
    faa = os.path.join(work, "q.faa")
    with open(faa, "w") as fh:
        fh.writelines(f">{q}\n{p}\n" for q, p in queries.items())
    alignments, parts = [], placement.chunks(target, work)
    for part in parts:
        run = subprocess.run([miniprot, "-t", str(threads)] + list(options or ["-I"]) + [
            "-T", str(table or 1), "--gff", "--trans", f"--outs={OUTS}", f"--outn={OUTN}", part, faa],
            capture_output=True, text=True, check=False)
        if run.returncode != 0:
            return [], "miniprot failed: " + (run.stderr.strip().splitlines() or ["?"])[-1][:200]
        alignments += parse_alignments(run.stdout)
    return (best_of_each(alignments) if len(parts) > 1 else alignments), None


def confirm(alignments, lengths, odb, hmm_db, hmmsearch_bin, threads, work):
    """The alignments' regions the HMMs confirm (steps 2-3) -> {locus: [each region's confirmed alignment of the best
    HMM score, with hmm_score and hmm_coverage]}; lengths: {query: its protein's length}."""
    found = regions(alignments, lengths)
    candidates = {f"c{i}": (j, a) for i, (j, a) in enumerate(
        ((locus, j), a) for locus, rs in found.items() for j, r in enumerate(rs) for a in r["top"])}
    scores = hmmsearch({n: a["protein"] for n, (_j, a) in candidates.items()}, hmm_db, hmmsearch_bin, threads, work)
    by_region = {}  # (locus, region) -> its confirmed alignment of the best HMM score
    for name, (key, r) in candidates.items():
        got = odb.by_locus(scores.get(name, {}))
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


def _records(lines, block=1 << 24):
    """A FASTA's records (header, its lines joined) from a stream (read in blocks) or an iterable of lines."""
    blocks = iter(lambda: lines.read(block), "") if hasattr(lines, "read") else ["".join(lines)]
    rest = ""
    for got in blocks:
        text = rest + got
        cut = text.rfind("\n>")
        if cut < 0:
            rest = text
            continue
        rest = text[cut + 1:]
        for rec in text[:cut].split("\n>"):
            yield rec
    if rest.strip():
        yield from rest.split("\n>")


def _field(desc, key):
    i = desc.find(key)
    if i < 0:
        return None
    j = desc.find(" ", i)
    return desc[i + len(key):j if j >= 0 else None]


def split_orfs(lines, out, prefix=""):
    """esl-translate's FASTA (lines: a stream or its lines) written to out as orf_pieces' pieces, each on one line:
    >{prefix}name.k source=.. coords=.. length=.. (since 2026-11-08 read in blocks, an ORF without an assembly gap
    written as it is)."""
    gap, write = "X" * GAP_X, out.write
    for rec in _records(lines):
        header, _, body = rec.partition("\n")
        name, _, desc = header.lstrip(">").strip().partition(" ")
        name = prefix + name
        if not name:
            continue
        seq = body.replace("\n", "").replace(" ", "").replace("\r", "")
        source, coords = _field(desc, "source="), _field(desc, "coords=")
        a, b = (int(v) for v in coords.split(".."))
        if gap not in seq:  # one piece, the ORF itself (orf_pieces')
            n = len(seq)
            if ORF_MIN <= n <= HMMSEARCH_MAX:
                write(f">{name}.0 source={source} coords={a}..{a + (3 * n - 1 if a <= b else 1 - 3 * n)} "
                      f"length={n}\n{seq}\n")
            continue
        for k, (p, x, y) in enumerate(orf_pieces(seq, a, b)):
            write(f">{name}.{k} source={source} coords={x}..{y} length={len(p)}\n{p}\n")


HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SPLIT = ("import sys; sys.path.insert(0, sys.argv[1]); from eukore.loci import split_orfs; "
         "split_orfs(sys.stdin, sys.stdout, sys.argv[2])")


def orf_file(genome, orfs, table, threads, work, translate="esl-translate"):
    """A genome's six-frame ORFs (esl-translate, split_orfs) into orfs, in `threads` pieces of its contigs at once,
    each esl-translate's and its split's own processes (since 2026-11-08: one esl-translate, its ORFs split by this
    process on one core, 11.6 s for Drosophila's 4.6 million), joined in the genome's order: the ORFs and their
    order those of one run, their names each piece's (p<k>.)."""
    pieces = [genome]
    if threads > 1:
        folder = os.path.join(work, "orf_pieces")
        os.makedirs(folder, exist_ok=True)
        pieces = placement.chunks(genome, folder, max(1, placement.genome_bases(genome) / threads))
    outs = [f"{orfs}.{k}" for k in range(len(pieces))]

    def one(k):
        err = f"{outs[k]}.err"
        with open(outs[k], "w") as fh, open(err, "w") as eh:
            tr = subprocess.Popen([translate, "-l", str(ORF_MIN), "-c", str(table or 1), pieces[k]],
                                  stdout=subprocess.PIPE, stderr=eh)
            split = subprocess.Popen([sys.executable, "-c", SPLIT, HERE, f"p{k}."], stdin=tr.stdout, stdout=fh,
                                     stderr=eh)
            tr.stdout.close()
            split.wait()
            tr.wait()
        if tr.returncode or split.returncode:
            with open(err) as eh:
                raise RuntimeError(f"esl-translate failed (exit {tr.returncode}, its split's {split.returncode}): "
                                   f"{eh.read().strip()[-200:]}")
        os.remove(err)
        if pieces[k] != genome:
            os.remove(pieces[k])
    with concurrent.futures.ThreadPoolExecutor(len(pieces)) as pool:
        list(pool.map(one, range(len(pieces))))
    with open(orfs, "wb") as fh:
        for path in outs:
            with open(path, "rb") as part:
                shutil.copyfileobj(part, fh, 1 << 24)
            os.remove(path)


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
        orf_file(part, orfs, table, threads, work, translate)
        found = collections.defaultdict(list)  # ORF -> [(locus, score, hmm_from, hmm_to, env_from, env_to)]
        for line in domtbl_lines(hmmsearch_runs(hmm, orfs, dom, threads, hmmsearch_bin,
                                                ["--domE", str(RESCUE_DOM_E)])):
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
    contigs = placement.read_spans(genome_path, {c: [(s - RESCUE_FLANK, e + RESCUE_FLANK) for _l, s, e in w]
                                                 for c, w in by_contig.items()})
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
    intron = int(3.6 * placement.genome_bases(genome_path) ** 0.5)
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


def annotated_transcripts(chosen, annotation):
    """The annotation's transcripts sharing coding bases with the calls' spans (annotation_mod.read_cds), chosen:
    {locus: alignment}; annotation: (its GFF, NCBI's assembly report or None, its source)."""
    from . import annotation as annotation_mod  # noqa: PLC0415
    gff, report = (tuple(annotation) + (None,))[:2]
    spans = collections.defaultdict(list)
    for a in chosen.values():
        spans[a["contig"]].append((a["start"], a["end"]))
    return annotation_mod.read_cds(gff, spans, annotation_mod.aliases(report, gff))


def annotated(chosen, contigs, annotation, odb, hmmsearch_bin, threads, work, hmm_db, table, transcripts=None):
    """Step 6: the genome's own annotated gene model (annotation: (its GFF, NCBI's assembly report or None, its
    source)) in place of a call's, where the GFF holds a transcript on the call's contig and strand sharing coding
    bases with it (the most: annotation_mod.best_transcript) whose translation the locus' own HMM confirms (its best
    HMM, at the cutoff) at ANNOTATION_SHARE of the call's HMM score or more. chosen: {locus: alignment} -> the same, an
    annotated model with model=the source (annotation.SOURCES'), the transcript, and miniprot's model kept as
    miniprot=dict(start, end, exons) (annotation.py's benchmark compares that). The annotation's exons are its CDS: its 5' piece's phase left out (a partial gene).
    transcripts: annotated_transcripts' (read here where not given)."""
    from . import annotation as annotation_mod  # noqa: PLC0415
    from .gencode import translate  # noqa: PLC0415 (until 2026-11-08 through benchmarks' alias: a worker loaded it)
    source = (tuple(annotation) + ("ncbi",) * 3)[2]
    if transcripts is None:
        transcripts = annotated_transcripts(chosen, annotation)
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
        cds = copies_mod.revcomp(cds) if t["strand"] == "-" else cds
        models[locus] = dict(exons=exons, start=exons[0][0], end=exons[-1][1], transcript=name)
        proteins[locus] = translate(cds, table).rstrip("*")
    scores = hmmsearch(proteins, hmm_db, hmmsearch_bin, threads, work)
    out = dict(chosen)
    for locus, m in models.items():
        got = odb.by_locus(scores.get(locus, {}))
        if not got or max(got, key=lambda h: got[h][0]) != locus:
            continue
        s, cov = got[locus]
        a = chosen[locus]
        if s >= odb.cutoff[locus] and s >= ANNOTATION_SHARE * a["hmm_score"]:
            out[locus] = dict(a, exons=m["exons"], start=m["start"], end=m["end"], shifts=[], hmm_score=s,
                              hmm_coverage=cov, protein=proteins[locus], model=source, transcript=m["transcript"],
                              miniprot=dict(start=a["start"], end=a["end"], exons=a["exons"]))
    return out


def call(genome_path, queries, odb, miniprot="miniprot", hmmsearch_bin="hmmsearch", threads=1, tmp=None,
         hmm_db=None, table=1, rescue_missing=True, translate="esl-translate", annotation=None):
    """The BUSCO genes of a genome (table: its NCBI genetic code) -> dict(loci {locus: dict(status, copy, spliced,
    features [(kind, start, end)],
    contig, start, end, strand, identity, coverage, score, hmm_score, hmm_coverage, cutoff, copies, query, protein,
    rescued, model (miniprot, or the annotation's source: ncbi, ensembl, toga2), miniprot (its model: start, end, exons),
    others (the coding sequences of the locus' other complete copies, since 2026-11-07: dup_identity))}, failed (why the genome could not be
    called, or None), windows {name: sequence} of the copies' flanks).
    rescue_missing: the rescue pass (step 5) of the loci left Missing, in a genome of RESCUE_MAX_BASES or less;
    annotation: (the genome's GFF, NCBI's assembly report or None, its source) for step 6 (annotated()), or None."""
    lengths = {q: len(p) for q, p in queries.items()}
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        hmm_db = hmm_db or odb.hmm_db(os.path.join(work, "loci.hmm"))
        alignments, failed = align(genome_path, queries, miniprot, threads, table, work)
        if failed:
            return dict(loci={}, failed=failed, windows={})
        confirmed = confirm(alignments, lengths, odb, hmm_db, hmmsearch_bin, threads, work)
        missing = [locus for locus in odb.loci if not confirmed.get(locus)]
        if (rescue_missing and missing and placement.genome_bases(genome_path) <= RESCUE_MAX_BASES
                and shutil.which(translate)):
            confirmed.update(rescue(genome_path, missing, queries, odb, miniprot, hmmsearch_bin, threads, table,
                                    hmm_db, work, translate))
        chosen, statuses = {}, {}
        for locus in odb.loci:
            rs = confirmed.get(locus, [])
            complete = [r for r in rs if r["hmm_coverage"] >= COMPLETE]
            if rs:
                chosen[locus] = max(complete or rs, key=lambda r: (r["hmm_score"], r["score"]))
                statuses[locus] = ("Single" if len(complete) == 1 else "Duplicated" if complete else "Fragmented",
                                   [r for r in complete if r is not chosen[locus]])
        transcripts = annotated_transcripts(chosen, annotation) if annotation and annotation[0] and chosen else None
        # the genome's bases read: each copy's and its other copies' spans, the annotated models', WINDOW each side
        # (since 2026-11-08 these alone: every contig with a call was held whole, a chromosome-level genome whole)
        need = collections.defaultdict(list)
        for locus, best in chosen.items():
            for r in [best] + statuses[locus][1]:
                need[r["contig"]].append((r["start"] - WINDOW, r["end"] + WINDOW))
        for t in (transcripts or {}).values():
            need[t["contig"]].append((t["cds"][0][0] - WINDOW, max(b for _a, b in t["cds"]) + WINDOW))
        contigs = placement.read_spans(genome_path, need)
        if transcripts is not None:
            chosen = annotated(chosen, contigs, annotation, odb, hmmsearch_bin, threads, work, hmm_db, table,
                               transcripts)
    loci, windows = {}, {}
    for locus in odb.loci:
        if locus not in chosen:
            loci[locus] = dict(status="Missing")
            continue
        best = chosen[locus]
        status, others = statuses[locus]
        seq = contigs[best["contig"]]
        fitted = copies_mod.whole(seq, best["start"], best["end"], best["strand"], best["exons"])
        raw = seq[best["start"] - 1:best["end"]]
        feats = list(fitted["features"])
        feats += frameshift_features([(s, e) for k, s, e in feats if k == "exon"], best.get("shifts"))
        feats += [("soft_masked", a, b) for a, b in on_copy(soft_runs(raw), best["start"], best["end"], best["strand"])]
        feats += [("low_complexity", a, b) for a, b in dust(fitted["copy"])]
        loci[locus] = dict(status=status, copy=fitted["copy"], spliced=in_frame(fitted["spliced"], best.get("shifts")),
                           features=sorted(
            feats, key=lambda f: (f[1], f[0])), contig=best["contig"], contig_length=len(seq), start=best["start"],
            end=best["end"], strand=best["strand"], identity=best["identity"], coverage=best["coverage"],
            score=best["score"], hmm_score=best["hmm_score"], hmm_coverage=best["hmm_coverage"],
            cutoff=odb.cutoff[locus], copies=len(others) + (status != "Fragmented"),
            others=[in_frame(copies_mod.whole(contigs[o["contig"]], o["start"], o["end"], o["strand"],
                                              o["exons"])["spliced"], o.get("shifts")) for o in others],
            query=best["query"].split("~", 1)[1],
            protein=best["protein"], exons=best["exons"], rescued=bool(best.get("rescued")),
            model=best.get("model", "miniprot"), transcript=best.get("transcript", "-"),
            miniprot=best.get("miniprot") or dict(start=best["start"], end=best["end"], exons=best["exons"]))
        w1, w2 = max(1, best["start"] - WINDOW), min(len(seq), best["end"] + WINDOW)
        windows[f"{locus}|{best['contig']}|{w1}-{w2}"] = seq[w1 - 1:w2].upper()
    return dict(loci=loci, failed=None, windows=windows)


DUP_MIN_ID = 0.5  # a duplicated copy's other copies compared down to this identity (below it: "<50")


def dup_identity(loci, vsearch, threads=1, tmp=None):
    """How alike a Duplicated locus' copies are (since 2026-11-07): the identity (vsearch, --iddef 2, %) of the nearest
    other complete copy's coding sequence to the chosen copy's -> {locus: "99.8", or "<50" where none is DUP_MIN_ID
    identical}; {} without vsearch. Uncollapsed haplotypes of a diploid's assembly are 99% and more alike, a whole
    genome duplication's copies less (package.duplication_of)."""
    dups = {locus: r for locus, r in loci.items() if r.get("others")}
    if not dups or not vsearch:
        return {}
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        query, db, out = (os.path.join(work, n) for n in ("others.fna", "chosen.fna", "hits.tsv"))
        with open(query, "w") as fq, open(db, "w") as fd:
            for locus, r in dups.items():
                fd.write(f">{locus}\n{r['spliced']}\n")
                fq.writelines(f">{locus}|{i}\n{s}\n" for i, s in enumerate(r["others"]) if s)
        run = subprocess.run([vsearch, "--usearch_global", query, "--db", db, "--id", str(DUP_MIN_ID), "--iddef", "2",
                              "--maxaccepts", "0", "--maxrejects", "0", "--userout", out,
                              "--userfields", "query+target+id", "--threads", str(max(1, threads)), "--quiet"],
                             capture_output=True, text=True, check=False)
        if run.returncode != 0 or not os.path.isfile(out):
            return {}  # not measured (not "unlike")
        best = {}
        with open(out) as fh:
            for line in fh:
                q, t, ident = line.rstrip("\n").split("\t")[:3]
                if q.split("|")[0] == t:
                    best[t] = max(best.get(t, 0.0), float(ident))
    return {locus: f"{best[locus]:.1f}" if locus in best else f"<{DUP_MIN_ID * 100:.0f}" for locus in dups}
