"""EUKDB's rRNA (since 2026-11-05; from protal's genome_18s.py and rdna_freeze.py, made complete and KSGP-free): every
18S span and complete rDNA operon of the build's genomes and of the published long reads (EUKARYOME 2.1, ROD v1.2),
labelled by NCBI taxid, checked without any KSGP.

- A genome's units: windows around conserved 18S sites (SEEDS) reaching DOWN bases past them on the 18S's strand;
  barrnap (--kingdom euk) finds the 18S, 5.8S, 28S and 5S; an 18S cut by introns is one span. A unit is the 18S span's
  start to its 28S's end, uncut (protal cut its records at 8,191 bases), with its 5.8S; ITSx gives the ITS1/5.8S/ITS2
  borders and the 28S's start where installed. A 5S in the same block (within GAP_5S of the 28S's end or of the 18S's
  start, on either strand, as in some fungi's and protists' repeats) extends the unit and is one of its regions; a 5S
  elsewhere in the genome is not taken. `lsu_complete`: barrnap's 28S hit is not partial.
- Published records: ROD's operon variants by their own NCBI taxid; EUKARYOME's long reads by their species name
  matched to NCBI's names and synonyms (else their genus's), its compromised-quality ids left out. Their regions as a
  genome's units'.
- Checks, each a column (a record failing one is kept and flagged): ambiguous bases, the longest homopolymer, the 99%
  18S type within its species, and own_species: a genome's 18S types against every other source's 18S (vsearch),
  agreeing to the family or deeper with its NCBI lineage (t), sharing no more than the phylum at FOREIGN_ID or more
  (f, a contaminant's rDNA), else undecided (-).
"""

import collections
import gzip
import os
import re
import shutil
import subprocess
import tempfile

from . import util

MOTIFS = ("GTGCCAGCAGCCGCGGTAA", "AAACTTAAAGGAATTGACGG", "GGTGGTGCATGGCCGTTCTTAGTT")  # 18S ~570, ~1,140, ~1,200
SEED = 12
UP = 3000           # a window's bases before a seed on its 18S's strand
DOWN = 15000        # and after it: the whole 28S (human's ends ~9.5 kb past the 18S's start) and an inline 5S
MAX_INTRON = 1500   # an 18S cut by introns: pieces this close are one span
MIN_LEN = 1200      # an 18S span this long at least
MAX_ITS = 3000      # the 18S's end to the 5.8S, the 5.8S's end to the 28S: at most
GAP_5S = 3000       # an inline 5S: this close to the unit's ends
SAME_TYPE = 0.99
NEAR = 0.5
FOREIGN_ID = 97.0
PREFIX = 1000       # a published read's first bases: its 18S, for the own-species pool
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ITSX_FIELD = re.compile(r"(SSU|ITS1|5\.8S|ITS2|LSU): (\d+)-(\d+)")
RECORD_COLUMNS = ("record", "kind", "taxid", "species_taxid", "accession", "source", "regions", "length",
                  "lsu_complete", "ambiguous", "homopolymer", "type", "own_species", "screen", "label_from")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


SEEDS = sorted({s[i:i + SEED] for m in MOTIFS for s in (m, revcomp(m)) for i in range(len(m) - SEED + 1)})
FORWARD = {m[i:i + SEED] for m in MOTIFS for i in range(len(m) - SEED + 1)}


class BarrnapError(RuntimeError):
    pass


def barrnap_gff(barrnap, fasta, threads=1):
    command = [barrnap, "--kingdom", "euk", "--threads", str(threads), "--quiet", fasta]
    run = subprocess.run(command, capture_output=True, text=True, check=False)
    if run.returncode != 0:
        raise BarrnapError(f"barrnap failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
    return run.stdout


def contigs(path):
    name, parts = None, []
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = (line[1:].split() or [""])[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def windows(seq):
    """Merged [start, end) windows around the seeds in an upper-case sequence."""
    sites = []
    for seed in SEEDS:
        before, after = (UP, DOWN) if seed in FORWARD else (DOWN, UP)
        i = seq.find(seed)
        while i >= 0:
            sites.append((max(0, i - before), min(len(seq), i + after)))
            i = seq.find(seed, i + 1)
    spans = []
    for a, b in sorted(sites):
        if spans and a <= spans[-1][1]:
            spans[-1][1] = max(spans[-1][1], b)
        else:
            spans.append([a, b])
    return spans


def parse_gff(gff):
    """barrnap's GFF -> {window: {feature name: [(start, end, strand, partial)]}} (1-based inclusive)."""
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for line in gff.splitlines():
        f = line.split("\t")
        if line.startswith("#") or len(f) < 9:
            continue
        name = f[8].split("Name=", 1)[1].split(";", 1)[0] if "Name=" in f[8] else ""
        out[f[0]][name].append((int(f[3]), int(f[4]), f[6], "aligned only" in f[8] or "partial" in f[8]))
    return out


def joined(features, gap=MAX_INTRON):
    out = []
    for strand in "+-":
        for start, end, _s, _p in sorted(f for f in features if f[2] == strand):
            if out and out[-1][2] == strand and start - out[-1][1] <= gap:
                out[-1][1] = max(out[-1][1], end)
            else:
                out.append([start, end, strand])
    return out


def units_of(seq, feats):
    """A window's (or a read's) 18S spans and their units -> [dict(span (start, end, strand), sequence (the span,
    its strand), exons (the 18S pieces joined), unit (sequence or None), regions [(name, start, end)] on the unit,
    lsu_complete, its stretch for ITSx and where that starts)]. Coordinates 1-based in `seq`."""
    out = []
    s18 = feats.get("18S_rRNA", [])
    for start, end, strand in joined(s18):
        pieces = [f for f in s18 if f[2] == strand and start <= f[0] and f[1] <= end]
        span = seq[start - 1:end]
        exons = "".join(seq[a - 1:b] for a, b, _s, _p in sorted(pieces))
        if strand == "-":
            span, exons = revcomp(span), revcomp(exons)
        item = dict(span=(start, end, strand), sequence=span, exons=exons, unit=None, regions=None,
                    lsu_complete=False)
        fwd = strand == "+"
        s58 = [f for f in feats.get("5_8S_rRNA", []) if f[2] == strand and
               (0 < f[0] - end <= MAX_ITS if fwd else 0 < start - f[1] <= MAX_ITS)]
        if s58:
            f58 = min(s58, key=lambda f: f[0]) if fwd else max(s58, key=lambda f: f[1])
            lsu = [f for f in feats.get("28S_rRNA", []) if f[2] == strand and
                   ((f[1] > f58[1] + 100 and f[0] <= f58[1] + MAX_ITS) if fwd else
                    (f[0] < f58[0] - 100 and f[1] >= f58[0] - MAX_ITS))]
            if lsu:
                f28 = min(lsu, key=lambda f: f[0]) if fwd else max(lsu, key=lambda f: f[1])
                a, b = (start, f28[1]) if fwd else (f28[0], end)
                five = [f for f in feats.get("5S_rRNA", [])  # inline: in the same block, either strand
                        if (fwd and (0 <= f[0] - b <= GAP_5S or 0 <= a - f[1] <= GAP_5S)) or
                        (not fwd and (0 <= a - f[1] <= GAP_5S or 0 <= f[0] - b <= GAP_5S))]
                a = min([a] + [f[0] for f in five])
                b = max([b] + [f[1] for f in five])
                unit = seq[a - 1:b]
                if not fwd:
                    unit = revcomp(unit)

                def at(x1, x2):  # a feature's [start, end) on the unit
                    return (x1 - a, x2 - a + 1) if fwd else (b - x2, b - x1 + 1)

                regions = [("18S",) + at(start, end), ("5.8S",) + at(f58[0], f58[1]), ("28S",) + at(f28[0], f28[1])]
                regions += [("5S",) + at(f[0], f[1]) for f in five]
                r18, r58 = regions[0], regions[1]
                regions.insert(1, ("ITS1", r18[2], r58[1]))
                regions.insert(3, ("ITS2", r58[2], regions[3][1]))
                item.update(unit=unit, regions=sorted(regions, key=lambda r: r[1]), lsu_complete=not f28[3],
                            its=unit[max(0, r18[2] - 300):min(len(unit), r58[2] + 3300)],
                            its_at=max(0, r18[2] - 300))
        out.append(item)
    return out


def run_itsx(stretches, itsx, threads=1):
    """{stretch: ITSx's positions {region: (begin, end) 1-based}} for each distinct stretch."""
    if not stretches or not itsx:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta, prefix = os.path.join(tmp, "its.fasta"), os.path.join(tmp, "itsx")
        names = {}
        with open(fasta, "w") as fh:
            for s in stretches:
                if s not in names:
                    names[s] = f"u{len(names)}"
                    fh.write(f">{names[s]}\n{s}\n")
        subprocess.run([itsx, "-i", fasta, "-o", prefix, "-t", "all", "--cpu", str(max(1, threads)), "--save_regions",
                        "none", "--preserve", "T", "--complement", "F", "--graphical", "F", "--summary", "F"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        found = {}
        if os.path.isfile(prefix + ".positions.txt"):
            with open(prefix + ".positions.txt") as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    found[f[0]] = {m.group(1): (int(m.group(2)), int(m.group(3)))
                                   for m in (ITSX_FIELD.match(x) for x in f[2:]) if m}
        return {s: found.get(n, {}) for s, n in names.items()}


def refine(item, found):
    """A unit's ITS1/5.8S/ITS2/28S borders from ITSx where it found ITS1 and ITS2 (barrnap's otherwise)."""
    if not item.get("unit") or not all(k in found for k in ("ITS1", "ITS2")):
        return item
    at = item["its_at"]
    its1 = (at + found["ITS1"][0] - 1, at + found["ITS1"][1])
    its2 = (at + found["ITS2"][0] - 1, at + found["ITS2"][1])
    if its1[1] >= its2[0]:
        return item
    regions = {r[0]: r for r in item["regions"]}
    end28 = regions["28S"][2]
    new = [r for r in item["regions"] if r[0] not in ("18S", "ITS1", "5.8S", "ITS2", "28S")]
    new += [("18S", regions["18S"][1], its1[0]), ("ITS1",) + its1, ("5.8S", its1[1], its2[0]), ("ITS2",) + its2,
            ("28S", its2[1], end28)]
    item["regions"] = sorted(new, key=lambda r: r[1])
    return item


def regions_text(regions):
    return ",".join(f"{n}:{a}-{b}" for n, a, b in regions) if regions else "-"


def extract_genome(path, barrnap, itsx=None, threads=1):
    """A genome's rRNA -> [dict(kind 18S or rDNA, contig, start, end, strand, sequence, exons, regions, lsu_complete)]:
    each 18S span (MIN_LEN at least), and each unit with its 5.8S and 28S as an rDNA record, uncut."""
    seqs, origin = {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        win = os.path.join(tmp, "windows.fna")
        with open(win, "w") as fh:
            for name, seq in contigs(path):
                seq = seq.upper()
                for a, b in windows(seq):
                    wid = f"w{len(seqs)}"
                    seqs[wid], origin[wid] = seq[a:b], (name, a)
                    fh.write(f">{wid}\n{seq[a:b]}\n")
        gff = barrnap_gff(barrnap, win, threads) if seqs else ""
    feats = parse_gff(gff)
    items = []
    for wid, f in feats.items():
        for item in units_of(seqs[wid], f):
            contig, offset = origin[wid]
            item.update(contig=contig, start=offset + item["span"][0], end=offset + item["span"][1],
                        strand=item["span"][2])
            items.append(item)
    found = run_itsx([i["its"] for i in items if i.get("unit")], itsx, threads)
    rows = []
    for item in items:
        if len(item["sequence"]) >= MIN_LEN:
            rows.append(dict(kind="18S", contig=item["contig"], start=item["start"], end=item["end"],
                             strand=item["strand"], sequence=item["sequence"], exons=item["exons"], regions=None,
                             lsu_complete=False))
        if item.get("unit"):
            refine(item, found.get(item["its"], {}))
            rows.append(dict(kind="rDNA", contig=item["contig"], start=item["start"], end=item["end"],
                             strand=item["strand"], sequence=item["unit"], exons=item["exons"],
                             regions=item["regions"], lsu_complete=item["lsu_complete"]))
    return rows


def published_regions(reads, barrnap, itsx=None, threads=1):
    """Published long reads [(name, sequence)] -> {name: (regions, lsu_complete, 18S exons)}: each read as a window."""
    if not reads:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "reads.fna")
        with open(fasta, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in reads)
        gff = barrnap_gff(barrnap, fasta, threads)
    feats = parse_gff(gff)
    seqs = dict(reads)
    items = {}
    for name, f in feats.items():
        units = [u for u in units_of(seqs[name], f) if u.get("unit")]
        if units:
            items[name] = max(units, key=lambda u: len(u["unit"]))
    found = run_itsx([u["its"] for u in items.values()], itsx, threads)
    return {n: (refine(u, found.get(u["its"], {}))["regions"], u["lsu_complete"], u["exons"]) for n, u in items.items()}


# --- labels by NCBI taxid ----------------------------------------------------------------------------------------

def name_index(ncbi):
    """{lower-case name: taxid} of NCBI's scientific names and synonyms (a name two taxa share: none)."""
    index, clash = {}, set()
    for taxid, name in ncbi.name.items():
        for n in [name] + ncbi.synonyms.get(taxid, []):
            key = n.lower()
            if key in index and index[key] != taxid:
                clash.add(key)
            index[key] = taxid
    for key in clash:
        index.pop(key, None)
    return index


def eukaryome_taxid(label, index, ncbi):
    """A EUKARYOME label (k__...;g__Genus;s__Genus_species) -> (taxid, how) by its species name, else its genus's;
    (None, None) where neither is NCBI's."""
    parts = {p[:1]: p[3:].replace("_", " ") for p in label.split(";") if "__" in p}
    sp = parts.get("s", "")
    if sp and (t := index.get(sp.lower())) and ncbi.rank.get(ncbi.species_of(t)) == "species":
        return ncbi.species_of(t), "eukaryome_name"
    genus = parts.get("g", "")
    if genus and (t := index.get(genus.lower())) and ncbi.rank.get(t) == "genus":
        return t, "eukaryome_genus"
    return None, None


# --- checks ------------------------------------------------------------------------------------------------------

def homopolymer(seq):
    return max((len(m.group(0)) for m in re.finditer(r"(.)\1*", seq)), default=0)


def ambiguous(seq):
    return sum(1 for c in seq if c not in "ACGT") / max(1, len(seq))


def identity(a, b):
    """1 - edits / the shorter's length, the shorter aligned to its best stretch of the longer (Myers' bit-vector)."""
    if len(a) > len(b):
        a, b = b, a
    if not a:
        return 0.0
    m = len(a)
    peq = {}
    for i, c in enumerate(a):
        peq[c] = peq.get(c, 0) | (1 << i)
    mask, high = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score, best = mask, 0, m, m
    for c in b:
        eq = peq.get(c, 0)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (mask & ~(xh | pv))
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph, mh = (ph << 1) & mask, (mh << 1) & mask
        pv, mv = mh | (mask & ~(xv | ph)), ph & xv
        best = min(best, score)
    return 1 - best / m


def types_of(exons):
    """Sequences -> {sequence: type number}, each the first type whose first sequence it is SAME_TYPE identical to,
    the commonest first."""
    reps, out = [], {}
    for seq, _n in sorted(collections.Counter(exons).items(), key=lambda kv: (-kv[1], -len(kv[0]), kv[0])):
        for i, rep in enumerate(reps):
            if identity(seq, rep) >= SAME_TYPE:
                out[seq] = i + 1
                break
        else:
            reps.append(seq)
            out[seq] = len(reps)
    return out


RANKS = ("phylum", "class", "order", "family", "genus")


def ssu_part(seq, row):
    """A record's 18S: an 18S record whole, an operon's 18S region, else (a read without one found) its first PREFIX."""
    if row["kind"] == "18S":
        return seq
    return next((seq[a:b] for n, a, b in row.get("regions_list") or () if n == "18S"), seq[:PREFIX])


def own_species(records, ncbi, vsearch, threads=1):
    """The own-species check (module docstring) -> {record: t, f or -}: each genome's 18S types against the 18S of
    every other genome and source (vsearch --usearch_global at 80%, the first 32 hits), their nearest (within NEAR
    points of the best) agreeing with the genome's NCBI lineage to the family or deeper: t; sharing no more than the
    phylum at FOREIGN_ID or more: f; else -. records: [(name, sequence, row)] (row: accession, source, species_taxid,
    regions); a published record is not checked (its label is its source's)."""
    if not vsearch or not records:
        return {}
    ssu = {name: ssu_part(seq, row)[:2500] for name, seq, row in records}
    owner = {name: row["accession"] for name, _s, row in records}
    queries = {name: ssu[name] for name, _s, row in records if row["source"] == "genome"}
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        db, q, hits = (os.path.join(tmp, f) for f in ("db.fna", "q.fna", "hits.tsv"))
        with open(db, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in ssu.items() if s)
        with open(q, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in queries.items() if s)
        if not queries:
            return {}
        subprocess.run([vsearch, "--usearch_global", q, "--db", db, "--id", "0.8", "--iddef", "2", "--maxaccepts",
                        "32", "--maxrejects", "256", "--userout", hits, "--userfields", "query+target+id",
                        "--threads", str(max(1, threads)), "--quiet"], check=True)
        found = collections.defaultdict(list)
        with open(hits) as fh:
            for line in fh:
                qn, tn, ident = line.rstrip("\n").split("\t")
                if owner.get(tn) != owner.get(qn):
                    found[qn].append((float(ident), tn))
    taxid_of = {name: int(row["species_taxid"]) for name, _s, row in records}
    for name in queries:
        hits = found.get(name, [])
        if not hits:
            out[name] = "-"
            continue
        top = max(i for i, _t in hits)
        mine = ncbi.ranks(taxid_of[name])
        depth = 0
        for ident, target in hits:
            if ident < top - NEAR:
                continue
            theirs = ncbi.ranks(taxid_of[target])
            depth = max([depth] + [i + 1 for i, r in enumerate(RANKS) if mine.get(r) and mine.get(r) == theirs.get(r)])
        out[name] = "t" if depth >= 4 else "f" if depth <= 1 and top >= FOREIGN_ID else "-"
    return out


def write_records(folder, records):
    """records [(name, sequence, row)] -> folder/records.fna.gz and records.tsv (RECORD_COLUMNS)."""
    os.makedirs(folder, exist_ok=True)
    fasta, table = os.path.join(folder, "records.fna.gz"), os.path.join(folder, "records.tsv")
    with gzip.open(util.part_of(fasta), "wt", compresslevel=6) as fa, open(util.part_of(table), "w", newline="\n") as tb:
        tb.write("\t".join(RECORD_COLUMNS) + "\n")
        for name, seq, row in records:
            fa.write(f">{name}\n{seq}\n")
            tb.write("\t".join(str(row.get(c, "-")) if row.get(c) not in (None, "") else "-" for c in RECORD_COLUMNS)
                     + "\n")
    os.replace(util.part_of(fasta), fasta)
    os.replace(util.part_of(table), table)


def finish(records, ncbi, vsearch=None, threads=1):
    """The checks over all records [(name, sequence, row)] (row: kind, accession, source, species_taxid, taxid,
    regions_list, lsu_complete, label_from) -> the same with their RECORD_COLUMNS filled, sorted by species and name."""
    by_species = collections.defaultdict(list)
    for name, seq, row in records:
        by_species[row["species_taxid"]].append((name, ssu_part(seq, row)))
    type_of = {}
    for sp, items in by_species.items():
        types = types_of([p for _n, p in items])
        for name, part in items:
            type_of[name] = types[part]
    own = own_species(records, ncbi, vsearch, threads)
    out = []
    for name, seq, row in records:
        row = dict(row, record=name, length=len(seq), regions=regions_text(row.get("regions_list")),
                   lsu_complete="t" if row.get("lsu_complete") else "f", ambiguous=f"{ambiguous(seq):.4f}",
                   homopolymer=homopolymer(seq), type=type_of.get(name, "-"), own_species=own.get(name, "-"),
                   screen=row.get("screen", "-"))
        out.append((name, seq, row))
    return sorted(out, key=lambda r: (int(r[2]["species_taxid"]), r[0]))


def genome_records(accession, rows, taxid, species_taxid, source="genome"):
    """A genome's extraction (extract_genome's rows) -> records [(name, sequence, row)]."""
    out, n = [], collections.Counter()
    for r in rows:
        n[r["kind"]] += 1
        name = f"{accession}_{r['kind']}_{n[r['kind']]}"
        out.append((name, r["sequence"], dict(kind=r["kind"], accession=accession, source=source, taxid=taxid,
                                              species_taxid=species_taxid, regions_list=r["regions"],
                                              lsu_complete=r["lsu_complete"], label_from="assembly")))
    return out


def save(folder, rows):
    """A genome's extraction into its done folder (rrna.tsv: kind, contig, start, end, strand, regions,
    lsu_complete, sequence)."""
    with open(os.path.join(folder, "rrna.tsv"), "w", newline="\n") as fh:
        fh.write("kind\tcontig\tstart\tend\tstrand\tregions\tlsu_complete\tsequence\n")
        for r in rows:
            fh.write(f"{r['kind']}\t{r['contig']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{regions_text(r['regions'])}"
                     f"\t{'t' if r['lsu_complete'] else 'f'}\t{r['sequence']}\n")
    return {"18S": sum(r["kind"] == "18S" for r in rows), "rDNA": sum(r["kind"] == "rDNA" for r in rows)}


def load(folder):
    """A genome's saved extraction -> extract_genome's rows ([] where there is none)."""
    path = os.path.join(folder, "rrna.tsv")
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            regions = None if f[5] == "-" else [(p.rpartition(":")[0], int(p.rpartition(":")[2].split("-")[0]),
                                                 int(p.rpartition(":")[2].split("-")[1])) for p in f[5].split(",")]
            rows.append(dict(kind=f[0], contig=f[1], start=int(f[2]), end=int(f[3]), strand=f[4], regions=regions,
                             lsu_complete=f[6] == "t", sequence=f[7]))
    return rows


def has_tools(barrnap="barrnap"):
    return bool(shutil.which(barrnap))
