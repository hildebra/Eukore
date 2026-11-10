"""The read benchmark (since 2026-11-06): how well reads of a eukaryote are told apart by the package, by any read
aligner, with no profiler in between (benchmarks.py's reads.tsv; docs/benchmarks.md).

Test genomes: up to PER_GROUP species' genomes per supergroup (choose: those with the most other species of their
genus in the package first, so that a held-out species has relatives to be placed by), kept on disk through the
build's prune. From each, simulated reads (READ_LEN, ERROR substitutions, either strand):
- marker reads: from its BUSCO copies' and rRNA records' spans in its assembly, introns and all, COVERAGE-fold: the
  reads a reference should place;
- genome reads: GENOME_READS from anywhere in its genome: the reads it should mostly leave alone (off-target);
and BACKGROUND_READS from each of a few bacterial and archaeal genomes (data/read_background.tsv), which any
assignment gets wrong.

The references, each aligned to by minimap2 (-x sr, secondary alignments kept):
- Eukore: the BUSCO copies whole and the rRNA records;
- Eukore exons: the copies' coding sequences alone (an exon-level reference, EukDetect's kind) and the rRNA records;
- Eukore 18S: the rRNA records' 18S alone (an 18S reference, PR2's kind);
- other references (compare: EukDetect2's database, MetaPhlAn 4's eukaryote markers, PR2's 18S; tools/
  compare_references.py makes the latter two), each's records labelled by its own taxid table, cut to the package's
  clades in a trial package.

A read's hits covering MIN_COVER of it, MIN_ID identical or more and within NEAR of its best hit name the lowest
common ancestor of their taxa: its assignment. Each read is assigned three times: with its own species in the
reference (`in`), with its own species' records left out (`held_out`: a species new to the reference, to be placed by
its relatives), and with its whole order's left out (`novel`, since 2026-11-09: a lineage no genome holds, right at its
class at best). Per rank (species .. phylum) an assignment is correct (its taxon there is the read's), wrong, or above
it (an assignment of a higher rank, which is right but says less). A background read is a prokaryote's (wrong wherever
placed) or, given its taxid (since 2026-11-09: human chromosome 21, the host), a eukaryote's, scored as a test's.
"""

import bisect
import collections
import os
import random
import shutil
import subprocess

from . import labels, package, qc, rrna, util
from . import odb as odb_mod
from .benchmarks import rows_of as benchmarks_rows  # a package table's rows, streamed
from .genomes import genome_file
from .placement import contig_lengths
from .rrna import revcomp

READ_LEN = 150
ERROR = 0.005          # substitutions per base
COVERAGE = 3           # marker reads: each locus' bases read this many times
GENOME_READS = 20_000  # reads of a test genome from anywhere in it
BACKGROUND_READS = 20_000
PER_GROUP = 3          # test genomes per supergroup (eukore_build.py --read_genomes)
MIN_ID = 0.90
MIN_COVER = 0.8
NEAR = 0.01
RANKS = ("species", "genus", "family", "order", "class", "phylum")
# novel (since 2026-11-09): the test species' whole order out of the reference, a lineage the reference holds none of
# (a metagenome's unknown eukaryote): its reads are right at its class at best, wrong at its order and below
NOVEL_RANK = "order"
MODES = ("in", "held_out", "novel")
REFERENCES = ("Eukore", "Eukore exons", "Eukore 18S")
COLUMNS = ("test", "kind", "species_taxid", "supergroup", "congeners", "reference", "mode", "reads", "n",
           "assigned") + tuple(f"{r}_{x}" for r in RANKS for x in ("correct", "wrong"))
# detection (detection.tsv): each test genome a sample, its reads from anywhere in it (the first `depth`) and as many
# background reads (depth / their genomes of each); a taxon is detected with MIN_READS reads placed to it at the rank
# scored (the species with it in the reference, `in`; its genus with it held out) and MIN_SHARE of the sample's reads
# placed at that rank, as profilers cut by relative abundance: a fixed count alone detects more congeners the deeper
# the sample (the 2026-11-06 trial: Eukore's false species 4, 7, 20 in 24 samples of 1,000, 5,000, 20,000 reads)
MIN_READS = 3
MIN_SHARE = 0.01
DEPTHS = (1_000, 5_000, 20_000)
DETECT_RANK = {"in": "species", "held_out": "genus", "novel": "class"}
DETECT_COLUMNS = ("test", "species_taxid", "supergroup", "congeners", "reference", "mode", "rank", "depth",
                  "in_reference", "true_reads", "detected", "false_taxa")
BASES = "ACGT"
# minimap2's secondary alignments kept a read; a species' rRNA records are in the references one per 99% type (its
# copies of an rDNA array, its 18S records and its units' 18S alike), so that a held-out species' own records do not
# take them all (until 2026-11-06 a species of 51 records or more left its held-out reads no other hit)
MAX_SECONDARY = 50


def genus_of(ncbi, taxid):
    return next((t for t in ncbi.ancestors(ncbi.merged.get(taxid, taxid)) if ncbi.rank.get(t) == "genus"), None)


def choose(ncbi, reps, per_group=PER_GROUP, exists=os.path.isfile):
    """The test genomes: reps {accession: (species taxid, genome path)} -> [accession], up to per_group per supergroup,
    those whose genus holds the most other species first, then by accession."""
    genus = {a: genus_of(ncbi, sp) for a, (sp, _p) in reps.items()}
    size = collections.Counter(g for g in genus.values() if g)
    by_group = collections.defaultdict(list)
    for a, (sp, path) in reps.items():
        if path and exists(path):
            by_group[labels.group(ncbi, sp)].append(a)
    out = []
    for group in sorted(by_group):
        ranked = sorted(by_group[group], key=lambda a: (-(size[genus[a]] - 1 if genus[a] else 0), a))
        out += ranked[:per_group]
    return out


# --- references --------------------------------------------------------------------------------------------------

def compared(compare):
    """--read_compare's references as a list [(name, FASTA, table)]: one (a tuple), several (a list) or none."""
    if not compare:
        return []
    return [tuple(compare)] if isinstance(compare[0], str) else [tuple(c) for c in compare]


def references(folder, out, ncbi, compare=None, inside=None):
    """The references' FASTAs written into out -> {name: (path, {record: taxid})}. compare: other references, each
    (name, FASTA, table of record<TAB>taxid) (compared); inside: a test of a taxid (a trial package's clades) the
    compared references are cut to."""
    os.makedirs(out, exist_ok=True)
    species_of = {r["accession"]: int(r["species_taxid"])
                  for r in package.read_table(folder, "genomes.tsv") if r["representative"] == "t"}
    exons = collections.defaultdict(list)
    for f in benchmarks_rows(os.path.join(folder, "busco", "features.tsv")):
        if f["kind"] in ("exon", "frameshift"):  # qc.spliced leaves the frameshifts' bases out
            exons[(f["accession"], f["locus"])].append((f["kind"], int(f["start"]), int(f["end"])))
    loci = [r["locus"] for r in benchmarks_rows(os.path.join(folder, "busco", "loci.tsv"))]
    paths = {name: os.path.join(out, name.replace(" ", "_") + ".fna") for name in REFERENCES}
    taxa = {name: {} for name in REFERENCES}
    with open(paths["Eukore"], "w") as whole, open(paths["Eukore exons"], "w") as cds, \
            open(paths["Eukore 18S"], "w") as ssu:
        for locus in loci:
            copies = os.path.join(folder, "busco", "copies", f"{locus}.fna.gz")
            if not os.path.isfile(copies):
                continue
            for acc, seq in package.read_copies(folder, locus).items():
                if acc not in species_of:
                    continue
                name = f"{acc}|{locus}"
                whole.write(f">{name}\n{seq}\n")
                taxa["Eukore"][name] = species_of[acc]
                coding = qc.spliced(seq, exons[(acc, locus)])
                if coding:
                    cds.write(f">{name}\n{coding}\n")
                    taxa["Eukore exons"][name] = species_of[acc]
        rows = {r["record"]: r for r in benchmarks_rows(os.path.join(folder, "rrna", "records.tsv"))}
        path = os.path.join(folder, "rrna", "records.fna.gz")
        typed = set()  # (species, 99% type): one record of each (MAX_SECONDARY)
        for name, seq in (odb_mod.read_fasta(path) if os.path.isfile(path) else ()):
            r = rows.get(name)
            if not r:
                continue
            if r.get("type", "-") != "-":
                if (r["species_taxid"], r["type"]) in typed:
                    continue
                typed.add((r["species_taxid"], r["type"]))
            taxid = int(r["taxid"])
            for ref, fh in (("Eukore", whole), ("Eukore exons", cds)):
                fh.write(f">{name}\n{seq}\n")
                taxa[ref][name] = taxid
            regions = rrna.parse_regions(r["regions"])
            part = seq if r["kind"] == "18S" else next((seq[a:b] for n, a, b in regions or () if n == "18S"), "")
            if part:
                ssu.write(f">{name}\n{part}\n")
                taxa["Eukore 18S"][name] = taxid
    out_refs = {name: (paths[name], taxa[name]) for name in REFERENCES}
    for name, fasta, table in compared(compare):
        wanted = {}
        with util.open_text(table) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1].isdigit() and (inside is None or inside(int(f[1]))):
                    wanted[f[0].split()[0]] = int(f[1])
        path = os.path.join(out, name.replace(" ", "_") + ".fna")
        with open(path, "w") as fh:
            for rec, seq in odb_mod.read_fasta(fasta):
                if rec in wanted:
                    fh.write(f">{rec}\n{seq}\n")
        out_refs[name] = (path, wanted)
    return out_refs


# --- reads -------------------------------------------------------------------------------------------------------

def read_of(contig, pos, rng):
    """READ_LEN bases of contig from pos (0-based), ERROR substitutions, on a random strand; None at an end or in Ns."""
    seq = contig[pos:pos + READ_LEN].upper()
    if len(seq) < READ_LEN or seq.count("N") > READ_LEN // 10:
        return None
    seq = "".join(rng.choice(BASES.replace(c, "")) if c in BASES and rng.random() < ERROR else c for c in seq)
    return revcomp(seq) if rng.random() < 0.5 else seq


def marker_spans(folder, accessions):
    """Genomes' BUSCO copies' and rRNA records' spans in their assemblies -> {accession: [(contig, start, end)]}
    (1-based), the tables read once for them all (until 2026-11-06 once a test genome: 3.7 million rows each)."""
    wanted, spans = set(accessions), collections.defaultdict(list)
    for r in benchmarks_rows(os.path.join(folder, "busco", "copies.tsv")):
        if r["accession"] in wanted and r["start"].isdigit():
            spans[r["accession"]].append((r["contig"], int(r["start"]), int(r["end"])))
    for r in benchmarks_rows(os.path.join(folder, "rrna", "records.tsv")):
        if r.get("accession") in wanted and r.get("start", "-").isdigit():
            spans[r["accession"]].append((r["contig"], int(r["start"]), int(r["end"])))
    return spans


def simulate(genome, spans, rng, n_genome=GENOME_READS):
    """A genome's reads -> (marker reads [sequence], genome reads [sequence], in random order: detection's samples are
    their first reads, until 2026-11-06 those of the first contigs in the file)."""
    lengths = contig_lengths(genome)
    ends, total = [], 0  # each contig's end in the genome joined (a pick by bisect: a genome of a million contigs)
    for _c, n in lengths:
        total += n
        ends.append(total)
    picks = collections.defaultdict(list)
    for _ in range(n_genome if total > READ_LEN else 0):
        x = rng.randrange(total)
        k = bisect.bisect_right(ends, x)
        picks[lengths[k][0]].append(x - (ends[k] - lengths[k][1]))
    by_contig = collections.defaultdict(list)
    for contig, start, end in spans:
        by_contig[contig].append((start, end))
    marker, anywhere = [], []
    with util.open_text(genome) as fh:
        name, parts = None, []

        def take(name, seq):
            for start, end in by_contig.get(name, ()):
                lo, hi = max(0, start - 1 - READ_LEN // 2), max(0, end - READ_LEN // 2)
                for _ in range(max(1, round((end - start + 1) * COVERAGE / READ_LEN))):
                    read = read_of(seq, rng.randint(lo, max(lo, hi)), rng)
                    if read:
                        marker.append(read)
            for pos in picks.get(name, ()):
                read = read_of(seq, min(pos, max(0, len(seq) - READ_LEN)), rng)
                if read:
                    anywhere.append(read)
        for line in fh:
            if line.startswith(">"):
                if name in by_contig or name in picks:
                    take(name, "".join(parts))
                name, parts = (line[1:].split() or [""])[0], []
            elif name in by_contig or name in picks:
                parts.append(line.strip())
        if name in by_contig or name in picks:
            take(name, "".join(parts))
    rng.shuffle(anywhere)
    return marker, anywhere


# --- assignment ----------------------------------------------------------------------------------------------------

def paf_groups(path):
    """A minimap2 PAF's alignments per read (consecutive) -> (read, [(identity, target)]) of those covering
    MIN_COVER of the read."""
    read, hits = None, []
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 11:
                continue
            if f[0] != read:
                if read is not None:
                    yield read, hits
                read, hits = f[0], []
            qlen, qs, qe, matches, block = int(f[1]), int(f[2]), int(f[3]), int(f[9]), int(f[10])
            if block and (qe - qs) >= MIN_COVER * qlen:
                hits.append((matches / block, f[5]))
    if read is not None:
        yield read, hits


def assign(hits, taxa, ncbi, leave_out=None, clade=None, cache=None):
    """A read's hits [(identity, target)] -> the taxon of their lowest common ancestor, or None; leave_out: a species
    whose records are out of the reference (held out); clade: (rank, taxid) whose records are out of it (novel)."""
    if leave_out is not None:
        hits = [(i, t) for i, t in hits if t in taxa and ncbi.species_of(taxa[t]) != leave_out]
    if clade is not None:
        cache = {} if cache is None else cache
        hits = [(i, t) for i, t in hits if t in taxa and ranks_of(ncbi, taxa[t], cache).get(clade[0]) != clade[1]]
    hits = [(i, t) for i, t in hits if t in taxa]
    if not hits:
        return None
    best = max(i for i, _t in hits)
    if best < MIN_ID:
        return None
    return qc.lca(ncbi, {taxa[t] for i, t in hits if i >= best - NEAR})


def ranks_of(ncbi, taxid, cache):
    """{rank: taxid} of a taxon's lineage at RANKS, and the deepest rank index it reaches."""
    if taxid not in cache:
        at = {}
        for t in ncbi.ancestors(ncbi.merged.get(taxid, taxid)):
            r = ncbi.rank.get(t)
            if r in RANKS and r not in at:
                at[r] = t
        cache[taxid] = at
    return cache[taxid]


def score(truth, taxon, ncbi, cache):
    """Per rank, an assignment (taxon) of a read of species truth: correct, wrong or above (None: unscorable)."""
    mine, theirs = ranks_of(ncbi, truth, cache), ranks_of(ncbi, taxon, cache)
    deepest = min((RANKS.index(r) for r in theirs), default=len(RANKS))
    out = {}
    for i, r in enumerate(RANKS):
        if r not in mine:
            out[r] = None
        elif i < deepest:
            out[r] = "above"
        else:
            out[r] = "correct" if theirs.get(r) == mine[r] else "wrong" if r in theirs else "above"
    return out


# --- detection -----------------------------------------------------------------------------------------------------

def holdings(taxa, ncbi, cache):
    """A reference's taxa {record: taxid} -> (its species, {genus: its species}, {class: its orders})."""
    species = {ncbi.species_of(t) for t in set(taxa.values())} - {None}
    genera, orders = collections.defaultdict(set), collections.defaultdict(set)
    for sp in species:
        at = ranks_of(ncbi, sp, cache)
        if at.get("genus"):
            genera[at["genus"]].add(sp)
        if at.get("class") and at.get(NOVEL_RANK):
            orders[at["class"]].add(at[NOVEL_RANK])
    return species, genera, orders


def in_reference(sp, mode, held, ncbi, cache):
    """Whether a reference (holdings) can place a test species right: holds it (in), another species of its genus
    (held out), another order of its class (novel)."""
    species, genera, orders = held
    if mode == "in":
        return sp in species
    at = ranks_of(ncbi, sp, cache)
    if mode == "novel":
        return bool(at.get("class") and orders.get(at["class"], set()) - {at.get(NOVEL_RANK)})
    g = at.get("genus")
    return bool(g and genera.get(g, set()) - {sp})


def detect(placed, test, sp, rank, depth, background, ncbi, cache, present=()):
    """A sample: the test's first `depth` reads placed ({source: {read index: taxon}}) and depth / len(background) of
    each background genome's -> (reads placed to its right taxon at rank, whether it is detected, other taxa
    detected): detected with MIN_READS reads and MIN_SHARE of the reads placed at rank. present: taxa at rank the
    sample truly holds besides the test (a eukaryote background's, the host: since 2026-11-09), no false detection."""
    counts = collections.Counter()
    per = depth // max(1, len(background))
    for source, upto in [(test, depth)] + [(b, per) for b in background]:
        for i, taxon in placed.get(source, {}).items():
            if i < upto:
                t = ranks_of(ncbi, taxon, cache).get(rank)
                if t:
                    counts[t] += 1
    floor = max(MIN_READS, MIN_SHARE * sum(counts.values()))
    truth = ranks_of(ncbi, sp, cache).get(rank)
    right = counts.pop(truth, 0) if truth else 0
    return right, right >= floor, sum(1 for t, n in counts.items() if n >= floor and t not in present)


def detection_summary(rows):
    """detection.tsv's headline: per reference, mode and depth over its samples (held out: those with congeners)
    the right taxa detected (tp), missed with it not in the reference (fn_absent) or in it (fn_missed), wrong taxa
    detected (fp), precision, recall and F1."""
    sums = collections.defaultdict(collections.Counter)
    for r in rows:
        if r["mode"] == "held_out" and not int(r["congeners"]):
            continue
        s = sums[(r["reference"], r["mode"], int(r["depth"]))]
        s["samples"] += 1
        if r["detected"] == "t":
            s["tp"] += 1
        else:
            s["fn_missed" if r["in_reference"] == "t" else "fn_absent"] += 1
        s["fp"] += int(r["false_taxa"])
    out = collections.defaultdict(dict)
    for (ref, mode, depth), s in sorted(sums.items()):
        tp, fn, fp = s["tp"], s["fn_missed"] + s["fn_absent"], s["fp"]
        p = tp / (tp + fp) if tp + fp else None
        rc = tp / (tp + fn) if tp + fn else None
        out[ref][f"{mode}_{depth}"] = dict(s, precision=round(p, 4) if p is not None else None,
                                           recall=round(rc, 4) if rc is not None else None,
                                           f1=round(2 * tp / (2 * tp + fp + fn), 4) if tp + fp + fn else None)
    return dict(out)


# --- all of it -----------------------------------------------------------------------------------------------------

def run(folder, ncbi, genomes_folder, tests, minimap2, threads, work, compare=None, background=(), inside=None,
        log=print, seed=1, n_genome=GENOME_READS, n_background=BACKGROUND_READS, detections=None):
    """The read benchmark of the package in folder (tests: accessions of genomes in genomes_folder; background:
    [(name, FASTA)]) in work (deleted after) -> reads.tsv's rows; detections: a list given takes detection.tsv's rows
    (each test genome a sample at DEPTHS, detect)."""
    os.makedirs(work, exist_ok=True)
    rng = random.Random(seed)
    genomes = {r["accession"]: r for r in package.read_table(folder, "genomes.tsv") if r["representative"] == "t"}
    species = {int(r["taxid"]): dict(r, supergroup=r.get("group") or r.get("supergroup", "-"))  # format 5: its group
               for r in package.read_table(folder, "species.tsv")}
    genus_size = collections.Counter(genus_of(ncbi, int(g["species_taxid"])) for g in genomes.values())
    refs = references(folder, os.path.join(work, "refs"), ncbi, compare, inside)
    reads_fa = os.path.join(work, "reads.fa")
    meta = {}  # test -> (kind, species taxid)
    counts = collections.Counter()  # (test, reads) -> reads simulated
    spans = marker_spans(folder, tests)
    with open(reads_fa, "w") as fh:
        for acc in tests:
            sp = int(genomes[acc]["species_taxid"])
            marker, anywhere = simulate(genome_file(genomes_folder, acc), spans.get(acc, []), rng, n_genome)
            meta[acc] = ("eukaryote", sp)
            for kind, reads in (("marker", marker), ("genome", anywhere)):
                counts[(acc, kind)] = len(reads)
                fh.writelines(f">{acc}|{kind}|{i}\n{s}\n" for i, s in enumerate(reads))
            log(f"  reads of {acc}: {len(marker):,} from its markers, {len(anywhere):,} from anywhere")
        for name, path, *taxid in background:  # (name, FASTA[, NCBI taxid: a eukaryote's, since 2026-11-09])
            _m, anywhere = simulate(path, [], rng, n_background)
            meta[name] = ("background", int(taxid[0]) if taxid and taxid[0] else 0)
            counts[(name, "genome")] = len(anywhere)
            fh.writelines(f">{name}|genome|{i}\n{s}\n" for i, s in enumerate(anywhere))
    cache = {}
    # novel: each test's order, out of the reference (a test without an order has no novel mode)
    novel = {acc: ranks_of(ncbi, sp, cache).get(NOVEL_RANK) for acc, (what, sp) in meta.items() if what == "eukaryote"}
    rows = []
    for ref, (fasta, taxa) in refs.items():
        if not taxa:
            continue
        paf, err = os.path.join(work, "hits.paf"), os.path.join(work, "minimap2.err")
        # --split-prefix: a reference over minimap2's index part (-I, 8 Gb: the full package's) is indexed in parts,
        # and a read's hits of all parts merged, consecutive (paf_groups; until 2026-11-06 a read was counted once a
        # part, its assignment the last part's)
        with open(paf, "w") as out, open(err, "w") as eh:
            run = subprocess.run([minimap2, "-x", "sr", "-c", "--secondary=yes", "-N", str(MAX_SECONDARY), "-p", "0.5",
                                  "--split-prefix", os.path.join(work, "mm2split"), "-t", str(max(1, threads)), fasta,
                                  reads_fa], stdout=out, stderr=eh)
        if run.returncode:
            with open(err) as eh:
                raise RuntimeError(f"minimap2 failed (exit {run.returncode}) on {ref}: {eh.read().strip()[-300:]}")
        os.remove(err)
        tally = collections.defaultdict(collections.Counter)
        placed = {mode: collections.defaultdict(dict) for mode in MODES}  # genome reads: {source: {index: taxon}}
        for read, hits in paf_groups(paf):
            test, kind, i = read.rsplit("|", 2)
            what, sp = meta[test]
            for mode in MODES:
                clade = None
                if mode == "novel" and what == "eukaryote":
                    if not novel.get(test):
                        continue
                    clade = (NOVEL_RANK, novel[test])
                taxon = assign(hits, taxa, ncbi, sp if mode == "held_out" and what == "eukaryote" else None, clade,
                               cache)
                if taxon is None:
                    continue
                if kind == "genome":
                    placed[mode][test][int(i)] = taxon
                t = tally[(test, kind, mode)]
                t["assigned"] += 1
                if sp:  # a test's read, or a eukaryote background's (the host's: right on its own lineage)
                    for r, verdict in score(sp, taxon, ncbi, cache).items():
                        if verdict in ("correct", "wrong"):
                            t[f"{r}_{verdict}"] += 1
                else:  # a bacterium's read: wrong at every rank its assignment reaches (Eukaryota reaches none)
                    reached = ranks_of(ncbi, taxon, cache)
                    for r in RANKS[min((RANKS.index(x) for x in reached), default=len(RANKS)):]:
                        t[f"{r}_wrong"] += 1
        os.remove(paf)
        for (test, kind), n in sorted(counts.items()):
            what, sp = meta[test]
            for mode in MODES:
                if (what == "background" and mode != "in") or (mode == "novel" and not novel.get(test)):
                    continue
                t = tally[(test, kind, mode)]
                g = genus_of(ncbi, sp) if sp else None
                rows.append(dict(test=test, kind=what, species_taxid=sp or "-",
                                 supergroup=species.get(sp, {}).get("supergroup", "-") if sp else "-",
                                 congeners=genus_size[g] - 1 if g else 0, reference=ref, mode=mode, reads=kind, n=n,
                                 assigned=t["assigned"], **{c: t[c] for c in COLUMNS[10:]}))
        if detections is not None:
            held = holdings(taxa, ncbi, cache)
            names = [b[0] for b in background]
            hosts = [meta[n][1] for n in names if meta[n][1]]  # the eukaryote background's: truly in every sample
            depths = [d for d in DEPTHS if d <= n_genome] or [n_genome]
            for test in tests:
                sp = meta[test][1]
                g = genus_of(ncbi, sp)
                for mode in MODES:
                    if mode == "novel" and not novel.get(test):
                        continue
                    present = {ranks_of(ncbi, h, cache).get(DETECT_RANK[mode]) for h in hosts} - {None}
                    for depth in depths:
                        right, found, wrong = detect(placed[mode], test, sp, DETECT_RANK[mode], depth, names, ncbi,
                                                     cache, present)
                        detections.append(dict(
                            test=test, species_taxid=sp, supergroup=species.get(sp, {}).get("supergroup", "-"),
                            congeners=genus_size[g] - 1 if g else 0, reference=ref, mode=mode,
                            rank=DETECT_RANK[mode], depth=depth,
                            in_reference="t" if in_reference(sp, mode, held, ncbi, cache) else "f",
                            true_reads=right, detected="t" if found else "f", false_taxa=wrong))
        log(f"  reads aligned to {ref} ({len(taxa):,} records)")
    for path in [p for p, _t in refs.values()] + [reads_fa]:
        if os.path.isfile(path):
            os.remove(path)
    return rows


def background_rows(table):
    """data/read_background.tsv's rows -> [(accession, assembly, organism, NCBI taxid or None)]: a eukaryote's taxid
    (since 2026-11-09) makes its reads right on its own lineage (the host's), a prokaryote's are wrong wherever
    placed."""
    out = []
    if not os.path.isfile(table):
        return out
    with open(table) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 3 or f[0] == "accession":
                continue
            out.append((f[0], f[1], f[2], int(f[3]) if len(f) > 3 and f[3].strip().isdigit() else None))
    return out


def read_background(table, folder):
    """data/read_background.tsv's genomes in folder (fetched by the fetch stage) -> [(name, path, taxid or None)] of
    those there."""
    out = []
    for acc, _assembly, organism, taxid in background_rows(table):
        path = os.path.join(folder, f"{acc}.fna.gz")
        if os.path.isfile(path):
            out.append((organism.replace(" ", "_"), path, taxid))
    return out


def background_url(accession, assembly):
    """An NCBI genome's FTP URL from its accession and assembly name (GCF_000005845.2, ASM584v2); a sequence's
    (NC_000021.9: one chromosome, since 2026-11-09) E-utilities FASTA."""
    if not accession.startswith(("GCF_", "GCA_")):
        return f"{EUTILS}/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    digits = accession.split("_")[1].split(".")[0]
    full = f"{accession}_{assembly}"
    return (f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[:3]}/{digits[:3]}/{digits[3:6]}/{digits[6:9]}/"
            f"{full}/{full}_genomic.fna.gz")


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_background(table, folder, log=print):
    """data/read_background.tsv's genomes into folder (kept when there); a sequence's plain FASTA compressed there."""
    os.makedirs(folder, exist_ok=True)
    for acc, assembly, organism, _taxid in background_rows(table):
        path = os.path.join(folder, f"{acc}.fna.gz")
        if os.path.isfile(path):
            continue
        plain = not acc.startswith(("GCF_", "GCA_"))
        got = path[:-len(".gz")] + ".part" if plain else path
        why = util.fetch_file(background_url(acc, assembly), got, gzip_check=not plain)
        if why:
            log(f"WARNING: read background {acc} ({organism}): {why}")
        elif plain:
            with open(got, "rb") as src, util.bgzf_open(util.part_of(path), "wb") as out:
                shutil.copyfileobj(src, out, 1 << 24)
            os.replace(util.part_of(path), path)
            os.remove(got)
