"""Eukore's benchmarks (since 2026-11-06): how complete, how right and how discriminating a package is. The package
stage measures them as it builds (into the package's benchmarks/, before its manifest, which keeps summary.json's
figures as `benchmarks`); any package can be measured again afterwards:

    python3 -m eukore.benchmarks OUT/Eukore_2026.11 --build_dir WORK -t 8    # into OUT/Eukore_2026.11/benchmarks/

docs/benchmarks.md describes them. Each is a table of benchmarks/:
- completeness.tsv: per species' genome its BUSCOs Single, Duplicated, Fragmented and Missing, its copies'
  confidence, its assembly's level and contigs, its rRNA in the assembly: BUSCO's completeness, on Eukore's calls;
- orthodb.tsv: per BUSCO of each package species that BUSCO's dataset holds a protein of (OrthoDB's gene model in
  that species), Eukore's call against it: missed, or the share of OrthoDB's protein's 5-mers in the call's
  translation (recall) and of the translation's in OrthoDB's (precision). The right gene model of the right paralog
  gives both near 1. The build's calls are queried with the nearest dataset species' proteins, for these species
  their own (`guided`), so that they find OrthoDB's own protein again; the `held_out_` columns compare the same
  species' calls made without its own proteins (genomes.HELD_OUT, made by the process stage since 2026-11-06, read
  from the work folder): the measure that is not circular, of a species whose nearest dataset relative is another;
- resolution.tsv: per marker (each BUSCO's coding sequence; the 18S of each species' best rRNA record) and the rank a
  copy shares with its nearest copy of another species (qc.py's vsearch search), their identity's quantiles and the
  copies under 97% and 99%: how well a marker tells species apart, and the identity a read can expect per rank;
- rrna.tsv: per rRNA source and kind its records holding each region, a whole 28S, an 18S intron, each confidence;
- labels.tsv: per supergroup and reference (PR2, EUKARYOME, KSGP) its species by how they were named;
- runtime.tsv: per genome its size, threads and seconds (WORK/done/<accession>/done.json, genomes run since
  2026-11-06), given the work folder;
- stages.tsv: the package stage's steps and their seconds (the build's own run);
- composition.tsv (since 2026-11-10): per supergroup what the package holds of its species (a genome, rRNA alone,
  complete BUSCOs, an operon, an organelle gene);
- summary.json: the headline figures; summary.md (since 2026-11-10) the same as Markdown tables.
The paper's tables of the package alone (paper(): selection, breadth, overlap, annotation, genes, holdout, ..) are
described in docs/benchmarks.md. Since 2026-11-12 the comparisons with other references (EukDetect2, MetaPhlAn 4,
PR2), the read benchmark and every figure are benchmarks/'s, run on a finished package (benchmarks/README.md): the
build makes the database and its own checks alone.
A benchmark whose inputs are missing (BUSCO's dataset, vsearch, the work folder) is left out, and summary.json says so.
"""

import argparse
import collections
import functools
import json
import os
import re
import shutil
import statistics
import sys
import time

from . import odb as odb_mod
from . import genomes as genomes_mod
from . import gencode, package, qc, rrna, store, util

FOLDER = "benchmarks"
KMER = 5
CONCORDANT = 0.9    # recall and precision both this high: the call is OrthoDB's gene model
SIMILAR = 0.5       # both this high: the same gene, another model (exon ends, an isoform, another assembly)
DEPTHS = ("none", "phylum", "class", "order", "family", "genus")  # qc.depth's 0 .. 5: the rank shared
QUANTILES = (("q05", 0.05), ("q25", 0.25), ("median", 0.5), ("q75", 0.75), ("q95", 0.95))
SSU_MIN = 1000      # an 18S this long at least, for resolution.tsv (a read's fragment would understate identity)
REGIONS = ("18S", "ITS1", "5.8S", "ITS2", "28S", "5S")
TIERS = package.TIERS  # high, unverified, medium, low (format 6; a format-5 package has no unverified: 0)
MATCHES = ("taxid", "species", "synonym", "genus", "-")  # taxid: MIDORI2's own (since format 5)
REFERENCES = package.REFERENCES
ALL_BUSCOS = "BUSCO (all)"

COMPLETENESS_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "source", "level", "genome_size",
                        "contigs", "loci", "single", "duplicated", "fragmented", "missing", "complete", "high",
                        "unverified",
                        "medium", "low", "rrna_18s", "rdna_units", "rdna_complete")
ORTHODB_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "locus", "status", "odb_proteins",
                   "odb_length", "protein_length", "recall", "precision", "guided", "verdict")
HELD_COLUMNS = ("held_out_status", "held_out_query", "held_out_recall", "held_out_precision", "held_out_verdict",
                "held_out_place", "held_out_span", "held_out_rescued")
ORTHODB_COLUMNS += HELD_COLUMNS
RESOLUTION_COLUMNS = ("marker", "shared_rank", "copies", "under_97", "under_99") + tuple(q for q, _f in QUANTILES)
RRNA_COLUMNS = ("source", "kind", "records", "median_length") + REGIONS + ("lsu_complete", "introns") + TIERS
LABEL_COLUMNS = ("supergroup", "reference") + tuple("none" if m == "-" else m for m in MATCHES)
RUNTIME_COLUMNS = ("accession", "species_taxid", "supergroup", "genome_size", "genome_bytes", "searched", "threads",
                   "seconds")

# NCBI's genetic codes: every table since 2026-11-07 (eukore/gencode.py; the nuclear ones alone until then)
translate = gencode.translate


def kmers(protein, k=KMER):
    """A protein's k-mers, but those across a stop or an unknown residue."""
    protein = protein.rstrip("*")
    return {w for w in (protein[i:i + k] for i in range(len(protein) - k + 1)) if "*" not in w and "X" not in w}


def shared(a, b):
    """The share of a's k-mers found in b."""
    return len(a & b) / len(a) if a else 0.0


def quantile(values, q):
    """values (sorted) at q, nearest rank."""
    return values[min(len(values) - 1, int(q * len(values)))]


def rows_of(path):
    """A package table's rows (dicts), streamed: copies.tsv and features.tsv of all eukaryotes are millions of rows."""
    if not os.path.isfile(path):
        return
    with util.open_text(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if line.strip():
                yield dict(zip(head, line.rstrip("\n").split("\t")))


def write(folder, name, columns, rows):
    package._tsv(os.path.join(folder, name), columns, rows)  # noqa: SLF001 (the package's own TSV writer)


# --- the benchmarks ----------------------------------------------------------------------------------------------

def completeness(genomes, species, statuses, tiers, n_loci):
    """completeness.tsv's rows: each species' genome (genomes: genomes.tsv's representatives) with its BUSCOs' calls
    (statuses, tiers: {accession: Counter})."""
    out = []
    for acc, g in sorted(genomes.items(), key=lambda kv: int(kv[1]["species_taxid"])):
        s, sp = statuses.get(acc, collections.Counter()), species.get(g["species_taxid"], {})
        called = s["Single"] + s["Duplicated"]
        out.append(dict(accession=acc, species_taxid=g["species_taxid"], name=sp.get("name", g.get("organism")),
                        supergroup=sp.get("supergroup", "-"), source=g.get("source"), level=g.get("level"),
                        genome_size=g.get("genome_size"), contigs=g.get("contigs"), loci=n_loci,
                        single=s["Single"], duplicated=s["Duplicated"], fragmented=s["Fragmented"],
                        missing=n_loci - called - s["Fragmented"], complete=f"{called / max(n_loci, 1):.3f}",
                        **{t: tiers.get(acc, collections.Counter())[t] for t in TIERS},
                        rrna_18s=g.get("rrna_18s"), rdna_units=g.get("rdna_units"),
                        rdna_complete=g.get("rdna_complete")))
    return out


def odb_species(odb, ncbi, held):
    """{package species: {locus: [(dataset taxid, protein)]}} of the dataset's proteins of the package's species
    (a dataset taxid of a strain is its species')."""
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for locus, proteins in odb.proteins.items():
        for taxid, _name, seq in proteins:
            sp = ncbi.species_of(taxid)
            if sp in held:
                out[sp][locus].append((taxid, seq))
    return out


def compare(cds, table, references):
    """A call's coding sequence against the dataset's proteins of its species [(taxid, protein)] -> (recall,
    precision, OrthoDB protein's length, call's protein length, verdict), against the protein that matches best."""
    protein = translate(cds, table)
    mine = kmers(protein)
    recall, precision, length = max((shared(kmers(p), mine), shared(mine, kmers(p)), len(p)) for _t, p in references)
    return recall, precision, length, len(protein.rstrip("*")), verdict(recall, precision)


compare_cds = compare  # run()'s parameter `compare` (another reference) hides the name there


def verdict(recall, precision):
    """A call's verdict from its k-mer recall and precision against its reference protein (the one rule of compare,
    marker_rows and tools/validate_markers.py)."""
    low = min(recall, precision)
    return "concordant" if low >= CONCORDANT else "similar" if low >= SIMILAR else "discordant"


def best_match(protein, refs):
    """A protein against reference proteins {name: protein} -> (the name of the one sharing the most k-mers both ways,
    recall, precision against it)."""
    got = kmers(protein)
    best = max(refs, key=lambda a: shared(kmers(refs[a]), got) + shared(got, kmers(refs[a])))
    ref = kmers(refs[best])
    return best, shared(ref, got), shared(got, ref)


def place_of(call, held):
    """Where a held-out call lies against the package's call of the locus (both dicts with contig, start, end) ->
    (held_out_place, held_out_span): `same` (on the package call's contig, overlapping it) and the two spans'
    intersection over union, `elsewhere` (another contig, or no overlap) and 0, or `-` where either has no place."""
    try:
        a, b = int(held["start"]), int(held["end"])
        c, d = int(call["start"]), int(call["end"])
    except (KeyError, TypeError, ValueError):
        return "-", "-"
    inter = min(b, d) - max(a, c) + 1
    if held.get("contig") != call.get("contig") or inter <= 0:
        return "elsewhere", "0.000"
    return "same", f"{inter / (max(b, d) - min(a, c) + 1):.3f}"


def orthodb_row(acc, species_taxid, locus, call, cds, references, ncbi, held=None):
    """orthodb.tsv's row of one BUSCO of a genome whose species the dataset holds: call (dict: status, query,
    contig, start, end) or None (Missing); cds: the call's coding sequence; references: the dataset's proteins of the
    species [(taxid, protein)]; held: the held-out call (genomes.HELD_OUT's row: status, query, contig, start, end,
    rescued, cds) or None where the genome has none."""
    row = dict(accession=acc, species_taxid=species_taxid, locus=locus, odb_proteins=len(references),
               odb_length=max(len(p) for _t, p in references))
    table = ncbi.code_of(species_taxid)
    if call is None:
        row.update(status="Missing", protein_length="-", recall="-", precision="-", guided="-", verdict="missed")
    else:
        status, query = call["status"], call["query"]
        recall, precision, length, plen, verdict = compare(cds, table, references)
        q_taxid = query.split("_")[1] if query.count("_") >= 2 else ""  # <busco>_<taxid>_<id>; not an ancestral's
        guided = "t" if q_taxid.isdigit() and ncbi.species_of(int(q_taxid)) == species_taxid else "f"
        row.update(status=status, odb_length=length, protein_length=plen, recall=f"{recall:.3f}",
                   precision=f"{precision:.3f}", guided=guided, verdict=verdict)
    if held is None:
        row.update({c: "-" for c in HELD_COLUMNS})
    elif held["status"] == "Missing":
        row.update({c: "-" for c in HELD_COLUMNS}, held_out_status="Missing", held_out_verdict="missed")
    else:
        recall, precision, _length, _plen, verdict = compare(held["cds"], table, references)
        place, span = place_of(call, held) if call else ("-", "-")
        row.update(held_out_status=held["status"], held_out_query=held["query"], held_out_recall=f"{recall:.3f}",
                   held_out_precision=f"{precision:.3f}", held_out_verdict=verdict, held_out_place=place,
                   held_out_span=span, held_out_rescued=held.get("rescued", "-"))
    return row


def held_out_calls(work, acc, rank=None):
    """A genome's held-out calls (WORK/done/<accession>/held_out.tsv; or, rank given, its clade hold-out's at that rank,
    held_out_<rank>.tsv) -> {locus: its row (status, query, cds; contig, start, end and rescued since a later run)},
    None without."""
    folder, name = genomes_mod.done_dir(work, acc) if work else None, genomes_mod.held_out_file(rank)
    if not folder or not store.has(folder, name):
        return None
    return {r["locus"]: r for r in store.table(folder, name)}


HOLDOUT_LEVELS = ("species",) + genomes_mod.HELD_OUT_RANKS
DEPTH_NAMES = ("none", "phylum", "class", "order", "family", "genus")  # qc.depth's 0 .. 5
HOLDOUT_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "level", "left_out", "nearest", "loci", "called",
                   "concordant", "similar", "discordant", "missed", "same_gene", "span_90")


def clade_levels(work, acc):
    """A dataset species' genome's held-out calls at each level (since 2026-11-09: its species', then its genus',
    family's, order's left out; genomes.clade_held_out) -> {level: (calls (held_out_calls'), record (left_out, nearest:
    the depth its lineage shares with the nearest dataset species kept))}, {} without."""
    summary = genomes_mod.done_summary(work, acc) if work else None
    held = (summary or {}).get("held_out")
    out = {}
    if not held:
        return out
    calls = held_out_calls(work, acc)
    if calls:
        out["species"] = (calls, dict(left_out=len(held.get("left_out") or ()), nearest=held.get("nearest")))
    for rank in genomes_mod.HELD_OUT_RANKS:
        rec = (held.get("ranks") or {}).get(rank)
        if not rec or rec.get("no_taxon"):
            continue
        source = rec.get("same_as") or rank
        calls = held_out_calls(work, acc, None if source == "species" else source)
        if calls:
            out[rank] = (calls, rec)
    return out


def holdout_rows(tallies, records, species):
    """holdout.tsv's rows: per genome and level, its held-out calls against OrthoDB's gene models (tallies:
    {(accession, species taxid, level): Counter of verdicts, same_gene, span_90, loci}; records: {(accession,
    level): clade_levels' record})."""
    out = []
    for (acc, sp, level), t in sorted(tallies.items(), key=lambda kv: (kv[0][0], HOLDOUT_LEVELS.index(kv[0][2]))):
        rec = records[(acc, level)]
        info = species.get(str(sp), {})
        nearest = rec.get("nearest")
        out.append(dict(accession=acc, species_taxid=sp, name=info.get("name", "-"),
                        supergroup=info.get("supergroup", "-"), level=level, left_out=rec.get("left_out", "-"),
                        nearest="-" if nearest is None else DEPTH_NAMES[int(nearest)], loci=t["loci"],
                        called=t["loci"] - t["missed"], concordant=t["concordant"], similar=t["similar"],
                        discordant=t["discordant"], missed=t["missed"], same_gene=t["same_gene"],
                        span_90=t["span_90"]))
    return out


def reach_summary(rows):
    """summary.json's orthodb.reach (since 2026-11-10; Figure 1b): holdout.tsv's rows pooled by the rank the genome
    shares with the nearest dataset species kept (a genome's levels of one rank once: the first, the fewest left out)
    -> {rank: genomes, loci, and the shares of the loci called, on the gene of the own-protein call, and OrthoDB's
    gene model}, nearest rank first."""
    pooled, seen = collections.defaultdict(collections.Counter), set()
    for r in rows:
        if r["nearest"] in ("-", None) or (r["accession"], r["nearest"]) in seen:
            continue
        seen.add((r["accession"], r["nearest"]))
        c = pooled[r["nearest"]]
        c["genomes"] += 1
        for k in ("loci", "called", "same_gene", "concordant"):
            c[k] += int(r[k])
    return {rank: dict(genomes=c["genomes"], loci=c["loci"], called=round(c["called"] / c["loci"], 4),
                       same_gene=round(c["same_gene"] / c["loci"], 4), concordant=round(c["concordant"] / c["loci"], 4))
            for rank, c in ((rank, pooled.get(rank)) for rank in reversed(DEPTH_NAMES)) if c and c["loci"]}


def holdout_summary(rows):
    """summary.json's orthodb.clade_hold_out: per level, its genomes and loci, the share called, the share of those
    OrthoDB's gene model (concordant), the share on the package's call's gene, and the nearest dataset species kept
    (the genomes per shared rank)."""
    out = {}
    for level in HOLDOUT_LEVELS:
        rs = [r for r in rows if r["level"] == level]
        if not rs:
            continue
        loci, called = sum(r["loci"] for r in rs), sum(r["called"] for r in rs)
        out[level] = dict(genomes=len(rs), loci=loci, called=round(called / loci, 4) if loci else None,
                          concordant_of_called=round(sum(r["concordant"] for r in rs) / called, 4) if called else None,
                          same_gene_of_called=round(sum(r["same_gene"] for r in rs) / called, 4) if called else None,
                          nearest=dict(collections.Counter(r["nearest"] for r in rs)))
    return out


def verdict_summary(rows, prefix=""):
    """summary.json's figures of orthodb.tsv's rows, of the build's calls or (prefix "held_out_") the held-out ones."""
    verdicts = collections.Counter(r[prefix + "verdict"] for r in rows)
    called = len(rows) - verdicts["missed"]
    return {"species": len({r["accession"] for r in rows}), "pairs": len(rows),
            "called": round(called / len(rows), 4) if rows else None,
            "concordant_of_called": round(verdicts["concordant"] / called, 4) if called else None,
            "verdicts": dict(verdicts),
            "median_recall": median_of([float(r[prefix + "recall"]) for r in rows if r[prefix + "recall"] != "-"])}


def place_summary(rows):
    """summary.json's figures of the held-out calls' places against the package's calls (orthodb.tsv's rows): of the
    loci both called, the share on the same gene, and those spans' median overlap and share overlapping 90% or more;
    the held-out calls the rescue pass made."""
    placed = [r for r in rows if r.get("held_out_place", "-") != "-"]
    spans = [float(r["held_out_span"]) for r in placed if r["held_out_place"] == "same"]
    return {"both_called": len(placed),
            "same_gene": round(len(spans) / len(placed), 4) if placed else None,
            "median_span": median_of(spans),
            "span_at_least_90": round(sum(s >= 0.9 for s in spans) / len(spans), 4) if spans else None,
            "rescued": sum(r.get("held_out_rescued") == "t" for r in rows)}


def resolution_rows(by_marker):
    """resolution.tsv's rows: by_marker {marker: [(identity or None, depth or None)]} (None: no copy of another
    species 70% identical) -> per marker and shared rank (and `all`, and `unmatched`)."""
    out = []
    for marker, found in by_marker.items():
        groups = collections.defaultdict(list)
        unmatched = 0
        for ident, depth in found:
            if ident is None:
                unmatched += 1
            else:
                groups[DEPTHS[depth]].append(ident)
        every = [i for vs in groups.values() for i in vs]
        for rank, values in [(r, groups.get(r, [])) for r in DEPTHS] + [("all", every)]:
            n_extra = unmatched if rank == "all" else 0
            if not values and not n_extra:
                continue
            values = sorted(values)
            row = dict(marker=marker, shared_rank=rank, copies=len(values) + n_extra,
                       under_97=sum(v < 97.0 for v in values) + n_extra,
                       under_99=sum(v < 99.0 for v in values) + n_extra)
            row.update({q: f"{quantile(values, f):.2f}" if values else "-" for q, f in QUANTILES})
            out.append(row)
        if unmatched:
            out.append(dict(marker=marker, shared_rank="unmatched", copies=unmatched, under_97=unmatched,
                            under_99=unmatched, **{q: "-" for q, _f in QUANTILES}))
    return out


def best_ssu(records):
    """Each species' best rRNA record for its 18S (records: records.tsv's rows) -> {record: (species, 18S start, end)}:
    the most confident, then a genome's, then the longest 18S of SSU_MIN bases or more. A record labelled by its genus
    alone (EUKARYOME's) is no species' (until 2026-11-06 the genus was a "species", the nearest of its own)."""
    best = {}
    for r in records:
        if r.get("label_from") == "eukaryome_genus":
            continue
        if r["kind"] == "18S":
            span = (0, int(r["length"]))
        else:
            span = next(((a, b) for n, a, b in rrna.parse_regions(r["regions"]) or () if n == "18S"), None)
        if not span or span[1] - span[0] < SSU_MIN:
            continue
        key = (TIERS.index(r["confidence"]) if r.get("confidence") in TIERS else len(TIERS),
               r["source"] not in ("genome", "curated_unit"), -(span[1] - span[0]), r["record"])
        sp = r["species_taxid"]
        if sp not in best or key < best[sp][0]:
            best[sp] = (key, r["record"], span)
    return {rec: (int(sp), a, b) for sp, (_k, rec, (a, b)) in best.items()}


def ssu_nearest(folder, ncbi, vsearch, threads):
    """[(identity, depth)] of each species' best 18S against the other species' (qc.check's search)."""
    chosen = best_ssu(rows_of(os.path.join(folder, "rrna", "records.tsv")))
    if len(chosen) < 2:
        return []
    seqs = {}
    for name, seq in odb_mod.read_fasta(os.path.join(folder, "rrna", "records.fna.gz")):
        if name in chosen:
            _sp, a, b = chosen[name]
            seqs[(name, "18S")] = seq[a:b]
    checks = qc.check(seqs, {name: chosen[name][0] for name, _l in seqs}, ncbi, vsearch, threads)
    return [(c["nearest_id"], c["nearest_depth"]) for c in checks.values()]


def rrna_rows(folder):
    """rrna.tsv's rows: per source and kind of rRNA record its regions, whole 28S, introns, confidence."""
    groups = collections.defaultdict(list)
    for r in rows_of(os.path.join(folder, "rrna", "records.tsv")):
        groups[(r["source"], r["kind"])].append(r)
    out = []
    for (source, kind), rs in sorted(groups.items()):
        names = [{n for n, _a, _b in rrna.parse_regions(r["regions"]) or ()} for r in rs]
        if kind == "18S":
            names = [n | {"18S"} for n in names]
        lengths = sorted(int(r["length"]) for r in rs if r["length"].isdigit())
        tiers = collections.Counter(r.get("confidence", "-") for r in rs)
        out.append(dict(source=source, kind=kind, records=len(rs),
                        median_length=lengths[len(lengths) // 2] if lengths else 0,
                        **{reg: sum(reg in n for n in names) for reg in REGIONS},
                        lsu_complete=sum(r["lsu_complete"] == "t" for r in rs),
                        introns=sum(r.get("introns", "-") not in ("-", "") for r in rs),
                        **{t: tiers[t] for t in TIERS}))
    return out


MARKER_COLUMNS = ("marker", "accession", "species_taxid", "status", "reference", "recall", "precision", "verdict")


def marker_rows(folder, ncbi, marker_folder):
    """markers.tsv's rows (since 2026-11-07): each marker call (markers/copies.tsv, its protein) against UniProt's
    protein of its species in the marker set (its query set: a species' reference proteome's member, markers.py) by
    shared 5-mers both ways, as orthodb.tsv against OrthoDB's: `concordant` (0.9 both ways: that gene, its model),
    `similar` (0.5: that gene, another model or allele), `discordant` (another gene: a paralog the decoys did not
    tell), `missed` (no call where UniProt has one). Circular in part: the species' own protein was a query."""
    path = os.path.join(folder, "markers", "copies.tsv")
    if not os.path.isfile(path) or not marker_folder:
        return []
    from .odb import read_fasta  # noqa: PLC0415
    refs = collections.defaultdict(dict)  # (marker, species) -> {UniProt accession: protein}
    for name, seq in read_fasta(os.path.join(marker_folder, "refseq_db.faa.gz")):
        marker, taxid, acc = name.split("_", 2)
        refs[(marker, ncbi.species_of(int(taxid)))][acc] = seq
    genomes = {r["accession"]: r for r in rows_of(os.path.join(folder, "genomes.tsv")) if r["representative"] == "t"}
    calls = {(r["locus"], r["accession"]): r for r in rows_of(path)}
    loci_ = [r["locus"] for r in rows_of(os.path.join(folder, "markers", "loci.tsv"))]
    out = []
    for marker in loci_:
        from .package import read_copies  # noqa: PLC0415
        proteins = read_copies(folder, marker, "proteins", "markers")
        for acc, g in sorted(genomes.items()):
            sp = int(g["species_taxid"])
            mine = refs.get((marker, sp))
            if not mine:
                continue
            call = calls.get((marker, acc))
            if call is None or acc not in proteins:
                out.append(dict(marker=marker, accession=acc, species_taxid=sp, status="Missing", reference="-",
                                recall="-", precision="-", verdict="missed"))
                continue
            best, recall, precision = best_match(proteins[acc], mine)
            out.append(dict(marker=marker, accession=acc, species_taxid=sp, status=call["status"], reference=best,
                            recall=f"{recall:.3f}", precision=f"{precision:.3f}", verdict=verdict(recall, precision)))
    return out


PRIMER_COLUMNS = ("pair", "supergroup", "species", "amplified", "share")


def primer_rows(folder, species):
    """primers.tsv's rows (since 2026-11-07): per primer pair and supergroup (and `all`), the species with a record the
    pair runs on, those with one it amplifies (rrna/amplicons.tsv, organelle/amplicons.tsv)."""
    species_of = {}
    for sub, table in (("rrna", "records.tsv"), ("organelle", "genes.tsv")):
        path = os.path.join(folder, sub, table)
        species_of.update({r["record"]: r["species_taxid"] for r in rows_of(path)} if os.path.isfile(path) else {})
    tried, amplified = collections.defaultdict(set), collections.defaultdict(set)
    for sub in ("rrna", "organelle"):
        path = os.path.join(folder, sub, "amplicons.tsv")
        for r in rows_of(path) if os.path.isfile(path) else ():
            sp = species_of.get(r["record"])
            if sp is None:
                continue
            for group in (species.get(sp, {}).get("supergroup", "-"), "all"):
                tried[(r["pair"], group)].add(sp)
                if r["amplifies"] == "t":
                    amplified[(r["pair"], group)].add(sp)
    return [dict(pair=p, supergroup=g, species=len(s), amplified=len(amplified[(p, g)]),
                 share=f"{len(amplified[(p, g)]) / len(s):.3f}")
            for (p, g), s in sorted(tried.items(), key=lambda kv: (kv[0][0], kv[0][1] == "all", kv[0][1]))]


def label_rows(species):
    """labels.tsv's rows: per supergroup (and `all`) and reference its species by how they were named; a reference
    that named none (not in the package: an older format, MIDORI2 not given) left out."""
    counts = collections.defaultdict(collections.Counter)
    held = [ref for ref in REFERENCES if ref in REFERENCES[:3] or
            any(r.get(f"{ref.lower()}_match", "-") != "-" for r in species.values())]
    for r in species.values():
        for ref in held:
            how = r.get(f"{ref.lower()}_match", "-")
            for group in (r.get("supergroup", "-"), "all"):
                counts[(group, ref)][how if how in MATCHES else "-"] += 1
    return [dict(supergroup=g, reference=ref, **{"none" if m == "-" else m: c[m] for m in MATCHES})
            for (g, ref), c in sorted(counts.items(), key=lambda kv: (kv[0][0] == "all", kv[0]))]


def runtime_rows(work, genomes, species):
    """runtime.tsv's rows: each species' genome's run (its done.json), those timed (run since 2026-11-06). The
    done.json files read READ_THREADS at once (genomes.prefetch; since 2026-11-11: one after another)."""
    out = []
    with genomes_mod.cached_summaries():
        genomes_mod.prefetch(work, genomes)
        summaries = {acc: genomes_mod.done_summary(work, acc) for acc in genomes}
    for acc, g in sorted(genomes.items()):
        done = summaries[acc]
        if done is None or done.get("seconds") is None:
            continue
        out.append(dict(accession=acc, species_taxid=g["species_taxid"],
                        supergroup=species.get(g["species_taxid"], {}).get("supergroup", "-"),
                        genome_size=g.get("genome_size"), genome_bytes=done.get("bytes"),
                        searched=len(done.get("searched") or ()) or done.get("loci"), threads=done.get("threads"),
                        seconds=done["seconds"]))
    return out


# --- all of them -------------------------------------------------------------------------------------------------

def median_of(values):
    return round(statistics.median(values), 4) if values else None


# --- the paper's tables (since 2026-11-06; docs/benchmarks.md, "The paper's figures") ------------------------------

SELECTION_COLUMNS = ("supergroup", "fate", "taxa")
BREADTH_COLUMNS = ("supergroup", "ncbi_species", "eukore", "eukore_rrna_only", "compared", "both", "eukore_only",
                   "compared_only")
GENE_COLUMNS = ("accession", "species_taxid", "supergroup", "odb_nearest", "odb_self", "exon_bases", "intron_bases")
ANNOTATION_COLUMNS = ("accession", "species_taxid", "supergroup", "odb_nearest", "odb_self", "locus", "status",
                      "transcript", "gene", "exons", "ann_exons", "exact_exons", "introns", "ann_introns",
                      "exact_introns", "bases", "ann_bases", "shared_bases",
                      "source")  # since 2026-11-10 (the last column): the annotation's (ncbi, ensembl, toga2, ..)
COMPOSITION_COLUMNS = ("supergroup", "species", "with_genome", "rrna_only", "complete_median", "complete_90",
                       "with_rrna", "with_operon", "with_organelle")
FLAG_COLUMNS = ("supergroup", "confidence", "flags", "copies")
REGION_COLUMNS = ("supergroup", "source", "region", "records") + tuple(q for q, _f in QUANTILES)
ORGANELLE_COLUMNS = ("supergroup", "species", "organelle", "gene", "with_gene")
OVERLAP_COLUMNS = ("supergroup", "database", "species", "both", "eukore_only", "database_only")
COPY_COLUMNS = ("accession", "species_taxid", "supergroup", "level", "contigs", "rrna_18s", "rdna_units",
                "rdna_complete")
ODB_NAME = "BUSCO (OrthoDB)"
NEAREST_RANKS = ("genus", "family", "order", "class", "phylum")


def supergroup_of(ncbi, taxid, cache):
    if taxid not in cache:
        from . import labels  # noqa: PLC0415
        cache[taxid] = labels.group(ncbi, taxid)
    return cache[taxid]


def selection_rows(fates, ncbi):
    """selection.tsv: per supergroup what became of NCBI's taxa with an assembly (genomes.fates)."""
    counts, cache = collections.Counter(), {}
    for taxid, fate in fates.items():
        counts[(supergroup_of(ncbi, taxid, cache), fate)] += 1
    return [dict(supergroup=g, fate=f, taxa=n) for (g, f), n in sorted(counts.items())]


def compared_species(table, ncbi, inside=None):
    """The species another reference holds (its record<TAB>taxid table), within inside's clades."""
    out = set()
    with util.open_text(table) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                sp = ncbi.species_of(int(f[1]))
                if ncbi.rank.get(sp) == "species" and (inside is None or inside(sp)):
                    out.add(sp)
    return out


def breadth_rows(species, fates, other, ncbi):
    """breadth.tsv: per supergroup (and `all`) NCBI's species with an assembly, Eukore's species (with a genome; with
    rRNA alone), another reference's (other: its species, or None) and the two's overlap."""
    cache = {}
    ncbi_species = collections.Counter(supergroup_of(ncbi, t, cache) for t, f in (fates or {}).items()
                                       if f not in ("not a species", "placeholder", "hybrid"))
    with_genome = {int(t): r["supergroup"] for t, r in species.items() if r["representative"] != "-"}
    rrna_only = collections.Counter(r["supergroup"] for r in species.values() if r["representative"] == "-")
    theirs = {sp: supergroup_of(ncbi, sp, cache) for sp in (other or ())}
    groups = sorted(set(ncbi_species) | set(with_genome.values()) | set(rrna_only) | set(theirs.values()))
    rows = []
    for g in groups + ["all"]:
        mine = {sp for sp, x in with_genome.items() if g in (x, "all")}
        their = {sp for sp, x in theirs.items() if g in (x, "all")}
        rows.append(dict(supergroup=g, ncbi_species=sum(ncbi_species.values()) if g == "all" else ncbi_species[g],
                         eukore=len(mine), eukore_rrna_only=sum(rrna_only.values()) if g == "all" else rrna_only[g],
                         compared=len(their) if other is not None else "-",
                         both=len(mine & their) if other is not None else "-",
                         eukore_only=len(mine - their) if other is not None else "-",
                         compared_only=len(their - mine) if other is not None else "-"))
    return rows


def reference_species(references, ncbi, inside=None, index=None):
    """{reference name: its NCBI species} of labels.Reference's (PR2, EUKARYOME's long reads, KSGP): each species name
    of the reference that is an NCBI species' name or synonym (rrna.name_index), within inside's clades."""
    index = index if index is not None else rrna.name_index(ncbi)
    out = {}
    for ref in references:
        found = set()
        for name in ref.species:
            t = index.get(name)
            sp = ncbi.species_of(t) if t else None
            if sp and ncbi.rank.get(sp) == "species" and (inside is None or inside(sp)):
                found.add(sp)
        out[ref.name] = found
    return out


def overlap_rows(species, databases, ncbi):
    """overlap.tsv: per supergroup (and `all`) and database ({name: its NCBI species}) the database's species, those
    in Eukore too, Eukore's alone and the database's alone. Eukore's: every species of species.tsv (genome or rRNA)."""
    cache = {}
    mine = {int(t): r.get("supergroup", "-") for t, r in species.items()}
    rows = []
    for name, theirs in databases.items():
        their = {sp: mine.get(sp) or supergroup_of(ncbi, sp, cache) for sp in theirs}
        for g in sorted(set(mine.values()) | set(their.values())) + ["all"]:
            a = {sp for sp, x in mine.items() if g in (x, "all")}
            b = {sp for sp, x in their.items() if g in (x, "all")}
            rows.append(dict(supergroup=g, database=name, species=len(b), both=len(a & b), eukore_only=len(a - b),
                             database_only=len(b - a)))
    return rows


def copy_rows(genomes, species):
    """rrna_copies.tsv: each species' genome's rRNA as extracted (genomes.tsv: its 18S spans, rDNA units, whole
    units). An assembly collapses a tandem rDNA array: these are the copies assembled, not the genome's."""
    return [dict(accession=a, species_taxid=g["species_taxid"],
                 supergroup=species.get(g["species_taxid"], {}).get("supergroup", "-"), level=g.get("level", "-"),
                 contigs=g.get("contigs", "-"), rrna_18s=g.get("rrna_18s", "-"), rdna_units=g.get("rdna_units", "-"),
                 rdna_complete=g.get("rdna_complete", "-"))
            for a, g in sorted(genomes.items()) if g.get("source") != "no_assembly"]


def odb_distance(odb, ncbi, held):
    """{species: (rank shared with the nearest other species of BUSCO's dataset: genus .. phylum, or none; whether
    the species is itself in the dataset)} of the package's species (held)."""
    theirs = {ncbi.species_of(t) for t in odb.species()}
    names = {r: collections.Counter() for r in NEAREST_RANKS}
    for sp in theirs:
        for r, name in ncbi.ranks(sp).items():
            if r in names and name:
                names[r][name] += 1
    out = {}
    for sp in held:
        mine, own = ncbi.ranks(sp), sp in theirs
        out[sp] = (next((r for r in NEAREST_RANKS if mine.get(r) and names[r][mine[r]] - (1 if own else 0) > 0),
                        "none"), "t" if own else "f")
    return out


def annotation_rows(work, genomes, species, distance):
    """annotation.tsv: each compared copy (WORK/done/<accession>/annotation.tsv, annotation.py) with its genome's
    supergroup and distance to BUSCO's dataset. The genomes read READ_THREADS at once (util.threaded; since 2026-11-11:
    one after another)."""
    from . import annotation as ann  # noqa: PLC0415
    from .genomes import done_dir  # noqa: PLC0415
    out = []
    accs = sorted(genomes)
    read = functools.partial(store.members, wanted=(ann.TABLE,))
    for acc, got in zip(accs, util.threaded(read, [done_dir(work, acc) for acc in accs])):
        g = genomes[acc]
        sp = int(g["species_taxid"])
        near, own = distance.get(sp, ("-", "-"))
        for r in store.table_of(store.decoded(ann.TABLE, (got or {}).get(ann.TABLE))):
            out.append(dict(r, accession=acc, species_taxid=sp, odb_nearest=near, odb_self=own,
                            supergroup=species.get(str(sp), {}).get("supergroup", "-")))
    return out


def region_rows(folder, species):
    """rrna_regions.tsv: per supergroup, source and region the region's length's quantiles; the organelle genes'
    (organelle/genes.tsv, those not partial) as source `mitochondrion` or `plastid`, region the gene (SSU: 12S in
    animals' mitochondria; COX1 from its first base to its last, introns and all)."""
    lengths = collections.defaultdict(list)
    for r in rows_of(os.path.join(folder, "rrna", "records.tsv")):
        g = species.get(r["species_taxid"], {}).get("supergroup", "-")
        for name, a, b in rrna.parse_regions(r["regions"]) or ():
            lengths[(g, r["source"], name)].append(b - a)
    for r in rows_of(os.path.join(folder, "organelle", "genes.tsv")):
        if r.get("partial") == "f" and str(r.get("length", "")).isdigit():
            g = species.get(r["species_taxid"], {}).get("supergroup", "-")
            lengths[(g, r["organelle"], r["gene"])].append(int(r["length"]))
    out = []
    for (g, source, name), values in sorted(lengths.items()):
        values.sort()
        out.append(dict(supergroup=g, source=source, region=name, records=len(values),
                        **{q: quantile(values, f) for q, f in QUANTILES}))
    return out


def organelle_rows(folder, species):
    """organelle.tsv: per supergroup its species and those with each organelle gene."""
    held = collections.defaultdict(set)
    for r in rows_of(os.path.join(folder, "organelle", "genes.tsv")):
        held[(r["organelle"], r["gene"])].add(r["species_taxid"])
    per_group = collections.Counter(r["supergroup"] for r in species.values())
    out = []
    for g, n in sorted(per_group.items()):
        for (organelle, gene), taxa in sorted(held.items()):
            out.append(dict(supergroup=g, species=n, organelle=organelle, gene=gene,
                            with_gene=sum(1 for t in taxa if species.get(t, {}).get("supergroup") == g)))
    return out


def annotation_figures(rows):
    """annotation.tsv's rows' headline: genomes, copies, the share on an annotated gene, of those copies' introns the
    share exact, their coding bases' recall and precision."""
    matched = [r for r in rows if r["transcript"] != "-"]
    introns = sum(int(r["introns"]) for r in matched)
    shared = sum(int(r["shared_bases"]) for r in matched)
    return {"genomes": len({r["accession"] for r in rows}), "copies": len(rows),
            "on_an_annotated_gene": round(len(matched) / len(rows), 4) if rows else None,
            "introns_exact": round(sum(int(r["exact_introns"]) for r in matched) / introns, 4) if introns else None,
            "bases_recall": round(shared / max(1, sum(int(r["ann_bases"]) for r in matched)), 4),
            "bases_precision": round(shared / max(1, sum(int(r["bases"]) for r in matched)), 4)}


def count_of(text):
    try:
        return int(text)
    except (TypeError, ValueError):
        return 0


def composition_rows(folder, species, complete):
    """composition.tsv (since 2026-11-10): per supergroup (and `all`) what the package holds of its species
    (species.tsv): with a genome (or a transcriptome, or records: a representative) and with rRNA records alone; their
    genomes' median share of complete BUSCOs and the share at 90% or more (complete: completeness.tsv's rows); the
    species with an rRNA record, with a whole operon (rrna/records.tsv's rDNA) and with an organelle gene."""
    operon = {r["species_taxid"] for r in rows_of(os.path.join(folder, "rrna", "records.tsv")) if r["kind"] == "rDNA"}
    shares = collections.defaultdict(list)
    for r in complete:
        for g in (r["supergroup"], "all"):
            shares[g].append(float(r["complete"]))
    counts = collections.defaultdict(collections.Counter)
    for taxid, r in species.items():
        for g in (r["supergroup"], "all"):
            c = counts[g]
            c["species"] += 1
            c["with_genome"] += r["representative"] != "-"
            c["rrna_only"] += r["representative"] == "-"
            c["with_rrna"] += count_of(r.get("rrna_records")) > 0
            c["with_operon"] += str(taxid) in operon
            c["with_organelle"] += count_of(r.get("organelle_genes")) > 0
    rows = []
    for g in sorted(g for g in counts if g != "all") + ["all"]:
        c, v = counts[g], shares.get(g, [])
        rows.append(dict(supergroup=g, **{k: c[k] for k in COMPOSITION_COLUMNS[1:4]},
                         complete_median=f"{statistics.median(v):.3f}" if v else "-",
                         complete_90=f"{sum(x >= 0.9 for x in v) / len(v):.3f}" if v else "-",
                         **{k: c[k] for k in COMPOSITION_COLUMNS[6:]}))
    return rows


def paper(folder, out, ncbi, genomes, species, odb, work, fates, inside, flags, gene_bases, summary, log,
          databases=None):
    """The tables of the paper's figures beyond the benchmarks' own, of the package alone (each where its inputs are
    there; the other references' columns and rows and the read benchmark: benchmarks/eukore_bench/compare.py, since
    2026-11-12). databases: {name: its NCBI species} of the other references the build read (overlap.tsv, with
    BUSCO's dataset)."""
    if fates:
        write(out, "selection.tsv", SELECTION_COLUMNS, selection_rows(fates, ncbi))
    rows = breadth_rows(species, fates, None, ncbi)
    write(out, "breadth.tsv", BREADTH_COLUMNS, rows)
    summary["breadth"] = next(({k: v for k, v in r.items() if k != "supergroup"} for r in rows
                               if r["supergroup"] == "all"), {})
    databases = dict(databases or {})
    if odb:
        databases[ODB_NAME] = {sp for sp in (ncbi.species_of(t) for t in odb.species())
                               if sp and ncbi.rank.get(sp) == "species" and (inside is None or inside(sp))}
    if databases:
        rows = overlap_rows(species, databases, ncbi)
        write(out, "overlap.tsv", OVERLAP_COLUMNS, rows)
        summary["overlap"] = {r["database"]: {k: r[k] for k in ("species", "both", "eukore_only", "database_only")}
                              for r in rows if r["supergroup"] == "all"}
    write(out, "rrna_copies.tsv", COPY_COLUMNS, copy_rows(genomes, species))
    held = {int(g["species_taxid"]) for g in genomes.values()}
    distance = odb_distance(odb, ncbi, held) if odb else {}
    rows = []
    for acc, g in sorted(genomes.items()):
        sp = int(g["species_taxid"])
        near, own = distance.get(sp, ("-", "-"))
        exon, intron = gene_bases.get(acc, (0, 0))
        rows.append(dict(accession=acc, species_taxid=sp, supergroup=species.get(str(sp), {}).get("supergroup", "-"),
                         odb_nearest=near, odb_self=own, exon_bases=exon, intron_bases=intron))
    write(out, "genes.tsv", GENE_COLUMNS, rows)
    by_group = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        by_group[r["supergroup"]][0] += r["exon_bases"]
        by_group[r["supergroup"]][1] += r["intron_bases"]
    summary["intron_share"] = {g: round(i / (e + i), 4) for g, (e, i) in sorted(by_group.items()) if e + i}
    if work:
        rows = annotation_rows(work, genomes, species, distance)
        if rows:
            write(out, "annotation.tsv", ANNOTATION_COLUMNS, rows)
            summary["annotation"] = dict(annotation_figures(rows),
                                         sources=dict(collections.Counter(r.get("source", "ncbi") for r in rows)))
            # since 2026-11-10: per annotation source, and per rank shared with BUSCO's dataset (Figure 1b)
            for key, of in (("by_source", lambda r: r.get("source") or "ncbi"),
                            ("by_distance", lambda r: "own" if r["odb_self"] == "t" else r["odb_nearest"])):
                groups = collections.defaultdict(list)
                for r in rows:
                    groups[of(r)].append(r)
                summary["annotation"][key] = {g: annotation_figures(rs) for g, rs in sorted(groups.items())}
            log(f"benchmarks: {len(rows)} copies against their genome's annotation "
                f"({ {k: v for k, v in summary['annotation'].items() if not k.startswith('by_')} })")
        else:
            summary["left_out"]["annotation"] = "no genome compared with its annotation (annotation.py)"
    write(out, "flags.tsv", FLAG_COLUMNS, [dict(supergroup=g, confidence=t, flags=f, copies=n)
                                           for (g, t, f), n in sorted(flags.items())])
    write(out, "rrna_regions.tsv", REGION_COLUMNS, region_rows(folder, species))
    write(out, "organelle.tsv", ORGANELLE_COLUMNS, organelle_rows(folder, species))


def run(folder, ncbi, odb_folder=None, vsearch=None, threads=1, work=None, nearest=None, stages=None, out=None,
        log=print, fates=None, inside=None, databases=None, marker_folder=None, clades=None):
    """The benchmarks of the package in `folder` (finished, or the build's unfinished one once its tables are written)
    into out (default folder/benchmarks), with summary.json and (since 2026-11-10) summary.md -> summary.json's
    figures. clades: a trial's number of clades (summary.md says so). ncbi: the taxonomy (the build's, or the package's
    taxonomy/); odb_folder: BUSCO's dataset (orthodb.tsv); vsearch: resolution.tsv; nearest: the build's own
    {(accession, locus): (identity, depth)} of its checks (qc.check), so that the BUSCOs are not searched again;
    work: the build's work folder (runtime.tsv, annotation.tsv); stages: [(step, seconds)] of the package stage
    (stages.tsv). The paper's (paper()): fates, genomes.fates' {taxid: fate} (selection.tsv); inside: a trial's
    clades; databases: {name: NCBI species} of the other references read (overlap.tsv: reference_species, ROD's).
    The other references and the read benchmark: benchmarks/eukore_bench/compare.py (since 2026-11-12)."""
    started = time.time()
    out = out or os.path.join(folder, FOLDER)
    os.makedirs(out, exist_ok=True)
    loci = [r["locus"] for r in rows_of(os.path.join(folder, "busco", "loci.tsv"))]
    # the benchmarks split species by their group (format 5's species.tsv `group`; a format-4 package's supergroup)
    species = {r["taxid"]: dict(r, supergroup=r.get("group") or r.get("supergroup", "-"))
               for r in rows_of(os.path.join(folder, "species.tsv"))}
    genomes = {r["accession"]: r for r in rows_of(os.path.join(folder, "genomes.tsv")) if r["representative"] == "t"}
    species_of = {a: int(g["species_taxid"]) for a, g in genomes.items()}
    summary = {"species": len(genomes), "busco_loci": len(loci), "left_out": {}}

    # the dataset's own gene models of the package's species: which genomes' calls orthodb.tsv compares
    odb, models = None, {}
    if odb_folder and os.path.isfile(os.path.join(odb_folder, "scores_cutoff")):
        held = set(loci)
        odb = odb_mod.Odb(odb_folder, [l for l in odb_mod.all_loci(odb_folder) if l in held])
        models = odb_species(odb, ncbi, set(species_of.values()))
    else:
        summary["left_out"]["orthodb"] = "no BUSCO dataset (--odb)"
    compared = {a for a, sp in species_of.items() if sp in models}
    search = vsearch and nearest is None  # the BUSCOs searched here: a package measured after its build
    if not vsearch:
        summary["left_out"]["resolution"] = "no vsearch"

    statuses, tiers, calls = collections.defaultdict(collections.Counter), collections.defaultdict(
        collections.Counter), {}
    flags = collections.Counter()  # (supergroup, confidence, flags) -> copies: which checks fire, together
    for r in rows_of(os.path.join(folder, "busco", "copies.tsv")):
        acc = r["accession"]
        if acc not in genomes:
            continue
        statuses[acc][r["status"]] += 1
        tiers[acc][r.get("confidence", "-")] += 1
        flags[(species.get(str(species_of[acc]), {}).get("supergroup", "-"), r.get("confidence", "-"),
               r.get("flags", "-"))] += 1
        if acc in compared or search:
            calls[(acc, r["locus"])] = {k: r.get(k, "-") for k in ("status", "query", "contig", "start", "end")}
    # of genomes: a species without an assembly holds EukDetect2's records alone (until 2026-11-06 counted, its
    # "missing" BUSCOs lowering small supergroups' medians)
    rows = completeness({a: g for a, g in genomes.items() if g.get("source") != "no_assembly"}, species, statuses,
                        tiers, len(loci))
    write(out, "completeness.tsv", COMPLETENESS_COLUMNS, rows)
    by_group = collections.defaultdict(list)
    for r in rows:
        by_group[r["supergroup"]].append(float(r["complete"]))
    every = [float(r["complete"]) for r in rows]
    summary["completeness"] = {"median": median_of(every),
                               "at_least_90": round(sum(v >= 0.9 for v in every) / len(every), 4) if every else None,
                               "by_supergroup": {g: median_of(v) for g, v in sorted(by_group.items())}}
    log(f"benchmarks: completeness of {len(rows)} genomes, median {summary['completeness']['median']}")
    composition = composition_rows(folder, species, rows)  # since 2026-11-10: what the package holds per supergroup
    write(out, "composition.tsv", COMPOSITION_COLUMNS, composition)
    summary["composition"] = {k: v for k, v in composition[-1].items() if k != "supergroup"}

    # the copies' coding sequences, where a benchmark needs them: exons from features.tsv, copies per locus
    exons = collections.defaultdict(list)
    gene_bases = collections.defaultdict(lambda: [0, 0])  # accession -> its copies' exon and intron bases
    for f in rows_of(os.path.join(folder, "busco", "features.tsv")):
        if f["kind"] in ("exon", "intron") and f["accession"] in genomes:
            gene_bases[f["accession"]][f["kind"] == "intron"] += int(f["end"]) - int(f["start"])
        wanted = f["accession"] in compared or (search and f["accession"] in genomes)
        if f["kind"] in ("exon", "frameshift") and wanted:  # qc.spliced leaves the frameshifts' bases out
            exons[(f["accession"], f["locus"])].append((f["kind"], int(f["start"]), int(f["end"])))
    odb_rows, busco_nearest = [], collections.defaultdict(list)
    if vsearch and nearest is not None:
        for (acc, locus), found in nearest.items():
            if acc in genomes:
                busco_nearest[locus].append(found)
    need = compared | (set(genomes) if search else set())
    held = {a: held_out_calls(work, a) for a in compared}
    # the clade hold-out (since 2026-11-09): each level's held-out calls against OrthoDB's gene models
    levels = {a: clade_levels(work, a) for a in compared}
    clade_tallies, clade_records = collections.defaultdict(collections.Counter), {}
    for locus in loci:
        cds = {}
        if need and os.path.isfile(os.path.join(folder, "busco", "copies", f"{locus}.fna.gz")):
            cds = {a: qc.spliced(seq, exons[(a, locus)]) for a, seq in package.read_copies(folder, locus).items()
                   if a in need and (a, locus) in calls}
        for acc in sorted(compared):
            references = models[species_of[acc]].get(locus)
            if references:
                odb_rows.append(orthodb_row(acc, species_of[acc], locus, calls.get((acc, locus)), cds.get(acc, ""),
                                            references, ncbi, held[acc].get(locus) if held[acc] else None))
                table = ncbi.code_of(species_of[acc])
                for level, (level_calls, record) in levels[acc].items():
                    t = clade_tallies[(acc, species_of[acc], level)]
                    clade_records[(acc, level)] = record
                    t["loci"] += 1
                    h = level_calls.get(locus)
                    if not h or h["status"] == "Missing":
                        t["missed"] += 1
                        continue
                    t[compare_cds(h["cds"], table, references)[4]] += 1
                    place, span = place_of(calls.get((acc, locus)), h) if calls.get((acc, locus)) else ("-", "-")
                    if place == "same":
                        t["same_gene"] += 1
                        t["span_90"] += float(span) >= 0.9
        if search:
            got = qc.check({(a, locus): s for a, s in cds.items()}, species_of, ncbi, vsearch, threads)
            busco_nearest[locus] = [(c["nearest_id"], c["nearest_depth"]) for c in got.values()]
    if odb:
        for r in odb_rows:
            sp = species.get(str(r["species_taxid"]), {})
            r.update(name=sp.get("name", "-"), supergroup=sp.get("supergroup", "-"))
        write(out, "orthodb.tsv", ORTHODB_COLUMNS, odb_rows)
        summary["orthodb"] = dict(verdict_summary(odb_rows), species=len(compared))
        log(f"benchmarks: {len(odb_rows)} BUSCOs of {len(compared)} species against OrthoDB's own "
            f"({summary['orthodb']['verdicts']})")
        rows = [r for r in odb_rows if r["held_out_verdict"] != "-"]
        if rows:  # the calls made without the species' own proteins: the measure that is not circular
            summary["orthodb"]["held_out"] = dict(verdict_summary(rows, "held_out_"), place=place_summary(rows))
            log(f"benchmarks: held-out calls of {summary['orthodb']['held_out']['species']} species "
                f"({summary['orthodb']['held_out']['verdicts']}; {summary['orthodb']['held_out']['place']})")
        elif compared:
            summary["left_out"]["orthodb_held_out"] = "no held-out calls (genomes run before 2026-11-06, or no --build_dir)"
        rows = holdout_rows(clade_tallies, clade_records, species)
        if any(r["level"] != "species" for r in rows):  # the clade hold-out (since 2026-11-09)
            write(out, "holdout.tsv", HOLDOUT_COLUMNS, rows)
            summary["orthodb"]["clade_hold_out"] = holdout_summary(rows)
            summary["orthodb"]["reach"] = reach_summary(rows)
            log(f"benchmarks: clade hold-out of {len({r['accession'] for r in rows})} genomes: "
                + ", ".join(f"{lv} {s['concordant_of_called']} of {s['called']} called"
                            for lv, s in summary["orthodb"]["clade_hold_out"].items()))
        elif compared:
            summary["left_out"]["clade_hold_out"] = "no clade hold-out calls (genomes run before 2026-11-09, or " \
                                                    "--held_out_ranks none)"

    if vsearch:
        markers = dict(busco_nearest)
        markers[ALL_BUSCOS] = [v for locus, vs in busco_nearest.items() if locus in set(loci) for v in vs]  # BUSCOs alone
        ssu = ssu_nearest(folder, ncbi, vsearch, threads)
        if ssu:
            markers["18S"] = ssu
        rows = resolution_rows({m: v for m, v in markers.items() if v})
        write(out, "resolution.tsv", RESOLUTION_COLUMNS, rows)
        summary["resolution"] = {}
        for r in rows:  # every copy, and (since 2026-11-10, `congeners`) those whose nearest other species is a congener
            if r["marker"] in (ALL_BUSCOS, "18S") and r["shared_rank"] in ("all", "genus"):
                figures = {"copies": r["copies"], "under_99": round(r["under_99"] / r["copies"], 4),
                           "under_97": round(r["under_97"] / r["copies"], 4)}
                entry = summary["resolution"].setdefault(r["marker"], {})
                if r["shared_rank"] == "all":
                    entry.update(figures)
                else:
                    entry["congeners"] = figures
        log(f"benchmarks: resolution of {len(busco_nearest)} BUSCOs and {len(ssu)} species' 18S")

    rows = rrna_rows(folder)
    if rows:
        write(out, "rrna.tsv", RRNA_COLUMNS, rows)
        operons = [r for r in rows if r["kind"] == "rDNA"]
        n_operons = sum(r["records"] for r in operons)
        summary["rrna"] = {"records": sum(r["records"] for r in rows), "operons": n_operons,
                           "lsu_complete_of_operons": round(sum(r["lsu_complete"] for r in operons) / n_operons, 4)
                           if n_operons else None}
    rows = marker_rows(folder, ncbi, marker_folder)
    if rows:  # Eukore's own markers against UniProt's proteins of the same species (since 2026-11-07)
        write(out, "markers.tsv", MARKER_COLUMNS, rows)
        summary["markers"] = {}
        for m in sorted({r["marker"] for r in rows}):
            v = collections.Counter(r["verdict"] for r in rows if r["marker"] == m)
            summary["markers"][m] = dict(species=sum(v.values()), **{k: v[k] for k in ("concordant", "similar",
                                                                                      "discordant", "missed")})
        log(f"benchmarks: markers against UniProt's same species: {summary['markers']}")
    rows = primer_rows(folder, species)
    if rows:
        write(out, "primers.tsv", PRIMER_COLUMNS, rows)
        summary["primers"] = {r["pair"]: float(r["share"]) for r in rows if r["supergroup"] == "all"}
    else:
        summary["left_out"]["primers"] = "no amplicons (a package before format 5)"
    gap = list(rows_of(os.path.join(folder, "intraspecific", "barcode_gap.tsv")) or ())
    if gap:  # the panel's barcode gap (intraspecific.py): per marker, the species resolved of those measured
        by_marker = collections.defaultdict(list)
        for r in gap:
            if r["resolved"] in ("t", "f"):
                by_marker["BUSCO (all)" if r["marker"].endswith("at2759") else r["marker"]].append(r["resolved"] == "t")
        summary["barcode_gap"] = {m: {"species_markers": len(v), "resolved": round(sum(v) / len(v), 4)}
                                  for m, v in sorted(by_marker.items())}
    else:
        summary["left_out"]["barcode_gap"] = "no intraspecific panel (--intraspecific 0, or a package before format 5)"
    rows = label_rows(species)
    write(out, "labels.tsv", LABEL_COLUMNS, rows)
    summary["labels"] = {r["reference"]: round((r["taxid"] + r["species"] + r["synonym"]) / max(1, sum(r[m] for m in (
        "taxid", "species", "synonym", "genus", "none"))), 4) for r in rows if r["supergroup"] == "all"}

    if work:
        rows = runtime_rows(work, genomes, species)
        if rows:
            write(out, "runtime.tsv", RUNTIME_COLUMNS, rows)
            sized = [r for r in rows if str(r["genome_size"]).isdigit() and int(r["genome_size"]) > 0]
            gb = sum(int(r["genome_size"]) for r in sized) / 1e9
            summary["runtime"] = {"genomes": len(rows),
                                  "threads": "/".join(sorted({str(r["threads"]) for r in rows})),
                                  "median_seconds": median_of([float(r["seconds"]) for r in rows]),
                                  "seconds": round(sum(float(r["seconds"]) for r in rows), 1),
                                  "seconds_per_gb": round(sum(float(r["seconds"]) for r in sized) / gb, 1)
                                  if gb else None}  # all genomes' time over all their bases: small ones' overhead
        else:
            summary["left_out"]["runtime"] = "no genome timed (genomes run before 2026-11-06)"
    else:
        summary["left_out"]["runtime"] = "no build folder (--build_dir)"
    paper(folder, out, ncbi, genomes, species, odb, work, fates, inside, flags,
          {a: tuple(v) for a, v in gene_bases.items()}, summary, log, databases)
    if stages is not None:
        write(out, "stages.tsv", ("step", "seconds"),
              [dict(step=s, seconds=f"{t:.1f}") for s, t in stages] +
              [dict(step="benchmarks", seconds=f"{time.time() - started:.1f}")])
    with open(os.path.join(out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    # since 2026-11-10: the headline figures as Markdown tables, for a reader without a JSON viewer
    name = re.sub(r"^\..+-\d+\.(.+)\.partial$", r"\1", os.path.basename(os.path.normpath(folder)))  # util.part_of's
    with open(os.path.join(out, "summary.md"), "w", newline="\n") as fh:
        fh.write(summary_md(summary, name, clades, composition))
    return summary


def pct(v):
    return "-" if v is None else f"{100 * v:.1f}%"


def summary_md(s, name, clades=None, composition=None):
    """summary.json's figures (s) as Markdown: the headline table (the README's), then (composition.tsv's rows given)
    what the package holds per supergroup."""
    lines = ["| benchmark | figure |", "|---|---|",
             f"| package | {name}{f' (a trial: {clades} clades)' if clades else ''}, {s['species']:,} species' "
             f"genomes, {s['busco_loci']} BUSCOs |"]
    c = s.get("composition")
    if c:
        lines.append(f"| species; with a genome; with rRNA records alone | {c['species']:,}; {c['with_genome']:,}; "
                     f"{c['rrna_only']:,} |")
        lines.append(f"| species with an rRNA record; a whole operon; an organelle gene | {c['with_rrna']:,}; "
                     f"{c['with_operon']:,}; {c['with_organelle']:,} |")
    c = s.get("completeness") or {}
    lines.append(f"| median complete BUSCOs per genome | {pct(c.get('median'))} "
                 f"({pct(c.get('at_least_90'))} of genomes at 90% or more) |")
    o = s.get("orthodb")
    if o:
        lines.append(f"| BUSCOs OrthoDB has a gene model of, called | {pct(o['called'])} of {o['pairs']:,} "
                     f"({o['species']} species) |")
        lines.append(f"| of those called, OrthoDB's gene model (5-mers ≥ 90% both ways) | "
                     f"{pct(o['concordant_of_called'])} |")
        h = o.get("held_out")
        if h:
            lines.append(f"| the same, called without the species' own proteins | {pct(h['called'])} called, "
                         f"{pct(h['concordant_of_called'])} of them OrthoDB's gene model ({h['species']} species) |")
            p = h.get("place")
            if p and p.get("both_called"):
                lines.append(f"| of those also called with its own proteins, on the same gene | {pct(p['same_gene'])} "
                             f"of {p['both_called']:,}, {pct(p['span_at_least_90'])} of them sharing 90% of the "
                             f"span; {p['rescued']} called by the rescue pass |")
        for rank, r in (o.get("reach") or {}).items():
            lines.append(f"| nearest dataset species kept shares the {rank} (clade hold-out, {r['genomes']} genomes): "
                         f"called; on the same gene; OrthoDB's gene model | {pct(r['called'])}; {pct(r['same_gene'])}; "
                         f"{pct(r['concordant'])} of {r['loci']:,} |")
    a = s.get("annotation")
    if a:
        lines.append(f"| calls laid on their genome's annotation ({a['genomes']:,} genomes): introns exact; coding "
                     f"bases' recall, precision | {pct(a['introns_exact'])}; {pct(a['bases_recall'])}, "
                     f"{pct(a['bases_precision'])} |")
        for source, r in (a.get("by_source") or {}).items():
            lines.append(f"| the same, {source}'s annotation ({r['genomes']:,} genomes) | {pct(r['introns_exact'])}; "
                         f"{pct(r['bases_recall'])}, {pct(r['bases_precision'])} |")
    for marker, r in (s.get("resolution") or {}).items():
        label = "BUSCO" if marker.startswith("BUSCO") else marker
        if "under_99" in r:
            lines.append(f"| {label} copies under 99% identical to any other species | {pct(r['under_99'])} of "
                         f"{r['copies']:,} |")
        if r.get("congeners"):
            lines.append(f"| {label} copies under 99% identical to their nearest congener | "
                         f"{pct(r['congeners']['under_99'])} of {r['congeners']['copies']:,} |")
    for marker, r in (s.get("barcode_gap") or {}).items():
        lines.append(f"| panel species told from their congeners by {marker} (barcode gap) | {pct(r['resolved'])} of "
                     f"{r['species_markers']:,} |")
    m = s.get("markers")
    if m:
        lines.append("| Eukore's own loci against UniProt's same species: concordant | " + ", ".join(
            f"{k} {v['concordant']}/{v['species']}" for k, v in m.items()) + " |")
    r = s.get("rrna")
    if r:
        lines.append(f"| rRNA records; operons with a whole 28S | {r['records']:,}; "
                     f"{pct(r['lsu_complete_of_operons'])} of {r['operons']:,} |")
    p = s.get("primers")
    if p:
        lines.append("| primer pairs: species amplified of those with a record it runs on | " + ", ".join(
            f"{k} {pct(v)}" for k, v in sorted(p.items(), key=lambda kv: -kv[1])) + " |")
    for ref, modes in (s.get("detection") or {}).items():
        cells = []
        for mode, rank in (("in", "species"), ("held_out", "genus"), ("novel", "class")):
            got = sorted((int(k.rsplit("_", 1)[1]), v.get("f1")) for k, v in modes.items()
                         if k.rsplit("_", 1)[0] == mode)
            if got:
                cells.append(f"{rank}{' held out' if mode != 'in' else ''} " + "/".join(
                    "-" if f is None else f"{f:.2f}" for _d, f in got))
        depths = "/".join(f"{d:,}" for d in sorted({int(k.rsplit('_', 1)[1]) for k in modes}))
        lines.append(f"| detection F1, {ref} ({depths} reads) | {'; '.join(cells)} |")
    for ref, r in (s.get("reads") or {}).items():
        if r.get("censored"):
            lines.append(f"| {ref}'s genome reads censored (all secondary hits left out) | " + ", ".join(
                f"{mode} {pct(v)}" for mode, v in r["censored"].items()) + " |")
    for ref, v in (s.get("labels") or {}).items():
        lines.append(f"| species named in {ref} by name or synonym | {pct(v)} |")
    t = s.get("runtime")
    if t:
        per_gb = f", {t['seconds_per_gb'] / 60:.0f} min per Gb of genome" if t.get("seconds_per_gb") else ""
        lines.append(f"| genome run (BUSCOs, rRNA, tiara), median | {t['median_seconds'] / 60:.1f} min on "
                     f"{t.get('threads', '?')} threads{per_gb} |")
    for what, why in (s.get("left_out") or {}).items():
        lines.append(f"| {what} | not measured: {why} |")
    if composition:
        lines += ["", "| supergroup | species | with a genome | rRNA alone | complete BUSCOs, median | genomes ≥ 90% "
                      "complete | with rRNA | with an operon | with an organelle gene |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for r in composition:
            lines.append("| " + " | ".join(
                [str(r["supergroup"])] + [f"{int(r[k]):,}" for k in ("species", "with_genome", "rrna_only")] +
                [pct(float(r[k])) if r[k] != "-" else "-" for k in ("complete_median", "complete_90")] +
                [f"{int(r[k]):,}" for k in ("with_rrna", "with_operon", "with_organelle")]) + " |")
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("package", help="a package folder (Eukore_<version>/)")
    p.add_argument("--build_dir", "--work", dest="work", metavar="WORK",
                   help="the build's folder (eukore_build.py --build_dir; --work before 2026-11-09): runtime.tsv, and "
                        "the defaults of --odb and --taxdump")
    p.add_argument("--sources", metavar="DIR", help="the build's sources (eukore_build.py --sources; default "
                                                    "WORK/sources): the defaults of --odb and --taxdump")
    p.add_argument("--odb", help="BUSCO's dataset the package was called with (default SOURCES/eukaryota_odb12): "
                                 "orthodb.tsv")
    p.add_argument("--taxdump", help="NCBI's taxdump (default SOURCES/taxdump.tar.gz, else the package's "
                                     "taxonomy/): the whole taxonomy finds the dataset's strains' species")
    p.add_argument("--vsearch", default="vsearch", help="vsearch (resolution.tsv)")
    p.add_argument("--out", help="where the tables go (default PACKAGE/benchmarks; the package's manifest does not "
                                 "list them)")
    p.add_argument("-t", "--threads", type=int, default=4)
    opts = p.parse_args(argv)
    from . import ncbi as ncbi_mod  # noqa: PLC0415
    src = opts.sources or (os.path.join(opts.work, "sources") if opts.work else None)
    odb_folder = opts.odb or (odb_mod.folder_of(src) if src else None)
    taxdump = opts.taxdump or (os.path.join(src, "taxdump.tar.gz") if src else None)
    ncbi = ncbi_mod.Ncbi(taxdump if taxdump and os.path.isfile(taxdump) else os.path.join(opts.package, "taxonomy"))
    summary = run(opts.package, ncbi, odb_folder, shutil.which(opts.vsearch), opts.threads, opts.work,
                  out=opts.out, log=lambda text: print(text, file=sys.stderr, flush=True))
    json.dump(summary, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
