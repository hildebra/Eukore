"""Eukore's benchmarks (since 2026-11-06): how complete, how right and how discriminating a package is. The package
stage measures them as it builds (into the package's benchmarks/, before its manifest, which keeps summary.json's
figures as `benchmarks`); any package can be measured again afterwards:

    python3 -m eukore.benchmarks OUT/Eukore_2026.11 --work WORK -t 8     # into OUT/Eukore_2026.11/benchmarks/

docs/benchmarks.md describes them, docs/benchmarks/plot_benchmarks.py draws them. Each is a table of benchmarks/:
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
- summary.json: the headline figures.
A benchmark whose inputs are missing (BUSCO's dataset, vsearch, the work folder) is left out, and summary.json says so.
"""

import argparse
import collections
import functools
import json
import os
import shutil
import statistics
import sys
import time

from . import odb as odb_mod
from . import genomes as genomes_mod
from . import package, qc, rrna, util

FOLDER = "benchmarks"
KMER = 5
CONCORDANT = 0.9    # recall and precision both this high: the call is OrthoDB's gene model
SIMILAR = 0.5       # both this high: the same gene, another model (exon ends, an isoform, another assembly)
DEPTHS = ("none", "phylum", "class", "order", "family", "genus")  # qc.depth's 0 .. 5: the rank shared
QUANTILES = (("q05", 0.05), ("q25", 0.25), ("median", 0.5), ("q75", 0.75), ("q95", 0.95))
SSU_MIN = 1000      # an 18S this long at least, for resolution.tsv (a read's fragment would understate identity)
REGIONS = ("18S", "ITS1", "5.8S", "ITS2", "28S", "5S")
TIERS = ("high", "medium", "low")
MATCHES = ("species", "synonym", "genus", "-")
REFERENCES = ("PR2", "EUKARYOME", "KSGP")
ALL_BUSCOS = "BUSCO (all)"

COMPLETENESS_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "source", "level", "genome_size",
                        "contigs", "loci", "single", "duplicated", "fragmented", "missing", "complete", "high",
                        "medium", "low", "rrna_18s", "rdna_units", "rdna_complete")
ORTHODB_COLUMNS = ("accession", "species_taxid", "name", "supergroup", "locus", "status", "odb_proteins",
                   "odb_length", "protein_length", "recall", "precision", "guided", "verdict")
HELD_COLUMNS = ("held_out_status", "held_out_query", "held_out_recall", "held_out_precision", "held_out_verdict")
ORTHODB_COLUMNS += HELD_COLUMNS
RESOLUTION_COLUMNS = ("marker", "shared_rank", "copies", "under_97", "under_99") + tuple(q for q, _f in QUANTILES)
RRNA_COLUMNS = ("source", "kind", "records", "median_length") + REGIONS + ("lsu_complete", "introns") + TIERS
LABEL_COLUMNS = ("supergroup", "reference") + tuple("none" if m == "-" else m for m in MATCHES)
RUNTIME_COLUMNS = ("accession", "species_taxid", "supergroup", "genome_size", "genome_bytes", "searched", "threads",
                   "seconds")

# NCBI's nuclear genetic codes (ncbieaa, codons in TCAG order); another table reads as the standard one
BASES = "TCAG"
CODES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}


@functools.lru_cache(maxsize=None)
def codons(table=1):
    aa = CODES.get(table, CODES[1])
    return {a + b + c: aa[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
            for k, c in enumerate(BASES)}


def translate(cds, table=1):
    """A coding sequence's protein (X for a codon with another base than ACGT)."""
    code, cds = codons(table), cds.upper()
    return "".join(code.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


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
    low = min(recall, precision)
    verdict = "concordant" if low >= CONCORDANT else "similar" if low >= SIMILAR else "discordant"
    return recall, precision, length, len(protein.rstrip("*")), verdict


def orthodb_row(acc, species_taxid, locus, call, cds, references, ncbi, held=None):
    """orthodb.tsv's row of one BUSCO of a genome whose species the dataset holds: call (status, query) or None
    (Missing); cds: the call's coding sequence; references: the dataset's proteins of the species [(taxid, protein)];
    held: the held-out call (status, query, cds; genomes.HELD_OUT) or None where the genome has none."""
    row = dict(accession=acc, species_taxid=species_taxid, locus=locus, odb_proteins=len(references),
               odb_length=max(len(p) for _t, p in references))
    table = ncbi.code_of(species_taxid)
    if call is None:
        row.update(status="Missing", protein_length="-", recall="-", precision="-", guided="-", verdict="missed")
    else:
        status, query = call
        recall, precision, length, plen, verdict = compare(cds, table, references)
        q_taxid = query.split("_")[1] if query.count("_") >= 2 else ""  # <busco>_<taxid>_<id>; not an ancestral's
        guided = "t" if q_taxid.isdigit() and ncbi.species_of(int(q_taxid)) == species_taxid else "f"
        row.update(status=status, odb_length=length, protein_length=plen, recall=f"{recall:.3f}",
                   precision=f"{precision:.3f}", guided=guided, verdict=verdict)
    if held is None:
        row.update({c: "-" for c in HELD_COLUMNS})
    elif held[0] == "Missing":
        row.update(held_out_status="Missing", held_out_query="-", held_out_recall="-", held_out_precision="-",
                   held_out_verdict="missed")
    else:
        recall, precision, _length, _plen, verdict = compare(held[2], table, references)
        row.update(held_out_status=held[0], held_out_query=held[1], held_out_recall=f"{recall:.3f}",
                   held_out_precision=f"{precision:.3f}", held_out_verdict=verdict)
    return row


def held_out_calls(work, acc):
    """A genome's held-out calls (WORK/done/<accession>/held_out.tsv) -> {locus: (status, query, cds)}, None without."""
    path = os.path.join(genomes_mod.done_dir(work, acc), genomes_mod.HELD_OUT) if work else None
    if not path or not os.path.isfile(path):
        return None
    return {r["locus"]: (r["status"], r["query"], r["cds"]) for r in rows_of(path)}


def verdict_summary(rows, prefix=""):
    """summary.json's figures of orthodb.tsv's rows, of the build's calls or (prefix "held_out_") the held-out ones."""
    verdicts = collections.Counter(r[prefix + "verdict"] for r in rows)
    called = len(rows) - verdicts["missed"]
    return {"species": len({r["accession"] for r in rows}), "pairs": len(rows),
            "called": round(called / len(rows), 4) if rows else None,
            "concordant_of_called": round(verdicts["concordant"] / called, 4) if called else None,
            "verdicts": dict(verdicts),
            "median_recall": median_of([float(r[prefix + "recall"]) for r in rows if r[prefix + "recall"] != "-"])}


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
    the most confident, then a genome's, then the longest 18S of SSU_MIN bases or more."""
    best = {}
    for r in records:
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


def label_rows(species):
    """labels.tsv's rows: per supergroup (and `all`) and reference its species by how they were named."""
    counts = collections.defaultdict(collections.Counter)
    for r in species.values():
        for ref in REFERENCES:
            how = r.get(f"{ref.lower()}_match", "-")
            for group in (r.get("supergroup", "-"), "all"):
                counts[(group, ref)][how if how in MATCHES else "-"] += 1
    return [dict(supergroup=g, reference=ref, **{"none" if m == "-" else m: c[m] for m in MATCHES})
            for (g, ref), c in sorted(counts.items(), key=lambda kv: (kv[0][0] == "all", kv[0]))]


def runtime_rows(work, genomes, species):
    """runtime.tsv's rows: each species' genome's run (its done.json), those timed (run since 2026-11-06)."""
    from .genomes import done_dir  # noqa: PLC0415
    out = []
    for acc, g in sorted(genomes.items()):
        try:
            with open(os.path.join(done_dir(work, acc), "done.json")) as fh:
                done = json.load(fh)
        except (OSError, ValueError):
            continue
        if done.get("seconds") is None:
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


def run(folder, ncbi, odb_folder=None, vsearch=None, threads=1, work=None, nearest=None, stages=None, out=None,
        log=print):
    """The benchmarks of the package in `folder` (finished, or the build's unfinished one once its tables are written)
    into out (default folder/benchmarks) -> summary.json's figures. ncbi: the taxonomy (the build's, or the package's
    taxonomy/); odb_folder: BUSCO's dataset (orthodb.tsv); vsearch: resolution.tsv; nearest: the build's own
    {(accession, locus): (identity, depth)} of its checks (qc.check), so that the BUSCOs are not searched again;
    work: the build's work folder (runtime.tsv); stages: [(step, seconds)] of the package stage (stages.tsv)."""
    started = time.time()
    out = out or os.path.join(folder, FOLDER)
    os.makedirs(out, exist_ok=True)
    loci = [r["locus"] for r in rows_of(os.path.join(folder, "busco", "loci.tsv"))]
    species = {r["taxid"]: r for r in rows_of(os.path.join(folder, "species.tsv"))}
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
    for r in rows_of(os.path.join(folder, "busco", "copies.tsv")):
        acc = r["accession"]
        if acc not in genomes:
            continue
        statuses[acc][r["status"]] += 1
        tiers[acc][r.get("confidence", "-")] += 1
        if acc in compared or search:
            calls[(acc, r["locus"])] = (r["status"], r.get("query", "-"))
    rows = completeness(genomes, species, statuses, tiers, len(loci))
    write(out, "completeness.tsv", COMPLETENESS_COLUMNS, rows)
    by_group = collections.defaultdict(list)
    for r in rows:
        by_group[r["supergroup"]].append(float(r["complete"]))
    every = [float(r["complete"]) for r in rows]
    summary["completeness"] = {"median": median_of(every),
                               "at_least_90": round(sum(v >= 0.9 for v in every) / len(every), 4) if every else None,
                               "by_supergroup": {g: median_of(v) for g, v in sorted(by_group.items())}}
    log(f"benchmarks: completeness of {len(rows)} genomes, median {summary['completeness']['median']}")

    # the copies' coding sequences, where a benchmark needs them: exons from features.tsv, copies per locus
    exons = collections.defaultdict(list)
    if compared or search:
        for f in rows_of(os.path.join(folder, "busco", "features.tsv")):
            if f["kind"] == "exon" and (f["accession"] in compared or (search and f["accession"] in genomes)):
                exons[(f["accession"], f["locus"])].append(("exon", int(f["start"]), int(f["end"])))
    odb_rows, busco_nearest = [], collections.defaultdict(list)
    if vsearch and nearest is not None:
        for (acc, locus), found in nearest.items():
            if acc in genomes:
                busco_nearest[locus].append(found)
    need = compared | (set(genomes) if search else set())
    held = {a: held_out_calls(work, a) for a in compared}
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
            summary["orthodb"]["held_out"] = verdict_summary(rows, "held_out_")
            log(f"benchmarks: held-out calls of {summary['orthodb']['held_out']['species']} species "
                f"({summary['orthodb']['held_out']['verdicts']})")
        elif compared:
            summary["left_out"]["orthodb_held_out"] = "no held-out calls (genomes run before 2026-11-06, or no --work)"

    if vsearch:
        markers = dict(busco_nearest)
        markers[ALL_BUSCOS] = [v for vs in busco_nearest.values() for v in vs]
        ssu = ssu_nearest(folder, ncbi, vsearch, threads)
        if ssu:
            markers["18S"] = ssu
        rows = resolution_rows({m: v for m, v in markers.items() if v})
        write(out, "resolution.tsv", RESOLUTION_COLUMNS, rows)
        summary["resolution"] = {}
        for r in rows:
            if r["marker"] in (ALL_BUSCOS, "18S") and r["shared_rank"] == "all":
                summary["resolution"][r["marker"]] = {"copies": r["copies"],
                                                      "under_99": round(r["under_99"] / r["copies"], 4),
                                                      "under_97": round(r["under_97"] / r["copies"], 4)}
        log(f"benchmarks: resolution of {len(busco_nearest)} BUSCOs and {len(ssu)} species' 18S")

    rows = rrna_rows(folder)
    if rows:
        write(out, "rrna.tsv", RRNA_COLUMNS, rows)
        operons = [r for r in rows if r["kind"] == "rDNA"]
        n_operons = sum(r["records"] for r in operons)
        summary["rrna"] = {"records": sum(r["records"] for r in rows), "operons": n_operons,
                           "lsu_complete_of_operons": round(sum(r["lsu_complete"] for r in operons) / n_operons, 4)
                           if n_operons else None}
    rows = label_rows(species)
    write(out, "labels.tsv", LABEL_COLUMNS, rows)
    summary["labels"] = {r["reference"]: round((r["species"] + r["synonym"]) / max(1, sum(r[m] for m in (
        "species", "synonym", "genus", "none"))), 4) for r in rows if r["supergroup"] == "all"}

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
        summary["left_out"]["runtime"] = "no work folder (--work)"
    if stages is not None:
        write(out, "stages.tsv", ("step", "seconds"),
              [dict(step=s, seconds=f"{t:.1f}") for s, t in stages] +
              [dict(step="benchmarks", seconds=f"{time.time() - started:.1f}")])
    with open(os.path.join(out, "summary.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("package", help="a package folder (Eukore_<version>/)")
    p.add_argument("--work", help="the build's work folder: runtime.tsv, and the defaults of --odb and --taxdump")
    p.add_argument("--odb", help="BUSCO's dataset the package was called with (default WORK/sources/eukaryota_odb12): "
                                 "orthodb.tsv")
    p.add_argument("--taxdump", help="NCBI's taxdump (default WORK/sources/taxdump.tar.gz, else the package's "
                                     "taxonomy/): the whole taxonomy finds the dataset's strains' species")
    p.add_argument("--vsearch", default="vsearch", help="vsearch (resolution.tsv)")
    p.add_argument("--out", help="where the tables go (default PACKAGE/benchmarks; the package's manifest does not "
                                 "list them)")
    p.add_argument("-t", "--threads", type=int, default=4)
    opts = p.parse_args(argv)
    from . import ncbi as ncbi_mod  # noqa: PLC0415
    src = os.path.join(opts.work, "sources") if opts.work else None
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
