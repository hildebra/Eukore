#!/usr/bin/env python3
"""Copy numbers by read depth (since 2026-11-06): how many rDNA units (and mitochondrial and plastid genomes) a genome
holds per haploid nuclear genome, from the depth of its own sequencing reads. An assembly collapses an rDNA array into
a few units (genomes.tsv's assembly counts are a lower bound); its reads do not.

For a genome of a package:
- runs: its sequencing runs at ENA, the first of: a curated run list (data/depth_runs.tsv: the old reference
  assemblies, whose BioSamples hold no reads), the runs of the assembly's BioSample (NCBI Datasets), its BioProjects'
  runs of its taxid or species, the runs of its species (or its own taxid) whose strain, isolate, cultivar or ecotype
  is the assembly's strain. Only WGS of genomic DNA at random (no cDNA, amplicons, capture, Hi-C); Illumina first,
  then PacBio, then Nanopore; the run with most bases (pick_run).
- reads: reads of the run's first file (read 1 of a pair) until the genome is covered `depth` times (genomes.tsv's
  genome_size), at most `max_reads` (stream_reads). A file of THIN_MAX_BASES or fewer is read whole and every nth
  read taken (a file sorted by position, a run submitted as an aligned BAM, holds the genome's start first); a
  larger one gives its first reads (a flow cell's order, as random over the genome as any).
- the reference (reference): the genome's rDNA units of `rrna/records.tsv`, one per 18S type (the best by confidence,
  a whole 28S, length), else the species' published operons; its Single BUSCO copies (`busco/copies/` or `euk30/`,
  confidence not low, no assembly gap); its species' organelle genes (`organelle/genes.fna.gz`).
- mapping: minimap2 (`sr` for short reads, `map-hifi`, `map-pb`, `map-ont`), primary and supplementary alignments
  (a long read over several units of an array maps to the unit once per unit). An alignment counts where it is
  `min_identity` identical and the read's unaligned ends lie beyond the target's ends (a read from a paralog, or a
  repeat that shares a domain with a target, is clipped inside the target: not counted) (accept).
- depth: per base, over a target's interior (`trim` bases from either end left out: a read over a target's end maps
  only where enough of it overlaps). A BUSCO copy's and an organelle gene's mean over its whole interior; an rDNA
  unit's median over its 18S, 5.8S and 28S (not its spacers, which vary between units; a median, not swayed by an
  insertion some units carry: Drosophila's R1 and R2).
- copies: the rDNA units' depths summed (reads of one unit map to one of the types) over the median depth of the
  single-copy BUSCOs (estimate). Mitochondrial and plastid copies alike, the median of the organelle's genes.

Both alleles of a diploid map to its single-copy genes, so a copy number is per haploid genome. GC-rich rDNA is
covered less by PCR-amplified libraries: an estimate is a lower bound there too, but much the nearer one.
"""

import collections
import gzip
import io
import json
import os
import re
import statistics
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request

from eukore import odb, rrna, util

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
DATASETS_REPORT = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}/dataset_report"
RUN_FIELDS = ("run_accession", "sample_accession", "study_accession", "tax_id", "instrument_platform",
              "instrument_model", "library_strategy", "library_source", "library_selection", "library_layout",
              "read_count", "base_count", "strain", "fastq_ftp")
STRATEGIES = {"WGS"}
SOURCES = {"GENOMIC"}
SELECTIONS = {"RANDOM", "RANDOM PCR", "unspecified", "other"}  # not cDNA, PCR amplicons, capture, restriction digests
PLATFORMS = ("ILLUMINA", "PACBIO_SMRT", "OXFORD_NANOPORE")
HIFI_MODELS = re.compile(r"Sequel II|Revio|Vega", re.I)
IDENTITY = {"sr": 0.95, "map-hifi": 0.95, "map-pb": 0.8, "map-ont": 0.8}
TRIM = {"sr": 100, "map-hifi": 500, "map-pb": 500, "map-ont": 500}
OVERHANG = 20         # a read's unaligned end tolerated inside a target (its adapters, a bad last base)
MIN_INTERIOR = 200    # a target's interior this long at least to count
DEPTH = 5.0           # the genome covered this many times by the reads staged
MAX_READS = 20_000_000
THIN_MAX_BASES = 5e9  # a run's file of this many bases or fewer is read whole and thinned (thinning)
UNEVEN = 0.25         # single-copy BUSCOs under a fifth of their mean depth: a larger share flags `uneven`
MIN_LOCI = 10         # single-copy BUSCOs measured: fewer flags `few_loci`
LOW_DEPTH = 3.0       # their median depth: lower flags `low_depth`
SPREAD = 0.5          # their interquartile range over their median: more flags `spread`
AT_RICH = 0.75        # an organelle's sequences this AT-rich flag `mt_at_rich`, `pt_at_rich`: short-read libraries
                      # cover them poorly (yeast's mtDNA, 83% AT: 7-60x in 5 kb windows where its nuclear genes had 8x)
DEPTH_COLUMNS = ("accession", "species_taxid", "runs", "run_source", "platform", "preset", "reads", "read_length",
                 "rdna_units_ref", "rdna_ref_source", "sc_loci", "sc_depth", "sc_spread", "rdna_depth",
                 "rdna_copies_reads", "rdna_18s_copies", "mt_genes", "mt_copies_reads", "pt_genes",
                 "pt_copies_reads", "flags")
GENE_REGIONS = ("18S", "5.8S", "28S")
CONFIDENCE_RANK = {"high": 0, "medium": 1, "low": 2, "-": 1}


# --- package tables -----------------------------------------------------------------------------------------------

def read_dicts(path):
    """A package table (header first) -> its rows as dicts; [] for no file."""
    if not path or not os.path.isfile(path):
        return []
    with util.open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def fasta(path, keep=None):
    """A FASTA (plain or gzipped; odb.read_fasta) -> {name: sequence}, only `keep`'s names where given."""
    if not os.path.isfile(path):
        return {}
    return {name: seq for name, seq in odb.read_fasta(path) if keep is None or name in keep}


def regions_of(text):
    """records.tsv's `regions` -> [(name, start, end)] (rrna.parse_regions', [] for none)."""
    return rrna.parse_regions(text) or []


def rdna_units(records, accession, species_taxid):
    """The rDNA units to map to -> ([record rows], source): the genome's own operons, one per 18S type (the best by
    confidence, a whole 28S, length; a low one only where the type has nothing better), else the species' published
    operons (EUKARYOME's, ROD's, a curated unit) alike; ([], "-") for none."""
    def best_per_type(rows):
        chosen = {}
        for row in rows:
            key = (CONFIDENCE_RANK.get(row.get("confidence", "-"), 1), row.get("lsu_complete") != "t",
                   -int(row["length"]))
            t = row.get("type", "-")
            if t not in chosen or key < chosen[t][0]:
                chosen[t] = (key, row)
        return [row for _, row in sorted(chosen.values(), key=lambda kv: kv[0])]

    own = [r for r in records if r["kind"] == "rDNA" and r["accession"] == accession
           and any(n == "18S" for n, _, _ in regions_of(r["regions"]))]
    if own:
        return best_per_type(own), "genome"
    published = [r for r in records if r["kind"] == "rDNA" and r["species_taxid"] == str(species_taxid)
                 and any(n == "18S" for n, _, _ in regions_of(r["regions"]))]
    if published:
        return best_per_type(published), "published"
    return [], "-"


def single_copies(copies, accession):
    """The genome's BUSCO copies fit to measure depth by -> [loci]: Single, confidence not low (format 4), no
    assembly gap in them."""
    out = []
    for row in copies:
        if row["accession"] != accession or row.get("status") != "Single":
            continue
        if row.get("confidence", "high") == "low" or "gap" in row.get("flags", "-").split(","):
            continue
        out.append(row["locus"])
    return out


def busco_folder(package):
    for name in ("busco", "euk30"):
        if os.path.isdir(os.path.join(package, name)):
            return os.path.join(package, name)
    return None


def reference(package, accession, species_taxid):
    """The targets of a genome -> ({target: sequence}, {target: (kind, [(start, end) measured])}, info): kind `rdna`,
    `busco`, `mitochondrion` or `plastid`; names `<kind>|<record or locus>`."""
    return references(package, [(accession, species_taxid)])[accession]


def references(package, wanted):
    """reference() of several genomes ([(accession, species taxid)]) in one pass over the package's tables and FASTA
    files (since 2026-11-08: a pass for each genome, minutes each: months for every genome) -> {accession: (seqs,
    meta, info)}."""
    wanted = [(a, str(sp)) for a, sp in wanted]
    accessions, species = {a for a, _sp in wanted}, {sp for _a, sp in wanted}
    # the rDNA records of these genomes and species alone (rdna_units reads its own among them)
    records = [r for r in read_dicts(os.path.join(package, "rrna", "records.tsv"))
               if r["kind"] == "rDNA" and (r["accession"] in accessions or r["species_taxid"] in species)]
    by_acc, by_sp = collections.defaultdict(list), collections.defaultdict(list)
    for r in records:
        by_acc[r["accession"]].append(r)
        by_sp[r["species_taxid"]].append(r)
    units = {a: rdna_units(by_acc[a] + by_sp[sp], a, sp) for a, sp in wanted}
    unit_seqs = fasta(os.path.join(package, "rrna", "records.fna.gz"),
                      {r["record"] for rows, _source in units.values() for r in rows})
    loci = {}
    folder = busco_folder(package)
    if folder:
        copies = collections.defaultdict(list)
        for row in read_dicts(os.path.join(folder, "copies.tsv")):
            if row["accession"] in accessions:
                copies[row["accession"]].append(row)
        loci = {a: single_copies(copies[a], a) for a in accessions}
    by_locus = collections.defaultdict(set)
    for a, ls in loci.items():
        for locus in ls:
            by_locus[locus].add(a)
    copy_seqs = {locus: fasta(os.path.join(folder, "copies", f"{locus}.fna.gz"), accs)
                 for locus, accs in by_locus.items()}
    genes = collections.defaultdict(list)
    for r in read_dicts(os.path.join(package, "organelle", "genes.tsv")):
        if r["species_taxid"] in species:
            genes[r["species_taxid"]].append(r)
    gene_seqs = fasta(os.path.join(package, "organelle", "genes.fna.gz"),
                      {r["record"] for rows in genes.values() for r in rows})
    out = {}
    for accession, sp in wanted:
        seqs, meta = {}, {}
        rows, unit_source = units[accession]
        for row in rows:
            if row["record"] in unit_seqs:
                name = f"rdna|{row['record']}"
                seqs[name] = unit_seqs[row["record"]]
                meta[name] = ("rdna", [(a, b) for n, a, b in regions_of(row["regions"]) if n in GENE_REGIONS],
                              [(a, b) for n, a, b in regions_of(row["regions"]) if n == "18S"])
        for locus in loci.get(accession, ()):
            seq = copy_seqs[locus].get(accession)
            if seq:
                name = f"busco|{locus}"
                seqs[name] = seq
                meta[name] = ("busco", [(0, len(seq))], [])
        for row in genes[sp]:
            if row["record"] in gene_seqs and row.get("partial") != "t":
                name = f"{row['organelle']}|{row['record']}"
                seqs[name] = gene_seqs[row["record"]]
                meta[name] = (row["organelle"], [(0, len(seqs[name]))], [])
        out[accession] = (seqs, meta, {"rdna_units_ref": ",".join(r["record"] for r in rows) or "-",
                                       "rdna_ref_source": unit_source})
    return out


# --- finding and staging the reads --------------------------------------------------------------------------------

def get_text(url, attempts=3):
    """A URL's text (ENA's and NCBI's APIs); "" where the server refuses it (HTTP 4xx but 429: no such record). Raises
    OSError where it was not answered (a timeout, 429, 5xx): find_runs' answer is then not kept, a rerun asks again
    (until 2026-11-08 "" too, kept as "no run" for good)."""
    why = None
    for attempt in range(attempts):
        if attempt:
            time.sleep(5 * attempt)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": util.USER_AGENT})
            with urllib.request.urlopen(request, timeout=120) as response:
                return response.read().decode()
        except urllib.error.HTTPError as e:
            e.close()
            if 400 <= e.code < 500 and e.code != 429:  # the server's answer: no such record
                return ""
            why = e
        except OSError as e:
            why = e
    raise OSError(f"{url}: {why}")


def parse_tsv(text):
    lines = [x for x in text.splitlines() if x.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, x.split("\t"))) for x in lines[1:]]


def ena_runs_of(accession, get=get_text):
    """ENA's runs of a sample, study or run accession -> [row]."""
    query = urllib.parse.urlencode({"accession": accession, "result": "read_run", "fields": ",".join(RUN_FIELDS)})
    return parse_tsv(get(f"{ENA_FILEREPORT}?{query}"))


def ena_runs_by_strain(taxids, strain, get=get_text):
    """ENA's runs of these taxids whose strain, isolate, cultivar or ecotype is `strain` -> [row]."""
    value = strain.replace('"', "")
    taxa = " OR ".join(f"tax_eq({t})" for t in taxids)
    names = " OR ".join(f'{f}="{value}"' for f in ("strain", "isolate", "cultivar", "ecotype"))
    query = urllib.parse.urlencode({
        "result": "read_run", "fields": ",".join(RUN_FIELDS), "limit": "500",
        "query": f'({taxa}) AND ({names}) AND library_strategy="WGS" AND library_source="GENOMIC"'})
    return parse_tsv(get(f"{ENA_SEARCH}?{query}"))


def links_of(accession, get=get_text):
    """An assembly's BioSample and BioProjects (NCBI Datasets) -> (biosample or "", [bioprojects])."""
    try:
        report = json.loads(get(DATASETS_REPORT.format(accession)) or "{}")["reports"][0]
    except (KeyError, IndexError, ValueError, TypeError):
        return "", []
    biosample = report.get("assembly_info", {}).get("biosample", {}).get("accession", "")
    projects = sorted(set(re.findall(r'"bioproject_accession":\s*"([A-Z]+\d+)"', json.dumps(report))))
    return biosample, projects


def preset_of(platform, model):
    if platform == "ILLUMINA":
        return "sr"
    if platform == "PACBIO_SMRT":
        return "map-hifi" if HIFI_MODELS.search(model or "") else "map-pb"
    return "map-ont"


def usable(row):
    return (row.get("library_strategy") in STRATEGIES and row.get("library_source") in SOURCES
            and row.get("library_selection") in SELECTIONS and row.get("instrument_platform") in PLATFORMS
            and bool(row.get("fastq_ftp")))


def pick_run(rows):
    """The run to take of a genome's runs: usable, Illumina before PacBio before Nanopore, then most bases -> row or
    None."""
    rows = [r for r in rows if usable(r)]
    if not rows:
        return None

    def key(r):
        bases = int(r["base_count"]) if r.get("base_count", "").isdigit() else 0
        return (PLATFORMS.index(r["instrument_platform"]), -bases, r["run_accession"])
    return min(rows, key=key)


def read_curated(path):
    """data/depth_runs.tsv -> {accession: [runs]}."""
    return {r["accession"]: [x for x in r["runs"].split(",") if x] for r in read_dicts(path)}


def find_runs(genome, curated, get=get_text):
    """A genome's run to take -> {"how": manual|biosample|bioproject|strain|none, "run": row or None, "biosample":
    ...}. `genome` is a genomes.tsv row. A BioProject's runs count only of the genome's taxid or species (a project
    may hold many species)."""
    accession = genome["accession"]
    if accession in curated:
        rows = [r for run in curated[accession] for r in ena_runs_of(run, get)]
        run = pick_run(rows) or (rows[0] if rows else None)  # a curated run is taken as curated
        return {"how": "manual", "run": run, "biosample": "-"}
    biosample, projects = links_of(accession, get)
    found = {"biosample": biosample or "-"}
    if biosample:
        run = pick_run(ena_runs_of(biosample, get))
        if run:
            return dict(found, how="biosample", run=run)
    taxids = sorted({genome["taxid"], genome["species_taxid"]})
    for project in projects:
        run = pick_run([r for r in ena_runs_of(project, get) if r.get("tax_id") in taxids])
        if run:
            return dict(found, how="bioproject", run=run)
    strain = genome.get("strain", "na")
    if strain and strain not in ("na", "-"):
        run = pick_run(ena_runs_by_strain(taxids, strain, get))
        if run:
            return dict(found, how="strain", run=run)
    return dict(found, how="none", run=None)


def first_file(fastq_ftp):
    """A run's file to read: read 1 of a pair (`_1.fastq.gz`), else its one file -> https URL, or None where it has
    none (a curated run without files; until 2026-11-08 an IndexError that stopped the stage)."""
    files = [f for f in (fastq_ftp or "").split(";") if f]
    if not files:
        return None
    one = next((f for f in files if f.endswith("_1.fastq.gz")), files[0])
    return one if one.startswith(("http://", "https://")) else "https://" + one


def thinning(spots, read_length, max_bases, thin_max_bases=THIN_MAX_BASES):
    """How many reads of a run's file go per read taken (1: its first reads). A file of thin_max_bases or fewer (its
    spots times its reads' length) is read whole and thinned, so that the reads taken lie over the whole genome
    even where the file is sorted by position (a run submitted as an aligned BAM)."""
    file_bases = spots * read_length
    if not spots or file_bases > thin_max_bases or file_bases <= max_bases:
        return 1
    return max(1, int(file_bases // max_bases))


def take_reads(stream, out, max_reads, max_bases, spots=0, thin_max_bases=THIN_MAX_BASES):
    """FASTQ records from a binary stream (gzip) into `out` (a binary writer) until max_reads or max_bases: the
    first ones, or every nth over the whole file where its `spots` (ENA's read_count) say it is small enough
    (thinning, by the first 1,000 reads' length) -> (reads, bases, n)."""
    reads = bases = 0
    with gzip.GzipFile(fileobj=stream) as gz:
        text = io.BufferedReader(gz, 1 << 20)

        def records():
            while True:
                record = [text.readline() for _ in range(4)]
                if not record[3]:
                    return
                yield record
        source = records()
        head = [r for _, r in zip(range(1000), source)]
        length = sum(len(r[1].rstrip()) for r in head) / len(head) if head else 0
        every = thinning(spots, length, max_bases, thin_max_bases)
        for i, record in enumerate(r for part in (head, source) for r in part):
            if reads >= max_reads or bases >= max_bases:
                break
            if i % every:
                continue
            out.write(b"".join(record))
            reads += 1
            bases += len(record[1].rstrip())
    return reads, bases, every


def stream_reads(url, dest, max_reads, max_bases, attempts=3, opener=None, spots=0):
    """A run's file's reads (take_reads) into dest (gzip; written whole or not at all) -> (reads, bases, n), or
    raises OSError."""
    opener = opener or (lambda u: urllib.request.urlopen(
        urllib.request.Request(u, headers={"User-Agent": util.USER_AGENT}), timeout=300))
    why = None
    for _ in range(attempts):
        part = util.part_of(dest)
        try:
            with opener(url) as response, gzip.open(part, "wb", compresslevel=6) as fh:
                got = take_reads(response, fh, max_reads, max_bases, spots)
            os.replace(part, dest)
            return got
        except (OSError, EOFError) as e:
            why = e
            if os.path.exists(part):
                os.remove(part)
    raise OSError(f"{url}: {why}")


# --- mapping and depth --------------------------------------------------------------------------------------------

def accept(fields, min_identity, overhang=OVERHANG):
    """A PAF line's fields -> whether the alignment counts: identical enough, and the read's unaligned ends beyond
    the target's ends (or within `overhang` bases)."""
    qlen, qs, qe, strand = int(fields[1]), int(fields[2]), int(fields[3]), fields[4]
    tlen, ts, te = int(fields[6]), int(fields[7]), int(fields[8])
    matches, block = int(fields[9]), int(fields[10])
    if block == 0 or matches / block < min_identity:
        return False
    left, right = (qs, qlen - qe) if strand == "+" else (qlen - qe, qs)
    return (left <= overhang or ts <= overhang) and (right <= overhang or tlen - te <= overhang)


def depths(paf_lines, lengths, min_identity):
    """PAF lines -> ({target: per-base depth list}, reads counted, mean length of the reads counted)."""
    diff = {t: [0] * (n + 1) for t, n in lengths.items()}
    counted, qlens = set(), []
    for line in paf_lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or f[5] not in diff or "tp:A:S" in f[12:]:
            continue
        if not accept(f, min_identity):
            continue
        d = diff[f[5]]
        d[int(f[7])] += 1
        d[int(f[8])] -= 1
        if f[0] not in counted:
            counted.add(f[0])
            qlens.append(int(f[1]))
    out = {}
    for t, d in diff.items():
        run, cov = 0, []
        for x in d[:-1]:
            run += x
            cov.append(run)
        out[t] = cov
    return out, len(counted), (sum(qlens) / len(qlens) if qlens else 0.0)


def interior_depth(cov, spans, trim, median=False):
    """The mean (or median) depth over the spans' bases at least `trim` from the target's ends, or None where under
    MIN_INTERIOR bases are left."""
    lo, hi = trim, len(cov) - trim
    values = [cov[i] for a, b in spans for i in range(max(a, lo), min(b, hi))]
    if len(values) < MIN_INTERIOR:
        return None
    return float(statistics.median(values)) if median else sum(values) / len(values)


def minimap2_paf(minimap2, ref, reads, preset, threads=1):
    """minimap2's PAF lines of the reads on the reference (secondary alignments left out)."""
    command = [minimap2, "-x", preset, "-t", str(threads), "--secondary=no", "-c", ref, reads]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with process.stdout:
        yield from process.stdout
    with process.stderr:
        error = process.stderr.read()
    if process.wait() != 0:
        raise OSError(f"minimap2 failed: {error[-300:]}")


def quartiles(values):
    if len(values) < 4:
        return None
    q = statistics.quantiles(values, n=4)
    return q[0], q[2]


def estimate(target_depth, meta):
    """Per target depths ({target: depth or None}) and their kinds -> the copy numbers and their flags."""
    sc = sorted(v for t, v in target_depth.items() if meta[t][0] == "busco" and v is not None)
    out = {"sc_loci": len(sc), "sc_depth": None, "sc_spread": None, "rdna_depth": None, "rdna_copies_reads": None,
           "rdna_18s_copies": None, "mt_genes": 0, "mt_copies_reads": None, "pt_genes": 0, "pt_copies_reads": None}
    flags = []
    if not sc:
        out["flags"] = "no_single_copy"
        return out
    sc_depth = statistics.median(sc)
    out["sc_depth"] = sc_depth
    q = quartiles(sc)
    if q and sc_depth > 0:
        out["sc_spread"] = (q[1] - q[0]) / sc_depth
        if out["sc_spread"] > SPREAD:
            flags.append("spread")
    mean = sum(sc) / len(sc)
    if sum(1 for v in sc if v < mean / 5) > UNEVEN * len(sc):
        flags.append("uneven")  # reads of a few regions: a run sorted by position, a capture
    if len(sc) < MIN_LOCI:
        flags.append("few_loci")
    if sc_depth < LOW_DEPTH:
        flags.append("low_depth")
    rdna = [t for t in target_depth if meta[t][0] == "rdna"]
    if not rdna:
        flags.append("no_rdna_unit")
    elif sc_depth > 0:
        out["rdna_depth"] = sum(target_depth[t] or 0 for t in rdna)
        out["rdna_copies_reads"] = out["rdna_depth"] / sc_depth
        s18 = [target_depth.get(t + "#18S") for t in rdna]
        if any(v is not None for v in s18):
            out["rdna_18s_copies"] = sum(v or 0 for v in s18) / sc_depth
    for kind, key in (("mitochondrion", "mt"), ("plastid", "pt")):
        genes = [v for t, v in target_depth.items() if meta[t][0] == kind and v is not None]
        out[f"{key}_genes"] = len(genes)
        if genes and sc_depth > 0:
            out[f"{key}_copies_reads"] = statistics.median(genes) / sc_depth
    out["flags"] = ",".join(flags) or "-"
    return out


def measure(package, genome, reads, preset, minimap2="minimap2", threads=1, folder=None, targets=None):
    """A genome's copy numbers from its staged reads -> (estimate dict, [(target, kind, length, depth)]). targets: its
    reference()'s, where read already (references(): several genomes' in one pass)."""
    accession, species_taxid = genome["accession"], genome["species_taxid"]
    seqs, meta, info = targets or reference(package, accession, species_taxid)
    folder = folder or os.path.dirname(reads)
    ref = os.path.join(folder, f"{accession}.ref.fna")
    with open(ref, "w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n{seq}\n")
    cov, counted, mean_len = depths(minimap2_paf(minimap2, ref, reads, preset, threads),
                                    {t: len(s) for t, s in seqs.items()}, IDENTITY[preset])
    os.remove(ref)
    trim = TRIM[preset]
    target_depth, rows = {}, []
    for t, (kind, spans, spans_18s) in meta.items():
        target_depth[t] = interior_depth(cov[t], spans, trim, median=kind == "rdna")
        if kind == "rdna":
            target_depth[t + "#18S"] = interior_depth(cov[t], spans_18s, trim, median=True)
        rows.append((t, kind, len(seqs[t]), gc(seqs[t]), target_depth[t]))
    meta_all = dict(meta)
    meta_all.update({t + "#18S": ("rdna18s", [], []) for t, m in meta.items() if m[0] == "rdna"})
    result = estimate(target_depth, meta_all)
    flags = [] if result["flags"] == "-" else result["flags"].split(",")
    for kind, key in (("mitochondrion", "mt"), ("plastid", "pt")):
        organelle = "".join(s for t, s in seqs.items() if meta[t][0] == kind)
        if organelle and 1 - gc(organelle) >= AT_RICH:
            flags.append(f"{key}_at_rich")
    result["flags"] = ",".join(flags) or "-"
    result.update(info)
    result["mapped_reads"] = counted
    result["mapped_read_length"] = round(mean_len, 1)
    return result, rows


def gc(seq):
    """The share of G and C among a sequence's ACGT bases."""
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / acgt if acgt else 0.0


def fmt(value, digits=1):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)
