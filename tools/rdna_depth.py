#!/usr/bin/env python3
"""rDNA (and mitochondrial and plastid) copy numbers of a package's genomes by read depth (since 2026-11-06; the
method in eukore/depth.py). Three steps, each resumable (a rerun does only what is not done), the first alone online:

    # on a node with internet: each genome's run found at ENA, its first reads staged (WORK/runs, WORK/reads)
    python3 tools/rdna_depth.py stage --package Eukore_<v> --work WORK [--accessions A,B | --species T,U | --all] -t 4
    # on any node (minimap2 on PATH or --minimap2): the reads mapped, the depths measured (WORK/depth)
    python3 tools/rdna_depth.py measure --package Eukore_<v> --work WORK -t 16
    # the table: WORK/depth.tsv; --into_package also writes it into the package as depth.tsv and adds
    # rdna_copies_reads, mt_copies_reads and pt_copies_reads to its genomes.tsv (the manifest's files updated)
    python3 tools/rdna_depth.py table --package Eukore_<v> --work WORK [--into_package]

The genomes are the package's representatives (genomes.tsv's `representative` t), those of --accessions or --species
(species taxids) where given; --all for every one (~29,000: the download is --depth times each genome's size, at most
--max_reads reads, plus whole files of up to 5 Gb where thinned; stage them in batches).
A run is curated in data/depth_runs.tsv (accession, runs, note) where the assembly's BioSample holds no reads and its
strain is named otherwise at ENA (the old reference assemblies).
"""

import argparse
import concurrent.futures
import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from eukore import depth, package as package_mod, util  # noqa: E402

CURATED = os.path.join(HERE, "data", "depth_runs.tsv")


def genomes_of(package, accessions=None, species=None, every=False):
    rows = [r for r in depth.read_dicts(os.path.join(package, "genomes.tsv"))
            if r.get("representative") == "t" and r["accession"].startswith("GC")]
    if accessions:
        rows = [r for r in rows if r["accession"] in accessions]
    elif species:
        rows = [r for r in rows if r["species_taxid"] in species]
    elif not every:
        sys.exit("name the genomes: --accessions, --species or --all")
    return rows


def paths(work, accession):
    return {"runs": os.path.join(work, "runs", f"{accession}.json"),
            "reads": os.path.join(work, "reads", f"{accession}.fq.gz"),
            "staged": os.path.join(work, "reads", f"{accession}.json"),
            "depth": os.path.join(work, "depth", f"{accession}.json"),
            "targets": os.path.join(work, "depth", f"{accession}.targets.tsv")}


def write_json(path, data):
    part = util.part_of(path)
    with open(part, "w", newline="\n") as fh:
        json.dump(data, fh, indent=1)
    os.replace(part, path)


def stage_one(genome, work, curated, target_depth, max_reads):
    """A genome's run found and its first reads staged -> a line of what was done."""
    acc = genome["accession"]
    p = paths(work, acc)
    if os.path.isfile(p["staged"]):
        return f"{acc}: staged"
    if os.path.isfile(p["runs"]):
        with open(p["runs"]) as fh:
            found = json.load(fh)
    else:
        found = depth.find_runs(genome, curated)
        write_json(p["runs"], found)
    run = found["run"]
    if run is None:
        return f"{acc}: no run ({found['how']})"
    size = int(genome["genome_size"]) if genome.get("genome_size", "").isdigit() else 0
    max_bases = int(target_depth * size) if size else 10 ** 12
    url = depth.first_file(run["fastq_ftp"])
    spots = int(run["read_count"]) if run.get("read_count", "").isdigit() else 0
    try:
        reads, bases, every = depth.stream_reads(url, p["reads"], max_reads, max_bases, spots=spots)
    except OSError as e:
        return f"{acc}: {run['run_accession']} not staged: {e}"
    write_json(p["staged"], {"accession": acc, "run": run["run_accession"], "how": found["how"], "url": url,
                             "platform": run["instrument_platform"], "model": run["instrument_model"],
                             "preset": depth.preset_of(run["instrument_platform"], run["instrument_model"]),
                             "reads": reads, "bases": bases, "every": every, "genome_size": size,
                             "staged": datetime.datetime.now().isoformat(timespec="seconds")})
    return (f"{acc}: {run['run_accession']} ({found['how']}), {reads} reads"
            f"{f' (1 in {every})' if every > 1 else ''}, {bases / max(size, 1):.1f}x")


def stage(args):
    genomes = genomes_of(args.package, args.accessions, args.species, args.all)
    for sub in ("runs", "reads"):
        os.makedirs(os.path.join(args.work, sub), exist_ok=True)
    curated = depth.read_curated(args.runs)
    with concurrent.futures.ThreadPoolExecutor(args.threads) as pool:
        for line in pool.map(lambda g: stage_one(g, args.work, curated, args.depth, args.max_reads), genomes):
            print(line, flush=True)


def measure(args):
    genomes = genomes_of(args.package, args.accessions, args.species, True)
    os.makedirs(os.path.join(args.work, "depth"), exist_ok=True)
    for genome in genomes:
        acc = genome["accession"]
        p = paths(args.work, acc)
        if os.path.isfile(p["depth"]) or not os.path.isfile(p["staged"]):
            continue
        with open(p["staged"]) as fh:
            staged = json.load(fh)
        result, rows = depth.measure(args.package, genome, p["reads"], staged["preset"], args.minimap2,
                                     args.threads, os.path.join(args.work, "depth"))
        result.update({k: staged[k] for k in ("run", "how", "platform", "preset", "reads", "bases")})
        part = util.part_of(p["targets"])
        with open(part, "w", newline="\n") as fh:
            fh.write("target\tkind\tlength\tgc\tdepth\n")
            for t, kind, length, gc, d in rows:
                fh.write(f"{t}\t{kind}\t{length}\t{gc:.3f}\t{depth.fmt(d)}\n")
        os.replace(part, p["targets"])
        write_json(p["depth"], result)
        print(f"{acc}: rDNA {depth.fmt(result['rdna_copies_reads'])} copies, single-copy depth "
              f"{depth.fmt(result['sc_depth'])} over {result['sc_loci']} loci, flags {result['flags']}", flush=True)


def row_of(genome, work):
    """A genome's depth.tsv row (depth.DEPTH_COLUMNS)."""
    p = paths(work, genome["accession"])
    row = {c: "-" for c in depth.DEPTH_COLUMNS}
    row.update(accession=genome["accession"], species_taxid=genome["species_taxid"])
    if os.path.isfile(p["depth"]):
        with open(p["depth"]) as fh:
            r = json.load(fh)
        row.update(runs=r["run"], run_source=r["how"], platform=r["platform"], preset=r["preset"],
                   reads=str(r["reads"]), read_length=depth.fmt(r["mapped_read_length"], 0),
                   rdna_units_ref=r["rdna_units_ref"], rdna_ref_source=r["rdna_ref_source"],
                   sc_loci=str(r["sc_loci"]), sc_depth=depth.fmt(r["sc_depth"]),
                   sc_spread=depth.fmt(r["sc_spread"], 2), rdna_depth=depth.fmt(r["rdna_depth"]),
                   rdna_copies_reads=depth.fmt(r["rdna_copies_reads"], 0),
                   rdna_18s_copies=depth.fmt(r["rdna_18s_copies"], 0), mt_genes=str(r["mt_genes"]),
                   mt_copies_reads=depth.fmt(r["mt_copies_reads"], 0), pt_genes=str(r["pt_genes"]),
                   pt_copies_reads=depth.fmt(r["pt_copies_reads"], 0), flags=r["flags"])
    elif os.path.isfile(p["runs"]):
        with open(p["runs"]) as fh:
            found = json.load(fh)
        row.update(run_source=found["how"], flags="no_runs" if found["run"] is None else "not_measured")
    else:
        row["flags"] = "not_staged"
    return row


GENOME_COLUMNS = ("rdna_copies_reads", "mt_copies_reads", "pt_copies_reads")


READS_ELEMENTS = (("rrna", "rDNA_unit", "rdna_copies_reads"), ("rrna", "18S", "rdna_18s_copies"),
                  ("mitochondrion", "genome", "mt_copies_reads"), ("plastid", "genome", "pt_copies_reads"))


def reads_rows(rows):
    """depth.tsv's rows -> copy_numbers.tsv.gz's of basis reads: rDNA units and 18S, mitochondrial and plastid genomes
    per haploid nuclear genome (partial -: not known), where measured."""
    return [dict(accession=r["accession"], species_taxid=r["species_taxid"], set=s, element=e, copies=r[c],
                 partial="-", basis="reads")
            for r in rows for s, e, c in READS_ELEMENTS if r.get(c, "-") not in ("-", "")]


def into_package(package, rows):
    """depth.tsv written into the package, genomes.tsv given GENOME_COLUMNS, and the manifest's files and
    inputs.read_depth updated."""
    by_acc = {r["accession"]: r for r in rows}
    with open(os.path.join(package, "depth.tsv"), "w", newline="\n") as fh:
        fh.write("\t".join(depth.DEPTH_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in depth.DEPTH_COLUMNS) + "\n")
    path = os.path.join(package, "genomes.tsv")
    with open(path) as fh:
        lines = fh.read().splitlines()
    header = lines[0].split("\t")
    keep = [i for i, c in enumerate(header) if c not in GENOME_COLUMNS]
    at = header.index("local") if "local" in header else len(header)
    at = sum(1 for i in keep if i < at)
    out = []
    for n, line in enumerate(lines):
        f = line.split("\t")
        f = [f[i] for i in keep]
        extra = list(GENOME_COLUMNS) if n == 0 else [by_acc.get(f[0], {}).get(c, "-") for c in GENOME_COLUMNS]
        out.append("\t".join(f[:at] + extra + f[at:]))
    part = util.part_of(path)
    with open(part, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    os.replace(part, path)
    written = ["genomes.tsv", "depth.tsv"]
    copies = os.path.join(package, package_mod.COPY_NUMBERS)
    if os.path.isfile(copies):  # since 2026-11-08: its rows of basis reads, replaced
        kept = [r for r in package_mod.read_gz(copies) if r["basis"] != "reads"]
        package_mod.write_gz(copies, package_mod.COPY_NUMBER_COLUMNS, kept + reads_rows(rows))
        written.append(package_mod.COPY_NUMBERS)
    manifest_path = os.path.join(package, "manifest.json")
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    for rel in written:
        full = os.path.join(package, rel)
        manifest["files"][rel] = [os.path.getsize(full), util.md5(full)]
    manifest["inputs"]["read_depth"] = {
        "tool": "tools/rdna_depth.py", "measured": sum(1 for r in rows if r["rdna_copies_reads"] != "-"),
        "date": datetime.date.today().isoformat()}
    write_json(manifest_path, manifest)


def table(args):
    genomes = genomes_of(args.package, args.accessions, args.species, True)
    rows = [row_of(g, args.work) for g in genomes]
    if not args.all and not (args.accessions or args.species):
        rows = [r for r in rows if r["flags"] != "not_staged"]
    out = os.path.join(args.work, "depth.tsv")
    with open(out, "w", newline="\n") as fh:
        fh.write("\t".join(depth.DEPTH_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in depth.DEPTH_COLUMNS) + "\n")
    print(f"{out}: {len(rows)} genomes, {sum(1 for r in rows if r['rdna_copies_reads'] != '-')} measured")
    if args.into_package:
        into_package(args.package, rows)
        print(f"{args.package}: depth.tsv, genomes.tsv's {', '.join(GENOME_COLUMNS)}, "
              f"{package_mod.COPY_NUMBERS}'s rows of basis reads (where the package has it), manifest updated")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("step", choices=("stage", "measure", "table"))
    parser.add_argument("--package", required=True, help="the Eukore package folder")
    parser.add_argument("--work", required=True)
    parser.add_argument("--accessions", type=lambda s: set(s.split(",")))
    parser.add_argument("--species", type=lambda s: set(s.split(",")), help="species taxids")
    parser.add_argument("--all", action="store_true", help="every representative genome")
    parser.add_argument("--runs", default=CURATED, help="curated runs (accession, runs, note)")
    parser.add_argument("--depth", type=float, default=depth.DEPTH, help="genome coverage of the reads staged")
    parser.add_argument("--max_reads", type=int, default=depth.MAX_READS)
    parser.add_argument("--minimap2", default="minimap2")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--into_package", action="store_true")
    args = parser.parse_args(argv)
    {"stage": stage, "measure": measure, "table": table}[args.step](args)


if __name__ == "__main__":
    main()
