#!/usr/bin/env python3
"""Build EUKDB, the standalone eukaryote reference (FORMAT.md; since 2026-11-05).

    python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT -t 16                 # all stages on a node with internet
    python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage fetch -t 8    # a cluster: the node with internet,
    python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage process -t 32 #   then any node (nothing fetched),
    python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage package       #   then the package

Stages:
- fetch: the sources into WORK/sources (EukDetect2 and EukDetect1, NCBI's GenBank and RefSeq assembly tables of the
  eukaryote groups, NCBI's taxdump, EUKARYOME 2.1, ROD v1.2), then the genome set (genomes.jobs) into OUT/genomes;
- process: each genome's euk30 loci placed and its rRNA extracted (WORK/done/<accession>/), EukDetect2's copies' exon
  models (WORK/exons), the published reads' regions (WORK/published); a genome not fetched yet waits for the next run;
- package: OUT/EUKDB_<version>/ written (package.py), then OUT/genomes cut back to --keep_gb (genomes.prune: the
  model species' and the smallest kept first).
`all` runs them in turn. Every stage is resumable: a rerun does only what is not done.
"""

import argparse
import datetime
import functools
import glob
import json
import multiprocessing
import os
import shutil
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eukdb import exons, genomes, ncbi as ncbi_mod, package, rrna, sources, species, util  # noqa: E402

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def say(text):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {text}", flush=True)


# --- sources -----------------------------------------------------------------------------------------------------

def fetch_sources(src_dir, refseq=True):
    """Every source into src_dir (kept when there and whole) -> {name: folder or file}."""
    out = {"eukdetect2": os.path.join(src_dir, "eukdetect2"), "eukdetect1": os.path.join(src_dir, "eukdetect1"),
           "assemblies": os.path.join(src_dir, "ncbi_assemblies"), "taxdump": os.path.join(src_dir, "taxdump.tar.gz")}
    for folder in (out["eukdetect2"], out["eukdetect1"], out["assemblies"]):
        os.makedirs(folder, exist_ok=True)
    for name, md5 in sources.EUKDETECT_MD5.items():
        dest = os.path.join(out["eukdetect2"], name)
        if not (os.path.isfile(dest) and util.md5(dest) == md5):
            say(f"EukDetect2: fetching {name}")
            sources.download(f"{sources.EUKDETECT_MIRROR.rstrip('/')}/{name}/content", dest, md5)
    if not all(os.path.isfile(os.path.join(out["eukdetect1"], f)) for f in sources.EUKDETECT1_FILES):
        archive = os.path.join(out["eukdetect1"], sources.EUKDETECT1_ARCHIVE[0])
        say("EukDetect1: fetching its archive (1.7 GB)")
        sources.download(sources.EUKDETECT1_URL, archive, sources.EUKDETECT1_ARCHIVE[1])
        with tarfile.open(archive, "r:gz") as tar:
            for member in tar:
                if os.path.basename(member.name) in sources.EUKDETECT1_FILES:
                    member.name = os.path.basename(member.name)
                    tar.extract(member, out["eukdetect1"])
        os.remove(archive)
    for group in sources.ASSEMBLY_GROUPS:
        for kind, base in (("genbank", sources.ASSEMBLIES_URL), ("refseq", sources.REFSEQ_URL)):
            if kind == "refseq" and not refseq:
                continue
            dest = os.path.join(out["assemblies"], f"{kind}_{group}.txt")
            if not os.path.isfile(dest):
                say(f"NCBI: fetching {kind}'s {group} assembly table")
                sources.download(f"{base}/{group}/assembly_summary.txt", dest)
    if not os.path.isfile(out["taxdump"]):
        say("NCBI: fetching the taxdump")
        sources.download(TAXDUMP_URL, out["taxdump"])
    for name in ("EUKARYOME", "ROD"):
        sources.fetch(src_dir, name, sources.download, shutil.which("7za") or shutil.which("7z") or "7za", say)
    return out


def source_paths(src_dir):
    return {"eukdetect2": os.path.join(src_dir, "eukdetect2"), "eukdetect1": os.path.join(src_dir, "eukdetect1"),
            "assemblies": os.path.join(src_dir, "ncbi_assemblies"), "taxdump": os.path.join(src_dir, "taxdump.tar.gz")}


def extra_rows(ncbi, src, assembly_dir, loci, refseq=True):
    """The model species' rows (why model), then the RefSeq species' their sources lack (why refseq)."""
    rows = [dict(r, why="model") for r in species.eukaryotes(species.read_table())]
    if refseq:
        have, _link = sources.species_loci(src["eukdetect2"], src["eukdetect1"], sources.read_locus_map(), ncbi,
                                           set(loci))
        refseq_rows = {}
        for path in sorted(glob.glob(os.path.join(assembly_dir, "refseq_*.txt"))):
            group = os.path.basename(path)[len("refseq_"):-len(".txt")]
            if group in species.REFSEQ_GROUPS:
                for taxid, row in species.refseq_genomes(path, group).items():
                    refseq_rows.setdefault(taxid, row)
        rows += [dict(r, why="refseq") for r in species.gap_rows(refseq_rows, have, ncbi,
                                                                  exclude=[r["ncbi_taxid"] for r in rows])]
    return rows


# --- the genomes' rRNA hook --------------------------------------------------------------------------------------

def rrna_hook(path, job, out_folder, barrnap=None, itsx=None, threads=1):
    """A genome's rRNA into its done folder (rrna.extract_genome; a model species' curated rDNA unit where the genome
    gives no operon) -> its counts."""
    sets = job.get("itsx", "all")  # ITSx's profile sets of the genome's lineage (plan)
    rows = rrna.extract_genome(path, barrnap, itsx, threads, sets)
    if job.get("rdna_unit", "-") != "-" and not any(r["kind"] == "rDNA" for r in rows):
        unit = species.unit_fasta(job["rdna_unit"], os.path.join(out_folder, job["rdna_unit"] + ".fna"))
        if unit:
            rows += [dict(r, contig=job["rdna_unit"])
                     for r in rrna.extract_genome(unit, barrnap, itsx, threads, sets)]
            os.remove(unit)
    return rrna.save(out_folder, rows)


# --- stages ------------------------------------------------------------------------------------------------------

def clades(spec):
    """--species: NCBI taxids, comma-separated or a file of one per line (# comments) -> a set, or None (all)."""
    if not spec:
        return None
    if os.path.isfile(spec):
        with open(spec) as fh:
            spec = ",".join(line.split("#")[0].split("\t")[0].strip() for line in fh)
    return {int(t) for t in spec.split(",") if t.strip()}


def inside_of(ncbi, wanted):
    """A test of a taxid: inside one of the clades `wanted` (None: every taxid is)."""
    if wanted is None:
        return lambda t: True

    def inside(t):
        if t in (None, "", "-"):
            return False
        t = int(t)
        return any(a in wanted for a in ncbi.ancestors(ncbi.merged.get(t, t)))
    return inside


def subset(collected, inside):
    """collect()'s result limited to the species inside the clades of --species."""
    kept = {a: r for a, r in collected["genomes"].items() if inside(r["species_taxid"])}
    return dict(genomes=kept, copies={a: c for a, c in collected["copies"].items() if a in kept},
                features=[f for f in collected["features"] if f["accession"] in kept],
                taxids={t for t in collected["taxids"] if inside(t)})


def plan(opts, src, ncbi, assemblies):
    loci = sources.read_loci_table()
    rows = [r for r in extra_rows(ncbi, src, src["assemblies"], loci, not opts.no_refseq)
            if opts.inside(r["ncbi_taxid"])]
    job_map = genomes.jobs(ncbi, src, loci, sources.read_locus_map(), assemblies, rows, [sources.PLACED_DIR])
    job_map = {a: j for a, j in job_map.items() if opts.inside(j["species_taxid"])}
    for j in job_map.values():
        j["itsx"] = rrna.itsx_groups(ncbi, j["species_taxid"])
    say(f"genome set: {len(job_map)} genomes ({dict(genomes.why_counts(job_map.values()))})")
    return job_map, rows


def run_genomes(opts, job_map, staged_only):
    settings = dict(work=opts.work, genomes=opts.genomes, miniprot=opts.miniprot, threads=opts.threads_per_genome,
                    staged_only=staged_only,
                    rrna=functools.partial(rrna_hook, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"),
                                           threads=opts.threads_per_genome) if shutil.which("barrnap") else None)
    todo = [(j, settings) for a, j in sorted(job_map.items()) if not genomes.is_done(opts.work, a)]
    say(f"genomes to run: {len(todo)} of {len(job_map)}")
    waiting, failed = [], []
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for i, (acc, summary, why) in enumerate(pool.imap_unordered(genomes.run, todo), 1):
            if why:
                (waiting if why == "not fetched yet" else failed).append((acc, why))
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                say(f"  {i}/{len(todo)} genomes run ({len(failed)} failed, {len(waiting)} not fetched yet)")
    return waiting, failed


def published(opts, src_dir, ncbi):
    """The published reads (EUKARYOME, ROD) with their regions and NCBI labels -> records [(name, sequence, row)],
    their regions kept in WORK/published/<source>.tsv (a rerun reads them)."""
    out = []
    index = rrna.name_index(ncbi)
    reads = []
    rod = sources.ready(src_dir, "ROD")
    if rod:
        for r in sources.rod_operons(rod):
            if r.get("taxid", "").isdigit() and r.get("sequence") and opts.inside(r["taxid"]):
                t = int(r["taxid"])
                reads.append((f"ROD_{r['seqid']}", r["sequence"].upper(), t, "ROD", "rod_taxid",
                              r.get("assembly_id") or "-"))
    euk = sources.ready(src_dir, "EUKARYOME")
    if euk:
        bad = sources.eukaryome_compromised(euk)
        for acc, label, seq in sources.eukaryome_long(euk):
            if acc in bad:
                continue
            t, how = rrna.eukaryome_taxid(label, index, ncbi)
            if t and opts.inside(t):
                reads.append((f"EUK_{acc}", seq, t, "EUKARYOME", how, "-"))
    cache_dir = os.path.join(opts.work, "published")
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, "regions.tsv")
    known = {}
    if os.path.isfile(cache):
        with open(cache) as fh:
            for line in fh:
                name, regions, complete = line.rstrip("\n").split("\t")
                known[name] = (regions, complete)
    # by ITSx's profile sets, so that a batch's reads share theirs (an ITSx run per distinct set)
    todo = sorted(((n, s, rrna.itsx_groups(ncbi, t)) for n, s, t, *_ in reads if n not in known),
                  key=lambda r: (r[2], r[0]))
    if todo and shutil.which("barrnap"):
        say(f"published reads: {len(todo)} to annotate (barrnap, ITSx)")
        batches = [todo[i:i + 2000] for i in range(0, len(todo), 2000)]
        work = functools.partial(_annotate_batch, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"))
        with multiprocessing.Pool(max(1, opts.jobs)) as pool, open(cache, "a", newline="\n") as fh:
            for got in pool.imap_unordered(work, batches):
                for name, (regions, complete) in got.items():
                    known[name] = (regions, complete)
                    fh.write(f"{name}\t{regions}\t{complete}\n")
    for name, seq, taxid, source, how, accession in reads:
        regions, complete = known.get(name, ("-", "f"))
        regions_list = None if regions == "-" else [(p.rpartition(":")[0], int(p.rpartition(":")[2].split("-")[0]),
                                                     int(p.rpartition(":")[2].split("-")[1])) for p in regions.split(",")]
        out.append((name, seq, dict(kind="rDNA" if regions_list else "18S", accession=accession, source=source,
                                    taxid=taxid, species_taxid=ncbi.species_of(taxid), regions_list=regions_list,
                                    lsu_complete=complete == "t", label_from=how)))
    return out


def _annotate_batch(batch, barrnap=None, itsx=None):
    """[(name, sequence, ITSx's sets)] -> {name: (regions text, lsu_complete t/f)}."""
    found = rrna.published_regions([(n, s) for n, s, _t in batch], barrnap, itsx, sets={n: t for n, _s, t in batch})
    return {name: (rrna.regions_text(found[name][0]) if name in found else "-",
                   "t" if name in found and found[name][1] else "f") for name, _s, _t in batch}


def genome_rrna_records(opts, job_map, ncbi, assemblies):
    """Every done genome's rRNA as records, labelled by its assembly's taxid."""
    taxid_of = {r["accession"]: (r["taxid"], r["species_taxid"]) for r in assemblies}
    out = []
    for acc, job in sorted(job_map.items()):
        rows = rrna.load(genomes.done_dir(opts.work, acc))
        if not rows:
            continue
        taxid, sp = taxid_of.get(acc, (job.get("species_taxid"), job.get("species_taxid")))
        sp = sp or job.get("species_taxid")
        if not sp:
            continue
        out += rrna.genome_records(acc, rows, taxid or sp, ncbi.species_of(int(sp)))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", required=True, help="the build's work folder (sources, placements, extractions)")
    p.add_argument("--out", required=True, help="where EUKDB_<version>/ and genomes/ go")
    p.add_argument("--stage", choices=("all", "fetch", "process", "package"), default="all")
    p.add_argument("--version", default=datetime.date.today().strftime("%Y.%m"), help="the package's version (YYYY.MM)")
    p.add_argument("--keep_gb", type=float, default=1000.0,
                   help="the genomes kept in OUT/genomes once the package is written, GB (default 1000; the build may "
                        "hold more while it works)")
    p.add_argument("--no_refseq", action="store_true",
                   help="leave out the RefSeq species EukDetect lacks or holds on fewer than 20 loci")
    p.add_argument("--species", help="a partial build (a trial run): only the species inside these NCBI taxids' clades, "
                                     "comma-separated or a file of one per line; the manifest records them")
    p.add_argument("--miniprot", default="miniprot")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("-j", "--jobs", type=int, default=None, help="genomes at a time (default threads / 4)")
    opts = p.parse_args(argv)
    opts.jobs = opts.jobs or max(1, opts.threads // 4)
    opts.threads_per_genome = max(1, opts.threads // opts.jobs)
    opts.genomes = os.path.join(opts.out, "genomes")
    os.makedirs(opts.work, exist_ok=True)
    os.makedirs(opts.genomes, exist_ok=True)
    src_dir = os.path.join(opts.work, "sources")
    if opts.stage in ("all", "fetch"):
        fetch_sources(src_dir, not opts.no_refseq)
    src = source_paths(src_dir)
    ncbi = ncbi_mod.Ncbi(src["taxdump"])
    wanted = clades(opts.species)
    opts.inside = inside_of(ncbi, wanted)
    assemblies = sources.read_assemblies(sorted(glob.glob(os.path.join(src["assemblies"], "genbank_*.txt"))) +
                                         sorted(glob.glob(os.path.join(src["assemblies"], "refseq_*.txt"))))
    job_map, rows = plan(opts, src, ncbi, assemblies)
    if opts.stage == "fetch":
        todo = [j for a, j in sorted(job_map.items())
                if not os.path.isfile(genomes.genome_file(opts.genomes, a)) and not genomes.is_done(opts.work, a)]
        say(f"fetching {len(todo)} genomes into {opts.genomes}")
        failed = 0
        for i, j in enumerate(todo, 1):
            why = util.fetch_file(j["url"], genomes.genome_file(opts.genomes, j["accession"]))
            failed += bool(why)
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                say(f"  {i}/{len(todo)} fetched ({failed} failed)")
        return 0
    if opts.stage in ("all", "process", "package"):
        if not shutil.which("barrnap"):
            say("WARNING: no barrnap: no rRNA records (install barrnap with nhmmer)")
        elif not shutil.which("ITSx"):
            say("WARNING: no ITSx: the ITS1/5.8S/ITS2/28S borders are barrnap's alone")
        else:
            why = rrna.check_itsx(shutil.which("ITSx"))
            if why:  # a damaged ITSx_db gives silently wrong borders: stop before any work
                sys.exit(f"{why}\n{shutil.which('ITSx')}: reinstall ITSx (bioconda's itsx) or give it its whole "
                         "ITSx_db/HMMs folder")
    if opts.stage in ("all", "process"):
        waiting, failed = run_genomes(opts, job_map, staged_only=opts.stage == "process")
        if shutil.which(opts.miniprot):
            only = []
            if wanted is not None:
                link, _records = sources.read_link(os.path.join(src["eukdetect2"], "busco_taxid_genome_link.txt"))
                only = ["--only", ",".join(sorted(g for g, (t, _grp) in link.items() if opts.inside(t))) or "-"]
            exons.main(["--eukdetect", src["eukdetect2"], "--eukdetect1", src["eukdetect1"], "--work",
                        os.path.join(opts.work, "exons"), "--out", os.path.join(opts.work, "exons"),
                        "-j", str(opts.threads), "--miniprot", opts.miniprot] + only)
        if waiting or failed:
            say(f"{len(waiting)} genomes not fetched yet, {len(failed)} failed (e.g. {(waiting + failed)[:3]}): the "
                "package waits for them (run again; --stage fetch first on a cluster)")
            return 0 if opts.stage == "process" else 1
    if opts.stage in ("all", "package"):
        collected = package.collect(ncbi, src, assemblies, sources.read_loci_table(), sources.read_locus_map(), rows,
                                    opts.work, [sources.PLACED_DIR],
                                    [os.path.join(opts.work, "exons", exons.EXONS)])
        if wanted is not None:
            collected = subset(collected, opts.inside)
        records = genome_rrna_records(opts, job_map, ncbi, assemblies) + published(opts, src_dir, ncbi)
        records = rrna.finish(records, ncbi, shutil.which("vsearch"), opts.threads)
        inputs = {"eukdetect2": sources.EUKDETECT_MD5, "eukdetect1": sources.EUKDETECT1_ARCHIVE[1],
                  "taxdump_md5": util.md5(src["taxdump"]),
                  "sources": {n: sources.state_of(src_dir, n).get("version") for n in ("EUKARYOME", "ROD")},
                  "assembly_tables": {os.path.basename(p): util.md5(p)
                                      for p in sorted(glob.glob(os.path.join(src["assemblies"], "*.txt")))}}
        if wanted is not None:
            inputs["partial"] = {"species_clades": sorted(wanted)}  # a trial package: these clades only
        # every genome is extracted: the kept ones cut back to the budget, then the package names those kept
        firsts = {a for r in species.read_table() for a in species.accessions(r)}
        rows_kept = genomes.prune(opts.genomes, opts.work, opts.keep_gb * 1e9, set(), firsts)
        kept = sum(s for _a, s, st in rows_kept if st != "deleted")
        say(f"genomes kept in {opts.genomes}: {kept / 1e9:.0f} GB (--keep_gb {opts.keep_gb:g}), "
            f"{sum(st == 'deleted' for *_x, st in rows_kept)} deleted")
        out = package.write(os.path.join(opts.out, f"EUKDB_{opts.version}"), opts.version, ncbi, collected,
                            sources.LOCI_TABLE, records, opts.genomes, inputs)
        say(f"EUKDB: {out} ({json.dumps(package.read(out)['counts'])})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
