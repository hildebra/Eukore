#!/usr/bin/env python3
"""Build EUKDB, the standalone eukaryote reference (FORMAT.md; since 2026-11-05).

    python3 eukdb_build.py --work WORK --out OUT -t 16                 # all stages on a node with internet
    python3 eukdb_build.py --work WORK --out OUT --stage fetch -t 8    # a cluster: the node with internet,
    python3 eukdb_build.py --work WORK --out OUT --stage process -t 32 #   then any node (nothing fetched),
    python3 eukdb_build.py --work WORK --out OUT --stage package       #   then the package

Stages:
- fetch: the sources into WORK/sources (NCBI's GenBank and RefSeq assembly tables of the --groups, NCBI's taxdump,
  BUSCO's eukaryota_odb12, EUKARYOME 2.1, ROD v1.2), then each species' genome (genomes.jobs) into OUT/genomes;
- process: each genome's euk30 loci called and its rRNA extracted (WORK/done/<accession>/); a genome not fetched yet
  waits for the next run;
- package: the copies' checks across species, the published rRNA reads' regions (WORK/published), OUT/genomes cut
  back to --keep_gb (genomes.prune: the model species' and the smallest kept first), then OUT/EUKDB_<version>/.
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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eukdb import genomes, ncbi as ncbi_mod, odb as odb_mod, package, rrna, sources, species, util  # noqa: E402

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
NO_ASSEMBLY = os.path.join(sources.DATA, "no_assembly.fna.gz")
PUBLISHED_BATCH = 250  # published reads a barrnap + ITSx run


def say(text):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {text}", flush=True)


# --- sources -----------------------------------------------------------------------------------------------------

def fetch_sources(src_dir, groups):
    """Every source into src_dir (kept when there and whole)."""
    os.makedirs(os.path.join(src_dir, "ncbi_assemblies"), exist_ok=True)
    for group in groups:
        for kind, base in (("genbank", sources.ASSEMBLIES_URL), ("refseq", sources.REFSEQ_URL)):
            dest = os.path.join(src_dir, "ncbi_assemblies", f"{kind}_{group}.txt")
            if not os.path.isfile(dest):
                say(f"NCBI: fetching {kind}'s {group} assembly table")
                sources.download(f"{base}/{group}/assembly_summary.txt", dest)
    taxdump = os.path.join(src_dir, "taxdump.tar.gz")
    if not os.path.isfile(taxdump):
        say("NCBI: fetching the taxdump")
        sources.download(TAXDUMP_URL, taxdump)
    odb_mod.fetch(src_dir, sources.download, say)
    for name in ("EUKARYOME", "ROD"):
        sources.fetch(src_dir, name, sources.download, shutil.which("7za") or shutil.which("7z") or "7za", say)


def tables(src_dir, kind, groups):
    return [p for g in groups for p in [os.path.join(src_dir, "ncbi_assemblies", f"{kind}_{g}.txt")]
            if os.path.isfile(p)]


# --- the genomes' rRNA hook --------------------------------------------------------------------------------------

def rrna_hook(path, job, out_folder, barrnap=None, itsx=None, threads=1):
    """A genome's rRNA into its done folder (rrna.extract_genome; a model species' curated rDNA unit where the genome
    gives no operon) -> its counts."""
    sets = job.get("itsx", "all")  # ITSx's profile sets of the genome's lineage (genomes.jobs)
    rows = rrna.extract_genome(path, barrnap, itsx, threads, sets)
    if job.get("rdna_unit", "-") != "-" and not any(r["kind"] == "rDNA" for r in rows):
        unit = species.unit_fasta(job["rdna_unit"], os.path.join(out_folder, job["rdna_unit"] + ".fna"))
        if unit:
            rows += [dict(r, contig=job["rdna_unit"])
                     for r in rrna.extract_genome(unit, barrnap, itsx, threads, sets)]
            os.remove(unit)
    return rrna.save(out_folder, rows)


# --- stages ------------------------------------------------------------------------------------------------------

def plan(opts, src_dir, ncbi):
    """The genome set -> (jobs, GenBank rows)."""
    assemblies = sources.read_assemblies(tables(src_dir, "genbank", opts.groups))
    refseq_rows = sources.read_assemblies(tables(src_dir, "refseq", opts.groups))
    no_assembly = genomes.no_assembly_genomes(NO_ASSEMBLY, os.path.join(opts.work, "no_assembly"), ncbi)
    job_map = genomes.jobs(ncbi, assemblies, refseq_rows, species.eukaryotes(species.read_table()), opts.inside,
                           no_assembly)
    say(f"genome set: {len(job_map)} genomes ({dict(genomes.why_counts(job_map.values()))}); "
        f"{sum(j['size'] or 0 for j in job_map.values()) / 1e9:.0f} Gb")
    return job_map, assemblies


def settings_of(opts, src_dir, ncbi, staged_only):
    folder = odb_mod.folder_of(src_dir)
    odb = odb_mod.Odb(folder, sources.read_loci_table())
    hmm_db = odb.hmm_db(os.path.join(opts.work, "euk30_loci.hmm"))
    tiara = None if opts.no_tiara else shutil.which(opts.tiara)
    return dict(work=opts.work, genomes=opts.genomes, miniprot=opts.miniprot, hmmsearch=opts.hmmsearch,
                threads=opts.threads_per_genome, staged_only=staged_only, odb=folder, loci=odb.loci, hmm_db=hmm_db,
                odb_lineages=odb_mod.species_lineages(odb, ncbi), tiara=tiara,
                rrna=functools.partial(rrna_hook, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"),
                                       threads=opts.threads_per_genome) if shutil.which("barrnap") else None)


def run_genomes(opts, job_map, settings):
    todo = [(j, settings) for a, j in sorted(job_map.items()) if not genomes.is_done(opts.work, a)]
    say(f"genomes to run: {len(todo)} of {len(job_map)}" + ("" if settings["tiara"] else " (no tiara: windows not "
                                                                                          "classified)"))
    waiting, failed = [], []
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for i, (acc, _summary, why) in enumerate(pool.imap_unordered(genomes.run, todo), 1):
            if why:
                (waiting if why == "not fetched yet" else failed).append((acc, why))
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                say(f"  {i}/{len(todo)} genomes run ({len(failed)} failed, {len(waiting)} not fetched yet)")
    return waiting, failed


def published(opts, src_dir, ncbi):
    """The published reads (EUKARYOME, ROD) with their regions and NCBI labels -> records [(name, sequence, row)],
    their regions kept in WORK/published/regions.tsv (a rerun reads them)."""
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
        # small batches over every thread: barrnap (nhmmer) runs one thread a batch and dominates (~0.5 s a read)
        batches = [todo[i:i + PUBLISHED_BATCH] for i in range(0, len(todo), PUBLISHED_BATCH)]
        work = functools.partial(_annotate_batch, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"))
        with multiprocessing.Pool(max(1, opts.threads)) as pool, open(cache, "a", newline="\n") as fh:
            for got in pool.imap_unordered(work, batches):
                for name, (regions, complete) in got.items():
                    known[name] = (regions, complete)
                    fh.write(f"{name}\t{regions}\t{complete}\n")
                fh.flush()
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


def genome_rrna_records(opts, job_map, ncbi):
    """Every done genome's rRNA as records, labelled by its assembly's taxid."""
    out = []
    for acc, job in sorted(job_map.items()):
        rows = rrna.load(genomes.done_dir(opts.work, acc))
        if rows:
            out += rrna.genome_records(acc, rows, job.get("taxid") or job["species_taxid"],
                                       ncbi.species_of(int(job["species_taxid"])))
    return out


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


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", required=True, help="the build's work folder (sources, calls, extractions)")
    p.add_argument("--out", required=True, help="where EUKDB_<version>/ and genomes/ go")
    p.add_argument("--stage", choices=("all", "fetch", "process", "package"), default="all")
    p.add_argument("--version", default=datetime.date.today().strftime("%Y.%m"), help="the package's version (YYYY.MM)")
    p.add_argument("--groups", default=",".join(sources.ASSEMBLY_GROUPS),
                   help=f"NCBI's assembly groups whose species the package holds (default all: "
                        f"{','.join(sources.ASSEMBLY_GROUPS)}; e.g. fungi,protozoa first)")
    p.add_argument("--keep_gb", type=float, default=1000.0,
                   help="the genomes kept in OUT/genomes once the package is written, GB (default 1000; the build may "
                        "hold more while it works)")
    p.add_argument("--species", help="a partial build (a trial run): only the species inside these NCBI taxids' clades, "
                                     "comma-separated or a file of one per line; the manifest records them")
    p.add_argument("--miniprot", default="miniprot")
    p.add_argument("--hmmsearch", default="hmmsearch")
    p.add_argument("--fcs_gx_db", help="NCBI FCS-GX's database (the gxdb folder's prefix, ~470 GB, best copied to "
                                       "/dev/shm; a node of ~512 GB): the package stage screens each genome's windows "
                                       "around its copies (copies.tsv's fcs)")
    p.add_argument("--run_gx", default="run_gx", help="FCS-GX's runner (bioconda's ncbi-fcs-gx)")
    p.add_argument("--tiara", default="tiara", help="tiara, which classifies the windows around the copies (default: "
                                                    "tiara on the PATH, else none; envs/tiara.yaml's, as it needs a "
                                                    "Python of 3.9 or older)")
    p.add_argument("--no_tiara", action="store_true", help="do not classify the windows around the copies")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("-j", "--jobs", type=int, default=None, help="genomes at a time (default threads / 4)")
    opts = p.parse_args(argv)
    opts.groups = [g for g in opts.groups.split(",") if g]
    unknown = set(opts.groups) - set(sources.ASSEMBLY_GROUPS)
    if unknown:
        p.error(f"--groups: {sorted(unknown)} (known: {', '.join(sources.ASSEMBLY_GROUPS)})")
    opts.jobs = opts.jobs or max(1, opts.threads // 4)
    opts.threads_per_genome = max(1, opts.threads // opts.jobs)
    opts.genomes = os.path.join(opts.out, "genomes")
    os.makedirs(opts.work, exist_ok=True)
    os.makedirs(opts.genomes, exist_ok=True)
    src_dir = os.path.join(opts.work, "sources")
    if opts.stage in ("all", "fetch"):
        fetch_sources(src_dir, opts.groups)
    ncbi = ncbi_mod.Ncbi(os.path.join(src_dir, "taxdump.tar.gz"))
    wanted = clades(opts.species)
    opts.inside = inside_of(ncbi, wanted)
    job_map, assemblies = plan(opts, src_dir, ncbi)
    if opts.stage == "fetch":
        todo = [j for a, j in sorted(job_map.items()) if j["url"] and not j.get("path")
                and not os.path.isfile(genomes.genome_file(opts.genomes, a)) and not genomes.is_done(opts.work, a)]
        say(f"fetching {len(todo)} genomes into {opts.genomes}")
        failed = 0
        for i, j in enumerate(todo, 1):
            why = util.fetch_file(j["url"], genomes.genome_file(opts.genomes, j["accession"]))
            failed += bool(why)
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                say(f"  {i}/{len(todo)} fetched ({failed} failed)")
        return 0
    if shutil.which("barrnap") and shutil.which("ITSx"):
        why = rrna.check_itsx(shutil.which("ITSx"))
        if why:  # a damaged ITSx_db gives silently wrong borders: stop before any work
            sys.exit(f"{why}\n{shutil.which('ITSx')}: reinstall ITSx (bioconda's itsx) or give it its whole "
                     "ITSx_db/HMMs folder")
    elif not shutil.which("barrnap"):
        say("WARNING: no barrnap: no rRNA records (install barrnap with nhmmer)")
    else:
        say("WARNING: no ITSx: the ITS1/5.8S/ITS2/28S borders are barrnap's alone")
    for tool in (opts.miniprot, opts.hmmsearch):
        if not shutil.which(tool):
            sys.exit(f"{tool}: not found (bioconda's miniprot and hmmer)")
    if opts.stage in ("all", "process"):
        settings = settings_of(opts, src_dir, ncbi, staged_only=opts.stage == "process")
        waiting, failed = run_genomes(opts, job_map, settings)
        if waiting or failed:
            say(f"{len(waiting)} genomes not fetched yet, {len(failed)} failed (e.g. {(waiting + failed)[:3]}): the "
                "package waits for them (run again; --stage fetch first on a cluster)")
            return 0 if opts.stage == "process" else 1
    if opts.stage in ("all", "package") and opts.fcs_gx_db:
        todo = [a for a in sorted(job_map) if genomes.is_done(opts.work, a)
                and not os.path.isfile(os.path.join(genomes.done_dir(opts.work, a), "fcs.tsv"))]
        say(f"FCS-GX: {len(todo)} genomes' windows to screen ({opts.fcs_gx_db})")
        failed = 0
        for i, acc in enumerate(todo, 1):
            got = genomes.fcs_screen(opts.work, acc, job_map[acc]["taxid"], opts.fcs_gx_db, opts.run_gx, opts.threads)
            failed += got is None
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                say(f"  {i}/{len(todo)} screened ({failed} failed)")
    if opts.stage in ("all", "package"):
        vsearch = shutil.which("vsearch")
        collected = package.collect(ncbi, job_map, assemblies, opts.work, vsearch, opts.threads)
        records = genome_rrna_records(opts, job_map, ncbi) + published(opts, src_dir, ncbi)
        records = rrna.finish(records, ncbi, vsearch, opts.threads)
        inputs = {"odb": odb_mod.state_of(src_dir), "taxdump_md5": util.md5(os.path.join(src_dir, "taxdump.tar.gz")),
                  "sources": {n: sources.state_of(src_dir, n).get("version") for n in ("EUKARYOME", "ROD")},
                  "assembly_tables": {os.path.basename(p): util.md5(p)
                                      for p in sorted(glob.glob(os.path.join(src_dir, "ncbi_assemblies", "*.txt")))},
                  "groups": opts.groups, "no_assembly": util.md5(NO_ASSEMBLY) if os.path.isfile(NO_ASSEMBLY) else None}
        if wanted is not None:
            inputs["partial"] = {"species_clades": sorted(wanted)}  # a trial package: these clades only
        # every genome is called and extracted: the kept ones cut back to the budget, then the package names those kept
        firsts = {a for r in species.read_table() for a in species.accessions(r)}
        rows_kept = genomes.prune(opts.genomes, opts.work, opts.keep_gb * 1e9, set(), firsts)
        kept = sum(s for _a, s, st in rows_kept if st != "deleted")
        say(f"genomes kept in {opts.genomes}: {kept / 1e9:.1f} GB (--keep_gb {opts.keep_gb:g}), "
            f"{sum(st == 'deleted' for *_x, st in rows_kept)} deleted")
        out = package.write(os.path.join(opts.out, f"EUKDB_{opts.version}"), opts.version, ncbi, collected,
                            sources.LOCI_TABLE, records, opts.genomes, inputs)
        say(f"EUKDB: {out} ({json.dumps(package.read(out)['counts'])})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
