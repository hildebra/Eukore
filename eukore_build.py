#!/usr/bin/env python3
"""Build Eukore, the standalone eukaryote reference (docs/FORMAT.md; since 2026-11-05).

    python3 eukore_build.py --build_dir WORK --out OUT -t 16                 # all stages on a node with internet
    python3 eukore_build.py --build_dir WORK --out OUT --stage fetch -t 8    # a cluster: the node with internet,
    python3 eukore_build.py --build_dir WORK --out OUT --stage process -t 32 #   then any node (nothing fetched),
    python3 eukore_build.py --build_dir WORK --out OUT --stage package       #   then the package
    python3 eukore_build.py --build_dir WORK --out OUT --stage process,package -t 32  # or both in one run

Stages:
- fetch: the sources into WORK/sources (NCBI's GenBank and RefSeq assembly tables of the --groups, NCBI's taxdump,
  BUSCO's eukaryota_odb12, EUKARYOME 2.1, ROD v1.2, RefSeq's organelle genomes, NCBI's FCS-GX and ANI reports, which
  rank the genomes: eukore/quality.py, PR2's and KSGP's taxonomies, Eukore's marker set: KOfam's profiles, Swiss-Prot's
  and NCBI's proteins, eukore/markers.py), then each species' genome (genomes.jobs) into OUT/genomes (--genomes);
- process: each genome's BUSCO genes called (a genome searched for fewer than those called is run again) and its rRNA
  extracted (WORK/done/<accession>/); a genome not fetched yet waits for the next run; then the organelle genes
  (mitochondrial and plastid rRNA, mitochondrial COX1) of RefSeq's organelle genomes, and of the assemblies' own
  organelle molecules where RefSeq has none (WORK/organelle);
- package: the copies' checks across species, the published rRNA reads' regions (WORK/published), the organelle genes
  of any assembly not searched yet, OUT/genomes (--genomes) cut back to --keep_gb where given (genomes.prune: the
  model species' and the smallest kept first; without it no genome file is deleted), then OUT/Eukore_<version>/ and,
  unless --no_release, OUT/release/: the website's data (website/) and the downloads, BGZF (eukore/release.py);
- release: OUT/release/ again of a finished OUT/Eukore_<version>/.
`all` runs fetch, process and package in turn; stages that follow one another may be given together (since
2026-11-10: `--stage process,package`, a cluster's last job after the process runs of a SLURM array; tools/eukore.sh
submits the chain). Every stage is resumable: a rerun does only what is not done. Every compressed file is BGZF
(gzip in blocks: bgzip -@ decompresses it on many threads), compressed on -t threads.
"""

import argparse
import bisect
import collections
import concurrent.futures
import contextlib
import datetime
import functools
import glob
import hashlib
import json
import multiprocessing
import os
import random
import shlex
import shutil
import socket
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eukore import genomes, ncbi as ncbi_mod, odb as odb_mod, package, rrna, sources, species, util  # noqa: E402
from eukore import annotation, benchmarks, labels, metadata as metadata_mod, organelle as organelle_mod  # noqa: E402
from eukore import groups as groups_mod, intraspecific as intraspecific_mod, markers as markers_mod, quality  # noqa: E402
from eukore import placement, primers as primers_mod, reads as reads_mod, regions as regions_mod  # noqa: E402
from eukore import refscreen, release as release_mod  # noqa: E402
from eukore import __version__  # noqa: E402
from eukore import store  # noqa: E402

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
NO_ASSEMBLY = os.path.join(sources.DATA, "no_assembly.fna.gz")
PUBLISHED_BATCH = 250  # published reads a barrnap + ITSx run
# the published reads' regions, kept in WORK/published: a new name where their borders change (2: ITSx's guessed 5.8S
# borders no longer taken, 2026-11-06; 3: regions on the record, not the unit, and the read's strand, 2026-11-06), so
# that a rerun annotates them again
PUBLISHED_REGIONS = "regions.3.tsv"
READ_BACKGROUND = os.path.join(sources.DATA, "read_background.tsv")  # the read benchmark's bacteria and archaea
FCS_THREADS = 4  # cores of each FCS-GX run (GX_NUM_CORES); -t / FCS_THREADS runs at once
PACK_THREADS = 8  # done folders packed at once (pack_pass: small files read, one written)
LEFTOVER_AGE = 86400  # a download's .part (or a TOGA2 .2bit) untouched this long (seconds) is a killed run's


def say(text):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {text}", flush=True)


# --- sources -----------------------------------------------------------------------------------------------------

def fetch_sources(src_dir, groups, ksgp_tax=None, refresh=False, threads=1):
    """Every source into src_dir (kept when there and whole). refresh (--refresh): those that choose the genomes and
    their gene models fetched again (NCBI's assembly tables, its FCS-GX and ANI reports, the annotation indexes; each
    old one kept where its download fails), so that the plan takes NCBI's newer assemblies and annotations. threads
    (-t): the marker set's hmmsearch (its ~280,000 candidates; the fetch node is often a shared login node)."""
    os.makedirs(os.path.join(src_dir, "ncbi_assemblies"), exist_ok=True)
    for group in groups:
        for kind, base in (("genbank", sources.ASSEMBLIES_URL), ("refseq", sources.REFSEQ_URL)):
            dest = os.path.join(src_dir, "ncbi_assemblies", f"{kind}_{group}.txt")
            if refresh or not os.path.isfile(dest):
                say(f"NCBI: fetching {kind}'s {group} assembly table")
                try:
                    sources.download(f"{base}/{group}/assembly_summary.txt", dest)
                except OSError as e:
                    if not os.path.isfile(dest):
                        raise
                    say(f"WARNING: {kind}'s {group} assembly table not fetched again ({e}): the old one kept")
    taxdump = os.path.join(src_dir, "taxdump.tar.gz")
    if not os.path.isfile(taxdump):
        say("NCBI: fetching the taxdump")
        sources.download(TAXDUMP_URL, taxdump)
    odb_mod.fetch(src_dir, sources.download, say)
    for name in ("EUKARYOME", "ROD"):
        sources.fetch(src_dir, name, sources.download, shutil.which("7za") or shutil.which("7z") or "7za", say)
    # the optional sources each on its own: one that fails (a site moved, down) is said, the others still fetched
    # (until 2026-11-06 one failure stopped the stage, NCBI's quality reports then missing from the genome ranking)
    for what, fetch in (("RefSeq's organelle release", lambda: organelle_mod.fetch(src_dir, sources.download, say)),
                        ("Rfam's SSU model", lambda: regions_mod.fetch(src_dir, sources.download, say)),
                        ("Eukore's marker set (KOfam, Swiss-Prot, NCBI's proteins)",
                         lambda: markers_mod.fetch(src_dir, sources.download, say, shutil.which("hmmsearch") or "hmmsearch",
                                                   threads)),
                        ("ENA's list of TSA transcriptomes",
                         lambda: genomes.fetch_tsa(src_dir, sources.download, say, refresh)),
                        ("NCBI's FCS-GX and ANI reports",
                         lambda: quality.fetch(src_dir, sources.download, say, refresh)),
                        ("PR2's and KSGP's taxonomies", lambda: labels.fetch(src_dir, sources.download, say, ksgp_tax)),
                        ("Ensembl's, VEuPathDB's and TOGA2's annotation indexes",
                         lambda: annotation.fetch_indexes(src_dir, sources.download, say, refresh)),
                        ("the read benchmark's background", lambda: reads_mod.fetch_background(
                            READ_BACKGROUND, os.path.join(src_dir, "read_background"), say)
                         if os.path.isfile(READ_BACKGROUND) else None)):
        try:
            fetch()
        except OSError as e:
            say(f"WARNING: {what} not fetched ({e}): the build goes on without it; a rerun of --stage fetch tries "
                "again")


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
    unknown = {r["species_taxid"] for r in assemblies} - set(ncbi.parent) - set(ncbi.merged)
    if unknown:
        say(f"WARNING: {len(unknown)} species of the assembly tables are not in the taxdump "
            f"(e.g. {sorted(unknown)[:3]}; "
            "a taxdump older than the tables?): left out (selection.tsv: not in the taxdump)")
    if quality.paths(src_dir):
        say(f"NCBI's FCS-GX and ANI reports: {quality.annotate(assemblies, refseq_rows, src_dir)} assemblies flagged")
    else:
        say("WARNING: no NCBI FCS-GX or ANI reports (--stage fetch fetches them): genomes ranked without them")
    no_assembly =genomes.no_assembly_genomes(NO_ASSEMBLY, os.path.join(opts.work, "no_assembly"), ncbi)
    job_map = genomes.jobs(ncbi, assemblies, refseq_rows, species.eukaryotes(species.read_table()), opts.inside,
                           no_assembly)
    # TOGA2's genomes of the vertebrate species NCBI has none of (DNA Zoo's, Bat1K's), with their gene models, and in
    # place of NCBI's where far more contiguous and as long
    toga2 = annotation.read_toga2(src_dir)
    job_map.update(genomes.toga2_jobs(ncbi, toga2, job_map, opts.inside, opts.groups))
    sizes_dir = os.path.join(src_dir, annotation.INDEX_DIR, "assembly_sizes")

    def kept(url, name):  # a small file of the genome comparison, kept in SRC (the fetch stage's: every stage the same)
        path = os.path.join(sizes_dir, name)
        if url and not os.path.isfile(path) and "fetch" in opts.stages:
            os.makedirs(sizes_dir, exist_ok=True)
            util.fetch_file(url, path, gzip_check=False)
        return path
    swaps = genomes.prefer_toga2(
        ncbi, toga2, job_map, lambda j: genomes.assembly_stats(kept(genomes.stats_url(j), f"{j['fetch']}.stats.txt")),
        lambda r: genomes.total_length(kept(genomes.chrom_sizes_url(r), f"{r['accession']}.chrom.sizes")))
    for old, new in swaps.items():
        del job_map[old]
        job_map[new["accession"]] = new
    if swaps:
        say(f"TOGA2: {len(swaps)} species' genomes TOGA2's in place of NCBI's (contig N50 {genomes.TOGA2_BETTER:g} "
            f"times or more, {genomes.TOGA2_LENGTH:.0%} of its length or more)")
    if not getattr(opts, "no_tsa", False):  # since 2026-11-07: transcriptomes of the species with no genome
        tsa = genomes.tsa_jobs(ncbi, genomes.read_tsa(genomes.tsa_path(src_dir)), job_map, opts.inside)
        job_map.update(tsa)
        if tsa:
            say(f"TSA: {len(tsa)} species without a genome have a transcriptome (ENA's TSA sets)")
    say(f"genome set: {len(job_map)} genomes ({dict(genomes.why_counts(job_map.values()))}); "
        f"{sum(j['size'] or 0 for j in job_map.values()) / 1e9:.0f} Gb")
    # a species' genome another assembly takes the place of (TOGA2's; NCBI's newer, after --refresh): its file
    # deleted once the new one is done, with --keep_gb (drop_superseded)
    changed = genomes.supersede(opts.work, job_map)
    if changed:
        say(f"genome set: {len(changed)} species' genomes changed since the last plan (e.g. {changed[:3]}); the old "
            "files are deleted once the new ones are done" + ("" if opts.keep_gb is not None else " (with --keep_gb)"))
    others = annotation.assign(job_map, annotation.read_index(src_dir), veupathdb=bool(opts.veupathdb_key))
    say(f"gene models: {sum(1 for j in job_map.values() if annotation.gff_url(j)):,} genomes may have NCBI's, "
        f"{others:,} VEuPathDB's, Ensembl's or TOGA2's (annotation.py"
        f"{'' if opts.veupathdb_key else '; no VEuPathDB API key: none of its'})")
    return job_map, assemblies


def loci_of(opts, src_dir):
    """The BUSCOs called: every one of the dataset (since 2026-11-05: the user picks among them, protal its euk30),
    or --loci's: a build stops where none of them is the dataset's, and says which are not (until 2026-11-06 a file of
    another dataset's ids gave a package of no BUSCOs)."""
    every = odb_mod.all_loci(odb_mod.folder_of(src_dir))
    if not opts.loci:
        return every
    with open(opts.loci) as fh:
        wanted = {line.split("\t")[0].strip() for line in fh if line.strip() and not line.startswith(("#", "locus"))}
    unknown = sorted(wanted - set(every))
    if len(unknown) == len(wanted):
        sys.exit(f"--loci {opts.loci}: none of its {len(wanted)} BUSCOs is in {odb_mod.LINEAGE} (e.g. {unknown[:3]})")
    if unknown and not getattr(opts, "loci_said", False):
        say(f"WARNING: --loci: {len(unknown)} of its BUSCOs are not in {odb_mod.LINEAGE}, left out: {unknown[:5]}")
        opts.loci_said = True
    return [locus for locus in every if locus in wanted]


def markers_of(opts, src_dir):
    """Eukore's own marker dataset (markers.py; since 2026-11-07), or None (--no_markers, or not built: the fetch
    stage builds it)."""
    return None if getattr(opts, "no_markers", False) else markers_mod.ready(src_dir)


def called_loci(opts, src_dir):
    """The loci called: the BUSCOs (loci_of) and the markers, as every stage counts a genome done for them (since
    2026-11-08: the fetch stage and the package's warning counted the BUSCOs alone, so that a genome pruned before the
    markers came was never fetched again for them)."""
    marker_folder = markers_of(opts, src_dir)
    return loci_of(opts, src_dir) + (odb_mod.all_loci(marker_folder) if marker_folder else [])


def odb_of(opts, src_dir):
    """The loci called (called_loci), in one dataset object."""
    marker_folder = markers_of(opts, src_dir)
    return odb_mod.Odb(odb_mod.folder_of(src_dir), called_loci(opts, src_dir), [marker_folder] if marker_folder else [])


def loci_items(opts, src_dir):
    """The loci called as package.write takes them: locus, its name, cutoff, the HMM's length, and its set (busco, or
    markers: Eukore's own, since 2026-11-07)."""
    odb = odb_of(opts, src_dir)
    return [dict(locus=l, name=odb.names.get(l, "-"), cutoff=odb.cutoff[l], hmm_length=odb.length[l],
                 set="markers" if l in odb.markers else "busco") for l in odb.loci]


def settings_of(opts, src_dir, ncbi, staged_only):
    folder = odb_mod.folder_of(src_dir)
    odb = odb_of(opts, src_dir)
    hmm_db = odb.hmm_db(os.path.join(opts.work, f"loci_{odb.digest()}.hmm"))  # a new locus or decoy: a new file
    tiara = None if opts.no_tiara else shutil.which(opts.tiara)
    return dict(work=opts.work, genomes=opts.genomes, miniprot=opts.miniprot, hmmsearch=opts.hmmsearch,
                threads=opts.threads_per_genome, staged_only=staged_only, odb=folder, loci=odb.loci, hmm_db=hmm_db,
                markers=markers_of(opts, src_dir),
                odb_lineages=odb_mod.species_lineages(odb, ncbi), odb_own=odb_mod.own_taxids(odb, ncbi), tiara=tiara,
                # the clade hold-out (since 2026-11-09): the dataset species' ranks, the ranks left out in turn
                odb_ranks={t: {r: x for r, x in ncbi.rank_ids(t).items() if r in genomes.DEPTH_RANKS}
                           for t in odb.species()},
                held_out_ranks=opts.held_out_ranks,
                rescue=not opts.no_rescue and bool(shutil.which("esl-translate")),
                esl_translate=shutil.which("esl-translate"), sourmash=sourmash_of(opts),
                veupathdb_key=opts.veupathdb_key, vsearch=shutil.which("vsearch"), scratch=opts.scratch,
                pack=not opts.done_folders,  # a genome's results one file (store.py; since 2026-11-10)
                rrna=functools.partial(rrna_hook, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"),
                                       threads=opts.threads_per_genome) if shutil.which("barrnap") else None)


def sourmash_of(opts):
    """sourmash, which sketches each genome (metadata.SKETCH), or None (--no_sketches, or not on the PATH)."""
    return None if opts.no_sketches else shutil.which(opts.sourmash)


def describe_pass(opts, job_map):
    """The metadata of each done genome still on disk that lacks it (metadata.py: a genome run before 2026-11-07, or
    whose measure failed), before the prune deletes the genomes -> how many were described, failed. Each under its
    genome's claim (describe_one; since 2026-11-10), so that several builds on a WORK (a SLURM array's tasks, after
    their genomes) share them out."""
    sourmash = sourmash_of(opts)
    todo = []
    for acc, j in sorted(job_map.items()):
        folder, path = genomes.done_dir(opts.work, acc), genomes.genome_file(opts.genomes, acc)
        if not j.get("path") and genomes.is_done(opts.work, acc) and os.path.isfile(path) \
                and metadata_mod.missing(folder, sourmash):
            todo.append((acc, path, opts.work, j.get("fetch", acc), sourmash))
    if not todo:
        return 0, 0
    say(f"assembly metadata: {len(todo)} genomes to describe" + ("" if sourmash else " (no sourmash: no sketches)"))
    failed = []
    # every --report_minutes (until 2026-11-08 at each 5% of the genomes: hours apart, or minutes)
    progress = progress_of(opts, "genomes described", [job_map[t[0]].get("size") or 0 for t in todo])
    # a genome's scan on one core and its sketch (sourmash) on another: -t / 2 at once (until 2026-11-08 -j, -t / 4)
    workers = max(1, opts.threads // (2 if sourmash else 1))
    with util.process_pool(workers) as pool:
        for acc, why in util.unordered(pool, describe_one, todo, progress,
                                       lambda: running(opts, progress, workers)):
            if why:
                failed.append((acc, why))
            progress.step(job_map[acc].get("size") or 0, failed=bool(why),
                          note=f" (e.g. {failed[0]})" if failed else "")
    return len(todo) - len(failed), len(failed)


def describe_one(args):
    """A pool's task: (accession, genome path, WORK, fetched accession, sourmash) -> (accession, why it failed or
    None): the genome's metadata it lacks (metadata.describe) made in a temporary folder and put with its results
    (store.put: a folder's or a packed genome's), under its claim. A genome another build holds is left to it; what
    is missing is looked at once the genome is claimed (another build's task may have made it since the plan)."""
    acc, path, work, fetched, sourmash = args
    folder = genomes.done_dir(work, acc)
    try:
        with genomes.claimed(work, acc) as held:
            want = metadata_mod.missing(folder, sourmash) if held else None
            if not want:
                return acc, None
            with tempfile.TemporaryDirectory() as tmp:
                why = metadata_mod.describe(path, tmp, acc, fetched, sourmash, want)
                made = {n: os.path.join(tmp, n) for n in os.listdir(tmp) if not n.startswith(".")}
                if made:
                    store.put(folder, made)
            return acc, why
    except Exception as e:  # noqa: BLE001  a genome's failure, never the pool's
        return acc, f"metadata: {type(e).__name__}: {e}"


def pack_pass(opts):
    """The done folders of the WORK packed, each into its one file (store.pack; since 2026-11-10: a build before then
    left a folder of 8-14 small files a genome), unless --done_folders. Safe beside other builds: a pack and a change
    of one genome wait for each other, the readers take both forms. In a random order, so that several builds'
    passes (a SLURM array's tasks) seldom meet at one genome."""
    if opts.done_folders:
        return
    folders = store.unpacked(os.path.join(opts.work, "done"))
    if not folders:
        return
    say(f"done folders: {len(folders):,} to pack, each into one file (WORK/done/<accession>.zip)")
    random.shuffle(folders)
    failed = []

    def pack(folder):
        try:
            store.pack(folder)
        except OSError as e:
            failed.append((os.path.basename(folder), str(e)[:100]))
    with concurrent.futures.ThreadPoolExecutor(max(1, min(PACK_THREADS, opts.threads))) as pool:
        list(pool.map(pack, folders))
    store.close()
    say(f"  {len(folders) - len(failed):,} packed"
        + (f", {len(failed):,} failed (left as folders, read as before; e.g. {failed[:2]})" if failed else ""))


def drop_superseded(opts, job_map):
    """The genome files another assembly took the place of, deleted once it is done (genomes.drop_superseded); kept
    without --keep_gb."""
    if opts.keep_gb is None:
        return
    gone = genomes.drop_superseded(opts.genomes, opts.work, job_map)
    if gone:
        say(f"{len(gone)} superseded genomes deleted, their replacements done (e.g. {gone[:3]})")


def run_one(item):
    """genomes.run in a worker, any error the genome's failure (until 2026-11-06 an error genomes.run did not expect,
    such as a FASTA header not UTF-8, stopped the pool and every other genome, on every rerun)."""
    try:
        return genomes.run(item)
    except Exception as e:  # noqa: BLE001
        return item[0]["accession"], None, f"{type(e).__name__}: {e}"


def genome_size(opts, job):
    """A genome's bases: its plan's size, else its file's (TOGA2's have none in the plan), else 0."""
    path = job.get("path") or genomes.genome_file(opts.genomes, job["accession"])
    return job.get("size") or (os.path.getsize(path) / placement.GZIP_BYTES_PER_BASE if os.path.isfile(path) else 0)


def genomes_to_run(opts, job_map, settings):
    """The jobs of job_map to run, each with its settings -> ([(job, settings)], how many again for gene models)."""
    # a done genome of a BUSCO dataset species without its held-out calls (run before 2026-11-06) is run again where
    # its genome is still there: never fetched again for a benchmark, its OrthoDB benchmark then left without them
    # and a done genome with better gene models than its run took (genomes.better_annotation: an index naming its
    # assembly since, a VEuPathDB key given since, its annotation fetched since), its genome and that annotation here
    # (since 2026-11-07: until then it kept miniprot's models, the annotation compared alone)
    # and (since 2026-11-08) a done genome run without the rescue pass where this run makes it (genomes.rescue_missing)
    def here(a, j):
        return os.path.isfile(j.get("path") or genomes.genome_file(opts.genomes, a))
    # each genome's done.json read once for all these questions (genomes.cached_summaries; since 2026-11-10)
    with genomes.cached_summaries():
        better = {a for a, j in job_map.items() if genomes.is_done(opts.work, a, settings["loci"])
                  and genomes.rerun_for_annotation(opts.work, j, opts.genomes, not settings["staged_only"])}
        # and (since 2026-11-09) a done genome of an older genome record (genomes.record_stale; --reversion stamps
        # them)
        whole = {a for a, j in job_map.items() if not genomes.is_done(opts.work, a, settings["loci"])
                 or ((genomes.held_out_missing(opts.work, j, settings)
                      or genomes.rescue_missing(opts.work, a, settings)
                      or genomes.record_stale(opts.work, a)) and here(a, j))}
        # each with its done.json as seen now: a genome another build on the WORK finishes meanwhile is left alone
        # (genomes.done_stamp; since 2026-11-08 it was run again whole)
        return [(dict(j, done_stamp=genomes.done_stamp(opts.work, a)) if a in whole
                 else dict(j, rerun="annotation", done_stamp=genomes.done_stamp(opts.work, a)), settings)
                for a, j in sorted(job_map.items()) if a in whole or a in better], len(better - whole)


def panel_settings(settings):
    """The intraspecific panel's further assemblies' settings: their BUSCOs and rRNA alone (no tiara, no held-out)."""
    return dict(settings, tiara=None, odb_own={})


def run_genomes(opts, job_map, settings, extras=None):
    """Every genome of job_map to run, and of extras (the intraspecific panel's further assemblies, panel_settings: in
    the same pool since 2026-11-08, until then a pool of their own after it, each pool's last genomes leaving the
    cores idle), run by schedule -> job_map's genomes (waiting: not fetched yet, failed: [(accession, why)]).
    A run that makes the passes after the genomes (no --genomes_only) waits for the genomes another build holds
    (others_end) and runs those that build left (since 2026-11-10: it went on without them, so that the last job
    after a SLURM array made its passes, and the package, without the genomes of a task still running or killed
    within genomes.CLAIM_STALE)."""
    extras = extras or {}
    todo, again = genomes_to_run(opts, job_map, settings)
    panel = genomes_to_run(opts, extras, panel_settings(settings))[0] if extras else []
    say(f"genomes to run: {len(todo)} of {len(job_map)}"
        + (f", {again} of them again for better gene models" if again else "")
        + (f", and {len(panel)} of the intraspecific panel's further assemblies" if panel else "")
        + ("" if settings["tiara"] else " (no tiara: windows not classified)"))
    waiting, failed = [], []
    # killed runs' partial folders and (a packed genome's half-written file; since 2026-11-10) files; not those of
    # another build's genomes going on (their claims live)
    for leftover in genomes.stale_partials(opts.work) + (genomes.stale_scratch(opts.work, opts.scratch)
                                                          if opts.scratch else []):
        if os.path.isdir(leftover):
            shutil.rmtree(leftover, ignore_errors=True)
        else:
            with contextlib.suppress(OSError):
                os.remove(leftover)
    todo += panel
    while todo:
        # largest first (since 2026-11-08; until then by accession, a 30 Gb genome begun last kept the other workers
        # idle for hours at the end)
        sizes = {j["accession"]: genome_size(opts, j) for j, _s in todo}
        todo = sorted(todo, key=lambda it: (-sizes[it[0]["accession"]], it[0]["accession"]))
        # every --report_minutes, whether a genome finished or not: genomes left and the time left, by their sizes (a
        # genome not fetched yet is no work), and the cores, memory and disk used (progress_of)
        progress = progress_of(opts, "genomes run", [sizes[j["accession"]] for j, _s in todo])
        held, elsewhere = [], 0

        def note():
            return (f", {len(waiting):,} not fetched yet" if waiting else "") \
                + (f", {elsewhere:,} run by another build" if elsewhere else "")
        for acc, _summary, why in schedule(opts, todo, sizes, progress, note):
            skipped = why in ("not fetched yet", genomes.RUNNING_ELSEWHERE, genomes.DONE_ELSEWHERE)
            if why in (genomes.RUNNING_ELSEWHERE, genomes.DONE_ELSEWHERE):
                elsewhere += 1
                if why == genomes.RUNNING_ELSEWHERE:
                    held.append(acc)
            elif why and acc in job_map:  # the panel's: a later run tries again, the package goes on without
                (waiting if why == "not fetched yet" else failed).append((acc, why))
            progress.step(sizes[acc], failed=bool(why) and not skipped, skipped=skipped, note=note())
        if opts.genomes_only or not held:  # --genomes_only: those genomes are the other builds', as are the passes
            break
        others_end(opts, held)
        # what the other builds left of them (a task killed at its time limit, a genome that failed there): run here
        todo = genomes_to_run(opts, {a: job_map[a] for a in held if a in job_map}, settings)[0] \
            + genomes_to_run(opts, {a: extras[a] for a in held if a in extras}, panel_settings(settings))[0]
        say(f"  the other builds' {len(held):,} genomes ended: {len(todo):,} of them left to run here")
    return waiting, failed


def others_end(opts, held, poll=60):
    """Waits until no other build holds a genome of held (accessions whose run met another build's claim): its claim
    gone (the genome done, or failed there), or no longer live (genomes.claim_live: a killed job's, its heartbeat
    stopped, after genomes.CLAIM_STALE)."""
    say(f"{len(held):,} genomes are run by another build (e.g. {sorted(held)[:3]}): waiting for them before the "
        "passes after the genomes")
    began = said = time.time()
    while True:
        live = [a for a in held if genomes.claim_live(genomes.claim_path(opts.work, a))]
        if not live:
            return
        if time.time() - said >= opts.report_minutes * 60:
            said = time.time()
            say(f"  {len(live):,} genomes still run by another build after {util.duration(said - began)} (e.g. "
                f"{sorted(live)[:3]})")
        time.sleep(poll)


MINIPROT_BYTES_PER_BASE = 7.5  # miniprot's index: ~22 GB for a 3 Gb genome (in placement.CHUNK_BASES pieces at most)
RUN_MEMORY = 3e9               # a genome run's other memory: its worker, tiara's model, the rRNA and metadata scans
MAX_GENOME_THREADS = 32        # a genome's threads at most (miniprot's index is made on one)
MEMORY_SHARE = 0.85            # of the memory found (SLURM's, the cgroup's, the node's): the estimates' room


def run_memory(bases):
    """A genome run's memory (bytes), estimated from its bases."""
    return MINIPROT_BYTES_PER_BASE * min(bases, placement.CHUNK_BASES) + RUN_MEMORY


def genome_threads(opts, memory):
    """A genome's threads: -t / -j, more for a genome whose memory takes a larger share of --memory_gb than -t / -j
    of -t (so that the cores the memory leaves are its: a node of 120 cores and 200 GB runs ~6 genomes of 3 Gb at
    once, ~20 threads each, not 30 at 4 each, ~700 GB)."""
    base = opts.threads_per_genome
    if not opts.memory:
        return base
    return max(base, min(MAX_GENOME_THREADS, opts.threads, round(opts.threads * memory / opts.memory)))


def schedule(opts, todo, sizes, progress, note, poll=60, task=None):
    """genomes.run (run_one; or task, a pool's function of (job, settings) that gives such a triple: backfill_one)
    of each (job, settings) in todo (largest first), at most -j at once, each begun when its memory
    (run_memory) fits in what --memory_gb leaves and -t / -j threads are free (since 2026-11-08: -j at once whatever
    their sizes, a node's memory over-run by its largest genomes) -> (accession, summary, why) as they finish. A
    genome whose memory does not fit lets a smaller one that does go first; one larger than --memory_gb runs alone.
    Its threads: genome_threads, at most those free. Progress is said on time between them (progress.tick, note()).
    A worker that dies (the OOM killer) stops the build (util.process_pool; until 2026-11-08 it waited forever)."""
    budget, cores, base = opts.memory, opts.threads, min(opts.threads_per_genome, opts.threads)
    # the genomes to begin, largest first, and their memory as negatives (ascending: bisect finds the first that fits
    # in what is free, not a scan of some 30,000 at each genome's end)
    pending = sorted(todo, key=lambda it: -run_memory(sizes[it[0]["accession"]]))
    negative = [-run_memory(sizes[it[0]["accession"]]) for it in pending]
    running, said, futures = {}, set(), {}

    def begin(item, memory, threads):
        job, settings = item
        settings = dict(settings, threads=threads,
                        rrna=functools.partial(settings["rrna"], threads=threads) if settings.get("rrna") else None)
        running[job["accession"]] = (memory, threads)
        futures[pool.submit(task or run_one, (job, settings))] = job["accession"]

    with util.process_pool(opts.jobs) as pool:
        while pending or running:
            free_threads = cores - sum(t for _m, t in running.values())
            free_memory = budget - sum(m for m, _t in running.values()) if budget else None
            while pending and (not running or (len(running) < opts.jobs and free_threads >= base)):
                # the largest that fits (the first, where nothing runs or no memory is known)
                i = bisect.bisect_left(negative, -free_memory) if running and budget else 0
                if i == len(pending):
                    break
                item, memory = pending.pop(i), -negative.pop(i)
                threads = max(1, min(genome_threads(opts, memory), free_threads if running else cores))
                if budget and memory > budget and "over" not in said:
                    said.add("over")
                    say(f"WARNING: {item[0]['accession']}'s run takes ~{util.size_text(memory)}, more than "
                        f"--memory_gb {budget / 1e9:,.0f}: run alone")
                begin(item, memory, threads)
                free_threads -= threads
                free_memory = free_memory - memory if budget else None
            done, _left = concurrent.futures.wait(futures, timeout=poll,
                                                  return_when=concurrent.futures.FIRST_COMPLETED)
            if not done:
                busy = sum(t for _m, t in running.values())
                progress.tick(note() + f", {len(running)} running on {busy} threads"
                              + (f", ~{util.size_text(sum(m for m, _t in running.values()))} of "
                                 f"{util.size_text(budget)} memory estimated" if budget else ""))
                continue
            for future in done:
                acc = futures.pop(future)
                running.pop(acc, None)
                try:
                    got = future.result()
                except concurrent.futures.BrokenExecutor:
                    raise
                except Exception as e:  # noqa: BLE001  the task not run (e.g. not pickled): the genome's failure
                    got = (acc, None, f"{type(e).__name__}: {e}")
                yield got


def progress_of(opts, what, weights, unit="Gb", scale=1e9):
    """A pass's progress (util.Progress): said every minute of its first 10, then every --report_minutes, each
    report with the cores, memory and disk the build uses (opts.usage; since 2026-11-08)."""
    every = opts.report_minutes * 60
    return util.Progress(what, weights, say, every, unit, scale, early=util.early_reports(every),
                         usage=opts.usage.line)


def running(opts, progress, workers=None):
    """', 4 running': a pool's items under way (at most its workers, default -j)."""
    return f", {min(max(1, workers or opts.jobs), progress.n - progress.done):,} running"


def backfill_pass(opts, job_map, settings):
    """The done genomes' steps a tool missing at their run left out (genomes.backfill: tiara's classes, the rRNA,
    a dataset species' clade hold-out), made now the tool is here (since 2026-11-08: never). Run by schedule (since
    2026-11-10: a pool of -j whatever the genomes' sizes): a clade hold-out is a call of the whole genome, with a
    genome run's memory, the other steps take a worker's alone."""
    with genomes.cached_summaries():
        needs = {a: genomes.backfill_needs(j, settings) for a, j in sorted(job_map.items())
                 if genomes.is_done(opts.work, a)}
    todo = [(job_map[a], settings) for a, steps in needs.items() if steps]
    if not todo:
        return
    say(f"backfill: {len(todo)} done genomes lack tiara's classes, their rRNA or their clade hold-out, made now")
    made, failed = collections.Counter(), []
    # schedule's sizes: the genome's bases where it is called again, 0 (a run's other memory alone) where not
    sizes = {j["accession"]: genome_size(opts, j) if "held_out" in needs[j["accession"]] else 0 for j, _s in todo}
    progress = progress_of(opts, "genomes backfilled", [1] * len(todo), "", 1)
    for acc, steps, why in schedule(opts, todo, sizes, progress, lambda: "", task=backfill_one):
        made.update(steps or [])
        if why:
            failed.append((acc, why))
        progress.step(1, failed=bool(why))
    say(f"  backfill: {dict(made)} made, {len(failed)} failed{' (e.g. ' + str(failed[:2]) + ')' if failed else ''}")


def backfill_one(item):
    try:
        return genomes.backfill(item)
    except Exception as e:  # noqa: BLE001
        return item[0]["accession"], [], f"{type(e).__name__}: {e}"


def published(opts, src_dir, ncbi):
    """The published reads (EUKARYOME, ROD) with their regions and NCBI labels -> records [(name, sequence, row)],
    their regions kept in WORK/published/PUBLISHED_REGIONS (a rerun reads them)."""
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
    # the regions cache: one file per process that annotated (since 2026-11-08: one file all runs appended to, the
    # lines of two package stages on one WORK interleaved), all read
    stem = os.path.join(cache_dir, PUBLISHED_REGIONS[:-len(".tsv")])
    known = {}
    for path in [stem + ".tsv"] + sorted(glob.glob(glob.escape(stem) + ".*.tsv")):
        known.update(read_regions(path))
    cache = f"{stem}.{socket.gethostname()}-{os.getpid()}.tsv"
    # a read is known where its line names it and (a line since 2026-11-08) its sequence's digest: a source's new
    # release with a read changed under its name has it annotated again
    known = {n: known[n] for n, s, *_ in reads if n in known and known[n][3] in (None, read_digest(s))}
    unread = set()  # the reads of a batch that failed (C13)
    # by ITSx's profile sets, so that a batch's reads share theirs (an ITSx run per distinct set)
    todo = sorted(((n, s, rrna.itsx_groups(ncbi, t)) for n, s, t, *_ in reads if n not in known),
                  key=lambda r: (r[2], r[0]))
    if todo and shutil.which("barrnap"):
        say(f"published reads: {len(todo)} to annotate (barrnap, ITSx)")
        # small batches over every thread: barrnap (nhmmer) runs one thread a batch and dominates (~0.5 s a read)
        batches = [todo[i:i + PUBLISHED_BATCH] for i in range(0, len(todo), PUBLISHED_BATCH)]
        work = functools.partial(_annotate_batch, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"))
        failed = []
        unread = {n for n, _s, _t in todo}
        with util.process_pool(opts.threads) as pool, open(cache, "a", newline="\n") as fh:
            if fh.tell() and not util.ends_in_newline(cache):  # a killed run's last line cut: the next starts anew
                fh.write("\n")
            digest = {n: read_digest(s) for n, s, _t in todo}
            for got, why in util.unordered(pool, work, batches):
                if why:
                    failed.append(why)
                fh.writelines(f"{name}\t{regions}\t{complete}\t{strand}\t{digest[name]}\n"
                              for name, (regions, complete, strand) in got.items())
                fh.flush()
                known.update({n: v + (digest[n],) for n, v in got.items()})
        unread -= set(known)
        if failed:
            # left out (since 2026-11-10: each went into the package as an 18S record of its whole length, a 4-6 kb
            # operon read among the 18S, on the strand it came in, until a rerun annotated it)
            say(f"WARNING: published reads: {len(failed)} batches of {PUBLISHED_BATCH} failed (e.g. {failed[0]}): "
                f"their {len(unread):,} reads left out of this package, annotated again by a rerun")
    for name, seq, taxid, source, how, accession in reads:
        if name in unread:
            continue
        regions, complete, strand, _digest = known.get(name, ("-", "f", "+", None))
        regions_list = rrna.parse_regions(regions)
        seq = rrna.revcomp(seq) if strand == "-" else seq  # the record on its unit's strand, as its regions
        out.append((name, seq, dict(kind="rDNA" if regions_list else "18S", accession=accession, source=source,
                                    taxid=taxid, species_taxid=ncbi.species_of(taxid), regions_list=regions_list,
                                    lsu_complete=complete == "t", label_from=how)))
    return out


def _annotate_batch(batch, barrnap=None, itsx=None):
    """[(name, sequence, ITSx's sets)] -> ({name: (regions text, lsu_complete t/f, strand)}, None), or ({}, why) where
    barrnap or ITSx failed on it (since 2026-11-08: the error stopped the package stage)."""
    try:
        found = rrna.published_regions([(n, s) for n, s, _t in batch], barrnap, itsx,
                                       sets={n: t for n, _s, t in batch})
    except Exception as e:  # noqa: BLE001  a batch's failure, never the stage's
        return {}, f"{type(e).__name__}: {e}"[:200]
    return {name: (rrna.regions_text(found[name][0]), "t" if found[name][1] else "f", found[name][3])
            if name in found else ("-", "f", "+") for name, _s, _t in batch}, None


def read_digest(seq):
    """A published read's digest in the regions cache: its sequence's MD5, 12 hex digits."""
    return hashlib.md5(seq.encode()).hexdigest()[:12]


def read_regions(path):
    """The published reads' regions cache -> {name: (regions text, lsu_complete t/f, strand, its read's digest, or None
    in a line before 2026-11-08)}; a line a killed run cut (not 4 or 5 fields, or no line end) left out, so that its
    read is annotated again (until 2026-11-06 it stopped every later run)."""
    known = {}
    if os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if line.endswith("\n") and len(f) in (4, 5) and f[3] in ("+", "-"):
                    known[f[0]] = (f[1], f[2], f[3], f[4] if len(f) == 5 else None)
    return known


def genome_rrna_records(opts, job_map, ncbi):
    """Every done genome's rRNA as records, labelled by its assembly's taxid."""
    out = []
    for acc, job in sorted(job_map.items()):
        rows = rrna.load(genomes.done_dir(opts.work, acc))
        if rows:
            out += rrna.genome_records(acc, rows, job.get("taxid") or job["species_taxid"],
                                       ncbi.species_of(int(job["species_taxid"])))
    return out


def annotate_one(work, genomes_folder, may_fetch, veupathdb_key, job):
    """annotation.annotate_genome in a worker -> (accession, copies compared or None, why it failed or None)."""
    try:
        return job["accession"], annotation.annotate_genome(work, job, genomes_folder, may_fetch, veupathdb_key), None
    except OSError as e:
        return job["accession"], None, str(e)[:100]
    except Exception as e:  # noqa: BLE001  one genome's failure, never the pass' (until 2026-11-08: OSError alone)
        return job["accession"], None, f"{type(e).__name__}: {e}"[:100]


def annotate_pass(opts, job_map, may_fetch):
    """The done genomes with an annotation (NCBI's, Ensembl's or TOGA2's) compared with it (annotation.py: its GFF
    beside the genome, from the fetch stage, or fetched here where may_fetch), the GFF deleted after."""
    todo = [j for _a, j in sorted(job_map.items()) if annotation.missing(opts.work, j)]
    if not todo:
        return
    say(f"gene models: {len(todo)} genomes' calls to compare with their annotation")
    compared, left, failed = 0, 0, []
    progress = progress_of(opts, "genomes compared", [1] * len(todo), "", 1)
    # -j at once, a GFF read on one core each (until 2026-11-08 one at a time, in this process)
    task = functools.partial(annotate_one, opts.work, opts.genomes, may_fetch, opts.veupathdb_key)
    with util.process_pool(opts.jobs) as pool:
        for acc, n, why in util.unordered(pool, task, todo, progress, lambda: running(opts, progress)):
            if why:
                failed.append((acc, why))
                progress.step(1, failed=True)
                continue
            compared, left = compared + (n is not None), left + (n is None)
            progress.step(1)
    say(f"  {compared} genomes compared, {left} without their annotation fetched yet (--stage fetch), "
        f"{len(failed)} failed{' (e.g. ' + str(failed[:2]) + ')' if failed else ''}")


def saved_molecules(opts, tmp, accession):
    """The organelle molecules a done genome's metadata scan kept (metadata.ORGANELLE, written with its stats.json),
    or None: the organelle pass reads them, not the genome again (since 2026-11-08). A packed genome's (store.py)
    written into tmp for the pass."""
    folder = genomes.done_dir(opts.work, accession)
    if not store.has(folder, metadata_mod.STATS):
        return None
    if not store.is_packed(folder):
        path = os.path.join(folder, metadata_mod.ORGANELLE)
        return path if os.path.isfile(path) else None
    path = os.path.join(tmp, f"{accession}.{metadata_mod.ORGANELLE}")
    return path if os.path.isfile(path) or store.extract(folder, metadata_mod.ORGANELLE, path) else None


def organelle_genes(opts, src_dir, job_map, ncbi):
    """The organelle genes' sources (eukore/organelle.py): RefSeq's organelle genomes, parsed once, and each species'
    own assembly where RefSeq has neither of its organelles' records (barrnap, nhmmer) -> RefSeq's records chosen,
    {(species taxid, organelle): rows}."""
    folder = os.path.join(opts.work, "organelle")
    chosen = organelle_mod.choose(organelle_mod.refseq_genes(src_dir, folder, ncbi, say), ncbi)
    if not chosen:
        say("WARNING: no RefSeq organelle release (--stage fetch fetches it): organelle genes of the assemblies "
            "alone, and no COX1 among them")
    barrnap, nhmmer = shutil.which("barrnap"), shutil.which("nhmmer")
    if barrnap and nhmmer:
        # an assembly a task on one core (its organelle molecules, small: barrnap, nhmmer and miniprot on them), -t at
        # once (until 2026-11-08 -j at once, -t / -j threads each, mostly idle)
        with tempfile.TemporaryDirectory() as tmp:  # the packed genomes' molecules, as files for the tools
            organelle_mod.assembly_pass(folder, job_map, lambda a: genomes.genome_file(opts.genomes, a), chosen, ncbi,
                                        barrnap, nhmmer, opts.threads, 1, say,
                                        shutil.which(opts.miniprot) if opts.organelle_models == "miniprot" else None,
                                        saved_of=functools.partial(saved_molecules, opts, tmp),
                                        progress=functools.partial(progress_of, opts, unit="", scale=1))
    else:
        say("WARNING: no barrnap or nhmmer: the assemblies' organelle molecules not searched")
    return chosen


def extras_of(opts, src_dir, ncbi, job_map, assemblies):
    """The intraspecific panel's further assemblies (genomes.intraspecific_jobs; since 2026-11-07) -> {accession: job},
    {} with --intraspecific 0."""
    if not opts.intraspecific:
        return {}
    refseq_rows = sources.read_assemblies(tables(src_dir, "refseq", opts.groups))
    extras = genomes.intraspecific_jobs(ncbi, assemblies, refseq_rows, job_map, genomes.read_panel(),
                                        opts.intraspecific, opts.inside)
    if extras:
        say(f"intraspecific panel: {len(extras)} further assemblies of "
            f"{len({j['species_taxid'] for j in extras.values()})} species (--intraspecific {opts.intraspecific})")
    return extras


def intraspecific_pass(opts, src_dir, ncbi, extras):
    """The panel's further assemblies, run (their BUSCOs and rRNA; no tiara, no held-out calls: panel_settings) with
    the plan's genomes by run_genomes: their organelle molecules searched for every organelle gene, then their genome
    files deleted: only their results are kept."""
    if not extras:
        return
    folder = os.path.join(opts.work, "organelle")
    chosen = organelle_mod.choose(organelle_mod.refseq_genes(src_dir, folder, ncbi, say), ncbi)
    barrnap, nhmmer = shutil.which("barrnap"), shutil.which("nhmmer")
    if barrnap and nhmmer:
        organelle_mod.assembly_pass(folder, extras, lambda a: genomes.genome_file(opts.genomes, a), chosen, ncbi,
                                    barrnap, nhmmer, opts.threads, 1, say,
                                    shutil.which(opts.miniprot) if opts.organelle_models == "miniprot" else None,
                                    every=True, saved_of=functools.partial(saved_molecules, opts),
                                    progress=functools.partial(progress_of, opts, unit="", scale=1))
    for acc in extras if opts.keep_gb is not None else ():  # without --keep_gb, kept
        path = genomes.genome_file(opts.genomes, acc)
        if genomes.is_done(opts.work, acc) and os.path.isfile(organelle_mod.assembly_file(folder, acc)) \
                and os.path.isfile(path):
            os.remove(path)


def intraspecific_writer(opts, ncbi, extras, collected, records, organelle, vsearch):
    """The barcode gap (intraspecific.py) -> a package writer of intraspecific/, the extras' genomes.tsv rows set."""
    rows, markers, whole, rrna_out, org_out = intraspecific_mod.collect(opts.work, extras,
                                                                        os.path.join(opts.work, "organelle"))
    reps = {a: int(g["species_taxid"]) for a, g in collected["genomes"].items() if g["representative"] == "t"}
    panel = {sp: a for a, sp in reps.items() if any(int(j["species_taxid"]) == sp for j in extras.values())}
    rep_markers = collections.defaultdict(dict)
    for (acc, locus), seq in collected.get("cds", {}).items():
        if acc in panel.values():
            rep_markers[acc][locus] = seq
    by_acc = collections.defaultdict(list)
    for _name, seq, row in records:
        by_acc[row.get("accession")].append(dict(kind=row["kind"], sequence=seq, regions=row.get("regions_list"),
                                                 lsu_complete=row.get("lsu_complete") in (True, "t")))
    for acc in panel.values():
        rep_markers[acc].update(intraspecific_mod.first_rrna(by_acc.get(acc, [])))
    org_by_sp = collections.defaultdict(list)
    for _name, seq, row in organelle:
        org_by_sp[int(row["species_taxid"])].append(dict(row, sequence=seq))
    for sp, acc in panel.items():
        rep_markers[acc].update(intraspecific_mod.organelle_seqs(org_by_sp.get(sp, [])))
    # the other species' copies the representatives' nearest is measured against (the BUSCOs': qc.check's already)
    others = collections.defaultdict(dict)
    for acc, sp in reps.items():
        for marker, seq in intraspecific_mod.first_rrna(by_acc.get(acc, [])).items():
            others[marker][f"{sp}|{acc}"] = seq
    for sp, rows_ in org_by_sp.items():
        for marker, seq in intraspecific_mod.organelle_seqs(rows_).items():
            others[marker][f"{sp}|org"] = seq
    # its nearest other species' identity and shared rank (qc.check's; the rank since 2026-11-09: resolved against a
    # congener alone)
    inter = {(reps[r["accession"]], r["locus"]): None if r.get("nearest_id") is None else
             (r["nearest_id"], r.get("nearest_depth")) for r in collected["copy_rows"]
             if r["accession"] in panel.values()}
    species_of = dict(reps, **{a: int(j["species_taxid"]) for a, j in extras.items()})
    gap = intraspecific_mod.barcode_gap(panel, rep_markers, markers, species_of, inter, others, ncbi, vsearch,
                                        opts.threads, {sp: groups_mod.group_of(ncbi, sp) for sp in panel})
    for acc, job in extras.items():
        if acc in whole:
            collected["genomes"].setdefault(acc, dict(accession=acc, species_taxid=job["species_taxid"],
                                                      taxid=job["taxid"], organism=job.get("organism", "")))
            collected["genomes"][acc].update(source="intraspecific", representative="f", loci=len(whole[acc]),
                                             fetched=job.get("fetch", acc), gencode=job.get("gencode", 1))
    say(f"intraspecific: {len(whole)} further assemblies of {len(panel)} species, {len(gap)} barcode gaps")
    return functools.partial(intraspecific_mod.write, rows=rows, whole=whole, rrna_out=rrna_out, org_out=org_out,
                             gap_rows=gap)


@contextlib.contextmanager
def timed(stages, step):
    """A step of the package stage, its seconds added to stages ([(step, seconds)]: benchmarks' stages.tsv)."""
    started = time.time()
    yield
    stages.append((step, time.time() - started))


def reference_screen(opts, src_dir, ncbi):
    """The reference screen of the user's MIDORI2 files (eukore/refscreen.py; since 2026-11-10): their records that
    are another organism's gene, a chimera of two organisms' or another gene, told by their windows against RefSeq's
    organelle genes and the other species' records -> (calls, counts), or None where it is not run (no --midori2,
    --no_reference_screen, no minimap2 or vsearch). Kept in WORK/reference_screen/<key>/ (calls.tsv, KSGP's lists
    drop.tsv, trim.tsv and review.tsv, counts.json), the key the files', RefSeq's release's, the taxdump's and the
    screen's version: a rerun reads it."""
    if not opts.midori2 or opts.no_reference_screen:
        return None
    minimap2, vsearch = shutil.which(opts.minimap2), shutil.which("vsearch")
    if not minimap2 or not vsearch:
        say(f"WARNING: no {opts.minimap2 if not minimap2 else 'vsearch'}: MIDORI2's records not screened "
            "(eukore/refscreen.py), every one names its species")
        return None
    release = os.path.join(organelle_mod.folder_of(src_dir), "release.txt")
    stamp = [refscreen.VERSION, util.md5(os.path.join(src_dir, "taxdump.tar.gz")),
             util.md5(release) if os.path.isfile(release) else "-"] + [util.md5(p) for p in opts.midori2]
    folder = os.path.join(opts.work, "reference_screen", hashlib.md5(",".join(stamp).encode()).hexdigest()[:10])
    kept = os.path.join(folder, "counts.json")
    if os.path.isfile(kept):
        with open(kept) as fh:
            counts = json.load(fh)
        calls = refscreen.read_calls(os.path.join(folder, "calls.tsv"))
        say(f"reference screen: {folder} is there ({counts.get('action_drop', 0):,} records dropped)")
        return calls, counts
    rows = organelle_mod.refseq_genes(src_dir, os.path.join(opts.work, "organelle"), ncbi, say)
    if not rows:
        say("WARNING: no RefSeq organelle release (--stage fetch fetches it): MIDORI2's records are screened against "
            "one another alone")
    os.makedirs(folder, exist_ok=True)
    calls, counts = refscreen.screen([refscreen.midori2_source(p) for p in opts.midori2],
                                     list(refscreen.trusted_records(rows)), ncbi, minimap2, vsearch, opts.threads,
                                     tmp=folder, log=say)
    refscreen.write_lists(folder, calls)
    for line in refscreen.summary(calls, counts, opts.midori2).splitlines():
        say(line)
    with open(util.part_of(kept), "w", newline="\n") as fh:
        json.dump(dict(counts), fh, indent=1)
    os.replace(util.part_of(kept), kept)
    return calls, dict(counts)


def databases_of(src_dir, ncbi, references, inside=None):
    """{name: NCBI species} of the other references the build read, for overlap.tsv: PR2, EUKARYOME's long reads and
    KSGP (labels.references, by their species names) and ROD (by its taxids)."""
    out = benchmarks.reference_species(references, ncbi, inside)
    rod = sources.ready(src_dir, "ROD")
    if rod:
        taxa = {ncbi.species_of(int(r["taxid"])) for r in sources.rod_operons(rod) if r.get("taxid", "").isdigit()}
        out["ROD"] = {sp for sp in taxa if sp and ncbi.rank.get(sp) == "species" and (inside is None or inside(sp))}
    return out


def measure(part, opts, src_dir, ncbi, collected, vsearch, stages, fates=None, tests=(), references=()):
    """The package's benchmarks (eukore/benchmarks.py) into its unfinished folder, and (since 2026-11-10, without
    --no_figures) their figures (eukore/figures.py) -> summary.json's figures, or None where they failed: a benchmark
    never stops the package. fates: genomes.fates (selection.tsv); tests: the read benchmark's genomes
    (reads.choose), kept through the prune; references: labels.references (overlap.tsv)."""
    nearest = {(r["accession"], r["locus"]): (r.get("nearest_id"), r.get("nearest_depth"))
               for r in collected["copy_rows"]} if vsearch else None
    minimap2 = shutil.which(opts.minimap2)
    if tests and not minimap2:
        say(f"WARNING: {opts.minimap2} not found: no read benchmark (bioconda's minimap2)")
    # its files in a folder of this run's (since 2026-11-08: WORK/read_benchmark, which two package stages on one WORK
    # wrote and deleted under each other); a killed run's removed by the next (util.dead_partials). Since 2026-11-10
    # in --local_tmp where given: the whole package's FASTAs and minimap2's index, written and read once each, until
    # then on WORK's shared file system
    reads = dict(genomes_folder=opts.genomes, tests=list(tests), minimap2=minimap2,
                 background=reads_mod.read_background(READ_BACKGROUND, os.path.join(src_dir, "read_background")),
                 work=util.part_of(os.path.join(getattr(opts, "scratch", None) or opts.work, "read_benchmark")))
    inside = opts.inside if opts.wanted is not None else None
    if not opts.read_compare:
        say("NOTE: no --read_compare: the read benchmark and breadth.tsv compare Eukore with no other reference "
            "(tools/compare_references.py makes them; tools/eukore.sh takes them from its COMPARE folder)")
    try:
        measured = benchmarks.run(part, ncbi, odb_mod.folder_of(src_dir), vsearch, opts.threads, opts.work, nearest,
                                  stages, log=say, fates=fates, compare=opts.read_compare, inside=inside, reads=reads,
                                  databases=databases_of(src_dir, ncbi, references, inside),
                                  marker_folder=markers_of(opts, src_dir),
                                  clades=len(opts.wanted) if opts.wanted is not None else None)
    except Exception as e:  # noqa: BLE001
        say(f"WARNING: the benchmarks failed ({type(e).__name__}: {e}): the package is written without them")
        shutil.rmtree(os.path.join(part, benchmarks.FOLDER), ignore_errors=True)
        return None
    if not opts.no_figures:  # since 2026-11-10: Figure 1 and the rest, before the manifest lists them
        draw_figures(part)
    return measured


def draw_figures(part):
    """The package's figures (eukore/figures.py) into its benchmarks/figures/; without matplotlib, or failed, left
    out: never the package's end."""
    started = time.time()
    try:
        from eukore import figures  # noqa: PLC0415 (matplotlib: only where figures are drawn)
    except ImportError as e:
        say(f"WARNING: no figures ({e}: envs/eukore.yaml's matplotlib-base): python3 -m eukore.figures draws them "
            "later")
        return
    try:
        made = figures.draw(part, log=say)
    except Exception as e:  # noqa: BLE001
        say(f"WARNING: the figures failed ({type(e).__name__}: {e}): the package is written without them")
        shutil.rmtree(os.path.join(part, benchmarks.FOLDER, figures.FOLDER), ignore_errors=True)
        return
    drawn = [n for n, p in made.items() if p and n not in ("index", "summary")]
    say(f"figures: {len(drawn)} drawn into {benchmarks.FOLDER}/{figures.FOLDER}/ ({time.time() - started:.0f} s); "
        f"left out: {', '.join(n for n, p in made.items() if not p) or 'none'}")


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


def at_a_time(work, items, n):
    """work(item) of each item, n at a time in threads (downloads: urllib and zlib let go of the GIL) -> their futures as
    they finish. A stop (Ctrl-C) drops those not started instead of waiting for them all."""
    pool = concurrent.futures.ThreadPoolExecutor(max(1, n))
    util.STOP.clear()
    try:
        yield from concurrent.futures.as_completed([pool.submit(work, item) for item in items])
    except BaseException:
        util.STOP.set()  # the running downloads end at their next try, not begun again
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)


def scratch_problem(path):
    """Why --local_tmp's folder can't be used (made where not there, then written in), or None. A SLURM job's
    $SLURM_LOCAL_SCRATCH may be named where the job asked for no local disk, its folder then not made."""
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        return f"{path}: {e.strerror} (a job's local scratch is there only where the job asks for it: --localscratch)"
    if not os.access(path, os.W_OK | os.X_OK):
        return f"{path}: not writable"
    return None


def command_of(argv):
    """The command line as a shell takes it, a VEuPathDB key on it hidden (a log is read by others)."""
    words, hide = [os.path.abspath(__file__)], False
    for w in argv:
        flag, eq, _key = w.partition("=")
        if hide:
            words.append("HIDDEN")
            hide = False
            continue
        if len(flag) > 4 and "--veupathdb_key".startswith(flag):  # argparse takes a flag's unique prefix too
            if eq:
                words.append(flag + "=HIDDEN")
                continue
            hide = True
        words.append(shlex.quote(w))
    return " ".join(words)


def separate_genomes(opts):
    """Whether --genomes puts the genomes elsewhere than OUT/genomes."""
    return os.path.abspath(opts.genomes) != os.path.abspath(os.path.join(opts.out, "genomes"))


def banner(opts, argv, whence):
    """The run's opening lines on STDOUT (since 2026-11-08): Eukore's version, where and when, the command, what each
    stage of the run will do and with what (the folders, the genome set, the cores, memory and downloads)."""
    stages = opts.stages
    package_dir = f"Eukore_{opts.version}"
    does = {
        "fetch": "the sources" + (" (their tables fetched again: --refresh)" if opts.refresh else "")
                 + ", then each species' genome and gene models into "
                 + ("--genomes" if separate_genomes(opts) else "OUT/genomes"),
        "process": "each genome's loci called and its rRNA extracted"
                   + (" alone (--genomes_only: the passes after, the organelle genes among them, left to a run without "
                      "it)" if opts.genomes_only else ", then the organelle genes"),
        "package": "the checks across species" + (", FCS-GX's screen" if opts.fcs_gx_db else "")
                   + ("" if opts.no_benchmarks else ", the benchmarks") + f", then OUT/{package_dir}/"
                   + ("" if opts.no_release else " and OUT/release/ (the website's data, the release's bundles)"),
        "release": f"OUT/release/ of OUT/{package_dir}/: the website's data and the release's bundles, BGZF"}
    if opts.wanted is None:
        species = "every species of the groups"
    else:
        taxids = sorted(opts.wanted)
        species = (f"a trial: the clades of {len(taxids)} taxids ({', '.join(map(str, taxids[:5]))}"
                   + (", ..." if len(taxids) > 5 else "") + ")")
    rows = [("stage", opts.stage + (": " + " -> ".join(stages) if len(stages) > 1 else ""))]
    rows += [(f"  {s}", does[s]) for s in stages]
    rows += [("groups", ", ".join(opts.groups)), ("species", species),
             ("loci", (f"the BUSCOs of --loci {opts.loci}" if opts.loci else f"every BUSCO of {odb_mod.LINEAGE}")
              + ("" if opts.no_markers else ", and Eukore's markers")),
             ("build_dir", os.path.abspath(opts.work)),
             ("out", os.path.abspath(opts.out) + (f"; the genomes in {os.path.abspath(opts.genomes)} (--genomes)"
                                                  if separate_genomes(opts) else ""))]
    if opts.scratch:
        rows.append(("local_tmp", os.path.abspath(opts.scratch)))
    cores = f"{opts.threads}"
    if "process" in stages:
        cores += f"; {opts.jobs} genomes at a time, {opts.threads_per_genome} threads each"
    if "fetch" in stages:
        cores += f"; {max(1, opts.downloads)} downloads at a time ({'curl' if util.CURL else 'urllib'})"
    rows.append(("threads", cores))
    if "process" in stages:
        rows.append(("memory", f"{util.size_text(opts.memory)} for the genomes' runs ("
                     + (whence if opts.memory_gb else f"{MEMORY_SHARE:.0%} of {whence}") + ")" if opts.memory else
                     "not known (--memory_gb): -j genomes at once whatever their sizes"))
    if "package" in stages:
        rows.append(("genomes", "every one kept (no --keep_gb)" if opts.keep_gb is None else
                     f"cut back to {opts.keep_gb:g} GB (--keep_gb)"))
    off = [what for what, flag in (("rescue pass", opts.no_rescue), ("transcriptomes", opts.no_tsa),
                                   ("tiara", opts.no_tiara), ("sketches", opts.no_sketches),
                                   ("benchmarks", opts.no_benchmarks),
                                   ("figures", opts.no_figures and not opts.no_benchmarks),
                                   ("read benchmark", not opts.read_genomes and not opts.no_benchmarks),
                                   ("other references in the read benchmark (--read_compare)",
                                    not opts.read_compare and not opts.no_benchmarks),
                                   ("intraspecific panel", not opts.intraspecific)) if flag]
    if off:
        rows.append(("off", ", ".join(off)))
    rows.append(("record", f"genome record {genomes.GENOME_RECORD}: "
                 + ("every done genome stamped with it, its calls kept (--reversion)" if opts.reversion else
                    "a done genome of an older one is run again")))
    rule = "=" * 100
    return "\n".join([rule, f" Eukore {__version__}: the standalone eukaryote reference's build",
                      f" {datetime.datetime.now():%Y-%m-%d %H:%M:%S} on {socket.gethostname()} (pid {os.getpid()}), "
                      f"Python {sys.version.split()[0]}",
                      " " + command_of(sys.argv[1:] if argv is None else argv), "-" * 100]
                     + [f" {k:<10}{v}" for k, v in rows] + [rule])


# the folders' flags before 2026-11-09, still taken as they were: --work is --build_dir, --scratch --local_tmp (so that
# a script's --scratch $SLURM_LOCAL_SCRATCH stays the runs' temporary files)
RENAMED = {"--work": "--build_dir", "--scratch": "--local_tmp"}


STAGES = ("fetch", "process", "package")  # in the build's order


def python_problem(version=None, filtered=None):
    """Why this Python can not run the build, or None (since 2026-11-10: said hours in, by the first use). The
    sources' archives are unpacked with tarfile's extraction filter (odb.py: filter="data"), which Python has from
    3.12, and from 3.11.4, 3.10.12 and 3.9.17 of the older lines."""
    import tarfile  # noqa: PLC0415
    version = sys.version_info if version is None else version
    filtered = hasattr(tarfile, "data_filter") if filtered is None else filtered
    if tuple(version[:2]) < (3, 9) or not filtered:
        return (f"Python {'.'.join(str(v) for v in version[:3])}: the build needs tarfile's extraction filter "
                "(Python 3.12, or the security releases 3.11.4, 3.10.12 and 3.9.17 on)")
    return None


def stages_of(text):
    """--stage's stages, in the build's order: `all` (every one), one stage, or stages that follow one another
    (`process,package`, `fetch,process`; since 2026-11-10), or `release` alone -> a tuple; None where it is none of
    these (`fetch,package`, `package,process`, a name unknown)."""
    names = tuple(s.strip() for s in text.split(","))
    if names == ("all",):
        return STAGES
    if names == ("release",):
        return names
    if names[0] in STAGES:
        first = STAGES.index(names[0])
        if names == STAGES[first:first + len(names)]:
            return names
    return None


def renamed(argv):
    """argv with the flags' old names (RENAMED) as their new ones, `--work=DIR` too."""
    out = []
    for a in argv:
        flag, eq, value = a.partition("=")
        out.append(RENAMED[flag] + eq + value if flag in RENAMED else a)
    return out


def main(argv=None):
    with contextlib.ExitStack() as held:  # the fetch stage's lock (util.single_run), let go however the run ends
        return _main(argv, held)


def _main(argv, held):
    p = argparse.ArgumentParser(description=f"Eukore {__version__}. " + __doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-V", action="version", version=f"Eukore {__version__}", help="Eukore's own version, and exit")
    p.add_argument("--build_dir", dest="work", metavar="WORK", required=True,
                   help="the build's own folder, kept for the whole build and every rerun: the sources the fetch stage "
                        "gets (sources/), each genome's results (done/<accession>/: a full build's weeks of runs), "
                        "the package stage's caches. Not a scratch folder: a file system that keeps its files "
                        "(--work before 2026-11-09)")
    p.add_argument("--out", required=True,
                   help="what the build makes: Eukore_<version>/ (the package), release/ (the website's data and the "
                        "downloads) and genomes/ (unless --genomes)")
    p.add_argument("--genomes", metavar="DIR",
                   help="the genomes' folder: the fetch stage downloads each species' genome, its gene models and "
                        "NCBI's assembly report into it, the process stage reads them from it, the package stage cuts "
                        "it back to --keep_gb (default OUT/genomes). Give every stage the same")
    p.add_argument("--stage", default="all", metavar="STAGE[,STAGE]",
                   help="fetch, process or package; all (the default): the three in turn; stages that follow one "
                        "another, comma-separated, in one run (process,package: a cluster's last job, one plan for "
                        "both; since 2026-11-10); release: OUT/release/ (the website's data and the downloads) of a "
                        "finished OUT/Eukore_<version>/ again, which the package stage writes last")
    p.add_argument("--genomes_only", action="store_true",
                   help="with --stage process: run the genomes, then the passes that take one genome at a time under "
                        "its claim (since 2026-11-10: the backfill with its clade hold-outs, the assembly metadata, "
                        "the done folders' packing), and stop: the passes over all genomes (organelle genes, "
                        "intraspecific panel, gene models) left to a process run without it. For several runs on one "
                        "WORK (a SLURM array: the genome claims share the genomes and those passes out); one run "
                        "without it after them all (--stage process,package: it waits for the genomes the others "
                        "still run, runs what they left, then makes the passes and the package; tools/eukore.sh "
                        "submits both)")
    p.add_argument("--done_folders", action="store_true",
                   help="keep each genome's results as a folder of 8-14 small files (WORK/done/<accession>/, as every "
                        "build before 2026-11-10) rather than one file (WORK/done/<accession>.zip, the default: far "
                        "fewer files on a shared file system; a WORK of folders is read as it is and packed by the "
                        "process stage)")
    p.add_argument("--no_release", action="store_true",
                   help="the package stage writes no OUT/release/ (the website's data and the downloads: a second "
                        "copy of the package, compressed)")
    p.add_argument("--version", default=datetime.date.today().strftime("%Y.%m"), help="the package's version (YYYY.MM)")
    p.add_argument("--groups", default=",".join(sources.ASSEMBLY_GROUPS),
                   help=f"NCBI's assembly groups whose species the package holds (default all: "
                        f"{','.join(sources.ASSEMBLY_GROUPS)}; e.g. fungi,protozoa first)")
    p.add_argument("--keep_gb", type=float,
                   help="the genomes kept in OUT/genomes (--genomes) by the package stage, GB, cut back to it before "
                        "the package is written (so that genomes.tsv says which are kept; the build may hold more "
                        "while it works); with it, a superseded assembly's and the intraspecific panel's genome files "
                        "are deleted once done too (default: no genome file deleted; until 2026-11-08 1000)")
    p.add_argument("--loci", help="call only these BUSCOs (a file, one per line or a loci table's first column; "
                                  "default: every BUSCO of the dataset)")
    p.add_argument("--species", help="a partial build (a trial run): only the species inside these NCBI taxids' clades, "
                                     "comma-separated or a file of one per line; the manifest records them")
    p.add_argument("--miniprot", default="miniprot")
    p.add_argument("--minimap2", default="minimap2", help="minimap2, the read benchmark's aligner")
    p.add_argument("--read_genomes", type=int, default=reads_mod.PER_GROUP,
                   help=f"the read benchmark's test genomes per supergroup, kept through the prune (default "
                        f"{reads_mod.PER_GROUP}; 0: no read benchmark)")
    p.add_argument("--read_compare", metavar="NAME=FASTA,TABLE", action="append",
                   help="another reference for the read benchmark and overlap.tsv (TABLE: record<TAB>NCBI taxid), "
                        "e.g. EukDetect2=eukdb.fasta,busco_taxid_genome_link.txt; again for each further one "
                        "(tools/compare_references.py: MetaPhlAn 4's eukaryote markers, PR2's 18S); the first is "
                        "breadth.tsv's")
    p.add_argument("--hmmsearch", default="hmmsearch")
    p.add_argument("--fcs_gx_db", help="NCBI FCS-GX's database (the gxdb folder's prefix, ~470 GB, best copied to "
                                       "/dev/shm; a node of ~512 GB): the package stage screens each genome's windows "
                                       "around its copies (copies.tsv's fcs)")
    p.add_argument("--run_gx", default="run_gx", help="FCS-GX's runner (bioconda's ncbi-fcs-gx)")
    p.add_argument("--tiara", default="tiara", help="tiara, which classifies the windows around the copies (default: "
                                                    "tiara on the PATH, else none; envs/tiara.yaml's, as it needs a "
                                                    "Python of 3.9 or older)")
    p.add_argument("--no_tiara", action="store_true", help="do not classify the windows around the copies")
    p.add_argument("--sourmash", default="sourmash", help="sourmash, which sketches each genome before the prune "
                                                          "(metadata/sketches.zip; default: sourmash on the PATH, else "
                                                          "none)")
    p.add_argument("--no_sketches", action="store_true", help="no sourmash sketches of the genomes")
    p.add_argument("--refresh", action="store_true",
                   help="fetch again what chooses the genomes and their gene models (NCBI's assembly tables, its "
                        "FCS-GX and ANI reports, Ensembl's, VEuPathDB's and TOGA2's indexes): a species' newer or "
                        "better genome is fetched and its old file deleted once the new one is done; a done genome "
                        "with better gene models than its run took is run again")
    p.add_argument("--veupathdb_key", default=os.environ.get(annotation.VEUPATHDB_KEY),
                   help="a VEuPathDB user's API key (its profile's Web Service Access), for VEuPathDB's gene models "
                        f"(default the environment's {annotation.VEUPATHDB_KEY}, the safer way: a key on the command "
                        "line is seen by every user of the node; without, none of its); sent to VEuPathDB's sites "
                        "alone, over https")
    p.add_argument("--no_rescue", action="store_true", help="no rescue pass of the BUSCOs miniprot missed (loci.py: "
                                                            "their HMMs in the genome's six-frame ORFs place them)")
    p.add_argument("--held_out_ranks", default=",".join(genomes.HELD_OUT_RANKS),
                   help="the clade hold-out of the OrthoDB benchmark: a BUSCO dataset species' genome called again "
                        "without its genus', family's, order's dataset species in turn, a call each (comma-separated; "
                        "`none`: the species' hold-out alone; default %(default)s)")
    p.add_argument("--reversion", action="store_true",
                   help="stamp every done genome of an older genome record (genomes.GENOME_RECORD: what made its "
                        "calls) with this Eukore's, its calls kept as they were made, so that it is not run again; "
                        "its done.json notes the record it had (`reversioned`). Without, such a genome is run again")
    p.add_argument("--no_markers", action="store_true",
                   help="no Eukore marker loci (RPB1, RPB2, EF-2, EF-1a, ...: data/markers.tsv), the BUSCOs alone")
    p.add_argument("--no_tsa", action="store_true",
                   help="no transcriptomes (ENA's TSA sets) for the species without a genome (default: their BUSCOs "
                        "called in them, source tsa)")
    p.add_argument("--intraspecific", type=int, default=3,
                   help="further assemblies of each species of data/intraspecific_panel.tsv run for the barcode gap "
                        "(intraspecific/; default 3; 0: none)")
    p.add_argument("--organelle_models", choices=organelle_mod.MODELS, default="miniprot",
                   help="the assemblies' organelle protein-coding genes: nhmmer's hits alone, or (default) miniprot's "
                        "exons by the organelle's code where nhmmer joined several hits (a gene with introns; "
                        "docs/methods.md)")
    p.add_argument("--ksgp_tax", help="a KSGP release's KSGP_v<V>_plus2.tax whose lineages species.tsv gives (protal's "
                                      f"KSGP, e.g. 4.4 or 5.2; default KSGP's served {labels.KSGP_VERSION}, fetched)")
    p.add_argument("--midori2", action="append",
                   help="a MIDORI2 RAW FASTA (its COI, e.g. MIDORI2_UNIQ_NUC_GB<release>_CO1_RAW.fasta.gz, from "
                        "reference-midori.info) whose lineages species.tsv gives, by NCBI taxid; given again for its "
                        "other genes' files (since 2026-11-10). Its records are screened (eukore/refscreen.py, with "
                        "minimap2 and vsearch): one that is another organism's gene or another gene names no species "
                        "(reference_screen.tsv)")
    p.add_argument("--no_reference_screen", action="store_true",
                   help="--midori2's records not screened: every one names its species (until 2026-11-10)")
    p.add_argument("--bold_tsv", help="a BOLD data package's TSV (its bin_uri and species columns): its BINs in "
                                      "crossrefs.tsv")
    p.add_argument("--no_benchmarks", action="store_true", help="write the package without its benchmarks/ "
                                                                  "(eukore/benchmarks.py; docs/benchmarks.md)")
    p.add_argument("--no_figures", action="store_true",
                   help="write the benchmarks without their figures (benchmarks/figures/: Figure 1 and the rest, "
                        "eukore/figures.py; since 2026-11-10; python3 -m eukore.figures draws them later)")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("-j", "--jobs", type=int, default=None,
                   help="genomes at a time at most (default threads / 4; fewer where --memory_gb holds fewer, each "
                        "then with more threads)")
    p.add_argument("--memory_gb", type=float,
                   help="the process stage's memory, GB: a genome is begun where its run's estimate (~7.5 bytes a "
                        "base of up to 4 Gb, and 3 GB) fits in what the running ones leave (default: SLURM's "
                        f"allocation, else the cgroup's limit, else the node's memory available; {MEMORY_SHARE:.0%}% "
                        "of it)")
    p.add_argument("--downloads", type=int, default=8,
                   help="the fetch stage's genomes and annotations downloaded at a time (default 8), by curl where on "
                        "the PATH (one NCBI connection: ~30-40 MB/s by curl, ~10 by Python's urllib)")
    p.add_argument("--local_tmp", dest="scratch", metavar="DIR",
                   help="a node's local disk (or a scratch file system without backups) for temporary files: the fetch "
                        "stage's unfinished downloads (.part files) and their gzip checks, each file copied into "
                        "OUT/genomes (--genomes) once whole; the process stage's (each genome run's copy of its "
                        "genome, miniprot's pieces of a large one, the rescue's ORFs, up to ~3 times a genome's size "
                        "each; and TMPDIR's, barrnap's and ITSx's), only the results written to --build_dir (default: "
                        "beside them, and TMPDIR; --scratch before 2026-11-09)")
    p.add_argument("--report_minutes", type=float, default=util.REPORT_EVERY / 60,
                   help="minutes between progress reports of the genomes' downloads and the process stage's passes: "
                        "items left, time left and, in the process stage, the cores, memory and disk used (default "
                        "10; a report every minute of a pass's first 10)")
    opts = p.parse_args(renamed(sys.argv[1:] if argv is None else argv))
    opts.groups = [g for g in opts.groups.split(",") if g]
    compare = []
    for spec in opts.read_compare or ():
        name, _eq, files = spec.partition("=")
        fasta, _comma, table = files.partition(",")
        if not (name and os.path.isfile(fasta) and os.path.isfile(table)) or name in reads_mod.REFERENCES or \
                name in {c[0] for c in compare}:
            p.error(f"--read_compare {spec}: NAME=FASTA,TABLE, both files there, a name of its own")
        compare.append((name, fasta, table))
    opts.read_compare = compare
    for path in opts.midori2 or ():
        if not os.path.isfile(path):
            p.error(f"--midori2 {path}: no such file")
    unknown = set(opts.groups) - set(sources.ASSEMBLY_GROUPS)
    if unknown:
        p.error(f"--groups: {sorted(unknown)} (known: {', '.join(sources.ASSEMBLY_GROUPS)})")
    why = python_problem()
    if why:
        p.error(why)
    opts.stages = stages_of(opts.stage)
    if opts.stages is None:
        p.error(f"--stage {opts.stage}: all, release, one of {', '.join(STAGES)}, or stages that follow one another "
                "(fetch,process or process,package)")
    if opts.genomes_only and opts.stages != ("process",):  # `all` would package without the passes
        p.error("--genomes_only: with --stage process alone")
    warning = None  # said after the banner
    if opts.scratch:  # checked now, not after the minutes of planning (since 2026-11-08)
        why = scratch_problem(opts.scratch)
        if why and opts.stages == ("fetch",):  # its .part files beside the genomes then, as without --local_tmp
            warning = f"WARNING: --local_tmp {why}: the .part files beside the genomes instead"
            opts.scratch = None
        elif why:  # the process stage's temporary files would load WORK's file system: stop
            p.error(f"--local_tmp {why}")
    opts.held_out_ranks = [] if opts.held_out_ranks in ("", "none") else opts.held_out_ranks.split(",")
    if set(opts.held_out_ranks) - set(genomes.HELD_OUT_RANKS):
        p.error(f"--held_out_ranks: {opts.held_out_ranks} (known: {','.join(genomes.HELD_OUT_RANKS)})")
    opts.held_out_ranks = [r for r in genomes.HELD_OUT_RANKS if r in opts.held_out_ranks]  # genus, family, order
    opts.jobs = opts.jobs or max(1, opts.threads // 4)
    opts.threads_per_genome = max(1, opts.threads // opts.jobs)
    opts.genomes = opts.genomes or os.path.join(opts.out, "genomes")  # --genomes since 2026-11-09
    wanted = opts.wanted = clades(opts.species)
    # the process stage's memory (since 2026-11-08): genomes begun while their estimates fit (schedule)
    memory, whence = (opts.memory_gb * 1e9, "--memory_gb") if opts.memory_gb else util.memory_budget()
    opts.memory = memory * (1 if opts.memory_gb else MEMORY_SHARE) if memory else None
    print(banner(opts, argv, whence), flush=True)  # the version and what the run will do (since 2026-11-08)
    if warning:
        say(warning)
    # every compressed file BGZF, its blocks compressed on the run's threads (since 2026-11-09)
    util.compression_threads(opts.threads)
    if opts.stages == ("release",):  # the website's data and the downloads of a finished package (since 2026-11-09)
        release_mod.write(os.path.join(opts.out, f"Eukore_{opts.version}"), os.path.join(opts.out, "release"), say,
                          threads=opts.threads)  # the downloads made at once (since 2026-11-10: one after another)
        return 0
    os.makedirs(opts.work, exist_ok=True)
    os.makedirs(opts.genomes, exist_ok=True)
    # the resources the build uses, said with the process stage's progress (since 2026-11-08)
    opts.usage = util.Usage([("WORK", opts.work), ("GENOMES" if separate_genomes(opts) else "OUT", opts.genomes)],
                            opts.threads)
    src_dir = os.path.join(opts.work, "sources")
    if "fetch" in opts.stages:
        # one fetch at a time on a WORK (since 2026-11-10: nothing stopped a second one, e.g. a job submitted twice;
        # the downloads have no claims of their own, and both wrote the same sources and genomes)
        try:
            held.enter_context(util.single_run(opts.work, "fetch"))
        except util.RunHeld as e:
            sys.exit(f"--stage fetch: {e}")
        fetch_sources(src_dir, opts.groups, opts.ksgp_tax, opts.refresh, opts.threads)
    ncbi = ncbi_mod.Ncbi(os.path.join(src_dir, "taxdump.tar.gz"))
    opts.inside = inside_of(ncbi, wanted)
    job_map, assemblies = plan(opts, src_dir, ncbi)
    extras = extras_of(opts, src_dir, ncbi, job_map, assemblies)  # the barcode gap's panel (since 2026-11-07)
    if opts.reversion:  # done genomes of an older genome record stamped with this one, kept (since 2026-11-09)
        olds = (genomes.reversion(opts.work, a) for a in sorted({**job_map, **extras}))
        stamped = collections.Counter(old for old in olds if old is not None)
        say(f"--reversion: {sum(stamped.values()):,} done genomes stamped with genome record {genomes.GENOME_RECORD}"
            + (f" (from {', '.join(f'{n:,} of record {r}' for r, n in sorted(stamped.items()))}), their calls kept"
               if stamped else ": none was of an older one"))
    # the genomes and their annotations fetched (--downloads at a time, through --local_tmp), by `all` too (until
    # 2026-11-06 `all` left them to the process stage's workers, --local_tmp and --downloads unused)
    if "fetch" in opts.stages:
        loci = called_loci(opts, src_dir)  # the BUSCOs and the markers: a genome not searched for one is fetched again
        # a done genome no longer here is fetched again where better gene models than its run took are there
        # (genomes.better_annotation): the process stage runs it again with them
        again = {a for a, j in job_map.items() if j["url"] and not j.get("path")
                 and not os.path.isfile(genomes.genome_file(opts.genomes, a)) and genomes.is_done(opts.work, a, loci)
                 and genomes.better_annotation(opts.work, j)}
        if again:
            say(f"{len(again):,} done genomes no longer here fetched again for better gene models (e.g. "
                f"{sorted(again)[:3]})")
        # and a done genome of an older genome record (since 2026-11-09): run again by this one, unless --reversion
        stale = {a for a, j in {**job_map, **extras}.items() if j["url"] and not j.get("path")
                 and not os.path.isfile(genomes.genome_file(opts.genomes, a)) and genomes.is_done(opts.work, a, loci)
                 and genomes.record_stale(opts.work, a)} - again
        if stale:
            say(f"{len(stale):,} done genomes of an older genome record than {genomes.GENOME_RECORD} no longer here "
                f"fetched again to be run again (e.g. {sorted(stale)[:3]}; --reversion keeps their calls instead)")
            again |= stale
        todo = [j for a, j in sorted({**job_map, **extras}.items()) if j["url"] and not j.get("path")
                and not os.path.isfile(genomes.genome_file(opts.genomes, a))
                and (not genomes.is_done(opts.work, a, loci) or a in again)]
        downloads = max(1, opts.downloads)
        say(f"fetching {len(todo):,} genomes ({sum(j.get('size') or 0 for j in todo) / 1e9:,.1f} Gb) into "
            f"{opts.genomes}, {downloads} at a time by {util.CURL or 'urllib'}, checked by "
            f"{util.GZIP_TEST or 'Python (no igzip or pigz)'}, .part files in {opts.scratch or 'there'}; "
            f"progress every {opts.report_minutes:g} minutes")
        # every --report_minutes: genomes left and the time left, by their sizes (a plant's genome is a thousand
        # yeasts'), and the download rate, all kept by this thread (the downloads' are at_a_time's)
        progress = util.Progress("genomes fetched", [j.get("size") or 0 for j in todo], say, opts.report_minutes * 60,
                                 early=util.early_reports(opts.report_minutes * 60))
        fetched = 0
        # each download's gzip check in a process of its own (a core apart from the downloads' threads; spawned, as
        # a process with threads is not forked), its .part file in --local_tmp
        checker = concurrent.futures.ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("spawn"))
        how = {"scratch": opts.scratch, "checker": checker}
        if opts.scratch:
            os.makedirs(opts.scratch, exist_ok=True)
        # killed runs' downloads (and TOGA2 .2bit files a killed conversion left), untouched for a day: not those of
        # another run downloading now (since 2026-11-08: every .part was deleted, a running download's too)
        for leftover in [p for folder in (opts.genomes, opts.scratch) if folder
                         for p in glob.glob(os.path.join(folder, "*.part")) + glob.glob(os.path.join(folder, "*.2bit"))]:
            with contextlib.suppress(OSError):
                if time.time() - os.path.getmtime(leftover) > LEFTOVER_AGE:
                    os.remove(leftover)

        def fetch_genome(j):
            # its annotation first (since 2026-11-08), so that the process stage never runs a genome whose gene
            # models are still to come (it ran them without, then again whole once they came)
            fetch_annotation(j)
            path, times = genomes.genome_file(opts.genomes, j["accession"]), collections.Counter()
            why = genomes.fetch_genome(j, path, times=times, **how)
            return j, why, 0 if why else os.path.getsize(path), times

        def fetch_annotation(j):
            """A job's annotation beside its genome (annotation.fetch) -> (path, source); where none has one, a note
            (annotation.note_none) and, for a done genome, the sources in its done.json. Raises OSError where one
            did not arrive (a later run asks again)."""
            acc = j["accession"]
            if not annotation.wanted(j) or annotation.found(opts.genomes, acc)[0] \
                    or annotation.none_noted(opts.genomes, j):
                return annotation.found(opts.genomes, acc)
            try:
                path, source = annotation.fetch(j, opts.genomes, opts.veupathdb_key, **how)
            except OSError:
                return None, None
            if path is None:
                annotation.note_none(opts.genomes, acc, annotation.none_before(j, None))
            # the sources asked before the one found (all where none had one): a done genome's not asked again
            genomes.note_annotation_none(opts.work, j, annotation.none_before(j, source))
            return path, source

        spent = collections.Counter()  # the downloads' seconds of transfer, check and copy (fetch_file's times)
        for done in at_a_time(fetch_genome, todo, downloads):
            j, why, size, times = done.result()
            fetched += size
            spent += times
            rate = fetched / max(1e-9, time.time() - progress.started) / 1e6
            busy = sum(spent[k] for k in ("transfer", "check", "copy")) or 1e-9
            each = spent["bytes"] / max(1e-9, spent["transfer"]) / 1e6  # one transfer's rate
            progress.step(j.get("size") or 0, failed=bool(why),
                          note=f"; {fetched / 1e9:,.1f} GB at {rate:,.1f} MB/s (a download's time: transfer "
                               f"{spent['transfer'] / busy:.0%} at {each:,.1f} MB/s, "
                               f"check {spent['check'] / busy:.0%}, copy {spent['copy'] / busy:.0%})")
        # the genomes' own gene models beside them, for the process stage's copies' models and their comparison
        # (annotation.py): NCBI's (RefSeq's, a GenBank submitter's; not where the assembly summary gives no
        # protein-coding genes), else Ensembl's, else TOGA2's
        # (and a done genome's better one than its run took: genomes.better_annotation)
        todo = [j for a, j in sorted(job_map.items()) if annotation.wanted(j)
                and not annotation.found(opts.genomes, a)[0] and not annotation.none_noted(opts.genomes, j)
                and (not store.has(genomes.done_dir(opts.work, a), annotation.TABLE)
                     or genomes.better_annotation(opts.work, j))]
        say(f"fetching {len(todo):,} genomes' annotations (the first of {', '.join(annotation.SOURCES)} that has one) "
            f"into {opts.genomes}")
        progress = util.Progress("annotations fetched", [1] * len(todo), say, opts.report_minutes * 60, "", 1,
                                 early=util.early_reports(opts.report_minutes * 60))
        def get(j):
            return j, annotation.fetch(j, opts.genomes, opts.veupathdb_key, **how)
        for done in at_a_time(get, todo, downloads):
            try:
                j, (path, source) = done.result()
                if path is None:  # none has one: noted beside the genome for its run (annotation.none_noted)
                    annotation.note_none(opts.genomes, j["accession"], annotation.none_before(j, None))
                # the sources asked before the one found (all where none had one): a done genome's not asked again
                genomes.note_annotation_none(opts.work, j, annotation.none_before(j, source))
                progress.step(1)
            except OSError as e:
                progress.step(1, failed=True, note=f" ({e})")
        checker.shutdown()
        held.close()  # the fetch is through: its lock let go before the stages after it
        if opts.stages == ("fetch",):
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
    if not opts.no_rescue and not shutil.which("esl-translate"):
        say("WARNING: no esl-translate (hmmer's easel tools): no rescue pass of the BUSCOs miniprot missed")
    # the tools' and Python's temporary files in --local_tmp too (the workers inherit it); the package stage's as
    # well (since 2026-11-10: the process stage's alone, a run of the package stage by itself wrote them to /tmp)
    if opts.scratch:
        os.makedirs(opts.scratch, exist_ok=True)
        os.environ["TMPDIR"] = tempfile.tempdir = os.path.abspath(opts.scratch)
    if "process" in opts.stages:
        # a run without the fetch stage fetches nothing: a genome not here waits for a later run
        settings = settings_of(opts, src_dir, ncbi, staged_only="fetch" not in opts.stages)
        began = time.time()
        say(f"process stage: {opts.jobs} genomes at a time, {opts.threads_per_genome} threads each, progress every "
            f"{opts.report_minutes:g} minutes" + opts.usage.line())
        waiting, failed = run_genomes(opts, job_map, settings, extras)
        # --genomes_only (since 2026-11-09): one of several runs on the WORK (a SLURM array), the passes over all
        # genomes left to one run after them all (until then each run made them over the genomes done so far, the
        # others still running). The passes that take a genome at a time under its claim are made by every run (since
        # 2026-11-10: by the last run alone, on one node: the backfill's clade hold-outs are whole genome calls):
        # the array's tasks share them out as they share the genomes, the last run makes what they left
        pack_pass(opts)
        backfill_pass(opts, job_map, settings)
        if opts.genomes_only:
            describe_pass(opts, job_map)
        else:
            intraspecific_pass(opts, src_dir, ncbi, extras)
            describe_pass(opts, job_map)
            drop_superseded(opts, job_map)
            organelle_genes(opts, src_dir, job_map, ncbi)
            annotate_pass(opts, job_map, may_fetch="fetch" in opts.stages)
        say(f"process stage: took {util.duration(time.time() - began)}"
            + (" (--genomes_only: the passes over all genomes left to a process run without it)"
               if opts.genomes_only else "") + opts.usage.line())
        # a rerun mends a genome not fetched yet or a download that may come; not one NCBI no longer serves, nor a
        # call that failed (it fails again): `all` packages without those (until 2026-11-06 it waited for them forever)
        retry = waiting + [(a, w) for a, w in failed if w.startswith("fetch failed: ") and not placement.given_up(
            {"status": "failed: " + w[len("fetch failed: "):]})]
        if failed:
            say(f"{len(failed)} genomes failed (e.g. {failed[:3]})")
        if retry:
            # a run that fetches (`all`) stops here, its rerun fetches them; one that does not (a cluster's
            # process,package, since 2026-11-10) can not mend them and packages without them, as --stage package does
            onward = "package" in opts.stages and "fetch" not in opts.stages
            say(f"{len(waiting)} genomes not fetched yet, {len(retry) - len(waiting)} downloads to try again (e.g. "
                f"{retry[:3]}): " + ("the package is written without them (--stage fetch, then this run again, "
                                     "adds them)" if onward else
                                     "run again (--stage fetch first on a cluster), or --stage package to write the "
                                     "package without them"))
            if not onward:
                return 0 if opts.stages == ("process",) else 1
    stages = []  # the package stage's steps and their seconds (benchmarks' stages.tsv)
    if opts.stages == ("package",):  # a tool missing at a genome's run, here now (else the process stage's pass)
        pack_pass(opts)
        backfill_pass(opts, job_map, settings_of(opts, src_dir, ncbi, staged_only=True))
    if "package" in opts.stages and opts.fcs_gx_db:
        todo = [a for a in sorted(job_map) if genomes.is_done(opts.work, a)
                and not store.has(genomes.done_dir(opts.work, a), "fcs.tsv")]
        # several genomes at once, FCS_THREADS cores each (since 2026-11-08: one at a time on every thread, each
        # run's ~13 Mb of windows mostly its start-up); the database in /dev/shm is mapped once for them all
        at_once = max(1, min(len(todo), opts.threads // FCS_THREADS))
        say(f"FCS-GX: {len(todo)} genomes' windows to screen ({opts.fcs_gx_db}), {at_once} at a time")
        failed = 0
        with timed(stages, "FCS-GX"):
            def screen(acc):
                return genomes.fcs_screen(opts.work, acc, job_map[acc]["taxid"], opts.fcs_gx_db, opts.run_gx,
                                          max(1, opts.threads // at_once), opts.scratch)
            for i, done in enumerate(at_a_time(screen, todo, at_once), 1):
                failed += done.result() is None
                if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                    say(f"  {i}/{len(todo)} screened ({failed} failed)")
    if "package" in opts.stages:
        vsearch = shutil.which("vsearch")
        if not vsearch:
            say("WARNING: no vsearch: no checks across species (own_species, same_as, the rRNA part screen), every "
                "record's own_species undecided")
        # killed runs' unfinished packages and read benchmark folders; not those of a package stage going on (since
        # 2026-11-08: every one was deleted)
        for leftover in (util.dead_partials(opts.out, f"Eukore_{opts.version}")
                         + util.dead_partials(opts.work, "read_benchmark")
                         + (util.dead_partials(opts.scratch, "read_benchmark") if opts.scratch else [])):
            shutil.rmtree(leftover, ignore_errors=True)
        loci = called_loci(opts, src_dir)
        short = [a for a in job_map if genomes.is_done(opts.work, a) and not genomes.is_done(opts.work, a, loci)]
        if short:
            say(f"WARNING: {len(short)} genomes were not searched for all {len(loci)} loci called (run before a locus "
                f"was, e.g. {short[:3]}): their call rates count them out; run --stage process again for them all")
        stale = sorted(a for a in job_map if genomes.is_done(opts.work, a) and genomes.record_stale(opts.work, a))
        if stale:  # since 2026-11-09: their calls go into the package as they are, genomes.tsv's genome_record says so
            say(f"WARNING: {len(stale):,} genomes' calls are of an older genome record than {genomes.GENOME_RECORD} "
                f"(e.g. {stale[:3]}): run --stage process again for them (their genomes fetched again where pruned), "
                "or --reversion to keep them; genomes.tsv's genome_record names each one's")
        undone = sorted(a for a in job_map if not genomes.is_done(opts.work, a))
        if undone:
            say(f"WARNING: {len(undone)} genomes of the plan are not done (e.g. {undone[:3]}): the package is written "
                "without them (manifest: inputs.genomes_not_done)")
        marker_folder = markers_of(opts, src_dir)
        if marker_folder and not markers_mod.current(src_dir):
            # since 2026-11-10 also a set built by other rules or from other input files (markers.ready's recipe and
            # inputs), not only one of before the decoys
            say("WARNING: the marker set was not built by this code's rules, tables and input files (a set of "
                "before data/marker_decoys.tsv's decoys, EFL for EF-1a, among them): run --stage fetch to build it "
                "again; until then the old set serves, a genome holding EFL alone may have it as its TEF1")
        with timed(stages, "BUSCO copies and their checks"):
            # the marker copies checked again by the set's profiles, decoys too (since 2026-11-09: markers.recheck)
            collected = package.collect(ncbi, job_map, assemblies, opts.work, vsearch, opts.threads,
                                        markers_mod.kinds(marker_folder),
                                        markers_mod.recheck(marker_folder, opts.hmmsearch, opts.threads))
        with timed(stages, "introns and genetic codes"):
            hmmalign = shutil.which("hmmalign")
            if not hmmalign:
                say("WARNING: no hmmalign: introns without their places in the genes, no genetic code check")
            hmm_of = odb_of(opts, src_dir).hmm
            package.genetics(collected, ncbi, hmm_of, hmmalign, opts.work, say, opts.threads)
        with timed(stages, "rRNA records and their checks"):
            records = genome_rrna_records(opts, job_map, ncbi) + published(opts, src_dir, ncbi)
            # each 18S's SSU model (since 2026-11-09): a bacterium's or an organelle's SSU is no 18S
            models, cmsearch = regions_mod.models_path(src_dir), shutil.which("cmsearch")
            ssu = functools.partial(regions_mod.ssu_models, cm=models, cmsearch=cmsearch, threads=opts.threads,
                                    tmp=opts.work) if cmsearch and regions_mod.is_cm(models) else None
            if ssu is None:
                say("WARNING: no cmsearch (infernal) or no Rfam SSU models (--stage fetch fetches them): the 18S not "
                    "told from a bacterium's or an organelle's SSU (records.tsv's ssu_model -)")
            records = rrna.finish(records, ncbi, vsearch, opts.threads, log=say, ssu_models=ssu)
        extra_tables = []
        with timed(stages, "18S variable regions and amplicons"):
            cm, cmalign = regions_mod.cm_path(src_dir), shutil.which("cmalign")
            if cmalign and os.path.isfile(cm):
                extra_tables.append(("rrna/regions.tsv", regions_mod.REGION_COLUMNS,
                                     regions_mod.run(records, cm, cmalign, opts.threads, opts.work, say)))
            else:
                say("WARNING: no cmalign (infernal) or no Rfam SSU model (--stage fetch fetches it): no 18S variable "
                    "regions")
            extra_tables.append(("rrna/amplicons.tsv", primers_mod.AMPLICON_COLUMNS,
                                 primers_mod.run(records, "rrna", threads=opts.threads)))
        with timed(stages, "organelle genes"):
            chosen = organelle_genes(opts, src_dir, job_map, ncbi)  # before the prune: the genomes searched while there
            held = {int(g["species_taxid"]) for g in collected["genomes"].values() if g["representative"] == "t"}
            organelle = organelle_mod.collect(ncbi, job_map, held | {int(r["species_taxid"]) for _n, _s, r in records},
                                              chosen, os.path.join(opts.work, "organelle"))
            organelle = organelle_mod.check(organelle, ncbi, vsearch, opts.threads)
            extra_tables.append(("organelle/amplicons.tsv", primers_mod.AMPLICON_COLUMNS,
                                 primers_mod.run(organelle, "organelle", threads=opts.threads)))
        writers = []
        if extras:
            with timed(stages, "intraspecific barcode gap"):
                writers.append(intraspecific_writer(opts, ncbi, extras, collected, records, organelle, vsearch))
        inputs = {"odb": odb_mod.state_of(src_dir), "taxdump_md5": util.md5(os.path.join(src_dir, "taxdump.tar.gz")),
                  "sources": {n: sources.state_of(src_dir, n).get("version") for n in ("EUKARYOME", "ROD")},
                  "assembly_tables": {os.path.basename(p): util.md5(p) for p in sorted(
                      tables(src_dir, "genbank", opts.groups) + tables(src_dir, "refseq", opts.groups))},
                  "genomes_not_done": len(undone),
                  "groups": opts.groups, "no_assembly": util.md5(NO_ASSEMBLY) if os.path.isfile(NO_ASSEMBLY) else None,
                  "ncbi_quality": {os.path.basename(p): util.md5(p) for p in quality.paths(src_dir)},
                  "refseq_organelle": {os.path.basename(p): util.md5(p)
                                       for p in organelle_mod.release_files(src_dir) if os.path.isfile(p)},
                  "markers": markers_mod.state_of(markers_of(opts, src_dir)),
                  "tsa": genomes.tsa_state(src_dir) if not getattr(opts, "no_tsa", False) else None}
        if wanted is not None:
            inputs["partial"] = {"species_clades": sorted(wanted)}  # a trial package: these clades only
        # every genome is called and extracted: the kept ones cut back to the budget, then the package names those kept
        firsts = {a for r in species.read_table() for a in species.accessions(r)}
        annotate_pass(opts, job_map, may_fetch="fetch" in opts.stages)
        # the read benchmark's genomes (reads.py), chosen while all are on disk and kept through the prune
        reps = {a: (int(g["species_taxid"]), genomes.genome_file(opts.genomes, a))
                for a, g in collected["genomes"].items() if g["representative"] == "t" and not job_map[a].get("path")}
        tests = reads_mod.choose(ncbi, reps, opts.read_genomes) if opts.read_genomes and not opts.no_benchmarks else []
        with timed(stages, "assembly metadata"):
            describe_pass(opts, job_map)  # before the prune: measured while the genomes are there
        drop_superseded(opts, job_map)
        # without --keep_gb, every genome kept (since 2026-11-08; the manifest written all the same)
        rows_kept = genomes.prune(opts.genomes, opts.work, float("inf") if opts.keep_gb is None else opts.keep_gb * 1e9,
                                  set(tests), firsts)
        kept = sum(s for _a, s, st in rows_kept if st != "deleted")
        say(f"genomes kept in {opts.genomes}: {kept / 1e9:.1f} GB ("
            + ("no --keep_gb: none deleted)" if opts.keep_gb is None else
               f"--keep_gb {opts.keep_gb:g}), {sum(st == 'deleted' for *_x, st in rows_kept)} deleted"))
        screened = None
        if opts.midori2 and not opts.no_reference_screen:  # since 2026-11-10: MIDORI2's misnamed, chimeric records
            with timed(stages, "reference screen"):
                screened = reference_screen(opts, src_dir, ncbi)
            if screened:
                extra_tables.append(("reference_screen.tsv", refscreen.TABLE_COLUMNS,
                                     refscreen.table_rows("MIDORI2", screened[0])))
        with timed(stages, "reference labels"):
            references = labels.references(
                src_dir, sources.ready(src_dir, "EUKARYOME"), opts.ksgp_tax, say, opts.midori2,
                (refscreen.dropped(screened[0]), refscreen.stamp(screened[1])) if screened else None)
            held_species = {int(g["species_taxid"]) for g in collected["genomes"].values()
                            if g["representative"] == "t"} | {int(r["species_taxid"]) for _n, _s, r in records}
            extra_tables.append(("crossrefs.tsv", labels.CROSSREF_COLUMNS, labels.crossrefs(
                ncbi, held_species, references, os.path.join(labels.folder_of(src_dir), labels.OTT_FILE),
                opts.bold_tsv, say)))
            groups_mod.check(ncbi, say)  # data/supergroups.tsv's clades against this taxdump
        hook = None if opts.no_benchmarks else functools.partial(
            measure, opts=opts, src_dir=src_dir, ncbi=ncbi, collected=collected, vsearch=vsearch, stages=stages,
            fates=genomes.fates(ncbi, assemblies, job_map, opts.inside if wanted is not None else None), tests=tests,
            references=references)
        meta = functools.partial(metadata_mod.write, work=opts.work, tables=tables(src_dir, "genbank", opts.groups)
                                 + tables(src_dir, "refseq", opts.groups), jobs=job_map)
        out = package.write(os.path.join(opts.out, f"Eukore_{opts.version}"), opts.version, ncbi, collected,
                            loci_items(opts, src_dir), records, opts.genomes, inputs, organelle, references, hook,
                            metadata=meta, tables=extra_tables, writers=writers, threads=opts.threads)
        say(f"Eukore: {out} ({json.dumps(package.read(out)['counts'])})")
        if not opts.no_release:  # the website's data and the downloads (since 2026-11-09; --stage release again)
            release_mod.write(out, os.path.join(opts.out, "release"), say, threads=opts.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
