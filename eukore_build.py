#!/usr/bin/env python3
"""Build Eukore, the standalone eukaryote reference (docs/FORMAT.md; since 2026-11-05).

    python3 eukore_build.py --work WORK --out OUT -t 16                 # all stages on a node with internet
    python3 eukore_build.py --work WORK --out OUT --stage fetch -t 8    # a cluster: the node with internet,
    python3 eukore_build.py --work WORK --out OUT --stage process -t 32 #   then any node (nothing fetched),
    python3 eukore_build.py --work WORK --out OUT --stage package       #   then the package

Stages:
- fetch: the sources into WORK/sources (NCBI's GenBank and RefSeq assembly tables of the --groups, NCBI's taxdump,
  BUSCO's eukaryota_odb12, EUKARYOME 2.1, ROD v1.2, RefSeq's organelle genomes, NCBI's FCS-GX and ANI reports, which
  rank the genomes: eukore/quality.py, PR2's and KSGP's taxonomies), then each species' genome (genomes.jobs) into
  OUT/genomes;
- process: each genome's BUSCO genes called (a genome searched for fewer than those called is run again) and its rRNA
  extracted (WORK/done/<accession>/); a genome not fetched yet waits for the next run; then the organelle genes
  (mitochondrial and plastid rRNA, mitochondrial COX1) of RefSeq's organelle genomes, and of the assemblies' own
  organelle molecules where RefSeq has none (WORK/organelle);
- package: the copies' checks across species, the published rRNA reads' regions (WORK/published), the organelle genes
  of any assembly not searched yet, OUT/genomes cut back to --keep_gb (genomes.prune: the model species' and the
  smallest kept first), then OUT/Eukore_<version>/.
`all` runs them in turn. Every stage is resumable: a rerun does only what is not done.
"""

import argparse
import collections
import concurrent.futures
import contextlib
import datetime
import functools
import glob
import json
import multiprocessing
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eukore import genomes, ncbi as ncbi_mod, odb as odb_mod, package, rrna, sources, species, util  # noqa: E402
from eukore import annotation, benchmarks, labels, metadata as metadata_mod, organelle as organelle_mod  # noqa: E402
from eukore import quality  # noqa: E402
from eukore import placement, reads as reads_mod  # noqa: E402
from eukore import __version__  # noqa: E402

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
NO_ASSEMBLY = os.path.join(sources.DATA, "no_assembly.fna.gz")
PUBLISHED_BATCH = 250  # published reads a barrnap + ITSx run
# the published reads' regions, kept in WORK/published: a new name where their borders change (2: ITSx's guessed 5.8S
# borders no longer taken, 2026-11-06; 3: regions on the record, not the unit, and the read's strand, 2026-11-06), so
# that a rerun annotates them again
PUBLISHED_REGIONS = "regions.3.tsv"
READ_BACKGROUND = os.path.join(sources.DATA, "read_background.tsv")  # the read benchmark's bacteria and archaea


def say(text):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {text}", flush=True)


# --- sources -----------------------------------------------------------------------------------------------------

def fetch_sources(src_dir, groups, ksgp_tax=None, refresh=False):
    """Every source into src_dir (kept when there and whole). refresh (--refresh): those that choose the genomes and
    their gene models fetched again (NCBI's assembly tables, its FCS-GX and ANI reports, the annotation indexes; each
    old one kept where its download fails), so that the plan takes NCBI's newer assemblies and annotations."""
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
        if url and not os.path.isfile(path) and opts.stage in ("all", "fetch"):
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
    say(f"genome set: {len(job_map)} genomes ({dict(genomes.why_counts(job_map.values()))}); "
        f"{sum(j['size'] or 0 for j in job_map.values()) / 1e9:.0f} Gb")
    # a species' genome another assembly takes the place of (TOGA2's; NCBI's newer, after --refresh): its file
    # deleted once the new one is done (drop_superseded)
    changed = genomes.supersede(opts.work, job_map)
    if changed:
        say(f"genome set: {len(changed)} species' genomes changed since the last plan (e.g. {changed[:3]}); the old "
            "files are deleted once the new ones are done")
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


def loci_items(opts, src_dir):
    """The BUSCOs called as package.write takes them: locus, OrthoDB name, cutoff, the HMM's length."""
    folder = odb_mod.folder_of(src_dir)
    odb = odb_mod.Odb(folder, loci_of(opts, src_dir))
    return [dict(locus=l, name=odb.names.get(l, "-"), cutoff=odb.cutoff[l], hmm_length=odb.length[l])
            for l in odb.loci]


def settings_of(opts, src_dir, ncbi, staged_only):
    folder = odb_mod.folder_of(src_dir)
    odb = odb_mod.Odb(folder, loci_of(opts, src_dir))
    hmm_db = odb.hmm_db(os.path.join(opts.work, "busco_loci.hmm"))
    tiara = None if opts.no_tiara else shutil.which(opts.tiara)
    return dict(work=opts.work, genomes=opts.genomes, miniprot=opts.miniprot, hmmsearch=opts.hmmsearch,
                threads=opts.threads_per_genome, staged_only=staged_only, odb=folder, loci=odb.loci, hmm_db=hmm_db,
                odb_lineages=odb_mod.species_lineages(odb, ncbi), odb_own=odb_mod.own_taxids(odb, ncbi), tiara=tiara,
                rescue=not opts.no_rescue and bool(shutil.which("esl-translate")),
                esl_translate=shutil.which("esl-translate"), sourmash=sourmash_of(opts),
                veupathdb_key=opts.veupathdb_key,
                rrna=functools.partial(rrna_hook, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"),
                                       threads=opts.threads_per_genome) if shutil.which("barrnap") else None)


def sourmash_of(opts):
    """sourmash, which sketches each genome (metadata.SKETCH), or None (--no_sketches, or not on the PATH)."""
    return None if opts.no_sketches else shutil.which(opts.sourmash)


def describe_pass(opts, job_map):
    """The metadata of each done genome still on disk that lacks it (metadata.py: a genome run before 2026-11-07, or
    whose measure failed), before the prune deletes the genomes -> how many were described, failed."""
    sourmash = sourmash_of(opts)
    todo = []
    for acc, j in sorted(job_map.items()):
        folder, path = genomes.done_dir(opts.work, acc), genomes.genome_file(opts.genomes, acc)
        if not j.get("path") and genomes.is_done(opts.work, acc) and os.path.isfile(path) \
                and metadata_mod.missing(folder, sourmash):
            todo.append((acc, path, folder, j.get("fetch", acc), sourmash))
    if not todo:
        return 0, 0
    say(f"assembly metadata: {len(todo)} genomes to describe" + ("" if sourmash else " (no sourmash: no sketches)"))
    failed = []
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for i, (acc, why) in enumerate(pool.imap_unordered(metadata_mod.describe_task, todo), 1):
            if why:
                failed.append((acc, why))
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                say(f"  {i}/{len(todo)} described ({len(failed)} failed{', e.g. ' + str(failed[0]) if failed else ''})")
    return len(todo) - len(failed), len(failed)


def drop_superseded(opts, job_map):
    """The genome files another assembly took the place of, deleted once it is done (genomes.drop_superseded)."""
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


def run_genomes(opts, job_map, settings):
    # a done genome of a BUSCO dataset species without its held-out calls (run before 2026-11-06) is run again where
    # its genome is still there: never fetched again for a benchmark, its OrthoDB benchmark then left without them
    # and a done genome with better gene models than its run took (genomes.better_annotation: an index naming its
    # assembly since, a VEuPathDB key given since, its annotation fetched since), its genome and that annotation here
    # (since 2026-11-07: until then it kept miniprot's models, the annotation compared alone)
    better = {a for a, j in job_map.items() if genomes.is_done(opts.work, a, settings["loci"])
              and genomes.rerun_for_annotation(opts.work, j, opts.genomes, not settings["staged_only"])}
    todo = [(j, settings) for a, j in sorted(job_map.items()) if not genomes.is_done(opts.work, a, settings["loci"])
            or a in better or (genomes.held_out_missing(opts.work, j, settings)
                               and os.path.isfile(j.get("path") or genomes.genome_file(opts.genomes, a)))]
    say(f"genomes to run: {len(todo)} of {len(job_map)}" + (f", {len(better)} of them again for better gene models"
                                                            if better else "")
        + ("" if settings["tiara"] else " (no tiara: windows not classified)"))
    waiting, failed = [], []
    for leftover in glob.glob(os.path.join(opts.work, "done", ".*.partial")):  # killed runs' (one build at a time)
        shutil.rmtree(leftover, ignore_errors=True)
    # every --report_minutes: genomes left and the time left, by their sizes (a genome not fetched yet is no work)
    progress = util.Progress("genomes run", [j.get("size") or 0 for j, _s in todo], say, opts.report_minutes * 60)
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for acc, _summary, why in pool.imap_unordered(run_one, todo):
            if why:
                (waiting if why == "not fetched yet" else failed).append((acc, why))
            progress.step(job_map[acc].get("size") or 0, failed=bool(why) and why != "not fetched yet",
                          skipped=why == "not fetched yet", note=f", {len(waiting):,} not fetched yet" if waiting else "")
    return waiting, failed


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
    cache = os.path.join(cache_dir, PUBLISHED_REGIONS)
    known = read_regions(cache)
    # by ITSx's profile sets, so that a batch's reads share theirs (an ITSx run per distinct set)
    todo = sorted(((n, s, rrna.itsx_groups(ncbi, t)) for n, s, t, *_ in reads if n not in known),
                  key=lambda r: (r[2], r[0]))
    if todo and shutil.which("barrnap"):
        say(f"published reads: {len(todo)} to annotate (barrnap, ITSx)")
        # small batches over every thread: barrnap (nhmmer) runs one thread a batch and dominates (~0.5 s a read)
        batches = [todo[i:i + PUBLISHED_BATCH] for i in range(0, len(todo), PUBLISHED_BATCH)]
        work = functools.partial(_annotate_batch, barrnap=shutil.which("barrnap"), itsx=shutil.which("ITSx"))
        with multiprocessing.Pool(max(1, opts.threads)) as pool, open(cache, "a", newline="\n") as fh:
            if fh.tell() and not ends_in_newline(cache):  # a killed run's last line cut: the next starts anew
                fh.write("\n")
            for got in pool.imap_unordered(work, batches):
                fh.writelines(f"{name}\t{regions}\t{complete}\t{strand}\n"
                              for name, (regions, complete, strand) in got.items())
                fh.flush()
                known.update(got)
    for name, seq, taxid, source, how, accession in reads:
        regions, complete, strand = known.get(name, ("-", "f", "+"))
        regions_list = rrna.parse_regions(regions)
        seq = rrna.revcomp(seq) if strand == "-" else seq  # the record on its unit's strand, as its regions
        out.append((name, seq, dict(kind="rDNA" if regions_list else "18S", accession=accession, source=source,
                                    taxid=taxid, species_taxid=ncbi.species_of(taxid), regions_list=regions_list,
                                    lsu_complete=complete == "t", label_from=how)))
    return out


def _annotate_batch(batch, barrnap=None, itsx=None):
    """[(name, sequence, ITSx's sets)] -> {name: (regions text, lsu_complete t/f, strand)}."""
    found = rrna.published_regions([(n, s) for n, s, _t in batch], barrnap, itsx, sets={n: t for n, _s, t in batch})
    return {name: (rrna.regions_text(found[name][0]), "t" if found[name][1] else "f", found[name][3])
            if name in found else ("-", "f", "+") for name, _s, _t in batch}


def ends_in_newline(path):
    with open(path, "rb") as fh:
        fh.seek(-1, os.SEEK_END)
        return fh.read(1) == b"\n"


def read_regions(path):
    """The published reads' regions cache -> {name: (regions text, lsu_complete t/f, strand)}; a line a killed run cut
    (not 4 fields, or no line end) left out, so that its read is annotated again (until 2026-11-06 it stopped every
    later run)."""
    known = {}
    if os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if line.endswith("\n") and len(f) == 4 and f[3] in ("+", "-"):
                    known[f[0]] = (f[1], f[2], f[3])
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


def annotate_pass(opts, job_map, may_fetch):
    """The done genomes with an annotation (NCBI's, Ensembl's or TOGA2's) compared with it (annotation.py: its GFF
    beside the genome, from the fetch stage, or fetched here where may_fetch), the GFF deleted after."""
    todo = [j for _a, j in sorted(job_map.items()) if annotation.missing(opts.work, j)]
    if not todo:
        return
    say(f"gene models: {len(todo)} genomes' calls to compare with their annotation")
    compared, left, failed = 0, 0, []
    for j in todo:
        try:
            n = annotation.annotate_genome(opts.work, j, opts.genomes, may_fetch, opts.veupathdb_key)
        except OSError as e:
            failed.append((j["accession"], str(e)[:100]))
            continue
        compared, left = compared + (n is not None), left + (n is None)
    say(f"  {compared} genomes compared, {left} without their annotation fetched yet (--stage fetch), "
        f"{len(failed)} failed{' (e.g. ' + str(failed[:2]) + ')' if failed else ''}")


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
        organelle_mod.assembly_pass(folder, job_map, lambda a: genomes.genome_file(opts.genomes, a), chosen, ncbi,
                                    barrnap, nhmmer, opts.jobs, opts.threads_per_genome, say)
    else:
        say("WARNING: no barrnap or nhmmer: the assemblies' organelle molecules not searched")
    return chosen


@contextlib.contextmanager
def timed(stages, step):
    """A step of the package stage, its seconds added to stages ([(step, seconds)]: benchmarks' stages.tsv)."""
    started = time.time()
    yield
    stages.append((step, time.time() - started))


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
    """The package's benchmarks (eukore/benchmarks.py) into its unfinished folder -> their figures, or None where
    they failed: a benchmark never stops the package. fates: genomes.fates (selection.tsv); tests: the read
    benchmark's genomes (reads.choose), kept through the prune; references: labels.references (overlap.tsv)."""
    nearest = {(r["accession"], r["locus"]): (r.get("nearest_id"), r.get("nearest_depth"))
               for r in collected["copy_rows"]} if vsearch else None
    minimap2 = shutil.which(opts.minimap2)
    if tests and not minimap2:
        say(f"WARNING: {opts.minimap2} not found: no read benchmark (bioconda's minimap2)")
    reads = dict(genomes_folder=opts.genomes, tests=list(tests), minimap2=minimap2,
                 background=reads_mod.read_background(READ_BACKGROUND, os.path.join(src_dir, "read_background")),
                 work=os.path.join(opts.work, "read_benchmark"))
    inside = opts.inside if opts.wanted is not None else None
    try:
        return benchmarks.run(part, ncbi, odb_mod.folder_of(src_dir), vsearch, opts.threads, opts.work, nearest,
                              stages, log=say, fates=fates, compare=opts.read_compare, inside=inside, reads=reads,
                              databases=databases_of(src_dir, ncbi, references, inside))
    except Exception as e:  # noqa: BLE001
        say(f"WARNING: the benchmarks failed ({type(e).__name__}: {e}): the package is written without them")
        shutil.rmtree(os.path.join(part, benchmarks.FOLDER), ignore_errors=True)
        return None


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


def main(argv=None):
    p = argparse.ArgumentParser(description=f"Eukore {__version__}. " + __doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-V", action="version", version=f"Eukore {__version__}", help="Eukore's own version, and exit")
    p.add_argument("--work", required=True, help="the build's work folder (sources, calls, extractions)")
    p.add_argument("--out", required=True, help="where Eukore_<version>/ and genomes/ go")
    p.add_argument("--stage", choices=("all", "fetch", "process", "package"), default="all")
    p.add_argument("--version", default=datetime.date.today().strftime("%Y.%m"), help="the package's version (YYYY.MM)")
    p.add_argument("--groups", default=",".join(sources.ASSEMBLY_GROUPS),
                   help=f"NCBI's assembly groups whose species the package holds (default all: "
                        f"{','.join(sources.ASSEMBLY_GROUPS)}; e.g. fungi,protozoa first)")
    p.add_argument("--keep_gb", type=float, default=1000.0,
                   help="the genomes kept in OUT/genomes by the package stage, GB, cut back to it before the package "
                        "is written (so that genomes.tsv says which are kept; default 1000; the build may hold more "
                        "while it works)")
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
                        f"(default the environment's {annotation.VEUPATHDB_KEY}; without, none of its)")
    p.add_argument("--no_rescue", action="store_true", help="no rescue pass of the BUSCOs miniprot missed (loci.py: "
                                                            "their HMMs in the genome's six-frame ORFs place them)")
    p.add_argument("--ksgp_tax", help="a KSGP release's KSGP_v<V>_plus2.tax whose lineages species.tsv gives (protal's "
                                      f"KSGP, e.g. 4.4 or 5.2; default KSGP's served {labels.KSGP_VERSION}, fetched)")
    p.add_argument("--no_benchmarks", action="store_true", help="write the package without its benchmarks/ "
                                                                  "(eukore/benchmarks.py; docs/benchmarks.md)")
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("-j", "--jobs", type=int, default=None, help="genomes at a time (default threads / 4)")
    p.add_argument("--downloads", type=int, default=8,
                   help="the fetch stage's genomes and annotations downloaded at a time (default 8), by curl where on "
                        "the PATH (one NCBI connection: ~30-40 MB/s by curl, ~10 by Python's urllib)")
    p.add_argument("--scratch", help="the fetch stage's unfinished downloads (.part files) and their gzip checks, e.g. "
                                     "a node's local disk or a scratch file system without backups; each file is "
                                     "copied into OUT/genomes once whole (default: beside it)")
    p.add_argument("--report_minutes", type=float, default=util.REPORT_EVERY / 60,
                   help="minutes between progress reports of the genomes' downloads and runs: genomes left, time left "
                        "(default 10)")
    opts = p.parse_args(argv)
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
        fetch_sources(src_dir, opts.groups, opts.ksgp_tax, opts.refresh)
    ncbi = ncbi_mod.Ncbi(os.path.join(src_dir, "taxdump.tar.gz"))
    wanted = opts.wanted = clades(opts.species)
    opts.inside = inside_of(ncbi, wanted)
    job_map, assemblies = plan(opts, src_dir, ncbi)
    # the genomes and their annotations fetched (--downloads at a time, through --scratch), by `all` too (until
    # 2026-11-06 `all` left them to the process stage's workers, --scratch and --downloads unused)
    if opts.stage in ("all", "fetch"):
        loci = loci_of(opts, src_dir)
        # a done genome no longer here is fetched again where better gene models than its run took are there
        # (genomes.better_annotation): the process stage runs it again with them
        again = {a for a, j in job_map.items() if j["url"] and not j.get("path")
                 and not os.path.isfile(genomes.genome_file(opts.genomes, a)) and genomes.is_done(opts.work, a, loci)
                 and genomes.better_annotation(opts.work, j)}
        if again:
            say(f"{len(again):,} done genomes no longer here fetched again for better gene models (e.g. "
                f"{sorted(again)[:3]})")
        todo = [j for a, j in sorted(job_map.items()) if j["url"] and not j.get("path")
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
                                 early=[s for s in util.REPORT_EARLY if s < opts.report_minutes * 60])
        fetched = 0
        # each download's gzip check in a process of its own (a core apart from the downloads' threads; spawned, as
        # a process with threads is not forked), its .part file in --scratch
        checker = concurrent.futures.ProcessPoolExecutor(1, mp_context=multiprocessing.get_context("spawn"))
        how = {"scratch": opts.scratch, "checker": checker}
        if opts.scratch:
            os.makedirs(opts.scratch, exist_ok=True)
        for leftover in glob.glob(os.path.join(opts.genomes, "*.part")) + (  # killed runs' downloads
                glob.glob(os.path.join(opts.scratch, "*.part")) if opts.scratch else []):
            os.remove(leftover)

        def fetch_genome(j):
            path, times = genomes.genome_file(opts.genomes, j["accession"]), collections.Counter()
            why = genomes.fetch_genome(j, path, times=times, **how)
            return j, why, 0 if why else os.path.getsize(path), times

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
                and not annotation.found(opts.genomes, a)[0]
                and (not os.path.isfile(os.path.join(genomes.done_dir(opts.work, a), annotation.TABLE))
                     or genomes.better_annotation(opts.work, j))]
        say(f"fetching {len(todo):,} genomes' annotations (the first of {', '.join(annotation.SOURCES)} that has one) "
            f"into {opts.genomes}")
        progress = util.Progress("annotations fetched", [1] * len(todo), say, opts.report_minutes * 60, "", 1)
        def get(j):
            return j, annotation.fetch(j, opts.genomes, opts.veupathdb_key, **how)
        for done in at_a_time(get, todo, downloads):
            try:
                j, (_path, source) = done.result()
                # the sources asked before the one found (all where none had one): a done genome's not asked again
                genomes.note_annotation_none(opts.work, j, annotation.none_before(j, source))
                progress.step(1)
            except OSError as e:
                progress.step(1, failed=True, note=f" ({e})")
        checker.shutdown()
        if opts.stage == "fetch":
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
    if opts.stage in ("all", "process"):
        settings = settings_of(opts, src_dir, ncbi, staged_only=opts.stage == "process")
        waiting, failed = run_genomes(opts, job_map, settings)
        describe_pass(opts, job_map)
        drop_superseded(opts, job_map)
        organelle_genes(opts, src_dir, job_map, ncbi)
        annotate_pass(opts, job_map, may_fetch=opts.stage == "all")
        # a rerun mends a genome not fetched yet or a download that may come; not one NCBI no longer serves, nor a
        # call that failed (it fails again): `all` packages without those (until 2026-11-06 it waited for them forever)
        retry = waiting + [(a, w) for a, w in failed if w.startswith("fetch failed: ") and not placement.given_up(
            {"status": "failed: " + w[len("fetch failed: "):]})]
        if failed:
            say(f"{len(failed)} genomes failed (e.g. {failed[:3]})")
        if retry:
            say(f"{len(waiting)} genomes not fetched yet, {len(retry) - len(waiting)} downloads to try again (e.g. "
                f"{retry[:3]}): run again (--stage fetch first on a cluster), or --stage package to write the "
                "package without them")
            return 0 if opts.stage == "process" else 1
    stages = []  # the package stage's steps and their seconds (benchmarks' stages.tsv)
    if opts.stage in ("all", "package") and opts.fcs_gx_db:
        todo = [a for a in sorted(job_map) if genomes.is_done(opts.work, a)
                and not os.path.isfile(os.path.join(genomes.done_dir(opts.work, a), "fcs.tsv"))]
        say(f"FCS-GX: {len(todo)} genomes' windows to screen ({opts.fcs_gx_db})")
        failed = 0
        with timed(stages, "FCS-GX"):
            for i, acc in enumerate(todo, 1):
                got = genomes.fcs_screen(opts.work, acc, job_map[acc]["taxid"], opts.fcs_gx_db, opts.run_gx,
                                         opts.threads)
                failed += got is None
                if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                    say(f"  {i}/{len(todo)} screened ({failed} failed)")
    if opts.stage in ("all", "package"):
        vsearch = shutil.which("vsearch")
        if not vsearch:
            say("WARNING: no vsearch: no checks across species (own_species, same_as, the rRNA part screen), every "
                "record's own_species undecided")
        for leftover in glob.glob(os.path.join(opts.out, f".*.Eukore_{opts.version}.partial")):  # killed runs'
            shutil.rmtree(leftover, ignore_errors=True)
        loci = loci_of(opts, src_dir)
        short = [a for a in job_map if genomes.is_done(opts.work, a) and not genomes.is_done(opts.work, a, loci)]
        if short:
            say(f"WARNING: {len(short)} genomes were not searched for all {len(loci)} BUSCOs called (a run before "
                f"format 4, e.g. {short[:3]}): run --stage process again for them all")
        undone = sorted(a for a in job_map if not genomes.is_done(opts.work, a))
        if undone:
            say(f"WARNING: {len(undone)} genomes of the plan are not done (e.g. {undone[:3]}): the package is written "
                "without them (manifest: inputs.genomes_not_done)")
        with timed(stages, "BUSCO copies and their checks"):
            collected = package.collect(ncbi, job_map, assemblies, opts.work, vsearch, opts.threads)
        with timed(stages, "rRNA records and their checks"):
            records = genome_rrna_records(opts, job_map, ncbi) + published(opts, src_dir, ncbi)
            records = rrna.finish(records, ncbi, vsearch, opts.threads, log=say)
        with timed(stages, "organelle genes"):
            chosen = organelle_genes(opts, src_dir, job_map, ncbi)  # before the prune: the genomes searched while there
            held = {int(g["species_taxid"]) for g in collected["genomes"].values() if g["representative"] == "t"}
            organelle = organelle_mod.collect(ncbi, job_map, held | {int(r["species_taxid"]) for _n, _s, r in records},
                                              chosen, os.path.join(opts.work, "organelle"))
        inputs = {"odb": odb_mod.state_of(src_dir), "taxdump_md5": util.md5(os.path.join(src_dir, "taxdump.tar.gz")),
                  "sources": {n: sources.state_of(src_dir, n).get("version") for n in ("EUKARYOME", "ROD")},
                  "assembly_tables": {os.path.basename(p): util.md5(p) for p in sorted(
                      tables(src_dir, "genbank", opts.groups) + tables(src_dir, "refseq", opts.groups))},
                  "genomes_not_done": len(undone),
                  "groups": opts.groups, "no_assembly": util.md5(NO_ASSEMBLY) if os.path.isfile(NO_ASSEMBLY) else None,
                  "ncbi_quality": {os.path.basename(p): util.md5(p) for p in quality.paths(src_dir)},
                  "refseq_organelle": {os.path.basename(p): util.md5(p)
                                       for p in organelle_mod.release_files(src_dir) if os.path.isfile(p)}}
        if wanted is not None:
            inputs["partial"] = {"species_clades": sorted(wanted)}  # a trial package: these clades only
        # every genome is called and extracted: the kept ones cut back to the budget, then the package names those kept
        firsts = {a for r in species.read_table() for a in species.accessions(r)}
        annotate_pass(opts, job_map, may_fetch=opts.stage == "all")
        # the read benchmark's genomes (reads.py), chosen while all are on disk and kept through the prune
        reps = {a: (int(g["species_taxid"]), genomes.genome_file(opts.genomes, a))
                for a, g in collected["genomes"].items() if g["representative"] == "t" and not job_map[a].get("path")}
        tests = reads_mod.choose(ncbi, reps, opts.read_genomes) if opts.read_genomes and not opts.no_benchmarks else []
        with timed(stages, "assembly metadata"):
            describe_pass(opts, job_map)  # before the prune: measured while the genomes are there
        drop_superseded(opts, job_map)
        rows_kept = genomes.prune(opts.genomes, opts.work, opts.keep_gb * 1e9, set(tests), firsts)
        kept = sum(s for _a, s, st in rows_kept if st != "deleted")
        say(f"genomes kept in {opts.genomes}: {kept / 1e9:.1f} GB (--keep_gb {opts.keep_gb:g}), "
            f"{sum(st == 'deleted' for *_x, st in rows_kept)} deleted")
        with timed(stages, "reference labels"):
            references = labels.references(src_dir, sources.ready(src_dir, "EUKARYOME"), opts.ksgp_tax, say)
        hook = None if opts.no_benchmarks else functools.partial(
            measure, opts=opts, src_dir=src_dir, ncbi=ncbi, collected=collected, vsearch=vsearch, stages=stages,
            fates=genomes.fates(ncbi, assemblies, job_map, opts.inside if wanted is not None else None), tests=tests,
            references=references)
        meta = functools.partial(metadata_mod.write, work=opts.work, tables=tables(src_dir, "genbank", opts.groups)
                                 + tables(src_dir, "refseq", opts.groups), jobs=job_map)
        out = package.write(os.path.join(opts.out, f"Eukore_{opts.version}"), opts.version, ncbi, collected,
                            loci_items(opts, src_dir), records, opts.genomes, inputs, organelle, references, hook,
                            metadata=meta)
        say(f"Eukore: {out} ({json.dumps(package.read(out)['counts'])})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
