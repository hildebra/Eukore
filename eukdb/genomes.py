"""EUKDB's genome set (since 2026-11-05): which eukaryote genomes the build fetches, one genome's run (fetched, its euk30
loci placed, its rRNA extracted, its results kept in the work folder), and the prune of the kept genomes to a budget.

A genome is taken for one or more reasons (`why`, comma-separated):
- eukdetect1: an EukDetect1 genome of a species EukDetect2 lacks; its coding copies placed with its own proteins;
- model: a model or host species (data/model_species.tsv); the loci its sources lack, from its relatives' proteins;
- refseq: a RefSeq reference or representative species EUKDB's sources lack or hold on fewer than species.GOOD_LOCI
  loci; as the model species;
- rrna: the species' representative assembly, for its rRNA (every species of the package has one where GenBank does).
Every genome gives both: its placed loci, where it has proteins to place, and its rRNA. Genomes over CHUNK_BASES go to
miniprot in pieces (placement.chunks). A genome is done when WORK/done/<accession>/done.json is written (last); a run
places and extracts only the genomes not done, so that one killed goes on where it stopped.
"""

import collections
import json
import os
import re
import shutil

from . import placement, species, util
from . import sources as src

PLACED_STATUS = ("genomic", "poor placement")


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def genome_file(folder, accession):
    return os.path.join(folder, f"{accession}.fna.gz")


def jobs(ncbi, sources, loci, locus_map, assemblies, extra_rows, placed_dirs=(), with_rrna=True):
    """The genome set -> {accession: job}: job = dict(accession, url, size, species_taxid, why [list], queries {query
    record: coding sequence}, own {query record: True where it is the genome's own coding sequence}, rdna_unit,
    name (the placements' genome name)). sources: {"eukdetect2": folder, "eukdetect1": folder or None}; extra_rows:
    the model and RefSeq species (species table rows, the model species first, each with a why: model or refseq);
    placed_dirs: folders of placed copies
    (data/placed) whose EukDetect1 genomes are placed already (their rRNA still taken)."""
    eukdetect, eukdetect1 = sources["eukdetect2"], sources.get("eukdetect1")
    out = {}

    def add(row, why, name, queries=None, own=False, unit="-"):
        url = placement.genome_url(row) if row else None
        if not url:
            return
        job = out.setdefault(row["accession"], dict(accession=row["accession"], url=url, size=row.get("size"),
                                                    species_taxid=row.get("species_taxid"), why=[], queries={}, own={},
                                                    rdna_unit="-", name=name))
        if why not in job["why"]:
            job["why"].append(why)
        for q, cds in (queries or {}).items():
            job["queries"][q] = cds
            job["own"][q] = own
        if unit and unit != "-":
            job["rdna_unit"] = unit

    placed_before = set()
    for folder in placed_dirs:
        for f in read_rows(os.path.join(folder, src.PLACED_TABLE)):
            if placement.final(f[5]):
                placed_before.add(f[1])
    # EukDetect1's genomes of species EukDetect2 lacks: their own coding sequences
    if eukdetect1 and locus_map:
        wanted = placement.candidates(eukdetect1, eukdetect, locus_map, ncbi)
        _g1, records1 = src.read_link1(os.path.join(eukdetect1, "busco_taxid_link.txt"), locus_map)
        matched = src.match_genomes({g: (t, grp) for g, (t, grp, _c) in wanted.items()}, assemblies, ncbi)
        for genome, (_taxid, _group, cds) in sorted(wanted.items()):
            row = matched[genome][0]
            # queries named by the euk30 locus (EukDetect1's records name its own loci)
            queries = {f"ed1-{genome}-{records1[r][1]}-S1": c for r, c in cds.items()}
            add(row, "eukdetect1", genome, {} if genome in placed_before else queries, own=True)
    # the model and RefSeq species: the loci their sources lack, from their relatives' proteins
    have, (genomes2, records2) = src.species_loci(eukdetect, eukdetect1, locus_map, ncbi, set(loci))
    lacking = {}
    for row in extra_rows:
        missing = set(loci) - have[ncbi.species_of(int(row["ncbi_taxid"]))]
        if missing and species.genome_name(row) not in placed_before:
            lacking[species.genome_name(row)] = (row, missing)
    if lacking:
        from . import exons  # noqa: PLC0415
        relatives, _link, _copies = exons.relatives_of(
            ncbi, eukdetect, eukdetect1, loci, locus_map,
            [os.path.join(d, src.PLACED_FASTA) for d in placed_dirs], link2=(genomes2, records2))
        for name, (row, missing) in sorted(lacking.items()):
            queries = relatives.queries(name, int(row["ncbi_taxid"]), dict.fromkeys(missing),
                                        species.record(row, "x").rsplit("-", 2)[0])
            add(species.assembly_row(row), row.get("why", "model"), name, queries, unit=row.get("rdna_unit", "-"))
    for row in extra_rows:  # placed before or nothing lacking: their genome for its rRNA
        add(species.assembly_row(row), "rrna", species.genome_name(row), unit=row.get("rdna_unit", "-"))
    # every EukDetect2 species' representative assembly, for its rRNA
    if with_rrna:
        copies = src.read_loci(os.path.join(eukdetect, "eukdb.fasta"), records2, set(loci), lengths_only=True)
        species_of = {g: ncbi.species_of(t) for g, (t, _grp) in genomes2.items() if g in copies}
        rep = {}
        for genome in sorted(species_of, key=lambda g: (-len(copies[g]), g)):
            rep.setdefault(species_of[genome], genome)
        matched = src.match_genomes({g: genomes2[g] for g in rep.values()}, assemblies, ncbi)
        for genome in sorted(rep.values()):
            row = matched[genome][0]
            if row and row.get("category") != "derived from metagenome":
                add(row, "rrna", genome)
    return out


def read_rows(path):
    """A placements table (plain or zstd; header first) -> its rows, split."""
    if not path or not os.path.isfile(path):
        return []
    with util.open_text(path) as fh:
        fh.readline()
        return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


RESULT_FILES = ("copies.fna", "placements.tsv", "features.tsv")
PLACEMENT_HEADER = "query\tlocus\tstatus\tcontig\tstart\tend\tstrand\texons\tcds_length\tlength\tidentity\tcoverage\t" \
                   "kmer_share\tmasked\tprotein"


def done_dir(work, accession):
    return os.path.join(work, "done", safe(accession))


def is_done(work, accession):
    return os.path.isfile(os.path.join(done_dir(work, accession), "done.json"))


def run(args):
    """One genome (a job of jobs(), with settings): fetched into settings["genomes"] unless there (or settings
    "staged_only": not fetched, left for a later run), its loci placed, its rRNA extracted (settings["rrna"]: a
    function of (genome path, job, out folder) -> its summary, rrna.extract_genome), its results written into
    WORK/done/<accession>/ (copies.fna, placements.tsv, features.tsv, then done.json) -> (accession, summary or None,
    why it is not done or None)."""
    job, settings = args
    acc = job["accession"]
    folder = done_dir(settings["work"], acc)
    path = genome_file(settings["genomes"], acc)
    if not os.path.isfile(path):
        if settings.get("staged_only"):
            return acc, None, "not fetched yet"
        os.makedirs(settings["genomes"], exist_ok=True)
        why = util.fetch_file(job["url"], path)
        if why:
            return acc, None, f"fetch failed: {why}"
    results = {}
    if job["queries"]:
        own = {q: c for q, c in job["queries"].items() if job["own"].get(q)}
        rel = {q: c for q, c in job["queries"].items() if not job["own"].get(q)}
        for queries, model in ((own, False), (rel, True)):
            if not queries:
                continue
            got = placement.place(queries, path, settings["miniprot"], settings["threads"], placement.MASK_MIN,
                                  placement.UNCAPPED, model=model)
            if any(f[0].startswith("miniprot failed") for _c, f, _x in got.values()):
                return acc, None, next(f[0] for _c, f, _x in got.values() if f[0].startswith("miniprot failed"))
            if model:
                got = {r: (c, f[:12], x, f[12] if len(f) > 12 else "-") for r, (c, f, x) in species.best(got).items()}
            else:
                got = {r: (c, f, x, job["name"]) for r, (c, f, x) in got.items()}
            results.update(got)
    tmp = util.part_of(folder)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    taken = 0
    with open(os.path.join(tmp, "copies.fna"), "w", newline="\n") as fa, \
            open(os.path.join(tmp, "placements.tsv"), "w", newline="\n") as tb, \
            open(os.path.join(tmp, "features.tsv"), "w", newline="\n") as ft:
        tb.write(PLACEMENT_HEADER + "\n")
        ft.write("locus\tkind\tstart\tend\tgenome_contig\tgenome_start\tgenome_end\tstrand\tmethod\tprotein\t"
                 "identity\tcoverage\n")
        for record in sorted(results):
            copy, fields, features, protein = results[record]
            m = src.RECORD_RE.match(record.split("~", 1)[0])
            locus = m.group("locus") if m else record
            tb.write("\t".join(map(str, [record, locus] + list(fields) + [protein])) + "\n")
            if copy is None or fields[0] != PLACED_STATUS[0]:
                continue
            taken += 1
            fa.write(f">{locus}\n{copy}\n")
            method = "placed_own" if job["own"].get(record) else "placed_relative"
            for r in features or ():
                if r[3] in ("exon", "intron"):
                    ft.write("\t".join([locus, r[3], r[4], r[5], r[7], r[8], r[9], r[10], method, protein,
                                        str(fields[8]), str(fields[9])]) + "\n")
    rrna = None
    if settings.get("rrna"):
        try:
            rrna = settings["rrna"](path, job, tmp)
        except (OSError, ValueError, EOFError) as e:
            shutil.rmtree(tmp, ignore_errors=True)
            return acc, None, f"rRNA extraction failed: {e}"
    summary = {"accession": acc, "why": job["why"], "copies_taken": taken, "rrna": rrna,
               "bytes": os.path.getsize(path)}
    with open(os.path.join(tmp, "done.json"), "w", newline="\n") as fh:
        json.dump(summary, fh, indent=1)
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    os.replace(tmp, folder)
    return acc, summary, None


def prune(genomes_folder, work, keep_bytes, protected=(), first=()):
    """The kept genomes cut back to keep_bytes once done (since 2026-11-05): `protected` accessions (a host's) kept
    whatever their size, then `first` (the model species'), then the smallest; a genome not done yet kept (pending);
    the rest deleted -> [(accession, bytes, status)], written as genomes_folder/manifest.tsv."""
    files = [f for f in os.listdir(genomes_folder) if f.endswith(".fna.gz")] if os.path.isdir(genomes_folder) else []
    order = []
    for f in files:
        acc = f[:-len(".fna.gz")]
        rank = 0 if acc in protected else 1 if acc in first else 2
        order.append((rank, os.path.getsize(os.path.join(genomes_folder, f)), acc))
    kept, rows = 0, []
    for rank, size, acc in sorted(order):
        if rank == 0 or not is_done(work, acc) or kept + size <= keep_bytes:
            kept += size
            rows.append((acc, size, "kept" if is_done(work, acc) or rank == 0 else "pending"))
        else:
            os.remove(genome_file(genomes_folder, acc))
            rows.append((acc, size, "deleted"))
    path = os.path.join(genomes_folder, "manifest.tsv")
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("accession\tbytes\tstatus\n" + "".join(f"{a}\t{s}\t{st}\n" for a, s, st in rows))
    os.replace(util.part_of(path), path)
    return rows


def why_counts(job_list):
    return collections.Counter(w for j in job_list for w in j["why"])
