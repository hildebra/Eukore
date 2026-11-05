"""EUKDB's package (FORMAT.md; format 2 since 2026-11-05: EUKDB's own locus calls in NCBI's genomes; format 4 since
2026-11-05: every BUSCO of the dataset, in busco/, with busco/loci.tsv), written from the build's genome set and work
folder, read and checked by its users (protal's converter).
"""

import collections
import datetime
import gzip
import json
import os
import shutil

from . import qc
from . import util

FORMAT = 4
MANIFEST = "manifest.json"
LOCI_COLUMNS = ("locus", "name", "cutoff", "hmm_length", "species", "called", "duplicated", "fragmented", "call_rate",
                "median_length", "median_cds_length", "median_protein_length")
GENOME_COLUMNS = ("accession", "fetched", "species_taxid", "taxid", "organism", "strain", "level", "genome_size",
                  "contigs", "category", "flags", "ftp", "source", "representative", "loci", "local")
COPY_COLUMNS = ("accession", "locus", "status", "contig", "contig_length", "start", "end", "strand", "length",
                "cds_length", "exons", "identity", "coverage", "hmm_score", "cutoff", "hmm_coverage", "copies",
                "query", "protein_length", "window_class", "fcs", "own_species", "same_as", "length_check")
FEATURE_COLUMNS = ("accession", "locus", "kind", "start", "end")
SPECIES_COLUMNS = ("taxid", "name", "representative", "source", "loci", "rrna_records")


def _tsv(path, columns, rows):
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "-")) if r.get(c, "-") not in (None, "") else "-" for c in columns) + "\n")
    os.replace(util.part_of(path), path)


def done_table(work, accession, name):
    """A done genome's table (WORK/done/<accession>/<name>, header first) -> [row dicts]."""
    from .genomes import done_dir  # noqa: PLC0415
    path = os.path.join(done_dir(work, accession), name)
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def done_copies(work, accession):
    """A genome's copies (WORK/done/<accession>/copies.fna) -> {locus: copy}."""
    from .genomes import done_dir  # noqa: PLC0415
    path = os.path.join(done_dir(work, accession), "copies.fna")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        lines = [line.strip() for line in fh]
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 2)}


def collect(ncbi, job_map, assemblies, work, vsearch=None, threads=1):
    """What the package holds of the BUSCO genes -> dict(genomes {accession: row}, copies {accession: {locus: sequence}},
    features [FEATURE_COLUMNS dicts], copy_rows [COPY_COLUMNS dicts], taxids). Each species' genome (its job) with
    the copies its run called; the species' other assemblies listed (source "assembly", no copies); the copies'
    checks across species (qc.check)."""
    by_accession = {r["accession"]: r for r in assemblies}
    genomes, copies, features, copy_rows = {}, {}, [], []
    for acc, job in sorted(job_map.items()):
        got = done_copies(work, acc)
        if not got:
            continue
        row = by_accession.get(acc, {})
        genomes[acc] = dict(accession=acc, fetched=job.get("fetch", acc), species_taxid=job["species_taxid"],
                            taxid=job.get("taxid") or job["species_taxid"],
                            organism=row.get("organism") or job.get("organism", ""),
                            strain=row.get("strain") or row.get("isolate") or "-", level=row.get("level") or job.get(
                                "level", ""), genome_size=row.get("size") or job.get("size"),
                            contigs=row.get("contigs"), category=row.get("category") or job.get("category", ""),
                            flags=",".join(job.get("flags") or ()) or "-", ftp=row.get("ftp"), source=job["why"][0],
                            representative="t", loci=len(got))
        copies[acc] = got
        for f in done_table(work, acc, "features.tsv"):
            features.append(dict(accession=acc, locus=f["locus"], kind=f["kind"], start=int(f["start"]),
                                 end=int(f["end"])))
        from .genomes import done_dir, fcs_of_copy  # noqa: PLC0415
        fcs = [dict(f, start=int(f["start"]), end=int(f["end"])) for f in done_table(work, acc, "fcs.tsv")]
        screened = os.path.isfile(os.path.join(done_dir(work, acc), "fcs.tsv"))
        for r in done_table(work, acc, "loci.tsv"):
            if r["status"] != "Missing":
                copy_rows.append(dict(r, accession=acc, fcs=fcs_of_copy(fcs, r["locus"], r["contig"], r["start"],
                                                                         r["end"]) if screened else "-"))
        for other in job.get("others", ()):
            o = by_accession.get(other)
            if o and other not in genomes:
                genomes[other] = dict(accession=other, fetched=other, species_taxid=job["species_taxid"],
                                      taxid=o["taxid"], organism=o["organism"],
                                      strain=o.get("strain") or o.get("isolate") or "-", level=o["level"],
                                      genome_size=o["size"], contigs=o["contigs"], category=o["category"],
                                      flags=",".join(o.get("flags") or ()) or "-", ftp=o["ftp"], source="assembly",
                                      representative="f", loci=0)
    exons = collections.defaultdict(list)
    for f in features:
        exons[(f["accession"], f["locus"])].append((f["kind"], f["start"], f["end"]))
    cds = {(r["accession"], r["locus"]): qc.spliced(copies[r["accession"]][r["locus"]],
                                                     exons[(r["accession"], r["locus"])]) for r in copy_rows}
    checks = qc.check(cds, {a: g["species_taxid"] for a, g in genomes.items()}, ncbi, vsearch, threads,
                      {(r["accession"], r["locus"]): int(r["protein_length"]) for r in copy_rows})
    for r in copy_rows:
        c = checks[(r["accession"], r["locus"])]
        r.update(own_species=c["own_species"], same_as=c["same_as"], length_check=c["length"])
    taxids = {g["taxid"] for g in genomes.values()} | {int(r["same_as"]) for r in copy_rows if r["same_as"] != "-"}
    return dict(genomes=genomes, copies=copies, features=features, copy_rows=copy_rows, taxids=taxids)


def loci_rows(loci, collected, ncbi):
    """busco/loci.tsv's rows: each BUSCO (loci: dicts of locus, name, cutoff, hmm_length) with how the package's
    species hold it: species with a genome, called (Single or Duplicated), duplicated, fragmented, the call rate over
    all and per supergroup (groups.GROUPS), the median copy, coding and protein lengths."""
    from . import groups  # noqa: PLC0415
    genomes = collected["genomes"]
    reps = {a: g for a, g in genomes.items() if g.get("representative") == "t"}
    group = {a: groups.group_of(ncbi, g["species_taxid"]) for a, g in reps.items()}
    sizes = collections.Counter(group.values())
    by_locus = collections.defaultdict(list)
    for r in collected["copy_rows"]:
        if r["accession"] in reps:
            by_locus[r["locus"]].append(r)
    out = []
    for item in loci:
        rs = by_locus.get(item["locus"], [])
        called = [r for r in rs if r["status"] in ("Single", "Duplicated")]
        per_group = collections.Counter(group[r["accession"]] for r in called)

        def median(key):
            values = sorted(int(r[key]) for r in called if str(r.get(key, "")).isdigit())
            return values[len(values) // 2] if values else 0
        row = dict(item, species=len(reps), called=len(called),
                   duplicated=sum(r["status"] == "Duplicated" for r in rs),
                   fragmented=sum(r["status"] == "Fragmented" for r in rs),
                   call_rate=f"{len(called) / max(len(reps), 1):.3f}", median_length=median("length"),
                   median_cds_length=median("cds_length"), median_protein_length=median("protein_length"))
        for name in sorted(sizes):
            row[f"call_{name}"] = f"{per_group[name] / sizes[name]:.3f}"
        out.append(row)
    return out, [f"call_{name}" for name in sorted(sizes)]


def write(out, version, ncbi, collected, loci, rrna=None, genomes_folder=None, inputs=None):
    """The package into `out` (written beside it and renamed into place; the manifest last) -> its folder.
    loci: the BUSCOs called ([{locus, name, cutoff, hmm_length}], in the dataset's order); rrna: (records [(name,
    sequence, row dict of records.tsv)]) or None."""
    part = util.part_of(out)
    shutil.rmtree(part, ignore_errors=True)
    for sub in ("taxonomy", "busco/copies", "rrna"):
        os.makedirs(os.path.join(part, sub))
    genomes = collected["genomes"]
    local = set()
    if genomes_folder and os.path.isdir(genomes_folder):
        local = {f[:-len(".fna.gz")] for f in os.listdir(genomes_folder) if f.endswith(".fna.gz")}
    for acc, r in genomes.items():
        r["local"] = f"genomes/{acc}.fna.gz" if acc in local else "-"
    _tsv(os.path.join(part, "genomes.tsv"), GENOME_COLUMNS,
         sorted(genomes.values(), key=lambda r: (str(r["species_taxid"]), r["representative"] != "t", r["accession"])))
    rows_of_loci, call_columns = loci_rows(loci, collected, ncbi)
    _tsv(os.path.join(part, "busco", "loci.tsv"), LOCI_COLUMNS + tuple(call_columns), rows_of_loci)
    order = {item["locus"]: i for i, item in enumerate(loci)}
    copies = collected["copies"]
    for item in loci:
        locus = item["locus"]
        path = os.path.join(part, "busco", "copies", f"{locus}.fna.gz")
        with gzip.open(path, "wt", compresslevel=6) as fh:
            for acc in sorted(copies):
                if copies[acc].get(locus):
                    fh.write(f">{acc}\n{copies[acc][locus].upper()}\n")
    _tsv(os.path.join(part, "busco", "copies.tsv"), COPY_COLUMNS,
         sorted(collected["copy_rows"], key=lambda r: (r["accession"], order.get(r["locus"], 1 << 20))))
    _tsv(os.path.join(part, "busco", "features.tsv"), FEATURE_COLUMNS,
         sorted(collected["features"], key=lambda f: (f["accession"], f["locus"], int(f["start"]), f["kind"])))
    n_records = collections.Counter()
    if rrna:
        from . import rrna as rrna_mod  # noqa: PLC0415
        rrna_mod.write_records(os.path.join(part, "rrna"), rrna)
        for _name, _seq, row in rrna:
            n_records[int(row["species_taxid"])] += 1
    rows = []
    for acc, r in sorted(genomes.items(), key=lambda kv: int(kv[1]["species_taxid"])):
        if r["representative"] != "t":
            continue
        sp = int(r["species_taxid"])
        rows.append(dict(taxid=sp, name=ncbi.name.get(sp, r["organism"]), representative=acc, source=r["source"],
                         loci=r["loci"], rrna_records=n_records.get(sp, 0)))
    held = {int(r["taxid"]) for r in rows}
    rows += [dict(taxid=sp, name=ncbi.name.get(sp, str(sp)), representative="-", source="rrna", loci=0,
                  rrna_records=n) for sp, n in sorted(n_records.items()) if sp not in held]
    _tsv(os.path.join(part, "species.tsv"), SPECIES_COLUMNS, rows)
    taxa = ncbi.write_subset(os.path.join(part, "taxonomy"),
                             set(collected["taxids"]) | {int(r["species_taxid"]) for r in genomes.values()} |
                             set(n_records) | ({int(row["taxid"]) for _n, _s, row in rrna} if rrna else set()))
    files = {}
    for dirpath, _dirs, names in os.walk(part):
        for name in sorted(names):
            path = os.path.join(dirpath, name)
            files[os.path.relpath(path, part).replace(os.sep, "/")] = [os.path.getsize(path), util.md5(path)]
    statuses = collections.Counter(r["status"] for r in collected["copy_rows"])
    manifest = {"format": FORMAT, "version": version, "built": datetime.datetime.now().isoformat(timespec="seconds"),
                "inputs": inputs or {},
                "counts": {"species": sum(1 for r in rows if r["representative"] != "-"),
                           "species_rrna_only": sum(1 for r in rows if r["representative"] == "-"),
                           "genomes": len(genomes), "copies": sum(len(c) for c in copies.values()),
                           "copy_status": dict(statuses), "rrna_records": sum(n_records.values()), "taxa": taxa},
                "files": dict(sorted(files.items()))}
    with open(os.path.join(part, MANIFEST), "w", newline="\n") as fh:
        json.dump(manifest, fh, indent=1)
    old = out + ".old"
    if os.path.isdir(out):
        os.replace(out, old)
    os.replace(part, out)
    shutil.rmtree(old, ignore_errors=True)
    return out


def read(folder):
    """A package's manifest, or None where there is none (unfinished)."""
    try:
        with open(os.path.join(folder, MANIFEST)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def check(folder, deep=False):
    """Why a package is not whole (no manifest, another format, a file missing or of another size; deep: or another
    MD5), or None."""
    manifest = read(folder)
    if manifest is None:
        return f"{folder}: no EUKDB package (no {MANIFEST})"
    if manifest.get("format") != FORMAT:
        return f"{folder}: format {manifest.get('format')}, this reader takes {FORMAT}"
    for rel, (size, digest) in manifest["files"].items():
        path = os.path.join(folder, rel)
        if not os.path.isfile(path):
            return f"{folder}: {rel} is missing"
        if os.path.getsize(path) != size or (deep and util.md5(path) != digest):
            return f"{folder}: {rel} is not the file its manifest names"
    return None


def read_table(folder, name):
    """A package table (genomes.tsv, species.tsv, busco/loci.tsv, busco/copies.tsv, busco/features.tsv,
    rrna/records.tsv) -> [row
    dicts]."""
    with util.open_text(os.path.join(folder, name)) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_copies(folder, locus):
    """A locus's copies -> {accession: sequence}."""
    out, name = {}, None
    with gzip.open(os.path.join(folder, "busco", "copies", f"{locus}.fna.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            elif name:
                out[name] = line.strip()
    return out
