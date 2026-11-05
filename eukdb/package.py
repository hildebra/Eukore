"""EUKDB's package (FORMAT.md, since 2026-11-05): written from the build's sources and work, read and checked by its
users (protal's converter).
"""

import collections
import datetime
import gzip
import json
import os
import shutil

from . import copies as copies_mod
from . import sources as src
from . import species as species_mod
from . import util

FORMAT = 1
MANIFEST = "manifest.json"
GENOME_COLUMNS = ("accession", "species_taxid", "taxid", "organism", "strain", "level", "genome_size", "contigs",
                  "category", "ftp", "source", "representative", "loci", "local")
FEATURE_COLUMNS = ("accession", "locus", "kind", "start", "end", "genome_contig", "genome_start", "genome_end", "strand",
                   "method", "protein", "identity", "coverage")
SPECIES_COLUMNS = ("taxid", "name", "representative", "loci", "sources", "rrna_records")
SOURCE_ORDER = {"eukdetect2": 0, "eukdetect1": 1, "model": 2, "refseq": 3, "assembly": 4}


def _tsv(path, columns, rows):
    with open(util.part_of(path), "w", newline="\n") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "-")) if r.get(c, "-") not in (None, "") else "-" for c in columns) + "\n")
    os.replace(util.part_of(path), path)


def collect(ncbi, sources, assemblies, loci, locus_map, extra_rows, work, placed_dirs=(), exon_tables=()):
    """What the package holds of euk30 -> dict(genomes {accession: row}, copies {accession: {locus: sequence}},
    features [FEATURE_COLUMNS dicts]). EukDetect2's genomes with their copies whole; EukDetect1's genomes of species
    EukDetect2 lacks, each locus its placed copy where one was taken (WORK/done, then placed_dirs), else its coding one;
    the model and RefSeq species' placed copies (WORK/done), as further genomes of their species; their exons from
    exon_tables (exons.py's, EukDetect2's copies), the placements' features and, for a coding copy, one exon."""
    eukdetect, eukdetect1 = sources["eukdetect2"], sources.get("eukdetect1")
    loci = list(loci)
    genomes2, records2 = src.read_link(os.path.join(eukdetect, "busco_taxid_genome_link.txt"))
    copies2 = src.read_loci(os.path.join(eukdetect, "eukdb.fasta"), records2, set(loci))
    genomes, copies, features = {}, {}, []
    coding = set()  # (EukDetect1 genome, locus) whose copy is its coding sequence
    taxid_of, source_of, row_of = {}, {}, {}
    for genome, by_locus in copies2.items():
        taxid_of[genome], source_of[genome] = genomes2[genome][0], "eukdetect2"
        copies[genome] = {l: s for l, s in by_locus.items() if s.strip("N")}
    # EukDetect1's genomes of species EukDetect2 lacks
    known = {ncbi.species_of(t) for t, _g in genomes2.values()}
    placed = collections.defaultdict(dict)  # EukDetect1 genome -> {euk30 locus: placed copy}
    record_of = {}
    if eukdetect1 and locus_map:
        genomes1, records1 = src.read_link1(os.path.join(eukdetect1, "busco_taxid_link.txt"), locus_map)
        copies1 = src.read_loci(os.path.join(eukdetect1, "ncbi_eukprot_met_arch_markers.fna"), records1)
        record_of = {gl: r for r, gl in records1.items()}
        for folder in placed_dirs:
            path = os.path.join(folder, src.PLACED_FASTA)
            if os.path.isfile(path):
                placed.update(src.read_loci(path, records1))
        for genome, by_locus in sorted(copies1.items()):
            taxid, _group = genomes1[genome]
            sp = ncbi.species_of(taxid)
            if sp in known or genome in genomes2 or ncbi.rank.get(sp) != "species" or ncbi.is_placeholder(sp):
                continue
            taxid_of[genome], source_of[genome] = taxid, "eukdetect1"
            copies[genome] = {}
            for locus, seq in by_locus.items():
                copies[genome][locus] = placed[genome].get(locus) or seq
                if locus not in placed[genome]:
                    coding.add((genome, locus))
    # the model and RefSeq species' placed copies, as further genomes of their species
    for row in extra_rows:
        name = species_mod.genome_name(row)
        got = done_copies(work, row["genbank"])
        if got:
            taxid_of[name], source_of[name] = int(row["ncbi_taxid"]), row.get("why", "model")
            row_of[name] = species_mod.assembly_row(row)
            copies[name] = got
    # their assemblies: EukDetect's genomes by name (sources.match_genomes), each one's own
    rep = {}
    by_species = collections.defaultdict(list)
    for g in copies:
        by_species[ncbi.species_of(taxid_of[g])].append(g)
    for sp, gs in by_species.items():
        rep[sp] = min(gs, key=lambda g: (-len(copies[g]), SOURCE_ORDER[source_of[g]], g))
    eukdetect_genomes = {g: (taxid_of[g], "-") for g in copies if g not in row_of}
    matched = src.own_accessions(src.match_genomes(eukdetect_genomes, assemblies, ncbi),
                                 {g: ncbi.species_of(taxid_of[g]) for g in eukdetect_genomes}, rep)
    accession_of = {}
    for g in copies:
        row = row_of.get(g) or matched.get(g, (None,))[0]
        acc = row["accession"] if row else "EUK_" + g.replace("/", "_")
        accession_of[g] = acc
        sp = ncbi.species_of(taxid_of[g])
        genomes[acc] = dict(accession=acc, species_taxid=sp, taxid=row["taxid"] if row else taxid_of[g],
                            organism=row["organism"] if row else ncbi.name.get(taxid_of[g], g),
                            strain=(row or {}).get("strain") or (row or {}).get("isolate") or "-",
                            level=(row or {}).get("level"), genome_size=(row or {}).get("size"),
                            contigs=(row or {}).get("contigs"), category=(row or {}).get("category"),
                            ftp=(row or {}).get("ftp"), source=source_of[g],
                            representative="t" if rep[sp] == g else "f", loci=len(copies[g]), name=g)
    # EukDetect1's copies this build's run placed (WORK/done of its assembly), in place of the coding ones
    for g in copies:
        if source_of[g] == "eukdetect1":
            for locus, seq in done_copies(work, accession_of[g]).items():
                copies[g][locus] = seq
                coding.discard((g, locus))
    for g, locus in sorted(coding):  # a coding copy: one exon over it
        features.append(dict(accession=g, locus=locus, kind="exon", start=0, end=len(copies[g][locus]),
                             method="coding", protein=g))
    # further GenBank assemblies of EukDetect's species (the model and RefSeq species: their reference alone)
    eukdetect_species = {ncbi.species_of(taxid_of[g]) for g in copies if g not in row_of}
    for row in sorted(assemblies, key=src.assembly_rank):
        if row["species_taxid"] in eukdetect_species and row["accession"] not in genomes \
                and not ncbi.is_placeholder(row["species_taxid"]):
            genomes[row["accession"]] = dict(accession=row["accession"], species_taxid=row["species_taxid"],
                                             taxid=row["taxid"], organism=row["organism"],
                                             strain=row["strain"] or row["isolate"] or "-", level=row["level"],
                                             genome_size=row["size"], contigs=row["contigs"],
                                             category=row["category"], ftp=row["ftp"], source="assembly",
                                             representative="f", loci=0)
    # the features: by accession
    features = [dict(f, accession=accession_of[f["accession"]]) for f in features]
    for g in copies:
        for f in done_features(work, row_of[g]["accession"] if g in row_of else accession_of[g]):
            features.append(dict(f, accession=accession_of[g]))
    by_genome = {g: accession_of[g] for g in copies}
    for table in exon_tables:
        for r in copies_mod.read_features(table):
            genome, _, locus = r[0].partition("~")
            if genome in by_genome and r[3] in ("exon", "intron") and source_of.get(genome) == "eukdetect2":
                features.append(dict(accession=by_genome[genome], locus=locus, kind=r[3], start=r[4], end=r[5],
                                     strand=r[10], method="annotated", protein="-"))
    return dict(genomes=genomes, copies={accession_of[g]: c for g, c in copies.items()}, features=features,
                taxids={taxid_of[g] for g in copies} | {r["taxid"] for r in genomes.values() if r.get("taxid")})


def done_copies(work, accession):
    """A genome's copies taken in this build's run (WORK/done/<accession>/copies.fna) -> {locus: copy}."""
    path = os.path.join(work or "", "done", accession.replace("/", "_"), "copies.fna")
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        lines = [line.strip() for line in fh]
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines) - 1, 2)}


def done_features(work, accession):
    path = os.path.join(work or "", "done", accession.replace("/", "_"), "features.tsv")
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def write(out, version, ncbi, collected, loci_table, rrna=None, genomes_folder=None, inputs=None):
    """The package into `out` (written beside it and renamed into place; the manifest last) -> its folder.
    rrna: (records [(name, sequence, row dict of records.tsv)]) or None (no rRNA half yet)."""
    part = util.part_of(out)
    shutil.rmtree(part, ignore_errors=True)
    for sub in ("taxonomy", "euk30/copies", "rrna"):
        os.makedirs(os.path.join(part, sub))
    genomes = collected["genomes"]
    local = set()
    if genomes_folder and os.path.isdir(genomes_folder):
        local = {f[:-len(".fna.gz")] for f in os.listdir(genomes_folder) if f.endswith(".fna.gz")}
    for acc, r in genomes.items():
        r["local"] = f"genomes/{acc}.fna.gz" if acc in local else "-"
    _tsv(os.path.join(part, "genomes.tsv"), GENOME_COLUMNS,
         sorted(genomes.values(), key=lambda r: (str(r["species_taxid"]), r["representative"] != "t", r["accession"])))
    loci = src.read_loci_table(loci_table)
    shutil.copyfile(loci_table, os.path.join(part, "euk30", "loci.tsv"))
    copies = collected["copies"]
    for locus in loci:
        path = os.path.join(part, "euk30", "copies", f"{locus}.fna.gz")
        with gzip.open(path, "wt", compresslevel=6) as fh:
            for acc in sorted(copies):
                if copies[acc].get(locus):
                    fh.write(f">{acc}\n{copies[acc][locus].upper()}\n")
    _tsv(os.path.join(part, "euk30", "features.tsv"), FEATURE_COLUMNS,
         sorted(collected["features"], key=lambda f: (f["accession"], f["locus"], int(f["start"]))))
    n_records = collections.Counter()
    if rrna:
        from . import rrna as rrna_mod  # noqa: PLC0415
        rrna_mod.write_records(os.path.join(part, "rrna"), rrna)
        for _name, _seq, row in rrna:
            n_records[int(row["species_taxid"])] += 1
    by_species = collections.defaultdict(list)
    for r in genomes.values():
        by_species[int(r["species_taxid"])].append(r)
    rows = []
    for sp, rs in sorted(by_species.items()):
        rep = next((r for r in rs if r["representative"] == "t"), None)
        if rep is None:
            continue
        rows.append(dict(taxid=sp, name=ncbi.name.get(sp, rep["organism"]), representative=rep["accession"],
                         loci=rep["loci"], rrna_records=n_records.get(sp, 0),
                         sources=",".join(sorted({r["source"] for r in rs if r["source"] != "assembly"},
                                                 key=SOURCE_ORDER.get))))
    _tsv(os.path.join(part, "species.tsv"), SPECIES_COLUMNS, rows)
    taxa = ncbi.write_subset(os.path.join(part, "taxonomy"),
                             set(collected["taxids"]) | {int(r["species_taxid"]) for r in genomes.values()})
    files = {}
    for dirpath, _dirs, names in os.walk(part):
        for name in sorted(names):
            path = os.path.join(dirpath, name)
            files[os.path.relpath(path, part).replace(os.sep, "/")] = [os.path.getsize(path), util.md5(path)]
    manifest = {"format": FORMAT, "version": version, "built": datetime.datetime.now().isoformat(timespec="seconds"),
                "inputs": inputs or {}, "counts": {"species": len(rows), "genomes": len(genomes),
                                                   "copies": sum(len(c) for c in copies.values()),
                                                   "rrna_records": sum(n_records.values()), "taxa": taxa},
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
    """A package table (genomes.tsv, species.tsv, euk30/features.tsv, rrna/records.tsv) -> [row dicts]."""
    with util.open_text(os.path.join(folder, name)) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def read_copies(folder, locus):
    """A locus's copies -> {accession: sequence}."""
    out, name = {}, None
    with gzip.open(os.path.join(folder, "euk30", "copies", f"{locus}.fna.gz"), "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
            elif name:
                out[name] = line.strip()
    return out
