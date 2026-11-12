"""A Eukore package beside the other databases, as collections of records (since 2026-11-12; compare.py runs it;
benchmarks/README.md):

- databases.tsv: the comparison table. Each database's words (benchmarks/databases.tsv: its markers, taxa, sequences,
  curation, updates, reference) with what the benchmarks count: the version read, its NCBI species, those in Eukore
  too, its records and bases;
- upset.tsv: the species by the databases that hold them (every combination, per group and over all): what several
  databases share, and what one holds alone;
- novel.tsv: per group and marker (18S, ITS, operon, mitochondrial genes) the species Eukore holds it of from their
  genome's own assembly, and those of them that no reference of that marker the build read names: the species a
  genome is the first sequence of that marker for (by name: an upper bound, a reference may name it otherwise);
- agreement.tsv: where Eukore and another database hold the same species' marker, whether the two sequences agree:
  the published reads' own-species check (the package's rrna/records.tsv), the reference screen's calls (the
  package's manifest) and another 18S reference's records aligned to Eukore's genomes' 18S (minimap2);
  agreement_records.tsv lists that reference's records that are another species' gene (`disagree`) or unlike the
  species' own genome's (`unlike`): one of the two databases holds another organism's sequence, or a divergent copy.

    python3 -m eukore_bench.databases --package OUT/Eukore_2026.11 --bench BENCH --out RESULTS [-t 8]
"""

import argparse
import collections
import json
import os
import shutil
import subprocess
import sys
import tempfile

from . import ROOT
from . import prep, resource
from eukore import benchmarks, ncbi as ncbi_mod, odb as odb_mod, rrna, util

EUKORE = "Eukore"
TABLE = os.path.join(ROOT, "benchmarks", "databases.tsv")
WORDS = ("markers", "taxa", "sequences", "curation", "updates", "reference")
DATABASE_COLUMNS = ("database",) + WORDS + ("version", "species", "both", "eukore_only", "database_only", "records",
                                           "bases")
UPSET_COLUMNS = ("supergroup", "databases", "n", "eukore", "species")
NOVEL_COLUMNS = ("supergroup", "marker", "references", "genome", "named", "new")
# a marker's references: those whose species names say whether a species' marker is published (species.tsv's
# <reference>_match columns; ROD's species of database_species.tsv). The build reads EUKARYOME's long reads alone
NOVEL = (("18S", "18S", ("PR2", "SILVA", "KSGP", "EUKARYOME")), ("ITS", "ITS", ("UNITE", "EUKARYOME")),
         ("operon", "operon", ("ROD", "EUKARYOME")), ("mitochondrial genes", "mitochondrion", ("MIDORI2",)))
NAMED = ("taxid", "species", "synonym")  # a <reference>_match that names the species itself
AGREEMENT_COLUMNS = ("database", "marker", "records", "agree", "disagree", "unlike", "undecided", "species",
                     "species_agree", "species_differ", "how")
RECORD_COLUMNS = ("database", "record", "verdict", "species_taxid", "name", "own_identity", "other_taxid",
                  "other_name", "other_identity", "shared_rank")
SSU_REFERENCES = ("PR2", "SILVA")  # another reference of BENCH/compare whose records are 18S: aligned to Eukore's
AGREE = 0.97        # a reference's 18S this alike to the species' own genome's: the same species' gene
MARGIN = 0.02       # this alike to another species' genome's, and nearer to it by this: that species' gene
MIN_COVER = 0.8     # of the record, or of Eukore's 18S (a record longer than it: the operon's other parts)
MIN_LENGTH = 500
RANKS = ("genus", "family", "order", "class", "phylum")


# --- the species each database holds -----------------------------------------------------------------------------------

def species_sets(out, species, theirs, ncbi):
    """{database: its NCBI species} of the build's database_species.tsv (PR2's names, EUKARYOME, KSGP, SILVA, UNITE,
    MIDORI2, ROD, BUSCO's dataset), compare.py's other references (theirs: of the same name, theirs) and Eukore's
    own (any record)."""
    sets = collections.defaultdict(set)
    for r in benchmarks.rows_of(os.path.join(out, "database_species.tsv")):
        sets[r["database"]].add(int(r["species_taxid"]))
    for name, held in (theirs or {}).items():
        sets[name] = set(held)
    sets[EUKORE] = {t for t in (ncbi.merged.get(int(x), int(x)) for x in species) if ncbi.rank.get(t) == "species"}
    return dict(sets)


def upset_rows(sets, species, ncbi):
    """upset.tsv's rows: per group (and `all`) and combination of databases (their names joined by " + ", Eukore
    first, then by name) the species exactly those hold."""
    cache = {}
    mine = {ncbi.merged.get(int(t), int(t)): r["supergroup"] for t, r in species.items()}
    holders = collections.defaultdict(list)
    for name in sorted(sets, key=lambda n: (n != EUKORE, n)):
        for sp in sets[name]:
            holders[sp].append(name)
    counts = collections.Counter()
    for sp, names in holders.items():
        group = mine.get(sp) or benchmarks.supergroup_of(ncbi, sp, cache)
        for g in (group, "all"):
            counts[(g, tuple(names))] += 1
    return [dict(supergroup=g, databases=" + ".join(names), n=len(names), eukore="t" if EUKORE in names else "f",
                 species=n)
            for (g, names), n in sorted(counts.items(), key=lambda kv: (kv[0][0] == "all", kv[0][0], -kv[1], kv[0][1]))]


# --- the species a genome is the first sequence of a marker for -------------------------------------------------------

def novel_rows(package, species, sets):
    """novel.tsv's rows: per group (and `all`) and marker the species holding it from their genome's own assembly
    (`genome`), those of them a reference of the marker names (`named`: by its name, a synonym or its taxid) and
    those none names (`new`); `references`: the marker's references the build read (none: `-`, nothing counted)."""
    held = resource.markers_held(package, species)
    read = {ref for ref in benchmarks.REFERENCES
            if any(r.get(f"{ref.lower()}_match", "-") != "-" for r in species.values())}
    group = {t: r["supergroup"] for t, r in species.items()}
    rows = []
    for label, marker, refs in NOVEL:
        used = [ref for ref in refs if ref in read or ref in sets]
        named = set()
        for t, r in species.items():
            if any(r.get(f"{ref.lower()}_match") in NAMED for ref in used if ref in read) or \
                    any(int(t) in sets[ref] for ref in used if ref not in read):
                named.add(t)
        for g in sorted(set(group.values())) + ["all"]:
            mine = {t for t in held[marker][1] if g in (group[t], "all")}
            rows.append(dict(supergroup=g, marker=label, references=", ".join(used) or "-", genome=len(mine),
                             named=len(mine & named) if used else "-", new=len(mine - named) if used else "-"))
    return rows


# --- the same species' marker in two databases -------------------------------------------------------------------------

def differing(by_species):
    """The species ({species: its records' verdicts}) none of whose records agrees and one or more disagrees or is
    unlike: the species whose marker the two databases hold otherwise."""
    return sum("agree" not in v and bool(v & {"disagree", "unlike"}) for v in by_species.values())


def published_agreement(package):
    """agreement.tsv's rows of the package's own checks: the published reads (rrna/records.tsv's EUKARYOME and ROD
    records) by their own-species check against the genomes' records, and the reference screen's calls of the
    MIDORI2 files (the manifest's inputs.references.MIDORI2.screen)."""
    counts = collections.defaultdict(collections.Counter)
    by_species = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in benchmarks.rows_of(os.path.join(package, "rrna", "records.tsv")):
        if r["source"] in ("genome", "curated_unit"):
            continue
        verdict = {"t": "agree", "f": "disagree"}.get(r.get("own_species", "-"), "undecided")
        counts[r["source"]][verdict] += 1
        by_species[r["source"]][r["species_taxid"]].add(verdict)
    rows = []
    for source, c in sorted(counts.items()):
        sp = by_species[source]
        rows.append(dict(database=source, marker="18S and operons", records=sum(c.values()), agree=c["agree"],
                         disagree=c["disagree"], unlike=0, undecided=c["undecided"], species=len(sp),
                         species_agree=sum("agree" in v for v in sp.values()), species_differ=differing(sp),
                         how="the build's own-species check against the genomes' records"))
    try:
        with open(os.path.join(package, "manifest.json")) as fh:
            screen = ((json.load(fh).get("inputs") or {}).get("references") or {}).get("MIDORI2", {}).get("screen")
    except (OSError, ValueError):
        screen = None
    if screen and screen.get("screened"):
        wrong = sum(int(screen.get(k, 0)) for k in ("misnamed", "chimera", "other_gene"))
        rows.append(dict(database="MIDORI2", marker="mitochondrial genes", records=screen["screened"],
                         agree=int(screen["screened"]) - int(screen.get("flagged", 0)), disagree=wrong,
                         unlike=0, undecided=int(screen.get("weak", 0)), species="-", species_agree="-",
                         species_differ="-",
                         how="the build's reference screen (reference_screen.tsv)"))
    return rows


def eukore_ssu(package, species, fasta):
    """The package's genomes' 18S (rrna/records of source genome, none of confidence low; one per species' 99% type)
    written to fasta -> {record: species taxid}."""
    rows = {}
    for r in benchmarks.rows_of(os.path.join(package, "rrna", "records.tsv")):
        if r["source"] in ("genome", "curated_unit") and r.get("confidence") != "low" and r["species_taxid"] in species:
            rows[r["record"]] = r
    taxa, typed = {}, set()
    path = os.path.join(package, "rrna", "records.fna.gz")
    with open(fasta, "w") as fh:
        for name, seq in (odb_mod.read_fasta(path) if os.path.isfile(path) else ()):
            r = rows.get(name)
            if not r:
                continue
            if r.get("type", "-") != "-":
                if (r["species_taxid"], r["type"]) in typed:
                    continue
                typed.add((r["species_taxid"], r["type"]))
            part = seq if r["kind"] == "18S" else next(
                (seq[a:b] for n, a, b in rrna.parse_regions(r["regions"]) or () if n == "18S"), "")
            if len(part) >= MIN_LENGTH:
                fh.write(f">{name}\n{part}\n")
                taxa[name] = int(r["species_taxid"])
    return taxa


def paf_hits(path):
    """A minimap2 PAF's alignments per query (consecutive) -> (query, [(identity, target)]) of those covering
    MIN_COVER of the query or of the target."""
    query, hits = None, []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            if f[0] != query:
                if query is not None:
                    yield query, hits
                query, hits = f[0], []
            cover = max((int(f[3]) - int(f[2])) / max(1, int(f[1])), (int(f[8]) - int(f[7])) / max(1, int(f[6])))
            if int(f[10]) and cover >= MIN_COVER:
                hits.append((int(f[9]) / int(f[10]), f[5]))
    if query is not None:
        yield query, hits


def verdict_of(hits, taxa, sp):
    """A reference's record of species sp against Eukore's 18S (hits [(identity, record)]) -> (verdict, its best
    identity to the species' own genome's or None, (identity, species) of its best other species or None): `agree`
    (AGREE or more to its own), `disagree` (AGREE or more to another species' and nearer to it by MARGIN: that
    species' gene under this name, or this species' genome another's), `unlike` (aligned, but under AGREE to its own
    genome's, or to others' alone: the two databases hold the species' 18S otherwise; a divergent copy, or one of
    the two is another organism's), `undecided` (not aligned to any)."""
    own = max((i for i, t in hits if taxa.get(t) == sp), default=None)
    other = max(((i, taxa[t]) for i, t in hits if t in taxa and taxa[t] != sp), default=None)
    if own is not None and own >= AGREE:
        return "agree", own, other
    if other and other[0] >= AGREE and other[0] >= (own or 0) + MARGIN:
        return "disagree", own, other
    return ("unlike" if own is not None or other else "undecided"), own, other


def shared_rank(ncbi, a, b):
    """The deepest rank of genus .. phylum two species share, or `none`."""
    mine, theirs = ncbi.rank_ids(a), ncbi.rank_ids(b)
    return next((r for r in RANKS if mine.get(r) and mine.get(r) == theirs.get(r)), "none")


def ssu_agreement(package, species, name, fasta, table, ncbi, minimap2, threads, work):
    """Another 18S reference's records (fasta, its record<TAB>taxid table) of the species Eukore holds a genome's 18S
    of, aligned to Eukore's genomes' 18S -> (agreement.tsv's row or None, agreement_records.tsv's rows)."""
    os.makedirs(work, exist_ok=True)
    ref = os.path.join(work, "eukore_18S.fna")
    taxa = eukore_ssu(package, species, ref)
    held = set(taxa.values())
    wanted = {}
    with util.open_text(table) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                sp = ncbi.merged.get(int(f[1]), int(f[1]))
                if ncbi.rank.get(sp) != "species":
                    sp = ncbi.species_of(sp)
                if sp in held and ncbi.rank.get(sp) == "species":
                    wanted[f[0].split()[0]] = sp
    queries = os.path.join(work, "queries.fna")
    with open(queries, "w") as fh:
        for rec, seq in odb_mod.read_fasta(fasta):
            if rec in wanted and len(seq) >= MIN_LENGTH:
                fh.write(f">{rec}\n{seq}\n")
            else:
                wanted.pop(rec, None)
    if not wanted or not taxa:
        return None, []
    paf, err = os.path.join(work, "hits.paf"), os.path.join(work, "minimap2.err")
    with open(paf, "w") as out, open(err, "w") as eh:
        run = subprocess.run([minimap2, "-x", "asm20", "-c", "--secondary=yes", "-N", "50", "-p", "0.5",
                              "-t", str(max(1, threads)), ref, queries], stdout=out, stderr=eh)
    if run.returncode:
        with open(err) as eh:
            raise RuntimeError(f"minimap2 failed (exit {run.returncode}) on {name}: {eh.read().strip()[-300:]}")
    counts, by_species, records, seen = collections.Counter(), collections.defaultdict(set), [], set()
    for rec, hits in paf_hits(paf):
        sp = wanted.get(rec)
        if sp is None or rec in seen:
            continue
        seen.add(rec)
        verdict, own, other = verdict_of(hits, taxa, sp)
        counts[verdict] += 1
        by_species[sp].add(verdict)
        if verdict in ("disagree", "unlike"):
            records.append(dict(database=name, record=rec, verdict=verdict, species_taxid=sp,
                                name=ncbi.name.get(sp, "-"),
                                own_identity=f"{100 * own:.1f}" if own is not None else "-",
                                other_taxid=other[1] if other else "-",
                                other_name=ncbi.name.get(other[1], "-") if other else "-",
                                other_identity=f"{100 * other[0]:.1f}" if other else "-",
                                shared_rank=shared_rank(ncbi, sp, other[1]) if other else "-"))
    for rec, sp in wanted.items():  # not aligned at all: undecided
        if rec not in seen:
            counts["undecided"] += 1
            by_species[sp].add("undecided")
    row = dict(database=name, marker="18S", records=sum(counts.values()), agree=counts["agree"],
               disagree=counts["disagree"], unlike=counts["unlike"], undecided=counts["undecided"],
               species=len(by_species), species_agree=sum("agree" in v for v in by_species.values()),
               species_differ=differing(by_species),
               how=f"its records aligned to Eukore's genomes' 18S (minimap2): {AGREE:.0%} or more alike to the "
                   f"species' own agree; as alike to another species' and nearer by {MARGIN:.0%} disagree; the "
                   "other aligned ones are unlike")
    return row, sorted(records, key=lambda r: (r["verdict"], r["name"], r["record"]))


# --- the comparison table ------------------------------------------------------------------------------------------

def read_words(path=TABLE):
    """benchmarks/databases.tsv -> {database: {column: its words}}, in the file's order."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    head = lines[0]
    for f in lines[1:]:
        row = dict(zip(head, f + [""] * (len(head) - len(f))))
        out[row["database"]] = {k: row.get(k) or "-" for k in WORDS}
    return out


def fasta_size(fasta, table):
    """(records, bases) of a reference's FASTA, its records of the taxid table alone."""
    wanted = set()
    with util.open_text(table) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                wanted.add(f[0].split()[0])
    records = bases = 0
    for rec, seq in odb_mod.read_fasta(fasta):
        if rec in wanted:
            records += 1
            bases += len(seq)
    return records, bases


def eukore_size(package):
    """(records, bases) of the package: its BUSCO and marker copies whole, rRNA records and organelle genes."""
    records = bases = 0
    for name in ("busco/copies.tsv", "markers/copies.tsv", "rrna/records.tsv", "organelle/genes.tsv"):
        for r in benchmarks.rows_of(os.path.join(package, *name.split("/"))):
            n = benchmarks.count_of(r.get("length"))
            if n:
                records += 1
                bases += n
    return records, bases


def versions(package, bench):
    """{database: the version read}: the build's references (the manifest's inputs), the benchmarks' downloads."""
    out = {}
    try:
        with open(os.path.join(package, "manifest.json")) as fh:
            manifest = json.load(fh)
    except (OSError, ValueError):
        manifest = {}
    inputs = manifest.get("inputs") or {}
    out[EUKORE] = f"{manifest.get('version', '-')} (Eukore {manifest.get('eukore', '-')})"
    for name, ref in (inputs.get("references") or {}).items():
        out[name] = str(ref.get("version", "-"))
    for name, v in (inputs.get("sources") or {}).items():  # EUKARYOME's and ROD's versions
        if isinstance(v, str):
            out.setdefault(name, v)
    odb = inputs.get("odb")
    if isinstance(odb, dict):
        out[benchmarks.ODB_NAME] = " ".join(str(odb[k]) for k in ("lineage", "date") if odb.get(k)) or "-"
    if bench:
        out.update(prep.versions(bench))
    return out


def database_rows(package, bench, sets, compare, sizes=True):
    """databases.tsv's rows: Eukore, then every database the benchmarks counted, then the listed ones not read."""
    words, version = read_words(), versions(package, bench)
    mine = sets[EUKORE]
    fastas = {name: (fasta, table) for name, fasta, table in compare}
    rows = []
    for name in [EUKORE] + sorted(n for n in sets if n != EUKORE) + [n for n in words if n not in sets]:
        held = sets.get(name)
        size = ("-", "-")
        if sizes and name == EUKORE:
            size = eukore_size(package)
        elif sizes and name in fastas:
            size = fasta_size(*fastas[name])
        rows.append(dict(database=name, **words.get(name, {k: "-" for k in WORDS}), version=version.get(name, "-"),
                         species=len(held) if held is not None else "-",
                         both=len(held & mine) if held is not None else "-",
                         eukore_only=len(mine - held) if held is not None and name != EUKORE else "-",
                         database_only=len(held - mine) if held is not None and name != EUKORE else "-",
                         records=size[0], bases=size[1]))
    return rows


# --- all of it -------------------------------------------------------------------------------------------------------

def run(package, bench, out, ncbi, species, summary, theirs=None, compare=(), minimap2=None, threads=1, tmp=None,
        log=print):
    """The databases' tables of package into out, their headline figures into summary (databases, upset, novel,
    agreement). theirs: {name: NCBI species} of compare.py's other references; compare: [(name, FASTA, table)]."""
    sets = species_sets(out, species, theirs, ncbi)
    if len(sets) == 1:
        summary.setdefault("left_out", {})["databases"] = \
            "no other database's species (database_species.tsv: a package built since 2026-11-12; or prep.py's)"
    rows = database_rows(package, bench, sets, compare)
    benchmarks.write(out, "databases.tsv", DATABASE_COLUMNS, rows)
    summary["databases"] = {r["database"]: {k: r[k] for k in ("version", "species", "both", "eukore_only",
                                                              "database_only", "records", "bases")} for r in rows}
    rows = upset_rows(sets, species, ncbi)
    benchmarks.write(out, "upset.tsv", UPSET_COLUMNS, rows)
    everywhere = [r for r in rows if r["supergroup"] == "all"]
    summary["upset"] = {"databases": sorted(sets), "species": sum(r["species"] for r in everywhere),
                        "eukore_alone": next((r["species"] for r in everywhere if r["databases"] == EUKORE), 0),
                        "in_every_database": next((r["species"] for r in everywhere if r["n"] == len(sets)), 0)}
    rows = novel_rows(package, species, sets)
    benchmarks.write(out, "novel.tsv", NOVEL_COLUMNS, rows)
    summary["novel"] = {r["marker"]: {k: r[k] for k in ("references", "genome", "named", "new")}
                        for r in rows if r["supergroup"] == "all"}
    rows, records = published_agreement(package), []
    minimap2 = minimap2 or shutil.which("minimap2")
    for name, fasta, table in compare:
        if name not in SSU_REFERENCES:
            continue
        if not minimap2:
            summary.setdefault("left_out", {})["agreement"] = "no minimap2 (the benchmarks' environment)"
            break
        work = tempfile.mkdtemp(prefix="agreement.", dir=tmp or out)
        try:
            row, found = ssu_agreement(package, species, name, fasta, table, ncbi, minimap2, threads, work)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        if row:
            rows.append(row)
            records += found
    if rows:
        benchmarks.write(out, "agreement.tsv", AGREEMENT_COLUMNS, rows)
        benchmarks.write(out, "agreement_records.tsv", RECORD_COLUMNS, records)
        summary["agreement"] = {f"{r['database']} ({r['marker']})": {k: r[k] for k in AGREEMENT_COLUMNS[2:10]}
                                for r in rows}
    log(f"databases: {len(sets) - 1} beside Eukore; {summary['upset']['eukore_alone']:,} species in Eukore alone; "
        + "; ".join(f"{m}: {v['new']} of {v['genome']} new" for m, v in summary["novel"].items()))
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--bench", required=True, help="prep.py's folder (BENCH)")
    p.add_argument("--out", required=True, help="the results' folder (compare.py's: its summary.json is added to)")
    p.add_argument("--taxdump", help="NCBI's whole taxdump (default: BENCH/downloads/taxdump.tar.gz)")
    p.add_argument("--minimap2", default=None, help="minimap2 (default: on PATH)")
    p.add_argument("--tmp", help="the alignments' temporary files (default: in --out)")
    p.add_argument("-t", "--threads", type=int, default=8)
    a = p.parse_args(argv)
    from . import compare as compare_mod  # noqa: PLC0415
    compare_mod.open_results(a.package, a.out)
    whole = prep.taxdump_of(a.bench, a.taxdump)
    ncbi = ncbi_mod.Ncbi(whole or os.path.join(a.package, "taxonomy"))
    summary = compare_mod.load_summary(a.out)
    compare = prep.read_references(a.bench)
    theirs = {name: benchmarks.compared_species(table, ncbi) for name, _fasta, table in compare}
    run(a.package, a.bench, a.out, ncbi, compare_mod.species_of(a.package), summary, theirs, compare, a.minimap2,
        a.threads, a.tmp, log=lambda text: print(text, file=sys.stderr, flush=True))
    compare_mod.write_summary(a.package, a.out, summary)
    json.dump({k: summary[k] for k in ("upset", "novel", "agreement") if k in summary}, sys.stdout, indent=1)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
