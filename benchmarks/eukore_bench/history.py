"""A Eukore release among releases (since 2026-11-12; compare.py runs it; benchmarks/README.md): what a reader
needs to reuse it and to follow it from one release to the next.

- sizes.tsv: the package's parts (its folders; busco/'s by kind of file), their files and bytes (the manifest's),
  with the build's cost beside in summary.json (the package stage's steps and the genomes' runs: stages.tsv,
  runtime.tsv);
- growth.tsv: per group (and `all`) and year, the species' genomes NCBI released that year (metadata/ncbi.tsv.gz's
  seq_rel_date): how fast the genomes a release can hold arrive, beside the other databases' release dates
  (benchmarks/databases.tsv's `released`, where given);
- changes.tsv, changes_species.tsv (with --previous, an earlier package): per group the species added, removed and
  kept, the genomes replaced by another assembly or by a newer version of the same, the species that gained or
  lost a genome; and those species one by one.

    python3 -m eukore_bench.history --package OUT/Eukore_2026.11 --out RESULTS [--previous OUT/Eukore_2026.05]
"""

import argparse
import collections
import json
import os
import sys

from . import databases
from eukore import benchmarks, package as package_mod, util

SIZE_COLUMNS = ("part", "files", "bytes")
GROWTH_COLUMNS = ("supergroup", "year", "genomes")
CHANGE_COLUMNS = ("supergroup", "previous", "species", "added", "removed", "kept", "replaced", "updated",
                  "gained_genome", "lost_genome")
SPECIES_COLUMNS = ("taxid", "name", "supergroup", "change", "previous_genome", "genome")
SPLIT = ("busco", "markers")  # folders whose parts are told apart: copies, cds, proteins and tables


def read_manifest(package):
    try:
        with open(os.path.join(package, "manifest.json")) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def size_rows(manifest):
    """sizes.tsv's rows: the manifest's files summed per part (a top folder, or a file at the top; busco/ and
    markers/ by their subfolders, their tables as `<folder>/ tables`), the largest first, `all` last."""
    files, sizes = collections.Counter(), collections.Counter()
    for path, value in (manifest.get("files") or {}).items():
        f = path.split("/")
        if len(f) == 1:
            part = path
        elif f[0] in SPLIT:
            part = f"{f[0]}/{f[1]}/" if len(f) > 2 else f"{f[0]}/ tables"
        else:
            part = f"{f[0]}/"
        files[part] += 1
        sizes[part] += int(value[0] if isinstance(value, (list, tuple)) else value)
    rows = [dict(part=p, files=files[p], bytes=sizes[p]) for p in sorted(sizes, key=lambda p: -sizes[p])]
    return rows + [dict(part="all", files=sum(files.values()), bytes=sum(sizes.values()))] if rows else []


def growth_rows(package, species):
    """growth.tsv's rows: per group (and `all`) and year the species' genomes (genomes.tsv's representatives) NCBI
    released that year (metadata/ncbi.tsv.gz's seq_rel_date of the accession, or of the one fetched)."""
    path = os.path.join(package, "metadata", "ncbi.tsv.gz")
    if not os.path.isfile(path):
        return []
    released = {}
    with util.open_text(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if "assembly_accession" not in head or "seq_rel_date" not in head:
            return []
        a, d = head.index("assembly_accession"), head.index("seq_rel_date")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(a, d) and f[d][:4].isdigit():
                released[f[a]] = int(f[d][:4])
    counts = collections.Counter()
    for r in package_mod.read_table(package, "genomes.tsv"):
        if r["representative"] != "t":
            continue
        year = released.get(r["accession"]) or released.get(r.get("fetched"))
        if year:
            for g in (species.get(r["species_taxid"], {}).get("supergroup", "-"), "all"):
                counts[(g, year)] += 1
    return [dict(supergroup=g, year=y, genomes=n)
            for (g, y), n in sorted(counts.items(), key=lambda kv: (kv[0][0] == "all", kv[0]))]


def assembly_of(accession):
    """An accession without its version and its GenBank or RefSeq prefix: GCA_000146045.2 and GCF_000146045.3 are one
    assembly."""
    stem = accession.split(".")[0]
    return stem[4:] if stem.startswith(("GCA_", "GCF_")) else stem


def change_rows(old, new):
    """(changes.tsv's rows, changes_species.tsv's rows) of two packages' species.tsv rows by taxid (old, new: each
    with its `supergroup`): per group (the new package's, a removed species' the old one's; and `all`) the species
    of the old package (`previous`) and the new one, those added, removed and kept; of the kept, the genomes
    `replaced` (another assembly), `updated` (the same assembly's newer version), and the species that gained or lost
    a genome."""
    counts, listed = collections.defaultdict(collections.Counter), []
    for taxid in sorted(set(old) | set(new), key=lambda t: (int(t) if t.isdigit() else 0, t)):
        a, b = old.get(taxid), new.get(taxid)
        row = b or a
        change = None
        before = a["representative"] if a else "-"
        after = b["representative"] if b else "-"
        if a is None:
            change = "added"
        elif b is None:
            change = "removed"
        elif before == "-" and after != "-":
            change = "gained_genome"
        elif before != "-" and after == "-":
            change = "lost_genome"
        elif before != after:
            change = "updated" if assembly_of(before) == assembly_of(after) else "replaced"
        for g in (row["supergroup"], "all"):
            c = counts[g]
            c["previous"] += a is not None
            c["species"] += b is not None
            c["kept"] += a is not None and b is not None
            if change:
                c[change] += 1
        if change:
            listed.append(dict(taxid=taxid, name=row.get("name", "-"), supergroup=row["supergroup"], change=change,
                               previous_genome=before, genome=after))
    rows = [dict(supergroup=g, **{k: c[k] for k in CHANGE_COLUMNS[1:]})
            for g, c in sorted(counts.items(), key=lambda kv: (kv[0] == "all", kv[0]))]
    return rows, listed


def species_of(package):
    return {r["taxid"]: dict(r, supergroup=r.get("group") or r.get("supergroup", "-"))
            for r in package_mod.read_table(package, "species.tsv")}


def run(package, out, species, summary, previous=None, log=print):
    """The release's tables of package into out, their headline figures into summary (sizes, cost, growth,
    changes)."""
    manifest = read_manifest(package)
    rows = size_rows(manifest)
    if rows:
        benchmarks.write(out, "sizes.tsv", SIZE_COLUMNS, rows)
        summary["sizes"] = {"bytes": rows[-1]["bytes"], "files": rows[-1]["files"],
                            "largest": {r["part"]: r["bytes"] for r in rows[:5] if r["part"] != "all"}}
    stages = list(benchmarks.rows_of(os.path.join(out, "stages.tsv")))
    cost = {}
    if stages:
        cost["package_stage_seconds"] = round(sum(float(r["seconds"]) for r in stages), 1)
    if (summary.get("runtime") or {}).get("seconds"):
        cost["genome_runs_seconds"] = summary["runtime"]["seconds"]
        cost["genome_runs_threads"] = summary["runtime"].get("threads")
    if cost:
        summary["cost"] = cost
    rows = growth_rows(package, species)
    if rows:
        benchmarks.write(out, "growth.tsv", GROWTH_COLUMNS, rows)
        years = {r["year"]: r["genomes"] for r in rows if r["supergroup"] == "all"}
        total = sum(years.values())
        last = max(years)
        summary["growth"] = {"genomes": total, "first_year": min(years), "last_year": last,
                             "released_in_last_two_years": round(sum(n for y, n in years.items() if y >= last - 1)
                                                                 / total, 4)}
    if previous:
        if not os.path.isfile(os.path.join(previous, "species.tsv")):
            raise OSError(f"{previous}: no species.tsv (--previous: an earlier package's folder)")
        rows, listed = change_rows(species_of(previous), species)
        benchmarks.write(out, "changes.tsv", CHANGE_COLUMNS, rows)
        benchmarks.write(out, "changes_species.tsv", SPECIES_COLUMNS, listed)
        before = read_manifest(previous)
        summary["changes"] = dict(next(({k: r[k] for k in CHANGE_COLUMNS[1:]} for r in rows
                                        if r["supergroup"] == "all"), {}),
                                  previous_version=before.get("version", "-"),
                                  previous_built=before.get("built", "-"))
        log(f"history: against {os.path.basename(os.path.normpath(previous))}: "
            f"{summary['changes'].get('added', 0):,} species added, {summary['changes'].get('removed', 0):,} removed, "
            f"{summary['changes'].get('replaced', 0):,} genomes replaced")
    return summary


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--package", required=True, help="a finished Eukore package (Eukore_<version>/)")
    p.add_argument("--out", required=True, help="the results' folder (compare.py's: its summary.json is added to)")
    p.add_argument("--previous", help="an earlier package: the changes from it (changes.tsv)")
    a = p.parse_args(argv)
    from . import compare  # noqa: PLC0415
    compare.open_results(a.package, a.out)
    summary = compare.load_summary(a.out)
    run(a.package, a.out, compare.species_of(a.package), summary, a.previous,
        log=lambda text: print(text, file=sys.stderr, flush=True))
    compare.write_summary(a.package, a.out, summary)
    json.dump({k: summary[k] for k in ("sizes", "cost", "growth", "changes") if k in summary}, sys.stdout, indent=1)
    print()
    return 0


RELEASED = "released"  # benchmarks/databases.tsv's optional column: a database's release date (YYYY-MM)


def release_dates(path=databases.TABLE):
    """{database: (year, month)} of benchmarks/databases.tsv's `released` column, where a row gives one."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    head = lines[0]
    if RELEASED not in head:
        return out
    for f in lines[1:]:
        row = dict(zip(head, f + [""] * (len(head) - len(f))))
        date = row.get(RELEASED, "")
        if date[:4].isdigit():
            out[row["database"]] = (int(date[:4]), int(date[5:7]) if date[5:7].isdigit() else 6)
    return out


if __name__ == "__main__":
    sys.exit(main())
