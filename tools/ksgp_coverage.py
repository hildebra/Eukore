#!/usr/bin/env python3
"""A package's eukaryote rRNA and organelle coverage against a KSGP release's (since 2026-11-14, format 7), for a user
that replaces KSGP's eukaryote records by Eukore's (protal):

    python3 tools/ksgp_coverage.py Eukore_2026.11 /path/to/KSGP_v5.2_plus2.tax --taxdump WORK/sources/taxdump.tar.gz \
        --out coverage.json

- the NCBI species and genera KSGP's nuclear 18S name (by KSGP's labels, eukore.rrna.label_taxon's names: species,
  genus) that hold an rRNA record in the package, by group (species.tsv's groups: Metazoa, Fungi, ...); a KSGP name
  NCBI does not hold is counted apart;
- the package's rRNA records by source, rank and confidence tier, its flags (the audit's among them) and the
  published records the audit cut; its organelle genes by source and tier;
- its species of organelle genes alone, by organelle; KSGP's plastid and mitochondrial species the package holds
  organelle genes of.
--taxdump: NCBI's whole taxonomy, to tell a KSGP species the package lacks from a name NCBI does not hold (without it
the package's own taxonomy/, which holds only its own taxa).
"""

import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import labels, links, ncbi as ncbi_mod, package, rrna  # noqa: E402


def ksgp_taxa(tax, ncbi, index):
    """KSGP's eukaryote records by compartment -> ({compartment: {NCBI species: records}}, {compartment: NCBI genera},
    {compartment: names NCBI does not hold})."""
    species = collections.defaultdict(collections.Counter)
    genera, unknown = collections.defaultdict(set), collections.defaultdict(set)
    with open(tax) as fh:
        for line in fh:
            rid, _, label = line.rstrip("\n").partition("\t")
            kind = links.kind_of(label)
            if not kind:
                continue
            parts = {p.split("__", 1)[0]: p.split("__", 1)[1] for p in rrna.ksgp_label(label).split(";")}
            genus = parts.get("g", "unclassified")
            if genus != "unclassified":
                g = index.get(rrna.re.sub(r"\(.*?\)", "", genus).strip().lower())
                if g and ncbi.rank.get(g) == "genus":
                    genera[kind].add(g)
            if parts.get("s", "unclassified") == "unclassified":
                continue
            name = rrna.species_name(genus, parts["s"].replace("_", " "))
            if not labels.is_species_name(name):
                continue
            t = index.get(name.lower())
            sp = ncbi.species_of(t) if t else None
            if sp and ncbi.rank.get(sp) == "species":
                species[kind][sp] += 1
            else:
                unknown[kind].add(labels.key(name))
    return species, genera, unknown


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("package")
    p.add_argument("ksgp", help="a KSGP release's plus2 tax file (KSGP_v<V>_plus2.tax)")
    p.add_argument("--taxdump", help="NCBI's taxdump.tar.gz (the build's WORK/sources/taxdump.tar.gz)")
    p.add_argument("--out", help="the report as JSON")
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(opts.taxdump or os.path.join(opts.package, "taxonomy"))
    index = rrna.name_index(ncbi)
    tax, version = opts.ksgp, labels.ksgp_version_of(opts.ksgp)
    k_species, k_genera, k_unknown = ksgp_taxa(tax, ncbi, index)
    records = list(package.read_table(opts.package, "rrna/records.tsv"))
    genes = list(package.read_table(opts.package, "organelle/genes.tsv"))
    rows = list(package.read_table(opts.package, "species.tsv"))

    def species_of(t):
        t = int(t)
        sp = ncbi.species_of(t)
        return sp if ncbi.rank.get(sp) == "species" else None
    held = {s for r in records for s in [species_of(r["species_taxid"])] if s}
    held_genera = {g for r in records for g in [ncbi.rank_ids(int(r["taxid"])).get("genus")] if g}
    by_group = collections.defaultdict(lambda: collections.Counter())
    for sp in set(k_species["nuclear"]) | held:
        g = labels.group(ncbi, sp)
        for key in (g, "all"):
            c = by_group[key]
            c["ksgp_species"] += sp in k_species["nuclear"]
            c["eukore_species"] += sp in held
            c["both"] += sp in k_species["nuclear"] and sp in held
    for gen in k_genera["nuclear"] | held_genera:
        g = labels.group(ncbi, gen)
        for key in (g, "all"):
            c = by_group[key]
            c["ksgp_genera"] += gen in k_genera["nuclear"]
            c["eukore_genera"] += gen in held_genera
            c["genera_both"] += gen in k_genera["nuclear"] and gen in held_genera
    coverage = {g: dict(c, species_missing=c["ksgp_species"] - c["both"],
                        genera_missing=c["ksgp_genera"] - c["genera_both"]) for g, c in sorted(by_group.items())}

    def tiers(rows_, column):
        out = collections.defaultdict(collections.Counter)
        for r in rows_:
            out[r.get(column, "-")][r.get("confidence", "-")] += 1
        return {k: dict(v) for k, v in sorted(out.items())}
    organelle_only = {int(r["taxid"]) for r in rows if r.get("source") == "organelle"}
    kinds = collections.defaultdict(set)
    for r in genes:
        kinds[r["organelle"]].add(int(r["species_taxid"]))
    report = {
        "package": os.path.basename(os.path.normpath(opts.package)), "ksgp": version,
        "ksgp_names_not_ncbi_species": {k: len(v) for k, v in k_unknown.items()},
        "coverage_nuclear": coverage,
        "rrna_by_source": tiers(records, "source"), "rrna_by_rank": tiers(records, "rank"),
        "rrna_by_label_from": dict(collections.Counter(r.get("label_from", "-") for r in records)),
        # the audit's (format 7): each flag's records, and the published records it cut
        "rrna_flags": dict(collections.Counter(f for r in records for f in r.get("flags", "-").split(",")
                                               if f not in ("-", "")).most_common()),
        "rrna_trimmed": dict(collections.Counter(r.get("source", "-") for r in records
                                                 if r.get("trimmed", "-") not in ("-", ""))),
        "organelle_by_source": tiers(genes, "source"),
        "species_by_source": dict(collections.Counter(r.get("source", "-") for r in rows)),
        "organelle_only_species": {k: len(organelle_only & s) for k, s in sorted(kinds.items())},
        "organelle_only_species_all": len(organelle_only),
        "ksgp_organelle_species": {k: {"ksgp": len(k_species[k]), "in_package": len(set(k_species[k]) & kinds[k])}
                                   for k in ("plastid", "mitochondrion", "nucleomorph")},
    }
    text = json.dumps(report, indent=1)
    if opts.out:
        with open(opts.out, "w") as fh:
            fh.write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
