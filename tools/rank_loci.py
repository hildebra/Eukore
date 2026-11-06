#!/usr/bin/env python3
"""How each BUSCO of eukaryota_odb12 suits euk30 (since 2026-11-05): its presence and conservation across all
eukaryotes, its presence in fungi and protists (EukDetect2's genomes, the measure euk30 was first chosen by), its copy
length, and whether Eukore's caller finds it in prokaryotes. -> a table ranked by `score`; data/euk30_loci.tsv's
columns come from it.

    python3 tools/rank_loci.py --odb WORK/sources/eukaryota_odb12 --taxdump WORK/sources/taxdump.tar.gz \\
        [--eukdetect2 DIR] [--panel DIR ...] --out ranking.tsv

- presence: of the dataset's 456 species (OrthoDB 12.1), the share in each supergroup (groups.GROUPS, by NCBI taxid)
  with a protein of the BUSCO (refseq_db.faa.gz); `min_presence` the weakest supergroup's (those of 5 species or
  more), `presence` all species'.
- conservation: each dataset protein's hmmsearch score against the BUSCO's own HMM, in bits per residue; the median
  per supergroup and over all (`conservation`); `min_conservation` the weakest supergroup's. A low one is a BUSCO
  whose sequence drifts far between lineages (hard to call, and to align reads to); a very high one tells species
  apart less.
- fungi_share, protist_share, median_length: as euk30 was first chosen (Fungi and Protist genomes of EukDetect2 with
  a copy of 8,000 bases or fewer; the median length of the copies), with --eukdetect2.
- prokaryotes: genomes of a prokaryote panel (tools/prokaryote_panel.py's output folders, --panel) in which Eukore's
  caller confirms the BUSCO at its cutoff, with its gates (`prokaryote_calls`) and without (`prokaryote_relaxed`);
  `-` where the BUSCO was not searched. Informational: a homologue in prokaryotes harms nothing where reads cannot
  cross between the domains (tools/prokaryote_separation.py measures it); broad eukaryote representation counts more
  (the user's decision, 2026-11-05).
- score: the mean of balanced_presence (the supergroups' presences averaged, each supergroup alike: the dataset is 57%
  animals), min_presence, min(1, conservation / 1.5) and, with EukDetect2, fungi_share and protist_share (from ~7,000
  real fungal and protist genomes, where the dataset holds few protists): how broadly and how callably the BUSCO
  represents eukaryotes.

--select N writes the N best-scoring BUSCOs whose median copy length is --max_median_length or less (a long copy is
fitted into protal's 16 kb by cutting introns) as a loci table (data/euk30_loci.tsv's format), in score order. With
--eukaryote_panel, a BUSCO is selected only where its call rate is --min_group_call or more in every panel supergroup
but those of FLOOR_EXEMPT: the score's mean lets a BUSCO missing from half a protist supergroup through on its strength
elsewhere (2026-11-06: SNF8 at rank 29 was called in 1 of 3 Haptista and 53% of SAR).
"""

import argparse
import collections
import glob
import os
import re
import statistics
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import groups, ncbi as ncbi_mod, odb as odb_mod  # noqa: E402

MIN_GROUP = 5
MIN_PANEL = 3      # a supergroup of the eukaryote panel counts with this many genomes
BROKEN_SHARE = 1 / 3  # a panel genome with fewer complete calls than this share of its supergroup's median: a
#                     broken assembly (a 4.8 Mb fern "reference genome", 3 of 129), left out
FLOOR_EXEMPT = ("Metamonada",)  # reduced, fast-evolving genomes (Giardia, Trichomonas): few BUSCOs reach half of them
CALLED = {"Single": 1.0, "Duplicated": 1.0, "Fragmented": 0.5}  # a call's weight in a call rate
LOCUS_RE = re.compile(r"-(\d+at\d+)-S\d+$")


def supergroup(ncbi, taxid):
    """groups.group_of: the package's supergroups (until 2026-11-06 rank_loci's own, without Vertebrata: its
    presence_ columns did not meet the panel's call_ ones)."""
    return groups.group_of(ncbi, taxid)


def conservation(folder, loci, hmmsearch="hmmsearch", threads=4):
    """{locus: {protein record: bits per residue}} of every dataset protein against its own BUSCO's HMM."""
    odb = odb_mod.Odb(folder, loci)
    out = collections.defaultdict(dict)
    with tempfile.TemporaryDirectory() as tmp:
        for locus in loci:
            faa = os.path.join(tmp, "p.faa")
            lengths = {}
            with open(faa, "w") as fh:
                for _t, name, seq in odb.proteins[locus]:
                    fh.write(f">{name}\n{seq}\n")
                    lengths[name] = len(seq)
            tbl = os.path.join(tmp, "t.tsv")
            subprocess.run([hmmsearch, "--cpu", str(threads), "--noali", "-o", os.devnull, "--tblout", tbl,
                            odb.hmm[locus], faa], check=True)
            with open(tbl) as fh:
                for line in fh:
                    if not line.startswith("#"):
                        f = line.split()
                        out[locus][f[0]] = float(f[5]) / max(lengths.get(f[0], 1), 1)
    return out


def eukdetect2_shares(folder, max_len=8000):
    """{locus: (fungi_share, protist_share, median_length)} from EukDetect2's link table and FASTA."""
    record_of, group_of = {}, {}
    with open(os.path.join(folder, "busco_taxid_genome_link.txt")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            m = LOCUS_RE.search(f[0]) if len(f) >= 3 else None
            if m:
                record_of[f[0]] = (f[2], m.group(1))
                group_of[f[2]] = "Fungi" if f[0].startswith("Fungi") else "Protist"
    longest = collections.defaultdict(dict)
    for name, seq in odb_mod.read_fasta(os.path.join(folder, "eukdb.fasta")):  # a wrapped FASTA's records whole
        if name in record_of:
            g, locus = record_of[name]
            longest[g][locus] = max(len(seq), longest[g].get(locus, -1))
    n = collections.Counter(group_of[g] for g in longest)
    out = {}
    for locus in {l for d in longest.values() for l in d}:
        have = {grp: [d[locus] for g, d in longest.items() if group_of[g] == grp and locus in d]
                for grp in ("Fungi", "Protist")}
        lengths = have["Fungi"] + have["Protist"]
        out[locus] = tuple(sum(x <= max_len for x in have[grp]) / max(n[grp], 1) for grp in ("Fungi", "Protist")) \
            + (int(statistics.median(lengths)) if lengths else 0,)
    return out


def panel_calls(folders):
    """{locus: (genomes confirmed with the caller's gates, without them)} of prokaryote panel folders."""
    gated, relaxed, searched = collections.defaultdict(set), collections.defaultdict(set), set()
    for folder in folders:
        hmm = os.path.join(folder, "euk30.hmm")
        if os.path.isfile(hmm):
            with open(hmm) as fh:
                searched.update(line.split()[1] for line in fh if line.startswith("NAME"))
        for path in glob.glob(os.path.join(folder, "done", "*.tsv")):
            with open(path) as fh:
                fh.readline()
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f[13] == "True":
                        (gated if f[1] == "gated" else relaxed)[f[2]].add(f[0])
    return gated, relaxed, searched


def cached(folder, name, make):
    """make()'s result, kept as JSON in folder/name (None: not kept)."""
    import json  # noqa: PLC0415
    path = os.path.join(folder, name) if folder else None
    if path and os.path.isfile(path):
        with open(path) as fh:
            return json.load(fh)
    got = make()
    if path:
        os.makedirs(folder, exist_ok=True)
        with open(path, "w") as fh:
            json.dump(got, fh)
    return got


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--odb", required=True)
    p.add_argument("--taxdump", required=True)
    p.add_argument("--eukdetect2", help="EukDetect2's files (busco_taxid_genome_link.txt, eukdb.fasta)")
    p.add_argument("--panel", nargs="*", default=[], help="tools/prokaryote_panel.py output folders")
    p.add_argument("--eukaryote_panel", help="tools/eukaryote_panel.py's output folder: the BUSCOs' call rates")
    p.add_argument("--hmmsearch", default="hmmsearch")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--out", required=True)
    p.add_argument("--select", type=int, help="write the N best BUSCOs as a loci table (--select_out)")
    p.add_argument("--select_out", default="euk30_loci.tsv")
    p.add_argument("--max_median_length", type=int, default=5000)
    p.add_argument("--min_group_call", type=float, default=0.5,
                   help="--select: the least call rate in each eukaryote panel supergroup (FLOOR_EXEMPT aside)")
    p.add_argument("--cache", help="a folder for the slow parts (conservation, EukDetect2's shares), reused by a rerun")
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(opts.taxdump)
    with open(os.path.join(opts.odb, "scores_cutoff")) as fh:
        loci = [line.split()[0] for line in fh if line.strip()]
    odb = odb_mod.Odb(opts.odb, loci)
    species = odb.species()
    group = {t: supergroup(ncbi, t) for t in species}
    sizes = collections.Counter(group.values())
    print(f"{len(species)} dataset species: {dict(sizes.most_common())}", flush=True)
    names = {}
    links = os.path.join(opts.odb, "links_to_ODB12.txt")
    if os.path.isfile(links):
        with open(links) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                names[f[0]] = f[1] if len(f) > 1 else "?"
    cons = cached(opts.cache, "conservation.json", lambda: conservation(opts.odb, loci, opts.hmmsearch, opts.threads))
    shares = cached(opts.cache, "eukdetect2_shares.json",
                    lambda: eukdetect2_shares(opts.eukdetect2)) if opts.eukdetect2 else {}
    gated, relaxed, searched = panel_calls(opts.panel)
    groups = [g for g, n in sizes.items() if n >= MIN_GROUP]
    called = {}  # supergroup -> [{locus: status}] of the eukaryote panel's genomes
    if opts.eukaryote_panel:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import eukaryote_panel  # noqa: PLC0415
        members = collections.defaultdict(list)
        for _acc, (grp, _order, calls) in eukaryote_panel.read(opts.eukaryote_panel).items():
            members[grp].append(calls)
        complete = {g: [sum(1 for s in c.values() if s in ("Single", "Duplicated")) for c in m]
                    for g, m in members.items()}
        broken = 0
        for g, m in members.items():
            floor = BROKEN_SHARE * statistics.median(complete[g])
            kept = [c for c, n in zip(m, complete[g]) if n >= floor]
            broken += len(m) - len(kept)
            members[g] = kept
        print(f"eukaryote panel: {broken} genomes under {BROKEN_SHARE:.2f} of their supergroup's median left out")
        called = {g: m for g, m in members.items() if len(m) >= MIN_PANEL}
        print(f"eukaryote panel: {sum(len(m) for m in members.values())} genomes; supergroups of {MIN_PANEL}+: "
              f"{', '.join(f'{g} {len(m)}' for g, m in sorted(called.items()))}", flush=True)
    rows = []
    for locus in loci:
        holders = {t for t, _n, _s in odb.proteins[locus]}
        presence = {g: sum(1 for t in species if group[t] == g and t in holders) / sizes[g] for g in groups}
        by_group = collections.defaultdict(list)
        for t, name, _s in odb.proteins[locus]:
            if name in cons[locus]:
                by_group[group[t]].append(cons[locus][name])
        allbits = [v for vs in by_group.values() for v in vs]
        con = {g: statistics.median(by_group[g]) if by_group[g] else 0.0 for g in groups}
        row = dict(locus=locus, name=names.get(locus, "?"), presence=len(holders & set(species)) / len(species),
                   balanced_presence=statistics.mean(presence.values()),
                   min_presence=min(presence.values()), weakest=min(presence, key=presence.get),
                   conservation=statistics.median(allbits) if allbits else 0.0, min_conservation=min(con.values()),
                   **{f"presence_{g}": v for g, v in presence.items()})
        hm = locus.split("at")[0]
        row["prokaryote_calls"] = len(gated.get(locus, ())) if hm in searched else "-"
        row["prokaryote_relaxed"] = len(relaxed.get(locus, ())) if hm in searched else "-"
        if locus in shares:
            row.update(fungi_share=shares[locus][0], protist_share=shares[locus][1], median_length=shares[locus][2])
        if called:
            rates = {g: statistics.mean(CALLED.get(calls.get(locus, "Missing"), 0.0) for calls in members)
                     for g, members in called.items()}
            row.update(call_rate=statistics.mean(rates.values()), min_call_rate=min(rates.values()),
                       weakest_call=min(rates, key=rates.get), **{f"call_{g}": v for g, v in rates.items()})
            parts = [row["call_rate"], row["min_call_rate"], min(1.0, row["conservation"] / 1.5)]
        else:
            parts = [row["balanced_presence"], row["min_presence"], min(1.0, row["conservation"] / 1.5)]
        if locus in shares:
            parts += [shares[locus][0], shares[locus][1]]
        row["score"] = sum(parts) / len(parts)
        rows.append(row)
    rows.sort(key=lambda r: (-r["score"], r["locus"]))
    columns = ["rank", "locus", "score", "call_rate", "min_call_rate", "weakest_call", "presence", "balanced_presence",
               "min_presence", "weakest", "conservation", "min_conservation", "fungi_share", "protist_share",
               "median_length", "prokaryote_calls", "prokaryote_relaxed", "name"] + \
        [f"call_{g}" for g in sorted(called)] + [f"presence_{g}" for g in groups]
    with open(opts.out, "w", newline="\n") as fh:
        fh.write("\t".join(columns) + "\n")
        for i, r in enumerate(rows, 1):
            r["rank"] = i
            fh.write("\t".join(f"{r[c]:.3f}" if isinstance(r.get(c), float) else str(r.get(c, "-")) for c in columns)
                     + "\n")
    print(f"-> {opts.out}: {len(rows)} BUSCOs; supergroups of {MIN_GROUP}+ species: {', '.join(groups)}")
    if opts.select:
        chosen = [r for r in rows if not r.get("median_length") or r["median_length"] <= opts.max_median_length]
        chosen = [r for r in chosen if all(r[f"call_{g}"] >= opts.min_group_call for g in called if g not in FLOOR_EXEMPT)]
        chosen = chosen[:opts.select]
        table = ("locus", "fungi_share", "protist_share", "median_length", "call_rate", "min_call_rate",
                 "weakest_call", "presence", "balanced_presence", "min_presence", "weakest", "conservation",
                 "min_conservation", "prokaryote_calls", "score", "name")
        with open(opts.select_out, "w", newline="\n") as fh:
            fh.write("\t".join(table) + "\n")
            for r in chosen:
                fh.write("\t".join(f"{r[c]:.3f}" if isinstance(r.get(c), float) else str(r.get(c, "-"))
                                   for c in table) + "\n")
        print(f"-> {opts.select_out}: the {len(chosen)} chosen")
    return 0


if __name__ == "__main__":
    sys.exit(main())
