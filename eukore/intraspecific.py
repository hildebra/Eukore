"""The barcode gap (since 2026-11-07; intraspecific/): how far each marker varies within a species, beside how far it is
from its nearest other species. The panel's species (data/intraspecific_panel.tsv) have further assemblies run
(genomes.intraspecific_jobs: other strains, isolates, BioSamples); their BUSCO copies, rRNA and organelle genes are
compared with each other and with the species' genome's.

- intra_min, intra_median: the identity (vsearch, --iddef 2, %) of the species' assemblies' copies of a marker to each
  other: their least alike pair, the median pair, of the pairs measured (- where none is); pairs_left_out (since
  2026-11-10): the pairs vsearch reports no alignment of (none over COVER of both copies at MIN_ID: a partial copy
  beside a whole one, or two copies under 50% alike), which are in neither figure. Until then such a pair was read
  as 50% identical, which "no alignment over the cover" does not say.
- inter_max: the identity of the species' genome's copy to its nearest copy of another species in the package (a
  BUSCO's: busco/copies.tsv's own check, qc.check's nearest_id; the 18S, the ITS, an organelle gene's: vsearch against
  every other species' in the package).
- gap: intra_min - inter_max; resolved: t where it is above 0 (every assembly of the species nearer to each other than
  to any other species) and the nearest other species is a congener, f where it is not above 0, - where the species
  has no other species' copy 50% alike, or (since 2026-11-09) none of its genus: a gap to a farther relative says
  nothing of telling the species from its congeners (until then t, where the package held no congener), or (since
  2026-11-10) none of its pairs is measured.
- inter_rank (since 2026-11-09): the deepest rank the nearest other species shares with the species (genus ..
  phylum, none; qc.depth), - where not known.
Pairs are aligned over COVER of both sequences at least (since 2026-11-10; of the query alone since 2026-11-09), the
species' own pairs and the genome's copy against the other species' alike; a pair's identity is the lower of its two
directions' (until 2026-11-10 the higher).

Markers: each BUSCO (its coding sequence), 18S, ITS (ITS1, 5.8S and ITS2 of an rDNA unit), and the organelle genes
(their coding sequence where known, else their record).
"""

import collections
import os
import statistics

from . import organelle as organelle_mod
from . import package, qc, rrna, util

MIN_ID = 0.5
COVER = 0.8  # a pair's alignment covers this share of both sequences (since 2026-11-10; of its query since 2026-11-09)
GAP_COLUMNS = ("marker", "species_taxid", "name", "group", "assemblies", "intra_min", "intra_median", "inter_max",
               "gap", "resolved",
               "inter_rank",  # since 2026-11-09 (format 6): the rank the nearest other species shares with it
               "pairs_left_out")  # since 2026-11-10: the species' pairs without an alignment, in no figure of the row
RANK_NAMES = ("none", "phylum", "class", "order", "family", "genus")  # qc.depth's 0 .. 5
COPY_COLUMNS = ("accession", "species_taxid", "locus", "status", "contig", "start", "end", "strand", "length",
                "cds_length", "exons", "hmm_score", "hmm_coverage", "copies", "dup_identity")


def first_rrna(rows):
    """A genome's 18S and ITS (its first rDNA unit holding them, else its first 18S record) -> {marker: sequence}."""
    out = {}
    for r in sorted(rows, key=lambda r: (r["kind"] != "rDNA", not r.get("lsu_complete"))):
        regions = {n: (a, b) for n, a, b in (r.get("regions") or ())}
        seq = r["sequence"]
        if "18S" not in out:
            if r["kind"] == "18S":
                out["18S"] = seq
            elif "18S" in regions:
                out["18S"] = seq[slice(*regions["18S"])]
        if "ITS" not in out and "ITS1" in regions and "ITS2" in regions:
            out["ITS"] = seq[regions["ITS1"][0]:regions["ITS2"][1]]
    return out


def organelle_seqs(rows):
    """Organelle gene rows (organelle.assembly_rows', or the package's) -> {marker: sequence}: a gene's coding sequence
    where known, else its record; an animal mitochondrion's SSU and LSU by their common names (12S, 16S)."""
    out = {}
    for r in rows:
        marker = r["gene"] if r["gene"] not in ("SSU", "LSU", "5S") else f"{r['organelle'][:4]}_{r['gene']}"
        seq = r.get("cds") or organelle_mod.coding_of(r, r["sequence"]) or r["sequence"]
        out.setdefault(marker, seq)
    return out


def collect(work, extras, organelle_folder):
    """The panel's further assemblies' results (their done folders, their organelle search files) -> (copy rows
    (COPY_COLUMNS), {accession: {marker: sequence}}, whole copies {accession: {locus: copy}}, rRNA {accession:
    [(name, sequence)]}, organelle genes {accession: [(name, sequence)]})."""
    from .genomes import done_dir  # noqa: PLC0415
    rows, markers, whole, rrna_out, org_out = [], {}, {}, {}, {}
    for acc, job in sorted(extras.items()):
        copies = package.done_copies(work, acc)
        if not copies:
            continue
        feats = collections.defaultdict(list)
        for f in package.done_table(work, acc, "features.tsv"):
            feats[f["locus"]].append((f["kind"], int(f["start"]), int(f["end"])))
        m = {}
        for r in package.done_table(work, acc, "loci.tsv"):
            if r["status"] == "Missing" or r["locus"] not in copies:
                continue
            rows.append(dict(r, accession=acc, species_taxid=job["species_taxid"]))
            m[r["locus"]] = qc.spliced(copies[r["locus"]], feats[r["locus"]])
        whole[acc] = copies
        rr = rrna.load(done_dir(work, acc))
        m.update(first_rrna(rr))
        rrna_out[acc] = [(f"{acc}_{r['kind']}_{i}", r["sequence"]) for i, r in enumerate(rr, 1)]
        path = organelle_mod.assembly_file(organelle_folder, acc)
        org = organelle_mod.assembly_rows(path) if os.path.isfile(path) else []
        m.update(organelle_seqs(org))
        org_out[acc] = [(f"{acc}|{r['gene']}", r["sequence"]) for r in org]
        markers[acc] = m
    return rows, markers, whole, rrna_out, org_out


def _pairs(query, db, vsearch, threads=1, maxaccepts=0):
    """vsearch of query {name: seq} against db {name: seq} -> [(query, target, identity)], or None where vsearch
    failed (since 2026-11-08: [], every pair then read as under MIN_ID and the gap called unresolved). A pair's
    alignment covers COVER of the query (since 2026-11-09: --iddef 2 reads identity over the aligned columns alone, and
    a partial record's few hundred bases matched made a nearest species or a strain look nearer than it is) and of
    the target (since 2026-11-10: the query's cover alone let a partial copy, a contig-edge 18S, inside a whole one
    count as a pair at its local identity, the whole one being the target)."""
    if not vsearch:
        return None
    if not query or not db:
        return []
    try:
        return util.vsearch_global(vsearch, query, db, MIN_ID, threads, maxaccepts, 256 if maxaccepts else 0,
                                   ["--query_cov", str(COVER), "--target_cov", str(COVER)])
    except OSError:
        return None


def barcode_gap(panel_species, rep_markers, extra_markers, species_of, inter, others, ncbi, vsearch, threads=1,
                group_of=None):
    """The barcode gap's rows (GAP_COLUMNS). panel_species: {species taxid: representative accession};
    rep_markers / extra_markers: {accession: {marker: sequence}}; species_of: {accession: species taxid}; inter: {(species
    taxid, marker): the representative's nearest other species' identity, or (identity, the depth its lineage shares
    with the species: qc.depth's)} known (the BUSCOs'); others: {marker: {name: sequence}} of every other species in the
    package for the markers whose inter is measured here, names `<species taxid>|...`."""
    by_marker = collections.defaultdict(dict)
    for acc, markers in list(rep_markers.items()) + list(extra_markers.items()):
        sp = species_of[acc]
        if sp not in panel_species:
            continue
        for marker, seq in markers.items():
            if seq:
                by_marker[marker][f"{sp}|{acc}"] = seq
    inter = dict(inter)
    for marker, db in others.items():  # the representatives' nearest other species, where not known
        query = {name: s for name, s in by_marker.get(marker, {}).items()
                 if name.split("|")[1] == panel_species.get(int(name.split("|")[0]))}
        best = {}
        for q, t, ident in _pairs(query, db, vsearch, threads, maxaccepts=64) or ():  # failed: not known
            if q.split("|")[0] != t.split("|")[0] and ident > best.get(q, (0.0, None))[0]:
                best[q] = (ident, int(t.split("|")[0]))
        cache = {}
        for name in query:
            sp = int(name.split("|")[0])
            ident, other = best.get(name, (None, None))
            inter[(sp, marker)] = None if ident is None else (
                ident, qc.depth(ncbi, sp, other, cache) if ncbi is not None else None)
    out = []
    for marker, seqs in sorted(by_marker.items()):
        hits = collections.defaultdict(dict)
        pairs = _pairs(seqs, seqs, vsearch, threads)
        if pairs is None:  # vsearch failed on this marker: no gap measured
            continue
        for q, t, ident in pairs:
            if q != t and q.split("|")[0] == t.split("|")[0]:
                key = tuple(sorted((q, t)))
                # one identity a pair: the lower of its two directions' (since 2026-11-10; until then the higher,
                # which with the query's cover alone was the partial copy's local one)
                hits[q.split("|")[0]][key] = min(hits[q.split("|")[0]].get(key, ident), ident)
        by_species = collections.defaultdict(list)
        for name in seqs:
            by_species[name.split("|")[0]].append(name)
        for sp, names in sorted(by_species.items()):
            if len(names) < 2:
                continue
            pairs = [(a, b) for i, a in enumerate(sorted(names)) for b in sorted(names)[i + 1:]]
            # the pairs measured; one vsearch reports no alignment of is left out and counted (since 2026-11-10; until
            # then read as MIN_ID, 50%: no alignment over the cover is not an identity, and with both sequences'
            # cover asked every pair of a partial copy would be the species' least alike one, at 50%)
            ids = [hits[sp][pair] for pair in pairs if pair in hits[sp]]
            near = inter.get((int(sp), marker))
            near, depth = near if isinstance(near, tuple) else (near, None)
            intra_min = min(ids) if ids else None
            # resolved against a congener alone (since 2026-11-09); nearer to a farther relative is unresolved all
            # the same; no pair measured: not measured
            resolved = "-" if near is None or intra_min is None else "f" if intra_min <= near else \
                "t" if depth == 5 else "-"
            out.append(dict(marker=marker, species_taxid=sp, name=ncbi.name.get(int(sp), "-") if ncbi else "-",
                            group=(group_of or {}).get(int(sp), "-"), assemblies=len(names),
                            intra_min="-" if intra_min is None else f"{intra_min:.1f}",
                            intra_median=f"{statistics.median(ids):.1f}" if ids else "-",
                            inter_max="-" if near is None else f"{near:.1f}",
                            gap="-" if near is None or intra_min is None else f"{intra_min - near:.1f}",
                            resolved=resolved, inter_rank="-" if depth is None else RANK_NAMES[depth],
                            pairs_left_out=len(pairs) - len(ids)))
    return out


def write(part, rows, whole, rrna_out, org_out, gap_rows):
    """intraspecific/ of the unfinished package: copies.tsv, copies.fna.gz (>accession|locus), rrna.fna.gz,
    organelle.fna.gz, barcode_gap.tsv."""
    folder = os.path.join(part, "intraspecific")
    os.makedirs(folder, exist_ok=True)
    package._tsv(os.path.join(folder, "copies.tsv"), COPY_COLUMNS, rows)
    package._tsv(os.path.join(folder, "barcode_gap.tsv"), GAP_COLUMNS, gap_rows)
    for name, records in (("copies.fna.gz", [(f"{a}|{lo}", s) for a, c in sorted(whole.items())
                                             for lo, s in sorted(c.items())]),
                          ("rrna.fna.gz", [r for a in sorted(rrna_out) for r in rrna_out[a]]),
                          ("organelle.fna.gz", [r for a in sorted(org_out) for r in org_out[a]])):
        path = os.path.join(folder, name)
        with util.bgzf_open(util.part_of(path)) as fh:
            fh.writelines(f">{n}\n{s.upper()}\n" for n, s in records)
        os.replace(util.part_of(path), path)
