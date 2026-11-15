"""Each clade's expected BUSCOs and each genome's tier (since 2026-11-14; the first full build's evaluation).

The 129 BUSCOs of eukaryota_odb12 are not all held by every eukaryote: the 2026-11 build's chromosome-level assemblies
of Metamonada held a median 24% of them, Discoba's 56%, Alveolata's 60% (at most 67%), animals' and plants' 100%. A
lineage's losses (reduced parasites: microsporidia, apicomplexans, Metamonada) read as missing data, and a flat
completeness cut would drop real genomes of those clades. So each group (groups.group_of) has an expected set: the
BUSCOs found (Single, Duplicated or Fragmented) in EXPECTED_SHARE or more of its reference genomes (whole assemblies:
REFERENCE_LEVELS, else contigs of REFERENCE_CONTIG mean length; REFERENCE_MIN of them, else every BUSCO). A genome's
`complete_expected` is the share of its group's expected BUSCOs it holds whole (Single or Duplicated).

A representative genome's tier (genomes.tsv `tier`):
- partial: NCBI's partial genome, or one small for its clade (genomes.SMALL), holding under POOR of its group's
  expected BUSCOs whole (or its group having none): a part of a genome (since 2026-11-15; until then whatever it held);
- poor: under POOR of its group's expected BUSCOs whole: a fragment, a skim, or another organism's genome with a
  few of the species' genes (the evaluation: 2,333 of the 2,664 NCBI genomes under half complete had contigs of a
  mean under 5 kb);
- draft: an assembly of contigs under DRAFT_CONTIG bases on average (a skim), but holding its genes;
- transcriptome: a TSA set;
- contaminated: one of the assemblies NCBI excludes (contaminated, unverified source organism, ...), taken where
  every assembly of its species is (genomes.ONLY_EXCLUDED; until 2026-11-14 such a species was left out);
- genome: the others.
A genome with FOREIGN_MIN or more copies another clade's (qc.cross_summary: a contamination FCS-GX did not report) is
flagged `foreign copies` in genomes.tsv (since 2026-11-14), its copies already low.
The package leaves out the copies of the DROPPED tiers' genomes (eukore_build.py --keep_poor_genomes keeps them, each
copy flagged poor_genome, confidence low); a draft's copies are flagged draft_genome (medium). The genome stays the
species' representative, its rRNA and organelle records kept.
"""

import collections

FOREIGN_MIN = 10
FOREIGN_FLAG = "foreign copies"
EXPECTED_SHARE = 0.5
REFERENCE_MIN = 3
REFERENCE_LEVELS = ("Complete Genome", "Chromosome")
REFERENCE_CONTIG = 100_000   # a scaffold-level assembly of contigs this long on average serves as a reference too
POOR = 0.5
DRAFT_CONTIG = 5_000
PARTIAL_FLAGS = {"partial", "partial genome", "small for its clade"}
EXCLUDED_FLAG = "excluded (the species' only assemblies)"  # genomes.ONLY_EXCLUDED
OLD_EXCLUDED_FLAG = "excluded, the species' only assemblies"  # its text until 2026-11-15 (a comma: read as two flags)
TIERS = ("genome", "draft", "transcriptome", "poor", "partial", "contaminated")
DROPPED = ("poor", "partial", "contaminated")
FOUND = ("Single", "Duplicated", "Fragmented")
WHOLE = ("Single", "Duplicated")
ALL = "all"  # the expected set of a group without REFERENCE_MIN references: every BUSCO


def mean_contig(genome):
    """A genome row's bases per contig, or None (no size or contig count)."""
    try:
        size, contigs = int(genome.get("genome_size") or 0), int(genome.get("contigs") or 0)
    except (TypeError, ValueError):
        return None
    return size / contigs if size and contigs else None


def is_transcriptome(genome):
    return genome.get("source") == "tsa" or genome.get("level") == "transcriptome"


def flags_of(genome):
    flags = genome.get("flags") or "-"
    if isinstance(flags, str):  # a package's genomes.tsv, or collect's row: the flags joined with commas
        flags = flags.replace(OLD_EXCLUDED_FLAG, EXCLUDED_FLAG).split(",")
    return {EXCLUDED_FLAG if f == OLD_EXCLUDED_FLAG else f for f in flags} - {"-"}


def is_reference(genome):
    """Whether a representative genome may tell its group's expected BUSCOs: a whole assembly, not partial."""
    if is_transcriptome(genome) or flags_of(genome) & (PARTIAL_FLAGS | {EXCLUDED_FLAG}) \
            or genome.get("source") == "no_assembly":
        return False
    if genome.get("level") in REFERENCE_LEVELS:
        return True
    m = mean_contig(genome)
    return m is not None and m >= REFERENCE_CONTIG


def expected_sets(genomes, statuses, group_of, loci):
    """{group: frozenset of expected BUSCOs} (ALL: every one, for the groups without REFERENCE_MIN references).
    genomes: {accession: row} of the representatives; statuses: {accession: {locus: status}}; group_of: {accession:
    group}; loci: the BUSCOs."""
    found = collections.defaultdict(collections.Counter)
    references = collections.Counter()
    for acc, g in genomes.items():
        if not is_reference(g):
            continue
        group = group_of.get(acc, "-")
        references[group] += 1
        for locus, status in statuses.get(acc, {}).items():
            if status in FOUND and locus in loci:
                found[group][locus] += 1
    out = {ALL: frozenset(loci)}
    for group, n in references.items():
        if n >= REFERENCE_MIN:
            out[group] = frozenset(l for l, c in found[group].items() if c >= EXPECTED_SHARE * n)
    return out


def complete_expected(statuses, expected):
    """A genome's share of its expected BUSCOs whole (Single or Duplicated), or None (an empty expected set)."""
    if not expected:
        return None
    return sum(1 for l in expected if statuses.get(l) in WHOLE) / len(expected)


def tier_of(genome, complete):
    """A representative genome's tier (module docstring) from its row and complete_expected."""
    flags = flags_of(genome)
    if EXCLUDED_FLAG in flags:  # NCBI excluded every assembly of its species (contaminated, ...)
        return "contaminated"
    # a partial genome holds few of its group's expected BUSCOs: one that holds POOR of them or more is not a part of
    # its genome, whatever its size says (since 2026-11-15: "small for its clade" took chromosome-level genomes of 100%,
    # Setaria italica and Spirodela polyrhiza among them, its median a family's or an order's)
    if flags & PARTIAL_FLAGS and (complete is None or complete < POOR):
        return "partial"
    if complete is not None and complete < POOR:
        return "poor"
    if is_transcriptome(genome):
        return "transcriptome"
    m = mean_contig(genome)
    if m is not None and m < DRAFT_CONTIG:
        return "draft"
    return "genome"


def assign(genomes, statuses, group_of, loci):
    """Each representative genome's expected set, complete_expected and tier, set in its row (`expected_loci`,
    `complete_expected`, `tier`) -> {group: expected set}."""
    sets = expected_sets(genomes, statuses, group_of, loci)
    for acc, g in genomes.items():
        expected = sets.get(group_of.get(acc, "-"), sets[ALL])
        c = complete_expected(statuses.get(acc, {}), expected)
        g["expected_loci"] = len(expected)
        g["complete_expected"] = "-" if c is None else f"{c:.3f}"
        g["tier"] = tier_of(g, c)
    return sets


# Eukore's own query proteins (since 2026-11-14): BUSCO's dataset holds few proteins of the divergent lineages (the
# 2026-11 build's held-out calls matched OrthoDB's gene model in 46% of Alveolata's BUSCOs, 62% of Discoba's, 85% of
# animals'), and miniprot misses genes whose nearest query is a phylum away. The package stage writes, of each group
# whose genomes' median complete share of every BUSCO is under QUERY_GROUPS_UNDER, its tier-genome representatives'
# clean copies' proteins (Single, confidence high, own species confirmed, not rescued, no stop inside) into
# WORK/eukore_queries.faa.gz, one per locus and genus (QUERIES_PER_GENUS); the next build's genomes of those groups take
# them as further queries (odb.Odb's extra), never a protein of their own species.
QUERY_GROUPS_UNDER = 0.85
QUERIES_PER_GENUS = 1
QUERIES = "eukore_queries.faa.gz"
QUERIES_STATE = "eukore_queries.json"


def query_groups(genomes, statuses, group_of, loci):
    """The groups whose tier-genome representatives hold a median under QUERY_GROUPS_UNDER of every BUSCO whole."""
    shares = collections.defaultdict(list)
    for acc, g in genomes.items():
        if g.get("tier") == "genome" and loci:
            s = statuses.get(acc, {})
            shares[group_of.get(acc, "-")].append(sum(1 for l in loci if s.get(l) in WHOLE) / len(loci))
    return sorted(g for g, v in shares.items() if v and sorted(v)[len(v) // 2] < QUERY_GROUPS_UNDER)


def query_proteins(copy_rows, proteins, genomes, groups, group_of, genus_of, kinds=()):
    """The query proteins (module comment) -> [(name <locus>_<species taxid>_eukore_<accession>, protein)]."""
    taken, out = set(), []
    for r in sorted(copy_rows, key=lambda r: (r["locus"], r["accession"])):
        acc, locus = r["accession"], r["locus"]
        g = genomes.get(acc, {})
        if locus in kinds or g.get("tier") != "genome" or group_of.get(acc) not in groups:
            continue
        if r.get("status") != "Single" or r.get("confidence") != "high" or r.get("own_species") != "t" \
                or r.get("rescued") == "t" or str(r.get("internal_stops") or 0) not in ("0", "-"):
            continue
        protein = (proteins.get((acc, locus)) or "").rstrip("*")
        genus = genus_of(int(g["species_taxid"]))
        if not protein or "*" in protein or (locus, genus) in taken:
            continue
        taken.add((locus, genus))
        out.append((f"{locus}_{int(g['species_taxid'])}_eukore_{acc}", protein))
    return out


def summary(genomes):
    """{tier: genomes} of the representatives."""
    return dict(collections.Counter(g.get("tier", "-") for g in genomes.values()))
