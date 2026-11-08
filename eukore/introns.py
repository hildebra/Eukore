"""The BUSCO copies' introns (since 2026-11-07; busco/introns.tsv, busco/intron_sites.tsv): each intron's phase, splice
sites, class and whether it looks like a minor-spliceosome (U12) intron, and its place in the gene as its locus' HMM
match state (align.py), so that an intron is compared across species: the sites a locus' introns share, how many
genomes of each group hold one there (intron gain and loss along the tree).

- phase: the coding bases before it (exons, frameshifts' bases out), mod 3: 0 between codons, 1 or 2 inside one.
- class: its first and last two bases: GT-AG (the major spliceosome's), GC-AG, AT-AC (most U12 introns), other.
- u12: t where its 5' end is U12's ([AG]TATCCTT, one difference allowed) and a U12 branch point (TCCTT[AG]A[CT], one
  difference) lies BRANCH_SPAN bases from its 3' end; else f.
- aa_position: the protein residue it falls in or before (0-based); hmm_position: that residue's match state, or -
  where the residue is in an insert (or the protein was not aligned).
- site: the locus' site it belongs to, `<hmm_position>.<phase>`: introns of the same phase at match states SITE_SLACK
  apart or less are one site (an intron slides a codon or two between species' alignments).
"""

import collections
import statistics

U12_DONOR = "ATATCCTT"
U12_BRANCH = "TCCTTAAC"
BRANCH_SPAN = (8, 40)
SITE_SLACK = 1
INTRON_COLUMNS = ("accession", "locus", "intron", "start", "end", "length", "phase", "donor", "acceptor", "class",
                  "u12", "aa_position", "hmm_position", "site")
SITE_COLUMNS = ("locus", "site", "hmm_position", "phase", "genomes", "share")


def _close(seq, motif, ok=None, most=1):
    """Whether seq differs from motif at most `most` places (ok {index: set of bases}: other bases allowed there)."""
    if len(seq) != len(motif):
        return False
    diff = sum(1 for i, (a, b) in enumerate(zip(seq, motif)) if a != b and a not in (ok or {}).get(i, ()))
    return diff <= most


def u12_like(intron):
    """Whether an intron's sequence looks like a minor-spliceosome (U12) intron's."""
    s = intron.upper()
    if not _close(s[:8], U12_DONOR, {0: {"G"}}):
        return False
    region = s[max(0, len(s) - BRANCH_SPAN[1]):len(s) - BRANCH_SPAN[0]]
    return any(_close(region[i:i + 8], U12_BRANCH, {5: {"G"}, 7: {"T"}}) for i in range(len(region) - 7))


def class_of(donor, acceptor):
    pair = f"{donor}-{acceptor}"
    return pair if pair in ("GT-AG", "GC-AG", "AT-AC") else "other"


def of_copy(copy, features, states=None):
    """A copy's introns ([(kind, start, end)] features on the copy; states: its residues' match states) -> rows
    without accession, locus and site."""
    shifted = sorted((int(s), int(e)) for k, s, e in features if k == "frameshift")
    exons = sorted((int(s), int(e)) for k, s, e in features if k == "exon")
    out, before = [], 0
    for i, ((s1, e1), (s2, _e2)) in enumerate(zip(exons, exons[1:])):
        before += e1 - s1 - sum(min(e, e1) - max(s, s1) for s, e in shifted if s < e1 and e > s1)
        a, b = e1, s2
        if b <= a:
            continue
        seq = copy[a:b].upper()
        aa = before // 3
        state = states[aa] if states and aa < len(states) else None
        out.append(dict(intron=i + 1, start=a, end=b, length=b - a, phase=before % 3, donor=seq[:2],
                        acceptor=seq[-2:], **{"class": class_of(seq[:2], seq[-2:])},
                        u12="t" if u12_like(seq) else "f", aa_position=aa,
                        hmm_position=state if state is not None else "-"))
    return out


def sites(rows):
    """Introns of one locus ([rows with hmm_position, phase]) -> {(hmm_position, phase): site name}, the site of each
    place: places of the same phase SITE_SLACK apart or less joined to the most common of them."""
    counts = collections.Counter((r["hmm_position"], r["phase"]) for r in rows if r["hmm_position"] != "-")
    named = {}
    centres = []
    for (pos, phase), _n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        centre = next((c for c in centres if c[1] == phase and abs(c[0] - pos) <= SITE_SLACK), None)
        if centre is None:
            centre = (pos, phase)
            centres.append(centre)
        named[(pos, phase)] = f"{centre[0]}.{centre[1]}"
    return named


def annotate(collected, states, group_of=None):
    """The package's introns (collected's intron_rows: of_copy's of each copy, accession and locus set, made by
    package.collect while it read the copies; states {(accession, locus): match states}) -> (intron rows
    (INTRON_COLUMNS), site rows (SITE_COLUMNS, and share_<group> with group_of {accession: group}), the share columns,
    per genome {accession: dict(introns, introns_per_kb_cds, median_intron_length)})."""
    rows, by_locus = [], collections.defaultdict(list)
    for row in collected["intron_rows"]:
        s = states.get((row["accession"], row["locus"]))
        aa = row["aa_position"]
        state = s[aa] if s and aa < len(s) else None  # of_copy's own rule
        row["hmm_position"] = state if state is not None else "-"
        rows.append(row)
        by_locus[row["locus"]].append(row)
    aligned = collections.defaultdict(set)  # the genomes each locus was aligned in: a site's share is of them
    for acc, locus in states:
        aligned[locus].add(acc)
    group_of = group_of or {}
    groups = sorted(set(group_of.values()))
    site_rows = []
    for locus, rs in sorted(by_locus.items()):
        named = sites(rs)
        holders = collections.defaultdict(set)
        for r in rs:
            r["site"] = named.get((r["hmm_position"], r["phase"]), "-")
            if r["site"] != "-":
                holders[r["site"]].add(r["accession"])
        n = len(aligned[locus]) or 1
        per_group = collections.Counter(group_of.get(a) for a in aligned[locus])
        for site, accs in sorted(holders.items(), key=lambda kv: tuple(map(int, kv[0].split(".")))):
            pos, phase = site.split(".")
            row = dict(locus=locus, site=site, hmm_position=pos, phase=phase, genomes=len(accs),
                       share=f"{len(accs) / n:.3f}")
            held = collections.Counter(group_of.get(a) for a in accs)
            for g in groups:
                row[f"share_{g}"] = f"{held[g] / per_group[g]:.3f}" if per_group[g] else "-"
            site_rows.append(row)
    cds = collections.Counter()
    for (acc, _locus), s in collected.get("cds", {}).items():
        cds[acc] += len(s)
    per_genome = collections.defaultdict(list)
    for r in rows:
        per_genome[r["accession"]].append(r["length"])
    genome = {acc: dict(introns=len(per_genome.get(acc, [])),
                        introns_per_kb_cds=f"{1000 * len(per_genome.get(acc, [])) / cds[acc]:.2f}",
                        median_intron_length=int(statistics.median(per_genome[acc])) if per_genome.get(acc) else "-")
              for acc in cds}
    return rows, site_rows, [f"share_{g}" for g in groups], genome
