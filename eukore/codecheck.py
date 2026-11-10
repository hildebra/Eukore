"""Each genome's genetic code checked against its BUSCO genes (since 2026-11-07; as Codetta does with Pfam): every codon
of its copies that lies at a well conserved position of its locus' HMM (a match state whose most likely amino acid has
probability INFORMATIVE or more) adds that position's log emission probabilities to the codon's tally. The amino acid
a codon codes is the best of its tally; one of TESTED (the stops, CTG: the codons nuclear codes are known to read
otherwise) differs from NCBI's table for the genome where it has MIN_OCCURRENCES such (a stop STOP_MIN)
positions and its tally favours another amino acid by MIN_LLR nats or more. A stop of the table read through (a
ciliate's TAA, given the standard code) is compared with a uniform background: no amino acid of its own.

miniprot placed the genes with NCBI's table; a gene whose code differs is still placed (miniprot aligns through an
in-frame stop, hmmsearch confirms it with X there), so its reassigned codons are seen at the positions their amino
acid fills in other species.

genomes.tsv's gencode_check: `agrees`, `differs:<codon>><amino acid>,...` or `undecided` (fewer than MIN_CODONS codons
seen often enough); best_table: the NCBI table that best fits the codons seen (NCBI's own where it fits as well).
metadata/codons.tsv.gz: per genome and codon, the table's amino acid, the one inferred, its positions and the
log-likelihood ratio.
"""

import collections
import math

from . import align as align_mod
from . import gencode

INFORMATIVE = 0.5
MIN_OCCURRENCES = 10
STOP_MIN = 5  # a stop codon at this many conserved positions is judged: under its own reading it is at none
STOPS = ("TAA", "TAG", "TGA")


def enough(codon, n):
    return n >= (STOP_MIN if codon in STOPS else MIN_OCCURRENCES)
MIN_LLR = 10.0
MIN_CODONS = 30
BACKGROUND = math.log(1 / 20)
CODON_COLUMNS = ("accession", "codon", "table_aa", "inferred_aa", "n", "llr")
NUCLEAR_TABLES = (1, 6, 10, 12, 15, 26, 27, 28, 29, 30, 31)
# the codons a nuclear code is known to read otherwise (the stops: Q, E, Y, W, C; CTG: S, A): the verdict's. Another
# codon's tally is biased toward the amino acids HMMs favour (on the 2026-11-07 trial ATG read as L in 30 of 33
# genomes: Met is rarely a conserved position's best), so it is reported (codons.tsv) but not judged
TESTED = ("TAA", "TAG", "TGA", "CTG")


def tally(pairs, emissions):
    """Codons at informative positions ([(cds, states, locus)]; emissions {locus: [{aa: p} per state]}) -> {codon:
    (n, {aa: summed log p})}."""
    sums = collections.defaultdict(lambda: collections.defaultdict(float))
    n = collections.Counter()
    for cds, states, locus in pairs:
        em = emissions.get(locus)
        if not em or not states:
            continue
        cds = cds.upper()
        for i, state in enumerate(states):
            if state is None or state >= len(em) or 3 * i + 3 > len(cds):
                continue
            probs = em[state]
            if max(probs.values()) < INFORMATIVE:
                continue
            codon = cds[3 * i:3 * i + 3]
            if len(codon) != 3 or set(codon) - set("ACGT"):
                continue
            n[codon] += 1
            for aa, p in probs.items():
                sums[codon][aa] += math.log(max(p, 1e-6))
    return {c: (n[c], dict(sums[c])) for c in n}


def score(n, sums, aa):
    return n * BACKGROUND if aa == "*" else sums.get(aa, n * math.log(1e-6))


def judged(tallies):
    """The TESTED codons seen often enough to judge (genomes.tsv's gencode_codons): a stop not seen is still a stop (a
    codon read as an amino acid would be seen), but a CTG not seen is not judged (CTG-clade yeasts avoid the codon)."""
    return ",".join(c for c in TESTED if enough(c, tallies.get(c, (0,))[0])) or "-"


def check(tallies, table):
    """A genome's codon tallies and NCBI's table for it -> (codon rows without accession, gencode_check, best_table)."""
    code = gencode.code(table)
    rows, differs = [], []
    seen = {c: v for c, v in tallies.items() if enough(c, v[0])}
    for codon in gencode.CODONS:
        if codon not in tallies:
            continue
        n, sums = tallies[codon]
        best = max(sums, key=sums.get)
        llr = sums[best] - score(n, sums, code[codon])
        rows.append(dict(codon=codon, table_aa=code[codon], inferred_aa=best, n=n, llr=f"{llr:.1f}"))
        if codon in seen and codon in TESTED and best != code[codon] and llr >= MIN_LLR:
            differs.append(f"{codon}>{best}")
    if len(seen) < MIN_CODONS:
        return rows, "undecided", "-"
    fit = {t: sum(score(n, s, gencode.code(t)[c]) for c, (n, s) in seen.items() if c in TESTED)
           for t in NUCLEAR_TABLES + (table,)}
    best_table = max(fit, key=lambda t: (fit[t], t == table))
    return rows, ("differs:" + ",".join(differs)) if differs else "agrees", best_table


def run(collected, states, hmm_of):
    """Every genome's code checked (collected's cds and genomes; states {(accession, locus): match states}; hmm_of
    {locus: HMM path}) -> (codon rows (CODON_COLUMNS), {accession: dict(gencode_check, best_table, gencode_codons)})."""
    emissions = {locus: align_mod.read_emissions(path) for locus, path in hmm_of.items() if path}
    # not a copy a check calls another organism's (since 2026-11-09: a contaminant's codons are its own code's)
    low = {(r["accession"], r["locus"]) for r in collected.get("copy_rows", ()) if r.get("confidence") == "low"}
    pairs = collections.defaultdict(list)
    for (acc, locus), cds in collected.get("cds", {}).items():
        if (acc, locus) in states and (acc, locus) not in low:
            pairs[acc].append((cds, states[(acc, locus)], locus))
    rows, per_genome = [], {}
    for acc, g in sorted(collected["genomes"].items()):
        if acc not in pairs:
            continue
        tallies = tally(pairs[acc], emissions)
        got, verdict, best = check(tallies, int(g.get("gencode") or 1))
        rows += [dict(r, accession=acc) for r in got]
        per_genome[acc] = dict(gencode_check=verdict, best_table=best, gencode_codons=judged(tallies))
    return rows, per_genome
