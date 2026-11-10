#!/usr/bin/env python3
"""The annotation swap measured against an independent truth (since 2026-11-09): where a genome's assembly has two
annotations of independent pipelines (NCBI's curated or its Eukaryotic Genome Annotation Pipeline's, and Ensembl's,
its BRAKER or Helixer runs', TOGA2's projections), its BUSCOs are called as the build calls them with the test
annotation given (loci.call: a copy takes the test's model where its HMM confirms it, loci.annotated), and both
miniprot's model and the copy shipped are laid on the truth annotation's coding sequences (annotation.compare): the
introns exact, the whole intron chain exact, the coding bases shared (F1). The build's own check of the swap (its
annotation.tsv) lays miniprot's model on the annotation that is swapped in: circular for the copies that take it.

    python3 tools/annotation_swap.py --sources WORK/sources --genome G.fna.gz --taxid 31033 --truth ncbi.gff.gz \
        --test ensembl.gff.gz [--test_source ensembl] [--report G_assembly_report.txt] [-t 8] [--out swap.tsv]

A copy is `swapped` where it took the test's model. Per copy: the truth's transcript, miniprot's and the shipped
model's exact introns of the truth's, whole chain exact, coding-base F1; `better`, `worse` or `same` where swapped (by
the chain, then F1). The summary: of the swapped copies, how many the swap made better or worse against the truth.
"""

import argparse
import collections
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import annotation, loci, ncbi as ncbi_mod, odb as odb_mod  # noqa: E402

COLUMNS = ("locus", "status", "model", "truth_transcript", "truth_introns", "miniprot_exact_introns",
           "miniprot_chain", "miniprot_f1", "shipped_exact_introns", "shipped_chain", "shipped_f1", "swap")


def measure(model, truth):
    """A model (dict contig, strand, status, exons) against the truth's transcripts -> (truth transcript, its introns,
    exact introns, whole chain exact, coding-base F1)."""
    row = annotation.compare(model, truth)
    if row["transcript"] == "-":
        return "-", 0, 0, False, 0.0
    f1 = 2 * row["shared_bases"] / max(1, row["bases"] + row["ann_bases"])
    chain = row["exact_introns"] == row["ann_introns"] == row["introns"] and row["exons"] == row["ann_exons"]
    return row["transcript"], row["ann_introns"], row["exact_introns"], chain, f1


def verdict(mini, shipped):
    """better, worse or same: the shipped model against miniprot's, by the whole chain, then F1 (0.01)."""
    if shipped[3] != mini[3]:
        return "better" if shipped[3] else "worse"
    if abs(shipped[4] - mini[4]) < 0.01:
        return "same"
    return "better" if shipped[4] > mini[4] else "worse"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sources", required=True, help="the build's sources (eukaryota_odb12, taxdump.tar.gz)")
    p.add_argument("--genome", required=True)
    p.add_argument("--taxid", type=int, required=True, help="the genome's NCBI taxid (its queries, its code)")
    p.add_argument("--truth", required=True, help="the truth annotation (GFF3 or GTF, gzipped or not)")
    p.add_argument("--test", required=True, help="the annotation the copies may take (the build's swap)")
    p.add_argument("--test_source", default="ensembl", help="its source's name (copies.tsv's model)")
    p.add_argument("--report", help="NCBI's assembly report: the sequences' other names (Ensembl's chromosomes, ...)")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--out", default="annotation_swap.tsv")
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(os.path.join(opts.sources, "taxdump.tar.gz"))
    folder = odb_mod.folder_of(opts.sources)
    odb = odb_mod.Odb(folder, odb_mod.all_loci(folder))
    queries = odb.queries(list(ncbi.ancestors(opts.taxid)), odb_mod.species_lineages(odb, ncbi))
    table = ncbi.code_of(opts.taxid)
    with tempfile.TemporaryDirectory() as tmp:
        got = loci.call(opts.genome, queries, odb, shutil.which("miniprot"), shutil.which("hmmsearch"), opts.threads,
                        tmp, None, table, annotation=(opts.test, opts.report, opts.test_source))
    called = {lo: r for lo, r in got["loci"].items() if r["status"] != "Missing"}
    spans = collections.defaultdict(list)
    for r in called.values():
        for model in (r, r["miniprot"]):
            spans[r["contig"]].append((model["start"], model["end"]))
    truth = annotation.read_cds(opts.truth, spans, annotation.aliases(opts.report, opts.truth))
    rows, tally = [], collections.Counter()
    for locus, r in sorted(called.items()):
        base = dict(contig=r["contig"], strand=r["strand"], status=r["status"])
        mini = measure(dict(base, exons=r["miniprot"]["exons"]), truth)
        shipped = measure(dict(base, exons=r["exons"]), truth)
        swapped = r["model"] != "miniprot"
        swap = verdict(mini, shipped) if swapped else "-"
        tally["copies"] += 1
        tally["swapped"] += swapped
        tally[f"swap_{swap}"] += swapped
        tally["truth_found"] += mini[0] != "-"
        tally["miniprot_chain"] += mini[3]
        tally["shipped_chain"] += shipped[3]
        rows.append(dict(locus=locus, status=r["status"], model=r["model"], truth_transcript=mini[0],
                         truth_introns=mini[1], miniprot_exact_introns=mini[2], miniprot_chain="t" if mini[3] else "f",
                         miniprot_f1=f"{mini[4]:.3f}", shipped_exact_introns=shipped[2],
                         shipped_chain="t" if shipped[3] else "f", shipped_f1=f"{shipped[4]:.3f}", swap=swap))
    with open(opts.out, "w", newline="\n") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in COLUMNS) + "\n" for r in rows)
    f1 = {k: [float(r[f"{k}_f1"]) for r in rows if r["truth_transcript"] != "-"] for k in ("miniprot", "shipped")}
    print(f"{tally['copies']} copies, {tally['truth_found']} on a truth transcript; {tally['swapped']} took the "
          f"{opts.test_source} model: {tally['swap_better']} better, {tally['swap_same']} the same, "
          f"{tally['swap_worse']} worse against the truth")
    print(f"whole intron chain exact: miniprot {tally['miniprot_chain']}, shipped {tally['shipped_chain']}; "
          + "; ".join(f"{k} mean F1 {sum(v) / len(v):.3f}" for k, v in f1.items() if v))


if __name__ == "__main__":
    main()
