#!/usr/bin/env python3
"""How well the assembly pass' organelle gene models find RefSeq's own exons (since 2026-11-07): a sample of RefSeq's
organelle records, each searched as if it were an assembly's molecule (its own species' genes left out of the
queries), by nhmmer (its hits joined) and by miniprot (the nearest species' proteins by the organelle's code). Per
protein-coding gene and model: the genes found, those whose exons are RefSeq's exactly, and the median F1 of their
exon bases against RefSeq's. Decides `--organelle_models`'s default (docs/methods.md).

    python3 tools/organelle_models.py --sources WORK/sources --taxdump WORK/sources/taxdump.tar.gz --sample 200
"""

import argparse
import collections
import os
import random
import shutil
import statistics
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import ncbi as ncbi_mod, organelle  # noqa: E402


def exon_bases(parts):
    return {(p, s) for a, b, s in parts for p in range(a, b + 1)}


def f1(truth, got):
    if not truth or not got:
        return 0.0
    tp = len(truth & got)
    return 2 * tp / (len(truth) + len(got))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sources", required=True, help="the build's sources folder (refseq_organelle/ fetched)")
    p.add_argument("--taxdump", required=True)
    p.add_argument("--sample", type=int, default=200, help="RefSeq records per organelle")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", default="organelle_models.tsv")
    p.add_argument("--work", help="a folder RefSeq's parsed genes are kept in (organelle.refseq_genes; default: temporary)")
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(opts.taxdump)
    with tempfile.TemporaryDirectory() as tmp:
        rows = organelle.refseq_genes(opts.sources, opts.work or tmp, ncbi)
    chosen = organelle.choose(rows, ncbi)
    nhmmer, miniprot = shutil.which("nhmmer"), shutil.which("miniprot")
    rng = random.Random(opts.seed)
    # records with an intron-holding gene first (fungi's, plants'), then others
    by_kind = collections.defaultdict(list)
    for (sp, kind), rs in chosen.items():
        if any(r["gene"] in organelle.CODING for r in rs):
            by_kind[kind].append((sp, rs))
    seqs = {}
    files = organelle.release_files(opts.sources)
    sample = {}
    for kind, items in by_kind.items():
        rng.shuffle(items)
        with_introns = [x for x in items if any(r["gene"] in organelle.CODING and r["exons"] != "1" for r in x[1])]
        rest = [x for x in items if x not in with_introns]
        sample[kind] = (with_introns[:opts.sample // 2] + rest)[:opts.sample]
    wanted = {r["seqid"] for kind in sample for _sp, rs in sample[kind] for r in rs}
    for path in files:
        for rec in organelle.gbff_records(path):
            if rec["seqid"] in wanted:
                seqs[rec["seqid"]] = rec["sequence"]
    stats = collections.defaultdict(lambda: dict(genes=0, found=0, exact=0, f1=[]))
    for kind, items in sample.items():
        for sp, rs in items:
            others = {k: v for k, v in chosen.items() if k[0] != sp}
            near = organelle.Nearest(others, ncbi)
            truth = {r["gene"]: r for r in rs if r["gene"] in organelle.CODING and r["trans_spliced"] == "f"
                     and r["partial"] == "f"}
            if not truth:
                continue
            molecules = [(s, kind, seqs[s]) for s in {r["seqid"] for r in truth.values()} if s in seqs]
            got = {"nhmmer": organelle.nhmmer_genes(molecules, [q for q in near.queries(sp, kind)
                                                                if q[0].split("|")[0] in truth], nhmmer, opts.threads),
                   "miniprot": organelle.miniprot_genes(molecules, near.proteins(sp, kind), miniprot,
                                                        ncbi.code_of(sp, organelle.COMPARTMENT[kind]), opts.threads)}
            for gene, t in truth.items():
                true_bases = exon_bases(organelle.parts_read(t["parts"]))
                for model, found in got.items():
                    s = stats[(kind, gene, model, t["exons"] != "1")]
                    s["genes"] += 1
                    g = found.get(gene)
                    if not g or g["seqid"] != t["seqid"]:
                        continue
                    s["found"] += 1
                    parts = organelle.parts_read(g.get("parts")) or [(int(g["start"]), int(g["end"]), g["strand"])]
                    score = f1(true_bases, exon_bases(parts))
                    s["f1"].append(score)
                    s["exact"] += score == 1.0
    with open(opts.out, "w", newline="\n") as fh:
        fh.write("organelle\tgene\tintrons\tmodel\tgenes\tfound\texact\tmedian_f1\n")
        for (kind, gene, model, introns), s in sorted(stats.items()):
            line = (f"{kind}\t{gene}\t{'t' if introns else 'f'}\t{model}\t{s['genes']}\t{s['found']}\t{s['exact']}\t"
                    f"{statistics.median(s['f1']) if s['f1'] else 0:.3f}")
            fh.write(line + "\n")
            print(line)


if __name__ == "__main__":
    main()
