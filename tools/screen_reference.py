#!/usr/bin/env python3
"""The reference screen on MIDORI2's gene files (since 2026-11-10; eukore/refscreen.py, as KSGP 4.5 and 5.2 were
screened): each record's windows against RefSeq's organelle genes and the other species' records; the records that
are another organism's gene (misnamed), that hold a piece of another organism's (chimera) or that are another gene
are listed, as KSGP's releases list theirs:

  OUT/calls.tsv    every record the first pass flagged: its call, its windows and the records behind it
  OUT/drop.tsv     the records to leave out: record, call, reason
  OUT/trim.tsv     the records to cut to their own part: record, start, end (1-based, inclusive), length, reason
  OUT/review.tsv   the records flagged but not settled, kept as they are: record, call, reason
  OUT/summary.txt  counts
  --fasta: OUT/<file>.screened.fasta.gz, each file without its dropped records, the cut ones cut (their headers kept,
  " screen:start-end" added)

Several gene files are screened together: a record of one gene's file that is another file's gene shows.

    python3 tools/screen_reference.py --sources WORK/sources --out midori2_screen -t 8 \\
        MIDORI2_UNIQ_NUC_GB273_CO1_RAW.fasta.gz [MIDORI2_UNIQ_NUC_GB273_Cytb_RAW.fasta.gz ...]
"""

import argparse
import os
import shutil
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import labels, ncbi as ncbi_mod, organelle, refscreen, util  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="MIDORI2 RAW FASTA files (their gene by their name)")
    p.add_argument("--sources", required=True, help="the build's sources folder (taxdump.tar.gz, refseq_organelle/)")
    p.add_argument("--taxdump", help="NCBI's taxdump.tar.gz (default: the sources folder's)")
    p.add_argument("--work", help="a folder RefSeq's parsed genes are kept in (the build's WORK/organelle; default: "
                                  "temporary, parsed again each run)")
    p.add_argument("--out", required=True, help="the lists' folder")
    p.add_argument("--fasta", action="store_true", help="write each file screened: without its dropped records, the "
                                                        "cut ones cut")
    p.add_argument("--keep", action="store_true", help="keep the searches' files (OUT/search)")
    p.add_argument("-t", "--threads", type=int, default=4)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--vsearch", default="vsearch")
    opts = p.parse_args(argv)
    minimap2, vsearch = shutil.which(opts.minimap2), shutil.which(opts.vsearch)
    if not minimap2 or not vsearch:
        sys.exit("the screen needs minimap2 and vsearch (envs/eukore.yaml)")
    began = time.time()
    os.makedirs(opts.out, exist_ok=True)
    ncbi = ncbi_mod.Ncbi(opts.taxdump or os.path.join(opts.sources, "taxdump.tar.gz"))
    with tempfile.TemporaryDirectory() as tmp:
        rows = organelle.refseq_genes(opts.sources, opts.work or tmp, ncbi, log=print)
    if not rows:
        print("WARNING: no RefSeq organelle genomes in the sources (--stage fetch): the records are screened against "
              "one another alone")
    calls, counts = refscreen.screen([refscreen.midori2_source(f) for f in opts.files],
                                     list(refscreen.trusted_records(rows)), ncbi, minimap2, vsearch, opts.threads,
                                     tmp=opts.out, log=print,
                                     keep=os.path.join(opts.out, "search") if opts.keep else None)
    refscreen.write_lists(opts.out, calls)
    with open(os.path.join(opts.out, "summary.txt"), "w", newline="\n") as fh:
        fh.write(refscreen.summary(calls, counts, opts.files) + f"runtime: {util.duration(time.time() - began)}\n")
    print(open(os.path.join(opts.out, "summary.txt")).read(), end="")
    if opts.fasta:
        for path in opts.files:
            base = os.path.basename(path)
            for ext in (".gz", ".fasta", ".fa", ".fna"):
                base = base[:-len(ext)] if base.endswith(ext) else base
            n = refscreen.write_screened(path, os.path.join(opts.out, base + ".screened.fasta.gz"), calls,
                                         labels.midori2_gene(path))
            print(f"{base}.screened.fasta.gz: {n:,} records")


if __name__ == "__main__":
    main()
