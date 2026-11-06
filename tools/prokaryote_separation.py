#!/usr/bin/env python3
"""Can reads cross between a locus' eukaryote copies and its prokaryote homologues? (since 2026-11-05.) Every region a
prokaryote panel (tools/prokaryote_panel.py) confirmed at a locus' cutoff, with or without the caller's gates, is cut
from its genome and aligned (vsearch, MIN_ID or more) to the coding sequences of the locus' copies in a Eukore package
(its exons joined). -> per locus: the prokaryote regions, those with a eukaryote copy at MIN_ID or more, the best
identity (`<50` where none), its pair. A read aligner that needs ~85% identity cannot place one domain's read on the
other's gene where the best identity is far below that.

    python3 tools/prokaryote_separation.py --panel PANEL [PANEL ...] --genome_dir DIR --package OUT/Eukore_<version>
"""

import argparse
import collections
import glob
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import package, placement, qc, rrna  # noqa: E402

MIN_ID = 0.5


def regions(panels):
    """{(accession, locus): [(contig, start, end)]} of the panels' confirmed regions (gated or relaxed)."""
    out = collections.defaultdict(set)
    for folder in panels:
        for path in glob.glob(os.path.join(folder, "done", "*.tsv")):
            with open(path) as fh:
                fh.readline()
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if f[13] == "True":
                        out[(f[0], f[2])].add((f[10], int(f[11]), int(f[12])))
    return out


def eukaryote_cds(folder, loci):
    """{locus: {accession: coding sequence}} of a package's copies (qc.spliced: their exons joined, their frameshifts'
    bases left out, as the package's own reader)."""
    feats = collections.defaultdict(list)
    for r in package.read_table(folder, "busco/features.tsv"):
        if r["kind"] in ("exon", "frameshift") and r["locus"] in loci:
            feats[(r["accession"], r["locus"])].append((r["kind"], int(r["start"]), int(r["end"])))
    out = collections.defaultdict(dict)
    for locus in loci:
        path = os.path.join(folder, "busco", "copies", f"{locus}.fna.gz")
        if not os.path.isfile(path):
            continue
        for acc, seq in package.read_copies(folder, locus).items():
            f = feats.get((acc, locus))
            out[locus][acc] = qc.spliced(seq, f) if f else seq
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--panel", nargs="+", required=True)
    p.add_argument("--genome_dir", required=True)
    p.add_argument("--package", required=True)
    p.add_argument("--vsearch", default="vsearch")
    p.add_argument("-t", "--threads", type=int, default=4)
    opts = p.parse_args(argv)
    found = regions(opts.panel)
    loci = sorted({locus for _a, locus in found})
    cds = eukaryote_cds(opts.package, set(loci))
    print("locus\tprokaryote_regions\twith_a_eukaryote_copy_at_50%\tbest_identity\tpair")
    for locus in loci:
        euk = cds.get(locus, {})
        prok = {}
        for (acc, l), spans in found.items():
            if l != locus:
                continue
            contigs = placement.read_contigs(os.path.join(opts.genome_dir, acc + ".fna.gz"), {c for c, _a, _b in spans})
            for i, (contig, a, b) in enumerate(sorted(spans)):
                seq = contigs.get(contig, "")[a - 1:b].upper()
                prok[f"{acc}|{contig}|{a}"] = seq
                prok[f"{acc}|{contig}|{a}|rc"] = rrna.revcomp(seq)
        if not euk:
            print(f"{locus}\t{len(prok) // 2}\t-\t-\tno eukaryote copies in the package")
            continue
        hits = []
        with tempfile.TemporaryDirectory() as tmp:
            q, db, out = (os.path.join(tmp, f) for f in ("q.fna", "db.fna", "hits.tsv"))
            with open(q, "w") as fh:
                fh.writelines(f">{n}\n{s}\n" for n, s in prok.items() if s)
            with open(db, "w") as fh:
                fh.writelines(f">{n}\n{s}\n" for n, s in euk.items() if s)
            subprocess.run([opts.vsearch, "--usearch_global", q, "--db", db, "--id", str(MIN_ID), "--iddef", "2",
                            "--maxaccepts", "4", "--maxrejects", "256", "--userout", out, "--userfields",
                            "query+target+id", "--threads", str(opts.threads), "--quiet"], check=True)
            with open(out) as fh:
                hits = [line.rstrip("\n").split("\t") for line in fh]
        by_region = collections.defaultdict(float)
        best = (0.0, "-")
        for qn, tn, ident in hits:
            region = qn.rsplit("|rc", 1)[0]
            by_region[region] = max(by_region[region], float(ident))
            best = max(best, (float(ident), f"{region} ~ {tn}"))
        n = len({k.rsplit("|rc", 1)[0] for k in prok})
        print(f"{locus}\t{n}\t{len(by_region)}\t{best[0]:.1f}" if hits else f"{locus}\t{n}\t0\t<{MIN_ID * 100:.0f}",
              f"\t{best[1]}" if hits else "\t-", sep="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
