#!/usr/bin/env python3
"""Are euk30's loci called in prokaryote genomes by EUKDB's caller? (since 2026-11-05; EukDetect2 dropped the BUSCOs it
found in 971 gut bacteria, which hold almost no archaea or alphaproteobacteria.) A panel of prokaryote genomes (one per
GTDB genus of a genome table, isolates first, then the most complete) is searched with ALL dataset proteins of the
loci and their ancestral variants (stricter than a build, which takes the 30 nearest species'), miniprot then
hmmsearch as loci.call. Every region scoring MIN_RATIO of its locus' cutoff or more is written, gated by the caller's
identity and coverage gates ("gated") and not ("relaxed"); confirmed = its best HMM is its locus at the cutoff.

    python3 tools/prokaryote_panel.py --table genomes.tsv --genome_dir DIR --odb WORK/sources/eukaryota_odb12 \\
        --out PANEL [--loci FILE] [-j 3 -t 2]

The table: accession, lineage (GTDB, d__ to s__), genome_category (mag or not), checkm2_completeness; a genome is
<genome_dir>/<accession>.fna.gz. Each genome's rows go into PANEL/done/<accession>.tsv (a rerun skips those done);
tools/rank_loci.py --panel PANEL reads them.
"""

import argparse
import csv
import functools
import multiprocessing
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukdb import loci, odb as odb_mod, sources  # noqa: E402

MIN_RATIO = 0.3
COLUMNS = ("accession", "mode", "locus", "score", "cutoff", "ratio", "best_hmm", "hmm_coverage", "identity",
           "coverage", "contig", "start", "end", "confirmed")


def panel(table, genome_dir, limit=None):
    """One genome per GTDB genus (isolates first, then the most complete) -> [(accession, lineage, path)]."""
    by_genus = {}
    with open(table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lin = r["lineage"]
            path = os.path.join(genome_dir, r["accession"] + ".fna.gz")
            if not lin.startswith(("d__Bacteria", "d__Archaea")) or not os.path.isfile(path):
                continue
            genus = lin.split(";g__")[1].split(";")[0] if ";g__" in lin else lin
            key = (r.get("genome_category") == "mag", -float(r.get("checkm2_completeness") or 0), r["accession"])
            if genus not in by_genus or key < by_genus[genus][0]:
                by_genus[genus] = (key, r["accession"], lin, path)
    return [v[1:] for v in sorted(by_genus.values())][:limit]


def search(item, odb_folder, wanted, out, threads=2, miniprot="miniprot", hmmsearch="hmmsearch"):
    acc, lin, path = item
    done = os.path.join(out, "done", acc + ".tsv")
    if os.path.isfile(done):
        return acc
    odb = odb_mod.Odb(odb_folder, wanted)
    queries = {f"{locus}~{name}": seq for locus in odb.loci for _t, name, seq in odb.proteins[locus]}
    queries.update({f"{locus}~{name}": seq for locus in odb.loci for name, seq in odb.ancestral[locus]})
    with tempfile.TemporaryDirectory() as tmp:
        faa = os.path.join(tmp, "q.faa")
        with open(faa, "w") as fh:
            fh.writelines(f">{q}\n{p}\n" for q, p in queries.items())
        run = subprocess.run([miniprot, "-t", str(threads), "-I", "--gff", "--trans", f"--outs={loci.OUTS}",
                              f"--outn={loci.OUTN}", path, faa], capture_output=True, text=True)
        al = loci.parse_alignments(run.stdout) if run.returncode == 0 else []
        lengths = {q: len(p) for q, p in queries.items()}
        cands = {}
        for mode, reg in (("gated", loci.regions(al, lengths)),
                          ("relaxed", loci.regions([dict(a, identity=1.0) for a in al], {q: 1 for q in lengths}))):
            for rs in reg.values():
                for r in rs:
                    cands[f"{mode}{len(cands)}"] = (mode, r)
        scores = loci.hmmsearch({n: r["protein"] for n, (_m, r) in cands.items()},
                                odb.hmm_db(os.path.join(tmp, "loci.hmm")), hmmsearch, threads, tmp)
    os.makedirs(os.path.dirname(done), exist_ok=True)
    with open(done + ".part", "w", newline="\n") as fh:
        fh.write(f"#{acc}\t{lin}\t{len(al)} alignments\n")
        for name, (mode, r) in cands.items():
            got = {odb.locus_of_hmm.get(h, h): v for h, v in scores.get(name, {}).items()}
            own, cov = got.get(r["locus"], (0.0, 0.0))
            best = max(got, key=lambda h: got[h][0]) if got else "-"
            cut = odb.cutoff[r["locus"]]
            if own >= MIN_RATIO * cut:
                fh.write("\t".join(map(str, [acc, mode, r["locus"], f"{own:.1f}", cut, f"{own / cut:.2f}", best,
                                             f"{cov:.2f}", r["identity"], f"{r.get('coverage', 0):.2f}", r["contig"],
                                             r["start"], r["end"], best == r["locus"] and own >= cut])) + "\n")
    os.replace(done + ".part", done)
    return acc


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--table", required=True)
    p.add_argument("--genome_dir", required=True)
    p.add_argument("--odb", required=True)
    p.add_argument("--loci", help="loci, one per line (default data/euk30_loci.tsv's)")
    p.add_argument("--out", required=True)
    p.add_argument("--limit", type=int)
    p.add_argument("-j", "--jobs", type=int, default=3)
    p.add_argument("-t", "--threads", type=int, default=2)
    opts = p.parse_args(argv)
    if opts.loci:
        with open(opts.loci) as fh:
            wanted = [line.split("\t")[0].strip() for line in fh if line.strip() and not line.startswith("locus")]
    else:
        wanted = sources.read_loci_table()
    os.makedirs(os.path.join(opts.out, "done"), exist_ok=True)
    odb_mod.Odb(opts.odb, wanted).hmm_db(os.path.join(opts.out, "euk30.hmm"))  # the loci searched (rank_loci reads it)
    items = panel(opts.table, opts.genome_dir, opts.limit)
    print(f"panel: {len(items)} genomes (one per genus), {len(wanted)} loci", flush=True)
    work = functools.partial(search, odb_folder=opts.odb, wanted=wanted, out=opts.out, threads=opts.threads)
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for i, _acc in enumerate(pool.imap_unordered(work, items), 1):
            if i % 20 == 0 or i == len(items):
                print(f"{i}/{len(items)} searched", flush=True)
    print("DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
