#!/usr/bin/env python3
"""Eukore's marker loci checked on genomes (since 2026-11-07): each genome given whose species has a protein in the
marker set (Swiss-Prot's or NCBI's, markers.py) has its markers called as the build calls them, once with every query
and once without its own species' proteins (held out: the queries are then other species' alone, the test that is not
circular), and each call's protein is compared with the species' protein by shared 5-mers both ways (benchmarks.py's
rule): concordant (0.9 both ways: that gene), similar (0.5: that gene, another model or allele; for a family marker,
kind family, also another member of the family), discordant (another gene: a paralog the decoys did not tell),
missed. Also says, per call, its HMM score and the best decoy's on the same protein (a call's HMM score over its best
decoy's: the family's, not a paralog's).

    python3 tools/validate_markers.py --sources WORK/sources --genomes OUT/genomes --out validation.tsv [ACCESSION ...]
"""

import argparse
import collections
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import benchmarks, gencode, genomes, loci, markers, ncbi as ncbi_mod, odb as odb_mod  # noqa: E402

COLUMNS = ("accession", "species_taxid", "organism", "marker", "kind", "queries", "status", "copies", "hmm_score",
           "best_decoy", "decoy_score", "reference", "recall", "precision", "verdict")


def verdict(protein, refs):
    got = benchmarks.kmers(protein)
    best = max(refs, key=lambda a: benchmarks.shared(benchmarks.kmers(refs[a]), got)
               + benchmarks.shared(got, benchmarks.kmers(refs[a])))
    ref = benchmarks.kmers(refs[best])
    recall, precision = benchmarks.shared(ref, got), benchmarks.shared(got, ref)
    v = ("concordant" if min(recall, precision) >= benchmarks.CONCORDANT else
         "similar" if min(recall, precision) >= benchmarks.SIMILAR else "discordant")
    return best, recall, precision, v


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sources", required=True, help="the build's sources (eukaryota_odb12, eukore_markers, taxdump)")
    p.add_argument("--genomes", required=True, help="the genomes' folder (<accession>.fna.gz)")
    p.add_argument("--out", default="marker_validation.tsv")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("accessions", nargs="*", help="genomes to check (default: every one there with a UniProt protein)")
    opts = p.parse_args(argv)
    ncbi = ncbi_mod.Ncbi(os.path.join(opts.sources, "taxdump.tar.gz"))
    folder = markers.folder_of(opts.sources)
    marker_loci = markers.all_loci(folder)
    kinds = markers.kinds(folder)
    odb = odb_mod.Odb(odb_mod.folder_of(opts.sources), marker_loci, [folder])
    lineages = odb_mod.species_lineages(odb, ncbi)
    refs = collections.defaultdict(dict)  # (marker, species) -> {accession: protein}
    taxids_of = collections.defaultdict(set)  # species -> its marker proteins' taxids
    for locus, ps in odb.proteins.items():
        for taxid, name, seq in ps:
            sp = ncbi.species_of(taxid)
            refs[(locus, sp)][name.split("_", 2)[2]] = seq
            taxids_of[sp].add(taxid)
    import glob  # noqa: PLC0415
    from eukore import sources  # noqa: PLC0415
    tables = sorted(glob.glob(os.path.join(opts.sources, "ncbi_assemblies", "genbank_*.txt")))
    rows = {r["accession"]: r for r in sources.read_assemblies(tables)}
    wanted = opts.accessions or sorted(f[:-len(".fna.gz")] for f in os.listdir(opts.genomes) if f.endswith(".fna.gz"))
    out = []
    with tempfile.TemporaryDirectory() as tmp:
        hmm_db = odb.hmm_db(os.path.join(tmp, "markers.hmm"))
        for acc in wanted:
            row = rows.get(acc)
            if row is None:
                continue
            sp = int(row["species_taxid"])
            if not any((m, sp) in refs for m in marker_loci):
                continue
            path = genomes.genome_file(opts.genomes, acc)
            table = ncbi.code_of(int(row["taxid"]))
            lineage = list(ncbi.ancestors(int(row["taxid"])))
            for label, leave in (("all", ()), ("held_out", taxids_of[sp])):
                queries = odb.queries(lineage, lineages, leave_out=leave)
                got = loci.call(path, queries, odb, shutil.which("miniprot"), shutil.which("hmmsearch"),
                                opts.threads, tmp, hmm_db, table, rescue_missing=False)
                for m in marker_loci:
                    if (m, sp) not in refs:
                        continue
                    r = got["loci"].get(m, {"status": "Missing"})
                    base = dict(accession=acc, species_taxid=sp, organism=row["organism"], marker=m,
                                kind=kinds.get(m, "clean"), queries=label, status=r["status"])
                    if r["status"] == "Missing":
                        out.append(dict(base, verdict="missed"))
                        continue
                    protein, _stops = gencode.protein_of(r["spliced"], table)
                    scores = loci.hmmsearch({"c": protein}, hmm_db, shutil.which("hmmsearch"), 1, tmp).get("c", {})
                    decoys = {h: v[0] for h, v in scores.items() if h.startswith("decoy-")}
                    best_decoy = max(decoys, key=decoys.get) if decoys else "-"
                    ref, recall, precision, v = verdict(protein, refs[(m, sp)])
                    out.append(dict(base, copies=r["copies"], hmm_score=f"{r['hmm_score']:.1f}", best_decoy=best_decoy,
                                    decoy_score=f"{decoys[best_decoy]:.1f}" if decoys else "-", reference=ref,
                                    recall=f"{recall:.3f}", precision=f"{precision:.3f}", verdict=v))
            print(f"{acc} {row['organism']}: " + ", ".join(f"{r['marker']}/{r['queries']}={r['verdict']}"
                                                          for r in out if r["accession"] == acc), flush=True)
    with open(opts.out, "w", newline="\n") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        fh.writelines("\t".join(str(r.get(c, "-")) for c in COLUMNS) + "\n" for r in out)
    counts = collections.Counter((r["marker"], r["queries"], r["verdict"]) for r in out)
    for key, n in sorted(counts.items()):
        print("\t".join(key), n)


if __name__ == "__main__":
    main()
