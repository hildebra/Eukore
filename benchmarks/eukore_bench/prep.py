"""The benchmarks' downloads (since 2026-11-12): every file the comparisons and the read benchmark need, into a folder
of their own (BENCH; the large references there, not in a home folder), each kept once there and whole (its size, a
gzip read through, an MD5 where the source gives one):

    BENCH/downloads/taxdump.tar.gz  NCBI's whole taxonomy (~65 MB; not with --taxdump, e.g. the build's own): a
                                package's taxonomy/ holds its own lineages alone, and a species another reference holds
                                alone would be unknown to it
    BENCH/downloads/pr2/        PR2's 18S (GitHub's release; ~60 MB)
    BENCH/downloads/metaphlan/  MetaPhlAn 4's database archive, its MD5 and marker info (~6 GB)
    BENCH/downloads/eukdetect2/ EukDetect2's database: its marker FASTA and record-to-taxid table alone (Zenodo; ~2 GB
                                of its ~7 GB, the rest bowtie2's indexes, which the comparison does not use)
    BENCH/read_background/      the read benchmark's bacteria and archaea (benchmarks/read_background.tsv; until
                                2026-11-12 the build's fetch stage fetched them)
    BENCH/read_genomes/         the read benchmark's test genomes (reads.choose of the package's genomes; linked
                                from --genomes where there, else downloaded from NCBI), tests.tsv listing them
    BENCH/compare/              the other references made comparable (references.py): <name>.fna.gz and
                                <name>_taxids.tsv, and references.tsv (name, FASTA, table) of all of them

    python3 -m eukore_bench.prep BENCH --package OUT/Eukore_2026.11 [--genomes DIR] [--only pr2,metaphlan,...]

benchmarks/eukore-bench-prep.sh runs it in the benchmarks' environment. --package gives the genomes the test genomes
are chosen from; the references are labelled by NCBI's whole taxonomy (taxdump(): --taxdump, else the download).
"""

import argparse
import json
import os
import sys
import urllib.request

from . import ROOT
from . import reads, references
from eukore import genomes as genomes_mod, ncbi as ncbi_mod, odb as odb_mod, package as package_mod, sources, util

PR2_VERSION = "5.1.1"
PR2_URL = "https://github.com/pr2database/pr2database/releases/download/v{v}/pr2_version_{v}_SSU_taxo_long.fasta.gz"
MPA_URL = "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases"
EUKDETECT_RECORD = "23479662"  # Zenodo: EukDetect2's database v1.1 (github.com/allind/EukDetect, 2026-10-28)
EUKDETECT_API = "https://zenodo.org/api/records/{record}"
EUKDETECT_FILES = ("eukdb.fasta", "busco_taxid_genome_link.txt")
READ_BACKGROUND = os.path.join(ROOT, "benchmarks", "read_background.tsv")  # data/ until 2026-11-12
REFERENCES = "references.tsv"
TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
STEPS = ("taxdump", "pr2", "metaphlan", "eukdetect2", "odb", "background", "tests", "convert")


def taxdump_of(bench, given=None):
    """The whole NCBI taxonomy the benchmarks label by: given (--taxdump), else BENCH's download, else None."""
    if given:
        return given
    path = os.path.join(bench, "downloads", "taxdump.tar.gz")
    return path if os.path.isfile(path) else None


def fetch(url, dest, md5=None, log=print):
    """url into dest unless there already (and of md5 where given) -> dest, or raises OSError."""
    if os.path.isfile(dest) and (md5 is None or util.md5(dest) == md5):
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    log(f"  {os.path.basename(dest)} <- {url}")
    why = util.fetch_file(url, dest, gzip_check=dest.endswith(".gz"))
    if why:
        raise OSError(f"{url}: {why}")
    if md5 is not None and util.md5(dest) != md5:
        os.remove(dest)
        raise OSError(f"{url}: its MD5 is not {md5} (removed; run again)")
    return dest


def text_at(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode()


def pr2(bench, log=print):
    name = os.path.basename(PR2_URL.format(v=PR2_VERSION))
    return fetch(PR2_URL.format(v=PR2_VERSION), os.path.join(bench, "downloads", "pr2", name), log=log)


def metaphlan(bench, version=None, log=print):
    """MetaPhlAn 4's database (version, else the server's mpa_latest) -> (archive, marker info)."""
    version = version or text_at(f"{MPA_URL}/mpa_latest").split()[0]
    folder = os.path.join(bench, "downloads", "metaphlan")
    md5_file = fetch(f"{MPA_URL}/{version}.md5", os.path.join(folder, f"{version}.md5"), log=log)
    with open(md5_file) as fh:
        md5 = fh.read().split()[0]
    archive = fetch(f"{MPA_URL}/{version}.tar", os.path.join(folder, f"{version}.tar"), md5, log)
    info = fetch(f"{MPA_URL}/{version}_marker_info.txt.bz2", os.path.join(folder, f"{version}_marker_info.txt.bz2"),
                 log=log)
    return archive, info


def eukdetect2(bench, record=EUKDETECT_RECORD, log=print):
    """EukDetect2's marker FASTA and record-to-taxid table (Zenodo's record, each file's MD5 checked) -> paths."""
    meta = json.loads(text_at(EUKDETECT_API.format(record=record)))
    files = {f["key"]: f for f in meta.get("files", [])}
    out = []
    for name in EUKDETECT_FILES:
        if name not in files:
            raise OSError(f"Zenodo record {record}: no {name}")
        f = files[name]
        md5 = f.get("checksum", "").partition("md5:")[2] or None
        out.append(fetch(f["links"]["self"], os.path.join(bench, "downloads", "eukdetect2", name), md5, log))
    return out


def background(bench, log=print):
    reads.fetch_background(READ_BACKGROUND, os.path.join(bench, "read_background"), log)


def test_genomes(bench, package, genomes=None, per_group=reads.PER_GROUP, log=print):
    """The read benchmark's test genomes (reads.choose of the package's NCBI genomes) into BENCH/read_genomes: linked
    from genomes (the build's folder) where there, else downloaded (genomes.tsv's ftp) -> their accessions."""
    folder = os.path.join(bench, "read_genomes")
    os.makedirs(folder, exist_ok=True)
    ncbi = ncbi_mod.Ncbi(os.path.join(package, "taxonomy"))
    rows = {r["accession"]: r for r in package_mod.read_table(package, "genomes.tsv")
            if r["representative"] == "t" and r["accession"].startswith(("GCA_", "GCF_"))
            and r.get("ftp") not in (None, "", "-")}
    reps = {a: (int(r["species_taxid"]), genomes_mod.genome_file(folder, a)) for a, r in rows.items()}
    chosen = reads.choose(ncbi, reps, per_group, exists=lambda _p: True)
    kept = []
    for acc in chosen:
        dest = genomes_mod.genome_file(folder, acc)
        there = genomes_mod.genome_file(genomes, acc) if genomes else None
        try:
            if os.path.isfile(dest):
                pass
            elif there and os.path.isfile(there):
                os.symlink(os.path.abspath(there), dest)
            else:
                ftp = rows[acc]["ftp"].rstrip("/").replace("ftp://", "https://")
                fetch(f"{ftp}/{os.path.basename(ftp)}_genomic.fna.gz", dest, log=log)
            kept.append(acc)
        except OSError as e:
            log(f"WARNING: test genome {acc} left out ({e})")
    with open(os.path.join(folder, "tests.tsv"), "w", newline="\n") as fh:
        fh.write("accession\tspecies_taxid\n" + "".join(f"{a}\t{rows[a]['species_taxid']}\n" for a in kept))
    log(f"  read benchmark: {len(kept)} test genomes in {folder}")
    return kept


def convert(bench, taxdump, log=print):
    """The downloads made comparable into BENCH/compare (references.py), references.tsv listing every one there."""
    out = os.path.join(bench, "compare")
    os.makedirs(out, exist_ok=True)
    down = os.path.join(bench, "downloads")
    pr2_fasta = os.path.join(down, "pr2", os.path.basename(PR2_URL.format(v=PR2_VERSION)))
    if os.path.isfile(pr2_fasta) and not os.path.isfile(os.path.join(out, "PR2.fna.gz")):
        references.main(["pr2", "--fasta", pr2_fasta, "--taxdump", taxdump, "--out", out])
    mpa = os.path.join(down, "metaphlan")
    tars = sorted(n for n in os.listdir(mpa) if n.endswith(".tar")) if os.path.isdir(mpa) else []
    if tars and not os.path.isfile(os.path.join(out, "MetaPhlAn4.fna.gz")):
        version = tars[-1][:-len(".tar")]
        references.main(["metaphlan", "--archive", os.path.join(mpa, tars[-1]), "--marker_info",
                         os.path.join(mpa, f"{version}_marker_info.txt.bz2"), "--out", out])
    listed = []
    for name in ("EukDetect2", "MetaPhlAn4", "PR2"):  # the first: breadth.tsv's compared reference
        if name == "EukDetect2":
            pair = tuple(os.path.join(down, "eukdetect2", f) for f in EUKDETECT_FILES)
        else:
            pair = (os.path.join(out, f"{name}.fna.gz"), os.path.join(out, f"{name}_taxids.tsv"))
        if all(os.path.isfile(p) for p in pair):
            listed.append((name,) + pair)
    with open(os.path.join(out, REFERENCES), "w", newline="\n") as fh:
        fh.write("name\tfasta\ttable\n" + "".join("\t".join(r) + "\n" for r in listed))
    log(f"  other references: {', '.join(r[0] for r in listed) or 'none'} ({os.path.join(out, REFERENCES)})")
    return listed


def odb_of(bench):
    """BUSCO's dataset in BENCH (the odb step's: downloaded, or the build's linked), or None."""
    folder = odb_mod.folder_of(os.path.join(bench, "downloads")) if bench else None
    return folder if folder and os.path.isdir(os.path.join(folder, "hmms")) else None


def odb(bench, given=None, log=print):
    """BUSCO's dataset for the tree's profiles (tree.py): the build's (given: SOURCES/eukaryota_odb12, the package's
    own version) linked into BENCH/downloads, else BUSCO's newest downloaded there."""
    down = os.path.join(bench, "downloads")
    os.makedirs(down, exist_ok=True)
    if odb_of(bench):
        return odb_of(bench)
    if given:
        if not os.path.isdir(os.path.join(given, "hmms")):
            raise OSError(f"{given}: no hmms/ in it (BUSCO's dataset folder, eukaryota_odb12/)")
        os.symlink(os.path.abspath(given), odb_mod.folder_of(down))
        return odb_mod.folder_of(down)
    return odb_mod.fetch(down, sources.download, log)


def versions(bench):
    """{name: the version downloaded} of the other references in BENCH (databases.tsv's `version`)."""
    out = {}
    down = os.path.join(bench, "downloads")
    if os.path.isfile(os.path.join(down, "pr2", os.path.basename(PR2_URL.format(v=PR2_VERSION)))):
        out["PR2"] = PR2_VERSION
    mpa = os.path.join(down, "metaphlan")
    tars = sorted(n for n in os.listdir(mpa) if n.endswith(".tar")) if os.path.isdir(mpa) else []
    if tars:
        out["MetaPhlAn4"] = tars[-1][:-len(".tar")]
    if os.path.isfile(os.path.join(down, "eukdetect2", EUKDETECT_FILES[0])):
        out["EukDetect2"] = f"Zenodo record {EUKDETECT_RECORD}"
    return out


def read_references(bench):
    """BENCH/compare/references.tsv -> [(name, FASTA, table)] (reads.compared's form)."""
    path = os.path.join(bench, "compare", REFERENCES)
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        next(fh, None)
        return [tuple(line.rstrip("\n").split("\t")[:3]) for line in fh if line.strip()]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bench", help="the benchmarks' folder (BENCH): downloads, references, test genomes")
    p.add_argument("--package", help="a finished Eukore package (Eukore_<version>/): its genomes give the test "
                                     "genomes")
    p.add_argument("--genomes", help="the build's genomes folder (eukore_build.py --genomes): test genomes there are "
                                     "linked, not downloaded")
    p.add_argument("--taxdump", help="NCBI's whole taxdump, e.g. the build's SOURCES/taxdump.tar.gz (default: "
                                     "downloaded into BENCH/downloads)")
    p.add_argument("--odb", help="BUSCO's dataset of the build (SOURCES/eukaryota_odb12: the package's own version), "
                                 "linked for the tree's profiles (default: BUSCO's newest, downloaded)")
    p.add_argument("--only", help=f"some of {','.join(STEPS)} (default: all)")
    p.add_argument("--read_genomes", type=int, default=reads.PER_GROUP, help="test genomes per supergroup")
    p.add_argument("--mpa_version", help="MetaPhlAn's database (default: the server's mpa_latest)")
    a = p.parse_args(argv)
    steps = a.only.split(",") if a.only else list(STEPS)
    unknown = set(steps) - set(STEPS)
    if unknown:
        p.error(f"--only: unknown {sorted(unknown)}")

    def log(text):
        print(text, flush=True)
    failed = []
    for step, run in (("taxdump", lambda: None if a.taxdump else fetch(
                          TAXDUMP_URL, os.path.join(a.bench, "downloads", "taxdump.tar.gz"), log=log)),
                      ("pr2", lambda: pr2(a.bench, log)),
                      ("metaphlan", lambda: metaphlan(a.bench, a.mpa_version, log)),
                      ("eukdetect2", lambda: eukdetect2(a.bench, log=log)),
                      ("odb", lambda: odb(a.bench, a.odb, log)),
                      ("background", lambda: background(a.bench, log)),
                      ("tests", lambda: test_genomes(a.bench, a.package, a.genomes, a.read_genomes, log)
                       if a.package else log("  tests: no --package, no test genomes")),
                      ("convert", lambda: convert(a.bench, taxdump_of(a.bench, a.taxdump), log)
                       if taxdump_of(a.bench, a.taxdump)
                       else log("  convert: no taxdump (the taxdump step, or --taxdump): nothing made comparable"))):
        if step not in steps:
            continue
        log(f"{step} ..")
        try:
            run()
        except OSError as e:
            failed.append(step)
            log(f"WARNING: {step}: {e}")
    if failed:
        log(f"failed: {', '.join(failed)} (run again: whole files are kept, the rest fetched)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
