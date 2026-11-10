"""Eukore's sources (since 2026-11-05; copied from protal's scripts/mini_db/eukaryotes.py, rdna_sources.py and
download_ref_data.py, so that Eukore imports nothing from them): the loci table, NCBI's GenBank and RefSeq assembly
tables, and the published eukaryote rDNA operons (EUKARYOME 2.1's long reads, ROD v1.2). EukDetect1's and EukDetect2's
copies were sources until 2026-11-05: Eukore now calls the loci in NCBI's genomes itself (odb.py, loci.py).
"""

import collections
import gzip
import json
import os
import re
import shutil
import subprocess
import zipfile

from . import util
from .util import open_text

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"
ASSEMBLIES_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank"
ASSEMBLY_GROUPS = ("fungi", "protozoa", "plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other")


def download(url, dest, md5=None):
    """A file fetched whole (util.fetch_file) and, with an md5, checked -> raises where it did not arrive."""
    why = util.fetch_file(url, dest, gzip_check=False)
    if why:
        raise OSError(f"{url}: {why}")
    if md5 and util.md5(dest) != md5:
        raise OSError(f"{dest}: MD5 {util.md5(dest)}, {md5} expected")


# No eukaryote's nuclear genome is smaller (microsporidia's, the smallest, are ~2.2 Mb): GenBank has records of 500-600
# bases as "complete genomes" of Aspergillus flavus, niger, fumigatus and others (a gene's sequence), which a build
# simulated as cells of 600 bases, and which were their species' best-ranked match. Nor is a species' genome a third of
# its other assemblies' median (an organelle's, a partial assembly).
MIN_GENOME = 2_000_000


PARTIAL_SHARE = 0.3


def read_assemblies(paths):
    """GenBank's assembly_summary.txt files -> [dict] of the latest assemblies (accession, taxid, species_taxid,
    organism, strain, isolate, level, asm_name, size, contigs, category, ftp, genes: its annotation's protein-coding
    genes, 0 for none, None where the table does not say), without those smaller than MIN_GENOME and,
    of a species with 3 assemblies or more, those smaller than PARTIAL_SHARE of their median size."""
    rows = []
    for path in paths:
        with open_text(path) as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    header = line.lstrip("#").strip().split("\t")
                    continue
                if header is None:
                    continue
                f = dict(zip(header, line.rstrip("\n").split("\t")))
                if f.get("version_status", "latest") != "latest":
                    continue
                excluded = f.get("excluded_from_refseq", "")
                flags = {x.strip() for x in excluded.split(";") if x.strip() and x.strip() != "na"}
                if f.get("genome_rep", "Full") == "Partial":  # NCBI's own: a part of the genome was assembled
                    flags.add("partial genome")
                flags = sorted(flags)
                rows.append(dict(accession=f["assembly_accession"], taxid=int(f["taxid"]),
                                 species_taxid=int(f["species_taxid"]), organism=f.get("organism_name", ""),
                                 strain=re.sub(r"^\w+=", "", f.get("infraspecific_name", "")), isolate=f.get("isolate", ""),
                                 level=f.get("assembly_level", ""), asm_name=f.get("asm_name", ""),
                                 size=int(f["genome_size"]) if f.get("genome_size", "").isdigit() else None,
                                 contigs=int(f["contig_count"]) if f.get("contig_count", "").isdigit() else None,
                                 category="derived from metagenome" if "metagenome" in excluded else "isolate",
                                 ftp=f.get("ftp_path", ""), flags=flags, refseq_category=f.get("refseq_category", ""),
                                 paired=f.get("gbrs_paired_asm", "na"), date=f.get("seq_rel_date", ""),
                                 submitter=f.get("submitter", ""), biosample=f.get("biosample", ""),
                                 type_material=f.get("relation_to_type_material", "") or "na",  # since 2026-11-07
                                 genes=(int(f["protein_coding_gene_count"]) if f["protein_coding_gene_count"].isdigit()
                                        else 0) if "protein_coding_gene_count" in f else None,
                                 group=re.sub(r"^(genbank|refseq|assembly_summary)_", "",
                                              os.path.basename(path)).rsplit(".", 1)[0]))
    rows = [r for r in rows if r["size"] is None or r["size"] >= MIN_GENOME]
    sizes = collections.defaultdict(list)
    for r in rows:
        if r["size"]:
            sizes[r["species_taxid"]].append(r["size"])
    median = {t: sorted(s)[len(s) // 2] for t, s in sizes.items() if len(s) >= 3}
    return [r for r in rows if not r["size"] or r["species_taxid"] not in median
            or r["size"] >= PARTIAL_SHARE * median[r["species_taxid"]]]


STATE = "source.json"


def folder_of(inputs, name):
    return os.path.join(inputs, SOURCES[name]["dir"])


def state_of(inputs, name):
    """A source's source.json ({} where there is none)."""
    path = os.path.join(folder_of(inputs, name), STATE)
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def ready(inputs, name):
    """Whether a source is downloaded and its checks passed (its folder, else None)."""
    state = state_of(inputs, name)
    spec = SOURCES[name]
    return folder_of(inputs, name) if state.get("ready") and state.get("version") == spec["version"] else None


def count_records(path):
    """FASTA records (gzipped or not), or the lines of a table but its header."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        if ".tab" in path or path.endswith(".tsv"):
            return sum(1 for _ in fh) - 1
        return sum(1 for line in fh if line.startswith(">"))


def unpack(name, folder, seven_zip="7za"):
    """A source's files as rdna_freeze.py reads them: ROD's tables out of its zip, EUKARYOME's long reads out of its zip
    and 7z (gzipped) and its compromised-quality ids as a table (openpyxl, where installed; else none, said: until
    2026-11-06 a missing openpyxl left all of EUKARYOME not ready). Every file is written under this process' own
    name (util.part_of) and EUKARYOME's archives are unpacked in this process' own folder, removed whatever happens
    (since 2026-11-10: fixed `.part` names and one `unpack` folder, which two runs fetching the sources at once wrote
    and removed under each other; a run that failed left the folder)."""
    if name == "ROD":
        with zipfile.ZipFile(os.path.join(folder, SOURCES[name]["files"][0][0])) as z:
            for member in z.namelist():
                base = os.path.basename(member)
                if base in SOURCES[name]["unpack"]:
                    part = util.part_of(os.path.join(folder, base))
                    with z.open(member) as src, open(part, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                    os.replace(part, os.path.join(folder, base))
    elif name == "EUKARYOME":
        for dead in util.dead_partials(folder, "unpack"):  # the folders of runs that were killed
            shutil.rmtree(dead, ignore_errors=True)
        work = util.part_of(os.path.join(folder, "unpack"))
        shutil.rmtree(work, ignore_errors=True)
        os.makedirs(work)
        dest = os.path.join(folder, "General_EUK_longread_v2.1.fasta.gz")
        part = util.part_of(dest)
        try:
            with zipfile.ZipFile(os.path.join(folder, "General_EUK_longread_v2.1.zip")) as z:
                z.extractall(work)
            archive = os.path.join(work, "General_EUK_longread_v2.1.7z")
            if not shutil.which(seven_zip):
                raise RuntimeError(f"{seven_zip}: not found (conda-forge's p7zip), for EUKARYOME's 7z archive")
            test = subprocess.run([seven_zip, "t", archive], capture_output=True, text=True)
            if test.returncode != 0 or "Everything is Ok" not in test.stdout:
                raise RuntimeError(f"{archive}: the 7z archive does not test whole")
            subprocess.run([seven_zip, "x", "-y", f"-o{work}", archive], check=True, stdout=subprocess.DEVNULL)
            with open(os.path.join(work, "General_EUK_longread_v2.1.fasta"), "rb") as src, \
                    util.bgzf_open(part, "wb", util.WORK_LEVEL) as dst:
                shutil.copyfileobj(src, dst, 1 << 22)
            os.replace(part, dest)
        finally:
            shutil.rmtree(work, ignore_errors=True)
            if os.path.exists(part):
                os.remove(part)
        compromised = os.path.join(folder, "EUK2_1_compromised_quality.tsv")
        from .labels import read_xlsx_rows  # noqa: PLC0415
        try:
            rows = read_xlsx_rows(os.path.join(folder, "EUK2_1_compromised_quality.xlsx"))
            at = rows[0].index("accession")
            ids = [r[at] for r in rows[1:] if len(r) > at and r[at]]
        except ImportError:
            print("WARNING: no openpyxl: EUKARYOME's compromised-quality ids not read (none left out)", flush=True)
            ids = []
        with open(util.part_of(compromised), "w", newline="\n") as fh:
            fh.write("accession\n" + "".join(f"{i}\n" for i in ids))
        os.replace(util.part_of(compromised), compromised)


def rod_operons(folder):
    """ROD's variants: dicts of its table's columns (seqid, reference, size, assembly_id, taxid, length, sequence,
    lineage, ...)."""
    with gzip.open(os.path.join(folder, "ROD_v1.2_operon_variants.tab.gz"), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(head, line.rstrip("\n").split("\t")))


def eukaryome_long(folder):
    """EUKARYOME's long reads: (EUK id, label k__ to s__, sequence)."""
    with gzip.open(os.path.join(folder, "General_EUK_longread_v2.1.fasta.gz"), "rt") as fh:
        name, parts = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    acc, _, label = name.partition(";")
                    yield acc, label, "".join(parts).upper()
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
        if name is not None:
            acc, _, label = name.partition(";")
            yield acc, label, "".join(parts).upper()


def eukaryome_compromised(folder):
    """The EUK ids EUKARYOME marks of compromised quality."""
    path = os.path.join(folder, "EUK2_1_compromised_quality.tsv")
    if not os.path.isfile(path):
        return set()
    with open(path) as fh:
        next(fh)
        return {line.strip() for line in fh if line.strip()}

SOURCES = {
    "ROD": {
        "dir": "ROD", "version": "v1.2",
        "files": [("ROD-v1.2-goldenrod.zip",
                   "https://zenodo.org/api/records/12558683/files/krabberod/ROD-v1.2-goldenrod.zip/content",
                   "0bd079ecf77fac32ba9a30b2785ccde1", 338556987)],
        "unpack": ["ROD_v1.2_operon_variants.tab.gz", "ROD_v1.2.NCBI_metadata.tab.gz",
                   "ROD_v1.2_genome_statistics.tab.gz"],
        "records": ("ROD_v1.2_operon_variants.tab.gz", 69480),
    },
    "EUKARYOME": {
        "dir": "EUKARYOME2.1", "version": "2.1",
        "files": [("General_EUK_longread_v2.1.zip",
                   "https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_longread_v2.1.zip", None, 83802082),
                  ("EUK2_1_compromised_quality.xlsx",
                   "https://sisu.ut.ee/wp-content/uploads/sites/643/EUK2_1_compromised_quality.xlsx", None, 8714361),
                  ("EUKARYOME2_1_classification-2.xlsx",
                   "https://sisu.ut.ee/wp-content/uploads/sites/643/EUKARYOME2_1_classification-2.xlsx", None, 4262950)],
        "unpack": ["General_EUK_longread_v2.1.fasta.gz", "EUK2_1_compromised_quality.tsv"],
        "records": ("General_EUK_longread_v2.1.fasta.gz", 400217),
    },
}


def fetch(inputs, name, download, seven_zip="7za", log=print):
    """A source into INPUTS/<its folder>, checked: download(url, dest, md5) fetches a file (download_ref_data.download:
    resuming, until the MD5 is right; without an MD5, until the length the server announces); a file there of the right
    bytes and MD5 is kept. Then unpacked and its records counted against the published number. -> its source.json, ready
    or with why not (a failed check stops nothing else: the build then lacks that source)."""
    spec = SOURCES[name]
    folder = folder_of(inputs, name)
    os.makedirs(folder, exist_ok=True)
    state = state_of(inputs, name)
    if state.get("ready") and state.get("version") == spec["version"]:
        log(f"{name} {spec['version']}: {folder} is there, checked {state.get('checked')}")
        return state
    files, problems = {}, []
    for file, url, md5, size in spec["files"]:
        dest = os.path.join(folder, file)
        if not (os.path.isfile(dest) and os.path.getsize(dest) == size and (md5 is None or util.md5(dest) == md5)):
            log(f"{name} {spec['version']}: fetching {file} ({size / 1e6:.0f} MB)")
            try:
                download(url, dest)  # its MD5 checked below, a problem of the source's, not a stop
            except OSError as e:  # the site moved or down: this source not ready (until 2026-11-06 the stage stopped)
                problems.append(f"{file}: {e}")
                continue
        got = util.md5(dest)
        files[file] = {"bytes": os.path.getsize(dest), "md5": got}
        if os.path.getsize(dest) != size:
            problems.append(f"{file}: {os.path.getsize(dest)} bytes, {size} published")
        if md5 is not None and got != md5:
            problems.append(f"{file}: MD5 {got}, {md5} published")
    if not problems:
        try:
            unpack(name, folder, seven_zip)
        except (RuntimeError, OSError, ImportError, zipfile.BadZipFile, subprocess.CalledProcessError) as e:
            problems.append(f"unpacking: {e}")
    counted = None
    if not problems:
        file, published = spec["records"]
        counted = count_records(os.path.join(folder, file))
        if counted != published:
            problems.append(f"{file}: {counted} records, {published} published")
    import datetime  # noqa: PLC0415
    state = {"source": name, "version": spec["version"], "files": files, "records": counted,
             "ready": not problems, "problems": problems,
             "checked": datetime.datetime.now().isoformat(timespec="seconds")}
    with open(util.part_of(os.path.join(folder, STATE)), "w") as fh:
        json.dump(state, fh, indent=1)
    os.replace(util.part_of(os.path.join(folder, STATE)), os.path.join(folder, STATE))
    log(f"{name} {spec['version']}: " + ("ready, " + f"{counted} records as published" if not problems else
                                          "NOT READY: " + "; ".join(problems)))
    return state
