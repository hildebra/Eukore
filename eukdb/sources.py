"""EUKDB's sources (since 2026-11-05; copied from protal's scripts/mini_db/eukaryotes.py, rdna_sources.py and
download_ref_data.py, so that EUKDBbuild imports nothing from them): EukDetect2's genomic euk30 copies (Zenodo
19056625) and EukDetect1's coding ones (figshare 12670856), the loci table and its EukDetect1 locus map, NCBI's GenBank
and RefSeq assembly tables, and the published eukaryote rDNA operons (EUKARYOME 2.1's long reads, ROD v1.2).
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
LOCI_TABLE = os.path.join(DATA, "euk30_loci.tsv")
PLACED_DIR = os.path.join(DATA, "placed")      # the copies placed so far (eukdetect1_genomic.py's, then placement.py's)
PLACED_FASTA = "ncbi_eukprot_met_arch_markers.genomic.fna.zst"
PLACED_TABLE = "eukdetect1_genomic.tsv.zst"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq"
ASSEMBLY_GROUPS = ("fungi", "protozoa", "plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other")


def is_placeholder(ncbi, taxid):
    return ncbi.is_placeholder(taxid)


def download(url, dest, md5=None):
    """A file fetched whole (util.fetch_file) and, with an md5, checked -> raises where it did not arrive."""
    why = util.fetch_file(url, dest, gzip_check=False)
    if why:
        raise OSError(f"{url}: {why}")
    if md5 and util.md5(dest) != md5:
        raise OSError(f"{dest}: MD5 {util.md5(dest)}, {md5} expected")


LOCI = 30


MAX_LOCUS = 16000


EUKDETECT_FILES = ("eukdb.fasta", "busco_taxid_genome_link.txt", "taxdump.tar.gz")


LOCUS_RE = re.compile(r"-(\d+at\d+)-S\d+$")


LEVELS = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}


def read_link(path):
    """busco_taxid_genome_link.txt -> ({genome: (taxid, group)}, {record: (genome, locus)}); group Fungi or Protist."""
    genomes, records = {}, {}
    with open_text(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            record, taxid, genome = fields[:3]
            m = LOCUS_RE.search(record)
            if not m:
                continue
            genomes[genome] = (int(taxid), "Fungi" if record.startswith("Fungi") else "Protist")
            records[record] = (genome, m.group(1))
    return genomes, records


def read_loci(path, records, loci=None, lengths_only=False):
    """eukdb.fasta -> {genome: {locus: sequence}}, a genome's longest sequence of a locus it has in several (-S1, -S2),
    of the loci asked for (all without `loci`); with lengths_only the lengths instead of the sequences."""
    result = collections.defaultdict(dict)
    keep = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                keep = records.get(line[1:].split(None, 1)[0].strip())
                if keep is not None and loci is not None and keep[1] not in loci:
                    keep = None
                continue
            if keep is None:
                continue
            seq = line.strip()
            genome, locus = keep
            old = result[genome].get(locus)
            old_length = (old if lengths_only else len(old)) if old is not None else -1
            if len(seq) > old_length:
                result[genome][locus] = len(seq) if lengths_only else seq.upper().replace("U", "T")
            keep = None  # one line per sequence
    return result


def read_loci_table(path=LOCI_TABLE):
    """euk30_loci.tsv -> its loci, in its order."""
    with open(path) as fh:
        return [line.split("\t")[0] for line in fh if line.strip() and not line.startswith(("#", "locus\t"))]


def read_locus_map(path=LOCI_TABLE):
    """euk30_loci.tsv -> {locus: its EukDetect1 locus} (the column eukdetect1_locus; "-" or none: no such locus)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "eukdetect1_locus" not in header:
            return {}
        col = header.index("eukdetect1_locus")
        rows = (line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#"))
        return {f[0]: f[col] for f in rows if len(f) > col and f[col] not in ("", "-")}


# EukDetect1 (Lind & Pollard 2021, Microbiome 9:58; figshare 10.6084/m9.figshare.12670856.v7): BUSCO eukaryota_odb10's
# genes of 3,713 eukaryotes, the non-vertebrate animals and the green algae among them, which EukDetect2 has none of.
# Its records are named as EukDetect2's (<group>-<genome>-<locus>-<copy>), its loci are another OrthoDB's (odb10):
# map_loci finds which of them each euk30 locus is. A copy S1 is single; D1, D2, ... a BUSCO found duplicated (left
# out); a record of several species' identical copies (CD-HIT, "...Collapse...") names no species (left out).
EUKDETECT1_FILES = ("busco_taxid_link.txt", "ncbi_eukprot_met_arch_markers.fna")


# EukDetect1's copies are BUSCO's coding sequences, spliced (no introns), where a read of the genome holds them:
# eukdetect1_genomic.py writes the copies it places in their genomes, introns and all, into this file beside them.
# A locus takes its genomic copy where there is one that clean_copy keeps, else its coding one (a species keeps its loci).
EUKDETECT1_GENOMIC = "ncbi_eukprot_met_arch_markers.genomic.fna"


EUKDETECT1_GROUPS = {"fungi": "Fungi", "protist": "Protist", "insect": "Metazoa", "worm": "Metazoa",
                     "metazoa": "Metazoa", "greenalgae": "Archaeplastida", "archaeplast": "Archaeplastida",
                     "archaeplasta": "Archaeplastida"}


# A EukDetect1 species with fewer of the euk30 loci EukDetect1 has (17 of 30) is left out: its reads would rarely reach
# enough loci to tell it from its relatives (EukDetect2's species have 24.8 of the 30 on average).
MIN_EUKDETECT1_LOCI = 6


RECORD_RE = re.compile(r"^(?P<group>[^-]+)-(?P<genome>.+)-(?P<locus>\d+at\d+)-(?P<copy>[SD]\d+)$")


def read_link1(path, locus_map):
    """EukDetect1's busco_taxid_link.txt -> ({genome: (taxid, group)}, {record: (genome, euk30 locus)}) of the single
    copies of the loci locus_map ({euk30 locus: EukDetect1 locus}) maps, a locus a genome has duplicated left out."""
    euk30_of = {v: k for k, v in locus_map.items()}
    genomes, single, duplicated = {}, {}, set()
    with open_text(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            m = RECORD_RE.match(fields[0]) if len(fields) >= 2 else None
            if not m or m.group("group") not in EUKDETECT1_GROUPS:
                continue
            genome = m.group("genome")
            genomes[genome] = (int(fields[1]), EUKDETECT1_GROUPS[m.group("group")])
            locus = euk30_of.get(m.group("locus"))
            if locus is None:
                continue
            if m.group("copy") == "S1":
                single[fields[0]] = (genome, locus)
            else:
                duplicated.add((genome, locus))
    return genomes, {r: gl for r, gl in single.items() if gl not in duplicated}


def genomic_copies(eukdetect1, records1, placed_dir=PLACED_DIR):
    """EukDetect1's genomic copies of the records read_link1 gave ({record: (genome, locus)}): placed_dir's (the
    checkout's), then those beside EukDetect1's files, a genome's there in place of its copies in the checkout.
    -> {genome: {locus: copy}}"""
    result = {}
    for path in (os.path.join(placed_dir, PLACED_FASTA) if placed_dir else "", os.path.join(eukdetect1, EUKDETECT1_GENOMIC)):
        if path and os.path.isfile(path):
            result.update(read_loci(path, records1))
    return result


def species_loci(eukdetect, eukdetect1, locus_map, ncbi, loci):
    """The euk30 loci each NCBI species has in euk30's sources (EukDetect2's genomes, EukDetect1's single copies of the
    loci locus_map maps; a download without EukDetect1: EukDetect2's alone) -> ({species taxid: {loci}}, EukDetect2's
    (genomes, records))."""
    genomes2, records2 = read_link(os.path.join(eukdetect, "busco_taxid_genome_link.txt"))
    link1 = os.path.join(eukdetect1 or "", "busco_taxid_link.txt")
    genomes1, records1 = read_link1(link1, locus_map) if locus_map and os.path.isfile(link1) else ({}, {})
    have = collections.defaultdict(set)
    for genomes, records in ((genomes2, records2), (genomes1, records1)):
        for genome, locus in records.values():
            if locus in loci:
                have[ncbi.species_of(genomes[genome][0])].add(locus)
    return have, (genomes2, records2)


# No eukaryote's nuclear genome is smaller (microsporidia's, the smallest, are ~2.2 Mb): GenBank has records of 500-600
# bases as "complete genomes" of Aspergillus flavus, niger, fumigatus and others (a gene's sequence), which a build
# simulated as cells of 600 bases, and which were their species' best-ranked match. Nor is a species' genome a third of
# its other assemblies' median (an organelle's, a partial assembly).
MIN_GENOME = 2_000_000


PARTIAL_SHARE = 0.3


def read_assemblies(paths):
    """GenBank's assembly_summary.txt files -> [dict] of the latest assemblies (accession, taxid, species_taxid,
    organism, strain, isolate, level, asm_name, size, contigs, category, ftp), without those smaller than MIN_GENOME and,
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
                rows.append(dict(accession=f["assembly_accession"], taxid=int(f["taxid"]),
                                 species_taxid=int(f["species_taxid"]), organism=f.get("organism_name", ""),
                                 strain=re.sub(r"^\w+=", "", f.get("infraspecific_name", "")), isolate=f.get("isolate", ""),
                                 level=f.get("assembly_level", ""), asm_name=f.get("asm_name", ""),
                                 size=int(f["genome_size"]) if f.get("genome_size", "").isdigit() else None,
                                 contigs=int(f["contig_count"]) if f.get("contig_count", "").isdigit() else None,
                                 category="derived from metagenome" if "metagenome" in excluded else "isolate",
                                 ftp=f.get("ftp_path", "")))
    rows = [r for r in rows if r["size"] is None or r["size"] >= MIN_GENOME]
    sizes = collections.defaultdict(list)
    for r in rows:
        if r["size"]:
            sizes[r["species_taxid"]].append(r["size"])
    median = {t: sorted(s)[len(s) // 2] for t, s in sizes.items() if len(s) >= 3}
    return [r for r in rows if not r["size"] or r["species_taxid"] not in median
            or r["size"] >= PARTIAL_SHARE * median[r["species_taxid"]]]


def norm(text):
    return re.sub(r"[^a-z0-9]", "", text.lower())


def assembly_rank(row):
    """Better assemblies first: complete genomes, then chromosomes, scaffolds, contigs; then fewer contigs."""
    return (LEVELS.get(row["level"], 4), row["contigs"] if row["contigs"] is not None else 1 << 30, row["accession"])


def match_genomes(genomes, assemblies, ncbi):
    """Each EukDetect2 genome's GenBank assembly: by its name (organism and strain or isolate, letters and digits alone)
    among its taxid's species' assemblies, else among those of the name, else its taxid's best assembly (not for an
    NCBI placeholder's taxid). -> {genome: (assembly row or None, how)}"""
    by_species, by_name = collections.defaultdict(list), collections.defaultdict(list)
    for row in assemblies:
        by_species[row["species_taxid"]].append(row)
        for key in {norm(row["organism"]), norm(row["organism"] + row["strain"]), norm(row["organism"] + row["isolate"]),
                    norm(row["organism"] + row["strain"] + row["isolate"])}:
            by_name[key].append(row)
    result = {}
    for genome, (taxid, _group) in genomes.items():
        species = ncbi.species_of(taxid)
        candidates = by_species.get(species, [])
        named = [r for r in by_name.get(norm(genome), [])]
        both = [r for r in named if r in candidates]
        for rows, how in ((both, "name and taxid"), (named, "name"),
                          ([] if ncbi.is_placeholder(species) else candidates, "taxid")):
            if rows:
                result[genome] = (min(rows, key=assembly_rank), how + (" (one)" if len(rows) == 1 else " (several)"))
                break
        else:
            result[genome] = (None, "none")
    return result


def own_accessions(matched, species_of, rep_of):
    """match_genomes' matches with each GenBank assembly a genome's alone (since 2026-11-03): of the genomes matched to
    one assembly the one matched by its name keeps it, else the species' representative, else the first by name; the
    others are matched "own accession" without an assembly, and so are written under an accession of their own
    (EUK_<genome>, as a genome without an assembly). The taxid fallback gives each of a species' genomes its name finds
    no assembly of the species' best one: r226's euk30 wrote up to 41 copies of a locus of Schizophyllum commune's
    genomes as GB_GCA_000143185.2, of which the converter kept the first (5,005 copies of 107 species missing from the
    full reference, 858 loci under a representative's accession another genome's copy; now 298 genomes of 109 species
    have their own). -> {genome: (assembly row or None, how)}"""
    claims = collections.defaultdict(list)
    for genome, (row, _how) in matched.items():
        if row:
            claims[row["accession"]].append(genome)
    result = dict(matched)
    for genomes in claims.values():
        keep = min(genomes, key=lambda g: (not matched[g][1].startswith("name"), rep_of.get(species_of[g]) != g, g))
        for genome in genomes:
            if genome != keep:
                result[genome] = (None, f"own accession ({matched[keep][0]['accession']} is {keep}'s)")
    return result


def base_name(name):
    """A name without PR2's "(Fungi)" (Cryptococcus(Fungi)) and NCBI's brackets ([Candida] boidinii)."""
    return re.sub(r"\([^()]*\)", "", name).replace("[", "").replace("]", "")


STATE = "source.json"


def folder_of(inputs, name):
    return os.path.join(inputs, SOURCES[name]["dir"])


def md5_of(path):
    import hashlib  # noqa: PLC0415
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


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
    and 7z (gzipped) and its compromised-quality ids as a table (openpyxl, where installed; else none)."""
    if name == "ROD":
        with zipfile.ZipFile(os.path.join(folder, SOURCES[name]["files"][0][0])) as z:
            for member in z.namelist():
                base = os.path.basename(member)
                if base in SOURCES[name]["unpack"]:
                    with z.open(member) as src, open(os.path.join(folder, base) + ".part", "wb") as dst:
                        shutil.copyfileobj(src, dst)
                    os.replace(os.path.join(folder, base) + ".part", os.path.join(folder, base))
    elif name == "EUKARYOME":
        work = os.path.join(folder, "unpack")
        os.makedirs(work, exist_ok=True)
        with zipfile.ZipFile(os.path.join(folder, "General_EUK_longread_v2.1.zip")) as z:
            z.extractall(work)
        archive = os.path.join(work, "General_EUK_longread_v2.1.7z")
        if not shutil.which(seven_zip):
            raise RuntimeError(f"{seven_zip}: not found (conda-forge's p7zip), for EUKARYOME's 7z archive")
        test = subprocess.run([seven_zip, "t", archive], capture_output=True, text=True)
        if test.returncode != 0 or "Everything is Ok" not in test.stdout:
            raise RuntimeError(f"{archive}: the 7z archive does not test whole")
        subprocess.run([seven_zip, "x", "-y", f"-o{work}", archive], check=True, stdout=subprocess.DEVNULL)
        dest = os.path.join(folder, "General_EUK_longread_v2.1.fasta.gz")
        with open(os.path.join(work, "General_EUK_longread_v2.1.fasta"), "rb") as src, \
                gzip.open(dest + ".part", "wb", compresslevel=6) as dst:
            shutil.copyfileobj(src, dst, 1 << 22)
        os.replace(dest + ".part", dest)
        shutil.rmtree(work)
        compromised = os.path.join(folder, "EUK2_1_compromised_quality.tsv")
        ids = read_xlsx_column(os.path.join(folder, "EUK2_1_compromised_quality.xlsx"), "accession")
        with open(compromised + ".part", "w", newline="\n") as fh:
            fh.write("accession\n" + "".join(f"{i}\n" for i in ids))
        os.replace(compromised + ".part", compromised)


def read_xlsx_column(path, column):
    """A column of an xlsx file's first sheet (openpyxl)."""
    import openpyxl  # noqa: PLC0415  (conda-forge's openpyxl, beside pandas)
    sheet = openpyxl.load_workbook(path, read_only=True).worksheets[0]
    rows = sheet.iter_rows(values_only=True)
    head = list(next(rows))
    at = head.index(column)
    return [str(r[at]) for r in rows if r and r[at]]


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


# EukDetect2's files (Zenodo record 19056625, CC BY 4.0) and their MD5s, into OUT/eukdetect2 (2.1 GB), and NCBI's
# GenBank assembly tables of the groups its genomes are in, into OUT/ncbi_assemblies (~16 MB).
EUKDETECT_MIRROR = "https://zenodo.org/api/records/19056625/files"


EUKDETECT_MD5 = {"eukdb.fasta": "91f3883f231ead6853f6883449cd331f",
                 "busco_taxid_genome_link.txt": "3eae92e5928da5af74eb262937ab96b8",
                 "taxdump.tar.gz": "5b89e986aa5dee81b9d54e6044d515ac"}


# EukDetect1's database (figshare 10.6084/m9.figshare.12670856.v7, 1.7 GB): its species EukDetect2 lacks (non-vertebrate
# animals, green algae, more fungi and protists) join euk30 (since 2026-10-26); the two files eukaryotes.py reads are
# kept in OUT/eukdetect1 (0.96 GB), the archive is removed.
EUKDETECT1_URL = "https://ndownloader.figshare.com/files/26173346"


EUKDETECT1_ARCHIVE = ("ncbi_eukprot_met_arch_database.tar.gz", "d2cd9ed430b63a1e1c52cb9866411b70")


ASSEMBLIES_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank"


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
        if not (os.path.isfile(dest) and os.path.getsize(dest) == size and (md5 is None or md5_of(dest) == md5)):
            log(f"{name} {spec['version']}: fetching {file} ({size / 1e6:.0f} MB)")
            download(url, dest, md5)
        got = md5_of(dest)
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
    with open(os.path.join(folder, STATE) + ".part", "w") as fh:
        json.dump(state, fh, indent=1)
    os.replace(os.path.join(folder, STATE) + ".part", os.path.join(folder, STATE))
    log(f"{name} {spec['version']}: " + ("ready, " + f"{counted} records as published" if not problems else
                                          "NOT READY: " + "; ".join(problems)))
    return state


