"""BUSCO's eukaryota_odb12 dataset (since 2026-11-05): its BUSCOs' proteins, HMMs and score cutoffs, from which
Eukore calls the loci in every genome itself (loci.py), in place of EukDetect1's and EukDetect2's copies.

The dataset (busco-data.ezlab.org; OrthoDB 12.1, 129 BUSCOs of 456 species, 2026-05-22) holds per BUSCO:
- refseq_db.faa.gz: its orthologs' proteins in the 456 species (~770 a BUSCO; records <busco>_<taxid>_<id>);
- ancestral_variants: 10 reconstructed ancestral proteins;
- hmms/<busco>.hmm and scores_cutoff: the profile and the bit score a protein reaches to be the BUSCO (BUSCO's rule).

A genome's queries are, per locus, the proteins of the NEAREST species to it in NCBI's taxonomy (the most shared
ancestors; then by name) and the ancestral variants: a close relative's protein places introns and exon ends best, the
ancestral ones reach where no relative is close.
"""

import collections
import gzip
import json
import os
import shutil
import tarfile

from . import util

URL = "https://busco-data.ezlab.org/v5/data"
LINEAGE = "eukaryota_odb12"
NEAREST = 30
STATE = "odb.json"


def latest(versions_path, lineage=LINEAGE):
    """BUSCO's file_versions.tsv -> (date, md5) of a lineage's newest dataset."""
    best = None
    with open(versions_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == lineage and (best is None or f[1] > best[0]):
                best = (f[1], f[2])
    if best is None:
        raise ValueError(f"{versions_path}: no {lineage}")
    return best


def folder_of(src_dir, lineage=LINEAGE):
    return os.path.join(src_dir, lineage)


def state_of(src_dir, lineage=LINEAGE):
    try:
        with open(os.path.join(src_dir, f"{lineage}.{STATE}")) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def fetch(src_dir, download, log=print, lineage=LINEAGE):
    """The newest dataset unpacked into src_dir/<lineage> (kept where its state names the same version) -> its folder."""
    os.makedirs(src_dir, exist_ok=True)
    versions = os.path.join(src_dir, "busco_file_versions.tsv")
    download(f"{URL}/file_versions.tsv", versions)
    date, md5 = latest(versions, lineage)
    folder = folder_of(src_dir, lineage)
    state = state_of(src_dir, lineage)
    if state.get("date") == date and os.path.isfile(os.path.join(folder, "scores_cutoff")):
        log(f"{lineage} {date}: {folder} is there")
        return folder
    archive = os.path.join(src_dir, f"{lineage}.{date}.tar.gz")
    log(f"{lineage}: fetching {date}")
    download(f"{URL}/lineages/{lineage}.{date}.tar.gz", archive, md5)
    shutil.rmtree(folder, ignore_errors=True)
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(src_dir, filter="data")
    os.remove(archive)
    with open(os.path.join(src_dir, f"{lineage}.{STATE}"), "w", newline="\n") as fh:
        json.dump({"lineage": lineage, "date": date, "md5": md5}, fh)
    return folder


def hmm_header(path):
    """An HMM file's (NAME, LENG): BUSCO names a profile by its OrthoDB group without the "at<taxid>" (10248812)."""
    name, length = None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("NAME"):
                name = line.split()[1]
            elif line.startswith("LENG"):
                length = int(line.split()[1])
            elif line.startswith("HMM "):
                break
    return name, length


def all_loci(folder):
    """Every BUSCO of a dataset (its scores_cutoff), in its order."""
    with open(os.path.join(folder, "scores_cutoff")) as fh:
        return [line.split()[0] for line in fh if line.strip()]


def names_of(folder):
    """{BUSCO: OrthoDB's name of it} (links_to_ODB12.txt; {} without one)."""
    path = os.path.join(folder, "links_to_ODB12.txt")
    out = {}
    if os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > 1:
                    out[f[0]] = f[1]
    return out


class Odb:
    """A dataset's loci: proteins (refseq_db and ancestral variants), HMMs, cutoffs. loci: the BUSCOs to keep."""

    def __init__(self, folder, loci):
        self.folder = folder
        self.loci = list(loci)
        wanted = set(self.loci)
        self.cutoff = {}
        with open(os.path.join(folder, "scores_cutoff")) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2 and f[0] in wanted:
                    self.cutoff[f[0]] = float(f[1])
        missing = wanted - set(self.cutoff)
        if missing:
            raise ValueError(f"{folder}: no {sorted(missing)[:5]} (another dataset?)")
        self.hmm = {locus: os.path.join(folder, "hmms", f"{locus}.hmm") for locus in self.loci}
        self.names = names_of(folder)
        headers = {locus: hmm_header(path) for locus, path in self.hmm.items()}
        self.length = {locus: h[1] for locus, h in headers.items()}
        self.locus_of_hmm = {h[0] or locus: locus for locus, h in headers.items()}  # hmmsearch's names -> loci
        self.proteins = collections.defaultdict(list)  # locus -> [(taxid, name, protein)]
        for name, seq in read_fasta(os.path.join(folder, "refseq_db.faa.gz")):
            locus, _, rest = name.partition("_")
            taxid = rest.split("_", 1)[0]
            if locus in wanted and taxid.isdigit():
                self.proteins[locus].append((int(taxid), name, seq))
        self.ancestral = collections.defaultdict(list)
        path = os.path.join(folder, "ancestral_variants")
        if os.path.isfile(path):
            for name, seq in read_fasta(path):
                locus = name.rsplit("_", 1)[0]
                if locus in wanted:
                    self.ancestral[locus].append((name, seq))
        self.median = {locus: sorted(len(s) for _t, _n, s in ps)[len(ps) // 2] if ps else 0
                       for locus, ps in self.proteins.items()}

    def species(self):
        return sorted({t for ps in self.proteins.values() for t, _n, _s in ps})

    def hmm_db(self, path):
        """The loci's HMMs in one file (hmmsearch's query) -> path."""
        if not os.path.isfile(path):
            with open(util.part_of(path), "w", newline="\n") as out:
                for locus in self.loci:
                    with open(self.hmm[locus]) as fh:
                        out.write(fh.read())
            os.replace(util.part_of(path), path)
        return path

    def queries(self, lineage, lineages, nearest=NEAREST):
        """A genome's queries -> {query name: protein}: per locus the proteins of the `nearest` species closest to it
        (lineage: the genome's ancestor taxids; lineages: {dataset species taxid: its ancestor taxids}), then the
        ancestral variants. Query names are <locus>~<record>, so that a hit names its locus."""
        mine = set(lineage)
        depth = {t: len(mine & set(anc)) for t, anc in lineages.items()}
        out = {}
        for locus in self.loci:
            ranked = sorted(self.proteins.get(locus, ()), key=lambda p: (-depth.get(p[0], 0), p[1]))
            for _t, name, seq in ranked[:nearest]:
                out[f"{locus}~{name}"] = seq
            for name, seq in self.ancestral.get(locus, ()):
                out[f"{locus}~{name}"] = seq
        return out


def read_fasta(path):
    """(name, sequence) of a FASTA file, plain or gzipped; names to their first space."""
    opener = gzip.open if path.endswith(".gz") else open
    name, parts = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name, parts = line[1:].split(None, 1)[0] if line[1:].strip() else "", []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


def species_lineages(odb, ncbi):
    """{dataset species taxid: its NCBI ancestor taxids} (a genome's queries' ranking; workers get it as data)."""
    return {t: list(ncbi.ancestors(ncbi.merged.get(t, t))) for t in odb.species()}
