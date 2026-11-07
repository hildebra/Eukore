"""BUSCO's eukaryota_odb12 dataset (since 2026-11-05): its BUSCOs' proteins, HMMs and score cutoffs, from which
Eukore calls the loci in every genome itself (loci.py), in place of EukDetect1's and EukDetect2's copies.

The dataset (busco-data.ezlab.org; OrthoDB 12.1, 129 BUSCOs of 456 species, 2026-05-22) holds per BUSCO:
- refseq_db.faa.gz: its orthologs' proteins in the 456 species (~770 a BUSCO; records <busco>_<taxid>_<id>);
- ancestral_variants: 10 reconstructed ancestral proteins;
- hmms/<busco>.hmm and scores_cutoff: the profile and the bit score a protein reaches to be the BUSCO (BUSCO's rule).

A genome's queries are, per locus, the NEAREST proteins to it: the dataset's proteins ranked by their species' nearness
in NCBI's taxonomy (the most shared ancestors; then by name), the first NEAREST of them (a species with paralogues
gives several: ~15-30 species), and the ancestral variants. A close relative's protein places introns and exon ends
best, the ancestral ones reach where no relative is close.
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


def hmm_names(path):
    """Every profile's NAME in an HMM file (a marker locus may hold two: markers.py)."""
    with open(path) as fh:
        return [line.split()[1] for line in fh if line.startswith("NAME ")]


class Odb:
    """A dataset's loci: proteins (refseq_db and ancestral variants), HMMs, cutoffs. loci: the BUSCOs to keep; extra:
    further datasets in the same layout holding the rest of loci (since 2026-11-07: Eukore's markers, markers.py),
    whose decoy profiles (decoys/*.hmm) join the HMMs a candidate is scored against but are no locus."""

    def __init__(self, folder, loci, extra=()):
        self.folder = folder
        self.loci = list(loci)
        wanted = set(self.loci)
        self.cutoff, self.home, self.names = {}, {}, {}
        folders = [folder] + [f for f in extra if f]
        for f in folders:
            with open(os.path.join(f, "scores_cutoff")) as fh:
                for line in fh:
                    x = line.split()
                    if len(x) >= 2 and x[0] in wanted and x[0] not in self.cutoff:
                        self.cutoff[x[0]], self.home[x[0]] = float(x[1]), f
            self.names.update({k: v for k, v in names_of(f).items() if k not in self.names})
        missing = wanted - set(self.cutoff)
        if missing:
            raise ValueError(f"{folder}: no {sorted(missing)[:5]} (another dataset?)")
        self.markers = {locus for locus in self.loci if self.home[locus] != folder}  # Eukore's own loci
        self.hmm = {locus: os.path.join(self.home[locus], "hmms", f"{locus}.hmm") for locus in self.loci}
        headers = {locus: hmm_header(path) for locus, path in self.hmm.items()}
        self.length = {locus: h[1] for locus, h in headers.items()}
        self.locus_of_hmm = {h[0] or locus: locus for locus, h in headers.items()}  # hmmsearch's names -> loci
        for locus in self.markers:
            self.locus_of_hmm.update({name: locus for name in hmm_names(self.hmm[locus])})
        self.span = {}  # a marker profile's span share (hmm_spans): KOfam's are longer than the proteins they model
        for f in folders[1:]:
            path = os.path.join(f, "hmm_spans")
            if os.path.isfile(path):
                with open(path) as fh:
                    self.span.update({x[0]: float(x[1]) for x in (line.split() for line in fh) if len(x) >= 2})
        self.decoys = sorted(os.path.join(f, "decoys", n) for f in folders[1:]
                             if os.path.isdir(os.path.join(f, "decoys"))
                             for n in os.listdir(os.path.join(f, "decoys")) if n.endswith(".hmm"))
        self.proteins = collections.defaultdict(list)  # locus -> [(taxid, name, protein)]
        self.ancestral = collections.defaultdict(list)
        for f in folders:
            mine = {locus for locus in wanted if self.home[locus] == f}
            for name, seq in read_fasta(os.path.join(f, "refseq_db.faa.gz")):
                locus, _, rest = name.partition("_")
                taxid = rest.split("_", 1)[0]
                if locus in mine and taxid.isdigit():
                    self.proteins[locus].append((int(taxid), name, seq))
            path = os.path.join(f, "ancestral_variants")
            if os.path.isfile(path):
                for name, seq in read_fasta(path):
                    locus = name.rsplit("_", 1)[0]
                    if locus in mine:
                        self.ancestral[locus].append((name, seq))

    def by_locus(self, hits):
        """A protein's hmmsearch hits {profile: (score, HMM coverage)} -> {locus: (score, coverage)}, per locus its best
        profile's (a marker may have two), a decoy's under its own name; a marker profile's coverage over its span
        share (its members' median coverage, hmm_spans), at most 1: complete reads as it does for a BUSCO."""
        got = {}
        for h, (s, cov) in hits.items():
            locus = self.locus_of_hmm.get(h, h)
            if locus not in got or s > got[locus][0]:
                got[locus] = (s, min(1.0, cov / self.span.get(h, 1.0)))
        return got

    def species(self):
        """The BUSCO dataset's species (its loci's proteins'; the markers' UniProt species left out): the OrthoDB
        benchmark's and the held-out calls'."""
        return sorted({t for locus, ps in self.proteins.items() if locus not in self.markers for t, _n, _s in ps})

    def all_species(self):
        """Every species whose proteins are queries (the markers' too): their lineages rank a genome's queries."""
        return sorted({t for ps in self.proteins.values() for t, _n, _s in ps})

    def hmm_db(self, path):
        """The loci's HMMs and the markers' decoys in one file (hmmsearch's query) -> path."""
        if not os.path.isfile(path):
            with open(util.part_of(path), "w", newline="\n") as out:
                for p in [self.hmm[locus] for locus in self.loci] + self.decoys:
                    with open(p) as fh:
                        out.write(fh.read())
            os.replace(util.part_of(path), path)
        return path

    def digest(self):
        """A short digest of the loci and decoys (an hmm_db's name: a new locus makes a new one)."""
        import hashlib  # noqa: PLC0415
        return hashlib.md5(",".join(self.loci + [os.path.basename(d) for d in self.decoys]).encode()).hexdigest()[:10]

    def queries(self, lineage, lineages, nearest=NEAREST, leave_out=()):
        """A genome's queries -> {query name: protein}: per locus its `nearest` proteins of the species closest to it
        (lineage: the genome's ancestor taxids; lineages: {dataset species taxid: its ancestor taxids}), then the
        ancestral variants; leave_out: dataset taxids whose proteins are not taken (the genome's own species', for the
        OrthoDB benchmark's held-out calls). Query names are <locus>~<record>, so that a hit names its locus."""
        mine = set(lineage)
        depth = {t: len(mine & set(anc)) for t, anc in lineages.items()}
        leave_out = set(leave_out)
        out = {}
        for locus in self.loci:
            ranked = sorted((p for p in self.proteins.get(locus, ()) if p[0] not in leave_out),
                            key=lambda p: (-depth.get(p[0], 0), p[1]))
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
    """{dataset species taxid: its NCBI ancestor taxids} (a genome's queries' ranking; workers get it as data; the
    markers' species too, since 2026-11-07)."""
    return {t: list(ncbi.ancestors(ncbi.merged.get(t, t))) for t in odb.all_species()}


def own_taxids(odb, ncbi):
    """{NCBI species taxid: the dataset taxids of that species (itself, its strains)}: the proteins a genome of the
    species leaves out of its held-out calls."""
    out = collections.defaultdict(list)
    for t in odb.species():
        out[ncbi.species_of(t)].append(t)
    return dict(out)
