"""Each species' lineage in the eukaryote 18S references beside its NCBI taxid (since 2026-11-05, format 4): PR2's,
EUKARYOME's and KSGP's, so that a user of the package (protal) names a species as the reference it reads its 18S with
names it, without matching names itself. Also NCBI's own lineage and the species' supergroup.

A species is matched by name, in each reference alone:
- `species`: the reference names the species by NCBI's scientific name (case, "_" for " " and a genus's "(Fungi)"
  aside); `synonym`: by one of NCBI's synonyms (Candida sphaerica, an anamorph of Kluyveromyces lactis);
- `genus`: by neither, but the reference names its genus: the lineage down to the genus, no species;
- `-`: none.
A name the reference files under several lineages (a homonym of another kingdom, a misfiled record) takes the lineage
that shares the most names with NCBI's lineage of the species, then the one of the most records.

The references:
- PR2 (pr2_version_<V>_SSU_mothur.tax.gz): its 9 ranks, domain;supergroup;division;subdivision;class;order;family;
  genus;species, as PR2 writes them (Saccharomyces_cerevisiae); records of an organelle (PR2's ":plas", ":mito",
  ":apic", ":nucl" names) left out.
- EUKARYOME (its classification table, superkingdom to genus, and its long reads' species): superkingdom;kingdom;
  subkingdom;phylum;subphylum;class;order;family;genus;species, as EUKARYOME writes them.
- KSGP (KSGP_v<V>_plus2.tax: k__ (PR2's supergroup) to s__): its eukaryote nuclear records (not k__Bacteria or
  k__Archaea, no chloroplast or "plas"/"mito" names, as protal's IsPlastidLabel/IsMitochondrionLabel), as KSGP writes
  them. Its accession-named placeholders (KSGP's LCA labels, a genome's 18S) name no species.
"""

import collections
import json
import os
import re

from . import util

# --- NCBI --------------------------------------------------------------------------------------------------------

NCBI_RANKS = ("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")


def supergroup(ncbi, taxid):
    """A species' supergroup as busco/loci.tsv splits them (groups.group_of)."""
    from . import groups  # noqa: PLC0415
    return groups.group_of(ncbi, taxid)


def ncbi_lineage(ncbi, taxid):
    """NCBI's lineage of a taxon at NCBI_RANKS (domain is NCBI's superkingdom before 2025), ";"-joined, a rank NCBI does
    not name empty."""
    ranks = ncbi.ranks(taxid)
    ranks.setdefault("domain", ranks.get("superkingdom"))
    return ";".join(ranks.get(r) or "" for r in NCBI_RANKS)


# --- names -------------------------------------------------------------------------------------------------------

UNNAMED = re.compile(r"(^|\s)(sp|spp|cf|aff|nr|unclassified|unidentified|uncultured|environmental|incertae|"
                     r"sedis|\?)(\.|\s|$)", re.IGNORECASE)
PLACEHOLDER = re.compile(r"\d|~|%|^[A-Z]{1,2}_?\d|^GC[AF][_ ]|^EUK\d|^OBEP|_X+$|\.(fam|ord|cl|gen)\.", re.IGNORECASE)


def key(name):
    """A name as references are compared: lower case, "_" as " ", a "(Fungi)" after it and spaces collapsed."""
    name = re.sub(r"\(.*?\)", " ", (name or "").replace("_", " "))
    return " ".join(name.lower().split())


def is_species_name(name):
    """Whether a reference's species name names a species (not Saccharomyces_sp., an accession, Genus_X)."""
    words = (name or "").replace("_", " ").split()
    return len(words) >= 2 and words[1][:1].islower() and not UNNAMED.search(" ".join(words)) \
        and not PLACEHOLDER.search(words[0]) and not any(ch.isdigit() for ch in words[1])


def is_genus_name(name):
    return bool(name) and not UNNAMED.search(name.replace("_", " ")) and not PLACEHOLDER.search(name) \
        and name[:1].isupper() and " " not in name.strip()


# the kingdoms a genus name can be a homonym across, by the names the references and NCBI give them
KINGDOMS = {"metazoa": "animals", "fungi": "fungi", "viridiplantae": "plants", "streptophyta": "plants",
            "chlorophyta": "plants"}


def other_kingdom(lineage, ncbi_names):
    """Whether a reference's lineage is of another kingdom than the NCBI lineage's names (a genus homonym: an animal
    genus the reference holds as a plant's alone; until 2026-11-06 matched as its genus)."""
    theirs = {KINGDOMS[k] for k in (key(x) for x in lineage) if k in KINGDOMS}
    mine = {KINGDOMS[k] for k in ncbi_names if k in KINGDOMS}
    return bool(theirs and mine and not theirs & mine)


class Reference:
    """A reference's lineages by species name and by genus name: {key: Counter(lineage tuple)}."""

    def __init__(self, name, version, ranks, path=None, md5=None):
        self.name, self.version, self.ranks, self.path, self.md5 = name, version, tuple(ranks), path, md5
        self.species = collections.defaultdict(collections.Counter)
        self.genus = collections.defaultdict(collections.Counter)

    def add(self, lineage, n=1):
        """A record's lineage (len(ranks) names; species last, genus before it; "" for none)."""
        lineage = tuple(lineage)
        genus, species = lineage[-2], lineage[-1]
        if is_genus_name(genus):
            self.genus[key(genus)][lineage[:-1]] += n
        if is_species_name(species):
            self.species[key(species)][lineage] += n

    @staticmethod
    def best(counter, ncbi_names):
        """Of a name's lineages, the one sharing the most names with NCBI's lineage, then of the most records."""
        return max(counter.items(), key=lambda kv: (len({key(x) for x in kv[0]} & ncbi_names), kv[1], kv[0]))[0]

    def label(self, ncbi, taxid):
        """(lineage tuple, match) of an NCBI species: match species, synonym or genus; (None, "-") for none."""
        names = {key(x) for x in ncbi.ranks(taxid).values() if x}
        sci = ncbi.name.get(taxid, "")
        for how, candidates in (("species", [sci]), ("synonym", ncbi.synonyms.get(taxid, []))):
            for name in candidates:
                found = self.species.get(key(name))
                if found:
                    return self.best(found, names), how
        genus = ncbi.ranks(taxid).get("genus") or (sci.split() or [""])[0]
        found = self.genus.get(key(genus))
        if found:
            lineage = self.best(found, names)
            if not other_kingdom(lineage, names):
                return lineage, "genus"
        return None, "-"

    def text(self, lineage):
        """A lineage as the reference writes it (KSGP's with its k__ .. s__; a genus's without the species)."""
        if self.name == "KSGP":
            return ";".join(f"{r[0]}__{x}" for r, x in zip(self.ranks, lineage))
        return ";".join(lineage)

    def stamp(self):
        return {"version": self.version, "file": os.path.basename(self.path or ""), "md5": self.md5,
                "species": len(self.species), "genera": len(self.genus)}


def label_species(ncbi, taxids, references):
    """{taxid: {reference name: (lineage text or "-", match)}} of NCBI species."""
    out = {}
    for t in taxids:
        out[t] = {}
        for ref in references:
            lineage, how = ref.label(ncbi, t)
            out[t][ref.name] = (ref.text(lineage) if lineage else "-", how)
    return out


# --- the references ----------------------------------------------------------------------------------------------

PR2_VERSION = "5.1.1"
PR2_FILE = f"pr2_version_{PR2_VERSION}_SSU_mothur.tax.gz"
PR2_URL = f"https://github.com/pr2database/pr2database/releases/download/v{PR2_VERSION}/{PR2_FILE}"
PR2_RANKS = ("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")
KSGP_VERSION = "4.0"  # the newest KSGP its site serves (2026-11-05); a later release by --ksgp_tax
KSGP_FILE = f"KSGP_v{KSGP_VERSION}_plus2.tax"
KSGP_URL = f"https://ksgp.earlham.ac.uk/downloads/v{KSGP_VERSION}/{KSGP_FILE}"
KSGP_RANKS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
EUKARYOME_RANKS = ("superkingdom", "kingdom", "subkingdom", "phylum", "subphylum", "class", "order", "family", "genus",
                   "species")
ORGANELLE = re.compile(r"[_:](plas|mito|apic|nucl)(\b|_|$)|chloroplast|mitochondri", re.IGNORECASE)
STATE = "labels.json"


def read_pr2(path, version=PR2_VERSION):
    """PR2's mothur taxonomy (id, then its 9 ranks ;-joined) -> Reference."""
    ref = Reference("PR2", version, PR2_RANKS, path, util.md5(path))
    with util.open_text(path) as fh:
        for line in fh:
            _id, _, tax = line.rstrip("\n").partition("\t")
            names = [x for x in tax.strip().rstrip(";").split(";")]
            if len(names) != len(PR2_RANKS) or names[0] != "Eukaryota" or ORGANELLE.search(tax):
                continue
            ref.add(names)
    return ref


def ksgp_version_of(path):
    m = re.search(r"KSGP_v(\d+(?:\.\d+)*)", os.path.basename(path))
    return m.group(1) if m else "?"


def read_ksgp(path):
    """KSGP's plus2 taxonomy (id, then k__ .. s__ ;-joined) -> Reference of its eukaryote nuclear records."""
    ref = Reference("KSGP", ksgp_version_of(path), KSGP_RANKS, path, util.md5(path))
    with util.open_text(path) as fh:
        for line in fh:
            _id, _, tax = line.rstrip("\n").partition("\t")
            parts = tax.strip().rstrip(";").split(";")
            if len(parts) != len(KSGP_RANKS) or parts[0] in ("k__Bacteria", "k__Archaea") or ORGANELLE.search(tax):
                continue
            ref.add([p.split("__", 1)[1] if "__" in p else p for p in parts])
    return ref


def read_xlsx_rows(path):
    """An xlsx file's first sheet -> [tuple] (openpyxl)."""
    import openpyxl  # noqa: PLC0415
    sheet = openpyxl.load_workbook(path, read_only=True).worksheets[0]
    return [tuple("" if v is None else str(v).strip() for v in r) for r in sheet.iter_rows(values_only=True)]


def read_eukaryome(classification_rows, long_reads=(), version="2.1", path=None):
    """EUKARYOME's classification table (rows: its header, then superkingdom .. genus[, the kingdom-specific genus]) and
    its long reads' labels [(id, "k__..;s__Genus_species")] -> Reference. A long read's species takes its genus's
    lineage of the classification (by its kingdom-specific name, Drosophila(Metazoa), where the table gives one; by the
    plain name where only one lineage holds it: until 2026-11-06 a homonym's first row's), else its own k__ to g__
    (superkingdom, subkingdom and subphylum empty)."""
    ref = Reference("EUKARYOME", version, EUKARYOME_RANKS, path, util.md5(path) if path else None)
    rows = list(classification_rows)
    exact, by_key = {}, collections.defaultdict(list)
    for r in rows[1:]:
        names = list(r[:9]) + [""] * (9 - len(r[:9]))
        if names[8]:
            ref.add(names + [""])
            if names not in by_key[key(names[8])]:
                by_key[key(names[8])].append(names)
            for raw in (names[8], r[9] if len(r) > 9 else ""):
                if "(" in raw:  # a kingdom-specific name; a plain one only by by_key where it is not a homonym
                    exact.setdefault(raw.strip().lower(), names)
    from .rrna import species_name  # noqa: PLC0415
    for _id, label in long_reads:
        parts = {p.split("__", 1)[0]: p.split("__", 1)[1] for p in label.split(";") if "__" in p}
        species = species_name(parts.get("g", ""), parts.get("s", "")).replace(" ", "_")  # s__ is the epithet
        if not is_species_name(species):
            continue
        genus = parts.get("g", "")
        names = exact.get(genus.strip().lower())
        if names is None and len(by_key.get(key(genus), ())) == 1:
            names = by_key[key(genus)][0]
        if names is None:
            names = ["", parts.get("k", ""), "", parts.get("p", ""), "", parts.get("c", ""), parts.get("o", ""),
                     parts.get("f", ""), parts.get("g", "")]
        ref.add(list(names) + [species])
    return ref


def folder_of(src_dir):
    return os.path.join(src_dir, "labels")


def fetch(src_dir, download, log=print, ksgp_tax=None):
    """PR2's taxonomy and, without ksgp_tax (a KSGP release's plus2 tax file of the user's), KSGP's served one into
    WORK/sources/labels (kept when there)."""
    folder = folder_of(src_dir)
    os.makedirs(folder, exist_ok=True)
    wanted = [(PR2_FILE, PR2_URL)] + ([] if ksgp_tax else [(KSGP_FILE, KSGP_URL)])
    for name, url in wanted:
        dest = os.path.join(folder, name)
        if os.path.isfile(dest) and os.path.getsize(dest) > 0:
            log(f"labels: {dest} is there")
            continue
        log(f"labels: fetching {url}")
        download(url, dest)
    with open(os.path.join(folder, STATE), "w", newline="\n") as fh:
        json.dump({"pr2": PR2_FILE, "ksgp": ksgp_tax or KSGP_FILE}, fh)


def references(src_dir, eukaryome_folder=None, ksgp_tax=None, log=print):
    """The references found (PR2, EUKARYOME, KSGP) -> [Reference]; one missing is said and left out."""
    folder = folder_of(src_dir)
    out = []
    pr2 = os.path.join(folder, PR2_FILE)
    if os.path.isfile(pr2):
        out.append(read_pr2(pr2))
    else:
        log(f"WARNING: no {pr2} (--stage fetch): species.tsv without PR2's lineages")
    if eukaryome_folder:
        xlsx = os.path.join(eukaryome_folder, "EUKARYOME2_1_classification-2.xlsx")
        if os.path.isfile(xlsx):
            from . import sources  # noqa: PLC0415
            long_reads = [(acc, label) for acc, label, _s in sources.eukaryome_long(eukaryome_folder)] \
                if os.path.isfile(os.path.join(eukaryome_folder, "General_EUK_longread_v2.1.fasta.gz")) else []
            try:
                out.append(read_eukaryome(read_xlsx_rows(xlsx), long_reads, path=xlsx))
            except ImportError:
                log("WARNING: no openpyxl: species.tsv without EUKARYOME's lineages")
        else:
            log(f"WARNING: no {xlsx}: species.tsv without EUKARYOME's lineages")
    ksgp = ksgp_tax or os.path.join(folder, KSGP_FILE)
    if os.path.isfile(ksgp):
        out.append(read_ksgp(ksgp))
    else:
        log(f"WARNING: no {ksgp} (--stage fetch, or --ksgp_tax): species.tsv without KSGP's lineages")
    return out
