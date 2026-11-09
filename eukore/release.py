"""A package's release (since 2026-11-09): what Eukore's website (website/) shows of a package and the files it offers
to download, written by the package stage into OUT/release/, laid out as the website is, so that copying the folder
onto the site publishes it (`rsync -a OUT/release/ <site>/`):

    release/
    ├── data/stats.json        the package's counts: its species by genome source and supergroup, its marker
    │                          categories (BUSCOs, Eukore's loci, rDNA, organelle genes), its headline benchmarks
    ├── data/taxonomy.json     its species as a tree (supergroup, group, phylum ... genus, species), each node's counts
    ├── data/downloads.json    each download: file, bytes, MD5, what it holds
    └── downloads/Eukore_<version>/
        ├── Eukore_<version>.tar.gz                   the whole package
        ├── Eukore_<version>_marker_genes.tar.gz      the BUSCOs' and Eukore's loci's copies whole (exons, introns)
        ├── Eukore_<version>_marker_cds.tar.gz        their coding sequences and proteins (no introns)
        ├── Eukore_<version>_rrna_organelle.tar.gz    rDNA records and their 18S regions; the organelle genes (12S,
        │                                             16S, COX1, rbcL, ...)
        ├── Eukore_<version>_sketches.sig.gz          the genomes' sourmash sketches, one collection
        ├── Eukore_<version>_<table>.gz               each metadata table alone (species, genomes, assembly, ...)
        └── MD5SUMS

Every file is BGZF at level 9 (util.bgzf_open): gzip to every reader, decompressed on many threads by `bgzip -d -@`.
A bundle's tar holds its files under Eukore_<version>/ as in the package, so that bundles unpacked side by side make
up the package's folders (package.read_copies reads them). A rerun keeps the downloads made from the same package
(its manifest's MD5) and makes those missing.
"""

import collections
import datetime
import gzip
import json
import os
import re
import shutil
import tarfile
import time
import zipfile

from . import metadata, package, util
from . import __version__

DATA, DOWNLOADS = "data", "downloads"
STATE = ".release.json"
SOURCES = ("NCBI RefSeq", "NCBI GenBank", "TOGA2", "ENA TSA", "EukProt", "rRNA records only")
SOURCE_TEXT = {"NCBI RefSeq": "NCBI's RefSeq assembly", "NCBI GenBank": "NCBI's GenBank assembly",
               "TOGA2": "TOGA2's assembly (vertebrates with no NCBI genome)",
               "ENA TSA": "an INSDC transcriptome (no genome)",
               "EukProt": "EukProt's genome or transcriptome records (no assembly at NCBI)",
               "rRNA records only": "published rRNA reads alone (EUKARYOME, ROD)"}
TREE_RANKS = ("supergroup", "group", "phylum", "class", "order", "family", "genus", "species")
SHARED = ("species.tsv", "genomes.tsv")  # in every bundle: what the accessions and taxids are
SKETCHES = "metadata/sketches.zip"
TABLES = (("species.tsv", "Species", "one row per species: supergroup, group, NCBI's lineage and six references'"),
          ("genomes.tsv", "Genomes", "each species' genome: accession, quality flags, rRNA copies, genetic code, "
                                     "duplication"),
          ("crossrefs.tsv", "Cross-references", "each species' Open Tree, GBIF, WoRMS, IRMNG, UNITE SH and BOLD ids"),
          ("copy_numbers.tsv.gz", "Copy numbers", "each genome's copies of every locus, rRNA gene and unit, "
                                                  "organelle molecule and gene"),
          ("metadata/assembly.tsv.gz", "Assemblies", "each genome's size, contiguity, gaps, GC, masking, telomeres, "
                                                     "organelle molecules, where it came from"),
          ("metadata/ncbi.tsv.gz", "NCBI's assembly summary", "NCBI's own row of every genome"),
          ("metadata/tetranucleotides.tsv.gz", "Tetranucleotides", "each genome's 136 canonical tetranucleotides"),
          ("metadata/codons.tsv.gz", "Codons", "the genetic-code check's evidence per genome"))


def _rows(folder, name, columns):
    """A package table's rows, streamed, as tuples of the columns asked ('-' where a column is not there): the copy
    tables of a full build have millions of rows, never held as dicts."""
    path = os.path.join(folder, name)
    if not os.path.isfile(path):
        return
    with util.open_text(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        at = [head.index(c) if c in head else None for c in columns]
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            yield tuple(f[i] if i is not None and i < len(f) else "-" for i in at)


def source_of(genome):
    """A genome row's source as the website groups them (SOURCES; TOGA2 whoever assembled it)."""
    s = metadata.genome_source(genome)
    return "TOGA2" if s.startswith("TOGA2") else s if s in SOURCES else "other"


def _path(species):
    """A species row's place in the tree: [(rank, name)] from its supergroup down to itself, a name its parent's own
    (a group that is also its phylum) once; an empty rank the name above it with _X (PR2's way, as
    ncbi_lineage_filled's own), empty ranks in a row one level."""
    lineage = (species.get("ncbi_lineage_filled") or species.get("ncbi_lineage") or "").split(";")
    lineage += [""] * (8 - len(lineage))
    names = [species.get("supergroup") or "incertae sedis", species.get("group") or "other eukaryotes"] + \
        lineage[2:7] + [species.get("name") or str(species.get("taxid"))]
    out, above = [], "-"
    for rank, name in zip(TREE_RANKS, names):
        name = name.strip() if name and name != "-" else ""
        if name:
            above = re.sub(r"(_X)+$", "", name) or name
        else:
            name = f"{above}_X"
        if out and name.lower() == out[-1][1].lower() and rank != "species":
            continue
        out.append((rank, name))
    return out


def tree(species_rows):
    """The species as a tree -> [name, rank, species, with a genome, children] nested, children by species, most first;
    a species is [name, rank, 1, 0 or 1, its genome's accession or '-']. rank: an index of TREE_RANKS (-1: the root)."""
    root = {"name": "Eukaryota", "rank": -1, "species": 0, "genomes": 0, "children": {}}
    for r in species_rows:
        has_genome = int(r.get("representative", "-") != "-")
        node = root
        node["species"] += 1
        node["genomes"] += has_genome
        for rank, name in _path(r):
            child = node["children"].get(name)
            if child is None:
                child = node["children"][name] = {"name": name, "rank": TREE_RANKS.index(rank), "species": 0,
                                                  "genomes": 0, "children": {}, "accession": "-"}
            child["species"] += 1
            child["genomes"] += has_genome
            if rank == "species":
                child["accession"] = r.get("representative", "-")
            node = child

    def packed(n):
        if n["rank"] == TREE_RANKS.index("species"):
            return [n["name"], n["rank"], n["species"], n["genomes"], n["accession"]]
        kids = sorted(n["children"].values(), key=lambda c: (-c["species"], c["name"].lower()))
        return [n["name"], n["rank"], n["species"], n["genomes"], [packed(c) for c in kids]]
    return packed(root)


def _rank_counts(node, out=None):
    out = collections.Counter() if out is None else out
    if node[1] >= 0:
        out[TREE_RANKS[node[1]]] += 1
    if isinstance(node[4], list):
        for c in node[4]:
            _rank_counts(c, out)
    return out


def _histogram(values, edges=(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0001)):
    out = []
    for lo, hi in zip(edges, edges[1:]):
        out.append(dict(name=f"{lo:.0%}-{min(hi, 1):.0%}", loci=sum(1 for v in values if lo <= v < hi)))
    return out


def _loci_category(folder, locus_set, species_of, title, what):
    """A locus set's (busco/, markers/) copies: loci, copies, species, confidence; per locus its copies and species."""
    loci = {}
    for locus, name, call_rate in _rows(folder, f"{locus_set}/loci.tsv", ("locus", "name", "call_rate")):
        loci[locus] = dict(name=name if name not in ("", "-") else locus, call_rate=_float(call_rate))
    if not loci:
        return None
    copies, confidence = collections.Counter(), collections.Counter()
    species, per_locus = set(), collections.defaultdict(set)
    for acc, locus, status, conf in _rows(folder, f"{locus_set}/copies.tsv",
                                          ("accession", "locus", "status", "confidence")):
        if status == "Missing":
            continue
        copies[locus] += 1
        confidence[conf] += 1
        sp = species_of.get(acc)
        if sp is not None:
            species.add(sp)
            per_locus[locus].add(sp)
    out = dict(key=locus_set, title=title, what=what, loci=len(loci), sequences=sum(copies.values()),
               species=len(species), confidence=_confidence(confidence))
    if locus_set == "busco":  # 129 loci: their call rates' spread, not one bar each
        out["items_title"] = "BUSCOs by the share of species they are called in"
        out["items"] = [dict(name=h["name"], sequences=h["loci"]) for h in
                        _histogram([v["call_rate"] for v in loci.values() if v["call_rate"] is not None])]
        out["items_unit"] = "BUSCOs"
    else:
        out["items_title"] = "copies and species per locus"
        out["items"] = [dict(name=v["name"], locus=k, sequences=copies.get(k, 0), species=len(per_locus.get(k, ())))
                        for k, v in loci.items()]
        out["items_unit"] = "copies"
    return out


def _float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _confidence(counter):
    return {k: counter.get(k, 0) for k in ("high", "medium", "low")}


RRNA_REGIONS = ("18S", "ITS1", "5.8S", "ITS2", "28S", "5S")
ORGANELLE_NAMES = {("mitochondrion", "SSU"): "SSU rRNA (12S)", ("mitochondrion", "LSU"): "LSU rRNA (16S)",
                   ("plastid", "SSU"): "SSU rRNA (16S)", ("plastid", "LSU"): "LSU rRNA (23S)",
                   ("mitochondrion", "5S"): "5S rRNA", ("plastid", "5S"): "5S rRNA"}
ORGANELLE_ORDER = ("SSU", "LSU", "5S", "COX1", "COX3", "CYTB", "nad1", "nad5", "atp6", "rbcL", "matK", "tufA", "psbA",
                   "ycf1", "trnH-psbA")


def _rrna_category(folder):
    records, species, confidence = 0, set(), collections.Counter()
    regions, region_species = collections.Counter(), collections.defaultdict(set)
    sources, units = collections.Counter(), 0
    for kind, sp, source, regs, conf in _rows(folder, "rrna/records.tsv",
                                              ("kind", "species_taxid", "source", "regions", "confidence")):
        records += 1
        species.add(sp)
        confidence[conf] += 1
        sources[source] += 1
        units += kind == "rDNA"
        parts = {p.split(":")[0] for p in regs.split(",") if p and p != "-"}
        for part in parts or ({"18S"} if kind == "18S" else ()):  # an 18S span is its own region (regions -)
            regions[part] += 1
            region_species[part].add(sp)
    if not records:
        return None
    items = [dict(name=r, sequences=regions[r], species=len(region_species[r])) for r in RRNA_REGIONS if regions[r]]
    items.append(dict(name="whole units (18S to 28S)", sequences=units, species=None))
    return dict(key="rrna", title="Nuclear rDNA", what="18S spans and whole operons (18S, ITS1, 5.8S, ITS2, 28S, "
                "5S) of the genomes, EUKARYOME's and ROD's long reads, with their 18S V1-V9",
                loci=sum(1 for r in RRNA_REGIONS if regions[r]), sequences=records, species=len(species),
                confidence=_confidence(confidence), items_title="records holding each region", items=items,
                items_unit="records", sources=dict(sources.most_common()))


def _organelle_categories(folder):
    out = []
    for organelle, title in (("mitochondrion", "Mitochondrial genes"), ("plastid", "Plastid genes")):
        records, species, confidence = 0, set(), collections.Counter()
        genes, gene_species = collections.Counter(), collections.defaultdict(set)
        for gene, org, sp, conf in _rows(folder, "organelle/genes.tsv",
                                         ("gene", "organelle", "species_taxid", "confidence")):
            if org != organelle:
                continue
            records += 1
            species.add(sp)
            confidence[conf] += 1
            genes[gene] += 1
            gene_species[gene].add(sp)
        if not records:
            continue
        order = {g: i for i, g in enumerate(ORGANELLE_ORDER)}
        items = [dict(name=ORGANELLE_NAMES.get((organelle, g), g), sequences=n, species=len(gene_species[g]))
                 for g, n in sorted(genes.items(), key=lambda kv: (order.get(kv[0], 99), kv[0]))]
        out.append(dict(key=organelle, title=title,
                        what=("rRNA (12S, 16S, 5S) and COX1, COX3, CYTB, nad1, nad5, atp6, with their proteins by "
                              "the organelle's code" if organelle == "mitochondrion" else
                              "rRNA (16S, 23S, 5S) and rbcL, matK, tufA, psbA, ycf1, the trnH-psbA spacer"),
                        loci=len(genes), sequences=records, species=len(species), confidence=_confidence(confidence),
                        items_title="records and species per gene", items=items, items_unit="records"))
    return out


def _benchmarks(folder):
    try:
        with open(os.path.join(folder, "benchmarks", "summary.json")) as fh:
            summary = json.load(fh)
    except (OSError, ValueError):
        return None
    out = {}
    if "completeness" in summary:
        out["completeness"] = {k: summary["completeness"].get(k) for k in ("median", "at_least_90", "by_supergroup")}
    if "resolution" in summary:
        out["resolution"] = summary["resolution"]
    if "rrna" in summary:
        out["rrna"] = summary["rrna"]
    if "primers" in summary:
        out["primers"] = summary["primers"]
    return out


def stats(folder):
    """A package's counts for the website -> (stats, taxonomy tree): its species per genome source and supergroup
    (and group), its marker categories, its headline benchmarks; the tree (tree)."""
    manifest = package.read(folder) or {}
    species_rows = list(package.read_table(folder, "species.tsv"))
    genomes = {}
    for acc, fetched, source, sp, rep in _rows(folder, "genomes.tsv", ("accession", "fetched", "source",
                                                                       "species_taxid", "representative")):
        if rep == "t":
            genomes[acc] = dict(accession=acc, fetched=fetched, source=source, species_taxid=sp)
    species_of = {acc: g["species_taxid"] for acc, g in genomes.items()}
    by_group = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    for r in species_rows:
        rep = r.get("representative", "-")
        source = "rRNA records only" if rep == "-" else source_of(genomes.get(rep, dict(source=r.get("source", ""),
                                                                                         accession=rep)))
        by_group[r.get("supergroup") or "incertae sedis"][r.get("group") or "other eukaryotes"][source] += 1
    columns = [s for s in SOURCES + ("other",) if any(c[s] for g in by_group.values() for c in g.values())]
    supergroups = []
    for sg, groups in by_group.items():
        total = collections.Counter()
        for c in groups.values():
            total.update(c)
        supergroups.append(dict(name=sg, counts={s: total.get(s, 0) for s in columns}, species=sum(total.values()),
                                groups=sorted(({"name": g, "counts": {s: c.get(s, 0) for s in columns},
                                                "species": sum(c.values())} for g, c in groups.items()),
                                              key=lambda x: (-x["species"], x["name"]))))
    supergroups.sort(key=lambda x: (-x["species"], x["name"]))
    taxonomy = tree(species_rows)
    ranks = _rank_counts(taxonomy)
    categories = [c for c in (
        _loci_category(folder, "busco", species_of, "BUSCO genes",
                       "the near-universal single-copy orthologs of BUSCO's eukaryota_odb12, each copy whole with its "
                       "introns, its coding sequence and protein"),
        _loci_category(folder, "markers", species_of, "Eukore's marker loci",
                       "the classic phylogenetic markers BUSCO does not hold (RPB1, RPB2, EF-2, TEF1, β-tubulin, actin, "
                       "HSP90, GAPDH, histone H3)"),
        _rrna_category(folder)) if c] + _organelle_categories(folder)
    counts = manifest.get("counts", {})
    out = dict(package=os.path.basename(os.path.normpath(folder)), version=manifest.get("version"),
               eukore=manifest.get("eukore"), format=manifest.get("format"), built=manifest.get("built"),
               released=datetime.datetime.now().isoformat(timespec="seconds"), release_by=__version__,
               partial=manifest.get("inputs", {}).get("partial"),
               totals=dict(species=len(species_rows),
                           species_with_genome=sum(1 for r in species_rows if r.get("representative", "-") != "-"),
                           species_rrna_only=sum(1 for r in species_rows if r.get("representative", "-") == "-"),
                           genomes=len(genomes), assemblies_listed=counts.get("genomes"),
                           supergroups=ranks.get("supergroup", 0), groups=ranks.get("group", 0),
                           phyla=ranks.get("phylum", 0), classes=ranks.get("class", 0), orders=ranks.get("order", 0),
                           families=ranks.get("family", 0), genera=ranks.get("genus", 0),
                           marker_categories=len(categories),
                           marker_loci=sum(c["loci"] for c in categories),
                           sequences=sum(c["sequences"] for c in categories)),
               sources=dict(columns=columns, text={s: SOURCE_TEXT.get(s, s) for s in columns},
                            supergroups=supergroups),
               markers=categories, benchmarks=_benchmarks(folder), ranks=list(TREE_RANKS))
    return out, dict(ranks=list(TREE_RANKS), tree=taxonomy)


# --- downloads -----------------------------------------------------------------------------------------------------

def _in_sets(rel, subfolders, tables):
    parts = rel.split("/")
    return parts[0] in package.LOCUS_SETS and len(parts) >= 2 and (
        (len(parts) > 2 and parts[1] in subfolders) or (len(parts) == 2 and parts[1] in tables))


def _marker_genes(rel):
    return rel in SHARED or _in_sets(rel, ("copies",), ("loci.tsv", "copies.tsv", "features.tsv", "introns.tsv",
                                                       "intron_sites.tsv"))


def _marker_cds(rel):
    return rel in SHARED or _in_sets(rel, ("cds", "proteins"), ("loci.tsv", "copies.tsv"))


def _rrna_organelle(rel):
    return rel in SHARED or rel.startswith("rrna/") or (rel.startswith("organelle/")
                                                         and rel != "organelle/molecules.fna.gz")


BUNDLES = (("", "Complete package", "Eukore, whole",
            "every file of the package: the marker genes, rDNA, organelle genes, metadata, taxonomy and benchmarks, "
            "with the manifest of their sizes and MD5s", lambda rel: True),
           ("_marker_genes", "Marker genes", "Marker genes with introns",
            "the BUSCO genes and Eukore's own marker loci, each copy whole from first to last codon (exons and "
            "introns), one FASTA per locus; with their tables (copies, loci, exons and introns) and species.tsv, "
            "genomes.tsv", _marker_genes),
           ("_marker_cds", "Marker genes", "Marker genes without introns",
            "the same copies' coding sequences (exons joined) and proteins, one FASTA per locus; with copies.tsv, "
            "loci.tsv, species.tsv, genomes.tsv", _marker_cds),
           ("_rrna_organelle", "rRNA and organelle markers", "rRNA operons and organelle markers",
            "nuclear rDNA records (18S spans and whole operons: 18S, ITS1, 5.8S, ITS2, 28S, 5S), their 18S V1-V9 and "
            "in-silico PCR; the mitochondrial and plastid genes (12S, 16S, 5S, COX1, CYTB, rbcL, matK, ...) with "
            "their coding sequences, proteins and PCR; species.tsv, genomes.tsv (the organelle molecules whole are in "
            "the complete package)", _rrna_organelle))


def _plain(info):
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    info.mode = 0o644 if info.isfile() else 0o755
    return info


def _tar(folder, rels, dest, top):
    """rels of folder as a BGZF tar (dest), each under top/."""
    part = util.part_of(dest)
    with util.bgzf_open(part, "wb") as raw, tarfile.open(fileobj=raw, mode="w|", format=tarfile.PAX_FORMAT) as tar:
        for rel in rels:
            tar.add(os.path.join(folder, rel), arcname=f"{top}/{rel}", recursive=False, filter=_plain)
    os.replace(part, dest)


def _sketches(src, dest):
    """The package's sketches.zip (sourmash signatures, one file a genome) as one signature file (a JSON list),
    BGZF -> the signatures. Streamed: a full build's sketches never held at once."""
    part, n = util.part_of(dest), 0
    with zipfile.ZipFile(src) as z, util.bgzf_open(part) as out:
        out.write("[")
        for name in sorted(z.namelist()):
            if not name.endswith((".sig", ".sig.gz", ".json", ".json.gz")):
                continue
            data = z.read(name)
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            got = json.loads(data)
            for sig in got if isinstance(got, list) else [got]:
                out.write((",\n" if n else "") + json.dumps(sig, separators=(",", ":")))
                n += 1
        out.write("]\n")
    os.replace(part, dest)
    return n


def _table(src, dest):
    """A package table as dest, BGZF at level 9: a BGZF table copied as it is, another compressed again."""
    if src.endswith(".gz") and util.is_bgzf(src):
        part = util.part_of(dest)
        shutil.copyfile(src, part)
        os.replace(part, dest)
    else:
        util.bgzf_copy(src, dest)


def package_files(folder):
    """Every file of a package (relative, "/"), the manifest first."""
    rels = []
    for dirpath, _dirs, names in os.walk(folder):
        for name in names:
            rels.append(os.path.relpath(os.path.join(dirpath, name), folder).replace(os.sep, "/"))
    return sorted(rels, key=lambda r: (r != package.MANIFEST, r))


def downloads(folder, dest, log=print):
    """The package's downloads into dest (a folder of its own) -> [entry: key, group, title, what, file, bytes, md5,
    files]. Kept from an earlier run where the package is the same (STATE: its manifest's MD5), else made again."""
    name = os.path.basename(os.path.normpath(folder))
    stamp = util.md5(os.path.join(folder, package.MANIFEST))
    try:
        with open(os.path.join(dest, STATE)) as fh:
            same = json.load(fh).get("manifest_md5") == stamp
    except (OSError, ValueError):
        same = False
    if not same and os.path.isdir(dest):
        shutil.rmtree(dest)
    os.makedirs(dest, exist_ok=True)
    with open(os.path.join(dest, STATE), "w", newline="\n") as fh:
        json.dump({"manifest_md5": stamp, "package": name}, fh)
    rels = package_files(folder)
    plan = []
    for suffix, group, title, what, chosen in BUNDLES:
        members = [r for r in rels if chosen(r)]
        if suffix and all(r in SHARED for r in members):  # a bundle of the shared tables alone: none
            continue
        plan.append(dict(key=suffix.strip("_") or "complete", group=group, title=title, what=what,
                         file=f"{name}{suffix}.tar.gz", files=len(members),
                         make=lambda path, m=members: _tar(folder, m, path, name)))
    if os.path.isfile(os.path.join(folder, SKETCHES)):
        plan.append(dict(key="sketches", group="Sourmash sketches", title="Sourmash sketches of the genomes",
                         what="each genome's sourmash signature (taken before the build deleted the genome), one "
                              "collection: sourmash search, gather and compare read it as it is",
                         file=f"{name}_sketches.sig.gz", files=None,
                         make=lambda path: _sketches(os.path.join(folder, SKETCHES), path)))
    for rel, title, what in TABLES:
        if os.path.isfile(os.path.join(folder, rel)):
            base = os.path.basename(rel)
            base = base[:-3] if base.endswith(".gz") else base
            plan.append(dict(key=base.split(".")[0], group="Metadata tables", title=title, what=what,
                             file=f"{name}_{base}.gz", files=None, source=rel,
                             make=lambda path, r=rel: _table(os.path.join(folder, r), path)))
    out = []
    for item in plan:
        path = os.path.join(dest, item["file"])
        if not os.path.isfile(path):
            started = time.time()
            got = item["make"](path)
            if item["key"] == "sketches":
                item["files"] = got
            log(f"release: {item['file']} ({util.size_text(os.path.getsize(path))}, "
                f"{util.duration(time.time() - started)})")
        elif item["key"] == "sketches":
            with gzip.open(path, "rt") as fh:
                item["files"] = len(json.load(fh))
        entry = {k: v for k, v in item.items() if k != "make"}
        entry.update(bytes=os.path.getsize(path), md5=util.md5(path))
        out.append(entry)
    with open(util.part_of(os.path.join(dest, "MD5SUMS")), "w", newline="\n") as fh:
        fh.writelines(f"{e['md5']}  {e['file']}\n" for e in out)
    os.replace(util.part_of(os.path.join(dest, "MD5SUMS")), os.path.join(dest, "MD5SUMS"))
    return out


def _json(path, value, compact=False):
    part = util.part_of(path)
    with open(part, "w", newline="\n", encoding="utf-8") as fh:
        if compact:
            json.dump(value, fh, separators=(",", ":"), ensure_ascii=False)
        else:
            json.dump(value, fh, indent=1, ensure_ascii=False)
        fh.write("\n")
    os.replace(part, path)


def write(folder, out, log=print, files=True):
    """The package's release into out (data/, downloads/<package>/: the module's doc) -> out. files: the downloads
    made too (False: the website's data alone, the downloads' table kept)."""
    if package.read(folder) is None:
        raise SystemExit(f"release: {folder} is no finished Eukore package (no {package.MANIFEST})")
    name = os.path.basename(os.path.normpath(folder))
    started = time.time()
    numbers, taxonomy = stats(folder)
    data = os.path.join(out, DATA)
    os.makedirs(data, exist_ok=True)
    _json(os.path.join(data, "stats.json"), numbers)
    _json(os.path.join(data, "taxonomy.json"), taxonomy, compact=True)
    if files:
        entries = downloads(folder, os.path.join(out, DOWNLOADS, name), log)
        _json(os.path.join(data, "downloads.json"),
              dict(package=name, version=numbers["version"], built=numbers["built"],
                   folder=f"{DOWNLOADS}/{name}/", compression="BGZF (gzip), level 9", files=entries))
    log(f"release: {out} ({numbers['totals']['species']:,} species, {len(numbers['markers'])} marker categories"
        + (f", {len(entries)} downloads, {util.size_text(sum(e['bytes'] for e in entries))}" if files else "")
        + f"; {util.duration(time.time() - started)})")
    return out
