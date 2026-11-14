"""Eukore's rRNA records linked to PR2's and KSGP's (since 2026-11-14, format 7): taxonomy/links.tsv.

Eukore's published 18S are EUKARYOME's alone (its long reads and its SSU set); PR2 and KSGP give no record or taxon.
Most of their eukaryote records are EUKARYOME's under the same id or accession (KSGP 5.2: all 342,515 of its
EUKARYOME ids, 127,675 of its 134,167 other nuclear 18S; 2026-11-14), so a user who held them (protal, whose KSGP now
holds prokaryote 16S alone) finds each Eukore record's PR2 and KSGP records and their labels here.

A record's accession is its EUKARYOME id or INSDC accession (EUK_<id>: the id), ROD's operon its NCBI contig's; PR2's
ids (AB353770.1.1740_U) and KSGP's (an EUKARYOME id, an accession, `_dup<n>` a further record of it, `<x>~<id>` an
addition named after another record) are matched by their accession without its version, span or suffix. One row a
record and counterpart; a record with none has no row.
"""

import collections
import os
import re

from . import util

LINK_COLUMNS = ("record", "accession", "reference", "id", "lineage")
PROKARYOTE = ("k__Bacteria", "k__Archaea")


def accession_key(rid):
    """An id's accession: KSGP's addition `<x>~<id>` its id, `_dup<n>`, PR2's `.start.end_<code>` and the version left
    out (AB353770.1.1740_U, AB353770_dup1, AB353770.1 -> AB353770)."""
    rid = rid.split("~")[-1]
    rid = re.sub(r"_dup\d+$", "", rid)
    rid = re.sub(r"\.\d+\.\d+_[A-Za-z]+$", "", rid)
    return rid.split(".")[0]


def kind_of(label):
    """A KSGP label's compartment: nuclear, plastid, mitochondrion, nucleomorph (its k__ rank's :plas, :mito, :nucl),
    or None for a prokaryote."""
    k = label.split(";", 1)[0]
    if k.startswith(PROKARYOTE):
        return None
    for tag, kind in ((":plas", "plastid"), (":mito", "mitochondrion"), (":nucl", "nucleomorph")):
        if tag in k:
            return kind
    return "nuclear"


def read_index(path, reference):
    """A PR2 mothur taxonomy or a KSGP plus2 tax file -> {accession: [(id, lineage)]} of its eukaryote records."""
    out = collections.defaultdict(list)
    with util.open_text(path) as fh:
        for line in fh:
            rid, _, lineage = line.rstrip("\n").partition("\t")
            lineage = lineage.strip().replace("; ", ";")
            if not rid or not lineage:
                continue
            if reference == "KSGP" and kind_of(lineage) is None:
                continue
            if reference == "PR2" and not lineage.startswith("Eukaryota"):
                continue
            out[accession_key(rid)].append((rid, lineage))
    return out


def record_accession(name, row):
    """A record's accession for the links: EUKARYOME's id or accession, ROD's contig, else None."""
    if row.get("source") == "EUKARYOME" and name.startswith("EUK_"):
        return accession_key(name[len("EUK_"):])
    acc = row.get("contig") or row.get("seqid")
    return accession_key(acc) if acc not in (None, "-", "") else None


def rows(records, indexes):
    """records [(name, sequence, row)], indexes {reference: read_index's} -> LINK_COLUMNS rows, and the counts
    {reference: records linked}."""
    out, counts = [], collections.Counter()
    for name, _seq, row in records:
        acc = record_accession(name, row)
        if not acc:
            continue
        for reference, index in indexes.items():
            hits = index.get(acc, ())
            if hits:
                counts[reference] += 1
            for rid, lineage in hits:
                out.append(dict(record=name, accession=acc, reference=reference, id=rid, lineage=lineage))
    return out, dict(counts)


def references(src_dir, ksgp_tax=None, log=print):
    """The link files of a build -> ({reference: (path, version)}): PR2's mothur taxonomy (labels.fetch's) and a KSGP
    tax file (--ksgp_tax, else KSGP's served one); one missing is said and left out."""
    from . import labels  # noqa: PLC0415
    folder = labels.folder_of(src_dir)
    out = {}
    pr2 = os.path.join(folder, labels.PR2_FILE)
    if os.path.isfile(pr2):
        out["PR2"] = (pr2, labels.PR2_VERSION)
    else:
        log(f"WARNING: no {pr2}: taxonomy/links.tsv without PR2's records")
    ksgp = ksgp_tax or os.path.join(folder, labels.KSGP_FILE)
    if os.path.isfile(ksgp):
        out["KSGP"] = (ksgp, labels.ksgp_version_of(ksgp))
    else:
        log(f"WARNING: no {ksgp}: taxonomy/links.tsv without KSGP's records")
    return out


def table(records, src_dir, ksgp_tax=None, log=print):
    """taxonomy/links.tsv's rows and its stamp ({reference: version, file, md5, records linked})."""
    refs = references(src_dir, ksgp_tax, log)
    indexes = {r: read_index(path, r) for r, (path, _v) in refs.items()}
    out, counts = rows(records, indexes)
    stamp = {r: {"version": v, "file": os.path.basename(path), "md5": util.md5(path),
                 "records_linked": counts.get(r, 0)} for r, (path, v) in refs.items()}
    if log:
        log("taxonomy/links.tsv: " + ", ".join(f"{r} {n:,} records" for r, n in sorted(counts.items())))
    return out, stamp
