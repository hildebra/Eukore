"""NCBI's own quality reports on its assemblies, as flags Eukore ranks a species' genomes by (since 2026-11-06).

- FCS-GX (NCBI's Foreign Contamination Screen, run on every GenBank and RefSeq assembly; fcs_summary_*.txt.gz): an
  assembly whose contaminant share is MOSTLY or more is another organism's genome under the species' name
  (FCS_FOREIGN: never taken; Azolla's genomes are its cyanobacterium's, 2026-11-06); one with CONTAMINATED of it, or
  EUK_BP of another eukaryote's sequence (what can give a foreign gene or 18S copy), is taken only without better
  (FCS_CONTAMINATED). The accession looked up is the one Eukore fetches: the RefSeq copy where there is one (RefSeq's
  copies are screened, and cleaned).
- ANI (NCBI's species check of Saccharomycotina genomes against their species' type-strain genomes;
  ANI_report_saccharomycotina.txt): "Failed" is another species' genome (ANI_FAILED: never taken), "Inconclusive"
  mostly one below 95% to its species' type (ANI_INCONCLUSIVE: taken only without better).
genomes.rank() then prefers, of otherwise equal assemblies, one with no contaminant at all (row["fcs_bp"]).
"""

import gzip
import os

REPORTS = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS"
FILES = ("fcs_summary_genbank.txt.gz", "fcs_summary_refseq.txt.gz", "fungi/ANI_report_saccharomycotina.txt")
MOSTLY = 0.5           # contaminant share of an assembly that is another organism's genome
CONTAMINATED = 0.01    # contaminant share that demotes it
EUK_BP = 100_000       # bases of another eukaryote's sequence that demote it
EUK_COLUMNS = ("primates", "anml_minus_primates", "fung", "plnt", "prst")
FCS_FOREIGN = "FCS-GX: mostly foreign"
FCS_CONTAMINATED = "FCS-GX: contaminated"
ANI_FAILED = "ANI: failed"
ANI_INCONCLUSIVE = "ANI: inconclusive"
EXCLUDE = {FCS_FOREIGN, ANI_FAILED}


def folder_of(src_dir):
    return os.path.join(src_dir, "ncbi_quality")


def fetch(src_dir, download, log=print):
    """NCBI's reports (FILES, ~70 MB) into src_dir/ncbi_quality, kept when there."""
    folder = folder_of(src_dir)
    os.makedirs(folder, exist_ok=True)
    for name in FILES:
        dest = os.path.join(folder, os.path.basename(name))
        if not os.path.isfile(dest):
            log(f"NCBI: fetching {os.path.basename(name)}")
            download(f"{REPORTS}/{name}", dest)


def paths(src_dir):
    """The reports fetched -> [path]."""
    return [p for p in (os.path.join(folder_of(src_dir), os.path.basename(n)) for n in FILES) if os.path.isfile(p)]


def number(f, k):
    """A report's cell as a number, 0 where empty or not one ("na", "-": until 2026-11-06 one stopped the plan)."""
    try:
        return float(f.get(k) or 0)
    except ValueError:
        return 0.0


def read_fcs(src_dir, wanted=None):
    """FCS-GX's summaries -> {accession: dict(fraction, bp (contaminant bases), euk_bp)}, of `wanted` accessions."""
    out = {}
    for name in FILES[:2]:
        path = os.path.join(folder_of(src_dir), name)
        if not os.path.isfile(path):
            continue
        with gzip.open(path, "rt") as fh:
            head = fh.readline().lstrip("#").rstrip("\n").split("\t")
            for line in fh:
                acc = line.split("\t", 1)[0]
                if wanted is not None and acc not in wanted:
                    continue
                f = dict(zip(head, line.rstrip("\n").split("\t")))
                out[acc] = dict(fraction=number(f, "fraction_contaminated"),
                                bp=int(number(f, "total_contaminant_length")),
                                euk_bp=int(sum(number(f, k) for k in EUK_COLUMNS)))
    return out


def read_ani(src_dir):
    """The Saccharomycotina ANI report -> {GenBank accession: taxonomy-check status}."""
    path = os.path.join(folder_of(src_dir), os.path.basename(FILES[2]))
    out = {}
    if os.path.isfile(path):
        with open(path) as fh:
            head = None
            for line in fh:
                if line.startswith("#"):
                    head = line.lstrip("#").strip().split("\t")
                    continue
                f = dict(zip(head or (), line.rstrip("\n").split("\t")))
                if f.get("genbank-accession"):
                    out[f["genbank-accession"]] = f.get("taxonomy-check-status", "")
    return out


def flags_of(fcs, ani_status):
    """An assembly's flags by its FCS-GX summary (or None) and ANI status (or None)."""
    flags = []
    if fcs:
        if fcs["fraction"] >= MOSTLY:
            flags.append(FCS_FOREIGN)
        elif fcs["fraction"] >= CONTAMINATED or fcs["euk_bp"] >= EUK_BP:
            flags.append(FCS_CONTAMINATED)
    if ani_status == "Failed":
        flags.append(ANI_FAILED)
    elif ani_status == "Inconclusive":
        flags.append(ANI_INCONCLUSIVE)
    return flags


def annotate(assemblies, refseq_rows, src_dir):
    """Each GenBank assembly row (sources.read_assemblies) given its flags by NCBI's reports (added to row["flags"])
    and row["fcs_bp"] (its contaminant bases, None where not screened), by the accession Eukore fetches (its RefSeq
    copy where it has one) -> {flag: assemblies flagged}."""
    copy_of = {r["paired"]: r["accession"] for r in refseq_rows if r.get("paired", "na").startswith("GCA_")}
    wanted = {r["accession"] for r in assemblies} | set(copy_of.values())
    fcs, ani = read_fcs(src_dir, wanted), read_ani(src_dir)
    counts = {}
    for row in assemblies:
        screen = fcs.get(copy_of.get(row["accession"], ""), fcs.get(row["accession"]))
        row["fcs_bp"] = screen["bp"] if screen else None
        new = flags_of(screen, ani.get(row["accession"]))
        if new:
            row["flags"] = sorted(set(row.get("flags") or ()) | set(new))
            for f in new:
                counts[f] = counts.get(f, 0) + 1
    return counts
