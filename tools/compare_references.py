#!/usr/bin/env python3
"""Other references for the read benchmark (eukore_build.py --read_compare NAME=FASTA,TABLE; eukore/reads.py), each as
a FASTA and its record<TAB>NCBI taxid table (since 2026-11-06):

- metaphlan: MetaPhlAn 4's eukaryote markers (Blanco-Miguez et al. 2023), from its database's archive
  (mpa_<version>.tar, or its marker FASTA, .fna.bz2) and its marker_info.txt.bz2: the markers whose clade is t__EUK<taxid>
  (an NCBI taxid), from https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/.
- pr2: PR2's 18S (Guillou et al. 2013; pr2_version_<V>_SSU_taxo_long.fasta.gz, its 9 ranks in each title): a record's
  species where NCBI names it (its name or a synonym, as species.tsv's labels find them), else its genus; organelle
  records (plastid, mitochondrion, nucleomorph) and those naming neither are left out.

    python3 tools/compare_references.py metaphlan --archive mpa_vJan26_CHOCOPhlAnSGB_202605.tar \\
        --marker_info mpa_vJan26_CHOCOPhlAnSGB_202605_marker_info.txt.bz2 --out DIR
    python3 tools/compare_references.py pr2 --fasta pr2_version_5.1.1_SSU_taxo_long.fasta.gz \\
        --taxdump WORK/sources/taxdump.tar.gz --out DIR

DIR gets <name>.fna.gz and <name>_taxids.tsv; the build then takes --read_compare NAME=DIR/<name>.fna.gz,
DIR/<name>_taxids.tsv.
"""

import argparse
import bz2
import gzip
import io
import os
import re
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import labels, ncbi as ncbi_mod, rrna  # noqa: E402

EUK_CLADE = re.compile(r"'clade': 't__EUK(\d+)'")


def text_of(path):
    """A text file, plain, gzipped or bzip2'd (by its name)."""
    if path.endswith(".bz2"):
        return bz2.open(path, "rt")
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def fasta_records(fh):
    """(title, sequence) of a FASTA stream; title without its ">"."""
    title, parts = None, []
    for line in fh:
        if line.startswith(">"):
            if title is not None:
                yield title, "".join(parts)
            title, parts = line[1:].rstrip("\n"), []
        else:
            parts.append(line.strip())
    if title is not None:
        yield title, "".join(parts)


def euk_markers(marker_info):
    """MetaPhlAn's marker_info lines (marker<TAB>its dict) -> {marker: NCBI taxid} of the eukaryote markers."""
    out = {}
    for line in marker_info:
        name, _tab, info = line.partition("\t")
        m = EUK_CLADE.search(info)
        if m:
            out[name] = int(m.group(1))
    return out


def marker_fasta(archive):
    """The markers' FASTA stream: the archive's .fna.bz2 member, or the .fna(.bz2|.gz) file itself."""
    if not tarfile.is_tarfile(archive):
        return text_of(archive)
    tar = tarfile.open(archive)
    member = next((m for m in tar.getmembers() if m.name.endswith((".fna.bz2", ".fna"))
                   and "VSG" not in os.path.basename(m.name)), None)
    if member is None:
        raise SystemExit(f"{archive}: no marker FASTA (.fna.bz2) in it")
    raw = tar.extractfile(member)
    return io.TextIOWrapper(bz2.open(raw) if member.name.endswith(".bz2") else raw)


def pr2_taxid(title, index, ncbi):
    """A PR2 title (record|domain|...|genus|species) -> (record, NCBI taxid of its species, else its genus) or None."""
    fields = title.split()[0].split("|")
    record, ranks = fields[0], fields[1:]
    if len(ranks) != len(labels.PR2_RANKS) or ranks[0] != "Eukaryota" or labels.ORGANELLE.search("|".join(ranks)):
        return None
    genus, species = ranks[-2], ranks[-1]
    if labels.is_species_name(species):
        t = index.get(labels.key(species))
        sp = ncbi.species_of(t) if t else None
        if sp and ncbi.rank.get(sp) == "species":
            return record, sp
    if labels.is_genus_name(genus):
        t = index.get(labels.key(genus))
        if t and ncbi.rank.get(t) == "genus":
            return record, t
    return None


def write(out, name, records):
    """records [(record, sequence, taxid)] -> out/<name>.fna.gz and out/<name>_taxids.tsv; how many."""
    os.makedirs(out, exist_ok=True)
    n = 0
    with gzip.open(os.path.join(out, f"{name}.fna.gz"), "wt") as fa, \
            open(os.path.join(out, f"{name}_taxids.tsv"), "w", newline="\n") as tab:
        for record, seq, taxid in records:
            fa.write(f">{record}\n{seq}\n")
            tab.write(f"{record}\t{taxid}\n")
            n += 1
    return n


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="what", required=True)
    m = sub.add_parser("metaphlan")
    m.add_argument("--archive", required=True, help="MetaPhlAn's database archive (.tar) or its marker FASTA")
    m.add_argument("--marker_info", required=True, help="its marker_info.txt.bz2")
    m.add_argument("--name", default="MetaPhlAn4")
    m.add_argument("--out", required=True)
    r = sub.add_parser("pr2")
    r.add_argument("--fasta", required=True, help="PR2's pr2_version_<V>_SSU_taxo_long.fasta.gz")
    r.add_argument("--taxdump", required=True)
    r.add_argument("--name", default="PR2")
    r.add_argument("--out", required=True)
    opts = p.parse_args(argv)
    if opts.what == "metaphlan":
        with text_of(opts.marker_info) as fh:
            taxa = euk_markers(fh)
        records = ((t.split()[0], s, taxa[t.split()[0]]) for t, s in fasta_records(marker_fasta(opts.archive))
                   if t.split()[0] in taxa)
        n = write(opts.out, opts.name, records)
        print(f"{opts.name}: {n:,} of MetaPhlAn's {len(taxa):,} eukaryote markers, {len(set(taxa.values())):,} taxa")
    else:
        ncbi = ncbi_mod.Ncbi(opts.taxdump)
        index = rrna.name_index(ncbi)
        kept = []
        with text_of(opts.fasta) as fh:
            for title, seq in fasta_records(fh):
                got = pr2_taxid(title, index, ncbi)
                if got:
                    kept.append((got[0], seq, got[1]))
        n = write(opts.out, opts.name, kept)
        print(f"{opts.name}: {n:,} records, {len({t for *_x, t in kept}):,} taxa")
    return 0


if __name__ == "__main__":
    sys.exit(main())
