#!/usr/bin/env python3
"""Other references for the comparisons and the read benchmark (compare.py, reads.py; until 2026-11-12
tools/compare_references.py, for the build's --read_compare), each as a FASTA and its record<TAB>NCBI taxid table
(since 2026-11-06). prep.py runs it on its downloads; EukDetect2's database needs none (its eukdb.fasta and
busco_taxid_genome_link.txt are such a pair):

- metaphlan: MetaPhlAn 4's eukaryote markers (Blanco-Miguez et al. 2023), from its database's archive
  (mpa_<version>.tar, or its marker FASTA, .fna.bz2) and its marker_info.txt.bz2: the markers whose clade is t__EUK<taxid>
  (an NCBI taxid), from https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/.
- pr2: PR2's 18S (Guillou et al. 2013; pr2_version_<V>_SSU_taxo_long.fasta.gz, its 9 ranks in each title): a record's
  species where NCBI names it (its name or a synonym, as species.tsv's labels find them), else its genus; organelle
  records (plastid, mitochondrion, nucleomorph) and those naming neither are left out.
- silva, unite (since 2026-11-17): SILVA's SSU Ref NR99 eukaryote 18S (its RNA written as DNA) and UNITE's ITS (its
  general release of all eukaryotes), the files the build's labels read (SOURCES/labels), each record's species where
  NCBI names it, else its genus: the agreement of their 18S and ITS with the genomes' (databases.py).

    python3 -m eukore_bench.references metaphlan --archive mpa_vJan26_CHOCOPhlAnSGB_202605.tar \\
        --marker_info mpa_vJan26_CHOCOPhlAnSGB_202605_marker_info.txt.bz2 --out DIR
    python3 -m eukore_bench.references pr2 --fasta pr2_version_5.1.1_SSU_taxo_long.fasta.gz \\
        --taxdump Eukore_2026.11/taxonomy --out DIR

DIR gets <name>.fna.gz and <name>_taxids.tsv (compare.py's references.tsv lists them).
"""

import argparse
import bz2
import gzip
import io
import os
import re
import sys
import tarfile

from . import ROOT  # noqa: F401  (the checkout on the path: the database's eukore package)
from eukore import labels, ncbi as ncbi_mod, rrna, util  # noqa: E402

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
    """A PR2 title (record|gene|organelle|...|domain|...|genus|species: its 9 ranks last, as taxo_long writes it) ->
    (record, NCBI taxid of its species, else its genus) or None (an organelle's record, or neither named)."""
    fields = title.split()[0].split("|")
    record, ranks = fields[0], fields[-len(labels.PR2_RANKS):]
    where = fields[2] if len(fields) > len(labels.PR2_RANKS) + 2 else "nucleus"
    if len(fields) <= len(labels.PR2_RANKS) or ranks[0] != "Eukaryota" or where != "nucleus" or \
            labels.ORGANELLE.search("|".join(ranks)):
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


def named_taxid(species, genus, index, ncbi):
    """A record's species name (else its genus) -> its NCBI taxid (the species', else the genus') or None."""
    if labels.is_species_name(species):
        t = index.get(labels.key(species))
        sp = ncbi.species_of(t) if t else None
        if sp and ncbi.rank.get(sp) == "species":
            return sp
    if labels.is_genus_name(genus):
        t = index.get(labels.key(genus))
        if t and ncbi.rank.get(t) == "genus":
            return t
    return None


def silva_records(fh, index, ncbi):
    """SILVA's SSU Ref NR99 FASTA (>accession.start.end Eukaryota;...;organism; RNA) -> (record, DNA, taxid) of its
    eukaryote records naming an NCBI species or genus (since 2026-11-17). Its plastids and mitochondria are under
    Bacteria: a Eukaryota record is a nuclear 18S."""
    for title, seq in fasta_records(fh):
        record, _sp, tax = title.partition(" ")
        parts = tax.split(";")
        if parts[0] != "Eukaryota" or len(parts) < 3:
            continue
        species = labels.species_of_name(parts[-1])
        genus = (parts[-1].split() or [""])[0]  # the organism's first word: "Saccharomyces sp. ABC" its genus
        t = named_taxid(species, genus, index, ncbi)
        if t:
            yield record, seq.upper().replace("U", "T"), t


def unite_records(lines, index, ncbi):
    """UNITE's general FASTA release (>name|accession|SH|refs or reps|k__..;s__Genus_species; ITSx's ITS1-5.8S-ITS2)
    -> (record, sequence, taxid) of those naming an NCBI species or genus (since 2026-11-17). A record's name is its
    accession and SH, unique where an accession is in two SHs."""
    for title, seq in fasta_records(lines):
        f = title.split("|")
        if len(f) < 5:
            continue
        ranks = {p.split("__", 1)[0]: p.split("__", 1)[1] for p in f[4].split(";") if "__" in p}
        t = named_taxid(ranks.get("s", "").replace("_", " "), ranks.get("g", ""), index, ncbi)
        if t:
            yield f"{f[1]}_{f[2]}", seq.upper(), t


def write(out, name, records):
    """records [(record, sequence, taxid)] -> out/<name>.fna.gz and out/<name>_taxids.tsv; how many."""
    os.makedirs(out, exist_ok=True)
    n = 0
    with util.bgzf_open(os.path.join(out, f"{name}.fna.gz")) as fa, \
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
    for what, name, given in (("silva", "SILVA", "SILVA's SSU Ref NR99 FASTA (labels.SILVA_FILE)"),
                              ("unite", "UNITE", "UNITE's general release of all eukaryotes (labels.UNITE_FILE, .tgz)")):
        s = sub.add_parser(what)
        s.add_argument("--fasta", required=True, help=given)
        s.add_argument("--taxdump", required=True)
        s.add_argument("--name", default=name)
        s.add_argument("--out", required=True)
    opts = p.parse_args(argv)
    if opts.what in ("silva", "unite"):  # since 2026-11-17: the agreement of their 18S and ITS with the genomes'
        ncbi = ncbi_mod.Ncbi(opts.taxdump)
        index = rrna.name_index(ncbi)
        if opts.what == "silva":
            with text_of(opts.fasta) as fh:
                kept = list(silva_records(fh, index, ncbi))
        else:
            kept = list(unite_records(labels.unite_lines(opts.fasta) if tarfile.is_tarfile(opts.fasta)
                                      else text_of(opts.fasta), index, ncbi))
        n = write(opts.out, opts.name, kept)
        print(f"{opts.name}: {n:,} records, {len({t for *_x, t in kept}):,} taxa")
        return 0
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
