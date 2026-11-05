#!/usr/bin/env python3
"""Exon models of EukDetect2's euk30 copies (since 2026-11-05; copied from protal's scripts/mini_db/euk30_exons.py).

EukDetect2's copies are genomic, introns and all, but come without their exons. This aligns proteins of each locus to a
genome's own copies with miniprot (placement.place: the copies are the "genome", one contig per locus) and takes each
copy's exons from its best alignment:

- the proteins (Relatives): the copy's own translation where it is one open frame (species.orf_cds), then the copies of
  the genome's nearest relatives by NCBI's lineage (the deepest shared ancestor first), QUERIES a locus: EukDetect1's
  coding sequences of the loci it maps, EukDetect2's copies that are one open frame, and the coding sequences of the
  model species' placed copies (their exons joined: mammals' for a mammal);
- a copy's model is the alignment of the highest identity times coverage on that copy (an alignment that lands on
  another locus's copy is not its model), taken at species.taken's thresholds; a copy without one has no model;
- EUKDB keeps every copy whole, of any length: its features are those of the whole copy (copies.fit uncapped), the
  parts outside the alignment as introns; protal fits a copy into its 16 kb when it converts the package.

Writes eukdetect2_exons.tsv (copies.FEATURES rows, record <genome>~<euk30 locus>, accession the EukDetect2 genome) into
--out. A genome's results are kept in --work as they are done, and a rerun annotates only the genomes not done yet.

    python3 -m eukdb.exons --eukdetect SOURCES/eukdetect2 --eukdetect1 SOURCES/eukdetect1 --work WORK/exons -j 16
"""

import argparse
import collections
import multiprocessing
import os
import re
import shutil
import sys
import tempfile

from . import copies as euk30_copies
from . import sources as eukaryotes
from . import species as model_species
from .ncbi import Ncbi
from .util import part_of

EXONS = "eukdetect2_exons.tsv"
FITTED = "eukdetect2_fitted.fna"
QUERIES = 4      # a locus's proteins aligned to a genome's copy
PER_NODE = 8     # copies of a locus kept per ancestor taxon, whose relatives are drawn from it


class Relatives:
    """Each locus's copies by every ancestor taxon of their genome (PER_NODE a taxon and locus), from which a genome's
    nearest relatives' proteins are drawn, the deepest shared ancestor first."""

    def __init__(self, ncbi):
        self.ncbi = ncbi
        self.at = collections.defaultdict(list)  # (ancestor taxid, locus) -> [(genome, kind, sequence)]
        self.orf = {}                              # (genome, locus) -> coding sequence or None

    def add(self, genome, taxid, kind, by_locus):
        """A genome's copies ({locus: sequence}), kind "cds" (a coding sequence) or "genomic" (taken where one open
        frame)."""
        for t in self.ncbi.ancestors(self.ncbi.merged.get(taxid, taxid)):
            for locus, seq in by_locus.items():
                lst = self.at[(t, locus)]
                if len(lst) < PER_NODE:
                    lst.append((genome, kind, seq))

    def cds(self, genome, kind, locus, seq):
        if kind == "cds":
            return seq
        key = (genome, locus)
        if key not in self.orf:
            self.orf[key] = model_species.orf_cds(seq)
        return self.orf[key]

    def queries(self, genome, taxid, copies, prefix):
        """{query record: coding sequence} for a genome's copies ({locus: sequence}; a locus without a copy: None):
        its own copy where one open frame, then its relatives', QUERIES a locus; a query's record is
        <prefix>-<locus>-S1~<n>~<source genome>."""
        out = {}
        lineage = list(self.ncbi.ancestors(self.ncbi.merged.get(taxid, taxid)))
        for locus, seq in sorted(copies.items()):
            taken, n = set(), 0
            own = [(genome, "genomic", seq)] if seq else []
            for g, kind, s in own + [c for t in lineage for c in self.at.get((t, locus), ())]:
                if g in taken:
                    continue
                cds = self.cds(g, kind, locus, s)
                if cds and len(cds) >= 90:
                    taken.add(g)
                    out[f"{prefix}-{locus}-S1~{n}~{g}"] = cds
                    n += 1
                    if n == QUERIES:
                        break
        return out


def relatives_of(ncbi, eukdetect, eukdetect1, loci, locus_map, placed_paths=(), copies2=None, link2=None):
    """The Relatives of every euk30 copy: EukDetect2's genomic copies (copies2/link2 where read already), EukDetect1's
    coding sequences of the loci locus_map maps (where its files are), the model species' placed copies' coding
    sequences (placed_paths) -> (Relatives, EukDetect2's (genomes, records), EukDetect2's copies)."""
    genomes, records = link2 or read_link2(eukdetect)
    if copies2 is None:
        copies2 = eukaryotes.read_loci(os.path.join(eukdetect, "eukdb.fasta"), records, set(loci))
    relatives = Relatives(ncbi)
    link1 = os.path.join(eukdetect1 or "", "busco_taxid_link.txt")
    if eukdetect1 and locus_map and os.path.isfile(link1):
        genomes1, records1 = eukaryotes.read_link1(link1, locus_map)
        for genome, by_locus in sorted(eukaryotes.read_loci(os.path.join(eukdetect1, "ncbi_eukprot_met_arch_markers.fna"),
                                                            records1).items()):
            relatives.add(genome, genomes1[genome][0], "cds", by_locus)
    for name, (taxid, by_locus) in sorted(model_sources([p for p in placed_paths if p]).items()):
        relatives.add(name, taxid, "cds", dict(by_locus))
    for genome in sorted(copies2):
        relatives.add(genome, genomes[genome][0], "genomic", copies2[genome])
    return relatives, (genomes, records), copies2


def read_link2(eukdetect):
    return eukaryotes.read_link(os.path.join(eukdetect, "busco_taxid_genome_link.txt"))


def annotate(job):
    """One genome: its copies' exon models and fitted copies -> (genome, [FEATURES rows], {locus: fitted copy over
    MAX_LOCUS}), written into WORK/done (genome.exons.tsv, genome.fna; the table last: it says the genome is done)."""
    genome, copies, queries, settings = job
    from . import placement as eukdetect1_genomic  # noqa: PLC0415
    stem = os.path.join(settings["work"], "done", re.sub(r"[^A-Za-z0-9_.-]", "_", genome))
    rows, fitted = [], {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "copies.fna")
        with open(fasta, "w") as fh:
            fh.writelines(f">{locus}\n{seq}\n" for locus, seq in sorted(copies.items()))
        results = eukdetect1_genomic.place(queries, fasta, settings["miniprot"], 1, 10 ** 9, 10 ** 9, model=True) \
            if queries else {}
    by_locus = collections.defaultdict(list)
    for query, (copy, fields, features) in results.items():
        locus = query.split("~", 1)[0].rsplit("-", 2)[1]
        if copy is not None and fields[1] == locus and \
                model_species.taken(float(fields[8]), float(fields[9]), 1.0):  # on its own locus's copy
            by_locus[locus].append((float(fields[8]) * float(fields[9]), query, features))
    for locus, hits in sorted(by_locus.items()):
        _score, _query, features = max(hits, key=lambda h: (h[0], h[1]))
        exons = [(int(r[8]), int(r[9])) for r in features if r[3] == "exon"]  # on the copy, 1-based inclusive
        seq = copies[locus]
        fit = euk30_copies.fit(seq, 1, len(seq), "+", exons, settings["max_len"], mask_min=10 ** 9)
        if fit["how"] != "whole":
            fitted[locus] = fit["copy"]
        rows += euk30_copies.feature_rows(f"{genome}~{locus}", genome, locus, "copy", "+", fit)
    os.makedirs(os.path.dirname(stem), exist_ok=True)
    with open(part_of(stem + ".fna"), "w", newline="\n") as fh:
        fh.writelines(f">{genome}~{locus}\n{seq}\n" for locus, seq in sorted(fitted.items()))
    os.replace(part_of(stem + ".fna"), stem + ".fna")
    euk30_copies.write_features(stem + ".exons.tsv", rows)
    return genome, rows, fitted


def model_sources(placed_paths, table=model_species.TABLE):
    """The model species' placed copies as coding sequences (their exons joined) -> {genome name: (taxid, {locus:
    cds})}."""
    rows = {model_species.genome_name(r): r for r in model_species.eukaryotes(model_species.read_table(table))}
    copies = model_species.copies_of(placed_paths)
    out = {}
    for path in placed_paths:
        for r in euk30_copies.read_features(euk30_copies.exons_of(path) if path else ""):
            m = eukaryotes.RECORD_RE.match(r[0])
            name = f"model_{m.group('genome')}" if m and m.group("group") == model_species.GROUP else None
            if name in rows and r[3] == "exon" and m.group("locus") in copies.get(name, {}):
                seq = copies[name][m.group("locus")]
                out.setdefault(name, (int(rows[name]["ncbi_taxid"]), collections.defaultdict(str)))
                out[name][1][m.group("locus")] += seq[int(r[4]):int(r[5])]
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eukdetect", required=True, help="EukDetect2's files (eukdb.fasta, busco_taxid_genome_link.txt, "
                                                      "taxdump.tar.gz)")
    p.add_argument("--eukdetect1", help="EukDetect1's files, whose coding sequences serve as relatives' proteins")
    p.add_argument("--loci", default=eukaryotes.LOCI_TABLE)
    p.add_argument("--placed", default=eukaryotes.PLACED_DIR,
                   help="the checkout's placed copies, whose model species serve as relatives (their exons joined)")
    p.add_argument("--work", required=True, help="each genome's results, as they are done")
    p.add_argument("--out", help=f"folder for {EXONS} and {FITTED} (default --eukdetect)")
    p.add_argument("--only", help="EukDetect2 genomes to annotate, comma-separated (default all)")
    p.add_argument("--miniprot", default="miniprot")
    p.add_argument("-j", "--jobs", type=int, default=4)
    opts = p.parse_args(argv)
    if not shutil.which(opts.miniprot):
        sys.exit(f"{opts.miniprot}: not found (bioconda's miniprot)")
    loci = eukaryotes.read_loci_table(opts.loci)
    ncbi = Ncbi(os.path.join(opts.eukdetect, "taxdump.tar.gz"))
    placed = [os.path.join(opts.placed, eukaryotes.PLACED_FASTA)] if opts.placed else []
    relatives, (genomes, _records), copies = relatives_of(ncbi, opts.eukdetect, opts.eukdetect1, loci,
                                                          eukaryotes.read_locus_map(opts.loci), placed)
    wanted = sorted(set(opts.only.split(",")) & set(copies)) if opts.only else sorted(copies)
    done = os.path.join(opts.work, "done")
    os.makedirs(done, exist_ok=True)
    settings = dict(work=opts.work, miniprot=opts.miniprot, max_len=10 ** 12)  # whole: EUKDB keeps every copy whole
    stems = {g: os.path.join(done, re.sub(r"[^A-Za-z0-9_.-]", "_", g)) for g in wanted}
    todo = [(g, copies[g], relatives.queries(g, genomes[g][0], copies[g], f"ed2-{g}"), settings)
            for g in wanted if not os.path.isfile(stems[g] + ".exons.tsv")]
    print(f"euk30 exon models: {len(wanted)} EukDetect2 genomes ({len(wanted) - len(todo)} done before)", flush=True)
    with multiprocessing.Pool(max(1, opts.jobs)) as pool:
        for i, (genome, rows, fitted) in enumerate(pool.imap_unordered(annotate, todo), 1):
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                print(f"  {i}/{len(todo)} genomes annotated", flush=True)
    out = opts.out or opts.eukdetect
    rows, fitted, models = [], [], collections.Counter()
    for g in wanted:
        got = euk30_copies.read_features(stems[g] + ".exons.tsv")
        rows += got
        models["copies with an exon model"] += len({r[0] for r in got})
        models["copies fitted into MAX_LOCUS"] += len({r[0] for r in got if r[11] != "whole"})
        with open(stems[g] + ".fna") as fh:
            fitted.append(fh.read())
    models["copies"] = sum(len(copies[g]) for g in wanted)
    models["copies over MAX_LOCUS"] = sum(len(s.strip("N")) > eukaryotes.MAX_LOCUS for g in wanted
                                          for s in copies[g].values())
    euk30_copies.write_features(os.path.join(out, EXONS), rows)
    with open(part_of(os.path.join(out, FITTED)), "w", newline="\n") as fh:
        fh.write("".join(fitted))
    os.replace(part_of(os.path.join(out, FITTED)), os.path.join(out, FITTED))
    print(f"-> {os.path.join(out, EXONS)}, {os.path.join(out, FITTED)}: " +
          ", ".join(f"{n} {k}" for k, n in models.items()), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
