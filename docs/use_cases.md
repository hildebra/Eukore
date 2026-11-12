# Using a package: four short cases

Eukore is a database: what it is good for shows in its uses. Each case below reads a finished package with Eukore's
own reader (`eukore/package.py`; the tables are plain TSV, so R, pandas or awk do as well). Run them from the
checkout, with `PKG` a package folder. The examples name the group `Alveolata` and the primer pair `V4_TAReuk`: put
yours.

## 1. Choose loci for a clade

`busco/loci.tsv` gives each BUSCO's call rate over all species and per group (`call_<group>`, `species.tsv`'s
`group`; [loci.md](loci.md)). The loci called in most species of the clade are the ones to build on:

```python
from eukore import package
PKG = "Eukore_2026.11"
loci = package.read_table(PKG, "busco/loci.tsv")
column = "call_Alveolata"
chosen = [r["locus"] for r in sorted(loci, key=lambda r: -float(r[column]))[:30]]
for r in sorted(loci, key=lambda r: -float(r[column]))[:5]:
    print(r["locus"], r["name"], r[column])
```

Which identity then tells the clade's species apart is in the benchmarks' `thresholds.tsv`
([benchmarks/README.md](../benchmarks/README.md#guidance-the-identity-that-parts-taxa)): one row per BUSCO and rank.

## 2. Check a primer pair before ordering it

`rrna/amplicons.tsv` holds the in-silico PCR of each primer pair of `data/primer_pairs.tsv` on every rRNA record
(`organelle/amplicons.tsv` the organelle barcodes'). The species of a group a pair amplifies:

```python
records = {r["record"]: r for r in package.read_table(PKG, "rrna/records.tsv")}
species = {r["taxid"]: r for r in package.read_table(PKG, "species.tsv")}
tried, amplified = set(), set()
for r in package.read_table(PKG, "rrna/amplicons.tsv"):
    if r["pair"] != "V4_TAReuk":
        continue
    taxid = records[r["record"]]["species_taxid"]
    if species.get(taxid, {}).get("group") == "Alveolata":
        tried.add(taxid)
        if r["amplifies"] == "t":
            amplified.add(taxid)
print(f"{len(amplified)} of {len(tried)} species amplified")
missed = sorted(species[t]["name"] for t in tried - amplified)
```

A row's mismatch counts (`forward_mismatches`, `reverse_3prime`, ..) let you apply your own rule. Every pair's
share per group is the package's `benchmarks/primers.tsv`.

## 3. A tree of the genomes

The benchmarks' tree step joins the BUSCO proteins of up to two species per genus and runs FastTree
([benchmarks/README.md](../benchmarks/README.md#how-good-the-records-are)). For your own selection, the proteins of
a locus are one call away, and `eukore/align.py` places them on BUSCO's profile as the build does:

```bash
PYTHONPATH=benchmarks python3 -m eukore_bench.tree --package $PKG --odb SOURCES/eukaryota_odb12 --out TREE -t 16
```

```python
proteins = package.read_copies(PKG, chosen[0], "proteins")     # {accession: protein}
```

`TREE/tree.nwk` is the tree (its tips the genomes' accessions, `tree_tips.tsv` names them), `tree_taxa.tsv` says
which NCBI taxa are one branch of it, `tree_strays.tsv` which genomes sit apart from their name's relatives.

## 4. A reference for mapping reads

A read mapper needs one FASTA and each record's species. Keep the copies whose checks found nothing against them
(`confidence` `high`), or also those that may be cut or gapped (`medium`); a `low` copy is another organism's:

```python
genomes = {r["accession"]: r["species_taxid"] for r in package.read_table(PKG, "genomes.tsv")
           if r["representative"] == "t"}
keep = {(r["accession"], r["locus"]) for r in package.read_table(PKG, "busco/copies.tsv")
        if r["confidence"] in ("high", "medium") and r["locus"] in chosen}
with open("reference.fna", "w") as fasta, open("reference_taxids.tsv", "w") as table:
    for locus in chosen:
        for accession, sequence in package.read_copies(PKG, locus).items():   # the copy, introns and all
            if (accession, locus) in keep and accession in genomes:
                fasta.write(f">{accession}|{locus}\n{sequence}\n")
                table.write(f"{accession}|{locus}\t{genomes[accession]}\n")
```

```bash
minimap2 -x sr -c --secondary=yes reference.fna reads.fastq.gz > reads.paf
```

`read_copies(PKG, locus, "cds")` gives the coding sequences instead (an exon-only reference); the rRNA records are
`rrna/records.fna.gz`. How a read's hits become a taxon is the user's own rule, or a profiler's (protal builds its
eukaryote database from a package this way): Eukore gives the sequences and what is known of each.
