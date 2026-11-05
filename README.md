# EUKDBbuild

This folder builds EUKDB, a standalone eukaryote reference. It holds the euk30 marker genes, complete with all their
introns, and the rDNA operons, complete and uncut, of every eukaryote species it can source. Both are cleaned and
labelled by NCBI's taxonomy, tied to no GTDB release or KSGP version. The package format is in `FORMAT.md`.

protal's database build converts a package to its own needs: KSGP names, copies fitted into 16 kb, records cut to
8,191 bases. Once the package is complete, protal downloads no eukaryote data of its own.

EUKDBbuild imports nothing from protal's `scripts/`. It started inside protal's repository on 2026-11-05, and is
meant to move into a repository of its own.

## Layout

| path | contents |
|---|---|
| `eukdb/util.py` | compressed files, downloads, MD5s |
| `eukdb/ncbi.py` | NCBI's taxonomy and the subset that ships in a package |
| `eukdb/sources.py` | EukDetect2 and EukDetect1, NCBI's assembly tables, EUKARYOME and ROD |
| `eukdb/species.py` | the model and host species (`data/model_species.tsv`) and the RefSeq species euk30's sources lack |
| `eukdb/placement.py` | euk30 loci placed in a genome with miniprot, from the species' own coding sequences or its relatives' proteins |
| `eukdb/exons.py` | exon models of EukDetect2's copies; the relatives index |
| `eukdb/copies.py` | a copy's exons and introns (and protal's 16 kb fitting, which EUKDB itself does not apply) |
| `eukdb/genomes.py` | the genome set: which genomes, one worker per genome (placement, features, rRNA), the prune |
| `eukdb/rrna.py` | a genome's 18S and whole rDNA units (inline 5S kept), published records by NCBI taxid, the checks |
| `eukdb/package.py` | collect, write, read and check a package |
| `eukdb_build.py` | the build |
| `data/` | the loci table, the model species table, curated rDNA units, the copies placed so far (`placed/`) |
| `tests/` | `python3 -m unittest discover -s EUKDBbuild/tests -t EUKDBbuild/tests` |

## Building

Three stages, so that the network-bound one runs where the network is:

```
python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage fetch                 # sources and genomes (network)
python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage process -t 64 -j 16   # placement, exons, rRNA; no network
python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage package --keep_gb 1000
```

`--stage all` (the default) runs the three in turn. The fetch stage gets EukDetect2 and EukDetect1, NCBI's assembly
tables and taxdump, EUKARYOME, ROD and the genomes. The process stage keeps every genome it was given; the package
stage prunes `OUT/genomes` to `--keep_gb` (the model species first, then the smallest) and writes
`OUT/EUKDB_<version>/` (`--version`, YYYY.MM by default). `--no_refseq` leaves out the RefSeq species EukDetect lacks
or holds on fewer than 20 loci (~1,700 genomes, ~735 GB gzipped). Each genome's results go into
`WORK/done/<accession>/`, so a rerun continues where the last one stopped.

`--species 4930,7227,...` (NCBI taxids, or a file of one per line) builds a trial package of the species inside those
clades alone: the genome set, the exon models, the published rRNA and the package are limited to them, and the
manifest's `inputs.partial` names the clades. The sources are fetched whole all the same.

A package is checked against its manifest by `eukdb.package.check(folder, deep=True)` (None when whole).

protal takes a package with `scripts/download_ref_data.py --eukdb OUT/EUKDB_<version>`.

## Not done yet
- The rRNA part screen (`records.tsv`'s `screen` column is `-`).

## Requirements
- Python 3.10 or later; 3.14 reads zstd itself, older versions need the `zstd` command.
- miniprot.
- barrnap with nhmmer.
- ITSx.
- vsearch.
