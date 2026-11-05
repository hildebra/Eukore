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
| `data/` | the loci table, the model species table, curated rDNA units, the copies placed so far (`placed/`) |
| `tests/` | `python3 -m unittest discover -s EUKDBbuild/tests -t EUKDBbuild/tests` |

The build (`eukdb_build.py`: the genome set, the rRNA records, the package writer) is being added.

## Requirements
- Python 3.10 or later; 3.14 reads zstd itself, older versions need the `zstd` command.
- miniprot.
- barrnap with nhmmer.
- ITSx.
- vsearch.
