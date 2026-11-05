# EUKDB package format, version 1

EUKDB is a standalone eukaryote reference. It holds the euk30 marker genes and the rDNA operons of every eukaryote
species it can source, complete (genes with all their introns, operons uncut) and cleaned. Its taxonomy is NCBI's,
tied to no GTDB release and no KSGP version. A user of the package (protal's build) takes what it needs and converts
it: naming, length limits, record cuts. Paths below are relative to the package folder, `EUKDB_<YYYY.MM>/`. A `.tsv`
has a header line, tab-separated columns and no quoting. `-` means "none". Coordinates are 0-based, end excluded,
on the sequence as the package holds it, unless named `genome_*` (1-based, inclusive, on the assembly's contig).

## manifest.json
`format` (1), `version` (`YYYY.MM`), `built` (ISO time), `inputs` (each source with its version, URL and MD5s, and the
date of NCBI's tables and taxdump), `counts` (species, genomes, copies, records), and `files`: every file of the
package with `[bytes, md5]`. The manifest is written last: a package without one is unfinished.

## taxonomy/
NCBI's taxonomy, cut to every taxon on the package's lineages: `nodes.tsv` (taxid, parent, rank), `names.tsv` (taxid,
name, class: scientific name or synonym), `merged.tsv` (old, new). `eukdb.ncbi.Ncbi(folder)` reads it.

## species.tsv
One row per NCBI species the package holds:

| column | |
|---|---|
| `taxid` | NCBI species taxid |
| `name` | NCBI scientific name |
| `representative` | the genome whose copies represent the species (an accession of `genomes.tsv`) |
| `loci` | euk30 loci its representative holds |
| `sources` | where its copies came from, comma-separated: `eukdetect2`, `eukdetect1`, `model`, `refseq` |
| `rrna_records` | its records in `rrna/records.tsv` |

## genomes.tsv
One row per genome a copy or a record comes from:

| column | |
|---|---|
| `accession` | GenBank or RefSeq assembly accession; `EUK_<name>` for an EukDetect genome without one |
| `species_taxid`, `taxid` | NCBI's |
| `organism`, `strain` | as NCBI names them |
| `level`, `genome_size`, `contigs`, `category` | the assembly's (`category`: isolate or derived from metagenome) |
| `ftp` | the assembly's folder on NCBI's FTP server, `-` without one |
| `source` | `eukdetect2`, `eukdetect1`, `model`, `refseq`, `assembly` (a further assembly of the species) |
| `representative` | `t` for the species' representative |
| `loci` | euk30 loci this genome has copies of |
| `local` | the genome's file in the build's kept genomes (`genomes/<accession>.fna.gz` beside the package), `-` where it is not kept |

## euk30/
- `loci.tsv`: `locus` (the euk30 locus, an OrthoDB eukaryota id), `eukdetect1_locus`, `median_length`.
- `copies/<locus>.fna.gz`: one record per genome, `>accession`, the copy whole: the genome from the gene's first codon
  to its last, all introns, on the gene's strand, upper case, soft-masked repeats of 100 bases or more outside the exons
  as N. EukDetect2's copies are as EukDetect2 gives them, flanks included.
- `features.tsv`: per copy its exons and introns.

| column | |
|---|---|
| `accession`, `locus` | the copy |
| `kind` | `exon` or `intron` |
| `start`, `end` | on the copy |
| `genome_contig`, `genome_start`, `genome_end` | where it lies in the assembly (placed copies), `-` otherwise |
| `strand` | the gene's strand on the contig |
| `method` | `placed_own` (its own coding sequence), `placed_relative` (a relative's protein), `annotated` (an EukDetect2 copy's model), `coding` (an EukDetect1 coding copy: one exon over it) |
| `protein` | the genome whose protein gave the model |
| `identity`, `coverage` | miniprot's |

A copy without a row has no exon model.

## rrna/
- `records.fna.gz`: every record, `>record`, uncut.
- `records.tsv`:

| column | |
|---|---|
| `record` | its name: `<accession>_rDNA_<n>`, `<accession>_18S_<n>`, `EUK_<EUKARYOME id>`, `ROD_<variant id>` |
| `kind` | `rDNA` (an operon) or `18S` (an 18S span alone) |
| `taxid`, `species_taxid` | NCBI's |
| `accession` | its genome (`-` for a published read without one) |
| `source` | `genome`, `curated_unit`, `EUKARYOME`, `ROD` |
| `regions` | `name:start-end,...` in order on the record: `5'ETS`, `18S`, `ITS1`, `5.8S`, `ITS2`, `28S`, and `5S` where a 5S lies in the same block (the unit as one stretch, on either strand: some protists' and fungi's repeats hold it); a region not found is absent |
| `length` | of the record |
| `lsu_complete` | `t` where barrnap's 28S hit covers its whole model |
| `ambiguous`, `homopolymer` | the share of non-ACGT bases; the longest one-base run |
| `type` | its 99% 18S type within the species (`1`, `2`, ...) |
| `own_species` | `t`, `f`, or `-`: whether its nearest records of other sources agree with its NCBI lineage to the family or deeper (a genome's contaminant rDNA fails) |
| `screen` | the part screen's verdict: `-`, or the rank a part's nearest records disagree at |
| `label_from` | how its taxid was found: `assembly`, `rod_taxid`, `eukaryome_name`, `eukaryome_genus` |

A record that fails a check is kept and flagged, so that a user chooses.

## genomes/ (beside the package, optional)
The eukaryote genomes the build kept, up to its `--keep_gb` (`<accession>.fna.gz`, as NCBI serves them). `genomes.tsv`'s
`local` column names them. A user simulating eukaryotes, or a host, takes them from here instead of downloading them.

## Changes
- 1 (2026-11-05): first format.
