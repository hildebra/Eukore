# EUKDB package format, version 4

EUKDB is a standalone eukaryote reference. It holds every BUSCO gene of eukaryota_odb12 and the rDNA operons of every eukaryote
species with an NCBI genome, complete (genes with all their introns, operons uncut) and checked. Its taxonomy is NCBI's,
tied to no GTDB release and no KSGP version. A user of the package (protal's build) takes what it needs (protal: the
30 BUSCOs of its euk30, chosen by `busco/loci.tsv`'s figures) and converts it: naming, length limits, masking,
record cuts. Paths below are relative to the package folder, `EUKDB_<YYYY.MM>/`. A
`.tsv` has a header line, tab-separated columns and no quoting. `-` means "none". Coordinates are 0-based, end
excluded, on the sequence as the package holds it, unless they are a contig's (1-based, inclusive).

## manifest.json
`format` (4), `version` (`YYYY.MM`), `built` (ISO time), `inputs` (BUSCO's dataset and its date, the taxdump's MD5,
EUKARYOME's and ROD's versions, NCBI's assembly tables' MD5s, the `groups` built, the no-assembly records' MD5;
`inputs.partial.species_clades` in a trial package limited to those clades, which a user should not take for EUKDB),
`counts` (species, species with rRNA alone, genomes, copies and their calls, rRNA records, taxa), and `files`: every
file of the package with `[bytes, md5]`. The manifest is written last: a package without one is unfinished.

## taxonomy/
NCBI's taxonomy, cut to every taxon on the package's lineages: `nodes.tsv` (taxid, parent, rank), `names.tsv` (taxid,
name, class: scientific name or synonym), `merged.tsv` (old, new). `eukdb.ncbi.Ncbi(folder)` reads it.

## species.tsv
One row per NCBI species the package holds:

| column | |
|---|---|
| `taxid` | NCBI species taxid |
| `name` | NCBI scientific name |
| `representative` | the genome its copies come from (an accession of `genomes.tsv`); `-` for a species with rRNA records alone |
| `source` | `ncbi` (NCBI's genome of the species by EUKDB's ranking), `model` (a model species' reference), `no_assembly` (records of a species NCBI has no assembly of), `rrna` (rRNA records alone) |
| `loci` | BUSCO genes called in its genome |
| `rrna_records` | its records in `rrna/records.tsv` |

## genomes.tsv
The species' genomes and their other assemblies:

| column | |
|---|---|
| `accession` | GenBank accession; `EUK_<taxid>` for a species' no-assembly records |
| `fetched` | the accession fetched: the RefSeq copy where the genome has one (RefSeq screens eukaryote genomes for contamination) |
| `species_taxid`, `taxid` | NCBI's |
| `organism`, `strain` | as NCBI names them |
| `level`, `genome_size`, `contigs`, `category` | the assembly's (`category`: isolate or derived from metagenome) |
| `flags` | NCBI's reasons it would not take the assembly into RefSeq (`excluded_from_refseq`), comma-separated |
| `ftp` | the assembly's folder on NCBI's FTP server |
| `source` | as `species.tsv`'s; `assembly` for a further assembly of the species (no copies) |
| `representative` | `t` for the species' genome |
| `loci` | BUSCO genes called in it |
| `local` | the genome's file in the build's kept genomes (`genomes/<accession>.fna.gz` beside the package), `-` where it is not kept |

A species' genome is NCBI's best by EUKDB's ranking: never one NCBI flags as contaminated, of an unverified source
organism, misassembled, chimeric or of low quality; another flagged one, or a metagenome's, only without better; then
a RefSeq reference, then a representative genome; then the better assembly level, the fewer contigs per megabase, the
newer.

## busco/
Every BUSCO gene of BUSCO's eukaryota_odb12 (129 at its 2026-05-22 release) called in every species' genome.
- `loci.tsv`: per BUSCO how widely the package's species hold it, for a user to choose its loci by:

| column | |
|---|---|
| `locus` | the BUSCO (an OrthoDB eukaryota id, `<group>at2759`) |
| `name` | OrthoDB's name of it |
| `cutoff`, `hmm_length` | BUSCO's score cutoff, the length of its HMM |
| `species` | species with a genome in the package |
| `called`, `duplicated`, `fragmented` | species whose genome holds it complete (Single or Duplicated); of them Duplicated; species with a fragment alone |
| `call_rate` | `called` over `species` |
| `call_<supergroup>` | the call rate within each supergroup of the package's species (Vertebrata, Metazoa, Fungi, Viridiplantae, Rhodophyta, SAR, Discoba, Metamonada, Amoebozoa, Haptista, Cryptophyceae, other eukaryotes) |
| `median_length`, `median_cds_length`, `median_protein_length` | of its complete copies |

- `copies/<locus>.fna.gz`: one record per genome, `>accession`, the copy whole: the genome from the gene's first codon
  to its last codon, all introns, on the gene's strand, upper case, nothing masked.
- `copies.tsv`: per copy how it was called and checked.

| column | |
|---|---|
| `accession`, `locus` | the copy |
| `status` | `Single` (one complete copy in the genome), `Duplicated` (several: the best is the copy), `Fragmented` (its HMM alignment covers less than 80% of the profile) |
| `contig`, `contig_length`, `start`, `end`, `strand` | where it lies in the assembly (1-based, inclusive) |
| `length`, `cds_length`, `exons` | of the copy; of its coding sequence; its exons |
| `identity`, `coverage` | miniprot's identity and the share of the query protein aligned |
| `hmm_score`, `cutoff`, `hmm_coverage` | its translation's score against the locus' BUSCO HMM, BUSCO's cutoff, the share of the profile aligned |
| `copies` | complete copies of the locus in the genome |
| `query` | the dataset protein that placed it (`<busco>_<taxid>_<id>`, or an ancestral variant) |
| `protein_length` | of its translation |
| `window_class` | tiara's class of the 50 kb around it (`eukarya`, `bacteria`, `archaea`, `organelle`, `unknown`), `-` where not classified |
| `fcs` | NCBI FCS-GX's call over the copy (`ACTION:division:top taxon`, its action `EXCLUDE`, `TRIM`, `FIX` or `REVIEW`), `-` for none or where not screened |
| `own_species` | `t`, `f` or `-`: its coding sequence's nearest copies of other species agree with its NCBI lineage to the family or deeper; or share no more than the phylum at 90% or more (another clade's gene: a contaminant, a misidentified genome) |
| `same_as` | the NCBI taxid of the lowest common ancestor of the species with a copy 99% identical to it, `-` for none |
| `length_check` | `short`, `long` (its protein in its locus' bottom or top 5%), `ok` |

- `features.tsv`: per copy its exons and introns, the assembly's soft-masked runs of 20 bases or more
  (`soft_masked`), and its low-complexity runs (`low_complexity`: DUST's score over 20 in 64-base windows).

| column | |
|---|---|
| `accession`, `locus` | the copy |
| `kind` | `exon`, `intron`, `soft_masked`, `low_complexity` |
| `start`, `end` | on the copy |

A copy failing a check is kept and flagged, so that a user chooses.

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
| `regions` | `name:start-end,...` in order on the record: `18S`, `ITS1`, `5.8S`, `ITS2`, `28S`, and `5S` where a 5S lies in the same block (the unit as one stretch, on either strand: some protists' and fungi's repeats hold it); a region not found is absent |
| `length` | of the record |
| `lsu_complete` | `t` where barrnap's 28S hit covers its whole model |
| `ambiguous`, `homopolymer` | the share of non-ACGT bases; the longest one-base run |
| `type` | its 99% 18S type within the species (`1`, `2`, ...) |
| `own_species` | `t`, `f`, or `-`: whether its nearest records of other sources agree with its NCBI lineage to the family or deeper (a genome's contaminant rDNA fails) |
| `screen` | the part screen's flags: `-`, or `part:rank,...`, each a part (`18S`, `junction`, `28S 5'`, `28S 3'`) whose nearest records of other genomes are another clade's at that rank (`kingdom` .. `family`) |
| `label_from` | how its taxid was found: `assembly`, `rod_taxid`, `eukaryome_name`, `eukaryome_genus` |

A record that fails a check is kept and flagged, so that a user chooses.

## genomes/ (beside the package, optional)
The eukaryote genomes the build kept, up to its `--keep_gb` (`<accession>.fna.gz`, as NCBI serves them). `genomes.tsv`'s
`local` column names them. A user simulating eukaryotes, or a host, takes them from here instead of downloading them.

## Changes
- 1 (2026-11-05): first format (EukDetect1's and EukDetect2's copies, placed copies of model and RefSeq species).
- 2 (2026-11-05): EUKDB's own locus calls in NCBI's genomes (BUSCO's eukaryota_odb12, miniprot, hmmsearch); one genome
  a species; `busco/copies.tsv` with the calls and checks; `features.tsv` without `method`, `protein`, `identity` and
  `coverage` (in `copies.tsv`), with `soft_masked` and `low_complexity`; copies unmasked; `genomes.tsv` with `fetched`
  and `flags`; `species.tsv`'s `source` in place of `sources`, and species with rRNA alone.
- 3 (2026-11-05): the organelle genes (`organelle/`; see that change's notes).
- 4 (2026-11-05): every BUSCO gene of the dataset, not the 30 of euk30 (the user picks its loci); `euk30/` renamed
  `busco/`; `busco/loci.tsv` of each BUSCO's calls over the package's species, in place of the euk30 loci table.
