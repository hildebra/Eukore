# Eukore package format, version 4

Eukore is a standalone eukaryote reference. It holds every BUSCO gene of eukaryota_odb12 and the rDNA operons of every
eukaryote species with an NCBI genome, complete (genes with all their introns, operons uncut) and checked, and the
species' mitochondrial and plastid rRNA and mitochondrial COX1 where NCBI has them. Its taxonomy is NCBI's, tied to no
GTDB release; each species also carries its lineage in PR2, EUKARYOME and KSGP (`species.tsv`), the versions named in
the manifest. Each copy and record carries its checks, a confidence tier and the reasons for it, and the hints a user
needs to spot a chimera (contig ends, assembly gaps, the screens). A user of the package (protal's build) takes what
it needs (protal: the 30 BUSCOs of its euk30, chosen by `busco/loci.tsv`'s figures) and converts it: naming, length
limits, masking, record cuts. Paths below are relative to the package folder, `Eukore_<YYYY.MM>/`. A
`.tsv` has a header line, tab-separated columns and no quoting. `-` means "none". Coordinates are 0-based, end
excluded, on the sequence as the package holds it, unless they are a contig's (1-based, inclusive).

## manifest.json
`format` (4), `version` (`YYYY.MM`), `eukore` (the version of Eukore that built it, `0.1`), `built` (ISO time),
`inputs` (BUSCO's dataset and its date, the taxdump's MD5,
EUKARYOME's and ROD's versions, NCBI's assembly tables' MD5s, the `groups` built, the no-assembly records' MD5, RefSeq's
organelle release files' MD5s (`refseq_organelle`), NCBI's FCS-GX and ANI reports' MD5s (`ncbi_quality`), the 18S
references species.tsv's lineages come from (`references`: PR2's, EUKARYOME's, KSGP's version, file, MD5 and the
species and genera it names); `inputs.partial.species_clades` in a trial package limited to those clades, which a user
should not take for Eukore), `counts` (species, species with rRNA alone, genomes, copies and their
calls, rRNA records, organelle genes by organelle and gene, species with organelle genes, taxa, `busco_loci`, the
copies' and records' confidence tiers, `labels`: per reference the species by how they matched), `benchmarks` (the
headline figures of `benchmarks/summary.json`, where the build measured them), and `files`: every file of the package
with `[bytes, md5]`. The manifest is written last: a package without one is unfinished.

## taxonomy/
NCBI's taxonomy, cut to every taxon on the package's lineages: `nodes.tsv` (taxid, parent, rank, gencode: NCBI's nuclear genetic code, the table a copy translates by), `names.tsv` (taxid,
name, class: scientific name or synonym), `merged.tsv` (old, new). `eukore.ncbi.Ncbi(folder)` reads it.

A taxid is NCBI's of the taxdump the manifest names (`inputs.taxdump_md5`). A taxid NCBI has merged into another since
is in `merged.tsv`: a user mapping older labels (or a later taxdump's) to the package goes through it.

## species.tsv
One row per NCBI species the package holds:

| column | |
|---|---|
| `taxid` | NCBI species taxid |
| `name` | NCBI scientific name |
| `representative` | the genome its copies come from (an accession of `genomes.tsv`); `-` for a species with rRNA records alone |
| `source` | `ncbi` (NCBI's genome of the species by Eukore's ranking), `model` (a model species' reference), `no_assembly` (records of a species NCBI has no assembly of), `rrna` (rRNA records alone) |
| `loci` | BUSCO genes called in its genome |
| `rrna_records` | its records in `rrna/records.tsv` |
| `organelle_genes` | its records in `organelle/genes.tsv` |
| `supergroup` | `Metazoa`, `Fungi`, `Viridiplantae`, `Rhodophyta`, `SAR`, `Discoba`, `Metamonada`, `Amoebozoa`, `Haptista`, `Cryptophyceae` (the first NCBI clade of these that holds it), else `other` |
| `ncbi_lineage` | NCBI's names at domain;kingdom;phylum;class;order;family;genus;species (a rank NCBI does not name empty) |
| `pr2`, `pr2_match` | its lineage in PR2 (domain;supergroup;division;subdivision;class;order;family;genus;species, as PR2 writes it: `Saccharomyces_cerevisiae`), and how it was found (below) |
| `eukaryome`, `eukaryome_match` | its lineage in EUKARYOME (superkingdom;kingdom;subkingdom;phylum;subphylum;class;order;family;genus;species: its classification table to the genus, its long reads' species) |
| `ksgp`, `ksgp_match` | its lineage in KSGP (`k__` .. `s__`, as KSGP writes it; its eukaryote nuclear records) |

A species is found in a reference by name: `species` (NCBI's scientific name), `synonym` (one of NCBI's synonyms, an
anamorph), `genus` (neither, but the reference names its genus: the lineage down to the genus, no species), `-` (none:
the lineage `-`). Names are compared without case, `_` as a space and a `(Fungi)` after a genus. A reference's records
of an organelle (PR2's `:plas`, `:mito`; KSGP's chloroplast, `plas` and `mito` names) and its names that name no species
(`sp.`, `Genus_X`, accessions, `?`) are not taken. A name a reference files under several lineages (a homonym of
another kingdom) takes the lineage that shares the most names with NCBI's lineage of the species, then the one of the
most records. The labels are a name match, not a placement: a user that names species by its own rules (protal's KSGP
rules) can compare its name with these. A row's taxid may be above species (a published read labelled by its genus).

## genomes.tsv
The species' genomes and their other assemblies:

| column | |
|---|---|
| `accession` | GenBank accession; `EUK_<taxid>` for a species' no-assembly records |
| `fetched` | the accession fetched: the RefSeq copy where the genome has one (RefSeq screens eukaryote genomes for contamination) |
| `species_taxid`, `taxid` | NCBI's |
| `organism`, `strain` | as NCBI names them |
| `level`, `genome_size`, `contigs`, `category` | the assembly's (`category`: isolate or derived from metagenome) |
| `flags` | NCBI's reasons it would not take the assembly into RefSeq (`excluded_from_refseq`), and NCBI's quality reports' (`FCS-GX: mostly foreign`, `FCS-GX: contaminated`, `ANI: failed`, `ANI: inconclusive`; see below), comma-separated |
| `ftp` | the assembly's folder on NCBI's FTP server |
| `source` | as `species.tsv`'s; `assembly` for a further assembly of the species (no copies) |
| `representative` | `t` for the species' genome |
| `loci` | BUSCO genes called in it |
| `rrna_18s`, `rdna_units`, `rdna_complete`, `rdna_contigs` | its 18S spans, rDNA units (from the 18S to the 28S), units with a whole 28S, and the contigs holding units, as `rrna/` found them in the assembly: copies in the assembly, a lower bound of the copies in the genome (an rDNA array is usually collapsed in an assembly); `-` where its rRNA was not searched |
| `local` | the genome's file in the build's kept genomes (`genomes/<accession>.fna.gz` beside the package), `-` where it is not kept |

A species' genome is NCBI's best by Eukore's ranking: never one NCBI flags as contaminated, of an unverified source
organism, misassembled, chimeric or of low quality, nor one NCBI's FCS-GX screen finds 50% or more another organism's
(`FCS-GX: mostly foreign`) or that fails NCBI's ANI species check (`ANI: failed`; Saccharomycotina); another flagged
one, or a metagenome's, only without better (`FCS-GX: contaminated`: 1% or more contaminant, or 100 kb or more of
another eukaryote's sequence; `ANI: inconclusive`); then a RefSeq reference, then a representative genome; then the
better assembly level, one FCS-GX finds no contaminant in, the fewer contigs per megabase, the newer. FCS-GX's summary
is the fetched assembly's (its RefSeq copy where it has one).

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
| `contig_edge` | bases from the copy to its contig's nearer end: a copy at an end may be cut, or joined to another sequence |
| `gaps` | runs of 10 Ns or more in the copy (assembly gaps: a misjoin's mark) |
| `confidence`, `flags` | `high`, `medium` or `low`, and why (comma-separated, `-` for none; below) |

A copy is `low` where a check calls it another organism's: `own_species` (f), `fcs` (FCS-GX's EXCLUDE, TRIM or FIX),
`window_prokaryote` (tiara's bacteria or archaea). It is `medium` where something may be wrong: `fragmented`,
`duplicated`, `length_short` or `length_long`, `gap`, `contig_edge` (under 1,000 bases), `own_species_undecided`
(`own_species` -), `window_organelle`, `fcs_review`. Else `high`.

- `features.tsv`: per copy its exons and introns, the assembly's soft-masked runs of 20 bases or more
  (`soft_masked`), its low-complexity runs (`low_complexity`: DUST's score over 20 in 64-base windows) and its
  assembly gaps (`gap`: runs of 10 Ns or more).

| column | |
|---|---|
| `accession`, `locus` | the copy |
| `kind` | `exon`, `intron`, `soft_masked`, `low_complexity`, `gap` |
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
| `contig`, `start`, `end`, `strand` | a genome's record: where it lies in the assembly (1-based, inclusive; the 18S span's, or the unit's from its first region to its last); `-` for a published read |
| `contig_edge` | a genome's record: bases to its contig's nearer end |
| `introns` | the 18S's introns on the record (`18S:start-end,...`: the gaps between barrnap's 18S pieces), `-` for none; the 18S region spans them |
| `gaps` | runs of 10 Ns or more in the record |
| `confidence`, `flags` | `high`, `medium` or `low`, and why (below) |

A record is `low` where a check calls it another organism's or of bad sequence: `own_species` (f), `screen`,
`ambiguous` (over 0.5%), `homopolymer` (36 or longer). It is `medium` where something may be wrong: `gap`, `contig_edge`
(under 1,000 bases), `own_species_undecided`, `lsu_partial` (an operon without a whole 28S), `genus_label` (a
published read labelled by its genus alone). Else `high`.

An operon's spacers are its ITS1 and ITS2 (`regions`); the intergenic spacer between units is not taken (a unit runs
from its 18S to its 28S, and its 5S where that lies in the same block). A 5S elsewhere in the genome (most animals' and
plants' 5S arrays) is not searched.

A record that fails a check is kept and flagged, so that a user chooses.

## organelle/
Each species' mitochondrial rRNA, plastid rRNA and mitochondrial COX1, where NCBI has them. Of each organelle, a
species' genes come from one source: its RefSeq organelle genome (the one with the most of these genes, then of the
species' own taxid, then the longest), else the organelle molecules of its own assembly (those NCBI's FASTA titles
name a mitochondrion, chloroplast, plastid or apicoplast).
- `genes.fna.gz`: one record per gene, `>record`, the gene from its first base to its last on its strand, introns and
  all, upper case.
- `genes.tsv`:

| column | |
|---|---|
| `record` | its name: `<seqid>_<gene>` |
| `gene` | `SSU` (12S in animals' mitochondria, 16S in plastids, 18S in plants' mitochondria, ...), `LSU` (16S in animals' mitochondria, 23S, 26S, ...), `5S`, `COX1` |
| `organelle` | `mitochondrion` or `plastid` (chloroplasts, apicoplasts and other plastids) |
| `taxid`, `species_taxid` | NCBI's |
| `source` | `refseq` (RefSeq's organelle genome, its annotated features) or `assembly` (the species' assembly: each gene by nhmmer with that gene of the 3 nearest RefSeq species, its exons joined; rRNA nhmmer misses by barrnap, `--kingdom mito` in animals' mitochondria, `bac` elsewhere) |
| `seqid` | the organelle sequence (RefSeq's accession, or the assembly's molecule) |
| `accession` | its assembly (`genomes.tsv`), `-` for RefSeq's |
| `start`, `end`, `strand` | on `seqid` (1-based, inclusive); `start` after `end` where the gene crosses a circular genome's origin |
| `length` | of the record |
| `exons` | RefSeq's feature's parts (a COX1's exons), or the nhmmer hits joined; 1 for barrnap's |
| `copies` | the gene's copies on the source (a plastid's inverted repeat: 2): the record is the first, or barrnap's longest whole hit |
| `partial` | `t` where RefSeq marks the feature partial, barrnap's hit is partial, or nhmmer's hits cover under 80% of the query's exons |
| `ambiguous` | the share of non-ACGT bases |

## benchmarks/ (optional)
How complete, right and discriminating the package is, measured by the build (`eukore/benchmarks.py`; absent with
`--no_benchmarks`). docs/benchmarks.md describes each table and `docs/benchmarks/plot_benchmarks.py` draws them. A
benchmark whose inputs the build lacked (BUSCO's dataset, vsearch) is left out, and `summary.json`'s `left_out` says
why. A reader needs none of them.

| file | one row per |
|---|---|
| `completeness.tsv` | species' genome: its BUSCOs Single, Duplicated, Fragmented, Missing, `complete` (Single and Duplicated over the BUSCOs called), its copies' confidence tiers, its assembly's level, size and contigs, its rRNA counts |
| `orthodb.tsv` | BUSCO of a species BUSCO's dataset holds a protein of: Eukore's `status`; `recall` and `precision`, the shares of OrthoDB's protein's 5-mers in the call's translation and the reverse; `guided` (`t` where the call's query was that species' own protein); `verdict`: `concordant` (both 0.9 or more), `similar` (both 0.5 or more), `discordant`, `missed`; `held_out_status`, `held_out_query`, `held_out_recall`, `held_out_precision`, `held_out_verdict`: the same for the species' BUSCOs called again without its own proteins (`-` where the build made no such calls) |
| `resolution.tsv` | marker (a BUSCO; `BUSCO (all)`; `18S`, each species' best record's) and `shared_rank` (the rank a copy shares with its nearest copy of another species: `genus` .. `phylum`, `none`; `all`; `unmatched`: none 70% identical): `copies`, those `under_97` and `under_99` identical, the identity's quantiles (`q05` .. `q95`) |
| `rrna.tsv` | rRNA source and kind: records, median length, records holding each region, a whole 28S, an 18S intron, each confidence tier |
| `labels.tsv` | supergroup (and `all`) and reference: its species by how they were named (`species`, `synonym`, `genus`, `none`) |
| `runtime.tsv` | genome timed (run since 2026-11-06): its size, bytes, BUSCOs searched, threads, seconds |
| `stages.tsv` | step of the package stage: its seconds |
| `summary.json` | the headline figures |

## genomes/ (beside the package, optional)
The eukaryote genomes the build kept, up to its `--keep_gb` (`<accession>.fna.gz`, as NCBI serves them). `genomes.tsv`'s
`local` column names them. A user simulating eukaryotes, or a host, takes them from here instead of downloading them.

## Changes
- 1 (2026-11-05): first format (EukDetect1's and EukDetect2's copies, placed copies of model and RefSeq species).
- 2 (2026-11-05): Eukore's own locus calls in NCBI's genomes (BUSCO's eukaryota_odb12, miniprot, hmmsearch); one genome
  a species; `busco/copies.tsv` with the calls and checks; `features.tsv` without `method`, `protein`, `identity` and
  `coverage` (in `copies.tsv`), with `soft_masked` and `low_complexity`; copies unmasked; `genomes.tsv` with `fetched`
  and `flags`; `species.tsv`'s `source` in place of `sources`, and species with rRNA alone.
- 3 (2026-11-05): `organelle/` (mitochondrial and plastid rRNA, mitochondrial COX1); `species.tsv`'s
  `organelle_genes`; the manifest's `inputs.refseq_organelle` and organelle counts.
- 4 (2026-11-05): every BUSCO gene of the dataset, not the 30 of euk30 (the user picks its loci); `euk30/` renamed
  `busco/`; `busco/loci.tsv` of each BUSCO's calls over the package's species, in place of the euk30 loci table;
  `taxonomy/nodes.tsv`'s `gencode`; `species.tsv`'s `supergroup`, `ncbi_lineage` and
  PR2, EUKARYOME and KSGP lineages; `genomes.tsv`'s rRNA counts; `busco/copies.tsv`'s `contig_edge`, `gaps`,
  `confidence` and `flags`; `features.tsv`'s `gap`; `rrna/records.tsv`'s `contig`, `start`, `end`, `strand`,
  `contig_edge`, `introns`, `gaps`, `confidence` and `flags`; the manifest's `inputs.references`. EUKARYOME's records
  are labelled by their species (its `s__` is the epithet alone: until 4 every one was labelled by its genus).
  Since 2026-11-06, still format 4: the reference renamed Eukore (was EUKDB), its folder `Eukore_<version>/` (was
  `EUKDB_<version>/`); the manifest's `eukore`, the version of Eukore that built it; `benchmarks/` and the manifest's
  `benchmarks`, the build's measures of the package (a reader of format 4 needs neither).
