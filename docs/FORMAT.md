# Eukore package format, version 6

Eukore is a standalone eukaryote reference. It holds every BUSCO gene of eukaryota_odb12 and the rDNA operons of every
eukaryote species with an NCBI genome, complete (genes with all their introns, operons uncut) and checked, and the
species' mitochondrial and plastid rRNA, mitochondrial COX1 and CYTB, and plastid rbcL, matK and tufA where NCBI has
them, and what the
build measured of each genome before deleting it (`metadata/`). Its taxonomy is NCBI's, tied to no
GTDB release; each species also carries its lineage in PR2, EUKARYOME and KSGP (`species.tsv`), the versions named in
the manifest. Each copy and record carries its checks, a confidence tier and the reasons for it, and the hints a user
needs to spot a chimera (contig ends, assembly gaps, the screens). A user of the package (protal's build) takes what
it needs (protal: the 30 BUSCOs of its euk30, chosen by `busco/loci.tsv`'s figures) and converts it: naming, length
limits, masking, record cuts. Paths below are relative to the package folder, `Eukore_<YYYY.MM>/`. A
`.gz` file is BGZF (blocked gzip, read by any gzip reader; since 2026-11-09). A
`.tsv` has a header line, tab-separated columns and no quoting. `-` means "none". Coordinates are 0-based, end
excluded, on the sequence as the package holds it, unless they are a contig's (1-based, inclusive).

Since format 5 (Eukore 0.2) it also holds the copies' coding sequences and proteins, their introns and the genomes'
genetic codes and duplication, more organelle genes with their checks, the 18S variable regions and in-silico PCR of
barcode primers, a two-level supergroup, three more references' lineages and other databases' ids, and a barcode gap
measured on a panel of species ([Changes](#changes)).

Since format 6 (Eukore 0.2.3, 2026-11-09) a confidence tier tells a sequence nothing could be checked against
(`unverified`) from one where something may be wrong (`medium`); the coding sequences are checked against the proteins
their HMMs confirmed; a Duplicated locus' copy is the one whose lineage fits; a sequence two clades' genomes share goes
to its owner; an 18S is checked to be a nucleus'; and each genome names what made its calls (`genome_record`).

### Confidence tiers (format 6)

Every copy, record and organelle gene has one of four tiers, and its reasons in `flags` (comma-separated, those of the
worse tiers first):
- `high`: a check across species confirms it (its nearest copies of other species agree with its lineage to the
  family or deeper), and nothing says something may be wrong;
- `unverified`: nothing says something may be wrong, but its own-species check could not decide
  (`own_species_undecided`: no other species of its family near enough, often a protist lineage the package holds
  few species of; for a published read, no genome of its family; until 2026-11-10 a published read's was not made);
- `medium`: something may be wrong (a fragment, a duplicate, a short span, an assembly gap, a contig's end, a stop
  inside its protein, a sequence shared with another clade's genome whose owner is undecided, ...);
- `low`: a check calls it another organism's, or not the gene (a contaminant, a bacterium's or an organelle's SSU as
  an 18S, a marker's paralog).

Until format 6 there were three: `unverified` was `medium`, `own_species_undecided` one of its reasons. A user who
drops `low` and weighs `medium` weighs `unverified` as it sees fit: it is a lack of evidence, not evidence.

## Every file

| file | what it holds | format |
|---|---|---|
| [`manifest.json`](#manifestjson) | versions, inputs, counts, every file's size and MD5 | |
| [`taxonomy/`](#taxonomy) | NCBI's taxonomy cut to the package; nuclear, mitochondrial and plastid genetic codes | codes: 5 |
| [`species.tsv`](#speciestsv) | one row per species: its genome, supergroup and group, NCBI lineage (and filled), lineages in PR2, EUKARYOME, KSGP, SILVA, UNITE, MIDORI2, type material | 4, more in 5 |
| [`crossrefs.tsv`](#crossrefstsv-format-5) | each species' ids in Open Tree, GBIF, WoRMS, IRMNG, Index Fungorum, UNITE (SH), BOLD (BIN) | 5 |
| [`reference_screen.tsv`](#reference_screentsv-optional) | the records of the MIDORI2 files the build was given that are another organism's gene, a chimera, another gene, or doubtful (optional) | |
| [`genomes.tsv`](#genomestsv) | each genome: accession, flags, rRNA copies, genetic code and its check, duplication, introns, type material, genome record, copies another clade's | 4, more in 5 and 6 |
| [`copy_numbers.tsv.gz`](#copy_numberstsvgz-format-5) | each genome's copies of every BUSCO and marker locus, rRNA gene and unit, organelle molecule and gene; by reads where measured | 5 |
| [`busco/loci.tsv`](#busco) | each BUSCO's call rate, overall and per group, and over the genomes alone | 4 |
| [`busco/copies/`, `cds/`, `proteins/`](#busco) | each BUSCO's copies whole (introns and all), their coding sequences, their proteins | cds, proteins: 5 |
| [`busco/copies.tsv`](#busco) | each copy's call, checks, confidence, internal stops, duplicated copies' identity, coding sequence check, how it was chosen | 4, more in 5 and 6 |
| [`busco/features.tsv`](#busco) | each copy's exons, introns, frameshifts, masked runs, gaps | 4 |
| [`busco/introns.tsv`, `intron_sites.tsv`](#busco) | each intron's phase, splice sites, class, U12 look, place on the HMM; the sites shared across species | 5 |
| [`markers/`](#markers-format-5) | Eukore's own marker loci (RPB1, RPB2, EF-2, EF-1α, ...): the same files as `busco/`, for genes BUSCO cannot hold | 5 |
| [`rrna/records.tsv`, `records.fna.gz`](#rrna) | 18S spans and rDNA units of genomes and published long reads, their parts, checks, confidence, SSU model | 4, more in 6 |
| [`rrna/regions.tsv`](#rrna) | each record's 18S V1-V5, V7-V9 | 5 |
| [`rrna/amplicons.tsv`](#rrna) | in-silico PCR of 10 18S, ITS and 28S primer pairs on each record | 5 |
| [`organelle/genes.tsv`, `genes.fna.gz`](#organelle) | mitochondrial and plastid rRNA and protein-coding genes, the trnH-psbA spacer, their checks and confidence | 4, more in 5 and 6 |
| [`organelle/cds.fna.gz`, `proteins.faa.gz`](#organelle) | the protein-coding genes' coding sequences and proteins by the organelle's code | 5 |
| [`organelle/amplicons.tsv`](#organelle) | in-silico PCR of 8 COI, rbcL, matK, psbA-trnH, 12S and 16S primer pairs | 5 |
| [`organelle/molecules.fna.gz`](#organelle) | the assemblies' organelle molecules, whole | 4 |
| [`intraspecific/`](#intraspecific-format-5) | the panel's further assemblies' copies, rRNA and organelle genes; the barcode gap | 5, more in 6 |
| [`metadata/`](#metadata) | each genome's measures (size, contiguity, GC, masking, telomeres, ...), tetranucleotides, sketches, NCBI's summary, codon evidence | codons: 5 |
| [`benchmarks/`](#benchmarks-optional) | the build's measures of the package | |

## manifest.json
`format` (6; a reader of 6 also takes 5 and 4), `version` (`YYYY.MM`), `eukore` (the version of Eukore that built it, e.g. `0.2.4`), `built` (ISO time),
`inputs` (BUSCO's dataset and its date, the taxdump's MD5,
EUKARYOME's and ROD's versions, NCBI's assembly tables' MD5s, the `groups` built, the no-assembly records' MD5, RefSeq's
organelle release files' MD5s (`refseq_organelle`), NCBI's FCS-GX and ANI reports' MD5s (`ncbi_quality`), the 18S
references species.tsv's lineages come from (`references`: PR2's, EUKARYOME's, KSGP's version, file, MD5 and the
species and genera it names; since 2026-11-10 MIDORI2's further `files` and its `screen`: the screen's version, the
records screened and flagged, those called `misnamed`, `chimera`, `other_gene` and `weak`, and those `dropped`,
`cut` and for `review`), Eukore's marker set (`markers`, since 2026-11-07: `data/markers.tsv`'s, KOfam's
`ko_list`'s and the members table's MD5s, Swiss-Prot's release, NCBI's query scope and proteins fetched per marker, when
it was built, and per marker its kind, members, cutoff and members dropped as foreign; `null` without one), ENA's list
of TSA transcriptome sets (`tsa`, since 2026-11-07: its file, MD5, date fetched and sets; `null` without one or with
`--no_tsa`);
`inputs.partial.species_clades` in a trial package limited to those clades, which a user
should not take for Eukore), `counts` (species, species with rRNA alone, genomes, copies and their
calls, rRNA records, organelle genes by organelle and gene, species with organelle genes, taxa, `busco_loci`, the
copies' and records' confidence tiers, `copy_unverified` (since 2026-11-10: per locus set the `unverified` copies and those
`alone_in_family`: of a species no other species of whose family, or genus where NCBI names no family, holds the
locus in the package), `labels`: per reference the species by how they matched; `metadata`: the rows of
`metadata/ncbi.tsv.gz`, the genomes, those described, sketched, the organelle molecules, and the genomes by
`genome_source` and `annotation_source`), `benchmarks` (the
headline figures of `benchmarks/summary.json`, where the build measured them), and `files`: every file of the package
with `[bytes, md5]`. The manifest is written last: a package without one is unfinished.

## taxonomy/
NCBI's taxonomy, cut to every taxon on the package's lineages: `nodes.tsv` (taxid, parent, rank, gencode: NCBI's nuclear genetic code, the table a BUSCO copy translates by; mgcode and pgcode: its mitochondrial and plastid codes, which the organelle genes translate by, since 2026-11-07; pgcode is 11 where the taxdump has no plastid column), `names.tsv` (taxid,
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
| `source` | `ncbi` (NCBI's genome of the species by Eukore's ranking), `model` (a model species' reference), `no_assembly` (records of a species NCBI has no assembly of), `toga2` (TOGA2's genome of a vertebrate species NCBI has no assembly of, accession `TOGA2_<assembly>`), `tsa` (format 5: an INSDC transcriptome of a species with no genome, its TSA set's accession; its copies flagged `transcript`), `rrna` (rRNA records alone) |
| `loci` | BUSCO genes called in its genome |
| `rrna_records` | its records in `rrna/records.tsv` |
| `organelle_genes` | its records in `organelle/genes.tsv` |
| `supergroup` | its deep lineage (format 5; `data/supergroups.tsv`, its most specific NCBI clade there): `Obazoa`, `Amoebozoa`, `Archaeplastida`, `TSAR`, `Haptista`, `Cryptista`, `Discoba`, `Metamonada`, `CRuMs`, `Ancyromonadida`, `Malawimonadida`, `Hemimastigophora`, `Provora`, else `incertae sedis`. Until format 4 the finer groups below (SAR one of them, `Haptista` all of it) |
| `group` | the finer level (format 5) `busco/loci.tsv`'s call rates and the benchmarks split species by: `Vertebrata`, `Metazoa` (the other animals), `Fungi`, `other Opisthokonta`, `Apusozoa and Breviatea`, `Amoebozoa`, `Viridiplantae`, `Rhodophyta`, `Glaucophyta`, `other Archaeplastida`, `Stramenopiles`, `Alveolata`, `Rhizaria`, `Telonemia`, `other TSAR`, `Haptophyta`, `Centroplasthelida`, `other Haptista`, `Cryptophyceae`, `other Cryptista`, `Discoba`, `Metamonada`, `other eukaryotes` |
| `ncbi_lineage` | NCBI's names at domain;kingdom;phylum;class;order;family;genus;species (a rank NCBI does not name empty) |
| `ncbi_lineage_filled` | the same, an empty rank filled (format 5): with the next unranked NCBI clade below the rank above it (a protist's `Sar`, `Stramenopiles`, ...), shallowest first; where none is left, the name above with `_X` (PR2's way: `Nematoda_X`, `Nematoda_XX`) |
| `pr2`, `pr2_match` | its lineage in PR2 (domain;supergroup;division;subdivision;class;order;family;genus;species, as PR2 writes it: `Saccharomyces_cerevisiae`), and how it was found (below) |
| `eukaryome`, `eukaryome_match` | its lineage in EUKARYOME (superkingdom;kingdom;subkingdom;phylum;subphylum;class;order;family;genus;species: its classification table to the genus, its long reads' species) |
| `ksgp`, `ksgp_match` | its lineage in KSGP (`k__` .. `s__`, as KSGP writes it; its eukaryote nuclear records) |
| `silva`, `silva_match` | its lineage in SILVA's SSU Ref NR99 (138.2; format 5): its eukaryote records' path, of any depth, then the species (`Saccharomyces_cerevisiae`; the organism's first two words where they name one) |
| `unite`, `unite_match` | its lineage in UNITE's general release of all eukaryotes (10.0, 2025-02-19; `k__` .. `s__`; format 5); its species hypotheses are in `crossrefs.tsv` |
| `midori2`, `midori2_match` | its lineage in a MIDORI2 release (root;domain;phylum;class;order;family;genus;species; a rank MIDORI2 does not give empty), where the build was given its COI file or several of its genes' (`--midori2`; format 5): MIDORI2 names each record's NCBI taxid, so a species is found by it (`taxid`) before its name. Since 2026-11-10 a record the reference screen drops ([`reference_screen.tsv`](#reference_screentsv-optional): another organism's gene, another gene) names no species |
| `type_material` | NCBI's relation of its genome to the species' type material (`assembly from type material`, `assembly from synonym type material`, `assembly designated as neotype`, ...; `genomes.tsv`'s), `-` for none (format 5) |

A species is found in a reference by name: `species` (NCBI's scientific name), `synonym` (one of NCBI's synonyms, an
anamorph), `genus` (neither, but the reference names its genus: the lineage down to the genus, no species), `taxid` (MIDORI2's
own NCBI taxid of the species; format 5), `-` (none:
the lineage `-`). Names are compared without case, `_` as a space and a `(Fungi)` after a genus. A reference's records
of an organelle (PR2's `:plas`, `:mito`; KSGP's chloroplast, `plas` and `mito` names) and its names that name no species
(`sp.`, `Genus_X`, accessions, `?`) are not taken. A name a reference files under several lineages (a homonym of
another kingdom) takes the lineage that shares the most names with NCBI's lineage of the species, then the one of the
most records. The labels are a name match, not a placement: a user that names species by its own rules (protal's KSGP
rules) can compare its name with these. A row's taxid may be above species (a published read labelled by its genus).

## crossrefs.tsv (format 5)
Each species' ids in other taxonomies and databases, one row per id: `taxid` (the species'), `reference`, `id`, `name`
(the name the reference gives it), `match` (how it was found):
- `OTT`: the Open Tree of Life's taxon (Open Tree Taxonomy 3.7.3), found by its NCBI taxid (`taxid`): no name matched.
- `GBIF`, `WoRMS`, `IRMNG`, `IndexFungorum`: the ids Open Tree gives the same taxon (`ott`): GBIF's backbone key, ...
- `UNITE`: UNITE's species hypotheses (SH codes) of the species' name or a synonym (`name`).
- `BOLD`: BOLD's BINs of the species' name or a synonym (`name`), where the build was given a BOLD data package
  (`--bold_tsv`; BOLD's bulk downloads need a login).

## reference_screen.tsv (optional)
Since 2026-11-10, where the build was given MIDORI2 files (`--midori2`) and has minimap2 and vsearch: the records of
those files the reference screen calls (`eukore/refscreen.py`; [methods.md](methods.md#the-reference-screen-since-2026-11-10)),
one row per record and gene file. A record with `action` `drop` names no species in `species.tsv`; the others do. A
reader takes a package with or without the table. `tools/screen_reference.py` writes the same calls as KSGP's lists
(`drop.tsv`, `trim.tsv`, `review.tsv`) and the files without the dropped records, for a user of MIDORI2 itself.

| column | meaning |
|---|---|
| `reference` | `MIDORI2` |
| `record` | the record's name in its file (`MG559732.1.<1.>690`: accession, version, span) |
| `gene` | its file's gene, by the organelle genes' names (`COX1`, `CYTB`, `SSU`, `LSU`, ...; `organelle/genes.tsv`'s) |
| `taxid`, `species_taxid` | the NCBI taxid MIDORI2 gives it; its species |
| `length` | its bases |
| `call` | `misnamed` (another clade's gene under its name: a contaminant, a misidentified specimen), `chimera` (a piece of another clade's gene in it), `other_gene` (another gene, whole or a piece), `weak` (doubtful: one of two records is misnamed and which is not told, or a piece that is not certain) |
| `action` | `drop` (left out), `trim` (cut to `keep_start`-`keep_end`: its own part), `review` (kept as it is) |
| `keep_start`, `keep_end` | the part kept of a cut record (1-based, inclusive), else `-` |
| `clade` | the clade its sequence is another's: NCBI's first rank of phylum, class and order at which it parts from the record's (`order Lepidoptera`); the other gene of an `other_gene` |
| `label_clade` | the record's own clade at that rank |
| `identity` | the identity (%) of its nearest record of that clade (or gene), over a window of 300 bases |
| `relatives`, `own_clade` | over the windows that clade claims: the mean identity of that clade's relatives (its second-best species under 97%), and of the record's own clade's nearest record; `-` for an `other_gene` |
| `windows` | its windows (300 bases every 50) by what they are: `foreign`, `undecided`, `own`, `gene`, `none` |
| `layout` | the stretches of windows in the record's order (`own:0-600 foreign:450-900`: 0-based starts and ends) |
| `targets` | the records nearest its foreign stretches (RefSeq's as `<accession>_<gene>`, else MIDORI2's names) |
| `reason` | the call in words |

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
| `source` | as `species.tsv`'s; `assembly` for a further assembly of the species (no copies), `intraspecific` for one the barcode gap's panel ran (its copies in `intraspecific/`; format 5) |
| `representative` | `t` for the species' genome |
| `loci` | BUSCO genes called in it |
| `rrna_18s`, `rdna_units`, `rdna_complete`, `rdna_contigs` | its 18S spans, rDNA units (from the 18S to the 28S), units with a whole 28S, and the contigs holding units, as `rrna/` found them in the assembly: copies in the assembly, a lower bound of the copies in the genome (an rDNA array is usually collapsed in an assembly); `-` where its rRNA was not searched |
| `rdna_copies_reads`, `mt_copies_reads`, `pt_copies_reads` | optional (a package `tools/rdna_depth.py` has measured): rDNA units, mitochondrial and plastid genomes per haploid nuclear genome by its reads' depth (`depth.tsv`, rounded); `-` where not measured |
| `local` | the genome's file in the build's kept genomes (`genomes/<accession>.fna.gz` beside the package), `-` where it is not kept |
| `gencode` | NCBI's nuclear genetic code for its taxid, the table its BUSCO copies were called and are translated by (format 5) |
| `duplicated_share`, `dup_identity_median`, `duplication` | of its BUSCOs called (Single or Duplicated) the share Duplicated; the median `dup_identity` of those; `haplotypes` (15% or more Duplicated, their copies 99.5% alike or more: an assembly holding both haplotypes of a diploid, or a recent hybrid), `haplotigs` (since 2026-11-10: 15% or more Duplicated, 95% to under 99.5% alike, and a copy of them on a contig under a quarter of the Single copies' median contig far more often than a Single copy lies on one: an assembly not purged of the second haplotype of an outbred diploid; `wgd_or_polyploid` until then), `wgd_or_polyploid` (15% or more Duplicated, less alike than 99.5% and not `haplotigs`: a whole-genome duplication, a polyploid, an older hybrid, or two haplotypes of an outbred diploid assembled equally long, which identity and contig lengths do not tell from a recent polyploid), `low` (under 15%), `-` (under 20 BUSCOs called, or a genome run before 0.2). A heuristic: read it with the genome's size and the species' biology (format 5) |
| `introns`, `introns_per_kb_cds`, `median_intron_length` | its BUSCO copies' introns (`busco/introns.tsv`): their number, per kb of the copies' coding sequence, their median length (format 5) |
| `genome_record`, `foreign_copies`, `foreign_source` | what made its calls (Eukore's genome record, `genomes.GENOME_RECORD`: 1 for a run before 2026-11-09, 2 since; a build runs a genome of an older one again unless `--reversion` keeps its calls); its copies another clade's (`own_species` f); the NCBI taxid of those clades' lowest common ancestor, `-` (format 6) |
| `gencode_check`, `gencode_codons`, `best_table` | its code checked against its BUSCO genes ([methods.md](methods.md)) on the codons nuclear codes are known to read otherwise (the stops TAA, TAG, TGA and CTG): `agrees`, `differs:<codon>><amino acid>,...` (NCBI's table reads them otherwise than the genes' conserved positions do; `<stop>>*` since 2026-11-11: a stop the table reads as an amino acid at 1% or less of that amino acid's other codons inside the copies, 200 or more: unused, a stop), `undecided` (too few codons at conserved positions); of those four, the ones seen at 10 conserved positions or more (a stop at 5: under its own reading it is at none; a stop not seen is a stop; a CTG not seen is not judged: CTG-clade yeasts avoid it), `-` for none; the NCBI table that fits them best (NCBI's own where it fits as well, a stop called unused read as a stop; `-` where undecided). Every codon's evidence is `metadata/codons.tsv.gz` (format 5) |

A species' genome is NCBI's best by Eukore's ranking: never one NCBI flags as contaminated, of an unverified source
organism, misassembled, chimeric or of low quality, nor one NCBI's FCS-GX screen finds 50% or more another organism's
(`FCS-GX: mostly foreign`) or that fails NCBI's ANI species check (`ANI: failed`; Saccharomycotina); another flagged
one, or a metagenome's, only without better (`FCS-GX: contaminated`: 1% or more contaminant, or 100 kb or more of
another eukaryote's sequence; `ANI: inconclusive`); then a RefSeq reference, then a representative genome; then the
better assembly level, one FCS-GX finds no contaminant in, the fewer contigs per megabase, the newer. FCS-GX's summary
is the fetched assembly's (its RefSeq copy where it has one).

## copy_numbers.tsv.gz (format 5)
Since 2026-11-08: every species' genome's copy numbers in one long table, a row per genome and element, for a reader
to pivot into a genome-by-element matrix. The same counts are in the files each element's own folder holds; this
table puts them side by side, with the zeros (a locus searched and not found).

| column | |
|---|---|
| `accession`, `species_taxid` | the genome (`genomes.tsv`'s; a species' genome only, not EukDetect2's records) |
| `set` | `busco`, `markers`, `rrna`, `mitochondrion`, `plastid` |
| `element` | a BUSCO or marker locus (each one the genome was searched for: a genome run before a locus was has no row of it); `rDNA_unit`, `18S`, `5.8S`, `28S`, `5S`; `molecules` (the assembly's sequences its FASTA titles name that organelle), an organelle gene (`organelle/genes.tsv`'s `gene`, of source `assembly`: found on the genome's own molecules); `genome` (by reads) |
| `copies` | basis `assembly`: a locus' complete copies (`copies.tsv`'s `copies`; 0 for a Fragmented or Missing one), the rRNA's (rDNA units; 18S spans alone or in units, and a unit's 18S shorter than an 18S span, at a contig's end, as partial; the units' 5.8S, 28S and inline 5S: 5S arrays apart from the units are not extracted), the organelle molecules, an organelle gene's copies on them (a plastid's inverted repeat: 2); basis `reads`: copies per haploid nuclear genome (`depth.tsv`'s `rdna_copies_reads`, `rdna_18s_copies`, `mt_copies_reads`, `pt_copies_reads`) |
| `partial` | a Fragmented locus' fragment (1), rDNA units and 28S without a whole 28S, a partial organelle gene; `-` by reads |
| `basis` | `assembly`; `transcripts` for a transcriptome (TSA: a locus' copies are its isoforms, not genes); `reads` (rows `tools/rdna_depth.py table --into_package` adds and replaces) |

Copies in an assembly are a lower bound of the genome's: an rDNA array collapses to a few units, haplotypes may be
collapsed or kept (`genomes.tsv`'s `duplication`), and a family marker's copies (actin, H3, ...) are those nearly as
alike to its queries as the best (miniprot's alignments of 95% of a query's best score, 10 of each): its gene
family's near members, not all of them. The manifest's `counts.copy_numbers` gives its rows per set.

## busco/
Every BUSCO gene of BUSCO's eukaryota_odb12 (129 at its 2026-05-22 release) called in every species' genome.
- `loci.tsv`: per BUSCO how widely the package's species hold it, for a user to choose its loci by:

| column | |
|---|---|
| `locus` | the BUSCO (an OrthoDB eukaryota id, `<group>at2759`) |
| `name` | OrthoDB's name of it |
| `cutoff`, `hmm_length` | BUSCO's score cutoff, the length of its HMM |
| `species` | species with a genome in the package searched for the locus (since 0.2.2: a genome run before a marker was added counts for it nowhere, not as one without it) |
| `called`, `duplicated`, `fragmented` | species whose genome holds it complete (Single or Duplicated); of them Duplicated; species with a fragment alone |
| `call_rate` | `called` over `species`: every record, transcriptomes too |
| `genomes`, `genome_call_rate` | since 2026-11-10: the species searched whose record is a genome assembly (not a transcriptome, `tsa`, nor EukDetect2's records, `no_assembly`), and the share of them holding the locus complete; `-` for none |
| `transcriptomes`, `transcriptome_call_rate` | since 2026-11-10: the species searched whose record is a transcriptome, and their share: a BUSCO a transcriptome misses is a gene not expressed, not one lost |
| `genome_call_<group>` | since 2026-11-10: `genome_call_rate` within each group that has such a species |
| `call_<group>` | the call rate within each group of the package's species (`species.tsv`'s `group`: Vertebrata, Metazoa, Fungi, Viridiplantae, Stramenopiles, Alveolata, Rhizaria, ...; until format 4 `call_<supergroup>` of its eleven supergroups, SAR one of them) |
| `median_length`, `median_cds_length`, `median_protein_length` | of its complete copies |

- `copies/<locus>.fna.gz`: one record per genome, `>accession`, the copy whole: the genome from the gene's first codon
  to its last codon, all introns, on the gene's strand, upper case, nothing masked.
- `cds/<locus>.fna.gz` (format 5): the same records' coding sequences: the exons joined, the bases of miniprot's
  frameshifts left out (`features.tsv`'s `frameshift`), the stop codon kept.
- `proteins/<locus>.faa.gz` (format 5): their proteins, translated by the genome's table (`genomes.tsv`'s `gencode`),
  the stop that ends them left out; an internal stop is `*` (`copies.tsv`'s `internal_stops`).
- `copies.tsv`: per copy how it was called and checked.

| column | |
|---|---|
| `accession`, `locus` | the copy |
| `status` | `Single` (one complete copy in the genome), `Duplicated` (several: the best is the copy), `Fragmented` (its HMM alignment covers less than 80% of the profile); the shipped model's: since 2026-11-10 a locus whose annotated or recovered model is complete is `Single` where miniprot's own was a fragment |
| `contig`, `contig_length`, `start`, `end`, `strand` | where it lies in the assembly (1-based, inclusive) |
| `length`, `cds_length`, `exons` | of the copy; of its coding sequence; its exons |
| `identity`, `coverage` | miniprot's identity and the share of the query protein aligned (of miniprot's own model also where `model` names an annotation: the annotated model is not aligned to the query) |
| `hmm_score`, `cutoff`, `hmm_coverage` | its translation's score against the locus' BUSCO HMM, BUSCO's cutoff, the share of the profile aligned |
| `copies` | complete copies of the locus in the genome |
| `query` | the dataset protein that placed it (`<busco>_<taxid>_<id>`, or an ancestral variant; `<busco>_self` for the rescue pass' self protein) |
| `rescued` | `t` where the rescue pass called it (its HMM in the genome's six-frame ORFs placed it; [methods.md](methods.md)), `f`, or `-` for a genome run before 2026-11-06 |
| `model` | where its exons come from: the genome's own annotation, its translation confirmed by the HMM (`ncbi`, `veupathdb`, `ensembl` or `toga2`; [methods.md](methods.md)), `miniprot`, or `-` for a genome run before 2026-11-06 |
| `protein_length` | of its translation |
| `window_class` | tiara's class of the 50 kb around it (`eukarya`, `bacteria`, `archaea`, `organelle`, `unknown`), `-` where not classified |
| `fcs` | NCBI FCS-GX's call over the copy (`ACTION:division:top taxon`, its action `EXCLUDE`, `TRIM`, `FIX` or `REVIEW`), `-` for none or where not screened |
| `own_species` | `t`, `f` or `-`: its coding sequence's nearest copies of other species agree with its NCBI lineage to the family or deeper; or, at 90% or more, are known to part from it at the class or above (a phylum or class both name and differ in, else two groups; since format 6 two lineages that name no such rank to compare are not), its owner told (below): another clade's gene, a contaminant, a misidentified genome |
| `same_as` | the NCBI taxid of the lowest common ancestor of the species with a copy 99% identical to it, `-` for none |
| `length_check` | `short`, `long` (its protein under 0.75 or over 1.33 times the median of its relatives', its 3 or more hits of other species; else under 0.5 or over 2 times its locus' median; until format 6 its locus' bottom or top 5%, a tenth of every locus whatever its lengths), `ok` |
| `contig_edge` | bases from the copy to its contig's nearer end: a copy at an end may be cut, or joined to another sequence |
| `gaps` | runs of 10 Ns or more in the copy (assembly gaps: a misjoin's mark) |
| `confidence`, `flags` | `high`, `unverified`, `medium` or `low` ([the tiers](#confidence-tiers-format-6); `unverified` since format 6), and why (comma-separated, `-` for none; below) |
| `internal_stops` | stop codons inside its protein by the genome's table (a pseudogene, a frameshift miniprot did not see, or a genetic code NCBI's table does not give: `genomes.tsv`'s `gencode_check`; format 5) |
| `dup_identity` | a Duplicated copy's nearest other complete copy in the genome: the identity of their coding sequences (vsearch, %), `<50` where none is 50% identical; `-` for a Single or Fragmented copy, or a genome run before 0.2 (format 5) |
| `cds_check` | its coding sequence (as `cds/` holds it) against the protein its HMM confirmed: `clean`; `recovered` (the best HMM score's alignment was not clean, another alignment of the gene, or one made again on its place, is and was confirmed: it is the copy); `stops` (a stop inside its translation); `mismatch` (its translation is not that protein: under 90% of its 5-mers in it); `-` for a genome of genome record 1 (format 6; [methods.md](methods.md)) |
| `cross` | the species of another clade whose genome holds this copy's sequence (90% or more): for a copy another clade's (`own_species` f), its source; for an undecided one, the other holder; else `-` (format 6) |
| `choice` | how the copy was chosen among a Duplicated locus' complete copies: `hmm` (the best HMM score's, as before format 6), `lineage` (another, whose own-species verdict fits its lineage where the best one's did not; format 6) |

A copy is `low` where a check calls it another organism's or not the locus: `own_species` (f), `fcs` (FCS-GX's EXCLUDE,
TRIM or FIX), `window_prokaryote` (tiara's bacteria or archaea), `paralog` (format 6: a marker copy whose best profile
is a decoy's, EFL called as EF-1a by a genome run before that decoy). It is `medium` where something may be wrong:
`fragmented`, `duplicated` (since 2026-11-10 not in a genome whose `duplication` is `haplotypes`, `haplotigs` or
`wgd_or_polyploid`: two copies are that genome's truth, its `status` and `copies` still say so), `rescued` (the rescue
pass': its span may fall short of the gene's ends), `length_short` or
`length_long`, `gap`, `contig_edge` (under 1,000 bases), `window_organelle`, `fcs_review`, `transcript`, and since
format 6 `internal_stops`, `cds_mismatch` and `cross_shared` (its sequence held by another clade's genome too, and
neither clade's other species tell whose it is). It is `unverified` where none of these holds but `own_species` is `-`
(`own_species_undecided`; `medium` until format 6). Else `high`.

A sequence two clades' genomes hold goes to its owner (format 6; until then both copies were `f`, the true owner's
too): the clade whose other species hold it nearer; where neither holds it near, not the clade that holds the locus as
another sequence; else undecided (`cross`, `cross_shared`). A genome with 3 or more copies another clade's of one
phylum has its undecided copies shared with that phylum called `f` too (`genomes.tsv`'s `foreign_copies`).

- `features.tsv`: per copy its exons and introns, the assembly's soft-masked runs of 20 bases or more
  (`soft_masked`), its low-complexity runs (`low_complexity`: DUST's score over 20 in 64-base windows) and its
  assembly gaps (`gap`: runs of 10 Ns or more).

| column | |
|---|---|
| `accession`, `locus` | the copy |
| `kind` | `exon`, `intron`, `soft_masked`, `low_complexity`, `gap`, `frameshift` (bases inside an exon that miniprot's frameshift leaves out of the coding sequence; since 2026-11-06) |
| `start`, `end` | on the copy |

A copy failing a check is kept and flagged, so that a user chooses.

- `introns.tsv` (format 5): every intron of the copies, with its place in the gene as its locus' HMM match state
  (each protein aligned to BUSCO's HMM by hmmalign), so that introns compare across species.

| column | |
|---|---|
| `accession`, `locus`, `intron` | the copy, the intron's number in it (from the gene's start) |
| `start`, `end`, `length` | on the copy, as `features.tsv`'s intron |
| `phase` | the coding bases before it, mod 3: 0 between codons, 1 or 2 inside one |
| `donor`, `acceptor`, `class` | its first and last two bases; `GT-AG` (the major spliceosome's), `GC-AG`, `AT-AC` (most minor-spliceosome introns), `other` |
| `u12` | `t` where it looks like a minor-spliceosome (U12) intron: its 5' end `[AG]TATCCTT` and a branch point `TCCTT[AG]A[CT]` 8 to 40 bases from its 3' end, one difference allowed in each; else `f` |
| `aa_position` | the protein residue it falls in or before (0-based) |
| `hmm_position` | that residue's match state on the locus' HMM, `-` where the residue is in an insert |
| `site` | its site of the locus (`intron_sites.tsv`): `<hmm_position>.<phase>` of the site's most common place; introns of one phase a match state apart are one site |

- `intron_sites.tsv` (format 5): per locus and site, how widely its genomes hold an intron there: `locus`, `site`,
  `hmm_position`, `phase`, `genomes` (holding one), `share` (of the genomes whose copy was aligned), `share_<group>`
  (within each group, as `loci.tsv`'s `call_<group>`).

## markers/ (format 5)
Eukore's own marker loci (`data/markers.tsv`, `eukore/markers.py`): classic phylogenetic and barcode genes that are not
BUSCOs, because OrthoDB groups each with its paralogs or lineages duplicated it ([loci.md](loci.md)). Two kinds:
- `clean`, an ortholog apart from its old paralogs: `RPB1` and `RPB2` (RNA polymerase II's two largest subunits),
  `EF2` (EF-2) and `TEF1` (EF-1α);
- `family`, a multigene family whose best copy is taken (its `Duplicated` is its nature, no `confidence` flag):
  `TUB2` (β-tubulin), `ACT` (actin), `HSP90` (cytosolic HSP90), `GAPDH` (glycolytic GAPDH) and `H3` (histone H3).

Each is defined by KEGG KOfam profiles of the ortholog and of its paralogs (decoys: RNA polymerases I, III, IV and V's
subunits for RPB1 and RPB2; the U5 snRNP GTPase, mitochondrial EF-G and RIA1 for EF-2; eRF3, HBS1 and EF-Tu for EF-1α;
α, γ, δ and ε tubulin for β-tubulin; the actin-related proteins ARP1 to ARP10 and the actin-likes for actin; GRP94 and
TRAP1 for HSP90, which also take the plastid's and mitochondrion's HSP90; the Calvin cycle's GAPA/B, the spermatogenic
and the NAD(P) GAPDHs for GAPDH; CENP-A for H3): a region whose best profile is a decoy's is not the marker. A copy's
`hmm_coverage` is read over its profile's span share (`hmm_spans`: KOfam's profiles are longer than the proteins they
model), so that `complete` means for a marker what it means for a BUSCO. Their copies are called as the BUSCOs' are,
in genomes (introns and all) and transcriptomes, and the folder holds the same files as `busco/`: `loci.tsv` (its
`name` the marker's, `cutoff` the set's), `copies/`, `cds/`, `proteins/<locus>.*.gz`, `copies.tsv`, `features.tsv`,
`introns.tsv`, `intron_sites.tsv`. `genomes.tsv`'s `loci` counts both sets, its `duplication` the BUSCOs alone; the
manifest's `counts.marker_loci` the markers.

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
| `regions` | `name:start-end,...` in order on the record: `18S`, `ITS1`, `5.8S`, `ITS2`, `28S`, and `5S` where a 5S lies in the same block (the unit as one stretch, on either strand: some protists' and fungi's repeats hold it; a tandem array's spacer 5S is one unit's, the one whose 28S it follows; of a unit's several the one after its 28S); a region not found is absent. A published read's record is the read on its unit's strand (reverse-complemented where its unit lies on -), its regions on it |
| `length` | of the record |
| `lsu_complete` | `t` where barrnap's 28S hit covers its whole model |
| `ambiguous`, `homopolymer` | the share of non-ACGT bases; the longest one-base run |
| `type` | its 99% 18S type within the species (`1`, `2`, ...) |
| `own_species` | `t`, `f`, or `-`: whether its nearest records of other sources agree with its NCBI lineage to the family or deeper (a genome's contaminant rDNA fails; a genome's copies searched as one per 99% type, so that its own do not fill the hits; its own assembly's records in another source, ROD's, are its own: format 6), `f` by the 18S' own rule (since 2026-11-10; until then `busco/copies.tsv`'s): 97% or more like species all of another class or phylum, its own clade's nearest record (the deepest of its family, order and class another record is of) 2 points or more under that, and the other clade's other species holding it nearer (or that species alone in its clade and its own clade holding nothing 90% like it); a lineage with no record of its class beside it only at 99.5% or more to another phylum's species whose clade holds it too. Of each species the best hit counts; the clades that tell the owner are searched alone where the first 32 hits held none; an 18S with introns is compared as its exons and as its span. A published read's (EUKARYOME's, ROD's; since 2026-11-10, until then `-`) against the genomes' records alone: its source's other reads, labelled by the same hand, vouch for nothing |
| `screen` | the part screen's flags: `-`, or `part:rank,...`, each a part (`18S`, `junction`, `28S 5'`, `28S 3'`) whose nearest records of other genomes are another clade's at that rank (`kingdom` .. `family`) |
| `label_from` | how its taxid was found: `assembly`, `rod_taxid`, `eukaryome_name`, `eukaryome_genus` |
| `contig`, `start`, `end`, `strand` | a genome's record: where it lies in the assembly (1-based, inclusive; the 18S span's, or the unit's from its first region to its last); `-` for a published read |
| `contig_edge` | a genome's record: bases to its contig's nearer end |
| `introns` | the 18S's introns on the record (`18S:start-end,...`: the gaps between barrnap's 18S pieces), `-` for none; the 18S region spans them |
| `gaps` | runs of 10 Ns or more in the record |
| `confidence`, `flags` | `high`, `unverified`, `medium` or `low`, and why (below) |
| `cross` | the species of another clade whose record holds its 18S: an `f`'s source; an undecided record's other holder where its own clade's nearest 18S is 2 points or more under it or the two are 99.5% alike (format 6; the 18S' own rule since 2026-11-10) |
| `ssu_model` | the SSU model of Rfam that scores its 18S best (cmsearch of the profile HMMs of RF01960 and RF02542, eukaryote, RF00177, bacteria, RF01959, archaea; a prokaryote model by 20 bits over the eukaryote's best): `eukarya`, `bacteria` or `archaea` (a bacterium's, a mitochondrion's or a plastid's SSU: `low`), `-` where none scored it or not measured (format 6) |

A record is `low` where a check calls it another organism's or of bad sequence: `own_species` (f), `screen`,
`ambiguous` (over 0.5%), `homopolymer` (36 or longer), `ssu_prokaryote` (format 6: its 18S scored best by Rfam's
bacterial or archaeal SSU model, a bacterium's or an organelle's SSU). It is `medium` where something may be wrong:
`gap`, `contig_edge` (under 1,000 bases), `lsu_partial` (an operon without a whole 28S), `genus_label` (a published read
labelled by its genus alone), `cross_shared` (format 6). It is `unverified` where none of these holds but its
own-species check is undecided (`own_species_undecided`; `medium` until format 6; until 2026-11-10 every published
read's, its check not made). Else `high`.

An operon's spacers are its ITS1 and ITS2 (`regions`); the intergenic spacer between units is not taken (a unit runs
from its 18S to its 28S, and its 5S where that lies in the same block). A 5S elsewhere in the genome (most animals' and
plants' 5S arrays) is not searched.

A record that fails a check is kept and flagged, so that a user chooses.

- `regions.tsv` (format 5): each record's 18S variable regions V1-V5 and V7-V9, as Hadziavdic et al. 2014 place them on
  *S. cerevisiae*'s 18S (Z75578; `data/rrna_regions.tsv`), carried to every record by Infernal's cmalign through Rfam's
  eukaryote SSU model (RF01960): `record`, `region`, `start`, `end` (on the record, 0-based, end excluded), `length`,
  `complete` (`t` where the record's 18S reaches both of the region's ends; a partial 18S holds a region in part). A
  record without an 18S has none. V6 is conserved in eukaryotes and is not given; the 28S's D domains are not given (no
  published coordinates of them on one reference were found; the D1-D2 barcode is the `NL1_NL4` amplicon).
- `amplicons.tsv` (format 5): in-silico PCR of the eukaryote barcode and amplicon-survey primers
  (`data/primers.tsv`, pairs in `data/primer_pairs.tsv`: 18S V4 `V4_TAReuk`, V9 `V9_1391F_EukB` and
  `V9_1380F_1510R`, `NS1_NS4`; ITS `ITS1F_ITS4`, `ITS3_ITS4`, `gITS7_ITS4`, `ITS86F_ITS4`; 28S `NL1_NL4`, `LR0R_LR5`)
  on every record:

| column | |
|---|---|
| `record`, `pair` | the record, the primer pair |
| `forward_start`, `forward_mismatches`, `forward_3prime` | where the forward primer matches best on the record (0-based), its mismatches (IUPAC; an N in the record matches nothing), those in its 3'-most 2 bases |
| `reverse_start`, `reverse_mismatches`, `reverse_3prime` | the reverse primer's, after the forward one (on the record's strand: its reverse complement's start) |
| `start`, `end`, `length` | the product, primers included (0-based, end excluded) |
| `amplifies` | `t` where each primer has 3 mismatches or fewer, none in its 3'-most 2 bases, and the product's length is in the pair's range; else `f`. The counts are given for a user's own rule |

## organelle/
Each species' mitochondrial rRNA, plastid rRNA, mitochondrial COX1 and (since 2026-11-07) CYTB (vertebrates' and
insects' barcode), and plastid rbcL and matK (land plants' core barcode; rbcL also algae's) and tufA (green algae's
barcode), and since format 5 mitochondrial COX3, nad1, nad5 and atp6 and plastid psbA, ycf1 and the trnH-psbA spacer
(land plants' barcode), where NCBI has them. Of each organelle, a
species' genes come from one source: its RefSeq organelle genome (the one with the most of these genes, then of the
species' own taxid, then the longest), else the organelle molecules of its own assembly (those NCBI's FASTA titles
name a mitochondrion, chloroplast, plastid or apicoplast; since 2026-11-10 also a kinetoplast, maxicircle,
minicircle, chromatophore or cyanelle).
- `genes.fna.gz`: one record per gene, `>record`, the gene from its first base to its last on its strand, introns and
  all, upper case. A trans-spliced gene (a plant mitochondrion's nad1, nad5: exons on both strands or out of order;
  `trans_spliced` t) is its exons joined. The trnH-psbA record runs from trnH-GUG's first base, across the spacer, to
  200 bases into psbA's 3' end, so that the barcode's primers sit on it (`regions`).
- `cds.fna.gz`, `proteins.faa.gz` (format 5): the protein-coding genes' coding sequences (their exons joined) and
  proteins (translated by the species' organelle code: `taxonomy/nodes.tsv`'s `mgcode`, `pgcode`; since 2026-11-10 a
  RefSeq gene by its own `/transl_table` and from its `/codon_start`; the stop that ends
  them left out), where the exons are known: RefSeq's, miniprot's, or a gene of one exon. nhmmer's hits joined are not
  a gene's exons (`flags` `exons_approximate`).
- `genes.tsv`:

| column | |
|---|---|
| `record` | its name: `<seqid>_<gene>` |
| `gene` | `SSU` (12S in animals' mitochondria, 9S in kinetoplastids' maxicircles, 16S in plastids, 18S in plants' mitochondria, ...), `LSU` (16S in animals' mitochondria, 12S in kinetoplastids' since 2026-11-10, 23S, 26S, ...), `5S`, `COX1`, `COX3`, `CYTB`, `nad1`, `nad5`, `atp6` (mitochondrion), `rbcL`, `matK`, `tufA`, `psbA`, `ycf1`, `trnH-psbA` (plastid) |
| `organelle` | `mitochondrion` (a kinetoplast too) or `plastid` (chloroplasts, apicoplasts, chromatophores, cyanelles and other plastids) |
| `taxid`, `species_taxid` | NCBI's |
| `source` | `refseq` (RefSeq's organelle genome, its annotated features) or `assembly` (the species' assembly: each gene by nhmmer with that gene of the 3 nearest RefSeq species, its exons joined; rRNA nhmmer misses by barrnap, `--kingdom mito` in animals' mitochondria, `bac` elsewhere) |
| `seqid` | the organelle sequence (RefSeq's accession, or the assembly's molecule) |
| `accession` | its assembly (`genomes.tsv`), `-` for RefSeq's |
| `start`, `end`, `strand` | on `seqid` (1-based, inclusive); `start` after `end` where the gene crosses a circular genome's origin |
| `length` | of the record |
| `exons` | RefSeq's feature's parts (a COX1's exons), or the nhmmer hits joined; 1 for barrnap's |
| `copies` | the gene's copies on the source (a plastid's inverted repeat: 2; since 2026-11-10 a RefSeq feature marked `/pseudo` is none): the record is a whole one before a partial one, then the longest, then the first (until 2026-11-10 the first), or barrnap's longest whole hit |
| `partial` | `t` where RefSeq marks the feature partial, barrnap's hit is partial, or nhmmer's hits cover under 80% of the query's exons |
| `ambiguous` | the share of non-ACGT bases |
| `name` | its common name (format 5): an animal mitochondrion's `12S` and `16S`, a land plant mitochondrion's `18S` and `26S`, a plastid's `16S` and `23S`, else the gene |
| `model` | where its exons come from (format 5): `refseq`, `miniprot` (the 3 nearest RefSeq species' proteins aligned by the organelle's code, where nhmmer joined several hits (a gene with introns) or, since 2026-11-10, its single hit is not the gene's coding sequence: it covers under 90% of its query, or holds more internal stops than miniprot's model; `--organelle_models`), `nhmmer`, `barrnap` |
| `parts` | its exons on `seqid` in the order transcribed, `start-end:strand,...` (format 5), `-` where not known (nhmmer's hits joined) |
| `trans_spliced` | `t` for a trans-spliced gene (its record its exons joined), else `f` (format 5) |
| `regions` | the trnH-psbA record's parts on it, `psbA3':0-200,spacer:200-...,trnH:...` (0-based, end excluded), else `-` (format 5) |
| `internal_stops` | stops inside a protein-coding gene's protein by the organelle's code: a NUMT, a pseudogene, RNA editing, or a code NCBI does not give (format 5) |
| `own_species`, `same_as`, `length_check` | as `busco/copies.tsv`'s, each gene of each organelle compared across species (its coding sequence where it has one, else its record; format 5) |
| `confidence`, `flags` | `low` (`own_species` f, `ambiguous` over 0.005), `medium` (`internal_stops`: RNA editing or a code NCBI does not give leaves stops in a working gene, `low` until format 6; `partial`, `length_short` or `length_long`, `trans_spliced`, `exons_approximate`, `cross_shared`), `unverified` (`own_species_undecided`: format 6, `medium` until then), else `high` (format 5) |
| `cross` | the species of another clade whose organelle holds its sequence (as `busco/copies.tsv`'s; format 6) |

- `amplicons.tsv` (format 5): as `rrna/amplicons.tsv`, the organelle genes' primer pairs on their genes: COX1 `COI_Folmer` (LCO1490/HCO2198), `COI_Leray` (mlCOIintF/jgHCO2198); rbcL `rbcLa`; matK `matK_KIM`; `psbA-trnH`; an animal mitochondrion's 12S `MiFish-U`, `12S_V5` and 16S `16S_Palumbi`.
- `molecules.fna.gz` (since 2026-11-07): the organelle molecules of each species' genome as its assembly holds them
  (those its FASTA titles name a mitochondrion, chloroplast, plastid, apicoplast, or since 2026-11-10 a kinetoplast,
  maxicircle, minicircle, chromatophore or cyanelle; up to 15 Mb each), whole,
  `>accession|seqid organelle`, upper case: for organelle genes a later build or a user wants (the nad genes, a
  plastid's psbA), whatever source `genes.tsv` took. RefSeq's organelle genomes are not copied.

## metadata/
What the build measured of each species' genome (the assembly file it read: `genomes.tsv`'s `fetched`) while the file
was on disk, and NCBI's own summaries of the genomes (since 2026-11-07). A build given `--keep_gb` deletes most
genomes after: these stand for them. A genome done by a build before 2026-11-07 is described where its file was still
there; one without a row was not.

- `assembly.tsv.gz`: one row per species' genome (`representative` `t` of `genomes.tsv`): where it and its gene models
  came from, and what `eukore/metadata.py` measured in one pass over its FASTA (`-` where it was not described):

| column | |
|---|---|
| `accession`, `fetched` | the genome (`genomes.tsv`) and the assembly file read (its RefSeq copy where it has one) |
| `genome_source` | the database it was downloaded from: `NCBI RefSeq`, `NCBI GenBank`, `TOGA2 (<who assembled it>)` (TOGA2's own copy of an assembly outside NCBI: DNA Zoo, Bat1K, ...), `EukProt` (a species without an assembly: its records, data/no_assembly.fna.gz) |
| `genome_url` | the file downloaded, `-` for EukProt's records |
| `annotation_source` | the database of the gene models its BUSCO calls took (`busco/copies.tsv`'s `model`): `NCBI RefSeq`, `NCBI GenBank` (a submitter's), `VEuPathDB`, `Ensembl`, `TOGA2`; `-` for none (miniprot's models alone) |
| `md5`, `bytes` | of the file read, as NCBI served it (gzipped) |
| `sequences`, `length`, `longest` | its sequences, their bases, the longest |
| `acgt`, `n_bases`, `other_bases` | its A, C, G and T, its N, and its other IUPAC bases (either case) |
| `n50`, `l50`, `n90`, `l90` | of its sequences (scaffolds): the length the longest reach 50% (90%) of its bases at, and how many |
| `contigs`, `contig_n50`, `contig_l50` | its sequences split at its gaps |
| `gaps`, `gap_bases` | its runs of 10 Ns or more (NCBI's gaps) and their bases |
| `short_1kb`, `short_10kb` | the share of its bases in sequences under 1 kb, under 10 kb |
| `gc` | G and C over A, C, G and T |
| `gc_windows`, `gc_window_sd`, `gc_window_q05`, `gc_window_q95`, `gc_outlier_windows` | the GC of its 10 kb windows (each sequence's, from its start; its last one where it holds 5 kb of ACGT or more): their number, spread, 5% and 95% quantiles, and the share 0.10 or more from their median (another organism's contigs, isochores) |
| `masked` | the soft-masked (lower-case) share of its A, C, G and T: 0 where its submitter did not mask it, not where it holds no repeats |
| `cpg_oe` | CpG observed over expected (CG dinucleotides × ACGT / (C × G)): low in genomes with methylated CpG (vertebrates, plants) |
| `duplicate_sequences`, `duplicate_bases` | sequences identical to an earlier one (as written: the reverse complement is not compared), their bases |
| `telomere_motif`, `telomeric_ends`, `telomere_both` | the telomeric repeat most of its sequences' ends hold, of TTAGGG, TTTAGGG, TTTTAGGG, TTAGG, TTAGGC, TTGGGG and TTTTGGGG (an end: the outer 5 kb of a sequence of 10 kb or more; holding a motif: its copies on either strand cover 500 bases or more, the motif covering most; yeasts' irregular TG repeats are none of them), `-` for none; the ends holding it; the sequences holding it at both ends (whole chromosomes, telomere to telomere) |
| `mitochondrion_molecules`, `mitochondrion_length`, `mitochondrion_gc`; `plastid_...` | its organelle molecules (`organelle/molecules.fna.gz`): their number, bases, GC |
| `sketch` | `t` where `sketches.zip` holds its sketch |

- `tetranucleotides.tsv.gz`: one row per genome of `assembly.tsv.gz`: `accession`, then the count of each of the 136
  canonical tetranucleotides (a 4-mer and its reverse complement counted as the first of the two in ACGT order; 4-mers
  with another base than ACGT left out), in either case. Divided by the row's sum: the genome's composition, as binning
  and contamination screens compare it.
- `sketches.zip` (where sourmash was there): one sourmash FracMinHash signature per genome (`signatures/<accession>.sig.gz`,
  named by its accession; k 21 and 31, scaled 10,000, no abundances), a collection sourmash reads whole: `sourmash
  search`, `compare`, `gather` and ANI against any genome or metagenome sketched alike.
- `ncbi.tsv.gz`: NCBI's assembly_summary row (GenBank's and RefSeq's tables the build read; `inputs.assembly_tables`) of
  every assembly of `genomes.tsv` and of the RefSeq copy fetched: `genome` (its `genomes.tsv` accession), `fetched`
  (`t` for the assembly file Eukore read), then every column of NCBI's tables as NCBI names them (`assembly_accession`,
  `bioproject`, `biosample`, `refseq_category`, `asm_submitter`, `seq_rel_date`, `assembly_type`,
  `relation_to_type_material`, `genome_size_ungapped`, `gc_percent`, `replicon_count`, `scaffold_count`,
  `contig_count`, `annotation_provider`, `total_gene_count`, `protein_coding_gene_count`, ...; their meanings:
  NCBI's README_assembly_summary.txt).
- `codons.tsv.gz` (format 5): the evidence of `genomes.tsv`'s `gencode_check`, per genome and codon seen at a
  conserved position of its BUSCO genes: `accession`, `codon`, `table_aa` (NCBI's table's amino acid, `*` a stop),
  `inferred_aa` (the amino acid those positions favour), `n` (the positions), `llr` (how much more likely the inferred
  amino acid is than the table's: natural log; a stop's against a uniform background).

## intraspecific/ (format 5)
The barcode gap: for a panel of well-sampled species (`data/intraspecific_panel.tsv`, 40 species), further assemblies
of each (`--intraspecific`, at most 3: other BioSamples, strains and isolates; `genomes.tsv` lists them with source
`intraspecific`, `representative` f) are run as the species' genome is, so that each marker's variation within the
species is set beside its distance to the nearest other species ([methods.md](methods.md)).
- `copies.tsv`: the further assemblies' BUSCO calls (`accession`, `species_taxid`, `locus`, `status`, `contig`,
  `start`, `end`, `strand`, `length`, `cds_length`, `exons`, `hmm_score`, `hmm_coverage`, `copies`, `dup_identity`, as
  `busco/copies.tsv`'s); `copies.fna.gz`: their copies whole, `>accession|locus`.
- `rrna.fna.gz`, `organelle.fna.gz`: their rRNA records (`>accession_kind_n`) and organelle genes (`>accession|gene`).
- `barcode_gap.tsv`: per marker (each BUSCO by its id, `18S`, `ITS`, each organelle gene; `mito_SSU`, `plas_LSU`, ...
  for the organelles' rRNA) and panel species:

| column | |
|---|---|
| `marker`, `species_taxid`, `name`, `group` | the marker and the species |
| `assemblies` | the species' assemblies holding the marker (its genome and its further ones) |
| `intra_min`, `intra_median` | the identity of their copies to each other (vsearch, `--iddef 2`, %): the least alike pair, the median pair. Since 2026-11-10 over the pairs vsearch aligns over 80% of both copies, a pair's identity the lower of its two directions; a pair it does not report is left out and counted (`pairs_left_out`), `-` where none is measured (until then such a pair read as 50) |
| `inter_max` | the identity of the species' genome's copy to its nearest copy of another species in the package (a BUSCO's: its check's, `qc.check`), `-` where none is 50% alike |
| `gap`, `resolved` | `intra_min - inter_max`; `t` where above 0 (its assemblies nearer to each other than to any other species) and its nearest other species is a congener (format 6: a gap to a farther relative says nothing of its congeners), `f` where not above 0, `-` where not measured or no congener is nearest |
| `inter_rank` | the deepest rank its nearest other species shares with it (`genus` .. `phylum`, `none`), `-` where not known (format 6). Pairs are aligned over 80% of their query at least (format 6), of both sequences since 2026-11-10 |
| `pairs_left_out` | the pairs of its assemblies' copies not measured (no alignment over 80% of both: a partial copy, or a pair under 50% alike), left out of `intra_min` and `intra_median` (last column, since 2026-11-10) |

## benchmarks/ (optional)
How complete, right and discriminating the package is, measured by the build (`eukore/benchmarks.py`; absent with
`--no_benchmarks`). docs/benchmarks.md describes each table. Since 2026-11-12 the build writes no `reads.tsv`,
`detection.tsv` or `figures/`, and `breadth.tsv` and `overlap.tsv` no other reference's: benchmarks/ makes them from a
finished package, into a folder of its own (benchmarks/README.md). A benchmark whose inputs the build lacked (BUSCO's dataset, vsearch) is left out, and
`summary.json`'s `left_out` says why. A reader needs none of them.

| file | one row per |
|---|---|
| `completeness.tsv` | species' genome: its BUSCOs Single, Duplicated, Fragmented, Missing, `complete` (Single and Duplicated over the BUSCOs searched, Missing included; species without an assembly left out), its copies' confidence tiers, its assembly's level, size and contigs, its rRNA counts |
| `orthodb.tsv` | BUSCO of a species BUSCO's dataset holds a protein of: Eukore's `status`; `recall` and `precision`, the shares of OrthoDB's protein's 5-mers in the call's translation and the reverse; `guided` (`t` where the call's query was that species' own protein); `verdict`: `concordant` (both 0.9 or more), `similar` (both 0.5 or more), `discordant`, `missed`; `held_out_status`, `held_out_query`, `held_out_recall`, `held_out_precision`, `held_out_verdict`: the same for the species' BUSCOs called again without its own proteins (`-` where the build made no such calls); `held_out_place` (`same`: on the package call's contig and overlapping it, `elsewhere`, or `-` where either is not called or the held-out calls predate their places) and `held_out_span` (the two spans' intersection over union); `held_out_rescued` (`t` where the rescue pass made the held-out call) |
| `resolution.tsv` | marker (a BUSCO; `BUSCO (all)`; `18S`, each species' best record's) and `shared_rank` (the rank a copy shares with its nearest copy of another species: `genus` .. `phylum`, `none`; `all`; `unmatched`: none 70% identical): `copies`, those `under_97` and `under_99` identical, the identity's quantiles (`q05` .. `q95`) |
| `resolution_hist.tsv` | (since 2026-11-12) `resolution.tsv`'s identities themselves: per marker, `shared_rank` and step of 0.5 points (`identity`: the step's lower bound; `-` for `unmatched`), the `copies`; benchmarks/ works the identity that parts two ranks out of it (`thresholds.tsv`) |
| `rrna.tsv` | rRNA source and kind: records, median length, records holding each region, a whole 28S, an 18S intron, each confidence tier |
| `labels.tsv` | supergroup (and `all`) and reference: its species by how they were named (`species`, `synonym`, `genus`, `none`) |
| `primers.tsv` | primer pair and supergroup (and `all`): the species with a record the pair runs on, those with one it amplifies, their share (`amplicons.tsv`; format 5) |
| `runtime.tsv` | genome timed (run since 2026-11-06): its size, bytes, BUSCOs searched, threads, seconds |
| `stages.tsv` | step of the package stage: its seconds |
| `selection.tsv` | supergroup and fate: NCBI's taxa with a GenBank assembly by what became of them (`chosen`, `chosen, flagged`, `model species`, `no assembly (records)`, `every assembly excluded`, `placeholder`, `hybrid`, `not a species`, `no download`) |
| `selection_taxa.tsv` | (since 2026-11-12) NCBI taxon with a GenBank assembly: its `taxid`, `name`, `rank`, `supergroup` (its group) and `fate` (`selection.tsv`'s), NCBI's `assemblies` of it and the `reasons`: a chosen genome's flags, or the flags that excluded an excluded taxon's assemblies, each with its assemblies where more than one (`contaminated:2`), `; ` between them. The list `selection.tsv` counts: for a reader looking a species up, the published list of the genomes left out and why, and the breadth at every rank (benchmarks/'s `breadth_ranks.tsv`) |
| `breadth.tsv` | supergroup (and `all`): NCBI's species with an assembly, Eukore's (with a genome; with rRNA alone), another reference's (`compared`: `-` in a package since 2026-11-12; benchmarks/'s comparisons fill them), and the two's overlap |
| `genes.tsv` | species' genome: the rank it shares with the nearest other species of BUSCO's dataset (`odb_nearest`), whether it is in the dataset itself (`odb_self`), its copies' exon and intron bases |
| `annotation.tsv` | copy of a genome with its own annotation (NCBI's, RefSeq's or a submitter's; else Ensembl's; else TOGA2's): the annotated transcript sharing most coding bases with miniprot's model of it (also where the copy took the annotation's model, `model` the source), both's exons and introns, those exact, the coding bases each holds and shares; `source`, the annotation's (`ncbi`, `ensembl`, `toga2`, ..; last column, since 2026-11-10: until then left out of this table, though the build had it) |
| `flags.tsv` | supergroup, confidence and flags (as `busco/copies.tsv` writes them): copies |
| `rrna_regions.tsv` | supergroup, source and region: records holding it, its length's quantiles; the organelle genes as source `mitochondrion` or `plastid` and region the gene (`organelle/genes.tsv`, those not partial; a COX1 with its introns) |
| `organelle.tsv` | supergroup and organelle gene: species, and those with the gene |
| `overlap.tsv` | supergroup (and `all`) and database (PR2, EUKARYOME's long reads and KSGP by their species names; ROD by its taxids; BUSCO's OrthoDB dataset; the other references: benchmarks/'s, since 2026-11-12): its NCBI species, those in Eukore too (any record), Eukore's alone, the database's alone |
| `database_species.tsv` | (since 2026-11-12) database of `overlap.tsv` and NCBI species it holds (`species_taxid`), in Eukore or not: the sets `overlap.tsv` counts two by two, from which benchmarks/ counts the species several databases share (`upset.tsv`) |
| `rrna_copies.tsv` | each species' genome: its assembly level and contigs, its 18S spans, rDNA units and whole units as extracted (an assembly collapses a tandem array: the copies assembled, not the genome's) |
| `reads.tsv` | (until 2026-11-12; since, benchmarks/'s) test genome (or background genome), reference (`Eukore`, `Eukore exons`, `Eukore 18S`, each `--read_compare` one), mode (`in`: its species in the reference; `held_out`: left out; `novel`: its whole order left out, since format 6), reads (`marker`: from its markers' loci; `genome`: from anywhere): reads simulated, assigned, and per rank (`species` .. `phylum`) correct and wrong (a prokaryote background read: wrong at every rank its assignment reaches; a eukaryote background read, the host's, format 6: right on its own lineage); `congeners`: the package's other species of its genus; `censored` (last column, since 2026-11-10; 0 in mode `in`): the unassigned reads whose alignment list was full (minimap2's 50 secondary alignments) of records left out in that mode, so that a hit to a kept record may have been cut off: the wrong assignments of `held_out` and `novel` are a lower bound by up to these |
| `detection.tsv` | (until 2026-11-12; since, benchmarks/'s) each test genome as a sample (its first `depth` reads from anywhere in it, 1,000, 5,000 and 20,000, and as many background reads), per reference and mode: the rank scored (`species` with its species in the reference, `genus` held out, `class` novel), whether the reference holds it (`in_reference`: the species; held out, another of its genus; novel, another order of its class), the reads placed to its right taxon (`true_reads`), `detected` (3 reads or more, and 1% or more of the sample's reads placed at that rank) and the other taxa detected (`false_taxa`; a eukaryote background's own taxon, the host truly in the sample, is none: format 6) |
| `holdout.tsv` | (format 6) per BUSCO dataset species' genome and level (`species`: its own proteins left out of the queries; `genus`, `family`, `order`: every dataset species of its genus, family, order too): `left_out` (dataset species left out), `nearest` (the rank its lineage shares with the nearest dataset species kept: `genus` .. `phylum`, `none`), the BUSCOs OrthoDB holds of it (`loci`), `called`, their verdicts against OrthoDB's gene model (`concordant`, `similar`, `discordant`, `missed`), those on the package call's gene (`same_gene`) and sharing 90% of its span (`span_90`) |
| `composition.tsv` | (since 2026-11-10) supergroup (and `all`): the package's species, those with a genome (or a transcriptome, or records: a representative) and with rRNA records alone, their genomes' median share of complete BUSCOs (`complete_median`) and the share at 90% or more (`complete_90`), the species with an rRNA record, a whole operon (`rrna/records.tsv`'s `rDNA`) and an organelle gene |
| `summary.json` | the headline figures; since 2026-11-10 also `composition` (its `all` row), `orthodb.reach` (the clade hold-out's levels pooled by the rank shared with the nearest dataset species kept: called, on the same gene, OrthoDB's gene model), `annotation.by_source` and `by_distance`, `resolution.<marker>.congeners` (the copies whose nearest other species is a congener), `reads.<reference>.censored` (per mode, the genome reads censored) |
| `summary.md` | (since 2026-11-10) the headline figures as Markdown tables, `composition.tsv` the second |
| `figures/` | (2026-11-10 to 2026-11-12; since, benchmarks/'s `eukore-figures.sh`) the figures (`eukore/figures.py`): `Fig1.png`, the overview (breadth, the calls' reach from BUSCO's dataset, the markers' resolution, read detection), then one per benchmark; `index.md` captions each and names its tables; `summary.md` as above |

## depth.tsv (optional)
Copy numbers by read depth, written into a package by `tools/rdna_depth.py table --into_package` ([read_depth.md](read_depth.md)); the manifest's `inputs.read_depth` names the tool, the date and the genomes measured. Since 2026-11-08 the tool also writes them into `copy_numbers.tsv.gz` (basis `reads`) where the package has it. One row per
representative genome it was asked for. A genome's sequencing reads (one ENA run, a sample of its reads) are mapped to
its rDNA units, its single-copy BUSCOs and its species' organelle genes; a copy number is a depth over the
single-copy BUSCOs' median depth. Both alleles of a diploid map to a single-copy gene: the numbers are per haploid
nuclear genome. The assembly's own counts (rDNA units found in it) are a lower bound: an assembly collapses an rDNA
array.

| column | |
|---|---|
| `accession`, `species_taxid` | the genome |
| `runs` | the ENA run its reads come from |
| `run_source` | how the run was found: `manual` (Eukore's `data/depth_runs.tsv`), `biosample` (the assembly's BioSample's), `bioproject` (its BioProject's, of its taxid or species), `strain` (its species' at ENA with the assembly's strain, isolate, cultivar or ecotype), `none` |
| `platform`, `preset` | the run's (`ILLUMINA`, `PACBIO_SMRT`, `OXFORD_NANOPORE`) and minimap2's preset (`sr`, `map-hifi`, `map-pb`, `map-ont`) |
| `reads`, `read_length` | reads mapped from, and the mean length of those that aligned |
| `rdna_units_ref`, `rdna_ref_source` | the rDNA units mapped to (`rrna/records.tsv`'s records, one per 18S type), from the genome itself (`genome`) or the species' published operons (`published`) |
| `sc_loci`, `sc_depth`, `sc_spread` | the single-copy BUSCOs measured (Single, confidence not low, no assembly gap), the median of their exons' mean depths (since 2026-11-10; until then of the whole copies', introns and all: a young repeat in the introns of many BUSCOs raised it, and lowered every copy number divided by it), and their interquartile range over that median |
| `rdna_depth` | the rDNA units' depths summed (each the median over its 18S, 5.8S and 28S) |
| `rdna_copies_reads` | `rdna_depth` over `sc_depth`: rDNA units per haploid genome |
| `rdna_18s_copies` | the same over the 18S alone |
| `mt_genes`, `mt_copies_reads`; `pt_genes`, `pt_copies_reads` | the species' mitochondrial (plastid) records of `organelle/` measured, and the median of their depths over `sc_depth`: mitochondrial (plastid) genomes per haploid nuclear genome, in the tissue sequenced |
| `flags` | `-`, or: `low_depth` (`sc_depth` under 3), `few_loci` (under 10), `spread` (`sc_spread` over 0.5), `uneven` (over a quarter of the BUSCOs under a fifth of their mean: reads of a few regions), `no_rdna_unit`, `no_single_copy`, `mt_at_rich`, `pt_at_rich` (the organelle's records 75% AT or more: short-read libraries cover them poorly, the number too low); `no_runs`, `not_measured`, `not_staged` for a genome without numbers; `exons_over_gene` (since 2026-11-10: the exons' depth over 1.25 times the whole copies': reads of processed pseudogenes, `sc_depth` too high) |

All values `-` where not measured. A number is one library's: libraries of one strain differed by up to 30% (README).

## genomes/ (beside the package, optional)
The eukaryote genomes the build kept, up to its `--keep_gb` where given, else all (`<accession>.fna.gz`, as NCBI serves
them). `genomes.tsv`'s
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
  `benchmarks`, the build's measures of the package (a reader of format 4 needs neither). Since 2026-11-07, still
  format 4: `organelle/genes.tsv`'s mitochondrial `CYTB` and plastid `rbcL`, `matK` and `tufA`; `organelle/molecules.fna.gz`; `metadata/` and the
  manifest's `counts.metadata` (a reader of format 4 needs none of them).
- 5 (2026-11-07, Eukore 0.2): `taxonomy/nodes.tsv`'s `mgcode` and `pgcode`; translations by every NCBI table
  (until 5 the nuclear tables 1, 6, 10, 12, 15 and 26-31 alone; another read as the standard code); `busco/cds/`,
  `busco/proteins/`, `busco/introns.tsv`, `busco/intron_sites.tsv`, `metadata/codons.tsv.gz`; `busco/copies.tsv`'s
  `internal_stops` and `dup_identity`; `genomes.tsv`'s `gencode`, `duplicated_share`, `dup_identity_median`,
  `duplication`, `introns`, `introns_per_kb_cds`, `median_intron_length`, `gencode_check` and `best_table`;
  `organelle/`'s COX3, nad1, nad5, atp6, psbA, ycf1 and trnH-psbA, trans-spliced genes (until 5 left out),
  `cds.fna.gz` and `proteins.faa.gz`, `genes.tsv`'s `name`, `model`, `parts`, `trans_spliced`, `regions`,
  `internal_stops`, `own_species`, `same_as`, `length_check`, `confidence` and `flags`;
  `rrna/regions.tsv` (the 18S variable regions), `rrna/amplicons.tsv` and `organelle/amplicons.tsv` (in-silico PCR),
  `benchmarks/primers.tsv`; `species.tsv`'s `supergroup` regrouped (the deep lineages: Obazoa, TSAR, ...), its new `group`
  (the finer level, SAR split) and `busco/loci.tsv`'s `call_<group>` in place of `call_<supergroup>`;
  `species.tsv`'s `ncbi_lineage_filled`, SILVA's, UNITE's and MIDORI2's lineages (`taxid`: MIDORI2's match by NCBI taxid) and
  `type_material`; `genomes.tsv`'s `type_material`; `crossrefs.tsv`; `intraspecific/` (the barcode gap); species of source
  `tsa` (a transcriptome; `copies.tsv`'s flag `transcript`). A reader that chose loci by format 4's `call_<supergroup>`
  (protal's euk30 choice) reads `call_<group>`.
  Since 2026-11-08, still format 5: `copy_numbers.tsv.gz` (every genome's copy numbers of its loci, rRNA and
  organelles in one table; the manifest's `counts.copy_numbers`), and `inputs.markers` and `inputs.tsa`.
  Since 2026-11-09, still format 5: every `.gz` file is BGZF (htslib's blocked gzip, at level 9; until then gzip's
  own stream, level 6). A gzip reader reads it as before; `bgzip -d -@` decompresses it on many threads and
  `samtools faidx` indexes a FASTA. The release beside the package (`OUT/release/`, eukore/release.py: the website's
  data and the downloads) is no part of the package.
- 6 (2026-11-09, Eukore 0.2.3): the confidence tier `unverified` (nothing against a sequence, its own-species check
  undecided: until 6 `medium`), in every `confidence` column and the manifest's counts, `benchmarks/completeness.tsv`'s
  and `rrna.tsv`'s `unverified`; `busco/copies.tsv`'s (and `markers/`') `cds_check`, `cross` and `choice`, and its flags
  `internal_stops`, `cds_mismatch`, `cross_shared`, `paralog`; `length_check` against the copy's relatives (until 6 its
  locus' outer 5%); `own_species` with the owner of a shared sequence told and unnamed ranks no longer another clade;
  `genomes.tsv`'s `genome_record`, `foreign_copies` and `foreign_source`; `rrna/records.tsv`'s `cross` and `ssu_model`
  and its flags `ssu_prokaryote` and `cross_shared`; `organelle/genes.tsv`'s `cross`, its `internal_stops` `medium`
  (until 6 `low`); `intraspecific/barcode_gap.tsv`'s `inter_rank`, `resolved` `t` against a congener alone. A reader
  of format 5 reads format 6 but for the tier `unverified`, which it may read as `medium`.
- 2026-11-10, no new format (optional): `reference_screen.tsv` and the manifest's `inputs.references.MIDORI2.screen`
  (and `files`), where the build was given MIDORI2 files; `species.tsv`'s `midori2` no longer names a species by a
  record the screen drops. MIDORI2's later header layout (GB273's: a tab, every NCBI rank as `rank_Name_taxid`) is
  read: until then such a file named no species (every `midori2` `-`). `rrna/records.tsv`'s `own_species` of a
  published read (EUKARYOME's, ROD's) is made, against the genomes' records (until then `-`, the read `unverified`
  at best), and an 18S relative under 90% no longer weighs for a clade. The 18S' `own_species` has its own rule (a lone
  lineage's 18S undecided, a 2-point margin to its own clade, the clades searched alone), an 18S with introns is
  compared as its exons too, an ITS1 of up to 4,500 bases is a unit's (3,000 until then), and a unit keeps the 5S
  after its 28S.
- 2026-11-06, no new format (optional, added to a built package): `depth.tsv`, `genomes.tsv`'s `rdna_copies_reads`,
  `mt_copies_reads` and `pt_copies_reads`, the manifest's `inputs.read_depth` (`tools/rdna_depth.py`). A reader
  takes a package with or without them.
