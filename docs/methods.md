# How Eukore is made

What the build does to each genome, rRNA record and organelle gene, and what the package tells its user. The tables
named here are described column by column in [FORMAT.md](FORMAT.md).

## BUSCO copies are made

Eukore calls every BUSCO of eukaryota_odb12 itself (OrthoDB 12.1, 129 BUSCOs; `--loci` restricts a trial), as
compleasm does. A user picks its loci among them: `busco/loci.tsv` gives each one's call rate over the package's
species and per supergroup (protal picks its euk30). A genome searched for fewer BUSCOs than the build calls (a run
before format 4) is run again. EukDetect1 and EukDetect2 were its sources of copies (30 loci) until 2026-11-05.
1. **The genome set.** For each NCBI species with a GenBank assembly in the `--groups`, the assembly that ranks best
   (`eukore/genomes.py`):
   - never one NCBI flags as contaminated, of an unverified source organism, misassembled, chimeric or of low
     quality;
   - never one that NCBI's own FCS-GX screen finds mostly another organism's (50% or more), or that fails NCBI's ANI
     species check (Saccharomycotina against their type strains; `eukore/quality.py`);
   - another flagged one, or a metagenome's, only without better. Eukore adds its flags: `partial genome` (NCBI's
     `genome_rep`); `small for its clade`, an assembly under a quarter of the typical genome of its genus (else
     family, else order): a part of the genome, whatever NCBI calls it; `FCS-GX: contaminated` (1% or more
     contaminant, or 100 kb or more of another eukaryote's sequence); `ANI: inconclusive`;
   - then a RefSeq reference, then a representative genome; then the better assembly level, one FCS-GX finds no
     contaminant in, fewer contigs per megabase, newer.

   On NCBI's tables and reports of 2026-11-05/06, NCBI's reports changed 249 species' genomes (189 contaminated, 28
   inconclusive ANI, 29 contaminant-free of two equal, 3 mostly foreign) and left out 24 species whose only genomes are
   another organism's (5 *Azolla* ferns' are their cyanobacterium's, *Setophaga vitellina*'s a bacterium's). 1,352
   species' genomes stay flagged contaminated and 218 inconclusive, having nothing better (`genomes.tsv`'s `flags`).

   Its RefSeq copy is fetched where it has one (RefSeq screens eukaryote genomes for contamination). A model species
   (`data/model_species.tsv`) takes its reference assembly. NCBI's placeholder species and interspecies hybrids
   ("A x B") are left out.

   A vertebrate species NCBI has no assembly of takes TOGA2's (since 2026-11-06; TOGA2's `assemblies_and_species.tsv`:
   DNA Zoo's, Bat1K's, the Hiller lab's), of its assemblies the one of the longest scaffold N50, as TOGA2 serves it
   (a UCSC .2bit, read by `eukore/twobit.py`), with TOGA2's gene models: `genomes.tsv`'s `source` `toga2`, its
   accession `TOGA2_<assembly>`. On 2026-11-06, 72 species (65 mammals; 27 of DNA Zoo's).

   A TOGA2 assembly outside NCBI also takes the place of the NCBI genome chosen for its species (a RefSeq reference
   genome too; never a model species') where it is far better: its contig N50 twice that genome's or more (NCBI's
   assembly statistics) and its total length 90% of it or more (TOGA2's `chrom.sizes`: not a part of the genome,
   however contiguous). The NCBI genome is listed among the species' other assemblies; NCBI's contamination screens
   of it and its NCBI gene models are given up (TOGA2's models come instead). On 2026-11-06, 7 species: *Rhinopithecus
   bieti*, *Hipposideros armiger*, *Neomonachus schauinslandi*, *Myocastor coypus*, *Mastomys natalensis*, *Pusa
   sibirica*, *Alouatta macconnelli* (contig N50 2.5-3,800 times NCBI's); *Lyroderma lyra*'s, 77% of the length, not.

   A build run again on genomes already done (since 2026-11-07) compares its plan with the last one
   (`WORK/genome_set.tsv`): a species whose genome changed (TOGA2's in its place; after `--refresh`, NCBI's newer
   assembly or version) gets the new one fetched and run, and the old file deleted once that one is done. A done
   genome is run again where better gene models than its run took are there (the source order below), so that its
   copies take them, not only the benchmark's comparison.
2. **The calls** (`eukore/odb.py`, `eukore/loci.py`).
   - miniprot aligns, per locus, the 30 dataset proteins of the species nearest to the genome in NCBI's taxonomy (a
     species with paralogues gives several, so ~15-30 species) and
     BUSCO's ancestral variants, with the genome's genetic code (NCBI's: ciliates and diplomonads read TAA and TAG as
     glutamine).
   - An alignment needs identity 0.2 and coverage 0.4 of its protein. The overlapping ones form one region.
   - hmmsearch scores the region's 5 best alignments; it is the locus when the locus' own HMM scores one best and at
     BUSCO's cutoff, and that alignment of the best HMM score is its gene.
   - compleasm's gates (0.4 and 0.6) lost up to half the BUSCOs of divergent genomes (of 129: *Cyanidioschyzon* 33,
     now 62; *Bodo saltans* 62, now 77; *Galdieria* 105, now 112). Until 2026-11-06 each alignment was also given
     the next one's translation (miniprot writes it first), so that the HMM checked a neighbour's protein: earlier
     calls, the trial packages' among them, are not to be trusted.
   - A locus is Single, Duplicated (the best copy taken), Fragmented (under 80% of the profile) or Missing.
   - **The rescue pass** (since 2026-11-06; `--no_rescue` leaves it out). In a genome far from every dataset species,
     miniprot finds no seed for many genes with the distant proteins: *Blastocystis hominis* called without its own
     proteins lost 27 of the 98 BUSCOs OrthoDB annotates in it (its genes compact, its introns ~30 bases), though its
     HMMs find 26 of them in its six-frame ORFs, on the right gene. So for each locus still Missing, its HMM is searched
     in the genome's ORFs (esl-translate, ORFs of 20 residues or more, split at runs of 10 or more X: an assembly
     gap's Ns, which do not end an ORF and would make it longer than hmmsearch's 100,000 residues; domains of E-value
     1e-3 or less); its best hit
     and the hits on that contig and strand within 20 kb of it, not overlapping on the HMM, place it and, joined in the
     HMM's order, make a "self" protein (nearly the genome's own). miniprot aligns the locus' queries and the self
     protein on that place and 5 kb each side again, with sensitive seeding (`-k5 -M0 -n1 -l4`, introns as long as the
     first pass allows), and the HMMs confirm as before. In the held-out *B. hominis*, the BUSCOs missed went from 32 to
     12 (all on the right gene). A rescued copy is `rescued` and at most `medium`: its span may fall short of the gene's
     ends (median overlap 0.65 with the call made with its own proteins), and its exon ends come from distant proteins.
     A genome over 1 Gb is not rescued (its ORFs' search too long; such genomes are mostly animals' and plants', near
     the dataset's species).
   - **Frameshifts** (since 2026-11-06). miniprot's GFF keeps a frameshift's bases inside its CDS, so that the exons
     joined read out of frame after it (17 of 69 held-out *B. hominis* copies had an internal stop). The copy's coding
     sequence leaves them out (the alignment's CIGAR places them; `frameshift` features of the copy): its held-out
     copies concordant with OrthoDB's went from 21 to 24, its own-protein ones from 52 to 55.
   - **The genome's own gene models** (since 2026-11-06). Where the genome's assembly (that accession and version) has
     an annotation, a copy takes the annotated transcript sharing the most coding bases with miniprot's model where
     the locus' HMM confirms its translation (its best HMM, at the cutoff) at 95% of miniprot's score or more
     (`model` in copies.tsv: its source); else miniprot's model stays. The annotation is the first there of:
     - NCBI's (about 37% of fungal and protist species, 11-29% of the other groups', RefSeq's or a submitter's);
     - VEuPathDB's, with a VEuPathDB user's API key (`--veupathdb_key`; its organism table: 817 genomes of protists,
       fungi, vectors and hosts, curated): 71 species NCBI annotates none of. Its own sequence names (`Pf3D7_01_v3`)
       are matched to the genome's by length, where one sequence alone has it;
     - Ensembl's (its species.json, 2026-10: ~5,900 assemblies; Ensembl's own pipeline, BRAKER, Helixer, VEuPathDB's
       and community annotations, the newest build): 2,037 species of NCBI's assemblies NCBI annotates none of,
       mostly invertebrates (1,551; Darwin Tree of Life's, ERGA's), 230 plants, 196 vertebrates, 58 fungi, 2
       protists;
     - TOGA2's (Malovichko et al. 2026; 2,297 vertebrate assemblies, genes projected from a reference species'):
       907 vertebrate species more.
     Ensembl names a chromosome by its number, TOGA2 an accession without its version: NCBI's assembly report gives
     the genome's names.
     Against NCBI's curated models miniprot's are exact (every intron) for 81% of *Drosophila*'s and *Arabidopsis*'
     BUSCOs, 38% of *C. elegans*' and 46% of *B. hominis*'; every curated model passed its HMM. Re-aligning doubtful
     models (other miniprot settings, GeneWise, exonerate on the gene's window) did not do better: it gained 2 exact
     models of 67 in held-out *B. hominis* and lost up to 43 of 129 in *Arabidopsis*.
3. **The copy.** The gene from first to last codon, introns and all, uncut and unmasked. Its features are its exons
   and introns, the assembly's soft-masked runs and its low-complexity runs (DUST); the user masks as it sees fit.
   - The 50 kb around each copy are kept for contamination screens.
   - tiara classifies them where installed (`--tiara`).
   - NCBI's FCS-GX screens them with `--fcs_gx_db` (its ~470 GB database; a node of ~512 GB). A copy overlapping a
     range FCS-GX calls another organism's gets an `fcs` column.
4. **The checks across species** (`eukore/qc.py`), each a column of `busco/copies.tsv`. They cover what EukDetect
   checked with its alien index, its 99% collapse and its length filter:
   - `own_species`: the copy's nearest copies of other species agree with its lineage, or are another clade's
     (a contaminant, a misidentified genome);
   - `same_as`: the common ancestor of species with 99% identical copies;
   - `length_check`: its protein in the bottom or top 5% of its locus.
5. **Species without an assembly.** 75 species NCBI has no assembly of (EukDetect2 took them from EukProt, and
   re-assembled *Dientamoeba fragilis*' transcriptome) come from `data/no_assembly.fna.gz`. That file is the one
   remnant of EukDetect2 (CC BY 4.0; Shih, Zhao, Pollard & Lind 2026), made by `tools/no_assembly_records.py`. Their
   loci are called in those records as in a genome.

The 30 loci themselves are EukDetect2's choice among the BUSCOs it kept after it dropped those found in gut bacteria
(see [loci.md](loci.md), "The loci in prokaryotes", for Eukore's own check).

## rRNA records are made

- **Genomes.** Each genome's 18S spans and whole rDNA units, from the 18S start to the 28S end (18S, ITS1, 5.8S,
  ITS2, 28S). A 5S is included only where it lies in the same block. barrnap finds the genes; ITSx sets the
  ITS1/5.8S/ITS2/28S borders, searching the profile sets of the genome's lineage. Where ITSx finds only one end
  of the 5.8S, it still reports both ITS spans, with the other border guessed (a 45-base 5.8S); barrnap's border
  is used there instead. Sequences titled as a
  mitochondrion, chloroplast, plastid or apicoplast are left out: barrnap's eukaryote models call a mitochondrial
  SSU an 18S, and that bacterial-type gene then takes bacterial reads. Organelle rRNA comes from the organelle
  pass below.
- **Published records.** EUKARYOME 2.1's long reads and ROD v1.2's operons, labelled by NCBI taxid.
- **Checks**, each a column of `rrna/records.tsv`:
  - ambiguous bases and homopolymers;
  - the 99% 18S type within the species;
  - `own_species`, against the package's other sources;
  - the part screen (`eukore/screen.py`, `screen`): a record whose 18S, junction or 28S half is another organism's,
    by the nearest records of other genomes on NCBI's lineages (a contaminant contig, a mosaic, a chimeric read).

## organelle genes are made

`eukore/organelle.py`, for each species of the package and each of its organelles, from one source:
- **RefSeq's organelle genomes** (its release's mitochondrion and plastid GenBank files, ~2.2 GB, parsed once into
  `WORK/organelle`). Their annotated rRNA features (SSU, LSU, 5S; the 16S is an animal mitochondrion's LSU, elsewhere
  the SSU), a mitochondrion's COX1 CDS and (since 2026-11-07) its CYTB, and a plastid's rbcL, matK and tufA CDS (by
  `/gene`, else `/product`; rbcS, the small subunit, is not rbcL, nor the plastid's cytochrome b6 or b559 CYTB). Of a species' records, the one with the most of these genes, then its own taxid's, then the
  longest.
- **Else its own assembly's organelle molecules**, those NCBI's FASTA titles name a mitochondrion, chloroplast, plastid
  or apicoplast. nhmmer finds each gene, queried with that gene of the 3 nearest RefSeq species (its exons joined);
  barrnap the rRNA nhmmer misses (`--kingdom mito` in animals' mitochondria, `bac` elsewhere). This runs after the
  genomes in the process stage, and again in the package stage before the prune, on the genomes still on disk. An
  assembly searched before a gene of its organelles was added (its `WORK/organelle/assembly/` file's `#searched`
  line, none before 2026-11-07) is searched again where its genome is still on disk.

A gene is kept whole, introns and all (fungi's and plants' COX1 hold group I and II introns). Checked on real genomes
(2026-11-05): human's mitochondrion (queried with the chimpanzee's) gave its 12S, 16S and COX1 within 6 bases of
RefSeq's annotation; yeast's (queried with *Kluyveromyces lactis*') its COX1 exactly (13,818-26,701, 8 exons), its 15S
and 21S within 90 bases. barrnap alone found no rRNA in yeast's mitochondrion, and Arabidopsis' plastid 16S and 23S
but not its 5S. A species with neither source has no organelle genes. The organelle genes have no checks across species
yet.

## assembly metadata is measured

`eukore/metadata.py`, for each species' genome, while its file is on disk: at the end of its run (`genomes.run`), and
for a genome done before 2026-11-07 (or whose measure failed) in a pass after the process stage's runs and again in
the package stage before the prune. One pass over the gzipped FASTA (decompressed by igzip or pigz where there;
~27 s a Gb on one core, 2026-11-07, WSL) measures its contiguity, gaps, composition, GC windows, soft-masking, CpG,
duplicate sequences, telomeric repeats at its sequences' ends and its 136 canonical tetranucleotides, and keeps its
organelle molecules whole; sourmash sketches it alongside (k 21 and 31, scaled 10,000; `--no_sketches` none). Into
`WORK/done/<accession>/` (`stats.json`, `organelle.fna.gz`, `sketch.sig.gz`); the package gathers them into
`metadata/` with NCBI's assembly_summary rows of every genome ([FORMAT.md](FORMAT.md)). Each measure was checked
against slow counts and in any chunking (`tests/test_metadata.py`).

## What else a package tells its user

So that a user (protal's build) has what it needs without computing it again:
- **Names.** `species.tsv` gives each species' NCBI lineage and supergroup, and its lineage in PR2 (5.1.1), EUKARYOME
  (2.1) and KSGP (`--ksgp_tax`, else the 4.0 KSGP's site serves), with how it matched: its name, a synonym, its genus,
  or none (`eukore/labels.py`). NCBI's taxids merged since are in `taxonomy/merged.tsv`.
- **Gene structure.** The BUSCO copies' exons and introns (`busco/features.tsv`), the 18S's introns on each rRNA record
  and its regions (18S, ITS1, 5.8S, ITS2, 28S, an inline 5S), the organelle genes' exon counts.
- **Copies in the genome.** Per BUSCO copy its complete copies in the genome (`copies`); per genome its 18S spans and
  rDNA units in the assembly, those with a whole 28S and the contigs holding them (`genomes.tsv`), a lower bound where
  the array is collapsed; per organelle gene its copies on the organelle genome.
- **Chimeras.** Each copy's and each genome record's distance to its contig's end and its assembly gaps, beside
  `own_species`, `same_as`, tiara's window, FCS-GX and the rRNA part screen.
- **Confidence.** `high`, `medium` or `low` per copy and per record, with the flags that set it ([FORMAT.md](FORMAT.md)): a user
  drops the `low` ones and weighs the `medium` ones as it sees fit.

## The copies' genetics (since 2026-11-07, format 5)

- **Coding sequences and proteins.** Each copy's exons joined (miniprot's frameshift bases left out) and translated by
  the genome's NCBI table (`eukore/gencode.py` holds every NCBI table; until 0.2 the nuclear ones alone):
  `busco/cds/`, `busco/proteins/`. A stop inside a protein stays `*` and is counted (`internal_stops`).
- **Places on the HMM.** hmmalign (HMMER) aligns each locus' proteins to BUSCO's HMM of it (`eukore/align.py`; cached
  under WORK/align/). A residue's match state is a position shared across species; the alignments are not shipped.
- **Introns** (`eukore/introns.py`). Each intron's phase, splice-site dinucleotides and class, a U12 look (the minor
  spliceosome's 5' end and branch point), and the match state of the residue it falls in or before. Introns of one
  locus and phase a match state apart are one site; `intron_sites.tsv` gives how widely each site is held, overall and
  per group: an intron conserved across eukaryotes, one gained in a clade, one lost.
- **Genetic code check** (`eukore/codecheck.py`, the idea of Codetta: Shulgina & Eddy 2021). Each codon of a genome's
  copies at a conserved match state (its most likely amino acid at 0.5 or more) adds the state's log emission
  probabilities to the codon's tally; the codon's amino acid is the tally's best. Of the codons nuclear codes are known
  to reassign (TAA, TAG, TGA: Q, E, Y, W, C; CTG: S, A), one NCBI's table reads otherwise, seen at 10 such positions or
  more (a stop at 5: under its own reading it is at none; on the trial no genome held more than 3) and 10 nats more
  likely, is reported (`gencode_check`: `differs:CTG>S`). The other codons' tallies are kept
  (`metadata/codons.tsv.gz`) but not judged: they lean toward the amino acids HMMs favour (on the 7-clade trial of
  2026-11-07 ATG read as Leu in 30 of 33 genomes; Met is rarely a conserved position's best). miniprot placed the
  genes with NCBI's table, through an in-frame stop where it had to, so a stop the genome reads as an amino acid is
  still seen. The table that fits the tested codons best is `best_table`. Checked on that trial: every genome agrees
  with its NCBI table; given table 12 instead, *S. cerevisiae*, *A. thaliana* and *D. melanogaster* are called
  `differs:CTG>L`, best table 1. Given table 1, *C. albicans* is not caught: it uses CTG at 3 conserved positions
  only (CTG-clade yeasts avoid the codon), so `gencode_codons` does not list CTG. The ciliate *Euplotes raikovi*'s
  transcriptome (TSA GKLH01; NCBI's table 10, TGA Cys) agrees with table 10; given table 1 it is called
  `differs:TGA>C`, best table 10 (TGA at 8 conserved positions).
- **Duplication** (`loci.dup_identity`, `package.duplication_of`). Each Duplicated BUSCO's other complete copies are
  kept by the run, and vsearch measures the nearest one's identity to the chosen copy (`dup_identity`). A genome with
  many duplicated BUSCOs whose copies are 99.5% alike holds both haplotypes of a diploid (or a recent hybrid); less
  alike, a whole-genome duplication or a polyploid (`genomes.tsv`'s `duplication`). The thresholds are a heuristic.

## Organelle genes, format 5 (since 2026-11-07)

- **More genes.** Mitochondrial COX3, nad1, nad5, atp6; plastid psbA, ycf1 and the trnH-psbA spacer (land plants'
  barcode: its record runs from trnH-GUG across the spacer to 200 bases into psbA's 3' end, where the barcode's psbA
  primer sits). The protein-coding genes are named by one table (`organelle.CDS_RULES`: a feature's /gene, else its
  /product). GenBank locations are read whole, so a trans-spliced gene (a plant mitochondrion's nad1 and nad5:
  exons on both strands, or out of order) is kept as its exons joined (until 0.2 it was left out).
- **Exons in assemblies** (`--organelle_models`, default `miniprot`). miniprot aligns the 3 nearest RefSeq species'
  proteins to the assembly's organelle molecules by the species' organelle code (`taxonomy/nodes.tsv`'s mgcode,
  pgcode). Where nhmmer joined several hits (a gene with introns, its exons otherwise unknown) and miniprot aligns the
  whole protein, miniprot's model is taken; nhmmer's single hit is kept (`organelle.with_miniprot`). The rule is
  measured (`tools/organelle_models.py`, 2026-11-07): 100 RefSeq records of each organelle, each searched as an
  assembly with its own species' records left out of the queries, against RefSeq's own exons:

| genes | nhmmer: found, exons exact, median F1 of exon bases | miniprot: the same |
|---|---|---|
| without introns (mitochondrial COX1, COX3, CYTB, atp6, nad1, nad5; plastid matK, psbA, rbcL, tufA, ycf1; 16 to 95 genes each) | all (67 of 68 ycf1); 53-96%; 1.00 | all (66 of 68 ycf1); 53-96%; 1.00 |
| with introns (1 to 36 genes each) | all (16 of 17 nad1); none; 0.27-0.60 (0.79-0.99 for tufA, ycf1, matK) | 88-100%; 0-20%; 0.955-0.997 |

  Without introns nhmmer's ends are exact as often or more (COX3 67 of 86, miniprot 46; CYTB 53 of 74, miniprot 41);
  with introns miniprot's exons are nearly right (its splice sites a few bases off: group I and II introns have no
  spliceosomal signals) where nhmmer's hits give none. A gene RefSeq annotates trans-spliced or partial is not
  measured.
- **Proteins and checks.** Each protein-coding gene's coding sequence translated by the organelle's code
  (`organelle/proteins.faa.gz`; an internal stop is a NUMT, a pseudogene, RNA editing or a code NCBI does not give);
  `qc.check` on each gene of each organelle across species (own-species, 99% collapses, length outliers); a tier
  (`organelle.confidence`).

## 18S variable regions and in-silico PCR (since 2026-11-07)

- **Regions** (`eukore/regions.py`). Hadziavdic et al. 2014 give V1-V5 and V7-V9 on *S. cerevisiae*'s 18S (Z75578).
  Infernal's cmalign aligns that 18S to Rfam's eukaryote SSU model (RF01960) once, which gives each region's model
  positions, then every record's 18S (each distinct one once, cached): a region is the record's bases from its first
  to its last model position. Checked: the yeast reference's regions come back exactly; *A. thaliana*'s V4 is 396
  bases. An 18S over 4,000 bases (a misjoined unit) is not aligned; cmalign's matrices are capped at 2 GB per thread.
  On the 7-clade trial, 1,327 distinct 18S took 36 CPU minutes (~1.6 CPU seconds each).
- **Primers** (`eukore/primers.py`). Each pair of `data/primer_pairs.tsv` on each record it applies to (rRNA records;
  an organelle gene, 12S and 16S on mitochondria): the best site of the forward primer, then of the reverse primer's
  reverse complement after it, IUPAC-aware (a record's N matches nothing). It amplifies with 3 mismatches or fewer per
  primer, none in its 3'-most 2 bases (ecoPCR's rule), and a product in the pair's length range. Every pair amplifies
  its reference record (`data/primer_references.fna`, `tests/test_primers.py`); TAReuk's V4 falls on yeast's 18S where
  Stoeck et al. 2010 place it (565-981). Choices made on the way: matK by the CBOL KIM primers (390F/1326R miss
  *Arabidopsis* by 5 mismatches); psbA-trnH tested on tobacco (*Arabidopsis*' trnH holds a one-base indel under
  trnHf_05); Folmer's COI on *Bombyx* (*Drosophila* and zebrafish mismatch LCO1490 by 4). On the trial ITS1F/ITS4
  amplifies all 20 fungi and none of the 17 stramenopiles, NS1/NS4 39 of 41 species (`benchmarks/primers.tsv`).

## Taxonomy, format 5 (since 2026-11-07)

- **Supergroups** (`data/supergroups.tsv`, `eukore/groups.py`). Two levels from NCBI clades (the most specific
  listed clade wins): the deep lineages of Burki et al. 2020 as EukProt v3 names them, and a finer group for call rates
  and benchmarks (SAR split into Stramenopiles, Alveolata, Rhizaria and Telonemia). NCBI has no Obazoa, Archaeplastida
  or Cryptista node: their clades are listed one by one. The build says where a listed taxid's name no longer matches
  its taxdump's.
- **Filled lineage.** An empty NCBI rank takes the next unranked clade below the rank above (a protist's Sar,
  Stramenopiles, ...), shallowest first, else the name above with `_X` (PR2's way).
- **References.** SILVA 138.2 (SSU Ref NR99: paths of any depth), UNITE 10.0 (all eukaryotes, with SH codes) and, given
  by the user, MIDORI2 (COI; its records carry NCBI taxids, so a species is found by taxid) join PR2, EUKARYOME and
  KSGP, by the same name match. `crossrefs.tsv`: Open Tree Taxonomy 3.7.3's taxon of each species found by its NCBI
  taxid (no name matched), and the GBIF, WoRMS, IRMNG and Index Fungorum ids Open Tree gives it; UNITE's SH codes and,
  given a BOLD data package, BOLD's BINs by name.
- **Type material.** NCBI's `relation_to_type_material` of each genome (`genomes.tsv`, `species.tsv`); the genome
  ranking does not use it.

## The barcode gap and transcriptomes (since 2026-11-07)

- **Panel** (`data/intraspecific_panel.tsv`, `eukore/intraspecific.py`). 40 well-sampled species (2 to 2,400
  assemblies each); up to `--intraspecific` (3) further assemblies of each are run like the species' genome (no tiara,
  no held-out calls), each of another BioSample, strain and isolate, unflagged first; their organelle molecules are
  searched for every gene; their genome files are then deleted. For each marker (each BUSCO's coding sequence, the
  18S, the ITS, each organelle gene) the identities within the species are set beside the species' genome's nearest
  other species (a BUSCO's from its own check; the others by vsearch against every other species in the package).
  On the 7-clade trial (2 further assemblies of each of its 6 panel species, 2026-11-07), the 18S and the ITS resolve
  *A. thaliana*, *C. albicans*, *C. elegans* and *D. melanogaster* (ITS gaps of 33 to 40 points), not *S. cerevisiae*
  (other *Saccharomyces* in the trial, the hybrid *S. pastorianus* among them, hold its 18S at 100%) nor "*Blastocystis
  hominis*" (its assemblies are of different subtypes: their 18S only 75% alike). A BUSCO's gap is `-` where no other
  species' copy is 70% alike (the trial holds no congener of most panel species). Those assemblies' FASTA files hold
  no organelle molecules, so their organelle genes are not compared.
- **Transcriptomes** (`genomes.tsa_jobs`). ENA's list of INSDC's eukaryote TSA sets names their NCBI taxids
  (NCBI's own lists do not); a species with no NCBI, model or TOGA2 genome takes its set of the most bases (not an
  environmental sample), in place of EukDetect2's records. Its BUSCOs are called as in a genome; no rRNA, tiara or
  annotation; its copies are `medium` at best (`transcript`: no introns, Duplicated its isoforms), `low` where the
  own-species check calls them another clade's (cultures' food, cross-contamination). On ENA's list of 2026-11-07,
  5,025 species without a genome have a transcriptome (341 Gb), most of them animals and plants, some 400 protists.
