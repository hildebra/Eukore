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
2. **The calls** (`eukore/odb.py`, `eukore/loci.py`).
   - miniprot aligns, per locus, the proteins of the 30 dataset species nearest to the genome in NCBI's taxonomy and
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
     in the genome's ORFs (esl-translate, ORFs of 20 residues or more; domains of E-value 1e-3 or less); its best hit
     and the hits on that contig and strand within 20 kb of it, not overlapping on the HMM, place it and, joined in the
     HMM's order, make a "self" protein (nearly the genome's own). miniprot aligns the locus' queries and the self
     protein on that place and 5 kb each side again, with sensitive seeding (`-k5 -M0 -n1 -l4`, introns as long as the
     first pass allows), and the HMMs confirm as before. In the held-out *B. hominis*, the BUSCOs missed went from 32 to
     12 (all on the right gene). A rescued copy is `rescued` and at most `medium`: its span may fall short of the gene's
     ends (median overlap 0.65 with the call made with its own proteins), and its exon ends come from distant proteins.
     A genome over 1 Gb is not rescued (its ORFs' search too long; such genomes are mostly animals' and plants', near
     the dataset's species).
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
  ITS1/5.8S/ITS2/28S borders, searching the profile sets of the genome's lineage.
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
  the SSU) and COX1 CDS. Of a species' records, the one with the most of these genes, then its own taxid's, then the
  longest.
- **Else its own assembly's organelle molecules**, those NCBI's FASTA titles name a mitochondrion, chloroplast, plastid
  or apicoplast. nhmmer finds each gene, queried with that gene of the 3 nearest RefSeq species (its exons joined);
  barrnap the rRNA nhmmer misses (`--kingdom mito` in animals' mitochondria, `bac` elsewhere). This runs after the
  genomes in the process stage, and again in the package stage before the prune, on the genomes still on disk.

A gene is kept whole, introns and all (fungi's and plants' COX1 hold group I and II introns). Checked on real genomes
(2026-11-05): human's mitochondrion (queried with the chimpanzee's) gave its 12S, 16S and COX1 within 6 bases of
RefSeq's annotation; yeast's (queried with *Kluyveromyces lactis*') its COX1 exactly (13,818-26,701, 8 exons), its 15S
and 21S within 90 bases. barrnap alone found no rRNA in yeast's mitochondrion, and Arabidopsis' plastid 16S and 23S
but not its 5S. A species with neither source has no organelle genes. The organelle genes have no checks across species
yet.

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
