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
   ("A x B") are left out; since 2026-11-15 so is a bare "<taxon> sp." of a clade above the genus ("Eukaryota sp."
   held 200 unrelated assemblies in 2026-11; a genus' "Fusarium sp." stays). Of a species with 3 assemblies or more,
   one under 30% of their median size is left out (a part, an organelle's); since 2026-11-15 the median is that of
   its unflagged assemblies where it has 3, and NCBI's reference or representative genome is never left out (nine
   oversized contig assemblies of *Cygnus atratus* had dropped its 1.12 Gb chromosome-level reference).

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
   - A locus is Single, Duplicated (the best copy taken), Fragmented (under 80% of the profile) or Missing. The
     share is the union of every domain hmmsearch reports: counting only domains of E-value 1e-3 or less was
     measured (2026-11-10) and called whole genes fragments (27 of 27 shifts judged by OrthoDB's own models).
   - **Complete against the query** (since 2026-11-15, genome record 4). A lineage whose proteins lack a stretch of the
     profile (an insertion or extension of the animals') had its whole genes Fragmented: in the first full build 1,410
     of 1,707 Fragmented calls of the dataset species held 90% or more of OrthoDB's own protein at 90% of its length
     (1,200 identical), Fragmented 13-15% of Alveolata's, Discoba's and Haptophyta's calls against under 1% of
     animals' and plants'. Each dataset protein's own coverage of its locus' profile is measured once
     (`odb.query_coverage`, `query_coverage.tsv` in the dataset's folder), and a call's coverage is read over its
     query's (at least 0.5), as the markers' is over their members' (`hmm_spans`). On the 7-clade trial
     (Eukore_2026.11-paper) 267 of 385 Fragmented calls turn Single (mostly *Blastocystis*'); of the 34 that OrthoDB
     annotates, 32 are its protein (5-mers 0.9 both ways) and 2 not; of the 4 that stay Fragmented, none is.
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
   - **Split copies** (since 2026-11-14, genome record 3). A genome skim's contigs (a few kb) cut many genes: of the
     first full build's 2,664 genomes under 50% complete, 2,333 were assemblies of mean contig under 5 kb, most with
     more Fragmented loci than Missing. A locus with no complete copy after the rescue pass is joined from a query's
     partial alignments on different contigs that tile its protein in order (10 residues of overlap at most) and reach
     their contigs' ends where the gene goes on (200 bases or less), 4 pieces at most; where the HMMs confirm their
     proteins joined, as any call, it is one copy: Single where they cover 80% of the profile, else Fragmented, its
     `model` joined, its pieces joined with 100 N between two (a `join` feature), at most `medium` (`joined`). For
     it miniprot keeps a query's alignments of 20% of its best score (a split gene's lesser piece scores under the
     95% the calls take, which are taken from these). Yeast, *Drosophila* and *B. hominis* cut into 1-3 kb contigs:
     Single 62 -> 75, 68 -> 87 and 37 -> 40 (of 108, 125 and 72 whole), every join on the gene the whole genome holds;
     without the contig-end rule 7 of 60 joins were a paralog's or another gene's pieces. The whole genomes' calls
     did not change.
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
     (a contaminant, a misidentified genome); since 2026-11-09 a sequence two clades' genomes share goes to its owner
     ([below](#checks-and-choices-format-6-since-2026-11-09));
   - `same_as`: the common ancestor of species with 99% identical copies;
   - `length_check`: its protein against its relatives' (until 2026-11-09 the bottom or top 5% of its locus).

   Since 2026-11-17 the search behind the first two compares a sequence's first 5,000 bases (`qc.CHECK_BASES`): only
   RPB1's coding sequences and the organelle genes compared as records with their introns (a fungus' mitochondrial
   LSU: 4.4 kb at the median, 43 kb the longest) are longer. vsearch aligns the whole of both sequences and leaves its
   fast aligner past some 5 kb each: a fungal mitochondrial LSU cost 39 core-seconds a query, 2.4 cut to 5,000 bases,
   with the same own-species verdict for 359 of 360 organelle genes drawn (RPB1: 11 and 0.8 core-seconds, the nearest
   copy the same for 27 of 28). The 2026-11-17 build's organelle step, 6 h 01 m of its final job's 9 h 31 m, was this
   search. `length_check` reads the whole sequences.
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
  mitochondrion, chloroplast, plastid or apicoplast (since 2026-11-10 also a kinetoplast, maxicircle, minicircle,
  chromatophore, cyanelle or nucleomorph) are left out: barrnap's eukaryote models call a mitochondrial
  SSU an 18S, and that bacterial-type gene then takes bacterial reads. Organelle rRNA comes from the organelle
  pass below. An SSU on any other contig (a contaminant's 16S, a symbiont's in a transcriptome, an organelle's on a
  contig not titled so) is told by Rfam's prokaryote models scoring it best (`ssu_prokaryote`; the checks below) and,
  since 2026-11-19, set apart in `rrna/prokaryote_ssu.*`: no species' record (until then a `low` one, which took
  every bacterial read the 2026-11-17 build's rRNA records placed in the read benchmark).
- **Published records.** EUKARYOME 2.1's long reads and ROD v1.2's operons, labelled by NCBI taxid.
- **Checks**, each a column of `rrna/records.tsv`:
  - ambiguous bases and homopolymers;
  - the 99% 18S type within the species;
  - `own_species`, against the package's other sources; a published read's against the genomes' records alone
    (since 2026-11-10: until then not made, its label taken as its source gave it);
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
  or apicoplast (since 2026-11-10 also a kinetoplast, maxicircle, minicircle, chromatophore or cyanelle). nhmmer
  finds each gene, queried with that gene of the 3 nearest RefSeq species (its exons joined);
  barrnap the rRNA nhmmer misses (`--kingdom mito` in animals' mitochondria, `bac` elsewhere). This runs after the
  genomes in the process stage, and again in the package stage before the prune, on the organelle molecules the
  metadata scan kept (`WORK/done/<accession>/organelle.fna.gz`; since 2026-11-08, else the genome where still on
  disk). An assembly searched before a gene of its organelles was added (its `WORK/organelle/assembly/` file's
  `#searched` line, none before 2026-11-07) is searched again. Since 2026-11-18 an assembly's molecule that is another
  species' RefSeq organelle genome (its sequence name a RefSeq record the package takes for that species) is not
  taken: NCBI's RefSeq assembly of a species can carry another species' RefSeq organelle genome as its own
  (*Salvelinus* sp. IW2-2015's chrMT is *S. alpinus*' `NC_000861.1`). Its genes went into the package under both
  species, one record name (`<seqid>_<gene>`) twice, and protal's mitochondrial records stopped on it.
- **Else its candidate contigs** (since 2026-11-14): of a genome of 1 Gb or less titling no molecule an organelle,
  contigs of an organelle genome's length (AT-rich ones for a mitochondrion), searched as titled molecules are. Since
  2026-11-15 a candidate gives genes only where nhmmer or miniprot find a protein-coding gene of the organelle on it
  (a nuclear contig of a plastome's length holds an rRNA-like stretch more easily than a COX1 or rbcL), barrnap's
  bacterial models are not run on candidates (a contaminant's 16S became an organelle's), a contig the
  mitochondrion's genes are on is no plastid candidate, and its genes are flagged `candidate_contig` (medium).
- **One molecule, one organelle** (since 2026-11-19): a molecule the genome titles an organelle is no candidate of
  the other, and in a species that can hold a plastid (not an animal or a fungus) a mitochondrion's candidate needs
  its COX1 or COX3: a plastid holds homologues of the mitochondrion's other protein-coding genes (petB and petD, the
  two halves of a cytochrome b; atpI of atp6; ndhA and ndhF of nad1 and nad5). Until then a genome titling its
  plastid and no mitochondrion had the plastid among its mitochondrion's candidates: *Galdieria partita*'s
  chloroplast `AP025529.1` gave the mitochondrion a partial CYTB (Eukore's search with the 3 nearest RefSeq
  mitochondria hits its petB and petD, 84,458-85,606 on the minus strand, 43% of the query) and its 16S and 23S as
  the SSU and LSU, beside its own plastid genes: one molecule under both organelles, its SSU and LSU one record name
  twice (protal's conversion stopped on 11 such records, 2026-11-19). An untitled plastome was filed as a
  mitochondrion alone, its plastid genes lost. Only searches of species that can hold a plastid change (their key's
  `one_organelle1`); a fungus' or an animal's is kept.

A gene read along its molecule spans its first base to its last, introns and all; one across a circular molecule's
origin (its parts in the order transcribed turn once) runs from its first part across the origin to its last. Until
2026-11-15 the largest space between a gene's parts was taken for the origin wherever the rest of the circle was
shorter: 129 RefSeq genes with a long intron, most of them hexacorals' nad5, shipped as their second exon, the genes
between and their first exon.

A gene is kept whole, introns and all (fungi's and plants' COX1 hold group I and II introns). Checked on real genomes
(2026-11-05): human's mitochondrion (queried with the chimpanzee's) gave its 12S, 16S and COX1 within 6 bases of
RefSeq's annotation; yeast's (queried with *Kluyveromyces lactis*') its COX1 exactly (13,818-26,701, 8 exons), its 15S
and 21S within 90 bases. barrnap alone found no rRNA in yeast's mitochondrion, and Arabidopsis' plastid 16S and 23S
but not its 5S. A species with neither source has no organelle genes. The organelle genes' checks across species
(since 2026-11-07) are [below](#organelle-genes-format-5-since-2026-11-07).

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
- **Confidence.** `high`, `unverified` (since 2026-11-09), `medium` or `low` per copy and per record, with the flags
  that set it ([FORMAT.md](FORMAT.md#confidence-tiers-format-6)): a user drops the `low` ones and weighs the `medium`
  ones (something may be wrong) and the `unverified` ones (nothing could be checked against them) as it sees fit.

## The copies' genetics (since 2026-11-07, format 5)

- **Coding sequences and proteins.** Each copy's exons joined (miniprot's frameshift bases left out) and translated by
  the genome's NCBI table (`eukore/gencode.py` holds every NCBI table; until 0.2 the nuclear ones alone):
  `busco/cds/`, `busco/proteins/`. A stop inside a protein stays `*` and is counted (`internal_stops`).
- **Places on the HMM.** hmmalign (HMMER) aligns each locus' proteins to BUSCO's HMM of it (`eukore/align.py`; each
  distinct protein's row on the HMM cached under WORK/align/). A residue's match state is a position shared across
  species; the alignments are not shipped.
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
  `differs:TGA>C`, best table 10 (TGA at 8 conserved positions). A stop NCBI's table reads as an amino acid but the
  genome does not use is at no conserved position, and neither is one used rarely: oligotrich ciliates read TAA and
  TAG as Gln at 1–4 of their 22–34 conserved Gln positions. Since 2026-11-11 such a stop is judged by every in-frame
  codon inside the copies (their last left out): where it is 1% or less of its amino acid's other codons there, and
  those 200 or more, it is called a stop (`differs:TGA>*`) and `best_table` reads it as one. Measured on 18
  ciliates of NCBI's tables 1, 4, 6, 10, 28 and 29 and the 34 genomes of the standard code of the 2026-11-10 build: a
  stop truly read as an amino acid was 7.1% or more of its amino acid's codons (*Strombidium stylifer*'s TAA), *Fabrea
  salina*'s TGA (NCBI's table 10, Cys) 0.4%: called a stop, best table 1. Copies mispredicted through an intron add
  stops (up to 41% of Cys in *Blastocystis*), so these counts can hide an unused stop but never call a used one so;
  for the same reason they are not used to call a stop read as an amino acid.
- **Duplication** (`loci.dup_identity`, `package.duplication_of`). Each Duplicated BUSCO's other complete copies are
  kept by the run, and vsearch measures the nearest one's identity to the chosen copy (`dup_identity`). A genome with
  many duplicated BUSCOs whose copies are 99.5% alike holds both haplotypes of a diploid (or a recent hybrid); less
  alike, an outbred diploid's second haplotype where the copies lie on contigs far shorter than the Single copies'
  (`haplotigs`, since 2026-11-10), else a whole-genome duplication or a polyploid, or two haplotypes assembled
  equally long (`genomes.tsv`'s `duplication`). The thresholds are a heuristic.

## Organelle genes, format 5 (since 2026-11-07)

- **More genes.** Mitochondrial COX3, nad1, nad5, atp6; plastid psbA, ycf1 and the trnH-psbA spacer (land plants'
  barcode: its record runs from trnH-GUG across the spacer to 200 bases into psbA's 3' end, where the barcode's psbA
  primer sits). The protein-coding genes are named by one table (`organelle.CDS_RULES`: a feature's /gene, else its
  /product). GenBank locations are read whole, so a trans-spliced gene (a plant mitochondrion's nad1 and nad5:
  exons on both strands, or out of order) is kept as its exons joined (until 0.2 it was left out).
- **Exons in assemblies** (`--organelle_models`, default `miniprot`). miniprot aligns the 3 nearest RefSeq species'
  proteins to the assembly's organelle molecules by the species' organelle code (`taxonomy/nodes.tsv`'s mgcode,
  pgcode). Where nhmmer joined several hits (a gene with introns, its exons otherwise unknown) and miniprot aligns the
  whole protein, miniprot's model is taken, and (since 2026-11-10) also where nhmmer's single hit is not the
  gene's coding sequence: it covers under 90% of its query (one exon of a gene with introns, reached into the
  intron), or, read as a coding sequence, it holds more internal stops by the organelle's code than miniprot's model
  (an end reached out of frame or past the stop); otherwise nhmmer's single hit is kept
  (`organelle.with_miniprot`). Measured on 190 RefSeq organelle genomes searched as assemblies (973 genes): genes
  without introns shipped with internal stops fell from 59 to 9, intron genes shipped as one hit from 24 to 14;
  capping how far an end is reached was measured too and made the intron-free control worse. The rule is
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
  KSGP, by the same name match. A species matched by its genus alone takes no lineage of another deep lineage than
  its own (a genus homonym): until 2026-11-15 only animals, fungi and plants were told apart, and the bird genus
  *Lessonia* took PR2's kelp lineage; now Stramenopiles, Alveolata, Rhizaria, Rhodophyta, Glaucophyta, Amoebozoa,
  Haptista, Cryptista, Discoba and Metamonada too, judged on every name of the species' NCBI lineage. `crossrefs.tsv`: Open Tree Taxonomy 3.7.3's taxon of each species found by its NCBI
  taxid (no name matched), and the GBIF, WoRMS, IRMNG and Index Fungorum ids Open Tree gives it; UNITE's SH codes and,
  given a BOLD data package, BOLD's BINs by name.
- **The reference screen** (since 2026-11-10): MIDORI2's records are screened before they name species; see
  [below](#the-reference-screen-since-2026-11-10).
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

## Eukore's marker loci

Since 2026-11-07 (0.2): nine classic phylogenetic and barcode genes that are not BUSCOs ([loci.md](loci.md)) are called
as loci of Eukore's own (`data/markers.tsv`, `eukore/markers.py`; the package's `markers/`).

- **Definition.** Each is one or two KEGG KOfam profiles of the ortholog (its targets) and the profiles of its
  paralogs (its decoys): RNA polymerases I, III, IV and V for RPB1 and RPB2; the U5 snRNP GTPase, mitochondrial EF-G and
  RIA1 for EF-2; eRF3, HBS1 and EF-Tu for EF-1α; α, γ, δ and ε tubulin for β-tubulin; ARP1 to ARP10 and the actin-likes
  for actin; GRP94 and TRAP1 for HSP90; the Calvin cycle's GAPA/B, the spermatogenic and the NAD(P) GAPDHs for GAPDH;
  CENP-A for H3. The decoys join the BUSCOs' HMMs: a region whose best profile is a decoy's is not the marker, the
  BUSCOs' own rule. Checked on model species' Swiss-Prot paralogs: Arabidopsis's plastid HSP90-5 and mitochondrial
  HSP90-6 score best against GRP94's profile, its chloroplast GAPA and GAPB against theirs, human's testis GAPDHS against
  its own, yeast's CSE4 against CENP-A's; only Arabidopsis's plastid glycolytic GAPCp scores as its cytosolic GAPC
  does (a duplicate of it). The skeletal muscle actin's profile is neither target nor decoy: it scores plants' and
  yeast's actins as well as the actin profiles do.
- **Query proteins.** Swiss-Prot's eukaryote proteins and NCBI's named as the marker (each row's Entrez query, with a
  length range: RefSeq's of every eukaryote, GenBank's of the protists, no environmental sample; of an organism the 3 of
  the length nearest the median). Until 2026-11-09 the scope's Entrez expression took GenBank's fungi, animals and
  plants too (Entrez's NOT joins two terms: TEF1's query held 32,708 records, every one of GenBank's 17,851 fungal;
  with the scope as meant, 3,109): the contaminants below were among them. With the scope as meant the members are
  RefSeq's mostly, and fewer species' (built 2026-11-09: RPB1 3,095 to 1,792 species, RPB2 2,432 to 1,154, EF-2 3,322
  to 1,476, EF-1α 3,789 to 861, β-tubulin 4,215 to 2,635, actin 2,559 to 1,656, HSP90 2,867 to 1,516, GAPDH 1,612 to
  1,113, H3 4,754 to 1,941). On the 7-clade trial (2026-11-09) the marker calls barely moved: of 285, 206 identical
  and 72 another exon end; no fungus' call changed (the trial's fungi have RefSeq members); in *Blastocystis* RPB1's
  four calls became clean models a little nearer *Blastocystis*' own RPB1 (each had 1-2 stops; their HMM score ~10%
  lower), and ST6's HSP90 a real HSP90 (the old call another gene). A fungus or plant with no RefSeq relative near is
  queried with farther proteins now: not measured. A name only gathers: a candidate is the marker's where its best profile of all the
  targets and decoys is one of the marker's targets and one of them reaches its KO's own threshold (KOfam's rule; a
  plant's actin scores best against the vertebrate actin profile, under its threshold, and reaches the other
  eukaryotes' one), which leaves out the misnamed and the paralogs. A member is of 0.6 to 1.6 times the median length
  and complete (its profile's coverage at 0.8 of the median's), both medians of whole records, Swiss-Prot's and
  RefSeq's: GenBank's protist records are many PCR fragments (EF-1α's: 351 residues on average against RefSeq's 460;
  16,000 of them are left out). Of a species one: the best checked record (Swiss-Prot's, then RefSeq's curated NP_,
  RefSeq's predicted XP_, GenBank's), then one of typical length (within 10% of the median: Arabidopsis's plastid
  GAPCp, a transit peptide longer, scores 1.7 bits over its cytosolic GAPC), then the best scoring.
- **Contamination.** A GenBank member (RefSeq annotates assemblies NCBI screened; Swiss-Prot is curated) is dropped
  as another organism's protein where its nearest member by shared 5-mers (both ways) is of another phylum, at 0.9 or
  more and more alike than the marker's members of two phyla are but for 1% of pairs (pairs of its kingdom where the
  nearest is of it: animals' β-tubulins are alike across phyla, fungi's are not; else pairs of two kingdoms), and no
  member of its own phylum is within 0.15 of it (that is the true owner of a contaminant's copy: *Meyerozyma* sp. JA9's
  RPB1, with *M. guilliermondii*'s at 0.88 behind *Naganishia cerealis*'s copy at 0.97); within its kingdom, only where
  its phylum has other members (a lone phylum's conserved β-tubulin or H3 is no evidence). On the build of 2026-11-07
  this drops 12 of ~32,000 members, each another organism's protein to look at: *Naganishia cerealis*'s RPB1, EF-2
  and β-tubulin (*Meyerozyma*'s), *Cerrena zonata*'s β-tubulin and H3 (*Hyphopichia*'s), *Anthurium amnicola*'s EF-2
  (a thrips's) and H3 (*Rhizophagus*'s, its mycorrhiza), an insect's actin filed as *Colletotrichum*'s, a vertebrate's
  as *Acacia*'s, β-tubulins of one fungal phylum filed under another; one doubtful (*Mnemiopsis*'s β-tubulin). The
  build's `members.tsv` gives every candidate assigned, its source, score, coverage and fate (member, foreign,
  incomplete, other).
- **Calls.** As the BUSCOs' (miniprot, then the HMMs), with the marker's cutoff half its members' 5th percentile score.
  KOfam's profiles are longer than the proteins they model (EF-1α's 772 positions for ~460 residues), so a call's HMM
  coverage is read over its profile's median coverage of the members (`hmm_spans`), and `complete` means what it
  means for a BUSCO. A family marker's call is its best copy (`Duplicated` is no flag); a genome's `duplication` is
  read on its BUSCOs alone.
- **Validation** (`tools/validate_markers.py`, 2026-11-07, on the trial's 47 genomes): each genome whose species has
  a member (17: the models, five *Saccharomyces*, two *Candida albicans*, two *Trichuris trichiura*, two
  *Blastocystis*, a second *D. melanogaster*) has its markers called as the build calls them, with every query and
  held out (its species' own proteins left out of the queries), and each call's protein compared with its species'
  member by shared 5-mers both ways. Held out, 115 of 124 calls are concordant (0.9 or more both ways), 7 similar
  (0.5), 1 discordant, 1 missed; held out changes no verdict (9 calls' scores by under 1.2 bits, an exon's end). RPB1
  and RPB2 are concordant in all 26. Every call scores over its best decoy (HSP90 against GRP94
  the nearest: 1,165 against 994 bits in *D. melanogaster*). The rest: the second *D. melanogaster* assembly
  (GCA_038396815) calls EF-1α, β-tubulin and GAPDH with lower scores than the reference assembly does (808 against
  1,056 bits; 618 against 709), which calls them all concordant: that assembly's errors; *S. pastorianus*'s GAPDH, a
  hybrid's other subgenome copy; *C. elegans*'s H3, a variant; *Blastocystis hominis*'s EF-2 and EF-1α (0.7), its
  Swiss-Prot proteins of another subtype. *Blastocystis* ST4's HSP90 is missed with its own protein as the query: the
  one region it finds in GCA_000963395 is 45% identical, a GRP94 gene, scoring 1,033 bits against GRP94's profile and
  782 against HSP90's: the decoy rule kept it out. Before the members' checks, KOfam's two actin profiles failed
  Arabidopsis's actin (its best, the vertebrates', under its threshold) and Arabidopsis's member was its plastid GAPCp:
  both fixed above.
- **EFL** (since 2026-11-09; `data/marker_decoys.tsv`, `markers.decoy_profile`). The elongation factor-like GTPase that
  replaced EF-1α in green algae, diatoms, haptophytes, chytrids and others (Keeling & Inagaki 2004) has no KEGG KO, and
  EF-1α's profile scores it best: of NCBI's 237 EFL records (2026-11-09, mostly PCR fragments) all scored best against
  TEF1's profile and 51 over TEF1's cutoff, so a genome holding EFL alone had it called as its TEF1. Its decoy is built
  of NCBI's EFL records that TEF1's target scores best, under the members' 5th percentile and unlike every member
  (shared 5-mers under half), aligned to EF-1α's profile (hmmalign) and built (hmmbuild, `decoy-EFL`). A member that
  scores the decoy over its target is the decoy's protein filed under the marker's name and is dropped
  (`members.tsv`'s `decoy`). Built on the trial's marker set (2026-11-09): 112 seeds of NCBI's 153 records; every one of
  the 237 EFL records then scores best against the decoy (none against TEF1; 26 against eRF3's KO, as before), 3,780
  of TEF1's 3,789 members against TEF1, and the other 9 were EFL: eight green algae's (*Volvox*, *Chlorella*,
  *Ostreococcus*, *Tetraselmis*, ...) and *Salpingoeca rosetta*'s "EF-1α", scoring ~830 bits against the decoy and
  470-500 against TEF1 (its cutoff 425). Those had been TEF1's queries for green algae. The package stage checks every
  marker copy again by the set's profiles and decoys (`markers.recheck`): a copy whose best profile is a decoy's, a
  genome's called before the decoy existed, is `low` (`paralog`).

## Checks and choices, format 6 (since 2026-11-09)

A critical review of the build (2026-11-09) found calls made by older code kept beside new ones, coding sequences out
of frame shipped as `high`, the true owner of a contaminant's sequence called another clade's, a tier that counted
"could not be checked" as "may be wrong", and bacterial SSU accepted as 18S. These are the changes.

- **The genome record** (`genomes.GENOME_RECORD`, `done.json`'s `genome_record`). A done genome was done for good
  whatever code made its calls; five fixes to the calls of 2026-11-06 reached none of the genomes done before them, and
  the trial package of 2026-11-07 held 32 genomes of 33 called before the frameshift fix (521 of 3,796 copies with a
  stop inside, 30-40% of each *Blastocystis* genome's; *S. pastorianus*' duplication not measured). Each run now
  stamps its genome record (2: this format's calls), and a change to what a genome's run makes that a package should
  not mix with older runs' raises it: a done genome of an older record is run again (fetched again where pruned).
  `--reversion` stamps every done genome with the current record instead, keeping its calls (the choice of a build
  whose genomes are done and good enough; `done.json`'s `reversioned` notes what it had). The package stage names
  genomes of an older record and `genomes.tsv`'s `genome_record` says each one's.
- **The coding sequence checked** (`loci.cds_state`, `copies.tsv`'s `cds_check`). The HMM confirmed miniprot's own
  translation of an alignment, the package shipped the exons joined: whatever made the two differ (each of the
  parser's past errors did) shipped unseen. Each confirmed alignment's coding sequence, as its copy would ship it, is
  now translated and set beside the protein the HMM confirmed: `clean`, `stops` or `mismatch` (under 90% of its 5-mers
  in that protein). A region's alignments are taken clean first (complete before a fragment where the best is complete,
  then by HMM score), a locus' regions too (a Duplicated locus' gene before a pseudogene's copy), and a copy still not
  clean is aligned again on its place (its span and 2 kb each side, its locus' queries, miniprot as it is, with a
  cheaper intron (`-J15`) and with the rescue pass' sensitive seeding): a clean model the HMMs confirm there, on the
  copy's contig and strand and overlapping it, takes its place (`recovered`). The rule changes no call that was clean:
  the exon-refinement trials of 2026-11-06 that hurt curated genomes chose by HMM score among every model, this one
  only replaces a broken one. Checked on four trial genomes (2026-11-09, all 129 BUSCOs, no annotation): the best HMM
  score's alignment was not clean in 3 of *B. hominis*' 121 confirmed regions, 26 of *Blastocystis* sp. RusS' 103, 5
  of *Trichuris trichiura*'s 126 and 2 of *C. elegans*' 127 (each a stop inside, none a mismatch). Of the copies
  shipped (110, 100, 121, 124), another alignment of the gene mended 1, 12, 3 and 1, aligning again 1, 10, 1 and 0;
  0, 4, 1 and 1 ship with a stop. (The trial package of 2026-11-07, its calls made before the frameshift fix, shipped
  22, 42, 28 and 10.) A copy that stays unclean ships as it is, `medium` (`internal_stops`, `cds_mismatch`);
  the organelle genes' stops are `medium` too (RNA editing leaves stops in a working plastid gene; until 2026-11-09
  `low`). An annotated model is not taken where its protein has a stop inside and miniprot's has none, nor where it is
  over 1.5 times miniprot's protein and the profile (a fusion of two genes: the HMM's score does not see the second).
- **A Duplicated locus' copy by lineage** (`package.collect`, `copies.tsv`'s `choice`). The best HMM score chose among
  a locus' complete copies before any check: in a divergent genome a contaminant's copy from a clade the HMMs know well
  scores best, and the genome's own was dropped (KOfam's GAPDH and HSP90 profiles, built with bacteria's proteins too,
  may give a bacterium's copy the same edge). A run now keeps the three best other complete copies whole (`alternates.*` in the done folder);
  the package checks them across species beside the copies, and takes the one whose own-species verdict ranks best
  (agrees, then undecided, then shared with another clade, then another clade's), the best HMM score's where none
  ranks above it (`choice` `hmm` or `lineage`); its copy, features, coding sequence and introns are written. A genome
  run before (genome record 1) has no alternates: its copies stay the best HMM score's.
- **Cross-contamination: the owner told** (`qc.verdict`, both `own_species` checks). A copy whose nearest copies (at
  90% or more) are another clade's was `f`, and so was that clade's own copy of the same sequence where its own
  relatives lay more than 1 point (the 18S: 0.5) behind the contaminated genome's: both genomes' copies `low`, the true
  owner's too. The unit test passed only because its congener lay 0.7 points behind. Now the sequence goes to the clade
  whose other species hold it nearer; where neither holds it near, not to the clade that holds the locus as another
  sequence (a species' family holds the gene, none near this copy: the copy is not the family's); else it is undecided,
  `cross_shared` (`medium`), with the other holder in `cross`. A genome with 3 or more copies another clade's of one
  phylum has its undecided copies shared with that phylum called `f` too; `genomes.tsv`'s `foreign_copies` and
  `foreign_source` sum a genome's. The exposure is large: on NCBI's tables of 2026-11-05, 1,352 species keep a genome
  NCBI calls contaminated, and the transcriptome-only species that the check compares against (some 5,000, cultures
  known for cross-contamination) are 27% of TSAR's species, 37% of Amoebozoa's, 61-63% of Haptista's and Cryptista's.
- **Unnamed ranks** (`qc.apart`). "Another clade" meant sharing no named rank down to the class; two sister protist
  genera whose phylum, class, order and family NCBI names none of shared none and were another clade's at 90%. Another
  clade is now one that a phylum or class both name tells apart, else (naming none to compare) another group
  (`groups.group_of`: Stramenopiles, Alveolata, Fungi, Metazoa, ...). The markers' member filter reads it the same way.
- **The 18S a nucleus'** (`regions.ssu_models`, `records.tsv`'s `ssu_model`). barrnap's eukaryote models alone call a
  bacterium's or a mitochondrion's SSU an 18S (methods above), the seeds' motifs are universal, and tiara and FCS-GX
  screen the BUSCO windows alone: a contaminant contig's 16S, a NUMT or an endosymbiont's insertion gave a species an
  "18S" (`medium`, nothing being near it). cmsearch scores each distinct 18S by the profile HMMs of Rfam's eukaryote
  (RF01960, and RF02542 for microsporidia) and prokaryote (RF00177, RF01959) SSU models: a prokaryote model's best by 20
  bits over the eukaryote's is `ssu_prokaryote`, `low`, and since 2026-11-19 set apart (`rrna/prokaryote_ssu.*`). An 18S over 10,000 bases (since 2026-11-13) is no one 18S:
  18S-like pieces each within 1,500 bases of the next chained into one span (up to 847 kb), or a published read whose
  regions were not found, taken whole (EUKARYOME's of 32,759 bases); it is `too_long`, `low`, and none of the checks
  run on it (the screen's vsearch, whose alignments cost the product of the lengths, waited an hour on one); a screen
  part over 10,000 bases is not screened. The chain's cause is fixed in the extraction (since 2026-11-13): barrnap's
  18S (and 28S) pieces within 1,500 bases are one span only while their shares of the model (barrnap's "aligned only N
  percent", a whole hit 1) add up to 1.3 at most, one molecule cut by introns covering it about once, and the span
  stays within 10,000 bases (15,000 for a 28S); two whole 18S hits are two 18S. The done genomes' saved records of
  joined pieces are extracted again from their own sequence, which holds every piece (barrnap on it alone, no genome
  needed), and kept in `WORK/rrna/`. ROD's operons are NCBI assemblies' own: a record's own
  assembly in ROD (another version, RefSeq's copy) is now its own genome, no confirmation of it (`screen.genome_of`).
- **The tiers** ([FORMAT.md](FORMAT.md#confidence-tiers-format-6)). `high` needed a relative of the copy's family
  near it, so a lineage the package holds few species of could not be `high`, and its copies were `medium` beside
  fragments and contig ends. On NCBI's tables of 2026-11-05, 5.9% of the species with a genome have no other species
  of their family (TSAR 12%, Amoebozoa 19%, Metamonada 37%, CRuMs all). `unverified` is now its own tier, between
  `high` and `medium`. The length check flagged each locus' outer 5% whatever the lengths (in the trial, with the
  undecided check, 650 of 1,769 `medium` copies were `medium` for these alone); it now compares a copy with its
  relatives' median (0.75-1.33), else its locus' (0.5-2).
- **Statistics of the right copies.** The intron sites' shares count whole genomic copies alone (a transcriptome's,
  which holds no introns, a fragment and a `low` copy no longer lower them); the genetic-code check leaves `low` copies
  out (a contaminant's codons are its own code's); the barcode gap's `resolved` is `t` against a congener alone
  (`inter_rank`; until 2026-11-09 `t` where the package held no congener of the species) and its pairs align over 80%
  of their query.

## The reference screen (since 2026-11-10)

KSGP 4.5 and 5.2 found records in their 16S and 18S that hold pieces of two phyla's genes, or another phylum's gene
under a name (`chimera_screen.py`: windows against GTDB genomes' 16S). MIDORI2, the mitochondrial reference the build
can name species by (`--midori2`), takes GenBank's records as they are named, so it carries such records too: a
contaminant's COI under its host's name, a specimen misidentified across orders, two genes joined. The same screen,
on NCBI's lineages and for any gene (`eukore/refscreen.py`), runs on MIDORI2's files before they name species; a
record it drops names none (`species.tsv`'s `midori2`), and `reference_screen.tsv` lists every call.
`tools/screen_reference.py` writes the calls as KSGP's lists (`drop.tsv`, `trim.tsv`, `review.tsv`) and the files
screened, for a user of MIDORI2 itself.

- **Windows and references.** Each record is cut into 300-base windows (every 150 bases in the first pass, every 50
  in the check). They are set beside RefSeq's organelle genes (trusted: curated genomes, every gene of them, so that
  another gene shows) and one record of each other species of the files (its longest).
- **First pass** (minimap2 `asm20`): a record is flagged where a window is 95% or more identical to a record of
  another gene or of a species NCBI's lineage parts from its own at the order or above (the first of phylum, class,
  order both name and differ in).
- **The check.** Each window's copies (97% or more; minimap2) and relatives (80% or more; vsearch), and, where a
  window looks another clade's, its own clade's records searched alone (a large clade's relatives fill a search's
  hits). A window is claimed by another clade where a record of it is a copy and 5 points nearer than the window's own
  clade's nearest record; nearer than that, a window tells the two clades apart no better than a conserved stretch
  does (plants' mitochondrial genes are 97-100% alike across orders).
- **Whose a sequence is.** A copy does not say: the copy may be the misnamed record, and so may several copies be (a
  batch of contaminated specimens are copies of their source under wrong names; KSGP's second-genome rule). A clade
  holds a sequence where it holds relatives of it: the other clade's second-best species under 97% (a RefSeq record
  counting as two) against the record's own clade's nearest record, averaged over the consecutive windows the clade
  claims (single windows' relatives at 80-90% differ by chance). 3 points nearer, and no third clade holding half the
  windows about as well, the stretch is foreign; 3 points the other way, it is the record's own (the copy is the
  misnamed one); else undecided.
- **Calls.** A foreign stretch of 2 windows or more (or one at 99%) is certain. With 2 windows of its own clade or
  more the record is a `chimera`, cut to its longest own part where 300 bases stay, else dropped; without, `misnamed`
  and dropped where the stretches are half its windows or more. Another gene's windows (RefSeq's other gene 5 points
  nearer than any of its own gene) make `other_gene` alike. Undecided stretches and uncertain pieces are `weak`: kept,
  listed for review.
- **What it does not see.** Clades NCBI does not part at the order or above (a misidentified congener, a chimera of
  two species of one order); a contaminant of a clade the references hold no relatives of; a batch of contaminated
  specimens under several names of the record's own clade.

On MIDORI2 GB273's LONGEST files of COI, CYTB, lrRNA and srRNA (511,981 records; RefSeq's organelle release of
2026-11-06; 2026-11-10, ~2 hours on 6 cores):

| gene | records | misnamed | chimera (cut) | other gene (cut) | for review |
|---|---|---|---|---|---|
| COX1 | 242,986 | 118 | 11 (9) | 6 (1) | 98 |
| CYTB | 71,358 | 29 | 12 (11) | 1 | 32 |
| lrRNA | 119,217 | 52 | 9 (7) | 3,768 (242) | 205 |
| srRNA | 78,420 | 75 | 19 (12) | 173 (95) | 100 |

3,896 records are dropped and 377 cut. Most of the drops are one fault of MIDORI2's: 9,185 records sit in both its
srRNA and its lrRNA file, and of their lrRNA copies 3,367 are SSU by RefSeq's genes (98-100% to a RefSeq SSU; their
srRNA copies pass), mostly fungal mitochondrial SSU, the barcode of Agaricomycetes (1,299) and Lecanoromycetes
(1,084), and plants'. The misnamed COI are of the kinds the screen is for: a Coleoptera record that is a hemipteran's
COI, Lepidoptera records that are dipterans', arthropods' that are vertebrates', a fish's that is a
cichlid's; the chimeras mostly long records whose one end is another phylum's or order's (an arthropod COI ending in a
vertebrate's, a rodent COI starting as an artiodactyl's).

Checked with records planted in the COI file (200 records' genes under the name of a species of another order, 200
records whose second half is another order's record's; 2026-11-10): of the misnamed, 170 dropped (162 `misnamed`, 8
`chimera`), 28 for review, 2 not flagged; of the chimeras, 104 caught (80 cut, the cut a median 29 bases from the
join; 24 dropped), 55 for review, 41 not (the other record's second half often covers another stretch of the gene,
which the record's windows do not meet). None of the 800 records the planted ones were made from was dropped, nor
any record the screen of the file alone keeps: the true owner of a sequence is not taken for its contaminant.

### The other sequences Eukore takes in

Which other outside sequences go into the package, and what checks them (2026-11-10):

| sequences | used as | checked by |
|---|---|---|
| EUKARYOME 2.1's long reads, ROD v1.2's operons | `rrna/records`, by their own labels | the part screen (`screen`) and, since 2026-11-10, `own_species` against the genomes' 18S |
| MIDORI2's files (`--midori2`) | names of species (`species.tsv`); since format 7 its srRNA file's mitochondrial SSU, one a species (`organelle/genes`), since 2026-11-15 its CO1 and Cytb files' COX1 and CYTB too | the reference screen above, protal's rules (length, ambiguous bases, a reading frame), then as every organelle gene |
| NCBI's and Swiss-Prot's proteins of the marker loci | the markers' members, their HMMs | `markers.foreign()`: a member another phylum's at near identity is dropped |
| the genomes, transcriptomes and `data/no_assembly.fna.gz`'s records | everything a genome gives | `own_species`, the part screen, FCS-GX, tiara, genetic code check |
| RefSeq's organelle genomes | organelle models; the reference screen's trusted records | none: curated whole genomes |
| PR2, KSGP, SILVA, UNITE, EUKARYOME's labels | names of species alone; PR2's and KSGP's record ids (`taxonomy/links.tsv`) | nothing to screen: Eukore reads their names, not their sequences |
| EUKARYOME's SSU set (since 2026-11-14, format 7) | `rrna/records` (nuclear 18S), `organelle/genes` (its organelle SSU) | the audit (`audit.py`, `chimeras.py`: ends, flanks, chimeras), then as its long reads: the part screen, `own_species` against the genomes' 18S, the SSU model, telomere repeats |
| BUSCO's `refseq_db.faa`, benchmarks/'s other references | miniprot's queries; the benchmarks' other databases | not screened: a seed protein names nothing, and a benchmark compares a database as it is |

- **Published reads' own species** (`rrna.own_species`). A published read's label was taken as its source gave it,
  `own_species` `-` and the read `unverified` at best; a misnamed read or a contaminant's operon in ROD's assembly was
  caught only where the part screen found another clade's copy at 99%. Now each read's 18S is set beside the genomes'
  18S alone (its source's other reads, labelled by the same hand, vouch for nothing), and the 18S' own rule decides
  (`rrna.verdict`, since 2026-11-10: not the coding sequences' `qc.verdict`; a lone lineage's 18S is undecided, not
  another clade's, unless 99.5% like another phylum's; its own clade must be 2 points under the other's; measured
  on 4,090 labelled 18S: true 18S of a species without a relative of its class called another clade's, 46 before
  and 2 after, with no loss of the contaminants caught).
- **Relatives under 90% hold nothing** (`qc.verdict`'s `related`, `rrna.RELATED`). Weighing two clades, a copy 99% like
  one family's single genome and 80% like its own family's 18S went to its own family (80 > 0), though at 80% an 18S is
  any eukaryote's. A relative under 90% now counts for neither clade (the BUSCO copies' check unchanged).
- **On the trial's records** (eukore-013: 2,350 ROD operons, 74 EUKARYOME reads): 2,397 of the 2,424 published
  reads `t`, 27 `-`, none `f`; no genome's record changed. 100 of its passing genomes' and reads' records renamed as
  a species of another phylum: 91 `f`, 9 `-`, none `t` (12 `t` without `related`: Arabidopsis' operons, the trial
  holding one Brassicaceae genome); no other record's verdict changed.

## Eukore as a user's only source of eukaryote rRNA and names, format 7 (since 2026-11-14)

protal asked that Eukore become its only source of eukaryote rRNA records and eukaryote names, so that KSGP holds
prokaryote 16S alone. Measured on KSGP 5.2 (817,884 records): 476,682 nuclear eukaryote 18S (342,515 of EUKARYOME,
94,791 of PR2 5.1.1's accessions, 7,803 of SILVA 138.2's eukaryote accessions, 20,931 of GenBank in neither, 8,805 of
NCBI's assemblies, 1,474 of RefSeq), 953 plastid 16S, 143 nucleomorph 18S and 58,806 mitochondrial SSU (MIDORI's).
Format 6 held none of PR2's, GenBank's, SILVA's or RefSeq's, and of EUKARYOME's 400,217 long reads only the 160,138
labelled by a species or genus NCBI names.

- **EUKARYOME alone** (`sources.EUKARYOME_SSU`, `eukore_build.eukaryome_reads`). KSGP 5's eukaryote 18S is EUKARYOME
  2.1's SSU set: all 342,515 of KSGP 5.2's EUKARYOME ids are in it, and 127,675 of its 134,167 other nuclear 18S under
  the same accession (PR2 94,020 of 95,113, KSGP's genome 18S all 8,805, SILVA 7,688 of 7,844, GenBank 16,020 of
  20,931, RefSeq 1,142 of 1,474; 125,651 with the same sequence, 1,474 one inside the other). The SSU set (558,632
  records) holds 407,694 EUK ids (386,813 the 18S of a long read) and 150,938 accessions (GenBank's 18S, `GCA_` genome
  18S, UNITE's). So the published 18S are EUKARYOME's alone: its long reads (but those it marks of compromised
  quality) and every record of its SSU set but the 18S of a long read taken (the read is taken whole); PR2 and KSGP
  give no record and no taxon. Their records are linked to Eukore's by accession, with their labels
  (`eukore/links.py`, `taxonomy/links.tsv`; PR2's taxonomy fetched, KSGP's `--ksgp_tax`). EUKARYOME tags a record of
  no nucleus by its kingdom `k___<tag>`: its mitochondrial (2,054), plastid (248, apicoplast 5, chromatophore 1) and
  nucleomorph (396) SSU are organelle genes; its bacterial (669) and archaeal (188) SSU and pseudogenes (7) are left
  out. Until 2026-11-14 its long reads of these tags named by an NCBI genus or species (2,010 mitochondrial, 374
  nucleomorph, 118 plastid reads among them) were nuclear operon reads.
- **Coverage** (no checks run; NCBI's taxdump of the trial build, KSGP 5.2's labels by `rrna.label_taxon`'s names):
  623,160 published records (EUKARYOME's long reads 396,176, its SSU set 157,504, ROD's operons 69,480). Of KSGP 5.2's
  59,026 NCBI species of its nuclear 18S, 57,962 hold an rRNA record (98.2%), 58,072 with the genomes' species (954
  not: Fungi 349, Metazoa 206, Alveolata 143, Viridiplantae 69, ...: KSGP's own additions of PR2, RefSeq and SILVA
  records EUKARYOME lacks); 4,065 species KSGP does not name. Of its 27,606 NCBI genera, 27,128. The records' label
  rank: EUKARYOME's species 135,638, genus 159,439, above the genus the rest (its environmental long reads mostly).
  Organelle genes of 18,663 species' mitochondria, 15,032 species' plastids and 69 species' nucleomorphs; organelle
  genes alone of 22,421 species (mitochondrion 10,972, plastid 11,759, nucleomorph 13). Of KSGP 5.2's plastid species
  513 of 527 and nucleomorph species 49 of 52 hold organelle genes; of its mitochondrial SSU species (MIDORI's) 13,084
  of 58,600 without MIDORI2's records.
- **MIDORI2's mitochondrial SSU** (`eukore_build.midori2_published`, its srRNA file among `--midori2`'s: by default the
  newest release's UNIQ files, fetched from Zenodo). One record a species whose mitochondrial genes no RefSeq genome
  gives: its longest the reference screen keeps (a cut record its kept part, a dropped one never, one left for review
  only where the species has no other, `screen_weak`), as KSGP took one MIDORI record a species (58,806 in KSGP 5.2); a
  species' every haplotype (MIDORI2's UNIQ files) would multiply the organelle checks' searches for no species more. The
  reference screen now runs before the organelle genes. On MIDORI2 GB273's LONGEST srRNA file with that release's screen
  (`~/refscreen/all.v5`: SSU misnamed 75, chimera 19, other gene 173, weak 100): 60,710 records taken (73 cut by the
  screen, 93 for review), 131 dropped, 17,402 of species with a RefSeq mitochondrial genome not taken, 166 of no NCBI
  taxon. Species with mitochondrial genes 18,663 -> 78,611, with organelle genes alone 22,421 -> 68,128 (mitochondrion
  57,033); of KSGP 5.2's mitochondrial SSU species 58,486 of 58,600 hold one. Since 2026-11-15 by protal's rules, with
  COX1 and CYTB ([below](#midori2s-cox1-and-cytb-and-one-taxonomic-framework-since-2026-11-15)).
- **The audit** (`eukore/audit.py`, `eukore/chimeras.py`; [FORMAT: rrna/](FORMAT.md#rrna)): KSGP 4.1-4.5's checks of
  their records (KSGP 5.2's and 4.5's READMEs) run on every published read at build time, calibrated on EUKARYOME
  2.1 (2026-11-14; WSL `~/euk-audit`):
  - *Primers.* EUKARYOME's records are trimmed of their primers already: 364,417 of its SSU set end 3 bases past
    EukB's site (the 18S's own last bases), 184,843 begin at one V4 site, 255,487 long reads end at one 28S site; 555
    long reads end exactly with RCA95m. A gene-end primer's site is the 18S's own end, so no primer is cut, only a
    barcode or adapter before one (bases of no rRNA, `data/primers.tsv`'s rDNA primers, 12 added).
  - *Vectors.* UniVec_Core 25-mers in 328 of the SSU set's records and 877 long reads, mostly single 25-mers of
    oligo-dT primers and poly-A vectors (SMART CDS III, the I.M.A.G.E. libraries' NotI oligo-dT, pSP64) on poly-A runs:
    K-mers of low complexity are left out, and a vector inside the rRNA counts from 50 bases (VecScreen's internal
    match: 30).
  - *Ends and flanks*, on 10,000 random long reads and 10,000 records of the SSU set without a long read (barrnap, then
    the audit with nhmmer): cut 312 long reads (an 18S read whose flank holds its 5.8S or 28S, an operon barrnap did not
    make whole, 166; a barcode before a primer 105; DNA beyond the gene 51; telomere repeats 2; vectors 2) and 289 SSU
    records (a barcode before a primer 172; DNA beyond a whole 18S 67; N tails 50; another gene beside the 18S 8;
    vectors 2); flagged `flank` 82 long reads (until the 28S was spared: 711, its 28S going on to the primer site the
    read ends at) and 646 SSU records (DNA of median 200 bases beside an 18S that stops short), `no_rrna` 32 and 189.
  - *The window screen* (KSGP 4.5), on every distinct 18S of the SSU set (484,839, its 5,763 `GCA_` genome 18S the
    references; the build adds Eukore's own genome 18S and ROD's): by phylum, as KSGP screens 16S, it called 4,987
    records mislabels and 2,128 chimeras, most of them Endogonomycetes (Mucoromycota) reads whose windows are 97-99%
    like Mortierellomycota genomes' 18S, a lineage without a genome of its own: a kingdom's phyla are as alike in the
    18S as a bacterial phylum's families in the 16S. By EUKARYOME's kingdom: 55,496 of the distinct 18S flagged by the
    first pass; of the records 392 chimeras cut to their own part, 350 chimeras flagged (their own part under 500
    bases), 144 mislabels (another kingdom's windows over half of it and none of its own: until KSGP 4.5's share rule
    was taken, 3,410, most of them Sarcomonadea (Cercozoa) reads with two windows at 98% to stramenopile genomes), 322
    recurrent, 8,496 uncertain (one window, or a conserved stretch of a lineage without a genome relative) and 4,978
    under a placeholder kingdom (no flag). Of KSGP 5.2's drops among them it calls (chimera or mislabel) 166 of 1,541
    chimeras (class Q: mostly of two genera of one kingdom, UCHIME's), 8 of 14 another organism's 18S (B), 5 of 13 genus
    homonyms (C), 6 of 80 of the wrong kingdom (U).
  - *UCHIME* (KSGP 4.2's rule), on 4,986 random records of the SSU set and 1,004 that KSGP 5.2 dropped as chimeras
    (classes Q and E), against all 555,012 of its nuclear 18S: of the random records 1 confirmed and 4 uncertain; of
    KSGP's chimeras 366 confirmed (224 cut to their own part) and 490 uncertain (their parts' best records of one genus,
    or of none named). Until the unnamed genera of EUKARYOME's environmental records were set aside, 682 of KSGP's
    chimeras were uncertain. Since 2026-11-15 a record's own sequence alone is left out of the database it is
    screened against (vsearch's `--self`): until then each part's every query was, so that with three parts a third
    of the published 18S could be no parent (a true chimera scored 35 alone, 2.5 with a parent in its part); and the
    distinct 18S are named in their names' order (a set's order, another in each process, renamed them every run and no
    kept part was found again). These figures were measured before both changes.
  - *The window screen's references* (since 2026-11-15): of each genome its 3 longest distinct 18S, and one genome for
    an assembly's every copy (its own records, ROD's, EUKARYOME's `EUK_GCA_` and `_dup` records). Every genome 18S and
    operon was one (4.43 M in the first full build, 8,140 of 22,740 genomes holding 31 or more): minimap2 skips a seed
    seen in over 500 references, and a window's 30 secondary hits fell to one genome, so the second genome the screen
    asks for never showed. The calibration above (5,763 references) has not been repeated on the build's references.
  - *Across nodes.* A read's ends and flanks are audited with its regions in the published batches, which the process
    stage's array tasks share (WORK/published/regions.4.tsv; the cache's name changed, so every read is annotated
    again once); the window screen's first pass and UCHIME's parts are searches of qc's protocol (`qc.MapParams`,
    `qc.UchimeParams`), which the `--qc_helpers` nodes take; the check of the flagged reads and UCHIME's confirmation
    (a few thousand reads) run on the final node, the confirmation's search shared too.
- **Labels above the genus** (`rrna.label_taxon`). A label's taxid is its species', else its genus' (as before), else
  the deepest of its family, order, class, phylum and kingdom NCBI names (a EUKARYOME name of a deep lineage by
  `data/eukaryome_names.tsv`'s NCBI name: Straminipila, Stramenopiles), else Eukaryota (a lineage of EUKARYOME's own,
  `Eukaryota.reg112b`, NCBI names none of). A taxon not of the eukaryotes is never taken (KSGP's `Paracoccus(Metazoa)`).
  The rank the record is labelled at is its label's (`rank`); a record labelled by a clade above its genus is
  `medium` (`higher_label`), as one labelled by its genus (`genus_label`).
- **EUKARYOME's taxonomy, one tree** (`eukore/eukaryome.py`; [FORMAT: taxonomy/](FORMAT.md#taxonomy)). On EUKARYOME 2.1's
  classification table (93,399 rows): 950 genera under several lineages, 1,014 names with several parents, rows a
  column off, a repeated header, ~40 subkingdoms the text of an Excel formula (`openpyxl` read the formulas, not their
  values, until 2026-11-14). Made one tree of 106,126 taxa, 52 of them written with a bracket above the genus
  (EUKARYOME's own oddities: `Nibbleridia(Provora)`), every name under one lineage. EUKARYOME's 400,217 reads' labels
  add 16,383 taxa (their species, and families and genera the table lacks). Of the NCBI species of KSGP 5.2's 18S and
  of NCBI's assemblies (78,329): 12,937 placed by their species, 727 by a synonym, 63,523 by their genus, 1,137 upward,
  4 nowhere (NCBI's `Eukaryota sp.`).
- **Organelle-only species**. The organelle genes of every species with a RefSeq organelle genome, not only those of
  a species with a genome or rRNA records: a species of organelle genes alone is a row of `species.tsv`, `source`
  `organelle`, with all its RefSeq genes (its rRNA, COX1, CYTB, ...), checked as every organelle gene. EUKARYOME's
  organelle SSU of a RefSeq genome whose genes the package holds are not taken again.
- **A telomere repeat** in an rRNA record (6 units or more of TTAGGG, TTAGG, TTTAGGG, TTGGGG or TTTTGGGG, in any phase,
  on either strand) is `low` (`telomere`): an rRNA joined to a chromosome's end (KSGP masked 133 bases of one in
  EUKARYOME's EUK1400052); one beside the rRNA is cut off by the audit (below).

## MIDORI2's COX1 and CYTB, and one taxonomic framework (since 2026-11-15)

protal's mitochondrial layer takes MIDORI2's COX1, CYTB and 12S (and 16S) of every species Eukore lacks
(`ribotal/scripts/mini_db/mito_records.py`), and KSGP 4.4 took MIDORI's mitochondrial SSU after screening it. Eukore
now ships MIDORI2's COX1 and CYTB beside its SSU, screened as both did, and every species, record and gene stands in
EUKARYOME's one tree, whose higher taxa are one name each.

- **MIDORI2's genes** (`organelle.midori2_records`, `eukore_build.midori2_published`). Of each named NCBI species
  without a RefSeq mitochondrial genome, one record of each of SSU, COX1 and CYTB: its longest that
  - the reference screen keeps (`refscreen.py`: another organism's gene, a chimera or another gene dropped; a foreign
    end cut; a doubtful record only where the species has no other, `screen_weak`);
  - passes protal's rules: the gene's length (SSU 250-3,000 bases, COX1 300-2,400, CYTB 300-1,600), at most 0.5%
    other bases than ACGT, a named species;
  - of a protein-coding gene, has a reading frame by the species' mitochondrial code (`organelle.reading_frame`): a
    record with a stop in every frame of both strands is a nuclear copy (a NUMT) or a pseudogene, left out; one stored
    reversed is turned; its coding sequence is its open frame, so the organelle check translates it without stops;
  - of the SSU, reads on its gene's strand by an 8-mer model of the SSU taken (KSGP 5.1's check, which turned 721 of
    MIDORI's mitochondrial SSU), else is turned.

  **One record a species screened** (since 2026-11-16, `organelle.midori2_choice`). The reference screen screened
  every record of MIDORI2's 15 UNIQ files: 4.0 million in the build of 2026-11-16 (it was measured on GB273's
  LONGEST files, 0.5 million), on the final node for hours. Now it screens, of each named NCBI species and gene, one
  record: the longest whose length is the gene's (SSU, COX1, CYTB, LSU as above; MIDORI2's other genes from about a
  third of their usual length in animals' mitochondria to about 1.6 times it: atp6 200-1,100 bases, atp8 100-400, COX2
  200-1,200, COX3 250-1,300, nad1 300-1,600, nad2 300-1,700, nad3 150-600, nad4 400-2,200, nad4L 150-500, nad5
  500-2,900, nad6 200-850: longer is several genes' span or an annotation's artifact), with at most 0.5% other bases
  than ACGT and, of a protein-coding gene, a reading frame by the species' mitochondrial code. A species' gene taken
  above is that record or none: where the screen drops it, the species has none of that gene (until then its next
  longest), and a species with no record of the gene's length has no record screened. The other records of a species
  name it in `species.tsv`'s `midori2` lineage as before, unscreened. The choice is kept beside the screen's calls
  (`WORK/reference_screen/<key>/chosen.tsv`). The length windows of the genes other than SSU, COX1, CYTB and LSU are
  set from animals' genes and not measured on MIDORI2's files.

  On MIDORI2 GB273's LONGEST files with that release's screen (`~/refscreen/all.v5`), of the species of no RefSeq
  mitochondrial genome:

  | | COX1 | CYTB | SSU |
  |---|---|---|---|
  | species taken | 213,812 | 50,844 | 56,311 |
  | of them cut by the screen, left for review | 8, 86 | 6, 29 | 72, 89 |
  | turned to the gene's strand | 2 | 0 | 1,067 |
  | left out: ambiguous bases | 8,792 | 1,901 | 1,706 |
  | left out: not the gene's length | 1,750 | 921 | 2,695 |
  | left out: no reading frame (NUMT, pseudogene) | 52 | 27 | - |
  | left out: dropped by the screen | 115 | 28 | 131 |
  | left out: no NCBI taxon, no named species | 806, 52 | 113, 3 | 166, 3 |

  The SSU's length floor is protal's: 2,732 of the srRNA file's records are shorter, 1,776 of them 150-199 bases
  (MiFish's 12S amplicon); until 2026-11-15 they were taken (60,710 SSU species then). Of the package's
  291,060 species (this measurement's) 262,590 hold mitochondrial genes (17,849 a RefSeq genome's), 254,536 hold
  organelle genes alone, 229,879 of them MIDORI2's alone. The organelle checks (`organelle.check`) then search
  COX1 as one locus of some 230,000 records; their run time on the cluster is not yet measured.
- **One taxonomic framework** (`eukore/eukaryome.py`; [FORMAT: taxonomy/](FORMAT.md#taxonomy)). EUKARYOME is
  Eukore's primary taxonomy: every species', record's and organelle gene's lineage is its place in EUKARYOME's one
  tree, NCBI's names only below the deepest taxon EUKARYOME holds. Measured on the package's 291,060 NCBI species
  (EUKARYOME's, ROD's, the genomes', RefSeq's organelle genomes' and MIDORI2's SSU, COX1 and CYTB; no build, the
  trial's taxdump and EUKARYOME 2.1; WSL `~/euk-audit/taxonomy.json`) as main placed them (0.2.6), a higher taxon
  was at times another lineage's or two names. Four faults, each fixed:
  - *A name EUKARYOME misfiles under another phylum.* EUKARYOME 2.1 files the golden moles as
    `Mammalia;Alcyonacea;Chrysochloridae`, and the tree, one Alcyonacea under the parent of most rows, took them into
    the soft corals' Anthozoa. A name below the subphylum filed under two phyla of a kingdom, each borne out by the
    name's parent there, is now two taxa: `Alcyonacea(Mammalia)`, the homonyms `Triadinium(Alveolata)` (a ciliate)
    and `Triadinium(Goniodomataceae)` (a dinoflagellate, which EUKARYOME's reads label `Triadinium(Dinoflagellata)`:
    a label goes to the taxon sharing the most of its names), Caudospora (a fungus, a microsporidian) and
    Bolidophyceae (under Bolidophyta for Bolidomonas, Bacillariophyta for Triparma): 4 of the table's 10 names under
    two phyla. A row whose phylum alone is another (Pezizomycotina under Basidiomycota,
    `Coniferophyta;Tracheophytina;...`: a column off) stays one taxon, its names under the parent of most rows.
  - *A genus name of another phylum or class.* A species placed by its genus took any genus of the name its kingdom
    allows: NCBI's golden mole Zebra under EUKARYOME's coral Zebra, insects under barnacle genera, bryozoans under a
    bug's, a slug under an aphid's. A genus whose phylum, or class in NCBI's phylum but of another subphylum, is not
    NCBI's, where EUKARYOME holds NCBI's, NCBI's taxon stands not above it and no name below it is NCBI's, is not
    taken (`Tree.conflicts`); the species is placed upward. A genus EUKARYOME classifies otherwise within NCBI's
    taxon keeps EUKARYOME's place: Galdieria in Cyanidiophyta (under EUKARYOME's subkingdom Rhodophyta; NCBI's
    Bangiophyceae), Neocystis in Trebouxiophyceae (NCBI's Chlorophyceae, both Chlorophytina), Candida in
    Saccharomycetes (NCBI's Pichiomycetes). A first rule without these exceptions moved 2,527 species, most of them
    out of EUKARYOME's classification into NCBI's.
  - *A family EUKARYOME writes with its order* (`Chilodontidae.Seguenziida`, the gastropods' beside the fishes'
    Chilodontidae) was not found by its name: a gastropod went under the fishes.
  - *An NCBI clade EUKARYOME names otherwise.* A species placed upward went under NCBI's own name where EUKARYOME
    lacks it: NCBI's Copepoda became a second class beside EUKARYOME's Hexanauplia, which holds the copepods
    EUKARYOME names. Each NCBI genus, family, order, class and phylum whose species placed by name stand under one
    EUKARYOME taxon of its rank (a taxon of under 1% of them a stray) is now that taxon (`Tree.learn`; 64,946 NCBI
    clades): Copepoda is Hexanauplia, NCBI's Thecostraca EUKARYOME's Hexanauplia too, which holds its barnacles. A
    taxon of the clade's name below it is taken before it (NCBI's phylum Microsporidia is EUKARYOME's Rozellomycota,
    its subphylum Microsporidia). A clade EUKARYOME splits (NCBI's Pinopsida into Pinopsida and Cupressopsida) is
    not learned: a taxon above its parts left the ranks between empty (the Triphoridae lost their Gastropoda), so it
    is found by its name as before. Since 2026-11-15 the build learns before the audit places the published reads, and
    learning forgets what earlier placements added under NCBI's names (`Tree.forget_placed`): the audit placed the
    reads before the package learned, a copepod genus EUKARYOME lacks stood under a class `Copepoda`, its species
    voted for it beside Hexanauplia, and Copepoda was learned as nothing. The figures below were measured without a
    build, so without that order.

  | NCBI taxa whose species stand under two or more EUKARYOME taxa of their rank (species not in the largest) | main (0.2.6) | now |
  |---|---|---|
  | phylum (of 77) | 14 (670) | 10 (641) |
  | class (of 288) | 48 (615) | 42 (574) |
  | order (of 1,277) | 183 (1,883) | 173 (1,849) |
  | family (of 7,456) | 556 (5,618) | 490 (5,305) |
  | genus (of 67,198) | 315 (463) | 315 (463) |

  391 species (0.13%) changed lineage, first at the subkingdom 18, phylum 26, subphylum 13, class 4, order 14,
  family 315, genus 1: the golden moles out of the corals (18), insects out of barnacle genera (10), Triparma
  under EUKARYOME's own row for it (10), bryozoans, a slug and molluscs out of arthropod genera, placeholders
  that EUKARYOME names (`Placozoa.phy.incertae.sedis` its Placozoa, `Fornicata.phy.incertae.sedis` its
  Fornicatia), and below the class species placed upward under their NCBI clade's EUKARYOME taxon.

  Every name stands on one lineage (0 names under two parents, before and after), every species' EUKARYOME kingdom
  holds its NCBI group, and of EUKARYOME's records labelled above species the taxon of their label and the one their
  NCBI taxon is placed at differ for 6,467 of 367,166 (7,003 in main; a record's lineage is its label's). What stays
  split is EUKARYOME's classification, which Eukore keeps: NCBI's Streptophyta are EUKARYOME's Tracheophyta, Bryophyta,
  Marchantiophyta and Zygnematophyta, its Euglenozoa Kinetoplastia, Euglenida and Diplonemia, its Rhodophyta also
  Cyanidiophyta, its Annelida also Myzostomida; NCBI's Pinopsida EUKARYOME's Pinopsida and Cupressopsida, its
  Eurotatoria Monogononta and Bdelloidea, its Clitellata also Branchiobdellae; Rhabditida and Strongylida, Noctuidae
  and Erebidae below. A few are EUKARYOME's own misfilings within a class, which Eukore cannot tell from a
  reclassification and keeps as EUKARYOME has them: the geometrid Phigalia under Elachistidae, the scarab Heteronyx
  under the wasps' Pompilidae, the skipper Psolos under Nymphalidae, the jumping spider Tasa under Anyphaenidae (one
  to six species each).
