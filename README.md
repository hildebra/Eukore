# EUKDB

EUKDB is a standalone eukaryote reference, and this repository builds it. It holds:
- every BUSCO gene of BUSCO's eukaryota_odb12 (129), complete with all its introns;
- the rDNA operons, complete and uncut;
- the mitochondrial and plastid rRNA and the mitochondrial COX1, where NCBI has them;

of every eukaryote species with an NCBI genome. The first two are called in NCBI's own assemblies and checked; the
organelle genes come from NCBI's RefSeq organelle genomes, else the assemblies. All are labelled by NCBI's taxonomy, tied
to no GTDB release. Each species also carries its lineage in PR2, EUKARYOME and KSGP, and each copy and record a
confidence tier, its checks and the hints of a chimera. The package format is in `FORMAT.md` (version 4).

A user converts a package to its own needs. protal's database build (`download_ref_data.py --eukdb`) does it with its
own reader of the format: the 30 BUSCOs of its euk30 (its choice, by `busco/loci.tsv` and "Choosing loci" below),
KSGP names, masking, copies fitted into 16 kb, records cut to 8,191 bases. EUKDB imports
nothing from protal. It started inside protal's repository (`EUKDBbuild/`, 2026-11-05) and keeps that history.

## Install

```
micromamba create -f envs/eukdb.yaml      # Python, miniprot, HMMER, barrnap, ITSx, vsearch, minimap2, p7zip, openpyxl, FCS-GX
micromamba activate eukdb
micromamba create -f envs/tiara.yaml      # optional: tiara apart, as it needs a Python of 3.9 or older
```

## How the BUSCO copies are made

EUKDB calls every BUSCO of eukaryota_odb12 itself (OrthoDB 12.1, 129 BUSCOs; `--loci` restricts a trial), as
compleasm does. A user picks its loci among them: `busco/loci.tsv` gives each one's call rate over the package's
species and per supergroup (protal picks its euk30). A genome searched for fewer BUSCOs than the build calls (a run
before format 4) is run again. EukDetect1 and EukDetect2 were its sources of copies (30 loci) until 2026-11-05.
1. **The genome set.** For each NCBI species with a GenBank assembly in the `--groups`, the assembly that ranks best
   (`eukdb/genomes.py`):
   - never one NCBI flags as contaminated, of an unverified source organism, misassembled, chimeric or of low
     quality;
   - never one that NCBI's own FCS-GX screen finds mostly another organism's (50% or more), or that fails NCBI's ANI
     species check (Saccharomycotina against their type strains; `eukdb/quality.py`);
   - another flagged one, or a metagenome's, only without better. EUKDB adds its flags: `partial genome` (NCBI's
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
2. **The calls** (`eukdb/odb.py`, `eukdb/loci.py`).
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
3. **The copy.** The gene from first to last codon, introns and all, uncut and unmasked. Its features are its exons
   and introns, the assembly's soft-masked runs and its low-complexity runs (DUST); the user masks as it sees fit.
   - The 50 kb around each copy are kept for contamination screens.
   - tiara classifies them where installed (`--tiara`).
   - NCBI's FCS-GX screens them with `--fcs_gx_db` (its ~470 GB database; a node of ~512 GB). A copy overlapping a
     range FCS-GX calls another organism's gets an `fcs` column.
4. **The checks across species** (`eukdb/qc.py`), each a column of `busco/copies.tsv`. They cover what EukDetect
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
(see "The loci in prokaryotes" below for EUKDB's own check).

## How the rRNA records are made

- **Genomes.** Each genome's 18S spans and whole rDNA units, from the 18S start to the 28S end (18S, ITS1, 5.8S,
  ITS2, 28S). A 5S is included only where it lies in the same block. barrnap finds the genes; ITSx sets the
  ITS1/5.8S/ITS2/28S borders, searching the profile sets of the genome's lineage.
- **Published records.** EUKARYOME 2.1's long reads and ROD v1.2's operons, labelled by NCBI taxid.
- **Checks**, each a column of `rrna/records.tsv`:
  - ambiguous bases and homopolymers;
  - the 99% 18S type within the species;
  - `own_species`, against the package's other sources;
  - the part screen (`eukdb/screen.py`, `screen`): a record whose 18S, junction or 28S half is another organism's,
    by the nearest records of other genomes on NCBI's lineages (a contaminant contig, a mosaic, a chimeric read).

## How the organelle genes are made

`eukdb/organelle.py`, for each species of the package and each of its organelles, from one source:
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

## What else a package tells its user (format 4)

So that a user (protal's build) has what it needs without computing it again:
- **Names.** `species.tsv` gives each species' NCBI lineage and supergroup, and its lineage in PR2 (5.1.1), EUKARYOME
  (2.1) and KSGP (`--ksgp_tax`, else the 4.0 KSGP's site serves), with how it matched: its name, a synonym, its genus,
  or none (`eukdb/labels.py`). NCBI's taxids merged since are in `taxonomy/merged.tsv`.
- **Gene structure.** The BUSCO copies' exons and introns (`busco/features.tsv`), the 18S's introns on each rRNA record
  and its regions (18S, ITS1, 5.8S, ITS2, 28S, an inline 5S), the organelle genes' exon counts.
- **Copies in the genome.** Per BUSCO copy its complete copies in the genome (`copies`); per genome its 18S spans and
  rDNA units in the assembly, those with a whole 28S and the contigs holding them (`genomes.tsv`), a lower bound where
  the array is collapsed; per organelle gene its copies on the organelle genome.
- **Chimeras.** Each copy's and each genome record's distance to its contig's end and its assembly gaps, beside
  `own_species`, `same_as`, tiara's window, FCS-GX and the rRNA part screen.
- **Confidence.** `high`, `medium` or `low` per copy and per record, with the flags that set it (FORMAT.md): a user
  drops the `low` ones and weighs the `medium` ones as it sees fit.

## Layout

| path | contents |
|---|---|
| `eukdb_build.py` | the build |
| `eukdb/odb.py`, `eukdb/loci.py` | BUSCO's dataset; the loci called in a genome |
| `eukdb/genomes.py` | the genome set, one worker per genome (loci, windows, tiara, rRNA), FCS-GX, the prune |
| `eukdb/qc.py` | the copies' checks across species |
| `eukdb/rrna.py`, `eukdb/screen.py` | the rRNA records and their checks; the part screen |
| `eukdb/organelle.py` | the mitochondrial and plastid rRNA and COX1 |
| `eukdb/quality.py` | NCBI's FCS-GX and ANI reports, as the genomes' flags |
| `eukdb/package.py` | collect, write, read and check a package; the copies' confidence |
| `eukdb/labels.py` | each species' PR2, EUKARYOME and KSGP lineage; NCBI's lineage and supergroup |
| `eukdb/ncbi.py`, `eukdb/sources.py`, `eukdb/species.py` | NCBI's taxonomy; the sources; the model species |
| `eukdb/placement.py`, `eukdb/copies.py`, `eukdb/util.py` | genome files; a copy's exons and introns; I/O |
| `data/` | the loci table, the model species table, curated rDNA units, the no-assembly records |
| `envs/` | the micromamba environments |
| `eukdb/depth.py`, `tools/rdna_depth.py`, `data/depth_runs.tsv` | rDNA and organelle copies by read depth; curated runs |
| `tools/` | `no_assembly_records.py` (run once) |
| `tests/` | `python3 -m unittest discover -s tests -t tests` |

## Building

Three stages, so that the network-bound one runs where the network is:

```
python3 eukdb_build.py --work WORK --out OUT --stage fetch --groups fungi,protozoa                 # network
python3 eukdb_build.py --work WORK --out OUT --stage process --groups fungi,protozoa -t 64 -j 16   # no network
python3 eukdb_build.py --work WORK --out OUT --stage package --groups fungi,protozoa -t 64 --keep_gb 1000 \
    [--fcs_gx_db /dev/shm/gxdb/all]
```

`--stage all` (the default) runs the three in turn.
- **fetch** gets the sources, then the genomes:
  - NCBI's GenBank and RefSeq assembly tables of the `--groups`;
  - NCBI's taxdump;
  - NCBI's FCS-GX summaries of its GenBank and RefSeq assemblies, and its Saccharomycotina ANI report (~70 MB);
  - BUSCO's eukaryota_odb12;
  - EUKARYOME 2.1 and ROD v1.2;
  - RefSeq's organelle genomes (mitochondrion and plastid);
  - PR2 5.1.1's taxonomy (5 MB) and KSGP 4.0's (243 MB; with `--ksgp_tax FILE` that KSGP release's tax file instead,
    protal's 4.4 or 5.2), for `species.tsv`'s lineages.
- **process** calls every genome it was given, then takes the organelle genes.
- **package**, in turn:
  - runs FCS-GX where asked;
  - checks the copies across species;
  - annotates the published rRNA reads (barrnap and ITSx over `-t` threads, 250 reads a batch);
  - runs the rRNA checks and part screen;
  - takes the organelle genes of any genome not searched yet;
  - prunes `OUT/genomes` to `--keep_gb` (the model species first, then the smallest);
  - writes `OUT/EUKDB_<version>/` (`--version`, YYYY.MM by default).

Each genome's results go into `WORK/done/<accession>/`, so a rerun continues where the last one stopped.

`--groups` (default all six: fungi, protozoa, plant, invertebrate, vertebrate_mammalian, vertebrate_other) chooses
NCBI's assembly groups. At NCBI's tables of 2026-11-05:

| groups | species | genome size (one per species) | gzipped download |
|---|---|---|---|
| fungi and protozoa | 8,772 | 445 Gb | ~135 GB |
| all six | 29,114 | 21.6 Tb | ~6.5 TB |

`--species 4930,7227,...` (NCBI taxids, or a file of one per line) builds a trial package of the species inside those
clades alone; the manifest's `inputs.partial` names the clades.

A package is checked against its manifest by `eukdb.package.check(folder, deep=True)` (None when whole).

## A trial (2026-11-05)

Seven clades: *Saccharomyces*, *Candida albicans*, *Blastocystis*, *Drosophila melanogaster*, *C. elegans*,
*Arabidopsis*, *Trichuris trichiura*. On 6 threads:
- 42 species against EukDetect2's 19, 1,150 loci against its 491. The genome stage took 8.6 minutes, the published
  rRNA 21.
- Of the 19 species both hold, EUKDB called 473 of EukDetect2's 491 loci, and 50 more.
- Where the genome is the same (the yeasts, *C. albicans*), its copies hold all of EukDetect2's copy's 21-mers, at
  the same length.
- The model species got 29-30 loci each.

## Choosing loci

The package holds every BUSCO; a user picks the ones it indexes. `busco/loci.tsv` gives each BUSCO's call rate over all
the package's species and per supergroup, its duplications and fragments, and its lengths. Three tools measure more:
- `tools/rank_loci.py` ranks the BUSCOs by:
  - presence in BUSCO's own dataset;
  - conservation: the dataset proteins' bits per residue against the BUSCO's HMM;
  - EukDetect2's fungus and protist shares;
  - call rates in a eukaryote panel (`tools/eukaryote_panel.py`: one genome per order per supergroup);
  - prokaryote calls (`tools/prokaryote_panel.py`).

  `--select N` writes the best N as a loci table.
- `tools/prokaryote_separation.py` measures how far a BUSCO's prokaryote homologues are from its eukaryote copies at
  the nucleotide level.

What the 2026-11-05 checks showed:
- **BUSCO's dataset overstates universality outside animals, fungi and plants.** It holds 837 species: 479 animals,
  201 fungi, 119 plants, 22 SAR, 6 Discoba, 6 Amoebozoa. Of the BUSCOs its presence and conservation ranked best, four
  of the seven not in euk30 were missing in *Blastocystis*' genomes. Choose by call rates in real genomes (the package's
  `busco/loci.tsv`, a eukaryote panel), not by the dataset's species list.
- **Prokaryote homologues do not matter for read profiling.** The best nucleotide identity of any prokaryote
  homologue to a eukaryote copy was 56% (the archaeal ribosome biogenesis protein); for the cytochrome c oxidase
  assembly proteins it was 50-53%, and most reached no eukaryote copy at 50%. A read aligner needing ~85% cannot cross
  between domains on them.

## rDNA copies by read depth (2026-11-06, a prototype)

An assembly collapses an rDNA array into a few units, so the units found in it are a lower bound. `tools/rdna_depth.py`
measures the copies from the genome's own sequencing reads instead (`eukdb/depth.py`; columns in FORMAT.md's
`depth.tsv`). It runs on a built package, apart from the build, in three resumable steps; only the first needs the
network:

```
python3 tools/rdna_depth.py stage   --package OUT/EUKDB_<v> --work DWORK --species 4932,6239 -t 4   # network
python3 tools/rdna_depth.py measure --package OUT/EUKDB_<v> --work DWORK -t 16                     # no network
python3 tools/rdna_depth.py table   --package OUT/EUKDB_<v> --work DWORK --into_package
```

1. **The run.** The genome's ENA run, found in this order:
   - a curated run (`data/depth_runs.tsv`: the old reference assemblies, whose BioSamples hold no reads);
   - the runs of the assembly's BioSample;
   - its BioProject's runs of its taxid or species;
   - its species' runs with the assembly's strain name.

   Only WGS of genomic DNA, randomly fragmented. Illumina first, then PacBio, then Nanopore; then the run with the
   most bases.
2. **The reads** (`stage`). Read 1 of the run, up to `--depth` (default 5) times the genome size or 20 M reads. A
   file of 5 Gb or less is read whole and thinned (every nth read). A larger one gives its first reads. Some runs were
   submitted as aligned BAMs and are sorted by position. One such run (*S. paradoxus*' ERR1938684) gave a median
   single-copy depth of 0 from its first reads.
3. **The mapping** (`measure`). minimap2 maps the reads to the genome's rDNA units (one per 18S type), its Single
   BUSCO copies and its species' organelle records. An alignment counts if it is 95% identical (80% for PacBio CLR
   and Nanopore), and the read's unaligned ends lie past the target's ends (a paralog's read is clipped inside).
4. **The copies.** The rDNA depth is the median over each unit's 18S, 5.8S and 28S, summed over the units. It is
   divided by the median of the single-copy BUSCOs' mean depths. Both alleles map, so the number is per haploid
   genome. Organelle genomes are counted the same way.
5. **The package** (`table --into_package`). `depth.tsv` is written into the package. `genomes.tsv` gets
   `rdna_copies_reads`, `mt_copies_reads` and `pt_copies_reads`, and the manifest's files are updated.

**Trial** (2026-11-06, the format-4 trial package, ~10x of reads per genome, Illumina):

| genome | run | rDNA copies (18S alone) | published |
|---|---|---|---|
| *S. cerevisiae* S288C | SRR1972546 (MiSeq) | 99 (100) | ~150 (100-200; Petes 1979) |
| | SRR2977196 (HiSeq), a second library | 131 (132) | |
| *C. albicans* SC5314 | SRR34991648 (NovaSeq X) | 67 (67) | varies by strain, ~20-180 (Rustchenko et al. 1993) |
| | SRR34991638 (NextSeq), a second library | 52 (51) | |
| *C. elegans* N2 | SRR19445835 (NovaSeq) | 58 (56) | ~55 (Ellis et al. 1986) |
| *A. thaliana* Col-0 | ERR13994817 (NovaSeq X) | 728 (860) | ~570-750 (Copenhaver & Pikaard 1996) |
| *D. melanogaster* iso-1 | SRR13070227 (NextSeq, whole flies) | 350 (334) | ~200-250 per X or Y array |

- **How close.** C. elegans and Arabidopsis fall on the published values. Yeast's lie at 0.7-0.9 of the classic ~150.
  *Drosophila*'s lies above its arrays' published size. Its three rDNA types are summed. A pooled count (all bases
  on all units over one unit's genes) gave 379, so the sum is not the cause.
- **Between libraries.** Two libraries of the same strain differed by up to 30% (yeast 99 and 131). The sampling
  error of the single-copy median is ~4% (~29 loci, spread 0.2-0.25). So a library's own bias, and real variation
  between a strain's stocks (rDNA arrays expand and contract), dominate. A number is therefore one library's.
  Format 4's ~100 single-copy BUSCOs narrow only the sampling part.
- **Organelle copies are not validated.** Short-read libraries cover AT-rich DNA poorly. Yeast's mtDNA (83% AT) ran
  from 7x to 60x in 5 kb windows, tracking its AT, where its nuclear genes had 8x. It gave 1-3 copies, against
  ~20-50 per cell. Such numbers are flagged `mt_at_rich`. They also depend on the tissue sequenced.
- **Automatic discovery.** Five more genomes of the trial package were tried without curated runs:
  - *S. paradoxus*: by its BioSample, thinned; 66 copies.
  - *T. trichiura*: by its BioProject; 425 copies.
  - *Blastocystis* ST1: by its BioSample; `low_depth`, as a xenic culture's reads are mostly bacterial.
  - Two found none. *S. kudriavzevii*'s assembly sample holds no reads (a Darwin Tree of Life-style sample whose runs
    lie on the specimen's own samples), and *Blastocystis* ST3 has no runs at ENA.
- **Cost.** The five genomes' reads at 10x were 2.4 GB gzipped and took ~5 minutes to stage. Mapping took 82 s on 5
  threads. A full run's cost is dominated by the download: `--depth` times the genome size per genome, capped at
  20 M reads, plus whole files of up to 5 Gb where thinned. With format 4's ~100 single-copy BUSCOs, a depth of 2-3
  should do.

## The loci in prokaryotes

EukDetect2 dropped the BUSCOs it found in 971 gut bacteria, and the 30 loci are among those it kept. On 2026-11-05 that
check was redone on EUKDB's own caller. The panel was 176 prokaryote genomes, one per GTDB genus of protal's local
inputs: 158 bacteria of 22 phyla and 18 archaea of 6. Each was searched with all ~23,000 dataset proteins of the 30
loci, which is stricter than a build (the 30 nearest species' proteins). The results:

| locus | OrthoDB name | called at the cutoff (compleasm's gates, the caller's until 2026-11-06) | without the identity/coverage gates |
|---|---|---|---|
| 10262062at2759 | ribosome biogenesis protein | 10 of 18 archaea (methanogens, haloarchaea, Thermoplasmata), 1.9-2.3x the cutoff, full length | 17 archaea |
| 1704689at2759 | cytochrome c oxidase assembly protein | 1 Pseudomonadota | 7 Pseudomonadota |
| 1726137at2759 | cytochrome c oxidase assembly | 1 Pseudomonadota | 4 Pseudomonadota |
| 1925287at2759 | tRNA m1A58 methyltransferase | 1 Desulfobacterota (half the profile) | 3 Desulfobacterota, 2 archaea |
| the other 26 | | none | 6 more calls, 5 loci, 0.3-0.8 of the profile |

These are the loci eukaryotes share with archaea, and the mitochondrion's alphaproteobacterial genes. A gut-bacteria
panel holds neither.
- **In a build**, a eukaryote's own copy outscores such a homologue. A prokaryote contig's copy is taken only where
  the eukaryote's own is missing; `own_species`, tiara and FCS-GX then catch it.
- **For a read profiler**, the nucleotide distance across domains is far beyond read alignment.
- **Still, these four loci** are the ones to watch, or to replace by BUSCOs without prokaryote homologues.

## Not done yet
- FCS-GX has been run with stand-ins only (its database is ~470 GB).
- Whether to replace the four loci with prokaryote homologues ("The loci in prokaryotes") by other odb12 BUSCOs.
- The organelle genes (format 3, 2026-11-05) have been built in tests only, not on NCBI's release; they have no checks
  across species (own_species, contaminants), no confidence and no exon positions (their count alone).
- Format 4 (2026-11-05) has been built in tests only. PR2's and KSGP's lineages have not been checked on their real
  files yet; EUKARYOME's were (13 of the 2026-11-05 trial's 43 species by name, the rest by genus).
- `genomes.tsv`'s rDNA counts are the assembly's (a lower bound); `tools/rdna_depth.py` measures a genome's own by its
  reads ("rDNA copies by read depth"). The intergenic spacer, 5S arrays away from the unit and the 28S's introns are
  not taken.
- rDNA copies by read depth (2026-11-06) has run on 10 genomes, not the release:
  - it is not yet a stage of `eukdb_build.py`;
  - PacBio and Nanopore runs are untested;
  - runs that lie on a specimen's other samples (Darwin Tree of Life) are not found;
  - a GC correction (rDNA against the BUSCOs' GC, in `depth/<accession>.targets.tsv`) is not made;
  - organelle copies are unvalidated.
