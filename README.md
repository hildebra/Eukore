# EUKDB

EUKDB is a standalone eukaryote reference, and this repository builds it. It holds:
- the euk30 marker genes, complete with all their introns;
- the rDNA operons, complete and uncut;

of every eukaryote species with an NCBI genome. Both are called in NCBI's own assemblies, checked, and labelled by
NCBI's taxonomy, tied to no GTDB release or KSGP version. The package format is in `FORMAT.md` (version 2).

A user converts a package to its own needs. protal's database build (`download_ref_data.py --eukdb`) does it with its
own reader of the format: KSGP names, masking, copies fitted into 16 kb, records cut to 8,191 bases. EUKDB imports
nothing from protal. It started inside protal's repository (`EUKDBbuild/`, 2026-11-05) and keeps that history.

## Install

```
micromamba create -f envs/eukdb.yaml      # Python, miniprot, HMMER, barrnap, ITSx, vsearch, p7zip, openpyxl, FCS-GX
micromamba activate eukdb
micromamba create -f envs/tiara.yaml      # optional: tiara apart, as it needs a Python of 3.9 or older
```

## How the euk30 copies are made

EUKDB calls the 30 loci itself, as compleasm does, on BUSCO's eukaryota_odb12 (OrthoDB 12.1; the 30 loci are among
its 129 BUSCOs). EukDetect1 and EukDetect2 were its sources of copies until 2026-11-05.
1. **The genome set.** For each NCBI species with a GenBank assembly in the `--groups`, the assembly that ranks best
   (`eukdb/genomes.py`):
   - never one NCBI flags as contaminated, of an unverified source organism, misassembled, chimeric or of low
     quality;
   - another flagged one, or a metagenome's, only without better;
   - then a RefSeq reference, then a representative genome; then the better assembly level, fewer contigs per
     megabase, newer.

   Its RefSeq copy is fetched where it has one (RefSeq screens eukaryote genomes for contamination). A model species
   (`data/model_species.tsv`) takes its reference assembly. NCBI's placeholder species and interspecies hybrids
   ("A x B") are left out.
2. **The calls** (`eukdb/odb.py`, `eukdb/loci.py`).
   - miniprot aligns, per locus, the proteins of the 30 dataset species nearest to the genome in NCBI's taxonomy, and
     BUSCO's ancestral variants.
   - An alignment needs identity 0.4 and coverage 0.6 of its protein. The overlapping ones form one region, whose best
     alignment is its gene.
   - hmmsearch confirms a region when the locus' own HMM scores it best and at BUSCO's cutoff.
   - A locus is Single, Duplicated (the best copy taken), Fragmented (under 80% of the profile) or Missing.
3. **The copy.** The gene from first to last codon, introns and all, uncut and unmasked. Its features are its exons
   and introns, the assembly's soft-masked runs and its low-complexity runs (DUST); the user masks as it sees fit.
   - The 50 kb around each copy are kept for contamination screens.
   - tiara classifies them where installed (`--tiara`).
   - NCBI's FCS-GX screens them with `--fcs_gx_db` (its ~470 GB database; a node of ~512 GB). A copy overlapping a
     range FCS-GX calls another organism's gets an `fcs` column.
4. **The checks across species** (`eukdb/qc.py`), each a column of `euk30/copies.tsv`. They cover what EukDetect
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

## Layout

| path | contents |
|---|---|
| `eukdb_build.py` | the build |
| `eukdb/odb.py`, `eukdb/loci.py` | BUSCO's dataset; the loci called in a genome |
| `eukdb/genomes.py` | the genome set, one worker per genome (loci, windows, tiara, rRNA), FCS-GX, the prune |
| `eukdb/qc.py` | the copies' checks across species |
| `eukdb/rrna.py`, `eukdb/screen.py` | the rRNA records and their checks; the part screen |
| `eukdb/package.py` | collect, write, read and check a package |
| `eukdb/ncbi.py`, `eukdb/sources.py`, `eukdb/species.py` | NCBI's taxonomy; the sources; the model species |
| `eukdb/placement.py`, `eukdb/copies.py`, `eukdb/util.py` | genome files; a copy's exons and introns; I/O |
| `data/` | the loci table, the model species table, curated rDNA units, the no-assembly records |
| `envs/` | the micromamba environments |
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
  - BUSCO's eukaryota_odb12;
  - EUKARYOME 2.1 and ROD v1.2.
- **process** calls every genome it was given.
- **package**, in turn:
  - runs FCS-GX where asked;
  - checks the copies across species;
  - annotates the published rRNA reads (barrnap and ITSx over `-t` threads, 250 reads a batch);
  - runs the rRNA checks and part screen;
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

## The loci in prokaryotes

EukDetect2 dropped the BUSCOs it found in 971 gut bacteria, and the 30 loci are among those it kept. On 2026-11-05 that
check was redone on EUKDB's own caller. The panel was 176 prokaryote genomes, one per GTDB genus of protal's local
inputs: 158 bacteria of 22 phyla and 18 archaea of 6. Each was searched with all ~23,000 dataset proteins of the 30
loci, which is stricter than a build (the 30 nearest species' proteins). The results:

| locus | OrthoDB name | called at the cutoff (the caller's gates) | without the identity/coverage gates |
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
