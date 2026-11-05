# EUKDBbuild

This folder builds EUKDB, a standalone eukaryote reference. It holds the euk30 marker genes, complete with all their
introns, and the rDNA operons, complete and uncut, of every eukaryote species with an NCBI genome. Both are called in
NCBI's own assemblies, checked, and labelled by NCBI's taxonomy, tied to no GTDB release or KSGP version. The package
format is in `FORMAT.md`.

protal's database build converts a package to its own needs: KSGP names, masking, copies fitted into 16 kb, records cut
to 8,191 bases. protal downloads no eukaryote data of its own when given a package.

EUKDBbuild imports nothing from protal's `scripts/`. It started inside protal's repository on 2026-11-05, and is meant
to move into a repository of its own.

## How the euk30 copies are made (since 2026-11-05)

EukDetect1 and EukDetect2 were EUKDB's sources of euk30 copies until 2026-11-05. EUKDB now calls the loci itself, as
compleasm does, on BUSCO's eukaryota_odb12 (OrthoDB 12.1; the 30 loci are among its 129 BUSCOs):
1. **The genome set.** For each NCBI species with a GenBank assembly in the `--groups`, the assembly that ranks best
   (`eukdb/genomes.py`):
   - never one NCBI flags as contaminated, of an unverified source organism, misassembled, chimeric or of low
     quality;
   - another flagged one, or a metagenome's, only without better;
   - then a RefSeq reference, then a representative genome; then the better assembly level, fewer contigs per
     megabase, newer.

   Its RefSeq copy is fetched where it has one. A model species (`data/model_species.tsv`) takes its reference
   assembly. NCBI's placeholder species are left out.
2. **The calls.** miniprot aligns, per locus, the proteins of the 30 dataset species nearest to the genome in NCBI's
   taxonomy, and BUSCO's ancestral variants. An alignment needs identity 0.4 and coverage 0.6 of its protein; the
   overlapping ones form one region, whose best alignment is its gene. hmmsearch confirms a region when the locus' own
   HMM scores it best and at BUSCO's cutoff. A locus is Single, Duplicated (the best copy taken), Fragmented (under
   80% of the profile) or Missing.
3. **The copy.** The gene from first to last codon, introns and all, uncut. Its features are its exons and introns,
   the assembly's soft-masked runs, and its low-complexity runs (DUST); the user masks as it sees fit. The 50 kb
   around it are classified by tiara where installed (`window_class`).
4. **The checks across species** (`eukdb/qc.py`), each a column of `euk30/copies.tsv`. They cover what EukDetect
   checked with its alien index, its 99% collapse and its length filter:
   - `own_species`: the copy's nearest copies of other species agree with its lineage, or are another clade's
     (a contaminant, a misidentified genome);
   - `same_as`: the common ancestor of species with 99% identical copies;
   - `length_check`: its protein in the bottom or top 5% of its locus.
5. **Species without an assembly.** 75 species NCBI has no assembly of (EukDetect2 took them from EukProt, and
   re-assembled *Dientamoeba fragilis*' transcriptome) come from `data/no_assembly.fna.gz`. That is the one remnant of
   EukDetect2 (CC BY 4.0), made by `tools/no_assembly_records.py`. Their loci are called in those records as in a
   genome.

The 30 loci themselves are EukDetect2's choice among the BUSCOs it kept, after it dropped those found in gut bacteria.

## Layout

| path | contents |
|---|---|
| `eukdb/util.py` | compressed files, downloads, MD5s |
| `eukdb/ncbi.py` | NCBI's taxonomy and the subset that ships in a package |
| `eukdb/sources.py` | the loci table, NCBI's assembly tables, EUKARYOME and ROD |
| `eukdb/odb.py` | BUSCO's eukaryota_odb12: proteins, HMMs, cutoffs, a genome's nearest queries |
| `eukdb/loci.py` | the loci called in a genome (miniprot, hmmsearch), the copies and their features |
| `eukdb/qc.py` | the copies' checks across species |
| `eukdb/species.py` | the model and host species (`data/model_species.tsv`), their curated rDNA units |
| `eukdb/genomes.py` | the genome set, one worker per genome (loci, windows, rRNA), the prune |
| `eukdb/placement.py` | a genome's URL, contigs and pieces for miniprot |
| `eukdb/copies.py` | a copy's exons and introns (and protal's 16 kb fitting, which EUKDB itself does not apply) |
| `eukdb/rrna.py` | a genome's 18S and whole rDNA units (inline 5S kept), published records by NCBI taxid, the checks |
| `eukdb/package.py` | collect, write, read and check a package |
| `eukdb_build.py` | the build |
| `data/` | the loci table, the model species table, curated rDNA units, the no-assembly records |
| `tools/` | `no_assembly_records.py` (the no-assembly records from EukDetect2, run once) |
| `tests/` | `python3 -m unittest discover -s EUKDBbuild/tests -t EUKDBbuild/tests` |

## Building

Three stages, so that the network-bound one runs where the network is:

```
python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage fetch                 # sources and genomes (network)
python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage process -t 64 -j 16   # loci, rRNA; no network
python3 EUKDBbuild/eukdb_build.py --work WORK --out OUT --stage package --keep_gb 1000
```

`--stage all` (the default) runs the three in turn. The fetch stage gets these, then the genomes:
- NCBI's GenBank and RefSeq assembly tables of the `--groups`;
- NCBI's taxdump;
- BUSCO's eukaryota_odb12;
- EUKARYOME 2.1 and ROD v1.2.

The process stage calls every genome it was given. The package stage runs the checks across species, prunes
`OUT/genomes` to `--keep_gb` (the model species first, then the smallest), and writes `OUT/EUKDB_<version>/`
(`--version`, YYYY.MM by default). Each genome's results go into `WORK/done/<accession>/`, so a rerun continues where
the last one stopped.

`--groups` (default all six: fungi, protozoa, plant, invertebrate, vertebrate_mammalian, vertebrate_other) chooses
NCBI's assembly groups. At NCBI's tables of 2026-11-05 the groups hold:

| groups | species | genome size (one per species) | gzipped download |
|---|---|---|---|
| fungi and protozoa | 8,772 | 445 Gb | ~135 GB |
| all six | 29,114 | 21.6 Tb | ~6.5 TB |

So `--groups fungi,protozoa` first is a sensible start.

`--species 4930,7227,...` (NCBI taxids, or a file of one per line) builds a trial package of the species inside those
clades alone. The genome set, the published rRNA and the package are limited to them, and the manifest's
`inputs.partial` names the clades. The sources are fetched whole all the same.

A package is checked against its manifest by `eukdb.package.check(folder, deep=True)` (None when whole).

protal takes a package with `scripts/download_ref_data.py --eukdb OUT/EUKDB_<version>`.

## Not done yet
- The rRNA part screen (`records.tsv`'s `screen` column is `-`).
- NCBI's FCS-GX on the windows around the copies (needs a ~512 GB node); tiara classifies them where installed.
- A one-off check that the 30 loci are called in no bacterial genome (EukDetect2's screen, redone on EUKDB's caller).

## Requirements
- Python 3.10 or later; 3.14 reads zstd itself, older versions need the `zstd` command.
- miniprot and HMMER (hmmsearch).
- barrnap with nhmmer, ITSx with its whole ITSx_db (the build checks it), vsearch.
- tiara (optional).
