# Eukore's benchmarks and figures

This folder holds what a paper about Eukore needs and the database build does not:
- what describes a finished Eukore package as a resource (since 2026-11-12): its breadth at every rank, the markers
  one species links, named lists' organisms, its records' quality, the identity that parts taxa, the changes between
  releases;
- the comparisons of a finished Eukore package with other databases, as collections of records (EukDetect2,
  MetaPhlAn 4, PR2, and the references the build read);
- the read benchmark;
- every figure.

The rest of the checkout makes the database and its website. The build keeps only its own checks of the package
(completeness, OrthoDB, resolution, rRNA, names, composition, run time; [docs/benchmarks.md](../docs/benchmarks.md)),
whose headline figures the website shows. Until 2026-11-12 the build also ran the comparisons, the read benchmark and
the figures (`--read_compare`, `--read_genomes`, `eukore/figures.py`).

Everything here runs on a finished package, in its own environment, with its large downloads in a folder you choose.

## Three steps

```bash
bash benchmarks/eukore-bench-prep.sh BENCH OUT/Eukore_2026.11 [--genomes OUT/genomes] [--taxdump SOURCES/taxdump.tar.gz] [--odb SOURCES/eukaryota_odb12]
bash benchmarks/eukore-benchmark.sh  BENCH OUT/Eukore_2026.11 RESULTS [-t 32] [--tmp /local/disk] [--no_tree] [--align_cache WORK/align] [--previous OUT/Eukore_2026.05]
bash benchmarks/eukore-figures.sh    BENCH OUT/Eukore_2026.11 RESULTS [--format pdf]
python3 -m eukore_bench.site RESULTS --release OUT/release        # onto the website (after eukore-figures.sh; PYTHONPATH=benchmarks)
```

- `BENCH` is a folder for large files (about 9 GB of downloads, plus the read benchmark's test genomes where they are not linked from the build's folder): the downloads, those genomes, the environment.
- `RESULTS` is where the tables, `summary.json` / `summary.md` and `figures/` go. It is never the package itself,
  whose manifest checks its files.

| Script | What it does | Where it runs |
|---|---|---|
| `eukore-bench-prep.sh` | Makes the environment (first run), then every download, then the other references made comparable. Run it again to finish what failed: whole files are kept. | A node with internet |
| `eukore-benchmark.sh` | The comparisons and the read benchmark (`eukore_bench/compare.py`). The read benchmark aligns reads of up to 3 test genomes per supergroup to the whole package and to each other reference with minimap2: tens of GB of memory. | A compute node (threads: `-t`, else SLURM's `--cpus-per-task`) |
| `eukore-figures.sh` | Figure 1 and every supplementary figure, `index.md` (captions and tables) and `summary.md`, into `RESULTS/figures` (`eukore_bench/figures.py`) | Anywhere |

On a SLURM cluster, for example:
```bash
sbatch -p download -c 2 --mem 8G --wrap "bash benchmarks/eukore-bench-prep.sh $BENCH $OUT/Eukore_2026.11 --genomes $GEN"
sbatch -c 32 --mem 200G --wrap "bash benchmarks/eukore-benchmark.sh $BENCH $OUT/Eukore_2026.11 $RESULTS --tmp \$SLURM_LOCAL_SCRATCH"
```

## The environment

[`env.yaml`](env.yaml) holds only what these steps run: Python, numpy, matplotlib, minimap2, hmmer, fasttree, curl and isa-l. The
databases compared with are downloads, not packages. `bench_env.sh` makes the environment at `BENCH/env`, with its
package cache at `BENCH/pkgs`, using micromamba, mamba or conda, whichever is installed. Nothing goes to your home
folder's conda directories, and nothing needs activating: the scripts put `BENCH/env/bin` first on `PATH`. Delete
`BENCH/pkgs` once the environment exists, if space matters.

## The downloads (`BENCH`)

| Folder | What | Size | From |
|---|---|---|---|
| `downloads/taxdump.tar.gz` | NCBI's whole taxonomy. A package's `taxonomy/` holds only its own lineages, so species that another reference alone holds would be unknown to it. `--taxdump` gives the build's own instead (`SOURCES/taxdump.tar.gz`: the same version as the package) | 65 MB | NCBI |
| `downloads/pr2/` | PR2 5.1.1's 18S (`taxo_long`) | 60 MB | PR2's GitHub release |
| `downloads/metaphlan/` | MetaPhlAn 4's database archive, its MD5 (checked) and marker info; the server's `mpa_latest` (2026-11: `mpa_vJan26_CHOCOPhlAnSGB_202605`) unless `--mpa_version` | 6 GB | the bioBakery server |
| `downloads/eukdetect2/` | EukDetect2's database v1.1: `eukdb.fasta` and `busco_taxid_genome_link.txt` alone, MD5s checked. The rest of its 7 GB is bowtie2 indexes, which the comparison does not use | 2 GB | Zenodo record 23479662 |
| `downloads/eukaryota_odb12/` | BUSCO's dataset, for the tree's profiles (`hmms/`): the build's own linked (`--odb SOURCES/eukaryota_odb12`: the package's version), else BUSCO's newest downloaded | 60 MB; nothing where linked | BUSCO's server |
| `read_background/` | The read benchmark's 7 bacteria, an archaeon and human chromosome 21 ([`read_background.tsv`](read_background.tsv)) | 40 MB | NCBI |
| `read_genomes/` | The read benchmark's test genomes (`tests.tsv`): up to `--read_genomes` (default 3) per supergroup, those with the most congeners in the package first. Linked from `--genomes` (the build's folder) where they are there, else downloaded | Depends on the genomes chosen (a vertebrate's is ~1 GB); nothing where linked | NCBI |
| `compare/` | The other references made comparable (`eukore_bench/references.py`): `<name>.fna.gz` and `<name>_taxids.tsv` (PR2, MetaPhlAn 4's eukaryote markers), and `references.tsv` listing every reference with its FASTA and record-to-taxid table. EukDetect2's pair is used as it is | Smaller than the downloads | |

`--only taxdump,pr2,metaphlan,eukdetect2,odb,background,tests,convert` runs some of the steps.

## What is compared

`eukore-benchmark.sh` copies the package's own `benchmarks/` tables into `RESULTS`, then adds:
- **`breadth.tsv`**: the build's counts per supergroup (NCBI's species with an assembly, Eukore's), with the first
  other reference's columns. The order is EukDetect2, MetaPhlAn 4, PR2: the first that is there.
- **`overlap.tsv`**: the build's rows (PR2's names, EUKARYOME, KSGP, ROD, BUSCO's dataset) and one per other reference,
  each against every Eukore species (genome or rRNA), per supergroup.
- **`breadth_ranks.tsv`, `content.tsv`, `panels.tsv`**: the package as a resource (`eukore_bench/resource.py`; below).
- **`databases.tsv`, `upset.tsv`, `novel.tsv`, `agreement.tsv`, `agreement_records.tsv`**: the package beside the
  other databases (`eukore_bench/databases.py`; below).
- **`reads.tsv`, `detection.tsv`**: the read benchmark (`eukore_bench/reads.py`; below).
- **`summary.json`, `summary.md`**: the package's headline figures with these added.

### The package as a resource

Eukore is a database: a paper describes it by what it holds, before any use of it. Three tables (since 2026-11-12;
`python3 -m eukore_bench.resource --package PACKAGE --bench BENCH --out RESULTS` makes them alone, in a minute):

- **`breadth_ranks.tsv`**: per group (and `all`) and rank (`species`, `genus`, `family`, `order`, `class`, `phylum`):
  - `ncbi_named`: NCBI's named taxa. A species counts where its name is a binomial (no `sp.`, `cf.`, placeholder or
    hybrid); a higher taxon where it holds such a species. It needs NCBI's whole taxonomy (`-` without it).
  - `ncbi_genome`: the taxa NCBI holds a genome of (the package's `benchmarks/selection_taxa.tsv`; `-` for a package
    built before 2026-11-12).
  - `eukore_genome`, `eukore`: Eukore's, with a genome and with any record.

  A taxon above the species counts in every group one of its species is of (Chordata in Vertebrata and Metazoa), once
  in `all`. NCBI's named species stand for the known diversity: the species described but never sequenced (Catalogue
  of Life's) are not counted.
- **`content.tsv`**: per group and marker kind, the species holding it from any source (`species`) and from their
  genome's own assembly (`same_assembly`), of the group's species (`of`). The kinds: `BUSCO` (a BUSCO gene called),
  `18S`, `ITS` (ITS1, 5.8S and ITS2 on one record), `28S`, `operon` (a record of kind `rDNA`), `mitochondrion`,
  `plastid`, and three combinations (`BUSCO + 18S`, `BUSCO + operon`, `BUSCO + operon + organelle`). The
  `same_assembly` counts are the species whose genes and barcodes are linked through one assembly; a RefSeq organelle
  genome or a published read of the species is another source.
- **`panels.tsv`**: the organisms of named lists ([`panels.tsv`](panels.tsv) here: WHO's fungal priority pathogens,
  a selection of human parasites, model organisms, the mock communities' eukaryotes; edit it, or give another with
  `--panels`). Each entry is found in NCBI's taxonomy by its name or a synonym (or the taxid the list gives), and the
  table says what the package holds of it: `status` (`genome`, `records`: rRNA records alone, `absent`, `not found`),
  the package's species of it, NCBI's species of it with a genome, its BUSCO genes called, rRNA records and organelle
  genes. An entry above the species (*Fusarium*, Mucorales) is held where a species of it is.

### Beside the other databases

The other databases are compared as collections of records, not as tools (since 2026-11-12;
`python3 -m eukore_bench.databases --package PACKAGE --bench BENCH --out RESULTS` makes these alone):

- **`databases.tsv`**: the comparison table, one row per database. Its words come from
  [`databases.tsv`](databases.tsv) here (markers, taxa, what a record is, curation, updates, reference), written by
  hand: check each cell against the database's own paper before printing it. The benchmarks add what they count: the
  `version` read, its NCBI `species`, those in Eukore too (`both`), Eukore's alone, its alone, and its `records` and
  `bases` where its sequences are in `BENCH` (Eukore's: its copies whole, rRNA records and organelle genes).
  `figures/index.md` shows the table.
- **`upset.tsv`**: the species by the databases that hold them. Per group (and `all`) and combination of databases
  (`databases`: their names joined by ` + `; `n` of them; `eukore` `t` where Eukore is one), the species exactly
  those hold. The sets are the build's (`benchmarks/database_species.tsv` of the package, since 2026-11-12: PR2,
  EUKARYOME's long reads, KSGP, SILVA, UNITE, MIDORI2, ROD, BUSCO's dataset) and the other references of `BENCH`. A
  package built before has the other references alone.
- **`novel.tsv`**: the species a genome is the first sequence of a marker for. Per group and marker (`18S`, `ITS`,
  `operon`, `mitochondrial genes`): the species holding it from their genome's own assembly (`genome`), those of them
  a reference of that marker names (`named`), and the rest (`new`). `references` lists the references asked: PR2,
  SILVA, KSGP and EUKARYOME for the 18S; UNITE and EUKARYOME for the ITS; ROD and EUKARYOME for the operon; MIDORI2
  for the mitochondrial genes; each only where the build read it. A species is named by its name, a synonym or its
  taxid (`species.tsv`'s `<reference>_match`), so `new` is an upper bound: a reference may hold the species under
  another name. EUKARYOME's species are its long reads' alone.
- **`agreement.tsv`**: where Eukore and another database hold the same species' marker, whether the sequences agree.
  One row per database and marker: its `records` compared, those that `agree`, `disagree` (another species' gene),
  are `unlike` and `undecided`; the `species` compared, those with a record that agrees (`species_agree`) and those
  with none that agrees and one that disagrees or is unlike (`species_differ`); `how`. Three sources:
  - EUKARYOME's and ROD's reads: the build's own-species check against the genomes' records (`rrna/records.tsv`).
  - MIDORI2's records: the build's reference screen (the manifest's counts; `reference_screen.tsv` lists them).
  - Another 18S reference of `BENCH` (PR2): its records of the species Eukore holds a genome's 18S of, aligned to
    those 18S by minimap2 (`-x asm20`; an alignment counts where it covers 80% of the record or of Eukore's 18S). A
    record 97% or more alike to its species' own genome agrees. One as alike to another species' genome, and nearer to
    it by 2 points, disagrees. The other aligned ones are unlike.
- **`agreement_records.tsv`**: that reference's records that disagree or are unlike, each with its identity to its
  own species' genome, the nearest other species and the rank the two share. A disagreement does not say which side is
  wrong: the reference's record may be misnamed, or the genome may be another organism's (a contaminant's rDNA, a
  misidentified strain). These are the records, and the genomes, to look at first.

### How good the records are

Headline numbers of the package's own checks (since 2026-11-12; `eukore_bench/checks.py`, from the package alone:
`python3 -m eukore_bench.checks --package PACKAGE --out RESULTS`), and a tree (`eukore_bench/tree.py`):

- **`quality.tsv`**: per group (and `all`), level and check, the genomes or records `flagged`, `of` the level's in the
  group. The levels:
  - `genome`: the species' genomes (those of an assembly). Checks: each of NCBI's flags and of its FCS-GX and ANI
    reports' (`genomes.tsv`'s `flags`), `any flag`, `copies another clade's` (`foreign_copies`), a duplicated genome
    (`duplication: haplotypes`, `haplotigs`, `wgd_or_polyploid`), `genetic code differs from NCBI's`.
  - `BUSCO copy`, `rRNA record`, `organelle gene`: each flag of the record's `flags`, and its confidence tier
    (`confidence: high`, ...).
- **`gencode.tsv`**: per group, the genomes whose genetic code, checked against their BUSCO genes
  (`genomes.tsv`'s `gencode_check`), agrees with NCBI's table, differs, is undecided or was not checked, and the
  tables that fit those that differ (`best_tables`: `6:3; 10:1`).
- **`rdna_published.tsv`**: the genomes' rDNA copies by read depth (`genomes.tsv`'s `rdna_copies_reads`, where
  `tools/rdna_depth.py` measured them) beside published numbers ([`rdna_published.tsv`](rdna_published.tsv) here:
  the five species of docs/read_depth.md's trial; add others with their source), and whether they lie `within` the
  published range.
- **The excluded genomes** are the package's own `benchmarks/selection_taxa.tsv` (since 2026-11-12): every NCBI taxon
  with an assembly, what became of it, NCBI's assemblies of it and the flags that excluded them. It is copied into
  `RESULTS` with the package's other tables.

**The tree** (`tree.nwk`, `tree_tips.tsv`, `tree_taxa.tsv`, `tree_strays.tsv`): whether the BUSCO proteins place each
genome where its name says. It runs where `BENCH` holds BUSCO's dataset (`eukore-bench-prep.sh`'s `odb` step; give
the build's with `--odb SOURCES/eukaryota_odb12`, the package's own version) and the environment has hmmalign and
FastTree; `--no_tree` leaves it out, `python3 -m eukore_bench.tree` runs it alone.

1. **The genomes.** Up to 2 species' genomes per genus (`--per_genus`), those with the most BUSCOs Single first: a
   genus is tested by two of its species.
2. **The loci.** The 50 BUSCOs (`--loci`) Single, and not of confidence low, in the most of those genomes. A genome
   with under half of them is left out.
3. **The alignment.** Each locus' proteins (`busco/proteins/`) are placed on BUSCO's profile by hmmalign
   (`eukore/align.py`, as the build places them): a residue on a match state is that column's. No alignment of the
   proteins to each other is made. `--align_cache WORK/align` reads the build's own rows. The loci's columns are
   joined; a column is kept where half of the genomes hold a residue.
4. **The tree.** FastTree (`-lg -gamma`).
5. **The taxa** (`tree_taxa.tsv`). Each NCBI taxon of genus to phylum with two genomes or more in the tree (and not
   all of them): `monophyletic` `t` where its genomes are one branch. The tree is unrooted, so a branch is one side
   of an edge. Else the branch that fits it best (`f1`: of the branch's genomes and the taxon's), the taxon's genomes
   outside it (`strays`) and other taxa's inside (`intruders`).
6. **The strays** (`tree_strays.tsv`). A genome outside its taxon's best branch, where the taxon has three genomes or
   more and a fifth of them at most (one at least) lie outside: the highest rank it strays at (`rank`, `taxon`; all
   of them in `ranks`) and the genomes of the branch it sits on (`beside`). These are the genomes to look at first: a
   misnamed genome, a contaminated one, or a taxon NCBI holds together and the genes do not.

The tree measures the genes' agreement with NCBI's taxonomy, not the true phylogeny: 50 proteins on profile columns
and FastTree resolve families and orders, not the deep branches of the eukaryotes.

### Guidance: the identity that parts taxa

A user of the package chooses a marker and an identity cutoff. **`thresholds.tsv`** (since 2026-11-12;
`eukore_bench/guidance.py`) gives, per marker (all BUSCOs together, each BUSCO, the 18S) and rank (`species`,
`genus`, `family`, `order`), the identity that best parts the sequences that should stay together (`within`) from
those that should be apart (`between`):

- A species: `within` is the intraspecific panel's species, each by the least alike pair of its own assemblies
  (`intraspecific/barcode_gap.tsv`'s `intra_min`). `between` is every copy's nearest copy of another species.
- A genus, family or order: `within` is the copies whose nearest other species shares that rank and no deeper one.
  `between` is those whose nearest shares only a higher rank, or none.

`threshold` is the cutoff with the best balanced accuracy (the highest of those as good), `within_kept` and
`between_apart` the shares on its right side, `accuracy` their mean. A marker and rank with under 5 sequences on a
side has no row. A threshold of 70 is the search's floor: any hit at all is of that rank. The identities come from
the package's `benchmarks/resolution_hist.tsv` (a package built since 2026-11-12, with vsearch).

The rows of `source` `published` are thresholds from the literature
([`thresholds_published.tsv`](thresholds_published.tsv): fungal ITS and 28S D1/D2, Vu et al. 2019), for context:
they are of other markers and taxa. Read the package's thresholds with its sampling in mind: a clade with few genomes
has few close relatives, and its thresholds come out low.

Short worked uses of a package (choosing loci for a clade, checking a primer pair, a tree, a mapping reference) are
in [docs/use_cases.md](../docs/use_cases.md).

### The release among releases

What a reader needs to reuse a release and to follow it to the next (since 2026-11-12; `eukore_bench/history.py`;
`python3 -m eukore_bench.history --package PACKAGE --out RESULTS [--previous OLDER_PACKAGE]` makes these alone):

- **`sizes.tsv`**: the package's parts (its folders; `busco/` and `markers/` by kind of file), their files and bytes,
  as the manifest lists them; `all` last. `summary.json`'s `cost` puts the build's time beside: the package stage's
  steps summed (`stages.tsv`) and the genomes' runs summed (`runtime.tsv`, on the threads each ran on). Peak memory
  is not recorded by the build: take it from the scheduler's accounting.
- **`growth.tsv`**: per group (and `all`) and year, the species' genomes NCBI released that year
  (`metadata/ncbi.tsv.gz`'s `seq_rel_date`). The figure sets the other databases' release dates beside
  (`databases.tsv`'s `released` column here, where filled): what a database frozen then could not hold.
- **`changes.tsv`, `changes_species.tsv`** (with `--previous`, an earlier package's folder): per group the species of
  the earlier package (`previous`) and of this one, those `added`, `removed` and `kept`; of the kept, the genomes
  `replaced` by another assembly or `updated` to a newer version of the same, and the species that gained or lost a
  genome; and those species one by one, with both genomes. This is the table a release note is written from.

Not made by code: a DOI per release. Deposit `OUT/release` (the build's release step) at Zenodo or a like archive
for each version, and cite that DOI with the version.

### The read benchmark

From each test genome, reads are simulated from its markers' loci (introns and all) and from anywhere in it, plus
reads from the background genomes. They are aligned by minimap2 to:
- Eukore;
- Eukore's exons alone;
- Eukore's 18S alone;
- each other reference.

A read's hits within 1 point of its best (90% identical at least, covering 80% of the read) name their lowest common
ancestor. No profiler is involved: the benchmark measures the reference. Each read is assigned three times:
- `in`: its species in the reference;
- `held_out`: its species left out;
- `novel`: its whole order left out.

Detection (`detection.tsv`): each test genome is a sample of 1,000, 5,000 and 20,000 of its reads, with as many
bacterial reads. A taxon is detected with 3 reads or more and 1% or more of the reads placed at its rank. The F1 per
reference is computed for the right species (in) and genus (held out). Human chromosome 21 in the background stands
for a gut or skin metagenome's host: its reads are right on its own lineage and wrong elsewhere.

## The figures

`RESULTS/figures/index.md` captions each figure and names its tables. Figure 1 (`Fig1.png`, at a journal's
two-column width) asks four questions:

| Panel | The question | What it shows | Its tables | Needs |
|---|---|---|---|---|
| a | How large? | Species per supergroup: NCBI's with a genome, Eukore's (with a genome; any record), each other reference's | `breadth.tsv`, `overlap.tsv` | The package; the other references from this folder |
| b | How far from BUSCO's dataset does a call stay right? | Left: the dataset's genomes called again without their genus', family's and order's dataset species. Right: each genome's calls laid on its own annotation, per annotation source | `holdout.tsv`, `annotation.tsv` | The package (the build's own tables) |
| c | Why genes, not the 18S? | Per marker kind: the copies under 99% identical to their nearest congener; the intraspecific panel's barcode gap | `resolution.tsv`, `intraspecific/barcode_gap.tsv` | The package |
| d | Does it place a metagenome's reads? | F1 of detecting a test genome's species and genus by the sample's depth, per reference | `detection.tsv`, `reads.tsv` | `eukore-benchmark.sh` |

A panel whose tables are missing says so in its place. The supplementary figures follow, one per benchmark. Drawn
from a package's own `benchmarks/` alone (`python3 -m eukore_bench.figures PACKAGE --benchmarks
PACKAGE/benchmarks --out DIR`), the figures that need no other reference are complete.

## What is here

| Path | |
|---|---|
| `eukore-bench-prep.sh`, `eukore-benchmark.sh`, `eukore-figures.sh` | The three steps |
| `bench_env.sh`, `env.yaml` | The environment, made at `BENCH/env` |
| `eukore_bench/` | The Python of the three steps: `prep.py`, `compare.py`, `resource.py`, `databases.py`, `checks.py`, `tree.py`, `guidance.py`, `history.py`, `reads.py`, `references.py`, `figures.py`, `site.py` (the results for the website). It reads packages with the database's own code (`eukore/` of this checkout) |
| `read_background.tsv` | The read benchmark's background genomes |
| `panels.tsv` | The named lists of organisms `panels.tsv` of the results looks up |
| `databases.tsv` | What each database compared is, in words: the comparison table's text |
| `rdna_published.tsv` | Published rDNA copy numbers, beside which the copies by read depth are set |
| `thresholds_published.tsv` | Published identity thresholds per rank, shown beside the package's own |
| `paper/` | The paper's own scripts and figures (not in git) |
