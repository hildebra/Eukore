# Building Eukore

How to install the build, run it (a trial, or all eukaryotes on a cluster), resume it, and read what it leaves.
The package it writes is described in [FORMAT.md](FORMAT.md), its methods in [methods.md](methods.md), its benchmarks in
[benchmarks.md](benchmarks.md).

## Install

Eukore runs on Linux (or WSL on Windows): its tools are bioconda's. It needs a conda-style package manager (micromamba,
mamba or conda), and nothing else from the system but `git`.

### 1. The code

```bash
git clone <this repository> Eukore      # or copy the folder: Eukore needs no build step
cd Eukore
```

Run everything from this folder (`eukore_build.py` finds its `eukore/` package and `data/` beside it).

### 2. The build environment

```bash
micromamba create -f envs/eukore.yaml
micromamba activate eukore
```

`envs/eukore.yaml` holds:
- Python 3.12 or later (the build stops at its start on a Python without `tarfile`'s extraction filter);
- miniprot and HMMER (hmmsearch, nhmmer, hmmalign; its easel tools' esl-translate for the rescue pass): the BUSCO calls,
  the introns' places and the genetic code check (hmmalign), and the
  organelle genes;
- barrnap and ITSx: the rRNA genes and the ITS borders;
- vsearch: the checks across species, the part screen, the duplicated copies' identity and the barcode gap;
- Infernal (cmalign): the 18S variable regions (since 0.2; without it, none);
- ISA-L's igzip: the downloads' gzip checks (else pigz, else Python's, ~3.5x slower);
- p7zip: EUKARYOME's archive;
- openpyxl: EUKARYOME's tables;
- ncbi-fcs-gx: optional, the FCS-GX screen;
- matplotlib: the benchmarks' figures, drawn by the package stage (`eukore/figures.py`; without it, left out).

Check that the tools are found:

```bash
for t in miniprot hmmsearch hmmalign esl-translate nhmmer barrnap ITSx vsearch cmalign 7za curl igzip; do printf '%-10s %s\n' $t "$(command -v $t || echo MISSING)"; done
python3 eukore_build.py --help
```

The build stops before any work if miniprot or hmmsearch is missing, or if ITSx's profile database is incomplete (a
damaged ITSx_db gives silently wrong borders). Without barrnap it makes no rRNA records, and without ITSx it uses
barrnap's borders alone: it warns either way.

### 3. Optional: tiara

tiara classifies the 50 kb around each BUSCO copy (eukaryote, prokaryote, organelle: `busco/copies.tsv`'s
`window_class`). It needs a Python of 3.9 or older, so it gets its own environment:

```bash
micromamba create -f envs/tiara.yaml
```

The build is then given it with `--tiara $MAMBA_ROOT_PREFIX/envs/eukore-tiara/bin/tiara` (or finds `tiara` on the
PATH). Without it, `window_class` is `-`, and `--no_tiara` says so explicitly.

### 4. Optional: FCS-GX's database

NCBI's FCS-GX screens the windows around the copies for other organisms' sequence (`busco/copies.tsv`'s `fcs`). Its
database is ~470 GB and wants a node of ~512 GB of memory, best read from `/dev/shm`:

```bash
sync_files.py get --mft https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir /data/gxdb
cp -r /data/gxdb /dev/shm/gxdb      # on the node that runs the package stage
```

Then give the package stage `--fcs_gx_db /dev/shm/gxdb/all`. NCBI's FCS-GX wiki describes the database and its
download. Without it, `fcs` is `-`. NCBI's own FCS-GX summaries of its assemblies, which rank the genomes, are fetched
either way.

### 5. Optional: protal's KSGP

`species.tsv` gives each species' KSGP lineage. The fetch stage takes KSGP 4.0, the newest the KSGP site serves. To
label species as a protal database does, give the build that KSGP release's tax file:
`--ksgp_tax /path/to/ksgp_4.4/KSGP_v4.4_plus2.tax`.

### 6. Check the install

```bash
python3 -m unittest discover -s tests -t tests      # ~80 tests with stand-in tools, a few seconds
```

### What it takes

| | |
|---|---|
| sources (`WORK/sources`) | ~4.4 GB: NCBI's assembly tables, taxdump, FCS-GX and ANI reports; BUSCO's eukaryota_odb12; EUKARYOME 2.1; ROD v1.2; RefSeq's organelle genomes (~2.2 GB); PR2; KSGP (243 MB); since 0.2 Rfam's SSU model (1 MB), SILVA 138.2 SSU NR99 (201 MB), UNITE 10.0 (50 MB), Open Tree Taxonomy 3.7.3 (111 MB), ENA's list of TSA sets (2 MB) |
| genomes, fungi and protozoa | 8,772 species, ~135 GB gzipped (445 Gb of sequence) |
| genomes, all six groups | 29,114 species, ~6.5 TB gzipped (21.6 Tb); fetched and pruned as the build goes |
| transcriptomes (since 0.2) | ~5,000 species without a genome, 341 Gb of transcripts (ENA's TSA sets, 2026-11); `--no_tsa` skips them |
| the intraspecific panel (since 0.2) | up to 3 further assemblies of 40 species (2 of the vertebrates, rice and maize): ~110, deleted once run |
| memory | the package stage: ~2 GB for a trial's; for all eukaryotes several GB more (the copies' rows, coding sequences, proteins and introns; since 2026-11-08 the copies whole and their features are read again from `WORK/done` as they are written, not held), measure on a first run; per genome job, miniprot's index on that genome (~7 bytes a base, in pieces of 4 Gb at most: a few GB, ~22 GB for a mammal's), the rest of the run holding only the copies' stretches (since 2026-11-08; until then a chromosome-level genome whole), ~3 GB besides; the process stage begins each genome where its estimate fits in `--memory_gb` (since 0.2.1), so that a node of 200 GB runs ~6 mammals at once, ~20 threads each, and dozens of fungi (simulated on the 29,114 species' genome sizes with 80 cores: 384–512 GB keeps the cores busy, 256 GB takes 8–18% longer, 200 GB 13–28%, 128 GB 25–57%: two thirds of the bases are in the 7,088 genomes over 1 Gb, whose memory then sets how many run); cmalign up to 2 GB per thread for the 18S regions; the fetch stage ~2 GB (the taxonomy) plus, per download of a TOGA2 `.2bit` genome, its longest chromosome (~250 MB for a mammal's; until 2026-11-08 ~6 times that: 12 downloads at once killed a build on a 10 GB node) |
| network | the fetch stage alone |

## Usage

### A first, trial build

A trial builds the package of a few clades on one machine, all stages in one run:

```bash
python3 eukore_build.py --build_dir work --out out --groups fungi --species 4930 -t 8
```

`--species` takes NCBI taxids (here *Saccharomyces*, 4930), comma-separated or a file of one per line. Only the species
inside those clades are built, and the manifest says so (`inputs.partial`): a trial package is not Eukore. The
2026-11-05 trial of seven clades (*Saccharomyces*, *Candida albicans*, *Blastocystis*, *Drosophila melanogaster*,
*C. elegans*, *Arabidopsis*, *Trichuris trichiura*; 42 genomes, 1 Gb) took on 6 threads:
- ~9 minutes for the genomes;
- ~21 minutes to annotate the published rRNA reads the first time (kept in `WORK/published` for later runs);
- ~7 minutes for the package stage after that.

### A full build, in three stages

On a cluster, the stage that needs the network runs where the network is, and the rest on compute nodes:

```bash
# a node with internet: the sources, then every species' genome
python3 eukore_build.py --build_dir WORK --out OUT --stage fetch --groups fungi,protozoa

# compute nodes, no network: each genome's BUSCOs and rRNA, then the organelle genes
python3 eukore_build.py --build_dir WORK --out OUT --stage process --groups fungi,protozoa -t 64 -j 16

# then the package: the checks across species, the published rRNA, the prune, the package, the release
python3 eukore_build.py --build_dir WORK --out OUT --stage package --groups fungi,protozoa -t 64 --keep_gb 1000 \
    [--fcs_gx_db /dev/shm/gxdb/all] [--tiara .../bin/tiara] [--ksgp_tax .../KSGP_v4.4_plus2.tax]

# only the release again (the website's data, the downloads) of a finished package
python3 eukore_build.py --build_dir WORK --out OUT --stage release --version 2026.11 -t 64
```

Give every stage the same `--build_dir`, `--out`, `--genomes`, `--sources`, `--groups` and `--species`: each stage plans the same
genome set from them. `--build_dir` is the build's own state, kept for every rerun (the sources, each genome's
results, the caches); `--out` what it makes (the package, the release); `--local_tmp` a node's throwaway temporary
files. Until 2026-11-09 `--build_dir` was `--work` and `--local_tmp` `--scratch`: both old names are still taken, as
they were. The genomes go into `OUT/genomes` unless `--genomes DIR` puts them in a folder of their own (since
2026-11-09), e.g. a download area apart from the package; to move a build's genomes there, move the folder's files
(the genomes, their gene models, assembly reports and `.annotation_none` notes, `manifest.tsv`) with it.
`--stage all` (the default) runs the three in turn on one machine with internet. Stages that follow one another may
be given together (since 2026-11-10): `--stage process,package` is one run with one plan, the passes after the genomes
made once (`fetch,process` likewise).

The process stage may run on several nodes at once, e.g. a SLURM array: each run claims the genomes it takes, so the
runs share them out without a split. With `--genomes_only` (since 2026-11-09) each runs its genomes, then the passes
that take one done genome at a time under its claim (since 2026-11-10: the backfill, whose clade hold-outs are whole
miniprot calls, the assembly metadata, the packing of older done folders; until then the last run made them alone, on
one node), and stops; one run without it, after them all, makes what they left and the passes over all genomes (the
intraspecific panel, the organelle genes, the gene models), then the package. No task list is written: the claims
are the queue, and a killed task's genomes and steps fall to the others. After the download, that is one chain of
two jobs, submitted together (since 2026-11-10):

```bash
B=(python3 eukore_build.py --build_dir WORK --out OUT --groups fungi,protozoa --version 2026.11)
# process.sh, 8 at once:  "${B[@]}" --stage process --genomes_only --local_tmp $TMP -t $SLURM_CPUS_PER_TASK
# final.sh, after them:   "${B[@]}" --stage process,package --local_tmp $TMP -t $SLURM_CPUS_PER_TASK --keep_gb 1000
jid=$(sbatch --parsable --array=1-8 process.sh)
sbatch --dependency=afterany:$jid final.sh
```

- `afterany`, not `afterok`: the final job starts once every process job has ended, however it ended. A task that ran
  into its time limit or lost a worker (the OOM killer) would leave an `afterok` job pending for good; the final run
  needs none of them to have ended well, since it runs the genomes they left itself. A genome whose claim is still
  live when it comes to it (a task killed less than 5 minutes ago, a process run still going) it waits for, then runs
  where that run left it, before the passes (since 2026-11-10: it went on without such a genome, and the package was
  written without it).
- A `process,package` run goes on to the package past genomes that are not fetched: it has no network to mend them,
  and a full build always has some assemblies NCBI no longer serves. It names them; `--stage fetch`, then the chain
  again, adds them. (`--stage process` alone still stops there with exit 0, `all` with 1.)
- Give `--version` on a cluster: its default is the month of the run's day, and a chain that crosses a month's end
  would otherwise write a second `Eukore_<version>`.
- `--local_tmp` serves the package stage too (since 2026-11-10): its tools' temporary files go to the node's disk.

`tools/eukore.sh` does this for a SLURM cluster: `bash tools/eukore.sh fetch` submits the download job,
`bash tools/eukore.sh build [N]` the chain (N process jobs, default 8, then the final job;
`AFTER=<job id>` makes the chain wait for the fetch job, and a failed fetch cancels it; `REVERSION=1` adds
`--reversion`). Its settings (the folders, the partitions, the build's options) are a small bash file,
`tools/eukore.conf.example`, copied to `local/eukore.conf` or named by `$EUKORE_CONF`.

Each run opens with a banner on STDOUT (since 2026-11-08): Eukore's version, the date, host and Python, the command
(a `--veupathdb_key` on it hidden), then what each of its stages will do and with what: the groups and species, the
loci, the folders, the threads, the process stage's memory, the package stage's `--keep_gb`, and what is switched off.

- **fetch** gets the sources into `WORK/sources` (or the folder `--sources` names), each kept once there and whole:
  - NCBI's GenBank and RefSeq assembly tables of the `--groups`;
  - NCBI's taxdump;
  - NCBI's FCS-GX summaries of its assemblies and its Saccharomycotina ANI report (~70 MB);
  - BUSCO's eukaryota_odb12;
  - EUKARYOME 2.1 and ROD v1.2;
  - RefSeq's organelle genomes (mitochondrion and plastid);
  - PR2 5.1.1's taxonomy (5 MB) and KSGP 4.0's (243 MB; not with `--ksgp_tax`);
  - since 0.2: SILVA 138.2's SSU Ref NR99, UNITE 10.0's release of all eukaryotes and Open Tree Taxonomy 3.7.3
    (`WORK/sources/labels`), Rfam's eukaryote SSU model (`WORK/sources/rfam`), ENA's list of eukaryote TSA sets
    (`WORK/sources/tsa`). MIDORI2 and BOLD are not fetched (their URLs change, BOLD's bulk data needs a login): give
    them with `--midori2` and `--bold_tsv`.
  - since 0.2, the inputs of Eukore's marker set (`WORK/sources/eukore_markers`, `eukore/markers.py`): KOfam's
    `ko_list` and its markers' and decoys' profiles (out of KOfam's 1.6 GB archive, fetched in 16 ranges at once and
    deleted), Swiss-Prot's eukaryote proteins (~90 MB, EBI's FTP) and NCBI's proteins named as each marker
    (`data/markers.tsv`'s Entrez query: RefSeq's of every eukaryote, GenBank's of the protists; 3 per organism, by
    E-utilities), and those of each decoy (`data/marker_decoys.tsv`). Most of an hour, NCBI's. The process stage
    builds the set of them (below; since 2026-11-11 the fetch stage did, its hmmsearch on the download node).

  The fetch stage downloads and unpacks; it computes nothing, so a cluster's download node (`tools/eukore.sh fetch`)
  needs few cores.

  It also fetches the transcriptomes of species without a genome and the intraspecific panel's further assemblies,
  as it fetches genomes.

  Then it fetches each species' genome into `OUT/genomes` (its RefSeq copy where it has one). Every minute for the
  first 10 minutes, then every 10 (`--report_minutes`; until 2026-11-09 the first after 1, 2, 3, 5 and 8 minutes), it
  says the genomes fetched, failed and left, the Gb left, the download rate and the time
  left, estimated by the genomes' sizes; the process stage reports its genomes alike. Beside each genome it fetches
  its gene models and NCBI's assembly report: NCBI's GFF (RefSeq's or a GenBank submitter's; not asked for where the
  assembly summary gives no protein-coding genes; a 404 or a GFF without coding sequences is none), else VEuPathDB's (with
  `--veupathdb_key`), else Ensembl's, else TOGA2's (their indexes fetched once into `SRC/annotation_sources`); the process stage takes a BUSCO copy's
  model from it where the HMM confirms it ([methods.md](methods.md)), compares miniprot's models with it
  (`annotation.tsv`) and deletes it.
- **process** first builds Eukore's marker set where it is not built from this code's tables, rules and inputs:
  the candidates scored with the KOfam profiles (hmmsearch, on `-t` threads) into the markers' query proteins, one
  per species, and the decoys' profiles (hmmalign, hmmbuild); the scores are kept for a rebuild; what it was made of is the package manifest's `inputs.markers`. One process run builds it at a time (a
  lock in its folder): a SLURM array's other tasks wait for it. An input missing (the fetch stage never fetched it)
  is said, and the genomes are called with the set of before, or without Eukore's markers. Then it
  calls each fetched genome's BUSCOs and extracts its rRNA into `WORK/done/<accession>.zip` (one file a
  genome since 2026-11-10; `--done_folders`: a folder of 8-14 small files, as before), `-j` genomes at a
  time with `-t / -j` threads each, the largest first (since 2026-11-08; the rRNA and the assembly metadata run beside
  miniprot). A genome not fetched yet is left for the next run; so is one whose annotation the fetch stage has not
  asked for yet (it fetches a genome's annotation before the genome itself). Several process runs may share one WORK,
  e.g. one per node (since 2026-11-08): each genome is claimed by the run that takes it (`WORK/done/.<accession>.claim`;
  touched every minute, a claim of a killed job or a node that died is taken over after 5 minutes untouched, so a
  restart runs the genomes its killed job held), a genome another run finished after this one planned is left alone
  (since 0.2.2: run again), and each run then makes the passes below too. Each run's plan, and its backfill and
  metadata passes, read the done genomes 16 at a time (since 2026-11-11: one after another, 5-10 calls of the shared
  file system each; with some 30,000 done genomes and 32 array tasks planning at once, half an hour with the cores
  idle before the first genome began), and a SLURM array's tasks share one plan of the genomes to run (since
  2026-11-11): the first task makes it under a lock (`WORK/plan/.<array job>.genomes.lock`) and writes
  `WORK/plan/<array job>.genomes.json`, the others wait for that file and read it (a killed planner's lock is taken
  over after 5 minutes, by a task that plans in its place; a plan that could not be written, e.g. on a full WORK,
  has each task make its own); with `--genomes_only` those that take a
  genome at a time under its claim (the backfill, the metadata), the others left to one run after them all (since
  2026-11-09; above). A failed rRNA extraction no longer fails its genome (since 2026-11-10: its calls were
  thrown away and made again by every run): the genome is done, `done.json` notes `rrna_failed`, and the backfill
  tries the rRNA again. A genome whose run fails has its reason noted in `WORK/failed/<accession>.txt` (with the
  time, host and SLURM job) until a later run does it, and each run's process stage ends by listing those of every
  run on the WORK in `WORK/failed.tsv` and its own failures by reason in its log (since 2026-11-11: the log named
  three of them; a full WORK, which notes nothing, still gets its reasons counted in the log). A dataset species' genome, called again for each hold-out, has miniprot's index built
  once in `--local_tmp` and read by the later calls (since 2026-11-10; genomes up to 4 Gb, some 4-8 bytes a
  base on the node's disk while the genome runs). A worker process that dies (the OOM killer,
  a crash) stops the stage with a message; a rerun goes on (since 0.2.2: the stage waited for it forever). A genome
  over 4 Gb is aligned in pieces of 4 Gb at most, a chromosome over that in windows of 4 Gb that reach 5 Mb
  into each other (since 2026-11-10; until then whole, alone in its piece, over the memory the scheduler reserved;
  since 0.2.2 a piece no longer takes one more whole chromosome past it). The rescue's gate and the intron bounds take
  NCBI's total length where the plan has it (since 2026-11-10: the gzip file's size over 0.3, up to 7% off, far more
  for an N-rich genome). Then it
  makes what a tool missing at a genome's run left out (tiara's classes, the rRNA, a clade hold-out: a backfill; since
  2026-11-10 begun as the genomes' memory fits, as the genome runs are: a hold-out is a miniprot call), takes the
  organelle genes (`WORK/organelle`; from the organelle molecules the metadata scan kept, not the genome again), and runs the
  intraspecific panel's further assemblies (their BUSCOs, rRNA and organelle genes; since 0.2.1 run in the same pool
  as the plan's genomes), deleting their genome files once done. Since 0.2.1 the genomes are begun as their memory fits
  (`--memory_gb`; a genome taking a larger share of it than `-t / -j` threads of `-t` gets that share of the threads,
  up to 32), each genome's hmmsearch runs are split by model over its threads, the rescue's six-frame ORFs are made in
  as many pieces at once, and the organelle and metadata passes run a genome a core. Each of these passes reports on
  the same clock as the fetch stage (every minute of its first 10, then every `--report_minutes`), whether or not an
  item finished meanwhile (since 0.2.1). Each report gives the items done,
  failed, running and left, the time so far and the time left, and the resources the build uses: the cores its
  processes kept busy since the last report (out of `-t`), their resident memory and its peak, the node's free memory
  and load, and the free disk of WORK and OUT. The stage's first and last lines give the same, the last with its time:

  ```
  [14:02:11]   12/3,410 genomes run (0 failed), 3,398 left (1,204.5 Gb); 2h 00m so far, ~3d 04h left (until ~2026-11-11 18:30), 8 running; 30.2 cores busy (of 32), 84.1 GB resident (peak 120.3 GB); node: 300.0 GB of 512.0 GB free, load 31.5; disk free: WORK, OUT 2.1 TB
  ```
- **package**, in turn:
  - makes what a tool missing at a genome's run left out, as the process stage does (tiara given here: its classes);
  - runs FCS-GX where asked (`-t / 4` genomes at once, 4 cores each);
  - checks the copies across species;
  - since 0.2: aligns the copies' proteins to BUSCO's HMMs (hmmalign, `-t` loci at once; cached in `WORK/align`), places their introns,
    checks each genome's genetic code;
  - annotates the published rRNA reads (barrnap and ITSx over `-t` threads, 250 reads a batch);
  - runs the rRNA checks and the part screen;
  - since 0.2: finds the 18S variable regions (cmalign; cached in `WORK/cmalign` batch by batch, an 18S cmalign
    failed on tried again by the next run) and runs the in-silico PCR;
  - takes the organelle genes of any genome not searched yet, checks them across species (since 0.2);
  - since 0.2: measures the intraspecific panel's barcode gap;
  - names each species in PR2, EUKARYOME, KSGP and (since 0.2) SILVA, UNITE and MIDORI2, and finds its ids in Open
    Tree's taxonomy (`crossrefs.tsv`);
  - describes each done genome still on disk that lacks its metadata (`eukore/metadata.py`: a genome done before
    2026-11-07), as the process stage does after its runs;
  - prunes `OUT/genomes` to `--keep_gb` where given (the model species first, then the smallest); without it, no
    genome file is deleted (since 2026-11-08);
  - writes `OUT/Eukore_<version>/`, with its benchmarks ([benchmarks.md](benchmarks.md)).

### Options

| option | default | |
|---|---|---|
| `--build_dir WORK` | required | the build's own folder, kept for the whole build and every rerun: the sources (`WORK/sources`, unless `--sources`), each genome's results (`WORK/done/<accession>/`: a full build's weeks of runs), the package stage's caches; not a scratch folder (`--work` before 2026-11-09, still taken) |
| `--out DIR` | required | where `Eukore_<version>/`, `release/` and `genomes/` (unless `--genomes`) go |
| `--genomes DIR` | `OUT/genomes` | the genomes' folder: the fetch stage downloads each species' genome, its gene models and NCBI's assembly report into it, the process stage reads them from it, the package stage cuts it back to `--keep_gb`; the same for every stage (since 2026-11-09) |
| `--sources DIR` | `WORK/sources` | the sources' folder: the fetch stage downloads into it what every other stage reads (NCBI's taxonomy and assembly tables, BUSCO's dataset, the marker set, the references' files); a folder of its own keeps them on a file system that keeps its files while WORK is elsewhere, e.g. on scratch. Builds may share it: a fetch takes its lock there too, one at a time. A moved folder serves as it was. A stage other than fetch stops at once without it (its `taxdump.tar.gz`); the same for every stage (since 2026-11-11) |
| `--stage` | `all` | `fetch`, `process`, `package` or `all` (the three in turn); stages that follow one another, with commas: `process,package`, `fetch,process` (one run, one plan; since 2026-11-10); `release`: `OUT/release/` again of a finished `OUT/Eukore_<version>/` (since 2026-11-09) |
| `--genomes_only` | | with `--stage process`: run the genomes, then the passes that take one genome at a time under its claim (since 2026-11-10: the backfill, the assembly metadata, the done folders' packing), and stop; the passes over all genomes (intraspecific panel, organelle genes, gene models) left to a run without it after them all (`--stage process,package`: it waits for the genomes other runs still hold, runs what they left, then makes the passes and the package); for several runs on one WORK, e.g. a SLURM array ([A full build](#a-full-build-in-three-stages), `tools/eukore.sh`; since 2026-11-09) |
| `--done_folders` | | keep each genome's results as a folder of 8-14 small files (`WORK/done/<accession>/`, as every build before 2026-11-10) rather than one file (`WORK/done/<accession>.zip`, the default: for all eukaryotes some 30,000 files in place of 300,000 in 30,000 folders on the shared file system). Either way a WORK holding the other form is read as it is; without this flag the process stage packs the folders it finds |
| `--no_release` | | the package stage writes no `OUT/release/` (the website's data and the downloads, a second copy of the package: since 2026-11-09 written last by default) |
| `--groups` | all six | NCBI's assembly groups: `fungi`, `protozoa`, `plant`, `invertebrate`, `vertebrate_mammalian`, `vertebrate_other` (comma-separated) |
| `--species` | | a trial: only the species inside these NCBI taxids' clades (comma-separated, or a file of one per line) |
| `--loci FILE` | every BUSCO | call only these BUSCOs (one per line, or a loci table's first column): a quicker trial |
| `--version` | `YYYY.MM` of today | the package's version, in its folder's name (Eukore's own version, 0.2.4, is in the manifest's `eukore`) |
| `-V` | | print Eukore's own version and exit |
| `--keep_gb` | none: no genome deleted | the genomes kept in `OUT/genomes`, GB: cut back to it before the package is written (its `genomes.tsv` says which are kept); with it, a superseded assembly's and the intraspecific panel's genome files are deleted once done too (until 2026-11-08 the default was 1000) |
| `-t`, `--threads` | 8 | threads in all; since 2026-11-09 also those that compress every `.gz` file the build writes (BGZF: the package's and the release's at level 9, the work folder's at 6) |
| `-j`, `--jobs` | threads / 4 | genomes at a time at most, `-t / -j` threads each; fewer where `--memory_gb` holds fewer, each then with more threads (since 0.2.1) |
| `--memory_gb` | SLURM's allocation, else the cgroup's limit, else the node's memory available; 85% of it | the process stage's memory: a genome is begun where its run's estimate (~7.5 bytes a base, of 4 Gb at most, and 3 GB) fits in what the running ones leave, the largest first, a smaller one in the room a larger one does not fit (since 0.2.1; until then `-j` at once whatever their sizes) |
| `--downloads` | 8 | the fetch stage's genomes and annotations downloaded at a time, by `curl` where it is on the PATH (one NCBI connection: ~30–40 MB/s by curl, ~10 by Python's urllib); their gzip checks run in one process of their own, a second core, by ISA-L's `igzip` where it is on the PATH (conda-forge's `isa-l`; ~3.5x Python's), else `pigz`, else Python |
| `--local_tmp DIR` | beside `OUT/genomes`; `TMPDIR` | a node's local disk (or a scratch file system without backups) for temporary files: the fetch stage's unfinished downloads (`.part` files) and their gzip checks, each file copied into `OUT/genomes` once whole; since 0.2.1 the process stage's too: each genome run's copy of its genome (read there by every step, not from `OUT` some six times), miniprot's pieces of a genome over 4 Gb, the rescue's ORFs, tiara's and vsearch's files, and `TMPDIR`'s (barrnap's, ITSx's); only the results are written to `WORK`; since 2026-11-10 `TMPDIR` for the package stage's tools as well. Room: up to ~3 times a genome's size per genome running (a few GB; ~10–20 GB for a mammal's, ~45 GB for a 30 Gb genome's pieces); `--scratch` before 2026-11-09, still taken |
| `--report_minutes` | 10 | minutes between progress reports of the genomes' downloads and the process stage's passes, after a report every minute of each one's first 10 (items left, time left; in the process stage the cores, memory and disk used) |
| `--tiara PATH`, `--no_tiara` | `tiara` on the PATH | tiara, which classifies the windows around the copies |
| `--sourmash PATH`, `--no_sketches` | `sourmash` on the PATH | sourmash, which sketches each genome before the prune (`metadata/sketches.zip`) |
| `--no_rescue` | the rescue pass on (with `esl-translate` on the PATH) | no rescue pass of the BUSCOs miniprot missed ([methods.md](methods.md)) |
| `--held_out_ranks` | `genus,family,order` | the OrthoDB benchmark's clade hold-out: a BUSCO dataset species' genome called again without its genus', family's and order's dataset species in turn, one more BUSCO search each ([benchmarks.md](benchmarks.md); `none`: its own species' hold-out alone; since 2026-11-09) |
| `--reversion` | off | stamp every done genome of an older genome record with the build's, its calls kept, so that it is not run again ([Resuming and rebuilding](#resuming-and-rebuilding); since 2026-11-09) |
| `--organelle_models` | `miniprot` | the assemblies' organelle protein-coding genes: `nhmmer`'s hits alone, or with `miniprot`'s exons by the organelle's code where nhmmer joined several hits (a gene with introns) or its single hit is not the coding sequence (under 90% of the query, or more internal stops than miniprot's model; since 2026-11-10; measured: [methods.md](methods.md); format 5) |
| `--intraspecific` | 3 | further assemblies of each species of `data/intraspecific_panel.tsv` run for the barcode gap (`intraspecific/`); 0: none |
| `--no_markers` | markers on | no Eukore marker loci (`data/markers.tsv`: RPB1, RPB2, EF-2, EF-1α, ...), the BUSCOs alone |
| `--no_tsa` | transcriptomes on | no INSDC transcriptomes (ENA's TSA sets) for the species without a genome |
| `--midori2` | none | a MIDORI2 RAW FASTA whose lineages `species.tsv` gives (its URLs change with each release: download it yourself); given again for its other genes' files. Its records are screened (since 2026-11-10; minimap2 and vsearch): one that is another organism's gene or another gene names no species (`reference_screen.tsv`; [methods](methods.md#the-reference-screen-since-2026-11-10)) |
| `--no_reference_screen` | screened | `--midori2`'s records not screened: every one names its species |
| `--bold_tsv` | none | a BOLD data package's TSV (`bin_uri`, `species`): its BINs in `crossrefs.tsv` (BOLD's bulk data needs a login) |
| `--refresh` | off | fetch again NCBI's assembly tables and quality reports and the annotation indexes: newer genomes and gene models taken, old genome files deleted once replaced ([Resuming and rebuilding](#resuming-and-rebuilding)) |
| `--veupathdb_key KEY` | the environment's `VEUPATHDB_API_KEY`, else none | a VEuPathDB user's API key (its profile's Web Service Access), for VEuPathDB's gene models; it is sent to VEuPathDB's sites alone, over https (since 0.2.2 a redirect elsewhere drops it), and written nowhere; the environment variable keeps it off the command line, which every user of the node can see |
| `--fcs_gx_db PREFIX`, `--run_gx` | none, `run_gx` | FCS-GX's database (e.g. `/dev/shm/gxdb/all`) and runner |
| `--ksgp_tax FILE` | KSGP 4.0, fetched | the KSGP release `species.tsv`'s KSGP lineages come from |
| `--miniprot`, `--hmmsearch` | on the PATH | the tools, where not on the PATH |
| `--no_benchmarks` | | write the package without its `benchmarks/` ([benchmarks.md](benchmarks.md)) |
| `--read_genomes N` | 3 | the read benchmark's test genomes per supergroup, kept through the prune (0: none) |
| `--read_compare NAME=FASTA,TABLE` | | another reference for the read benchmark and `overlap.tsv` (TABLE: record, tab, NCBI taxid), given again for each further one; the first is also `breadth.tsv`'s. E.g. EukDetect2's `eukdb.fasta` and `busco_taxid_genome_link.txt`; `tools/compare_references.py` makes MetaPhlAn 4's eukaryote markers and PR2's 18S into such pairs |
| `--minimap2` | on the PATH | the read benchmark's aligner |

### Resuming and rebuilding

Every stage is resumable: a rerun does only what is not done.
- A genome's results are in `WORK/done/<accession>.zip` (since 2026-11-10; `unzip -l` lists them, `unzip -p
  <file> done.json` prints one; a build before then, or `--done_folders`: the folder `WORK/done/<accession>/`).
  The file is written whole and renamed into place, so a killed job is redone whole, and a run again replaces
  the old results in one step. A WORK of folders needs nothing done to it: every stage reads both forms, and
  the process stage packs the folders (`--done_folders` keeps them).
- A genome run without a tool (tiara, barrnap) gets that step once the tool is there; one run without the rescue
  pass (`--no_rescue`, no `esl-translate`) is run again where it is made now (since 2026-11-08).
- A genome searched for fewer BUSCOs than the build calls (a run before format 4, or one with `--loci`) is run again.
- A fetch that stopped is taken up where it stopped. A source is downloaded again only where its version changed.
- The published rRNA reads' annotation is kept in `WORK/published`.
- A done genome is run again where better gene models than its run took are there (since 2026-11-07): the first of
  NCBI's, VEuPathDB's, Ensembl's and TOGA2's that its sources' indexes name for its assembly and that was not found
  to have none, where it comes before the one its run took (`done.json`'s `annotation` and `annotation_none`; a run
  before 2026-11-07 took NCBI's or none). The fetch stage fetches that annotation (and the genome again, where the
  prune deleted it); the process stage runs the genome again with it, its metadata, rRNA and held-out calls kept. This catches a genome run
  before an index named its assembly, before `--veupathdb_key` was given, or before its annotation was fetched.
- A species whose genome changed since the last plan (`WORK/genome_set.tsv`: TOGA2's assembly in its place; NCBI's
  newer assembly or version, after `--refresh`) gets the new one fetched and run; the old file (and its annotation)
  is deleted once the new one is done (`WORK/superseded.tsv`), its done folder kept.
- A done genome of an older genome record than the build's is run again (since 2026-11-09; the fetch stage fetches
  it again where the prune deleted it). The genome record (`genomes.GENOME_RECORD`, `done.json`'s `genome_record`)
  names what made a genome's calls: it is raised only where a change to Eukore's calls should not be mixed with older
  runs' in one package, not at every version (1: runs before 2026-11-09; 2: the coding sequences checked and recovered,
  a Duplicated locus' other copies kept for the choice by lineage: [methods.md](methods.md#checks-and-choices-format-6-since-2026-11-09)).
  `--reversion` stamps every done genome of an older record with the build's instead and keeps its calls, so that
  a build whose genomes were run by an earlier Eukore goes on without running them again: their copies keep the best
  HMM score's choice and their coding sequences unchecked (the package stage still checks each one's protein for
  stops and its lineage across species). Its `done.json` notes `reversioned` (the record it had, when, by which
  Eukore). The package stage says how many genomes are of an older record; `genomes.tsv`'s `genome_record` says each
  one's.

On a build begun before 2026-11-09 (every done genome of record 1), run the next stage with `--reversion` once to keep
its genomes' calls, or without it to run them all again with the checked coding sequences and the alternates:

```bash
python3 eukore_build.py --build_dir WORK --out OUT --stage process --reversion -t 120 -j 30 --memory_gb 200
```

`--refresh` (with the fetch stage, or `all`) fetches again what chooses the genomes and their gene models: NCBI's
assembly tables, its FCS-GX and ANI reports, and Ensembl's, VEuPathDB's and TOGA2's indexes (each old one kept where
its download fails). The next plan then takes NCBI's newer assemblies, and the rules above fetch and run what changed.
Without it, these sources stay as first fetched. For a newer taxonomy, delete `WORK/sources/taxdump.tar.gz` too.

### What a build leaves

| path | |
|---|---|
| `OUT/Eukore_<version>/` | the package ([FORMAT.md](FORMAT.md)), with its `benchmarks/` ([benchmarks.md](benchmarks.md)); its `manifest.json` is written last, so a folder without one is unfinished |
| `OUT/release/` | since 2026-11-09, the package's release, laid out as Eukore's website (`website/`) is, so that `rsync -a OUT/release/ <site>/` publishes it: `data/stats.json` (species by genome source and supergroup, the marker categories, the headline benchmarks), `data/taxonomy.json` (the species as a tree, each taxon's counts), `data/downloads.json`; `downloads/Eukore_<version>/`: the complete package (`Eukore_<version>.tar.gz`), the marker genes with introns (`_marker_genes.tar.gz`) and without (`_marker_cds.tar.gz`: coding sequences, proteins), the rDNA and organelle markers (`_rrna_organelle.tar.gz`), the genomes' sourmash sketches as one collection (`_sketches.sig.gz`), each metadata table alone (`_species.tsv.gz`, `_genomes.tsv.gz`, `_assembly.tsv.gz`, ...) and `MD5SUMS`. A bundle's files sit under `Eukore_<version>/` as in the package. A rerun keeps the downloads of the same package (its manifest's MD5). `--stage release` writes it again; `--no_release`, none |
| `OUT/genomes/` | the genomes kept after the prune (`<accession>.fna.gz`, as NCBI serves them; `manifest.tsv` lists kept and deleted). `genomes.tsv`'s `local` column names them |
| `WORK/sources/` (or `--sources`) | the sources, with their versions and checks (`source.json`, `labels.json`, ...) |
| `WORK/done/<accession>.zip` (`WORK/done/<accession>/` before 2026-11-10 and with `--done_folders`) | each genome's calls (`loci.tsv`, `copies.fna`, `features.tsv`), windows, rRNA (`rrna.tsv`), FCS-GX's ranges; for a species of BUSCO's dataset, its calls made without its own proteins (`held_out.tsv`, the OrthoDB benchmark's) |
| `WORK/organelle/`, `WORK/published/` | the organelle genes' searches; the published rRNA reads' annotation |

### Checking and reading a package

```bash
python3 -c "from eukore import package; print(package.check('OUT/Eukore_2026.11', deep=True) or 'whole')"
```

`check` compares every file with the manifest's size (and, with `deep`, its MD5), and says what is missing or
changed. The tables are plain TSV ([FORMAT.md](FORMAT.md)):

```python
from eukore import package, ncbi
pkg = "OUT/Eukore_2026.11"
species = package.read_table(pkg, "species.tsv")       # taxid, lineages in NCBI, PR2, EUKARYOME, KSGP, ...
copies = package.read_table(pkg, "busco/copies.tsv")   # each BUSCO copy's call, checks and confidence
loci = package.read_table(pkg, "busco/loci.tsv")       # each BUSCO's call rates, to choose loci by
seqs = package.read_copies(pkg, "10248812at2759")      # {accession: the copy, introns and all}
tax = ncbi.Ncbi(f"{pkg}/taxonomy")                     # the NCBI taxonomy the package was built on
```

A reader that is not Eukore's own (protal's) reads the tables by their header, so that new columns do not break it.

### Into protal

protal's database build converts a package into a GTDB release folder's eukaryote files:

```bash
python3 scripts/download_ref_data.py ... --eukdb OUT/Eukore_2026.11 [--eukdb_loci TABLE]
```

It takes protal's 30 BUSCOs (`--eukdb_loci` others), names the species by its KSGP and masks and fits the copies.
protal's `docs/databases.md` describes it.
