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
- Python 3.12 or later;
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
- matplotlib: the benchmarks' figures (`docs/benchmarks/plot_benchmarks.py`).

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
| memory | ~2 GB for a trial's package stage; per genome job, miniprot's on that genome (a few GB, more for the largest plant and animal genomes); cmalign up to 2 GB per thread for the 18S regions; the fetch stage ~2 GB (the taxonomy) plus, per download of a TOGA2 `.2bit` genome, its longest chromosome (~250 MB for a mammal's; until 2026-11-08 ~6 times that: 12 downloads at once killed a build on a 10 GB node) |
| network | the fetch stage alone |

## Usage

### A first, trial build

A trial builds the package of a few clades on one machine, all stages in one run:

```bash
python3 eukore_build.py --work work --out out --groups fungi --species 4930 -t 8
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
python3 eukore_build.py --work WORK --out OUT --stage fetch --groups fungi,protozoa

# compute nodes, no network: each genome's BUSCOs and rRNA, then the organelle genes
python3 eukore_build.py --work WORK --out OUT --stage process --groups fungi,protozoa -t 64 -j 16

# then the package: the checks across species, the published rRNA, the prune, the package
python3 eukore_build.py --work WORK --out OUT --stage package --groups fungi,protozoa -t 64 --keep_gb 1000 \
    [--fcs_gx_db /dev/shm/gxdb/all] [--tiara .../bin/tiara] [--ksgp_tax .../KSGP_v4.4_plus2.tax]
```

Give every stage the same `--work`, `--out`, `--groups` and `--species`: each stage plans the same genome set from
them. `--stage all` (the default) runs the three in turn on one machine with internet.

- **fetch** gets the sources into `WORK/sources`, each kept once there and whole:
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
  - since 0.2, Eukore's marker set (`WORK/sources/eukore_markers`, `eukore/markers.py`): KOfam's `ko_list` and its
    markers' and decoys' profiles (out of KOfam's 1.6 GB archive, fetched in 16 ranges at once and deleted),
    Swiss-Prot's eukaryote proteins (~90 MB, EBI's FTP) and NCBI's proteins named as each marker (`data/markers.tsv`'s
    Entrez query: RefSeq's of every eukaryote, GenBank's of the protists; 3 per organism, by E-utilities), scored with
    the profiles (hmmsearch, on `-t` threads) into the markers' query proteins, one per species. About an hour, most
    of it NCBI's; what it was made of is the package manifest's `inputs.markers`.

  It also fetches the transcriptomes of species without a genome and the intraspecific panel's further assemblies,
  as it fetches genomes.

  Then it fetches each species' genome into `OUT/genomes` (its RefSeq copy where it has one). After 1, 2, 3, 5 and
  8 minutes, then every 10 (`--report_minutes`), it says the genomes fetched, failed and left, the Gb left, the download rate and the time
  left, estimated by the genomes' sizes; the process stage reports its genomes alike. Beside each genome it fetches
  its gene models and NCBI's assembly report: NCBI's GFF (RefSeq's or a GenBank submitter's; not asked for where the
  assembly summary gives no protein-coding genes; a 404 or a GFF without coding sequences is none), else VEuPathDB's (with
  `--veupathdb_key`), else Ensembl's, else TOGA2's (their indexes fetched once into `SRC/annotation_sources`); the process stage takes a BUSCO copy's
  model from it where the HMM confirms it ([methods.md](methods.md)), compares miniprot's models with it
  (`annotation.tsv`) and deletes it.
- **process** calls each fetched genome's BUSCOs and extracts its rRNA into `WORK/done/<accession>/`, `-j` genomes at a
  time with `-t / -j` threads each. A genome not fetched yet is left for the next run. Then it takes the organelle
  genes (`WORK/organelle`), and runs the intraspecific panel's further assemblies (their BUSCOs, rRNA and organelle
  genes), deleting their genome files once done.
- **package**, in turn:
  - runs FCS-GX where asked;
  - checks the copies across species;
  - since 0.2: aligns the copies' proteins to BUSCO's HMMs (hmmalign; cached in `WORK/align`), places their introns,
    checks each genome's genetic code;
  - annotates the published rRNA reads (barrnap and ITSx over `-t` threads, 250 reads a batch);
  - runs the rRNA checks and the part screen;
  - since 0.2: finds the 18S variable regions (cmalign; cached in `WORK/cmalign`) and runs the in-silico PCR;
  - takes the organelle genes of any genome not searched yet, checks them across species (since 0.2);
  - since 0.2: measures the intraspecific panel's barcode gap;
  - names each species in PR2, EUKARYOME, KSGP and (since 0.2) SILVA, UNITE and MIDORI2, and finds its ids in Open
    Tree's taxonomy (`crossrefs.tsv`);
  - describes each done genome still on disk that lacks its metadata (`eukore/metadata.py`: a genome done before
    2026-11-07), as the process stage does after its runs;
  - prunes `OUT/genomes` to `--keep_gb` (the model species first, then the smallest);
  - writes `OUT/Eukore_<version>/`, with its benchmarks ([benchmarks.md](benchmarks.md)).

### Options

| option | default | |
|---|---|---|
| `--work DIR` | required | the build's work folder: sources, each genome's results, caches |
| `--out DIR` | required | where `Eukore_<version>/` and `genomes/` go |
| `--stage` | `all` | `fetch`, `process`, `package` or `all` |
| `--groups` | all six | NCBI's assembly groups: `fungi`, `protozoa`, `plant`, `invertebrate`, `vertebrate_mammalian`, `vertebrate_other` (comma-separated) |
| `--species` | | a trial: only the species inside these NCBI taxids' clades (comma-separated, or a file of one per line) |
| `--loci FILE` | every BUSCO | call only these BUSCOs (one per line, or a loci table's first column): a quicker trial |
| `--version` | `YYYY.MM` of today | the package's version, in its folder's name (Eukore's own version, 0.2, is in the manifest's `eukore`) |
| `-V` | | print Eukore's own version and exit |
| `--keep_gb` | 1000 | the genomes kept in `OUT/genomes`, GB: cut back to it before the package is written (its `genomes.tsv` says which are kept) |
| `-t`, `--threads` | 8 | threads in all |
| `-j`, `--jobs` | threads / 4 | genomes at a time |
| `--downloads` | 8 | the fetch stage's genomes and annotations downloaded at a time, by `curl` where it is on the PATH (one NCBI connection: ~30–40 MB/s by curl, ~10 by Python's urllib); their gzip checks run in one process of their own, a second core, by ISA-L's `igzip` where it is on the PATH (conda-forge's `isa-l`; ~3.5x Python's), else `pigz`, else Python |
| `--scratch DIR` | beside `OUT/genomes` | the fetch stage's unfinished downloads (`.part` files) and their gzip checks, e.g. a node's local disk or a scratch file system without backups; each file is copied into `OUT/genomes` once whole |
| `--report_minutes` | 10 | minutes between progress reports of the genomes' downloads and runs (genomes left, time left) |
| `--tiara PATH`, `--no_tiara` | `tiara` on the PATH | tiara, which classifies the windows around the copies |
| `--sourmash PATH`, `--no_sketches` | `sourmash` on the PATH | sourmash, which sketches each genome before the prune (`metadata/sketches.zip`) |
| `--no_rescue` | the rescue pass on (with `esl-translate` on the PATH) | no rescue pass of the BUSCOs miniprot missed ([methods.md](methods.md)) |
| `--organelle_models` | `miniprot` | the assemblies' organelle protein-coding genes: `nhmmer`'s hits alone, or with `miniprot`'s exons by the organelle's code where nhmmer joined several hits (a gene with introns; measured: [methods.md](methods.md); format 5) |
| `--intraspecific` | 3 | further assemblies of each species of `data/intraspecific_panel.tsv` run for the barcode gap (`intraspecific/`); 0: none |
| `--no_markers` | markers on | no Eukore marker loci (`data/markers.tsv`: RPB1, RPB2, EF-2, EF-1α, ...), the BUSCOs alone |
| `--no_tsa` | transcriptomes on | no INSDC transcriptomes (ENA's TSA sets) for the species without a genome |
| `--midori2` | none | a MIDORI2 RAW FASTA whose lineages `species.tsv` gives (its URLs change with each release: download it yourself) |
| `--bold_tsv` | none | a BOLD data package's TSV (`bin_uri`, `species`): its BINs in `crossrefs.tsv` (BOLD's bulk data needs a login) |
| `--refresh` | off | fetch again NCBI's assembly tables and quality reports and the annotation indexes: newer genomes and gene models taken, old genome files deleted once replaced ([Resuming and rebuilding](#resuming-and-rebuilding)) |
| `--veupathdb_key KEY` | the environment's `VEUPATHDB_API_KEY`, else none | a VEuPathDB user's API key (its profile's Web Service Access), for VEuPathDB's gene models; it is sent to VEuPathDB alone and written nowhere |
| `--fcs_gx_db PREFIX`, `--run_gx` | none, `run_gx` | FCS-GX's database (e.g. `/dev/shm/gxdb/all`) and runner |
| `--ksgp_tax FILE` | KSGP 4.0, fetched | the KSGP release `species.tsv`'s KSGP lineages come from |
| `--miniprot`, `--hmmsearch` | on the PATH | the tools, where not on the PATH |
| `--no_benchmarks` | | write the package without its `benchmarks/` ([benchmarks.md](benchmarks.md)) |
| `--read_genomes N` | 3 | the read benchmark's test genomes per supergroup, kept through the prune (0: none) |
| `--read_compare NAME=FASTA,TABLE` | | another reference for the read benchmark and `overlap.tsv` (TABLE: record, tab, NCBI taxid), given again for each further one; the first is also `breadth.tsv`'s. E.g. EukDetect2's `eukdb.fasta` and `busco_taxid_genome_link.txt`; `tools/compare_references.py` makes MetaPhlAn 4's eukaryote markers and PR2's 18S into such pairs |
| `--minimap2` | on the PATH | the read benchmark's aligner |

### Resuming and rebuilding

Every stage is resumable: a rerun does only what is not done.
- A genome's results are in `WORK/done/<accession>/`. `done.json` is written last, so a killed job is redone whole.
- A genome searched for fewer BUSCOs than the build calls (a run before format 4, or one with `--loci`) is run again.
- A fetch that stopped is taken up where it stopped. A source is downloaded again only where its version changed.
- The published rRNA reads' annotation is kept in `WORK/published`.
- A done genome is run again where better gene models than its run took are there (since 2026-11-07): the first of
  NCBI's, VEuPathDB's, Ensembl's and TOGA2's that its sources' indexes name for its assembly and that was not found
  to have none, where it comes before the one its run took (`done.json`'s `annotation` and `annotation_none`; a run
  before 2026-11-07 took NCBI's or none). The fetch stage fetches that annotation (and the genome again, where the
  prune deleted it); the process stage runs the genome again with it, its metadata kept. This catches a genome run
  before an index named its assembly, before `--veupathdb_key` was given, or before its annotation was fetched.
- A species whose genome changed since the last plan (`WORK/genome_set.tsv`: TOGA2's assembly in its place; NCBI's
  newer assembly or version, after `--refresh`) gets the new one fetched and run; the old file (and its annotation)
  is deleted once the new one is done (`WORK/superseded.tsv`), its done folder kept.

`--refresh` (with the fetch stage, or `all`) fetches again what chooses the genomes and their gene models: NCBI's
assembly tables, its FCS-GX and ANI reports, and Ensembl's, VEuPathDB's and TOGA2's indexes (each old one kept where
its download fails). The next plan then takes NCBI's newer assemblies, and the rules above fetch and run what changed.
Without it, these sources stay as first fetched. For a newer taxonomy, delete `WORK/sources/taxdump.tar.gz` too.

### What a build leaves

| path | |
|---|---|
| `OUT/Eukore_<version>/` | the package ([FORMAT.md](FORMAT.md)), with its `benchmarks/` ([benchmarks.md](benchmarks.md)); its `manifest.json` is written last, so a folder without one is unfinished |
| `OUT/genomes/` | the genomes kept after the prune (`<accession>.fna.gz`, as NCBI serves them; `manifest.tsv` lists kept and deleted). `genomes.tsv`'s `local` column names them |
| `WORK/sources/` | the sources, with their versions and checks (`source.json`, `labels.json`, ...) |
| `WORK/done/<accession>/` | each genome's calls (`loci.tsv`, `copies.fna`, `features.tsv`), windows, rRNA (`rrna.tsv`), FCS-GX's ranges; for a species of BUSCO's dataset, its calls made without its own proteins (`held_out.tsv`, the OrthoDB benchmark's) |
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
