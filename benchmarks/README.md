# Eukore's benchmarks and figures

This folder holds what a paper about Eukore needs and the database build does not:
- the comparisons of a finished Eukore package with other references (EukDetect2, MetaPhlAn 4, PR2);
- the read benchmark;
- every figure.

The rest of the checkout makes the database and its website. The build keeps only its own checks of the package
(completeness, OrthoDB, resolution, rRNA, names, composition, run time; [docs/benchmarks.md](../docs/benchmarks.md)),
whose headline figures the website shows. Until 2026-11-12 the build also ran the comparisons, the read benchmark and
the figures (`--read_compare`, `--read_genomes`, `eukore/figures.py`).

Everything here runs on a finished package, in its own environment, with its large downloads in a folder you choose.

## Three steps

```bash
bash benchmarks/eukore-bench-prep.sh BENCH OUT/Eukore_2026.11 [--genomes OUT/genomes] [--taxdump SOURCES/taxdump.tar.gz]
bash benchmarks/eukore-benchmark.sh  BENCH OUT/Eukore_2026.11 RESULTS [-t 32] [--tmp /local/disk]
bash benchmarks/eukore-figures.sh    BENCH OUT/Eukore_2026.11 RESULTS [--format pdf]
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

[`env.yaml`](env.yaml) holds only what these steps run: Python, numpy, matplotlib, minimap2, curl and isa-l. The
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
| `read_background/` | The read benchmark's 7 bacteria, an archaeon and human chromosome 21 ([`read_background.tsv`](read_background.tsv)) | 40 MB | NCBI |
| `read_genomes/` | The read benchmark's test genomes (`tests.tsv`): up to `--read_genomes` (default 3) per supergroup, those with the most congeners in the package first. Linked from `--genomes` (the build's folder) where they are there, else downloaded | Depends on the genomes chosen (a vertebrate's is ~1 GB); nothing where linked | NCBI |
| `compare/` | The other references made comparable (`eukore_bench/references.py`): `<name>.fna.gz` and `<name>_taxids.tsv` (PR2, MetaPhlAn 4's eukaryote markers), and `references.tsv` listing every reference with its FASTA and record-to-taxid table. EukDetect2's pair is used as it is | Smaller than the downloads | |

`--only taxdump,pr2,metaphlan,eukdetect2,background,tests,convert` runs some of the steps.

## What is compared

`eukore-benchmark.sh` copies the package's own `benchmarks/` tables into `RESULTS`, then adds:
- **`breadth.tsv`**: the build's counts per supergroup (NCBI's species with an assembly, Eukore's), with the first
  other reference's columns. The order is EukDetect2, MetaPhlAn 4, PR2: the first that is there.
- **`overlap.tsv`**: the build's rows (PR2's names, EUKARYOME, KSGP, ROD, BUSCO's dataset) and one per other reference,
  each against every Eukore species (genome or rRNA), per supergroup.
- **`reads.tsv`, `detection.tsv`**: the read benchmark (`eukore_bench/reads.py`; below).
- **`summary.json`, `summary.md`**: the package's headline figures with these added.

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
| `eukore_bench/` | The Python of the three steps: `prep.py`, `compare.py`, `reads.py`, `references.py`, `figures.py`. It reads packages with the database's own code (`eukore/` of this checkout) |
| `read_background.tsv` | The read benchmark's background genomes |
| `paper/` | The paper's own scripts and figures (not in git) |
