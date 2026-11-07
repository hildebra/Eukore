<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/logo/eukore_logo_dark.png">
    <img src="docs/logo/eukore_logo.png" alt="Eukore" width="360">
  </picture>
</p>

<p align="center">
  <b>A checked, complete reference of eukaryote marker genes, built from every eukaryote genome NCBI holds.</b><br>
  BUSCO genes with all their introns · whole rDNA operons · organelle rRNA, COX1, CYTB, rbcL, matK and tufA · assembly metadata
</p>

<p align="center">
  <a href="#what-is-in-eukore">Contents</a> ·
  <a href="#get-eukore">Get it</a> ·
  <a href="#use-a-package">Use it</a> ·
  <a href="#build-eukore">Build it</a> ·
  <a href="#benchmarks">Benchmarks</a> ·
  <a href="#documentation">Docs</a>
</p>

---

Eukore is two things:
- **a database**: for each eukaryote species with a genome at NCBI, its genes of BUSCO's 129 near-universal
  eukaryote orthologs, its rDNA operons and its organelle rRNA, each sequence complete, labelled by NCBI's taxonomy,
  checked across species and given a confidence tier;
- **the software that builds it** (this repository): from NCBI's assemblies to a versioned, self-describing package,
  resumable on a cluster, with benchmarks measured on every build.

It is made for tools that detect and profile eukaryotes in metagenomes (protal reads it as its eukaryote reference),
and for anyone needing clean eukaryote marker genes: phylogenetics, primer design, contamination screens. Eukore's
version is 0.12; it was named EUKDB until 2026-11-06.

## What is in Eukore

| | per species | from | complete as |
|---|---|---|---|
| **BUSCO genes** | each of the 129 BUSCOs of eukaryota_odb12 called in its genome, with its exons and introns | NCBI's best assembly of the species, called by Eukore (miniprot, then BUSCO's HMMs) | the gene from first to last codon, introns and all, uncut and unmasked |
| **rDNA operons** | its 18S spans and whole units (18S, ITS1, 5.8S, ITS2, 28S, an inline 5S) | its assembly (barrnap, ITSx); EUKARYOME's and ROD's long reads | from the 18S's start to the 28S's end |
| **organelle genes** | mitochondrial and plastid SSU, LSU, 5S rRNA; mitochondrial COX1 and CYTB; plastid rbcL, matK and tufA | RefSeq's organelle genomes, else its assembly's organelle molecules (kept whole too) | the gene with its introns |
| **assembly metadata** | its genome's size, contiguity, gaps, GC and its spread, soft-masking, CpG, duplicates, telomeres, organelle molecules, tetranucleotides, a sourmash sketch; NCBI's own summary | measured while the genome is on disk, before the build deletes it | in the package's `metadata/` |
| **taxonomy** | NCBI's lineage and supergroup; its lineage in PR2, EUKARYOME and KSGP | NCBI's taxdump, cut to the package | tied to no GTDB release |
| **checks** | per copy and record: own-species agreement, 99% collapses, length outliers, contig ends, assembly gaps, tiara's and FCS-GX's contamination screens, an rRNA chimera screen; a `high`/`medium`/`low` tier with its reasons | Eukore's checks across species | flagged, never silently dropped |
| **benchmarks** | completeness, agreement with OrthoDB's gene models, marker resolution, rRNA completeness, naming, run time | measured by each build | in the package's `benchmarks/` |

All eukaryote groups of NCBI are covered: fungi, protists, plants, invertebrates and vertebrates (on NCBI's tables
of 2026-11, about 29,000 species). The genome of a species is NCBI's best by Eukore's ranking, never one NCBI or its
FCS-GX and ANI screens find contaminated or misidentified where a better one exists. 75 species with no assembly
(protists from EukProt) come with their records.

A package is plain TSV and gzipped FASTA with a `manifest.json` of checksums; its format is
[docs/FORMAT.md](docs/FORMAT.md) (format 4).

```
Eukore_2026.11/
├── manifest.json            versions, inputs, counts, every file's size and MD5 (written last)
├── species.tsv              one row per species: NCBI, PR2, EUKARYOME and KSGP lineages, what it holds
├── genomes.tsv              the species' genomes: accession, quality flags, rRNA copies in the assembly
├── taxonomy/                NCBI's taxonomy, cut to the package (nodes, names, merged)
├── busco/                   loci.tsv (call rates per BUSCO), copies.tsv (calls, checks, confidence),
│                            features.tsv (exons, introns, masked runs), copies/<locus>.fna.gz
├── rrna/                    records.tsv, records.fna.gz
├── organelle/               genes.tsv, genes.fna.gz, molecules.fna.gz (the assemblies' organelle molecules)
├── metadata/                assembly.tsv.gz, tetranucleotides.tsv.gz, sketches.zip (sourmash), ncbi.tsv.gz
└── benchmarks/              the build's measures of the package
```

## Get Eukore

**Prebuilt releases**: the first full release (built on NCBI's genomes of 2026-11) is being built; trial packages
of a few clades exist so far. Until a release is published, build it yourself (below): a trial of seven clades took
about an hour on 6 cores, the full reference needs a cluster.

## Use a package

The tables are plain TSV, read by any tool. With Eukore's own reader:

```python
from eukore import package, ncbi
pkg = "Eukore_2026.11"
print(package.check(pkg, deep=True) or "whole")           # every file against the manifest's size and MD5

species = package.read_table(pkg, "species.tsv")           # taxid, lineages in NCBI, PR2, EUKARYOME, KSGP, ...
loci = package.read_table(pkg, "busco/loci.tsv")           # each BUSCO's call rates, to choose loci by
copies = package.read_table(pkg, "busco/copies.tsv")       # each copy's call, checks and confidence
seqs = package.read_copies(pkg, "10248812at2759")          # {accession: the copy, introns and all}
tax = ncbi.Ncbi(f"{pkg}/taxonomy")                         # the NCBI taxonomy the package was built on
```

A typical use keeps the `high` and `medium` copies of the BUSCOs it chooses (by `busco/loci.tsv`'s call rates in the
supergroups it cares about; [docs/loci.md](docs/loci.md)) and masks, cuts or renames them to its own needs. protal's
database build does this with its own reader (`download_ref_data.py --eukdb Eukore_2026.11`): the 30 BUSCOs of its
euk30, KSGP names, copies fitted into 16 kb. A reader other than Eukore's should read the tables by their header, so
that new columns do not break it.

## Build Eukore

Eukore runs on Linux (or WSL) with a conda-style package manager; its tools are bioconda's.

```bash
git clone https://github.com/hildebra/Eukore.git && cd Eukore
micromamba create -f envs/eukore.yaml && micromamba activate eukore
python3 -m unittest discover -s tests -t tests            # ~80 tests with stand-in tools, seconds
```

A **trial** builds a few clades on one machine (here *Saccharomyces*, taxid 4930); its manifest marks it partial:

```bash
python3 eukore_build.py --work work --out out --groups fungi --species 4930 -t 8
```

The **full reference** is built in three stages, the first where the network is, the rest on compute nodes:

```bash
python3 eukore_build.py --work WORK --out OUT --stage fetch                  # sources and genomes (~6.5 TB, pruned as it goes)
python3 eukore_build.py --work WORK --out OUT --stage process -t 64 -j 16    # each genome's genes and rRNA
python3 eukore_build.py --work WORK --out OUT --stage package -t 64 \
    [--fcs_gx_db /dev/shm/gxdb/all] [--tiara .../bin/tiara]                  # checks, rRNA, labels, package, benchmarks
```

Every stage is resumable: a rerun does only what is not done. Optional tools: tiara (the windows around each copy
classified) and NCBI's FCS-GX with its ~470 GB database (contamination screen). Installation, options, resources and
outputs are in [docs/building.md](docs/building.md).

### How it works

1. **Genome set.** Per NCBI species, its best assembly: NCBI's contamination, misassembly and ANI flags excluded,
   FCS-GX's contaminated assemblies demoted, then RefSeq reference, representative, assembly level, contiguity.
2. **BUSCO calls.** miniprot aligns the 30 proteins of BUSCO's dataset nearest to the genome (and its ancestral
   variants) with the genome's genetic code; hmmsearch confirms each region against all BUSCO HMMs at BUSCO's
   cutoffs; the best copy is kept whole, its exons and introns annotated.
3. **rRNA.** barrnap finds each genome's rRNA genes and ITSx the ITS borders (with the lineage's profiles); whole
   units are kept uncut. EUKARYOME's and ROD's long reads are labelled by NCBI taxid and annotated alike.
4. **Organelle genes** (rRNA, COX1, CYTB, rbcL, matK, tufA) from RefSeq's organelle genomes, else from the assembly's organelle
   molecules (nhmmer, barrnap).
5. **Assembly metadata.** One pass over each genome while it is on disk (`eukore/metadata.py`), and a sourmash
   sketch: what a user needs of the genome once the build has deleted it.
6. **Checks across species.** Each copy against the same gene of every other species (own-species agreement,
   identical copies, length outliers), an rRNA part screen for chimeras, tiara and FCS-GX on each copy's 50 kb
   window; each copy and record gets a confidence tier with its reasons.
7. **Package and benchmarks**, the manifest written last.

Details, and what was checked on real genomes, are in [docs/methods.md](docs/methods.md).

## Benchmarks

Each build measures its package ([docs/benchmarks.md](docs/benchmarks.md)): BUSCO completeness per genome and
supergroup; its calls against OrthoDB's own gene models of the same species; how far each marker tells species apart
(the identity of each copy's nearest other species, by the rank they share); rRNA operon completeness; naming in
PR2, EUKARYOME and KSGP; run time. `docs/benchmarks/plot_benchmarks.py` draws them for any package.

On a trial package (seven clades, 2026-11-06, with the rescue pass; not a release):

| benchmark | figure |
|---|---|
| package | Eukore_2026.11-rescue (a trial: 7 clades), 33 species' genomes, 129 BUSCOs |
| median complete BUSCOs per genome | 79.8% (6.1% of genomes at 90% or more) |
| BUSCOs OrthoDB has a gene model of, called | 96.6% of 321 (4 species) |
| of those called, OrthoDB's gene model (5-mers ≥ 90% both ways) | 83.5% |
| the same, called without the species' own proteins | 93.8% called, 74.8% of them OrthoDB's gene model (4 species) |
| of those also called with its own proteins, on the same gene | 100.0% of 301, 91.4% of them sharing 90% of the span; 20 called by the rescue pass |
| BUSCO copies under 99% identical to any other species | 65.6% of 3,676 |
| 18S copies under 99% identical to any other species | 36.6% of 41 |
| rRNA records; operons with a whole 28S | 2,719; 99.6% of 2,551 |
| species named in EUKARYOME by name or synonym | 30.2% |
| species named in KSGP by name or synonym | 34.9% |
| species named in PR2 by name or synonym | 30.2% |
| genome run (BUSCOs, rRNA, tiara; before the rescue pass), median | 2.3 min on 2 threads, 93 min per Gb of genome |

<p align="center">
  <img src="docs/benchmarks/figures/completeness.png" alt="BUSCO completeness per supergroup" width="860"><br>
  <img src="docs/benchmarks/figures/resolution.png" alt="Identity of each copy's nearest other species, by shared rank" width="860">
</p>

Among congeners the BUSCO genes tell species apart where the 18S cannot (65% of copies under 99% to every other species,
against 37% of 18S). All figures and what they show: [docs/benchmarks.md](docs/benchmarks.md).

## Documentation

| | |
|---|---|
| [docs/building.md](docs/building.md) | install, a trial, a full build on a cluster, options, resources, resuming, outputs |
| [docs/methods.md](docs/methods.md) | how the genomes are chosen and the genes, operons and organelle genes made and checked |
| [docs/FORMAT.md](docs/FORMAT.md) | the package format, file by file and column by column |
| [docs/benchmarks.md](docs/benchmarks.md) | what each benchmark measures, how to run and draw them, the figures |
| [docs/loci.md](docs/loci.md) | choosing BUSCOs for a use; the BUSCOs with prokaryote homologues |
| [docs/read_depth.md](docs/read_depth.md) | rDNA and organelle copies by read depth (`tools/rdna_depth.py`, a prototype) |

## Repository layout

| path | |
|---|---|
| `eukore_build.py` | the build (stages fetch, process, package) |
| `eukore/` | the package: `genomes.py` (genome set, one genome's run), `odb.py`, `loci.py` (BUSCO calls), `rrna.py`, `screen.py` (rRNA and its chimera screen), `organelle.py`, `qc.py` (checks across species), `quality.py` (NCBI's FCS-GX and ANI reports), `labels.py` (PR2, EUKARYOME, KSGP), `package.py` (write, read, check), `benchmarks.py` |
| `data/` | the model species, curated rDNA units, the records of species without an assembly |
| `envs/` | the micromamba environments (`eukore.yaml`; `tiara.yaml`, tiara needing an older Python) |
| `tools/` | choosing loci: `rank_loci.py`, `eukaryote_panel.py`, `prokaryote_panel.py`, `prokaryote_separation.py` |
| `tests/` | unit and end-to-end tests with stand-in tools |
| `docs/` | documentation, the format, the benchmarks' plot script and figures, the logo |
| `LICENSE` | GPL-2.0 |

## Status and limitations

- Version 0.12. The full reference has not been built on NCBI's release yet; format 4 has been built on trial clades.
- FCS-GX has been run with stand-ins only (its database is ~470 GB).
- The organelle genes have been built in tests only; they have no checks across species and no confidence yet.
- `genomes.tsv`'s rDNA copy numbers are those of the assembly, a lower bound (arrays collapse);
  `tools/rdna_depth.py` measures a genome's own from its reads, a prototype run on 10 genomes
  ([docs/read_depth.md](docs/read_depth.md)). The intergenic spacer, 5S arrays away from the unit and the 28S's
  introns are not taken.
- Four BUSCOs have prokaryote homologues (archaeal ribosome biogenesis, cytochrome c oxidase assembly; see
  [docs/loci.md](docs/loci.md)); whether to replace them in protal's euk30 is open.

## Sources and credits

Eukore is built from public data, each under its own terms; check them before redistributing a package:
NCBI's GenBank and RefSeq assemblies, taxonomy, FCS-GX and ANI reports and RefSeq organelle genomes; BUSCO's
eukaryota_odb12 (OrthoDB 12.1); EUKARYOME 2.1; the rRNA Operon Database (ROD) 1.2; PR2 5.1.1; KSGP. The records of
species without an assembly come from EukDetect2 (CC BY 4.0; Shih, Zhao, Pollard & Lind 2026). It runs miniprot,
HMMER, barrnap, ITSx, vsearch, and optionally tiara and NCBI FCS-GX: please cite them with Eukore.

Eukore started inside protal's repository (`EUKDBbuild/`, 2026-11-05) and keeps that history. Questions and bug
reports: [GitHub issues](https://github.com/hildebra/Eukore/issues).

## License

Eukore's code is licensed under the GNU General Public License, version 2.0 ([LICENSE](LICENSE)). The data a package
holds keep their sources' terms (above).
