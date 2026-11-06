# Benchmarks

Every Eukore build measures the package it writes: how complete its genomes are, whether its gene calls are the right
genes, how well its markers tell species apart, how whole its rRNA operons are, how its species are named in the 18S
references, and what the build cost. The tables go into the package (`Eukore_<version>/benchmarks/`,
[FORMAT.md](FORMAT.md)) and its headline figures into `manifest.json`, so a release carries its own report.
`docs/benchmarks/plot_benchmarks.py` draws them.

## Running them

During a build, nothing to do: the package stage measures the package once its tables are written, before its
manifest (`--no_benchmarks` skips them). A benchmark whose inputs are missing is left out, and
`benchmarks/summary.json`'s `left_out` says why; a failing benchmark never stops the package.

Any package can be measured again afterwards, and the figures drawn:

```bash
python3 -m eukore.benchmarks OUT/Eukore_2026.11 --work WORK -t 8      # into OUT/Eukore_2026.11/benchmarks/
python3 docs/benchmarks/plot_benchmarks.py OUT/Eukore_2026.11         # into docs/benchmarks/figures/
```

`--work` gives BUSCO's dataset (`WORK/sources/eukaryota_odb12`, for the OrthoDB benchmark), NCBI's whole taxonomy
(`WORK/sources/taxdump.tar.gz`, else the package's own `taxonomy/`) and the genomes' run times (`WORK/done/`); each
can also be given alone (`--odb`, `--taxdump`). `--out` writes the tables elsewhere (a measured-again package's
manifest does not list them). The plots need matplotlib (in `envs/eukore.yaml`); `--format svg` or `pdf` for print.

## What is measured

### Completeness (`completeness.tsv`)

BUSCO's own measure, on Eukore's calls: per species' genome, its BUSCOs of eukaryota_odb12 that are Single,
Duplicated, Fragmented or Missing, and `complete`, the Single and Duplicated over all BUSCOs called. It is read per
supergroup: BUSCO's dataset is mostly animals, fungi and plants, and a caller that reaches only those would show it
here. Each genome's copies' confidence tiers, assembly level and contigs are beside it, so that a poor genome is told
from a poor call.

![BUSCO completeness per supergroup](benchmarks/figures/completeness.png)

![Completeness against contiguity](benchmarks/figures/completeness_vs_assembly.png)

### The loci a user chooses among (`busco/loci.tsv`)

Not a benchmark table but the package's own: each BUSCO's call rate over all species and per supergroup, which a
user (protal's euk30) chooses its loci by ([loci.md](loci.md)).

![Each BUSCO's call rate per supergroup](benchmarks/figures/loci_call_rates.png)

### Gene models against OrthoDB (`orthodb.tsv`)

BUSCO's dataset holds, for each BUSCO, OrthoDB's protein in each of its ~450 species: that species' annotated gene.
Where a package species is one of them, its call is compared with it. The call's coding sequence (its exons, from
`busco/features.tsv`) is translated by the species' genetic code, and the two proteins' 5-mers are compared both
ways: `recall`, the share of OrthoDB's protein found in the call, and `precision`, the share of the call found in
OrthoDB's. The right gene model of the right gene gives both near 1:
- `concordant`: both 0.9 or more;
- `similar`: both 0.5 or more (another exon end, an isoform, another assembly of the species);
- `discordant`: another gene (a paralog, or an ohnolog of a duplicated genome) or a broken model;
- `missed`: OrthoDB annotates the gene, Eukore called nothing.

The package's calls are queried with the proteins of the dataset species nearest to the genome, for these species
their own (`guided`): they find OrthoDB's own protein again, which is circular. So the process stage calls a dataset
species' genome a second time with its own proteins (its species' and strains') left out of the queries, its nearest
other dataset species' and the ancestral variants kept (`WORK/done/<accession>/held_out.tsv`, since 2026-11-06; the
package keeps its first calls). The `held_out_` columns compare these calls with OrthoDB's gene model: the measure to
quote, of a genome whose nearest dataset relative is another species. It still measures the gene structure and the
copy chosen near the dataset species, not the reach into species far from every one (completeness per supergroup
shows that). OrthoDB's gene may also be annotated on another assembly of the species.

The held-out calls cost one more BUSCO search for each of the ~450 dataset species' genomes (their seconds in
`done.json`'s `held_out`, not in `runtime.tsv`). A genome run before 2026-11-06 is run again by `--stage process`
where its genome is still on disk (never fetched again for it); a package measured without them says so in
`summary.json`'s `left_out`.

![Calls against OrthoDB's own gene models](benchmarks/figures/orthodb.png)

### Resolution (`resolution.tsv`)

How well a marker tells species apart, the question a read profiler asks of it. For each copy, its nearest copy of
another species (vsearch over the copies' coding sequences, the search of the `own_species` check; one 18S per
species, its best record's of 1,000 bases or more): their identity and the deepest rank the two species share
(genus .. phylum, or none). Per marker (each BUSCO, all BUSCOs, the 18S) and shared rank, the table gives the
identity's quantiles and the copies under 97% and 99%. A copy under 99% to every other species is one that a read
assignment at 99% keeps to its species.

Read it with the package's sampling in mind: a species whose genus holds no other genome can only meet a farther
relative.

![Nearest other species by shared rank](benchmarks/figures/resolution.png)

### Confidence (`completeness.tsv`, `rrna.tsv`)

The copies' and records' tiers (`high`, `medium`, `low`), the sum of the checks that set them
([FORMAT.md](FORMAT.md)): another organism's sequence (`own_species`, FCS-GX, tiara), a duplicated or fragmented
call, a length in its locus' outer 5%, a contig end, an assembly gap; for rRNA, ambiguous bases, homopolymers, the part
screen, a partial 28S, a read labelled by its genus. Published reads are mostly `medium`: their own-species check has
little to compare them with in a trial, and a read labelled by its genus is never `high`.

![Confidence tiers](benchmarks/figures/confidence.png)

### rRNA operons (`rrna.tsv`)

Per source (the genomes, curated units, EUKARYOME, ROD) and kind (an operon, an 18S alone): the records holding each
region (18S, ITS1, 5.8S, ITS2, 28S, an inline 5S), a whole 28S, an 18S intron.

![rRNA records holding each region](benchmarks/figures/rrna.png)

### Names in the 18S references (`labels.tsv`)

Per supergroup and reference (PR2, EUKARYOME, KSGP), the species named by their name, an NCBI synonym, their genus
alone, or not at all: how far a user can carry a package species into those references' taxonomies.

![Species named in each reference](benchmarks/figures/labels.png)

### Run time (`runtime.tsv`, `stages.tsv`)

Each genome's run (its BUSCO calls, windows, tiara, rRNA) in seconds against its size, from `WORK/done/*/done.json`
(genomes run since 2026-11-06), and the package stage's steps.

![Run time](benchmarks/figures/runtime.png)

## The figures here

The figures above are those of a trial build, not of Eukore: the seven clades of the 2026-11-05 trial (*Saccharomyces*,
*Candida albicans*, *Blastocystis*, *Drosophila melanogaster*, *C. elegans*, *Arabidopsis*, *Trichuris trichiura*;
33 genomes, 1 Gb), all 129 BUSCOs, built on 2026-11-06 on 6 threads (26 minutes for the genomes, 6 for the package
stage; the published rRNA reads annotated in an earlier run). They show what the benchmarks look like; a release's
own are in its `benchmarks/` (`python3 docs/benchmarks/plot_benchmarks.py` redraws these from
`Eukore_2026.11-trial`).

What the trial shows:
- **Completeness** follows BUSCO's dataset: the plant 100%, the fungi and animals 81-85%, *Blastocystis* (SAR, 16
  genomes) 54%, with the most Fragmented BUSCOs (13%).
- **OrthoDB**: of the 321 BUSCOs OrthoDB annotates in four of the trial's species, Eukore called 96%, 84% of them
  OrthoDB's own gene model; in *S. mikatae* and *S. paradoxus* 93-94%. *Blastocystis hominis*' calls are more often
  `similar` (25 of 98) or `discordant` (16): where OrthoDB's gene model is of another assembly or isolate, such
  differences are expected, and these are the calls to look at first. BUSCO's dataset holds only 3 *S. cerevisiae*
  proteins; one (10249111at2759, 84 residues in OrthoDB, 386 in Eukore's call) shares nothing with the call.
- **OrthoDB, own proteins left out** (the held-out calls): the yeasts hardly change (*S. mikatae* 105 concordant of
  114, *S. paradoxus* 97 of 106: another *Saccharomyces* is in the dataset), but *Blastocystis hominis*, with no
  close relative in BUSCO's dataset, falls from 52 concordant and 5 missed of 98 to 21 and 32. Its own proteins made
  its calls look better than a genome far from the dataset gets: the held-out figure (80% concordant of the 88%
  called, over all four) is the one to quote, and the divergent lineages are where the calls need work.
- **Resolution**: 65% of the BUSCO copies are under 99% identical to every other species, against 37% of the species'
  18S: among congeners the BUSCOs tell species apart where the 18S cannot.

| benchmark | figure |
|---|---|
| package | Eukore_2026.11-trial (a trial: 7 clades), 33 species' genomes, 129 BUSCOs |
| median complete BUSCOs per genome | 72.9% (6.1% of genomes at 90% or more) |
| BUSCOs OrthoDB has a gene model of, called | 96.0% of 321 (4 species) |
| of those called, OrthoDB's gene model (5-mers ≥ 90% both ways) | 84.1% |
| the same, called without the species' own proteins | 87.5% called, 79.7% of them OrthoDB's gene model (4 species) |
| BUSCO copies under 99% identical to any other species | 65.4% of 3,387 |
| 18S copies under 99% identical to any other species | 36.6% of 41 |
| rRNA records; operons with a whole 28S | 2,719; 99.6% of 2,551 |
| species named in EUKARYOME by name or synonym | 30.2% |
| species named in KSGP by name or synonym | 34.9% |
| species named in PR2 by name or synonym | 30.2% |
| genome run (BUSCOs, rRNA, tiara), median | 2.3 min on 2 threads, 93 min per Gb of genome |

## Earlier: against EukDetect2 (2026-11-05)

Before the benchmarks, the first trial (the same seven clades, 30 BUSCOs) was compared with EukDetect2's database:
- 42 species against EukDetect2's 19, 1,150 loci against its 491. The genome stage took 8.6 minutes, the published
  rRNA 21, on 6 threads.
- Of the 19 species both hold, Eukore called 473 of EukDetect2's 491 loci, and 50 more.
- Where the genome is the same (the yeasts, *C. albicans*), its copies hold all of EukDetect2's copy's 21-mers, at
  the same length.
- The model species got 29-30 loci each.

## Not measured yet

- **Read assignment.** Reads simulated from held-out genomes, profiled against the package by a read profiler
  (protal): the end-to-end measure, which needs a profiler and so is not part of Eukore's build.
- **Divergent genomes.** OrthoDB's gene models exist only for the dataset's species; how far from them a call stays
  right is seen indirectly, in completeness per supergroup.
- **The organelle genes.** They have no checks across species yet, and so no benchmark.
