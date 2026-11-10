# rDNA copies by read depth (2026-11-06, a prototype)

An assembly collapses an rDNA array into a few units, so the units found in it are a lower bound. `tools/rdna_depth.py`
measures the copies from the genome's own sequencing reads instead (`eukore/depth.py`; columns in
[FORMAT.md](FORMAT.md)'s `depth.tsv`). It runs on a built package, apart from the build, in three resumable steps;
only the first needs the network:

```
python3 tools/rdna_depth.py stage   --package OUT/Eukore_<v> --work DWORK --species 4932,6239 -t 4   # network
python3 tools/rdna_depth.py measure --package OUT/Eukore_<v> --work DWORK -t 16                     # no network
python3 tools/rdna_depth.py table   --package OUT/Eukore_<v> --work DWORK --into_package
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
   divided by the median of the single-copy BUSCOs' exons' mean depths (since 2026-11-10; until then the whole
   copies', whose introns' repeats draw reads from all their genomic copies). Both alleles map, so the number is per haploid
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

## Not done yet
It has run on 10 genomes, not the release:
- it is not yet a stage of `eukore_build.py`;
- PacBio and Nanopore runs are untested;
- runs that lie on a specimen's other samples (Darwin Tree of Life) are not found;
- a GC correction (rDNA against the BUSCOs' GC, in `depth/<accession>.targets.tsv`) is not made;
- organelle copies are unvalidated.
