# Choosing loci

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

## The loci in prokaryotes

EukDetect2 dropped the BUSCOs it found in 971 gut bacteria, and the 30 loci are among those it kept. On 2026-11-05 that
check was redone on Eukore's own caller. The panel was 176 prokaryote genomes, one per GTDB genus of protal's local
inputs: 158 bacteria of 22 phyla and 18 archaea of 6. Each was searched with all ~23,000 dataset proteins of the 30
loci, which is stricter than a build (the 30 proteins of the nearest species). The results:

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

## Classic markers among the BUSCOs

`tools/classic_markers.py` searches the proteins of ten classic phylogenetic and barcode markers (yeast's and human's,
`data/classic_markers.faa`: RPB1, RPB2, EF-1α/TEF1, β-tubulin/TUB2, actin, HSP90, histone H3, GAPDH, EF-2,
calmodulin) with every BUSCO HMM and writes `data/classic_markers.tsv`. On eukaryota_odb12 (2026-11-07) none of them is
a BUSCO: no HMM scores any of their proteins near its cutoff (the best, HSP90 against 189968at2759, -1.9 bits). BUSCO
keeps single-copy genes, and these are gene families (actin, tubulin, EF-1α, HSP90, histones) or were not chosen
(RPB1, RPB2). Fungal multilocus phylogenies' RPB1, RPB2 and TEF1 would be new loci, not a name for existing ones.
