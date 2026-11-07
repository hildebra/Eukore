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
a BUSCO: no HMM scores any of their proteins near its cutoff (the best, HSP90 against 189968at2759, -1.9 bits). Nor
was any of them in eukaryota_odb10 (255 BUSCOs): none of its HMMs scores them at all.

Why (checked 2026-11-07): BUSCO takes an OrthoDB group of a lineage only where it is single-copy in 93% or more of the
lineage's species (90% in small lineages; Tegenfeldt et al. 2025, Nucleic Acids Res, doi:10.1093/nar/gkae987). The 129 BUSCOs are single-copy in 92 to
99% of the dataset's 456 species (median 95%, `info/ogs.id.info`). The classic markers' Eukaryota-level groups in
OrthoDB v12 (found by their yeast and human UniProt proteins; `phyletic_profile` of data.orthodb.org's API, of 5,952
species) are nearly universal but single-copy in few species:

| marker | OrthoDB group (Eukaryota) | species present | single-copy in | genes per species | likely why several copies |
|---|---|---|---|---|---|
| RPB1 | 4999407at2759 "DNA-directed RNA polymerase" | 99.6% | 6.3% | 2.2 | grouped with its paralogs, the largest subunits of RNA polymerases I and III (and plants' IV and V) |
| RPB2 | 4971508at2759 "DNA-directed RNA polymerase I" | 99.8% | 0.7% | 3.2 | grouped with polymerases I and III's second subunits (and plants' IV, V) |
| TEF1 (EF-1α) | 5014569at2759 | 95.8% | 2.5% | 3.3 | duplicated in many lineages (yeast's TEF1 and TEF2) |
| TUB2 (β-tubulin) | 4581638at2759 | 90.9% | 3.1% | 6.3 | a multigene family (plants, animals) |
| actin | 5012414at2759 | 97.4% | 5.6% | 4.5 | a multigene family, with actin-related proteins |
| HSP90 | 4969507at2759, 5015266at2759 | 46%, 44% | 28%, 1.6% | 1.5, 3.1 | cytosolic, ER (GRP94), mitochondrial (TRAP1) and plastid forms; split over groups, each in under half the species |
| histone H3 | 3474295at2759 | 4% | 0.5% | 3.4 | dozens of copies and variants; split into many lineage groups, none near-universal |
| GAPDH | 436619at2759 | 92.4% | 48.9% | 1.9 | cytosolic and plastid forms, duplicates; the nearest of all, still half of 93% |
| EF-2 | 4858627at2759 | 99.8% | 0.7% | 4.0 | duplicated in many lineages (yeast's EFT1 and EFT2) |
| calmodulin | 1901499at2759 | 59% | 6.3% | 4.5 | a family of EF-hand proteins |

So none fails for being rare: they fail BUSCO's single-copy rule, mostly because OrthoDB groups them with paralogs
older than eukaryotes' last common ancestor or because lineages duplicated them. Phylogenies that use them pick the
right paralog by hand (RPB1 against RPA190 and RPC160), which a single-copy screen cannot. Fungal multilocus
phylogenies' RPB1, RPB2 and TEF1 would be new loci of Eukore's own, with a paralog-aware caller, not a name for existing
ones. (`tools/classic_markers.py` makes the first part of this check; the OrthoDB figures were read from its API.)
