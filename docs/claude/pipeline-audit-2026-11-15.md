# Eukore pipeline audit, 2026-11-15

Audit of Eukore 0.2.6 (package format 7) at commit `e1ee769` on `main`, from the genome set to the package and its
benchmarks, and an answer to the question: is this now the best available representation of eukaryotic marker genes
in one unified, validated framework?

## Verdict

**In scope and integration, yes. As a validated resource, not yet.**

- **Nothing else public does what Eukore sets out to do.** A survey of the other references (section 2) found none
  that links genome-derived nuclear protein markers, rDNA operons and organelle barcode genes for the same species,
  with a check result on every record, across all eukaryote groups. EukDetect2, ROD and the PR2 interface each cover
  one part.
- **The claim cannot be made today, for three reasons.**
  1. No release is published. The one full build is Eukore 0.2.4, format 6. The current code (0.2.6, format 7) has
     never been built in full, and most of what this audit found is in code written after that build.
  2. Validation is uneven. It is strong for the BUSCO calls of animals, fungi and plants. It is thin for protists,
     and there is no accuracy measure at all for organelle genes, rRNA records or format 7's imported records
     (section 3).
  3. The audit found 31 defects that would put wrong, missing or mislabelled records into the next build, and 19
     more that bear on a published claim or on robustness (section 4). I re-checked about half of them myself;
     the rest stand as the stage audits reported them.
- **The software itself is in good order.** All 486 tests pass, coordinates and genetic codes are handled
  consistently, the package is written atomically with a checksummed manifest, and hold-out benchmarks are truly
  held out (section 5).

What the claim needs is listed in section 7. Which findings have since been fixed, on branch `audit-fixes`, is in
section 9.

## 1. What was audited, and how

| | |
|---|---|
| Code | `eukore_build.py`, every module of `eukore/`, `benchmarks/eukore_bench/`, `tools/eukore.sh`, at `e1ee769` |
| Tests | `python3 -m unittest discover -s tests -t tests` in WSL: 486 pass on Python 3.12 (5 skipped, tools absent) and on the build environment's Python 3.14 with real miniprot, HMMER, barrnap and vsearch (2 skipped: no minimap2, no cmalign there) |
| Evidence | the first full build's evaluation bundle, `local/Eukore_2026.11_evaluation.tar.gz` (Eukore 0.2.4, format 6, 33,664 genomes): manifest, `species.tsv`, `genomes.tsv`, `benchmarks/`, the final job's log |
| Other references | a web survey on 2026-11-15 of rRNA, organelle, protein-marker and ortholog resources |
| Earlier reports | `local/audit-2026-11-10.md`, `local/fix-plan-2026-11-10.md`, `local/evaluation-2026-11-14.md`: what they record as fixed or as known and open is not repeated here |

Method: six read-only audits ran in parallel, one per stage (orchestration; genome set; gene calls; rRNA; organelle
genes and taxonomy; checks, package and benchmarks). Each read its modules in full and ran small probes. I then
re-checked the main findings against the code and against the first build's tables. Each finding below says how it
stands:

- **confirmed**: I re-checked it myself, by reading the code, by a probe, or by recomputing from the first build's
  tables;
- **reported**: a stage audit verified it by trace or probe, and I did not repeat that;
- **suspected**: the mechanism is in the code, its effect on a real build is not measured.

Limits. Nothing was built: no finding in code newer than the first full build has been seen in a package. The HPC
work folder was not available. Figures about other databases were read by a web agent and are not re-checked against
the papers. While the audit ran, another session changed tracked files in this checkout (`README.md`, `website/`,
figures, `eukore/release.py`); those working-tree changes are not audited, and this report changes nothing but
itself.

## 2. Is it the best available?

### What exists

Figures for the other resources are from the survey of 2026-11-15 (their own sites and papers), not re-checked.

| resource | nuclear protein markers | rDNA | organelle genes | checks on its records | links the three per species |
|---|---|---|---|---|---|
| **Eukore** (first full build, format 6) | 129 BUSCOs and 9 classic markers as gene, CDS and protein; 33,664 species (5,028 of them from transcriptomes) | 4.66 M records; 38,190 species; 33,930 with a whole operon | 13,868 species (format 7 adds MIDORI2: about 214,000 species with COX1, measured on the inputs) | four tiers with reasons on every record | yes: 11,296 species hold all three, 24,367 hold BUSCO genes and rRNA |
| EukDetect2 (preprint 2026-06) | BUSCO genes of eukaryota_odb12; 6,594 species, fungi and protists only | none | none | BUSCOs with bacterial hits dropped, 99% collapse | no |
| ROD 1.2 | none | 69,480 operon variants from 11,935 assemblies (NCBI of mid-2023) | none | complete ends only | no |
| PR2 5.1.2 | none | about 240,000 18S, expert-curated | plastid 16S | manual curation | no |
| SILVA 144 | none | 266,197 eukaryote SSU (NR99); LSU still 138.2 | none | alignment quality; its new curation not yet applied to eukaryotes | no |
| EUKARYOME | none | 18S, ITS, 28S and long operon reads | none | curated | no |
| UNITE 10.0 | none | ITS, species hypotheses | none | reference sequences per hypothesis | no |
| MIDORI2 | none | none | mitochondrial genes of GenBank, as named there | none described | no |
| BOLD v5 | none | none | COI of about 1.7 M species | BIN clustering | no |
| RefSeq organelle | none | none | whole organelle genomes, about 18,000 mitochondrial taxa | NCBI curation | no |
| OrthoDB 12 / BUSCO | marker definitions; no per-genome calls released | none | none | n/a | no |
| NCBI Datasets, BlobToolKit | BUSCO scores per assembly; no sequences | none | none | n/a | no |
| MetaPhlAn 4 | about 489 eukaryote species | none | none | n/a | no |

The survey searched specifically for a competitor: a unified protein, rDNA and organelle reference, or a released
collection of BUSCO gene sequences across NCBI's eukaryote genomes. It found neither.

### Where Eukore leads

- **Breadth of genome-derived protein markers.** Every eukaryote group, animals and plants included: 33,664 species
  against EukDetect2's 6,594.
- **Gene structure.** Whole genes with introns, CDS and protein, the intron table and the genetic-code check. No
  other marker reference ships these.
- **Linking.** One species row ties its nuclear genes, operons, organelle genes, primer results and its names in
  seven taxonomies.
- **A stated trust level per record**, with the reasons, and benchmarks measured by every build, including hold-outs
  by clade distance.
- **Rebuildable**, with inputs checksummed in the manifest.

### Where others remain better

- **Curated protist 18S taxonomy: PR2.** Eukore's genome 18S have no manual curation, and its published 18S carry
  EUKARYOME's labels.
- **ITS species hypotheses: UNITE.** Eukore gives UNITE's codes, it does not replace them.
- **COI breadth: BOLD.** About 1.7 M species against the roughly 214,000 Eukore would take from MIDORI2.
- **Hand-curated phylogenomic sets: PhyloFisher** (240 genes with their paralogs told apart by hand). Eukore's
  family markers (actin, β-tubulin, HSP90, GAPDH, H3) are a genome's best-scoring member, not a chosen ortholog, and
  TEF1 and EF-2 behave the same way in practice (finding L7).
- **A public record.** Every resource above is published, versioned and citable. Eukore has no release, and its
  public repository's indexed description looked stale to the survey.

### One thing to check

The survey reports that BUSCO now ships eukaryote datasets as "odb12.2" (OrthoDB 12.2). Eukore pins
`eukaryota_odb12` (`eukore/odb.py:24-25`; the first build's manifest dates it 2026-05-22). Whether the newer set
differs for the eukaryote root is not verified.

## 3. What "validated" rests on

Measured values are from the first full build's `benchmarks/` (format 6).

| records | check or benchmark | measured | what it does not cover |
|---|---|---|---|
| BUSCO genes | calls against OrthoDB's own models, 836 dataset species | 99.1% called, 92.0% concordant | circular: the species' own proteins were queries |
| BUSCO genes | the same with the species' proteins held out | 99.0% called; 82.9% concordant at the benchmark's threshold, 58.7% at 0.98, 37.4% identical (recomputed here) | congeners' proteins and the HMM stay; 39 protist genomes |
| BUSCO genes | clade hold-out, by the nearest dataset species kept | concordant: genus 83.7%, family 86.9%, order 83.0%, class 80.7%, phylum 62.4%, none 34.1% | good assemblies of dataset species only: no skims, no transcriptomes |
| BUSCO genes | calls laid on the genome's own annotation, 7,541 genomes | introns exact 91.5%; coding bases recall 94.0%, precision 96.6% | measures miniprot's model, not the annotated model that is shipped where one exists |
| BUSCO genes | tiers | `high` 61.7%, `unverified` 6.1%, `medium` 32.0%, `low` 0.1% | no truth set; the own-species check cannot see a contaminant from the same class |
| classic markers | calls against UniProt's protein of the species | RPB1 1,738 of 1,776; RPB2 1,137 of 1,151; EF-2 1,353 of 1,452; actin 953 of 1,638; β-tubulin 1,607 of 2,410 | partly circular (members are queries); family markers are not orthologs |
| rRNA | completeness; tiers | 95.1% of 2.39 M operons with a whole 28S; `high` 80.8%, `low` 1.3% | no external truth; the agreement test against PR2 was not run |
| organelle genes | presence | 13,868 species | no accuracy measure in a build; the nhmmer and miniprot models were measured once on 190 RefSeq genomes |
| taxonomy | names matched in other references | EUKARYOME 31.3%, KSGP 35.2%, PR2 19.5% by name or synonym | matching, not correctness |
| MIDORI2 records (format 7) | the reference screen, on planted errors | 170 of 200 misnamed records dropped; 104 of 200 chimeras caught | within-order errors are invisible to it; never run in a build |
| contamination screens | tiara, FCS-GX on each copy's window | **not run** in the first build (finding B1) | everything |

Two things stand out.

**The `high` tier is where the evidence is, and protists have little of it.** Share of BUSCO copies that are `high`,
from `benchmarks/flags.tsv`: Fungi 77%, Viridiplantae 68%, Vertebrata 59%, Metazoa 54%, Discoba 51%, Stramenopiles
48%, Amoebozoa 31%, Alveolata 25%, Rhodophyta 22%, Haptophyta 18%, Metamonada 10%, Cryptophyceae 9%, Rhizaria 2%.
The held-out concordance tells the same story: Alveolata 46%, Amoebozoa 49%, Stramenopiles 59%, on 6 to 15 benchmark
genomes each.

**Several record types have no measurement.** The 99% collapses, the length check, the owner decisions for shared
sequences, the genome tiers, Eukore's own query proteins, joined split copies, and all of format 7's imports
(EUKARYOME's SSU set, MIDORI2's genes, the EUKARYOME tree) have unit tests and one-off trials, but no figure from a
build.

## 4. Findings

File and line numbers are at `e1ee769`. "First build" is the Eukore_2026.11 build of 2026-11-14.

### A. Would put wrong, missing or mislabelled records into the next build

Most are in code added after the first build (genome tiers, candidate organelle contigs, split-copy joins, format 7).

**Genome set and tiers**

| id | finding | where | status |
|---|---|---|---|
| G1 | **The genome tier `contaminated` is never assigned.** The flag is the text `excluded, the species' only assemblies`, which contains a comma. The package joins a genome's flags with commas and the tier code splits them on commas, so the flag never matches. A species whose every assembly NCBI excludes (256 in the first build) keeps its best one; its BUSCO copies ship as an ordinary genome's, it counts as a reference for its group's expected BUSCOs, and it can supply Eukore's own query proteins. | `eukore/genomes.py:97`, `eukore/package.py:446`, `eukore/tiers.py:60-62`, `:106` | confirmed by probe: `tier_of` returns `genome` |
| G2 | **The tier `partial` drops whole genomes.** A genome flagged `small for its clade` loses all its BUSCO copies whatever its completeness, because the flag is tested first. In the first build *Setaria italica* and *Spirodela polyrhiza* (both chromosome-level) carry the flag and hold 138 of 138 loci. The stage audit counts 52 of 302 such representatives at 90% complete or more; most flags come from a family or order median (Poaceae: anything under about 500 Mb). | `eukore/tiers.py:108`, `eukore/genomes.py:133-145` | confirmed for the two species; the count is reported |
| G3 | **A size filter removes reference genomes before the ranking.** An assembly under 30% of its species' median assembly size is dropped when the species has three or more. *Cygnus atratus*: the first build's representative is an 8.8 Gb contig-level assembly of 659,215 contigs, flagged `genome length too large`; the stage audit reports that the species' 1.12 Gb chromosome-level RefSeq reference was removed because nine oversized assemblies set the median. *Saccharum spontaneum* likewise. Dropped rows also vanish from the species' other assemblies and from the selection tables. | `eukore/sources.py:39-42`, `:83-90` | representative confirmed; the removed reference is reported |
| G4 | **A bare "<taxon> sp." taxid is treated as a species.** "Eukaryota sp." holds 200 unrelated assemblies under one representative in the first build; 164 representatives have a name ending in "sp.". The placeholder test covers only uncultured, unidentified, unclassified, environmental and metagenome names. | `eukore/ncbi.py:17` | confirmed in `genomes.tsv` |
| G5 | **Dropped-tier genomes read as empty genomes.** They stay in the list of genomes searched, so `copy_numbers.tsv.gz` gives 0 copies of every locus, `busco/loci.tsv` call rates keep them in the denominator, and `completeness.tsv` gives them 0, while the headline median leaves them out. | `eukore/package.py:220-233`, `:686-698`; `eukore/benchmarks.py:1036-1066` | confirmed for the call rates; the rest reported |
| G6 | **rRNA and organelle records of dropped, draft or excluded genomes carry no reason for it.** The genome tier reaches the BUSCO copies only. Such a record can be `high`, and it is the reference that published reads are checked against. | `eukore/package.py:130-133`; `eukore/rrna.py:1566-1616` | reported |
| G7 | **Transcriptome gaps.** Crosses of two species are not excluded (8 taken); 37 sets with no base count pass the 10 Mb floor; a species whose genome is dropped by its tier cannot fall back on its transcriptome. | `eukore/genomes.py:332-338` | reported |

**Gene calls**

| id | finding | where | status |
|---|---|---|---|
| L1 | **Most `Fragmented` calls in protists are whole genes.** A call is `Fragmented` under 80% of the profile's length. Where a lineage's protein is simply shorter than the profile, the whole gene fails that. In the first build's OrthoDB table, 1,410 of 1,707 Fragmented calls hold at least 90% of OrthoDB's own protein at 90% of its length, and 1,200 are identical to it. Share of called loci that are Fragmented: Alveolata 14.7%, Haptophyta 13.3%, Discoba 13.1%, against Metazoa and Viridiplantae 0.8%. Effect: protist completeness is undercounted (part of the "ceiling" of the first run's evaluation), and these copies are `medium`, left out of the intron statistics and of the query proteins. The classic markers already normalise by their members' coverage (`hmm_spans`); the BUSCOs do not. This is not the domain-threshold proposal that the audit of 2026-11-10 raised and simulation refuted (that one moved copies toward Fragmented; this moves them away). It would be settled by judging a copy against the profile coverage of its nearest dataset proteins. | `eukore/loci.py:79`, `:702-707` | confirmed by recomputing from `orthodb.tsv` |
| L2 | **A joined split copy ships residues twice where its pieces overlap.** Two contigs sharing 27 bases give a 309-residue protein from a 300-residue gene, `cds_check` clean, status Single. Overlapping contig ends are common in short-read assemblies. | `eukore/loci.py:285`, `:316`, `:660`, `:1168` | reported (probe with real miniprot) |
| L3 | **A joined copy's 100 N spacer is reported as an intron.** `introns.of_copy` takes every space between two exons as an intron and knows nothing of a `join`. The row (length 100, sites NN, class other) enters `introns.tsv`, the conserved intron sites and the genome's intron counts. | `eukore/introns.py:52-67` | confirmed by reading |
| L4 | **`join_split` cannot join pieces on contigs that are neighbours in the FASTA.** miniprot reports one piece when two contigs lie close in the file, so scaffold-ordered skims benefit least, and the cut-genome trial (cut in order) understates what a join can do. | `eukore/loci.py:603-605` | reported on synthetic genomes; real effect suspected |
| L5 | **The rescue pass reads seven nuclear codes as the standard code.** esl-translate accepts tables 15 and 26 to 31 without error and translates by table 1. No wrong sequence results (miniprot handles all eleven tables), but the rescue is weaker in *Blepharisma*, *Condylostoma*, *Mesodinium*, peritrichs, karyorelictids, *Pachysolen* and *Blastocrithidia*. | `eukore/loci.py:837` | reported |
| L6 | **`protein_length` is one too long for miniprot's models** (the stop's `*` is counted), exact for annotated models. It feeds the length check and `median_protein_length`. | `eukore/genomes.py:1274`, `eukore/loci.py:187` | reported |
| L7 | **TEF1 and EF-2 are labelled `clean` but behave as families.** TEF1 is Duplicated in 356 of 745 benchmark genomes, EF-2 in 334 of 1,452, so their copy is a best member, as for actin and β-tubulin. | `data/markers.tsv` | reported |

**Organelle genes**

| id | finding | where | status |
|---|---|---|---|
| O1 | **Candidate contigs can ship other sequences as organelle genes.** New since the first build: in an assembly with no titled organelle molecule, contigs are nominated by size and GC. A plastid candidate is any contig of 100 to 250 kb, with no GC test. barrnap with bacterial models then fills every rRNA that nhmmer missed, so a bacterial contaminant's 16S can become a plastid's or a mitochondrion's. One contig can be nominated for both organelles and give two records of one name. `genes.tsv` has no column telling a candidate contig from a titled molecule. | `eukore/organelle.py:671-746`, `:1051-1065` | plastid rule confirmed by reading; the double record reported; frequency suspected |
| O2 | **A gene with one long intron is shipped as if it crossed the origin.** `span_of` takes the largest space between a feature's parts as the genome's origin whenever the rest of the circle is shorter. A hexacoral nad5 annotated `join(2434..3153,15087..16202)` on an 18 kb genome ships as 15087 to the end plus 1 to 3153: exon 2, the genes between, exon 1. 129 RefSeq genes are affected (114 nad5). The CDS and protein are right; the gene record is not. | `eukore/organelle.py:226-236` | confirmed by reading; the count is reported |
| O3 | **Organelle genes are judged by the nuclear coding rule.** The 18S got its own own-species rule on 2026-11-10; organelle genes still go through `qc.verdict`, where a species alone in its family whose gene is 90% like another class's is called foreign and `low`. Plastid and mitochondrial rRNA and psbA are that conserved across classes. | `eukore/organelle.py:1429`, `eukore/qc.py:169-180` | rule reported by probe; extent not measured |
| O4 | **Fragmented organelle rRNA ships as a whole gene.** Where RefSeq annotates an rRNA in many pieces (apicomplexans, green algae), the longest piece is taken as the gene: 196 SSU or LSU have three or more features and are not marked partial; some are under 300 bases. They then serve as queries for related assemblies. | `eukore/organelle.py:462`, `:757` | reported on the real RefSeq parse |
| O5 | **MIDORI2 coding records take the first open frame.** 7,860 COX1 and 235 CYTB have several; by protein similarity the first is wrong in about 110. The nucleotide record is unaffected. | `eukore/organelle.py:1271` | reported on GB273 |

**Taxonomy and labels**

| id | finding | where | status |
|---|---|---|---|
| T1 | **The EUKARYOME lineage depends on the order of work.** The audit of published reads places taxa in the tree before `Tree.learn` runs in the package step. A clade placed early under its NCBI name is then not learned: copepods stay under a second class `Copepoda` beside EUKARYOME's `Hexanauplia`. The figures in `docs/methods.md` for the one-tree taxonomy were measured without a build, so they do not include this. | `eukore/audit.py:545-551` before `eukore/package.py:844-845` | order confirmed by reading; effect reported at unit level |
| T2 | **The genus homonym guard knows three kingdoms.** A genus match is refused only when both lineages name animals, fungi or plants. The bird genus *Lessonia* takes the kelp's PR2 lineage. In the first build: PR2 32 species, SILVA 24, KSGP 18, UNITE 14. | `eukore/labels.py:137-146`, `:190-195` | confirmed by reading; counts reported |

**rRNA and the audit of published reads**

| id | finding | where | status |
|---|---|---|---|
| R1 | **The window chimera screen's references are not reduced to one per genome.** Every genome 18S and operon record is a reference: 4.43 M in the first build, where 8,140 of 22,740 genomes hold 31 or more. A probe with real minimap2 shows two effects: a window identical to a sequence present some 600 times gets no hit at all (the preset skips seeds seen over 500 times), and at 40 copies all 30 secondary hits are one genome's, so the second genome the rule asks for never appears. A chimera or mislabel toward a well-sampled genome then passes unflagged. The calibration in `docs/methods.md` used 5,763 references. | `eukore/audit.py:578-581`; `eukore/chimeras.py:85-86`, `:156-211` | reference set confirmed by reading; minimap2's behaviour reported by probe; size in a build suspected |
| R2 | **UCHIME's result depends on how the work is split, and its kept parts are never reused.** A part's queries are removed from its database, so records screened in the same part cannot be each other's parents: with three parts a third of the published 18S cannot be. A probe with vsearch: a true chimera scores 35 alone and 2.5 with one parent in its part. The queries' names (`u<N>`) follow the iteration order of a Python set of strings, which differs in every process, so the cache keys differ in every run. | `eukore/chimeras.py:327-343`, `:403-413`; `eukore/audit.py:576`, `:595` | removal and naming confirmed by reading; scores reported by probe |
| R3 | **One assembly counts as several genomes.** In the window screen a genome's own record (`GCA_x.1`), EUKARYOME's copy (`EUK_GCA_x`) and ROD's (`GCF_x.2`) get three ids, so "seen in a second genome" is met by one assembly. In the part screen, EUKARYOME's further copies (`GCA_x_dup1`) fail the assembly pattern, so a contaminant operon is confirmed by its own assembly. `screen.genome_of` already has the right rule; `audit.genome_id` does not use it. | `eukore/audit.py:563-566`; `eukore/screen.py:221-233` | `audit.genome_id` confirmed by reading; the `_dup` case reported by probe |
| R4 | **A read labelled above every group can make a genome's 18S foreign.** Format 7 labels a read of no named lineage as Eukaryota. A genome 18S whose nearest hit at 97% is such a read, with its class-mates under 90%, is called another clade's and `low`. | `eukore/rrna.py:1196-1207`, `:1273-1282`; `eukore/qc.py:116-120` | reported on a toy taxonomy; frequency unknown |
| R5 | **Re-extraction of joined 18S can duplicate records**, and its cache does not know whether ITSx was present: without it, units with introns silently keep barrnap's ITS borders. | `eukore/rrna.py:591-593`, `:646-685` | reported with a stand-in barrnap; the ITSx case suspected |
| R6 | **A chimera's kept part can hold foreign bases.** The cut keeps what lies outside windows that stand foreign at 97%; foreign windows at 95 to 97% count as the record's own. | `eukore/chimeras.py:228-236`, `:271-281` | suspected |
| R7 | **In-silico PCR does not tell "primer site not in the record" from "primer does not match".** A record that starts after a site is `amplifies f`. NS1/NS4 at 50% against V9 at 92% probably reflects truncated records. This bears on the primer coverage figures. | `eukore/primers.py:54-90` | suspected |

**Build history**

| id | finding | where | status |
|---|---|---|---|
| H1 | **Eukore's own query proteins make a package depend on how it was built.** The package step writes the query set; any later process run reruns every genome of the served groups that is still on disk and was called with another set. Genomes already pruned keep their old calls, so the mix follows what happened to be on disk. Each rerun changes proteins, hence the set's digest, hence the next round of reruns. Neither the manifest nor `genomes.tsv` records which set made a genome's calls. | `eukore_build.py:303-331`, `:516-522`, `:2090-2108`; `eukore/genomes.py:1359` | mechanism confirmed by reading; that it never settles is suspected |
| H2 | **A species whose planned assembly changed is absent until the new one is done.** After a refresh or a TOGA2 swap, the package reads the new plan only; the old results in `WORK/done` are not used. Only the manifest's `genomes_not_done` count shows it. | `eukore/genomes.py:1688-1713`, `eukore/package.py:433-436` | reported |
| H3 | **A changed published read can be masked by a stale cache line.** Region files merge by file name and the digest filter runs after the merge. The read then ships whole as an 18S with no regions. Latent while EUKARYOME and ROD stay pinned. | `eukore_build.py:995-1004` | reported by probe |

### B. Bears on the truth of a published claim

| id | finding | where | status |
|---|---|---|---|
| B1 | **The contamination screens did not run in the first full build, and nothing in the package says so.** The log reads "no tiara: windows not classified"; no copy among 4.33 M carries a tiara or FCS-GX flag. A copy with no screen result gets no reason and can be `high`. The README's checks row and `benchmarks/databases.tsv` list both screens as Eukore's. | `eukore/package.py:123-127`, `:153-155` | confirmed in the log and `flags.tsv` |
| B2 | **The manifest cannot show an absent optional input.** No field for tiara, FCS-GX, the VEuPathDB key, which tools were present, `--loci`, `--keep_poor_genomes`, or published batches left out. A missing reference is only a missing key. A run report section on missing inputs, as the first run's evaluation proposed, cannot be built from it. | `eukore_build.py:2462-2477`, `eukore/package.py:916` | reported |
| B3 | **"OrthoDB's gene model" overstates the agreement.** Concordant means 90% of 5-mers shared both ways; a 400-residue protein can lose a 35-residue exon and pass. Held-out: 82.9% at that threshold, 58.7% at 0.98, 37.4% identical. The honest wording is "agrees with OrthoDB's protein at 90%". | `eukore/benchmarks.py`; `docs/benchmarks.md` | confirmed by recomputing |
| B4 | **"Completeness" is not BUSCO's.** The status is neither BUSCO's length rule nor compleasm's coverage (see L1), so the figures are not comparable with a BUSCO C% reported elsewhere. | `docs/methods.md:8`, README | reported |
| B5 | **"With a genome" counts transcriptomes.** 5,028 of the 33,664 species "with a genome" are transcriptome sets; the selection table calls them "chosen". | `eukore/genomes.py:612-615` | confirmed in `species.tsv` |
| B6 | **The comparisons with other databases favour Eukore.** Every `species.tsv` row counts as Eukore's, so in format 7 the roughly 291,000 species, mostly MIDORI2 and EUKARYOME imports, inflate "Eukore only" against EukDetect2 and MetaPhlAn, and overlap with EUKARYOME and MIDORI2 is by construction. Marker reads are simulated from Eukore's own record spans. | `benchmarks/eukore_bench/databases.py:69`, `compare.py:62-66`, `reads.py:411-416`; `eukore/benchmarks.py:785` | reported |
| B7 | **The share of species "named by name or synonym" omits one class from its denominator** (`upward`, new in format 7), which inflates EUKARYOME's figure. | `eukore/benchmarks.py:1209-1210` | reported |
| B8 | **The 18S agreement test against PR2 can call a correct record a disagreement.** minimap2 keeps 50 hits; in a species-rich clade the species' own 18S is not among them. Its 97% threshold also cannot tell species apart (nearest congener's 18S: median 99.7%). | `benchmarks/eukore_bench/databases.py:212-255` | suspected |
| B9 | **FCS-GX's absolute threshold demotes good references.** 100 kb of another eukaryote flags an assembly whatever its size; 176 species lose their RefSeq reference for this flag alone, 83 of them with under 0.1% contaminant. It is the documented rule; at genome scale it is a questionable one. | `eukore/quality.py:22`, `:102` | reported |
| B10 | **Re-measuring a package breaks its check.** `python3 -m eukore.benchmarks PACKAGE`, which the run report advises, rewrites files the manifest lists; `package.check` then fails. `package.check` is called by neither the build nor the release. | `eukore/package.py:909-914`, `eukore/run_report.py:316` | reported |
| B11 | **The website's "rRNA records alone" tile** counts every species without a representative genome, which in format 7 includes the organelle-only species. | `eukore/release.py:334` | reported |

### C. Robustness and speed

| id | finding | where | status |
|---|---|---|---|
| C1 | **The straggler takeover added after the first build cannot fire.** A claim's age is read from the file's change time, which the 60-second heartbeat resets, so the 900-second threshold is never reached. The 1 h 27 m stall of the first build is not fixed. If it did fire on a minimap2 or UCHIME part it would raise an uncaught error. | `eukore/qc.py:538-546`, `:565`, `:661-676`; `eukore/genomes.py:995-999` | confirmed by reading |
| C2 | **Benchmark searches are not kept across reruns.** Three searches share one folder and each call's cleanup deletes the previous call's hits, so a rerun repeats the 18S search behind the 3 h 36 m step. | `eukore/benchmarks.py:1146-1154`, `eukore/qc.py:693-699` | reported by probe |
| C3 | **The final job's backfill skips genomes whose claim still looks live.** Array tasks killed at the time limit leave claims that read as live for 300 seconds; the backfill returns and is not repeated. Those genomes miss their retried rRNA, tiara classes and hold-outs. | `eukore/genomes.py:1407-1409`, `eukore_build.py:2256` | logic reported; timing suspected |
| C4 | **Prune deletes genomes that still need a retry** (`rrna_failed`, `metadata_failed`), so the retry never happens. None in the first build. | `eukore/genomes.py:1829-1836` | reported |
| C5 | **Caches are not keyed by tool presence or by their inputs.** Published reads' regions and audits, and the joined rRNA, are kept whether or not ITSx, nhmmer or UniVec were there on the node that made them, and whatever `data/primers.tsv` was then. The manifest stamps the packaging node's tools, not the auditing nodes'. | `eukore_build.py:1004-1106`; `eukore/rrna.py:651-658`; `eukore/audit.py:367`, `:435-446` | suspected |
| C6 | **The store's change lock has the takeover race the genome claims had**, and no heartbeat. One change can be lost; rare. | `eukore/store.py:266-275` | reported |
| C7 | **A source list that fails to parse is silent.** A transcriptome list fetched on 2026-11-10 has 10,444 rows and no FASTA path: no set is read and no warning given. Entries in `WORK/unavailable.tsv` never expire. | `eukore/genomes.py` (`read_tsa`, `download_of`) | reported |
| C8 | **The reference screen on all 15 MIDORI2 files is unmeasured**, and opens one file per own-clade group at once (at least 223): at risk under a limit of 1,024 open files. So is the organelle check's COX1 search over some 230,000 records. | `eukore/refscreen.py:551`, `eukore/labels.py:437` | suspected |

## 5. What is sound

Checked and found right, with the bearing each has on trusting a package.

- **Coordinates and codes.** A probe with real miniprot on both strands, with introns of every phase and a
  frameshift, in all eleven nuclear tables: the copy is the gene, the package's CDS is the caller's, intron phases
  are right, the stop codon is inside the copy. One genetic code goes to miniprot, to the CDS check and to the
  package, and the shipped protein is always the translation of the shipped CDS.
- **Best-profile rule.** Enforced in every path that makes a call (first pass, rescue, recovery, join, annotated
  model); the markers' decoys compete.
- **Hold-outs.** A species' own dataset proteins, its strains' and Eukore's own queries are all left out; no
  annotated model is swapped in. Ancestral variants stay, as the docs say.
- **Organelle parsing.** On the real RefSeq release, whole CDS with a stop inside: animals 0 of 93,198, fungi 0 of
  4,243; land plants 346 of 58,624, which is RNA editing and rightly `medium`. So `join`, `complement`,
  `/codon_start`, `/transl_table` and `/transl_except` are handled.
- **rRNA.** Units on the plus and minus strand give identical records, regions and introns, and the slices match
  their source. The 99% typing's two pre-filters are strict necessary conditions: 60 random trials equal the
  exhaustive comparison. In-silico PCR handles IUPAC codes on both sides, the 3' rule and the reverse primer's
  orientation correctly. The V regions are computed after all trimming, so a cut cannot shift them. The search
  caches are keyed by parameters, queries and database, so a part made for another record set cannot be reused.
- **One job per species**, ranking as documented with ties broken by accession, merged taxids followed, RefSeq and
  GenBank copies paired without loss.
- **Assembly metadata.** Identical results for plain and gzipped input, any chunk size, lower case and N runs across
  chunk bounds.
- **Done genomes and claims.** A genome's results are written whole and renamed; a damaged file reads as not done;
  claims are exclusive with a locked takeover.
- **The package.** Built in a partial folder, manifest last, then renamed. Manifest counts equal the tables (4,329,328
  copies by status and by tier). `check(deep=True)` covers every file and BGZF end block. `docs/FORMAT.md`'s columns
  equal the code's for the eight main tables.
- **Failures surface.** Tool exit codes, download status and gzip ends are checked; a dead worker stops the build.
- **Tests.** 486 pass. With real tools the format 7 audit tests (flanks, UCHIME) pass too.

## 6. Documentation against the code

- `docs/methods.md`, genome set: assemblies NCBI flags as contaminated, mostly foreign or failing ANI are "never"
  taken; since the first build's fixes a species' only such assembly is kept. Genome tiers, rRNA donors and the
  transcriptome floor are not described; transcriptomes are said to give no rRNA, and now do.
- `docs/methods.md:118`: "from first to last codon" is the chosen alignment's first to last aligned codon; there is
  no search for a start codon.
- `docs/methods.md:317`: MIDORI2 is still "COI, given by the user".
- `docs/methods.md` and `docs/FORMAT.md`, organelle genes: described as from titled molecules; candidate contigs
  appear only under Changes.
- `docs/FORMAT.md`, first paragraph: "its taxonomy is NCBI's"; `docs/methods.md` makes EUKARYOME the primary
  taxonomy. The same paragraph lists the organelle genes of format 4.
- `docs/FORMAT.md:45`: `high` means agreement "to the family or deeper"; the code also passes a lone lineage at its
  order (`eukore/qc.py:162-168`).
- `docs/FORMAT.md`: tier `contaminated` is described and never assigned (G1); `protein_length` (L6).
- `docs/FORMAT.md`, `rrna/`: UCHIME is said to screen each record "against all others" (it is against the others
  less the record's own part, R2); a chimera's `kept a-b` is in the trimmed record's coordinates, not the published
  record's as stated; `lsu_complete` is "the whole model" (the code takes 80%); the flags `mislabel` and
  `chimera_weak` are described by phylum (the code uses EUKARYOME's kingdom).
- README (working tree), "Left out: organelle rRNA among the nuclear records; prokaryotes": only contigs titled as
  an organelle are left out; other such 18S are kept and flagged `ssu_prokaryote`.
- `docs/building.md:92` "about 80 tests", README "about 140": 486.
- `docs/building.md`: only `--species` marks a package partial, `--loci` does not; the options every stage must
  share omit several that change the plan; "no task list is written" although `WORK/plan/<array>.genomes.json` is.
- README: states terms for EukDetect2's records only. MIDORI2's, EUKARYOME's and ROD's sequences are shipped under
  "their own terms", unnamed; the survey found CC BY 4.0 for EUKARYOME and MIDORI2 and no stated licence for ROD.
  The manifest has no licence field.
- `data/model_species.tsv`: NCBI's zebrafish reference is now another assembly; the pig's listed accession is no
  longer in the GenBank table.

## 7. What the claim needs

In order. The first two are prerequisites for the rest.

1. **Fix what would corrupt the next build**: G1 to G5, L1 to L3, O1 to O3, R1 to R4, T1 and T2. Most are small in
   code. L1 needs one measurement first (status judged against the nearest dataset proteins' profile coverage, on
   the first build's OrthoDB table); R1 needs the references reduced to one 18S type per genome and the screen's
   calibration repeated on that set.
2. **Build format 7 in full**, with the inputs the first build lacked (VEuPathDB key, a current taxdump, the 805
   annotations) and with tiara and FCS-GX actually run. Record every optional input and tool in the manifest (B2),
   so that a package states what checked it.
3. **Measure what has no measurement** (section 3):
   - organelle genes against RefSeq's annotation, held out by species, in the build;
   - genome 18S against the same species' PR2 or EUKARYOME record, with the 50-hit limit fixed (B8);
   - the genome tiers and joined copies against the dataset species;
   - the tiers themselves against planted errors, as the reference screen was tested: contaminants at several
     taxonomic distances, to state what share `low` catches and what `high` wrongly keeps.
4. **Widen the protist evidence.** VEuPathDB's annotations would raise the 39 protist benchmark genomes; report the
   hold-out by supergroup beside the overall figure.
5. **Word the claims as measured**: "agrees with OrthoDB's protein at 90%" (B3), completeness as Eukore's own
   measure (B4), genomes and transcriptomes counted apart (B5), comparisons with other databases on genome-derived
   species only (B6).
6. **Publish a release** with a version, a DOI and the sources' licences named.

With 1 to 3 done, the statement "the most complete unified reference of eukaryote marker genes, with every record
checked and its trust level stated" is defensible. "Best" for a single marker stays with the specialist databases
named in section 2, which Eukore links to rather than replaces.

## 8. Not covered

- No build was run; findings in code newer than the first build are from reading and probes.
- `benchmarks/eukore_bench/figures.py`, the website's PHP, `tools/` other than `eukore.sh`, and the fetchers'
  network paths were not audited.
- The thresholds of the reference screen and of the window chimera screen were read, not re-derived.
- Figures about other databases are as the survey read them on 2026-11-15; several could not be confirmed from a
  primary source (MIDORI2's current release, MetaPhlAn's eukaryote count from its own files, BUSCO's odb12.2 size,
  licences of SILVA, UNITE, PR2 and BOLD).

## 9. Fixes made (branch `audit-fixes`, 2026-11-15)

The findings I had re-checked (marked confirmed in section 4) and that would change the next build were fixed the
same day, on branch `audit-fixes` from `e1ee769`, each with a test (a test of the fix's mechanism where the build
path itself is not under test). All 495 tests pass on Python 3.12 and on the build environment's 3.14 with real
miniprot, HMMER, barrnap and vsearch. Nothing has been built with them.

| id | fix | how it was checked |
|---|---|---|
| G1 | the flag `excluded (the species' only assemblies)` without a comma; the old text still read | tier `contaminated` given (test) |
| G2 | `partial` only for a genome under half its group's expected BUSCOs | a "small" genome of 100% is `genome` (test) |
| G3 | the size median of a species' unflagged assemblies; NCBI's reference or representative genome never left out | the *Cygnus atratus* case (test) |
| G4 | a bare "<taxon> sp." under a clade above the genus is a placeholder; a genus' "X sp." stays | test |
| G5 | a dropped genome is searched for nothing: out of `copy_numbers.tsv.gz`, the call rates and `composition.tsv`'s medians | test |
| L1 | a call complete against its query protein's own HMM coverage (at least 0.5); genome record 4 | measured first: on the 7-clade trial 267 of 385 Fragmented calls turn Single, of the 34 OrthoDB annotates 32 are its protein, 2 not; none of the 4 that stay is |
| L3 | the space between a joined copy's pieces is no intron | test |
| O1 | a candidate contig gives genes only beside a protein-coding gene of its organelle; no barrnap on candidates; one organelle per contig; flag `candidate_contig` | test with stand-in searches |
| O2 | a gene read along its molecule spans its first base to its last; organelle cache 5 | the hexacoral nad5 case and an origin crossing (test) |
| T1 | the tree learns before the audit places the reads, and learning forgets earlier placements' NCBI-named nodes | the copepod case; the test fails without the fix |
| T2 | the genus homonym guard covers ten deep protist lineages beside the three kingdoms, on the species' whole NCBI lineage | the *Lessonia* case (test); not re-counted on the first build |
| R1 | the window screen's references: 3 distinct 18S per genome | test of the selection; the screen's calibration on the build's references is not repeated |
| R2 | UCHIME against the whole database with `--self`; names in a fixed order; `UCHIME_RECORD` 2 | vsearch 2.32 honours `--self` in `--uchime_ref` (probe); a chimera of two co-screened records is found (test) |
| R3 | one id for an assembly's every copy, in both screens | test |
| H1 | own queries taken once, recorded in the manifest (`inputs.queries`) | the rule itself has no test |
| C1 | the straggler timed from when this node first saw the claim; `search_too` refuses window-screen and UCHIME parts | a claim kept alive by a heartbeat is searched here too (test) |

Consequences for the next build: every done genome is run again (genome record 4; record 3 already required it),
RefSeq's organelle release is parsed again and every assembly's organelle search made again (organelle cache 5), and
UCHIME's kept parts are searched again. `docs/FORMAT.md` (Changes), `docs/methods.md` and `docs/building.md` describe
the changes.

**Still open.** These findings are not fixed:

- G6, G7, L2, L4 to L7, O3 to O5, R4 to R7, H2, H3: reported by a stage audit and not re-checked by me, or needing a
  decision first:
  - O3 needs an organelle rule of its own;
  - L2 needs a decision on how overlapping pieces are trimmed;
  - R7 needs primer sites judged on the record's span.
- Section B (the wording of published claims, the manifest's missing-input fields) and section C other than C1.
- Section 7, items 2 to 6.
