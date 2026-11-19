<?php
require __DIR__ . '/config.php';
require __DIR__ . '/inc/data.php';
$page = 'home';
$stats = release_json('stats.json');
$totals = $stats['totals'] ?? null;
require __DIR__ . '/inc/header.php';
?>
<section class="cover">
  <div class="container cover-row">
    <div class="cover-text">
      <h1>Every eukaryote genome's marker genes, checked and complete</h1>
      <p>BUSCO genes with all their introns, whole rDNA operons and organelle barcode genes for each eukaryote species
        with a genome at NCBI, beside the audited 18S of EUKARYOME and the screened mitochondrial barcodes of MIDORI2:
        every sequence in NCBI's and EUKARYOME's taxonomy, with a confidence tier.</p>
    </div>
    <div class="cover-actions">
      <a class="button-boxed" href="main.php?site=tutorial">Get started!</a>
      <p class="cover-or">or</p>
      <div class="cover-pair">
        <a class="button-dashed" href="main.php?site=downloads">Download Eukore</a>
        <a class="button-dashed" href="<?= h(GITHUB_URL) ?>">Go to GitHub</a>
      </div>
    </div>
  </div>
</section>

<div class="container content">
<?php if ($totals): ?>
  <section class="tiles" aria-label="This release">
    <div class="tile"><span class="tile-label">Species</span><span class="tile-value"><?= num($totals['species']) ?></span></div>
    <div class="tile"><span class="tile-label">Genomes and transcriptomes</span><span class="tile-value"><?= num($totals['genomes']) ?></span></div>
    <div class="tile"><span class="tile-label">Marker loci and genes</span><span class="tile-value"><?= num($totals['marker_loci']) ?></span></div>
    <div class="tile"><span class="tile-label">Sequences</span><span class="tile-value"><?= num($totals['sequences']) ?></span></div>
  </section>
  <p class="tiles-note">In <?= release_badge($stats) ?>, built <?= h(date_text($stats['built'] ?? null)) ?>:
    <a href="main.php?site=stats">its statistics</a>.<?php if (!empty($stats['bundle'])): ?>
    This build (Eukore <?= h($stats['eukore'] ?? '') ?>, format <?= h($stats['format'] ?? '') ?>) is under
    evaluation, not a published release: what is described below is the current code's, of which a release will be
    built.<?php endif; ?></p>
<?php endif; ?>

  <div class="split">
    <h2 class="split-title">What?</h2>
    <div class="split-body">
      <p class="lead">Eukore is a database of eukaryote marker genes, and the software that builds it from NCBI's
        assemblies. For each eukaryote species with a genome it holds:</p>
      <ul>
        <li>its copies of BUSCO's 129 near-universal eukaryote genes and of nine classic phylogenetic markers (RPB1,
          RPB2, EF-2, TEF1, β-tubulin, actin, HSP90, GAPDH, H3), <b>each whole from first to last codon, introns and
          all</b>, with its coding sequence and protein;</li>
        <li>its rDNA: 18S spans and whole operons (18S, ITS1, 5.8S, ITS2, 28S, 5S) with the 18S's V1–V9;</li>
        <li>its mitochondrial and plastid genes: 12S, 16S, COX1, COX3, CYTB, nad1, nad5, atp6; rbcL, matK, tufA,
          psbA, ycf1 and the trnH-psbA spacer; in-silico PCR of 18 barcode primer pairs;</li>
        <li>its genome's metadata (size, contiguity, GC, masking, telomeres, a sourmash sketch), its genetic code
          checked, its introns' phases and splice sites;</li>
        <li>its lineage in NCBI's taxonomy and in EUKARYOME's, made one tree free of homonyms; its names in PR2, KSGP,
          SILVA, UNITE and MIDORI2, with Open Tree, GBIF, WoRMS and BOLD ids;</li>
        <li>checks across species on every copy (own-species agreement, contamination screens, length outliers) and a
          <b>high / unverified / medium / low</b> confidence tier with its reasons: flagged, never silently
          dropped.</li>
      </ul>
      <p>Beside the genomes, for the species that have none:</p>
      <ul>
        <li><b>published rRNA</b>: EUKARYOME's long reads and 18S set and ROD's operons, labelled down to the rank
          their label reaches and audited at build time (vectors, barcodes, DNA beside the gene and chimeras cut off
          or flagged);</li>
        <li><b>organelle genes</b> of every species with a RefSeq organelle genome, and MIDORI2's mitochondrial 12S,
          COX1 and CYTB, one record a species, screened for misnamed and chimeric records and for nuclear copies.</li>
      </ul>
    </div>
  </div>

  <div class="split">
    <h2 class="split-title">What not?</h2>
    <div class="split-body">
      <ul>
        <li>An assembly NCBI, its FCS-GX screen or its ANI check finds contaminated or misidentified, where the
          species has a better one; NCBI's placeholder species and crosses.</li>
        <li>The gene copies of a genome that is a part of one or holds under half the BUSCOs its clade should: the
          species stays, with its rRNA and organelle genes.</li>
        <li>The sequences of PR2, KSGP and SILVA: their records are linked to Eukore's and their names given, never
          taken.</li>
        <li>Prokaryote and organelle SSU among the 18S (a contaminant's 16S, a symbiont's in a transcriptome: set
          apart in their own file of the package); MIDORI2's records of the wrong length, with ambiguous bases or
          without a reading frame.</li>
      </ul>
    </div>
  </div>

  <div class="split">
    <h2 class="split-title">Why?</h2>
    <div class="split-body">
      <ul>
        <li>Metagenome profilers need eukaryote marker genes that are <b>complete</b>: protal reads Eukore as its
          only source of eukaryote marker genes, rRNA and names.</li>
<?php
// the build's own figures (stats.json's benchmarks.resolution: each marker's copies under 99% identical to their
// nearest congener's), not numbers written here: they change with every build
$res = $stats['benchmarks']['resolution'] ?? [];
$busco99 = $res['BUSCO (all)']['congeners']['under_99'] ?? null;
$ssu99 = $res['18S']['congeners']['under_99'] ?? null;
$cox99 = $res['COX1']['congeners']['under_99'] ?? null;
?>
        <li>Among congeners, BUSCO genes tell species apart where the 18S cannot<?php if ($busco99 !== null && $ssu99 !== null): ?>:
          in this build, <?= pct($busco99, 0) ?> of their copies are under 99% identical to the nearest congener's,
          against <?= pct($ssu99, 0) ?> of the 18S copies<?php if ($cox99 !== null): ?> (and <?= pct($cox99, 0) ?> of
          the mitochondrial COX1)<?php endif; ?><?php endif; ?>.</li>
        <li>One species' nuclear genes, rRNA and organelle barcodes stand side by side under one name: a link between
          metagenomes, amplicon surveys and barcoding.</li>
        <li>Phylogenetics, primer design and contamination screens get clean, labelled genes of every eukaryote group:
          fungi, protists, plants, invertebrates and vertebrates.</li>
        <li>Every package is plain TSV and BGZF-compressed FASTA with a manifest of checksums, versioned and measured
          by its own benchmarks.</li>
      </ul>
    </div>
  </div>

  <div class="split">
    <h2 class="split-title">Workflow</h2>
    <div class="split-body">
      <ol class="flow">
        <li><span class="flow-step">Sources</span>NCBI's assemblies and taxonomy, INSDC transcriptomes, TOGA2,
          EukProt; BUSCO's eukaryota_odb12; EUKARYOME and ROD; RefSeq's organelle genomes and MIDORI2; six taxonomic
          references</li>
        <li><span class="flow-step">Fetch</span>one best assembly per species: NCBI's contamination and ANI flags
          excluded, then RefSeq, level, contiguity</li>
        <li><span class="flow-step">Process</span>each genome once: BUSCO calls (miniprot, then BUSCO's HMMs, a rescue
          pass, split genes joined), rDNA (barrnap, ITSx), organelle genes, metadata and a sketch; the published
          reads annotated and audited</li>
        <li><span class="flow-step">Package</span>checks across species, confidence tiers and genome tiers, the
          reference screen, introns and genetic codes, 18S regions and in-silico PCR, names and ids, benchmarks</li>
        <li><span class="flow-step">Release</span>this site's statistics and the downloads, BGZF-compressed</li>
      </ol>
      <figure class="figure">
        <a href="images/overview_slide.png"><img src="images/overview_slide.png" alt="Eukore at a glance: its 25 sources, the build's fetch, process and package stages, and the package's parts" loading="lazy"></a>
        <figcaption>Eukore at a glance: its sources, the build and what the package holds. Every source and step in
          detail: <a href="images/overview.png">the full figure</a>.</figcaption>
      </figure>
    </div>
  </div>
</div>
<?php require __DIR__ . '/inc/footer.php'; ?>
