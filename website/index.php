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
        with a genome at NCBI, labelled by its taxonomy and given a confidence tier.</p>
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
    <a href="main.php?site=stats">its statistics</a>.</p>
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
        <li>its rDNA: 18S spans and whole operons (18S, ITS1, 5.8S, ITS2, 28S, 5S) with the 18S's V1–V9, beside
          EUKARYOME's and ROD's long reads;</li>
        <li>its mitochondrial and plastid genes: 12S, 16S, COX1, CYTB, nad1, nad5, atp6; rbcL, matK, tufA, psbA,
          ycf1 and the trnH-psbA spacer; in-silico PCR of 18 barcode primer pairs;</li>
        <li>its genome's metadata (size, contiguity, GC, masking, telomeres, a sourmash sketch), its genetic code
          checked, its introns' phases and splice sites;</li>
        <li>its lineage in NCBI and in PR2, EUKARYOME, KSGP, SILVA, UNITE and MIDORI2, with Open Tree, GBIF, WoRMS
          and BOLD ids;</li>
        <li>checks across species on every copy (own-species agreement, contamination screens, length outliers) and a
          <b>high / medium / low</b> confidence tier with its reasons: flagged, never silently dropped.</li>
      </ul>
    </div>
  </div>

  <div class="split">
    <h2 class="split-title">Why?</h2>
    <div class="split-body">
      <ul>
        <li>Metagenome profilers need eukaryote marker genes that are <b>complete</b>: protal reads Eukore as its
          eukaryote reference.</li>
        <li>Among congeners, BUSCO genes tell species apart where the 18S cannot: on a trial build, two thirds of their
          copies were under 99% identical to every other species, against a third of the 18S copies.</li>
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
          EukProt; BUSCO's eukaryota_odb12; EUKARYOME and ROD; RefSeq's organelle genomes; six taxonomic references</li>
        <li><span class="flow-step">Fetch</span>one best assembly per species: NCBI's contamination and ANI flags
          excluded, then RefSeq, level, contiguity</li>
        <li><span class="flow-step">Process</span>each genome once: BUSCO calls (miniprot, then BUSCO's HMMs, a rescue
          pass), rDNA (barrnap, ITSx), organelle genes, metadata and a sketch</li>
        <li><span class="flow-step">Package</span>checks across species and confidence tiers, introns and genetic
          codes, 18S regions and in-silico PCR, names and ids, benchmarks</li>
        <li><span class="flow-step">Release</span>this site's statistics and the downloads, BGZF-compressed</li>
      </ol>
      <figure class="figure">
        <a href="images/overview_slide.png"><img src="images/overview_slide.png" alt="Eukore at a glance: its 21 sources, the build's fetch, process and package stages, and the package's parts" loading="lazy"></a>
        <figcaption>Eukore at a glance: its sources, the build and what the package holds. Every source and step in
          detail: <a href="images/overview.png">the full figure</a>.</figcaption>
      </figure>
    </div>
  </div>
</div>
<?php require __DIR__ . '/inc/footer.php'; ?>
