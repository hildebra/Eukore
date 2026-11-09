<?php
// The release's statistics (data/stats.json, taxonomy.json: release.py's): its species by genome source and
// supergroup, its taxonomy as a tree, its marker categories, its headline benchmarks.
$stats = release_json('stats.json');
if (!$stats) {
    echo no_release('the statistics');
    return;
}
$t = $stats['totals'];
$sources = $stats['sources'];
// each source's colour slot, fixed by the source (not by its place among those present)
$slot = ['NCBI RefSeq' => 1, 'NCBI GenBank' => 2, 'TOGA2' => 3, 'ENA TSA' => 4, 'EukProt' => 5,
         'rRNA records only' => 6, 'other' => 0];

function stacked_bar(array $counts, int $scale, string $who, array $slot, array $text): string
{
    $out = '<div class="stack" role="img" aria-label="' . h($who . ': ' . implode(', ', array_map(
        fn($k, $v) => "$k $v", array_keys($counts), array_values($counts)))) . '">';
    foreach ($counts as $source => $n) {
        if ($n <= 0) {
            continue;
        }
        $w = 100 * $n / max(1, $scale);
        $out .= '<span class="seg s' . ($slot[$source] ?? 0) . '" style="width:' . number_format($w, 3, '.', '')
            . '%" data-tip="' . h($who . ' · ' . $source . ': ' . num($n) . ' species') . '"></span>';
    }
    return $out . '</div>';
}

$scale = 0;
foreach ($sources['supergroups'] as $sg) {
    $scale = max($scale, $sg['species']);
}
$categories = $stats['markers'] ?? [];
?>
<p class="page-intro">The statistics of <?= release_badge($stats) ?>, built <?= h(date_text($stats['built'])) ?> by
  Eukore <?= h($stats['eukore'] ?? '') ?> (package format <?= h($stats['format'] ?? '') ?>).
<?php if (!empty($stats['partial'])): ?>
  This is a <b>trial package</b>: it holds the species of a few clades only, not every eukaryote.
<?php endif; ?>
  Every number here is written by the build itself (<code>OUT/release/data/</code>).</p>

<section class="tiles" aria-label="Totals">
  <div class="tile"><span class="tile-label">Species</span><span class="tile-value"><?= num($t['species']) ?></span>
    <span class="tile-note"><?= num($t['species_with_genome']) ?> with a genome or transcriptome</span></div>
  <div class="tile"><span class="tile-label">rRNA records alone</span><span class="tile-value"><?= num($t['species_rrna_only']) ?></span>
    <span class="tile-note">species known by published rRNA reads</span></div>
  <div class="tile"><span class="tile-label">Taxa</span><span class="tile-value"><?= num($t['genera']) ?></span>
    <span class="tile-note">genera in <?= num($t['families']) ?> families, <?= num($t['phyla']) ?> phyla, <?= num($t['supergroups']) ?> supergroups</span></div>
  <div class="tile"><span class="tile-label">Marker categories</span><span class="tile-value"><?= num($t['marker_categories']) ?></span>
    <span class="tile-note"><?= num($t['marker_loci']) ?> loci and genes, <?= num($t['sequences']) ?> sequences</span></div>
</section>

<h2 id="sources">Species by genome source and supergroup</h2>
<hr>
<p>Each species' sequences come from one genome: NCBI's best assembly of it (RefSeq's copy where there is one), else
  an assembly of TOGA2's, an INSDC transcriptome, EukProt's records; species with no genome are known by their
  published rRNA reads alone. Open a supergroup for its groups.</p>
<ul class="legend" aria-label="Sources">
<?php foreach ($sources['columns'] as $s): ?>
  <li title="<?= h($sources['text'][$s] ?? $s) ?>"><span class="swatch s<?= $slot[$s] ?? 0 ?>"></span><?= h($s) ?></li>
<?php endforeach; ?>
</ul>
<div class="bars">
<?php foreach ($sources['supergroups'] as $sg): ?>
  <details class="bar-group">
    <summary class="bar-row">
      <span class="bar-name"><?= h($sg['name']) ?></span>
      <?= stacked_bar($sg['counts'], $scale, $sg['name'], $slot, $sources['text']) ?>
      <span class="bar-value"><?= num($sg['species']) ?></span>
    </summary>
<?php foreach ($sg['groups'] as $g): ?>
    <div class="bar-row bar-sub">
      <span class="bar-name"><?= h($g['name']) ?></span>
      <?= stacked_bar($g['counts'], $scale, $g['name'], $slot, $sources['text']) ?>
      <span class="bar-value"><?= num($g['species']) ?></span>
    </div>
<?php endforeach; ?>
  </details>
<?php endforeach; ?>
</div>
<details class="table-view">
  <summary>Show as a table</summary>
  <div class="table-scroll">
  <table>
    <thead><tr><th scope="col">Supergroup / group</th>
<?php foreach ($sources['columns'] as $s): ?>
      <th scope="col" class="num"><?= h($s) ?></th>
<?php endforeach; ?>
      <th scope="col" class="num">Species</th></tr></thead>
    <tbody>
<?php foreach ($sources['supergroups'] as $sg): ?>
      <tr class="row-major"><th scope="row"><?= h($sg['name']) ?></th>
<?php foreach ($sources['columns'] as $s): ?>
        <td class="num"><?= num($sg['counts'][$s] ?? 0) ?></td>
<?php endforeach; ?>
        <td class="num"><?= num($sg['species']) ?></td></tr>
<?php foreach ($sg['groups'] as $g): ?>
      <tr><th scope="row" class="indent"><?= h($g['name']) ?></th>
<?php foreach ($sources['columns'] as $s): ?>
        <td class="num"><?= num($g['counts'][$s] ?? 0) ?></td>
<?php endforeach; ?>
        <td class="num"><?= num($g['species']) ?></td></tr>
<?php endforeach; ?>
<?php endforeach; ?>
    </tbody>
  </table>
  </div>
</details>

<h2 id="taxonomy">Taxonomy</h2>
<hr>
<p>The species by their lineage: supergroup, group, then NCBI's phylum, class, order, family and genus (an empty rank
  takes the name above it with <code>_X</code>). Open a taxon for its subtaxa: each bar is the share of the opened
  taxon's species in that subtaxon.</p>
<div class="tree-tools">
  <ul class="legend" aria-label="Species">
    <li><span class="swatch tree-genome"></span>with a genome or transcriptome</li>
    <li><span class="swatch tree-rrna"></span>rRNA records alone</li>
  </ul>
  <label class="tree-search">Find a taxon <input type="search" id="tree-search" placeholder="e.g. Saccharomyces" autocomplete="off"></label>
</div>
<div id="tree-results" class="tree-results" hidden></div>
<div id="taxonomy-tree" class="tree" data-src="data/taxonomy.json" aria-live="polite">
  <noscript>
    <ul class="tree-plain">
<?php foreach ($sources['supergroups'] as $sg): ?>
      <li><b><?= h($sg['name']) ?></b> (<?= num($sg['species']) ?> species):
        <?= h(implode(', ', array_map(fn($g) => $g['name'] . ' ' . num($g['species']), $sg['groups']))) ?></li>
<?php endforeach; ?>
    </ul>
  </noscript>
</div>

<h2 id="markers">Marker gene categories</h2>
<hr>
<p>What the release holds of each kind of marker: its loci or genes, its sequences and the species they come from,
  and its sequences' confidence tiers (high: checked against its relatives, nothing against it; unverified: nothing
  against it, but no relative near enough to check it by; medium: something may be wrong, such as a fragment,
  a gap or a contig's end; low: a check calls it another organism's).</p>
<ul class="legend" aria-label="Confidence">
  <li><span class="swatch c-high"></span>high</li>
  <li><span class="swatch c-unverified"></span>unverified</li>
  <li><span class="swatch c-medium"></span>medium</li>
  <li><span class="swatch c-low"></span>low</li>
</ul>
<div class="cards">
<?php foreach ($categories as $c):
    $conf = $c['confidence'] + ['unverified' => 0];  // a format-5 release's has none
    $all = max(1, array_sum($conf));
    $max = 1;
    foreach ($c['items'] as $i) {
        $max = max($max, $i['sequences']);
    }
?>
  <article class="card">
    <h3><?= h($c['title']) ?></h3>
    <p class="card-what"><?= h($c['what']) ?></p>
    <dl class="card-numbers">
      <div><dt><?= $c['key'] === 'busco' || $c['key'] === 'markers' ? 'Loci' : ($c['key'] === 'rrna' ? 'Regions' : 'Genes') ?></dt><dd><?= num($c['loci']) ?></dd></div>
      <div><dt><?= $c['key'] === 'busco' || $c['key'] === 'markers' ? 'Copies' : 'Records' ?></dt><dd><?= num($c['sequences']) ?></dd></div>
      <div><dt>Species</dt><dd><?= num($c['species']) ?></dd></div>
    </dl>
    <div class="stack stack-thin" role="img" aria-label="confidence: high <?= num($conf['high']) ?>, unverified <?= num($conf['unverified']) ?>, medium <?= num($conf['medium']) ?>, low <?= num($conf['low']) ?>">
<?php foreach (['high', 'unverified', 'medium', 'low'] as $tier): if ($conf[$tier] > 0): ?>
      <span class="seg c-<?= $tier ?>" style="width:<?= number_format(100 * $conf[$tier] / $all, 3, '.', '') ?>%" data-tip="<?= h($c['title'] . ' · ' . $tier . ': ' . num($conf[$tier]) . ' (' . pct($conf[$tier] / $all) . ')') ?>"></span>
<?php endif; endforeach; ?>
    </div>
    <p class="card-conf"><?= pct($conf['high'] / $all, 0) ?> high · <?= pct($conf['unverified'] / $all, 0) ?> unverified · <?= pct($conf['medium'] / $all, 0) ?> medium · <?= pct($conf['low'] / $all, 0) ?> low</p>
    <h4><?= h(ucfirst($c['items_title'])) ?></h4>
    <ul class="items">
<?php foreach ($c['items'] as $i): ?>
      <li data-tip="<?= h($i['name'] . ': ' . num($i['sequences']) . ' ' . $c['items_unit'] . (isset($i['species']) && $i['species'] !== null ? ', ' . num($i['species']) . ' species' : '')) ?>">
        <span class="item-name"><?= h($i['name']) ?></span>
        <span class="item-bar"><span style="width:<?= number_format(100 * $i['sequences'] / $max, 3, '.', '') ?>%"></span></span>
        <span class="item-value"><?= num($i['sequences']) ?></span>
      </li>
<?php endforeach; ?>
    </ul>
  </article>
<?php endforeach; ?>
</div>

<?php $b = $stats['benchmarks'] ?? null; if ($b): ?>
<h2 id="benchmarks">Benchmarks</h2>
<hr>
<p>Each build measures its own package (<a href="<?= h(DOCS_URL) ?>/benchmarks.md">what each benchmark measures</a>):</p>
<section class="tiles">
<?php if (isset($b['completeness'])): ?>
  <div class="tile"><span class="tile-label">Median complete BUSCOs per genome</span><span class="tile-value"><?= pct($b['completeness']['median']) ?></span>
    <span class="tile-note"><?= pct($b['completeness']['at_least_90']) ?> of genomes at 90% or more</span></div>
<?php endif; ?>
<?php if (isset($b['resolution']['BUSCO (all)'])): ?>
  <div class="tile"><span class="tile-label">BUSCO copies under 99% identical to any other species</span><span class="tile-value"><?= pct($b['resolution']['BUSCO (all)']['under_99']) ?></span>
    <span class="tile-note">of <?= num($b['resolution']['BUSCO (all)']['copies']) ?> copies<?php if (isset($b['resolution']['18S'])): ?>; 18S: <?= pct($b['resolution']['18S']['under_99']) ?><?php endif; ?></span></div>
<?php endif; ?>
<?php if (isset($b['rrna']['lsu_complete_of_operons'])): ?>
  <div class="tile"><span class="tile-label">rDNA operons with a whole 28S</span><span class="tile-value"><?= pct($b['rrna']['lsu_complete_of_operons']) ?></span>
    <span class="tile-note">of <?= num($b['rrna']['operons']) ?> operons</span></div>
<?php endif; ?>
</section>
<?php if (!empty($b['completeness']['by_supergroup'])): ?>
<h3>Median BUSCO completeness by group</h3>
<ul class="items items-wide">
<?php foreach ($b['completeness']['by_supergroup'] as $name => $v): ?>
  <li data-tip="<?= h($name . ': ' . pct($v)) ?>"><span class="item-name"><?= h($name) ?></span>
    <span class="item-bar"><span style="width:<?= number_format(100 * $v, 3, '.', '') ?>%"></span></span>
    <span class="item-value"><?= pct($v) ?></span></li>
<?php endforeach; ?>
</ul>
<?php endif; ?>
<?php endif; ?>
<div id="tip" class="tip" role="tooltip" hidden></div>
