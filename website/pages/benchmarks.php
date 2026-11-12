<?php
// The benchmarks of a release (data/benchmarks.json and data/figures/: benchmarks/eukore_bench/site.py writes them
// from the benchmarks' results, after the build's release): what each database beside Eukore is and holds, the
// breadth at every rank, the markers one species links, the named lists' organisms, the checks, the same species'
// markers in two databases, the identity that parts taxa, the tree, the changes since the release before.
$b = release_json('benchmarks.json');
if (!$b) {
    echo '<div class="notice">No benchmarks have been published for this release yet: they appear once '
        . '<code>python3 -m eukore_bench.site RESULTS --release OUT/release</code> has written them beside the '
        . 'release\'s data (<a href="' . h(GITHUB_URL) . '/tree/main/benchmarks">benchmarks/</a>).</div>';
    return;
}
$t = $b['tables'] ?? [];
$s = $b['summary'] ?? [];
$figures = [];
foreach ($b['figures'] ?? [] as $f) {
    $figures[$f['name']] = $f;
}

// A table cell: a count with its thousands, anything else as it is.
function cell($v): string
{
    return is_numeric($v) && !preg_match('/^0\d/', (string) $v) && strpos((string) $v, '.') === false
        ? num($v) : h($v);
}

// A figure of the release (data/figures/<name>.png) with its caption; nothing where it was not drawn.
function figure_of(array $figures, string $name): string
{
    if (!isset($figures[$name])) {
        return '';
    }
    $f = $figures[$name];
    return '<figure class="figure"><img src="data/' . h($f['file']) . '" alt="' . h($f['caption']) . '" loading="lazy">'
        . '<figcaption>' . h($f['caption']) . '</figcaption></figure>';
}

// Rows (arrays by column) as a table: $columns [key => heading]; the numeric columns right-aligned.
function rows_table(array $rows, array $columns, array $numeric = []): string
{
    $out = '<div class="table-scroll"><table><thead><tr>';
    foreach ($columns as $key => $label) {
        $out .= '<th scope="col"' . (in_array($key, $numeric, true) ? ' class="num"' : '') . '>' . h($label) . '</th>';
    }
    $out .= '</tr></thead><tbody>';
    foreach ($rows as $r) {
        $out .= '<tr>';
        foreach ($columns as $key => $label) {
            $out .= '<td' . (in_array($key, $numeric, true) ? ' class="num"' : '') . '>' . cell($r[$key] ?? '–') . '</td>';
        }
        $out .= '</tr>';
    }
    return $out . '</tbody></table></div>';
}
?>
<p class="page-intro">What <?= release_badge(['package' => $b['package'] ?? 'Eukore'] + (release_json('stats.json') ?? [])) ?>
  holds and how good its records are, measured on the finished release
  (<a href="<?= h(GITHUB_URL) ?>/tree/main/benchmarks">benchmarks/</a>, <?= h(date_text($b['made'] ?? null)) ?>).
  Eukore is a database: the other databases are set beside it as collections of records.</p>

<?php if (!empty($t['databases'])): ?>
<h2 id="databases">Eukore beside other databases</h2>
<hr>
<p>What each database is, and what was counted of the version read: its species as NCBI names them, those Eukore
  holds too, and its records.</p>
<?= rows_table($t['databases'], ['database' => 'Database', 'markers' => 'Markers', 'taxa' => 'Taxa',
                                 'sequences' => 'A record is', 'curation' => 'Checked', 'reference' => 'Reference']) ?>
<?= rows_table($t['databases'], ['database' => 'Database', 'version' => 'Version', 'species' => 'NCBI species',
                                 'both' => 'In Eukore too', 'eukore_only' => "Eukore's alone",
                                 'database_only' => 'Its alone', 'records' => 'Records', 'bases' => 'Bases'],
               ['species', 'both', 'eukore_only', 'database_only', 'records', 'bases']) ?>
<?php if (isset($s['upset']['eukore_alone'])): ?>
<p><b><?= num($s['upset']['eukore_alone']) ?></b> species are in Eukore and in none of the other
  <?= num(count($s['upset']['databases'] ?? []) - 1) ?> databases.</p>
<?php endif; ?>
<?= figure_of($figures, 'upset') ?>
<?php endif; ?>

<?php if (!empty($t['breadth_ranks'])): ?>
<h2 id="breadth">Breadth at every rank</h2>
<hr>
<p>The eukaryote taxa NCBI names, those it holds a genome of, and Eukore's: with a genome, and with any record
  (rRNA reads alone too).</p>
<?= rows_table($t['breadth_ranks'], ['rank' => 'Rank', 'ncbi_named' => 'NCBI: named', 'ncbi_genome' => 'NCBI: with a genome',
                                     'eukore_genome' => 'Eukore: with a genome', 'eukore' => 'Eukore: any record'],
               ['ncbi_named', 'ncbi_genome', 'eukore_genome', 'eukore']) ?>
<?= figure_of($figures, 'breadth_ranks') ?>
<?php endif; ?>

<?php if (!empty($t['content'])): ?>
<h2 id="content">The markers one species links</h2>
<hr>
<p>The species holding each kind of marker, and those holding it from their genome's own assembly: genes and
  barcodes of one organism, linked through one assembly.</p>
<?= rows_table($t['content'], ['marker' => 'Marker', 'species' => 'Species, any source',
                               'same_assembly' => 'From their own assembly', 'of' => 'Of species'],
               ['species', 'same_assembly', 'of']) ?>
<?= figure_of($figures, 'content') ?>
<?php endif; ?>

<?php if (!empty($t['panels'])):
    $panels = [];
    foreach ($t['panels'] as $r) {
        $panels[$r['panel']][] = $r;
    }
?>
<h2 id="lists">Is my organism in it?</h2>
<hr>
<p>The <a href="main.php?site=stats#taxonomy">taxonomy tree</a> finds any taxon. The organisms of a few named lists,
  by what the release holds of them (open a list):</p>
<?php foreach ($panels as $name => $rows):
    $with = count(array_filter($rows, fn($r) => $r['status'] === 'genome'));
?>
<details class="table-view">
  <summary><b><?= h($name) ?></b>: <?= num($with) ?> of <?= num(count($rows)) ?> with a genome</summary>
  <?= rows_table($rows, ['entry' => 'Organism', 'ncbi_name' => "NCBI's name", 'rank' => 'Rank', 'status' => 'In Eukore',
                         'eukore_genome' => 'Species with a genome', 'loci' => 'BUSCO genes',
                         'rrna_records' => 'rRNA records', 'organelle_genes' => 'Organelle genes', 'note' => 'Note'],
                 ['eukore_genome', 'loci', 'rrna_records', 'organelle_genes']) ?>
</details>
<?php endforeach; ?>
<?php endif; ?>

<?php if (!empty($s['quality']['genome']) || !empty($t['agreement']) || !empty($t['gencode'])): ?>
<h2 id="quality">How good the records are</h2>
<hr>
<?php if (!empty($s['quality']['genome']['of'])):
    $q = $s['quality']['genome'];
?>
<p>Of the <?= num($q['of']) ?> species' genomes, those each check flags (a flagged genome is kept and marked, so
  that a user chooses):</p>
<ul class="items items-wide">
<?php foreach ($q['checks'] as $check => $n): ?>
  <li data-tip="<?= h($check . ': ' . num($n) . ' of ' . num($q['of'])) ?>"><span class="item-name"><?= h($check) ?></span>
    <span class="item-bar"><span style="width:<?= number_format(100 * $n / max(1, $q['of']), 3, '.', '') ?>%"></span></span>
    <span class="item-value"><?= pct($n / max(1, $q['of'])) ?></span></li>
<?php endforeach; ?>
</ul>
<?php endif; ?>
<?= figure_of($figures, 'quality') ?>
<?php if (!empty($t['gencode'])): ?>
<h3>Genetic codes</h3>
<p>Each genome's genetic code is checked against its own BUSCO genes:</p>
<?= rows_table($t['gencode'], ['genomes' => 'Genomes', 'agrees' => "Agree with NCBI's", 'differs' => 'Differ',
                               'undecided' => 'Undecided', 'unchecked' => 'Not checked',
                               'best_tables' => 'Tables that fit those that differ'],
               ['genomes', 'agrees', 'differs', 'undecided', 'unchecked']) ?>
<?php endif; ?>
<?php if (!empty($t['agreement'])): ?>
<h3>The same species in two databases</h3>
<p>Where Eukore and another database hold the same species' marker: the records that agree with Eukore's genomes,
  those that are another species' gene, and those unlike the species' own. A disagreement does not say which side is
  wrong.</p>
<?= rows_table($t['agreement'], ['database' => 'Database', 'marker' => 'Marker', 'records' => 'Records',
                                 'agree' => 'Agree', 'disagree' => "Another species'", 'unlike' => 'Unlike',
                                 'undecided' => 'Undecided', 'species_differ' => 'Species that differ'],
               ['records', 'agree', 'disagree', 'unlike', 'undecided', 'species_differ']) ?>
<?= figure_of($figures, 'agreement') ?>
<?php endif; ?>
<?php endif; ?>

<?php if (!empty($s['tree'])): ?>
<h2 id="tree">The genes' tree and NCBI's taxonomy</h2>
<hr>
<p>A tree of <?= num($s['tree']['genomes']) ?> genomes' BUSCO proteins (<?= num($s['tree']['loci']) ?> genes).
  NCBI's taxa with two genomes or more in it that are one branch of the tree:</p>
<ul class="items items-wide">
<?php foreach ($s['tree']['ranks'] as $rank => $r): ?>
  <li data-tip="<?= h($rank . ': ' . num($r['monophyletic']) . ' of ' . num($r['taxa'])) ?>"><span class="item-name"><?= h($rank) ?> (<?= num($r['taxa']) ?>)</span>
    <span class="item-bar"><span style="width:<?= number_format(100 * $r['share'], 3, '.', '') ?>%"></span></span>
    <span class="item-value"><?= pct($r['share']) ?></span></li>
<?php endforeach; ?>
</ul>
<p><?= num($s['tree']['strays']) ?> genomes sit apart from their name's relatives (the benchmarks'
  <code>tree_strays.tsv</code> lists them).</p>
<?= figure_of($figures, 'tree') ?>
<?php endif; ?>

<?php if (!empty($t['thresholds'])): ?>
<h2 id="thresholds">Which identity parts taxa</h2>
<hr>
<p>For a user choosing a marker and a cutoff: the identity at or above which two sequences are of one taxon of a rank,
  with the shares it keeps together and apart. Published thresholds of other markers and taxa are given for context.</p>
<?= rows_table($t['thresholds'], ['marker' => 'Marker', 'rank' => 'Rank', 'threshold' => 'Identity (%)',
                                  'within_kept' => 'Kept together', 'between_apart' => 'Kept apart',
                                  'source' => 'Source'], ['threshold', 'within_kept', 'between_apart']) ?>
<?= figure_of($figures, 'thresholds') ?>
<?php endif; ?>

<?php if (!empty($t['novel'])): ?>
<h2 id="novel">Species a genome is the first sequence for</h2>
<hr>
<p>The species whose marker Eukore takes from their genome, and those of them that none of the marker's references
  names (by name: an upper bound).</p>
<?= rows_table($t['novel'], ['marker' => 'Marker', 'references' => 'References asked', 'genome' => 'From a genome',
                             'named' => 'Named in a reference', 'new' => 'In none'], ['genome', 'named', 'new']) ?>
<?= figure_of($figures, 'novel') ?>
<?php endif; ?>

<?php if (!empty($s['changes']) || !empty($s['growth'])): ?>
<h2 id="changes">This release among releases</h2>
<hr>
<?php if (!empty($s['changes'])):
    $c = $s['changes'];
?>
<p>Since <?= h($c['previous_version'] ?? 'the release before') ?>: <?= num($c['added'] ?? 0) ?> species added,
  <?= num($c['removed'] ?? 0) ?> removed, <?= num($c['gained_genome'] ?? 0) ?> gained a genome;
  <?= num($c['replaced'] ?? 0) ?> genomes replaced by another assembly, <?= num($c['updated'] ?? 0) ?> by a newer
  version of the same (<?= num($c['previous'] ?? 0) ?> to <?= num($c['species'] ?? 0) ?> species).</p>
<?= figure_of($figures, 'changes') ?>
<?php endif; ?>
<?php if (!empty($s['growth'])): ?>
<p><?= pct($s['growth']['released_in_last_two_years']) ?> of the release's <?= num($s['growth']['genomes']) ?> genomes
  were released by NCBI in <?= h($s['growth']['last_year'] - 1) ?> or <?= h($s['growth']['last_year']) ?>.</p>
<?= figure_of($figures, 'growth') ?>
<?php endif; ?>
<?php endif; ?>
<div id="tip" class="tip" role="tooltip" hidden></div>
