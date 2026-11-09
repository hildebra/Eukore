<?php
// The release's downloads (data/downloads.json: release.py's), grouped as the release groups them.
$listing = release_json('downloads.json');
$stats = release_json('stats.json');
?>
<p class="page-intro">Eukore's releases are plain TSV tables and FASTA files, every one compressed as
  <b>BGZF</b> (blocked gzip) at the highest level: <code>gzip</code>, <code>zcat</code> and <code>tar -xzf</code>
  read them as any <code>.gz</code>, and <code>bgzip -@</code> decompresses them on many threads. The format of each
  file is described in <a href="<?= h(DOCS_URL) ?>/FORMAT.md">docs/FORMAT.md</a>.</p>

<?php if (!$listing || empty($listing['files'])): ?>
<?= no_release('the downloads') ?>
<?php else:
    $base = (DOWNLOAD_URL !== '' ? rtrim(DOWNLOAD_URL, '/') . '/' . basename($listing['folder']) . '/' : $listing['folder']);
    $groups = [];
    foreach ($listing['files'] as $f) {
        $groups[$f['group']][] = $f;
    }
?>
<p>Release <?= release_badge($stats ?: ['package' => $listing['package']]) ?>, built
  <?= h(date_text($listing['built'] ?? null)) ?>.
  Checksums of every file: <a href="<?= h($base) ?>MD5SUMS">MD5SUMS</a>.</p>
<?php foreach ($groups as $group => $files): ?>
<h2><?= h($group) ?></h2>
<hr>
<div class="downloads">
<?php foreach ($files as $f): ?>
  <div class="download">
    <div class="download-main">
      <a class="download-link" href="<?= h($base . $f['file']) ?>" download><?= h($f['title']) ?></a>
      <p><?= h(ucfirst($f['what'])) ?>.</p>
      <p class="download-file"><code><?= h($f['file']) ?></code></p>
    </div>
    <dl class="download-meta">
      <div><dt>Size</dt><dd><?= h(bytes_text($f['bytes'])) ?></dd></div>
<?php if (!empty($f['files'])): ?>
      <div><dt><?= $f['key'] === 'sketches' ? 'Signatures' : 'Files' ?></dt><dd><?= num($f['files']) ?></dd></div>
<?php endif; ?>
      <div class="md5"><dt>MD5</dt><dd><code><?= h($f['md5']) ?></code></dd></div>
    </dl>
  </div>
<?php endforeach; ?>
</div>
<?php endforeach; ?>
<?php endif; ?>

<h2>Unpacking</h2>
<hr>
<p>A bundle unpacks into the package's own folder, <code><?= h($listing['package'] ?? 'Eukore_<version>') ?>/</code>, its files where
  the complete package has them: bundles unpacked side by side make up those parts of the package.</p>
<?php $p = h($listing['package'] ?? 'Eukore_2026.11'); ?>
<pre><code># each file as published
md5sum -c MD5SUMS
# any tar unpacks a bundle; htslib's bgzip does it on 8 threads
tar -xzf <?= $p ?>_marker_genes.tar.gz
tar -I "bgzip -@ 8" -xf <?= $p ?>.tar.gz
# a table, or a FASTA file inside a bundle, read as any .gz
bgzip -d -@ 8 -c <?= $p ?>_assembly.tsv.gz | head
zcat <?= $p ?>/busco/copies/10248812at2759.fna.gz | head -2
# the genomes' sketches, as one sourmash collection
sourmash search query.sig <?= $p ?>_sketches.sig.gz</code></pre>

<h2>The software</h2>
<hr>
<p>Eukore's build, which makes these files from NCBI's assemblies, is on GitHub:
  <a href="<?= h(GITHUB_URL) ?>"><?= h(GITHUB_URL) ?></a> (GPL-2.0). Building a package yourself:
  <a href="main.php?site=tutorial">the tutorial</a>.</p>
<h3>License</h3>
<p>The data a release holds keep their sources' terms (<a href="main.php?site=acknowledgements#sources">the
  sources</a>): check them before redistributing a package.</p>
