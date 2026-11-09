<?php
// A short tutorial: install, run, use. The GitHub README and docs/ hold the rest.
$listing = release_json('downloads.json');
$pkg = h($listing['package'] ?? 'Eukore_2026.11');
?>
<p class="page-intro">This page is the short way through Eukore: getting a release, using it, and building one
  yourself. The <a href="<?= h(README_URL) ?>">README on GitHub</a> and the documentation in
  <a href="<?= h(DOCS_URL) ?>">docs/</a> tell the whole of it: <a href="<?= h(DOCS_URL) ?>/building.md">building.md</a>
  (installation, options, resources), <a href="<?= h(DOCS_URL) ?>/FORMAT.md">FORMAT.md</a> (every file and column),
  <a href="<?= h(DOCS_URL) ?>/methods.md">methods.md</a> and <a href="<?= h(DOCS_URL) ?>/loci.md">loci.md</a>.</p>

<nav class="toc" aria-label="On this page">
  <a href="#get">1. Get a release</a>
  <a href="#use">2. Use it</a>
  <a href="#install">3. Install the build</a>
  <a href="#run">4. Build Eukore</a>
  <a href="#publish">5. Publish a release</a>
</nav>

<h2 id="get">1. Get a release</h2>
<hr>
<p>Take the complete package, or only the part you need, from the <a href="main.php?site=downloads">downloads</a>
  (the marker genes with or without introns, the rRNA and organelle markers, the sketches, the metadata tables):</p>
<pre><code>md5sum -c MD5SUMS
tar -xzf <?= $pkg ?>.tar.gz          # or: tar -I "bgzip -@ 8" -xf <?= $pkg ?>.tar.gz</code></pre>
<p>A package is a folder of plain TSV tables and gzipped FASTA files:</p>
<pre><code><?= $pkg ?>/
├── manifest.json     versions, inputs, counts, every file's size and MD5 (written last)
├── species.tsv       one row per species: supergroup, group, NCBI and six references' lineages
├── genomes.tsv       the species' genomes: accession, quality flags, rRNA copies, genetic code
├── busco/            loci.tsv, copies.tsv, introns.tsv; copies/, cds/, proteins/&lt;locus&gt;.*.gz
├── markers/          Eukore's own loci (RPB1, RPB2, EF-2, TEF1, ...), as busco/
├── rrna/             records.tsv, records.fna.gz, regions.tsv (18S V1-V9), amplicons.tsv
├── organelle/        genes.tsv, genes.fna.gz, cds.fna.gz, proteins.faa.gz, molecules.fna.gz
├── metadata/         assembly.tsv.gz, ncbi.tsv.gz, codons.tsv.gz, sketches.zip
├── taxonomy/         NCBI's taxonomy, cut to the package
└── benchmarks/       the build's measures of the package</code></pre>

<h2 id="use">2. Use it</h2>
<hr>
<p>Every table is read by its header by any tool (R, pandas, awk), every FASTA by any reader of gzip. From the shell:</p>
<pre><code># the BUSCOs called in at least 90% of the species, and their names
awk -F'\t' 'NR == 1 || $9 >= 0.9 {print $1, $2, $9}' <?= $pkg ?>/busco/loci.tsv | head
# a locus's high-confidence copies (introns and all), by their accessions in copies.tsv
awk -F'\t' '$2 == "10248812at2759" && $29 == "high" {print $1}' <?= $pkg ?>/busco/copies.tsv > high.txt
seqkit grep -f high.txt <?= $pkg ?>/busco/copies/10248812at2759.fna.gz > high.fna
# the V4 region's coordinates on each 18S
awk -F'\t' '$2 == "V4"' <?= $pkg ?>/rrna/regions.tsv | head</code></pre>
<p>With Eukore's own reader (Python, the repository's <code>eukore/</code>):</p>
<pre><code>from eukore import package, ncbi
pkg = "<?= $pkg ?>"
print(package.check(pkg, deep=True) or "whole")           # every file against the manifest's size and MD5
species = package.read_table(pkg, "species.tsv")           # taxid, supergroup, group, lineages in NCBI, PR2, ...
loci = package.read_table(pkg, "busco/loci.tsv")           # each BUSCO's call rates, to choose loci by
copies = package.read_table(pkg, "busco/copies.tsv")       # each copy's call, checks and confidence
seqs = package.read_copies(pkg, "10248812at2759")          # {accession: the copy, introns and all}
cds = package.read_copies(pkg, "10248812at2759", "cds")    # {accession: its coding sequence}; "proteins" too
tax = ncbi.Ncbi(f"{pkg}/taxonomy")                         # the NCBI taxonomy the package was built on</code></pre>
<p>A typical use keeps the <code>high</code> and <code>medium</code> copies of the BUSCOs it chooses (by
  <code>busco/loci.tsv</code>'s call rates in the groups it cares about, <code>call_&lt;group&gt;</code>;
  <a href="<?= h(DOCS_URL) ?>/loci.md">loci.md</a>) and masks, cuts or renames them to its own needs. protal builds
  its eukaryote database this way (<code>download_ref_data.py --eukdb <?= $pkg ?></code>).</p>

<h2 id="install">3. Install the build</h2>
<hr>
<p>Eukore runs on Linux (or WSL) with a conda-style package manager; its tools (miniprot, HMMER, barrnap, ITSx,
  vsearch, Infernal) are bioconda's.</p>
<pre><code>git clone <?= h(GITHUB_URL) ?>.git && cd Eukore
micromamba create -f envs/eukore.yaml && micromamba activate eukore
python3 -m unittest discover -s tests -t tests            # the tests, with stand-in tools: seconds</code></pre>
<p>Optional: tiara (<code>envs/tiara.yaml</code>) classifies the windows around each copy, NCBI's FCS-GX with its
  ~470 GB database screens them for contamination (<a href="<?= h(DOCS_URL) ?>/building.md#install">building.md</a>).</p>

<h2 id="run">4. Build Eukore</h2>
<hr>
<p>A <b>trial</b> builds a few clades on one machine, here the genus <i>Saccharomyces</i> (NCBI taxid 4930); its
  manifest marks it partial:</p>
<pre><code>python3 eukore_build.py --work work --out out --groups fungi --species 4930 -t 8</code></pre>
<p>The <b>full reference</b> is built in three stages, the first where the network is, the rest on compute nodes;
  every stage is resumable, a rerun does only what is not done:</p>
<pre><code>python3 eukore_build.py --work WORK --out OUT --stage fetch                 # sources and genomes
python3 eukore_build.py --work WORK --out OUT --stage process -t 64 \
    [--memory_gb 200] [--scratch /local/disk]                               # each genome's genes and rRNA
python3 eukore_build.py --work WORK --out OUT --stage package -t 64 \
    [--fcs_gx_db /dev/shm/gxdb/all] [--tiara .../bin/tiara]                 # checks, package, benchmarks, release</code></pre>
<p>The package stage writes <code>OUT/Eukore_&lt;version&gt;/</code>, then <code>OUT/release/</code>: this
  website's statistics and every download, BGZF-compressed on <code>-t</code> threads (<code>--no_release</code>:
  none). <code>--stage release</code> writes <code>OUT/release/</code> again from a finished package.</p>

<h2 id="publish">5. Publish a release</h2>
<hr>
<p><code>OUT/release/</code> is laid out as this website is: copying it onto the site publishes the statistics and
  the downloads.</p>
<pre><code>OUT/release/
├── data/stats.json, taxonomy.json, downloads.json    read by the Statistics and Downloads pages
└── downloads/<?= $pkg ?>/                           the files, MD5SUMS

rsync -a OUT/release/ user@server:/var/www/eukore/</code></pre>
<p>The site itself (these PHP pages) is the repository's <code>website/</code> folder.</p>
