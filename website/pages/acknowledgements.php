<?php
// Who funds and makes Eukore, and the data and tools it is built from.
?>
<div class="split">
  <h2 class="split-title">Funding</h2>
  <div class="split-body">
    <p class="lead">Eukore is funded by the <a href="https://www.ukri.org/councils/bbsrc/">Biotechnology and Biological
      Sciences Research Council (BBSRC)</a>, part of UK Research and Innovation.</p>
<?php if (BBSRC_GRANTS): ?>
    <ul>
<?php foreach (BBSRC_GRANTS as [$ref, $grant_title]): ?>
      <li><b><?= h($ref) ?></b><?= $grant_title ? ': ' . h($grant_title) : '' ?></li>
<?php endforeach; ?>
    </ul>
<?php endif; ?>
  </div>
</div>

<div class="split">
  <h2 class="split-title">The lab</h2>
  <div class="split-body">
    <p class="lead">Eukore is developed by the Hildebrand lab: <a href="<?= h(LAB_URL) ?>"><?= h(LAB_URL) ?></a>.</p>
    <p>Questions, ideas and bug reports are welcome as <a href="<?= h(ISSUES_URL) ?>">GitHub issues</a>.</p>
  </div>
</div>

<div class="split" id="sources">
  <h2 class="split-title">Data</h2>
  <div class="split-body">
    <p>Eukore is built from public data, each under its own terms; please check them before redistributing a
      package, and cite them with Eukore:</p>
    <ul>
      <li>NCBI's GenBank and RefSeq assemblies, taxonomy, FCS-GX and ANI reports, RefSeq's organelle genomes;
        INSDC's transcriptome shotgun assemblies (via ENA)</li>
      <li>BUSCO's eukaryota_odb12 (OrthoDB 12.1); KOfam's profiles, Swiss-Prot's and NCBI's proteins (Eukore's own
        marker loci)</li>
      <li>EUKARYOME 2.1; the rRNA Operon Database (ROD) 1.2; Rfam's eukaryote SSU model</li>
      <li>PR2 5.1.1, KSGP, SILVA, UNITE, MIDORI2; Open Tree of Life, GBIF, WoRMS, IRMNG, BOLD</li>
      <li>TOGA2's vertebrate assemblies and gene models; Ensembl's and VEuPathDB's gene models</li>
      <li>the records of species without an assembly from EukDetect2 (CC BY 4.0; Shih, Zhao, Pollard &amp; Lind
        2026), of EukProt's genomes and transcriptomes</li>
    </ul>
  </div>
</div>

<div class="split">
  <h2 class="split-title">Tools</h2>
  <div class="split-body">
    <p>Eukore runs miniprot, HMMER, barrnap, ITSx, vsearch, Infernal, sourmash, and optionally tiara and NCBI's
      FCS-GX: please cite them with Eukore.</p>
    <p>Eukore began inside protal's repository and keeps that history.
      Its code is licensed under the GNU General Public License, version 2.0.</p>
  </div>
</div>
