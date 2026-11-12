<?php
// Eukore's website: its settings. The pages read the release's data from DATA_DIR (data/stats.json, taxonomy.json,
// downloads.json: eukore_build.py's OUT/release/data/) and link its files under DOWNLOAD_URL (OUT/release/downloads/).
// Publishing a release: rsync -a OUT/release/ <this folder>/  (README.md).

const SITE_NAME = 'Eukore';
const SITE_TAGLINE = 'Eukaryote marker genes from every genome NCBI holds';
const GITHUB_URL = 'https://github.com/hildebra/Eukore';
const README_URL = 'https://github.com/hildebra/Eukore#readme';
const DOCS_URL = 'https://github.com/hildebra/Eukore/tree/main/docs';
const ISSUES_URL = 'https://github.com/hildebra/Eukore/issues';
const LAB_URL = 'https://hildebra.github.io/';
const DATA_DIR = __DIR__ . '/data';
const DOWNLOAD_URL = '';  // the downloads' base URL: '' for this site's own folder (downloads.json's "folder")

// The BBSRC grants that fund Eukore, shown on the acknowledgements page: add each award's reference and title, e.g.
// ['BB/X000000/1', 'Its title']. None listed: the page names BBSRC alone.
const BBSRC_GRANTS = [];

// The site's pages: main.php?site=<key> (index.php is Home).
const PAGES = [
    'stats' => 'Statistics',
    'benchmarks' => 'Benchmarks',
    'downloads' => 'Downloads',
    'tutorial' => 'Tutorial',
    'acknowledgements' => 'Acknowledgements',
];
