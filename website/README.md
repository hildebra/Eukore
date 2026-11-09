# Eukore's website

A small PHP site in the style of [LotuS2's](https://lotus2.earlham.ac.uk/): the logo at the top left, the menu on the
right, a cover on the home page, a dark footer. It needs PHP 7.4 or later and nothing else: no database, no build
step. Its font (Mukta) comes from Google Fonts: host it with the site to send visitors to no third party.

| page | file | |
|---|---|---|
| Home | `index.php` | what Eukore is, why, its workflow (`images/overview_slide.png`, from `docs/figures/`), the release's totals |
| Statistics | `main.php?site=stats` (`pages/stats.php`) | the release's species by genome source and supergroup (open a supergroup for its groups), its taxonomy as a tree (open a taxon for its subtaxa, each with a bar of its share of the taxon's species; a search), its marker categories, its headline benchmarks |
| Downloads | `main.php?site=downloads` (`pages/downloads.php`) | every file of the release with its size and MD5, grouped |
| Tutorial | `main.php?site=tutorial` (`pages/tutorial.php`) | getting, using and building a release, pointing to the GitHub README |
| Acknowledgements | `main.php?site=acknowledgements` (`pages/acknowledgements.php`) | the BBSRC's funding, the lab (https://hildebra.github.io/), the data and tools |

`config.php` holds the settings: the GitHub links, the lab's page, the BBSRC grants to name (`BBSRC_GRANTS`, empty:
BBSRC alone), where the data and the downloads are.

## Its data: a build's release

The pages read nothing but three JSON files in `data/`, which the build writes with the package (`eukore_build.py`'s
package stage, or `--stage release` of a finished package; eukore/release.py) into `OUT/release/`, laid out as this
site is:

```
OUT/release/
├── data/
│   ├── stats.json        the Statistics page's numbers (and Home's totals)
│   ├── taxonomy.json     the species as a tree: [name, rank, species, with a genome, children | accession]
│   └── downloads.json    the Downloads page's files: title, what, file, bytes, MD5
└── downloads/Eukore_<version>/
    ├── Eukore_<version>.tar.gz, _marker_genes.tar.gz, _marker_cds.tar.gz, _rrna_organelle.tar.gz
    ├── Eukore_<version>_sketches.sig.gz, _species.tsv.gz, _genomes.tsv.gz, _assembly.tsv.gz, ...
    └── MD5SUMS
```

Publishing a release is copying that folder onto the site:

```bash
rsync -a website/ user@server:/var/www/eukore/              # the pages (once, and after a change)
rsync -a OUT/release/ user@server:/var/www/eukore/          # a release: its data and its downloads
```

The downloads can live elsewhere (another disk, a file server): copy `OUT/release/downloads/` there and set
`DOWNLOAD_URL` in `config.php` to its URL. `data/` holds the JSON of a trial package (7 clades, 43 species) until a
release replaces it; `downloads/` holds nothing in the repository.

Every download is BGZF at gzip's level 9: any gzip reader takes it, `bgzip -d -@ N` decompresses it on N threads. The
server should send the `.gz` files as they are (no `Content-Encoding: gzip` of its own): Apache's and nginx's
defaults do.

## Trying it locally

```bash
php -S 127.0.0.1:8099 -t website        # then http://127.0.0.1:8099/
```
