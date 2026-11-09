<?php
// The page's head and the site's header: the logo on the left, the menu on the right (as LotuS2's site). $page: the
// page's key (PAGES, or 'home'); $title: its title.
$page = $page ?? 'home';
$title = $title ?? SITE_NAME;
$nav = ['home' => ['index.php', 'Home']];
foreach (PAGES as $key => $label) {
    $nav[$key] = ['main.php?site=' . $key, $label];
}
?><!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title><?= h($page === 'home' ? SITE_NAME . ' — ' . SITE_TAGLINE : $title . ' — ' . SITE_NAME) ?></title>
  <meta name="description" content="Eukore: a checked, complete reference of eukaryote marker genes (BUSCO genes with their introns, rDNA operons, organelle genes) built from every eukaryote genome NCBI holds.">
  <link rel="icon" type="image/svg+xml" href="images/eukore_icon.svg">
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Mukta:wght@200;300;400;500;700&display=swap">
  <link rel="stylesheet" href="css/style.css">
</head>
<body class="page-<?= h($page) ?>">
<header class="site-header">
  <div class="container header-row">
    <a class="site-logo" href="index.php"><img src="images/eukore_header.svg" alt="Eukore: eukaryote marker genes" width="252" height="52"></a>
    <button class="menu-toggle" type="button" aria-controls="site-menu" aria-expanded="false">
      <span class="menu-bars" aria-hidden="true"></span><span class="sr-only">Menu</span>
    </button>
    <nav class="site-nav" id="site-menu" aria-label="Site">
      <ul>
<?php foreach ($nav as $key => [$href, $label]): ?>
        <li<?= $key === $page ? ' class="active"' : '' ?>><a href="<?= h($href) ?>"<?= $key === $page ? ' aria-current="page"' : '' ?>><?= h($label) ?></a></li>
<?php endforeach; ?>
        <li class="nav-github"><a href="<?= h(GITHUB_URL) ?>">GitHub <span aria-hidden="true">↗</span></a></li>
      </ul>
    </nav>
  </div>
</header>
<main>
