<?php
// The site's pages but Home: main.php?site=<key> (config.php's PAGES), as LotuS2's site.
require __DIR__ . '/config.php';
require __DIR__ . '/inc/data.php';
$page = $_GET['site'] ?? '';
if (!array_key_exists($page, PAGES)) {
    http_response_code(404);
    $page = 'missing';
    $title = 'Not found';
} else {
    $title = PAGES[$page];
}
require __DIR__ . '/inc/header.php';
?>
<div class="page-banner"><div class="container"><h1><?= h($title) ?></h1></div></div>
<div class="container content">
<?php
if ($page === 'missing') {
    echo '<p>There is no such page. <a href="index.php">Home</a></p>';
} else {
    require __DIR__ . '/pages/' . $page . '.php';
}
?>
</div>
<?php require __DIR__ . '/inc/footer.php'; ?>
