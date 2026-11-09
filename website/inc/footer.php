<?php
// The site's footer (as LotuS2's: who makes it, the pages), then the scripts.
?>
</main>
<footer class="site-footer">
  <div class="container footer-row">
    <div class="footer-about">
      <h3>About</h3>
      <p>Eukore is developed by the <a href="<?= h(LAB_URL) ?>">Hildebrand lab</a> and funded by the
        Biotechnology and Biological Sciences Research Council (BBSRC). Its code is free software (GPL-2.0) on
        <a href="<?= h(GITHUB_URL) ?>">GitHub</a>; questions and bug reports are welcome as
        <a href="<?= h(ISSUES_URL) ?>">issues</a>.</p>
    </div>
    <div class="footer-nav">
      <h3>Navigation</h3>
      <ul>
        <li><a href="index.php">Home</a></li>
<?php foreach (PAGES as $key => $label): ?>
        <li><a href="main.php?site=<?= h($key) ?>"><?= h($label) ?></a></li>
<?php endforeach; ?>
      </ul>
    </div>
  </div>
  <div class="container footer-base">
    <p>© <?= date('Y') ?> The Hildebrand lab · Eukore's data keep their sources' terms
      (<a href="main.php?site=acknowledgements#sources">sources</a>)</p>
  </div>
</footer>
<script src="js/site.js" defer></script>
</body>
</html>
