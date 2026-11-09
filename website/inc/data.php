<?php
// The release's data (config.php's DATA_DIR) and the pages' small helpers.

function release_json(string $name): ?array
{
    static $cache = [];
    if (!array_key_exists($name, $cache)) {
        $path = DATA_DIR . '/' . $name;
        $text = is_file($path) ? file_get_contents($path) : false;
        $cache[$name] = $text === false ? null : json_decode($text, true);
    }
    return $cache[$name];
}

function h($text): string
{
    return htmlspecialchars((string) $text, ENT_QUOTES | ENT_SUBSTITUTE, 'UTF-8');
}

// 1234567 -> "1,234,567"
function num($n): string
{
    return $n === null ? '–' : number_format((float) $n);
}

// 0.8175 -> "81.8%"
function pct($x, int $digits = 1): string
{
    return $x === null ? '–' : number_format(100 * (float) $x, $digits) . '%';
}

function bytes_text($n): string
{
    $units = ['bytes', 'KB', 'MB', 'GB', 'TB'];
    $n = (float) $n;
    $i = 0;
    while ($n >= 1000 && $i < count($units) - 1) {
        $n /= 1000;
        $i++;
    }
    return $i === 0 ? number_format($n) . ' bytes' : number_format($n, $n < 10 ? 2 : 1) . ' ' . $units[$i];
}

// "2026-11-09T07:34:41" -> "9 November 2026"
function date_text(?string $iso): string
{
    $t = $iso ? strtotime($iso) : false;
    return $t ? date('j F Y', $t) : '–';
}

// The release's name and a note where it is a trial (a package of some clades alone).
function release_badge(?array $stats): string
{
    if (!$stats) {
        return '';
    }
    $out = '<span class="badge">' . h($stats['package'] ?? 'Eukore') . '</span>';
    if (!empty($stats['partial'])) {
        $clades = $stats['partial']['species_clades'] ?? [];
        $out .= ' <span class="badge badge-trial" title="a trial package: the species of these NCBI taxids\' clades only: '
            . h(implode(', ', $clades)) . '">trial: ' . count($clades) . ' clades</span>';
    }
    return $out;
}

function no_release(string $what): string
{
    return '<div class="notice">No release has been published here yet: ' . h($what) . ' appear once a build\'s '
        . '<code>OUT/release/</code> is copied onto the site.</div>';
}
