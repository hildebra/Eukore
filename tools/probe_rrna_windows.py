"""Why a genome's rRNA windows come out empty (backfill: nhmmer "Unable to guess alphabet"; 2026-11-11).

Read-only. Repeats rrna.extract_genome's two passes (the seed scan, then placement.read_spans) on one genome and
reports each window whose bases come back empty or short, with what the two passes saw of its contig.

    python3 tools/probe_rrna_windows.py GENOME.fna.gz [more genomes...]
"""
import collections
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from eukore import metadata, placement, rrna  # noqa: E402


def titles_of(path):
    """Every FASTA title's first word, and the bases under it, read plainly (gzip text, as read_spans reads)."""
    out, name = [], None
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                name = (line[1:].split() or [""])[0]
                out.append([name, 0, line.rstrip("\n")[:100]])
            elif out:
                out[-1][1] += len(line.strip())
    return out


def probe(path):
    print(f"== {path} ({os.path.getsize(path) / 1e6:.0f} MB)")
    scan, keep, skipped = rrna.SeedScan(), False, []
    for title, data in metadata._pieces_of(path):  # noqa: SLF001
        if title is not None:
            text = title.decode("utf-8", "replace")
            keep = not rrna.NOT_NUCLEAR.search(text)
            if keep:
                scan.start((text.split() or [""])[0])
            else:
                skipped.append(text[:80])
        elif keep:
            scan.feed(data.translate(None, b"\n\r \t"))
    found, lengths = scan.finish()
    print(f"  scan: {len(lengths)} nuclear sequences, {len(skipped)} organelle ones left out, "
          f"{sum(len(w) for w in found.values())} windows on {len(found)} sequences")
    plain = titles_of(path)
    count = collections.Counter(n for n, _l, _t in plain)
    dups = {n: c for n, c in count.items() if c > 1}
    print(f"  file: {len(plain)} titles, {len(dups)} names more than once" + (f" (e.g. {list(dups.items())[:5]})"
                                                                               if dups else ""))
    plain_len = {}
    for n, length, _t in plain:
        plain_len.setdefault(n, []).append(length)
    spans = placement.read_spans(path, {n: [(a + 1, b) for a, b in w] for n, w in found.items()}) if found else {}
    bad = 0
    for name, windows in found.items():
        for a, b in windows:
            got = spans[name][a:b] if name in spans else None
            if got is None or len(got) < b - a:
                bad += 1
                if bad <= 10:
                    print(f"  window {name}:{a}-{b} ({b - a} bp): {len(got) if got is not None else 'no contig'} bp read"
                          f" back; scan's length {lengths.get(name)}, the file's lengths under that name "
                          f"{plain_len.get(name)}")
    print(f"  {bad} windows read back short or empty" + (" (first 10 above)" if bad > 10 else ""))
    if bad:
        named = [t for n, _l, t in plain if n in {k for k, w in found.items()}][:3]
        print(f"  titles of windowed sequences, e.g.: {named}")


for p in sys.argv[1:]:
    probe(p)
