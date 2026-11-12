"""Eukore's benchmarks against other references, and the paper's figures (since 2026-11-12): run on a finished
Eukore package, apart from its build (benchmarks/README.md).

- prep.py: the downloads (PR2, MetaPhlAn 4's database, EukDetect2's database, the read benchmark's genomes and its
  bacterial background) into a folder of their own, and the other references made comparable (references.py);
- compare.py: the package against them (breadth.tsv, overlap.tsv) and the read benchmark (reads.py: reads.tsv,
  detection.tsv), with the package's own benchmarks, into an output folder;
- figures.py: every figure, from that folder's tables.

The modules import the database's code (the eukore package of this checkout): this file puts the checkout first on
the path, so that `python3 -m eukore_bench.<module>` runs from anywhere with benchmarks/ on PYTHONPATH (the scripts
in benchmarks/ do so).
"""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
