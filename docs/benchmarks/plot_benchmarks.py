#!/usr/bin/env python3
"""Draws a package's figures into docs/benchmarks/figures/ (those of docs/benchmarks.md and the README). Since
2026-11-10 eukore/figures.py draws them (Figure 1 and the rest; the package stage into the package's
benchmarks/figures/); this script only sets the folder.

    python3 docs/benchmarks/plot_benchmarks.py OUT/Eukore_2026.11                    # into docs/benchmarks/figures
    python3 docs/benchmarks/plot_benchmarks.py PKG --benchmarks DIR --format svg
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from eukore import figures  # noqa: E402

if __name__ == "__main__":
    sys.exit(figures.main(default_out=os.path.join(HERE, "figures")))
