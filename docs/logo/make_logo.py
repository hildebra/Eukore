"""The Eukore logos from docs/logo/drafts/eukore_roots_text.svg (the source, its text live): the text as outlines (Lato,
SIL OFL), and the dark-background variant.

    python3 docs/logo/make_logo.py docs/logo /mnt/c/Windows/Fonts     # FONT_DIR holds LatoWeb-Bold.ttf, LatoWeb-Semibold.ttf

Needs fonttools. The PNGs are exported from the SVGs by a browser (2x, transparent).
"""
import re
import sys

from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont

logo_dir, font_dir = sys.argv[1], sys.argv[2]


def outline(font_file, runs, size, spacing, cx, baseline):
    """runs: [(text, fill)] -> <path> elements of the line, centred on cx, its baseline at baseline."""
    font = TTFont(font_file)
    cmap, glyphs, hmtx = font.getBestCmap(), font.getGlyphSet(), font["hmtx"]
    scale = size / font["head"].unitsPerEm
    chars = [(ch, fill) for text, fill in runs for ch in text]
    width = sum(hmtx[cmap[ord(ch)]][0] * scale for ch, _ in chars) + spacing * (len(chars) - 1)
    x, out, current, d = cx - width / 2, [], None, []
    for ch, fill in chars:
        name = cmap[ord(ch)]
        if fill != current and d:
            out.append(f'<path fill="{current}" d="{" ".join(d)}"/>')
            d = []
        current = fill
        pen = SVGPathPen(glyphs)
        glyphs[name].draw(TransformPen(pen, (scale, 0, 0, -scale, x, baseline)))
        if pen.getCommands():
            d.append(pen.getCommands())
        x += hmtx[name][0] * scale + spacing
    out.append(f'<path fill="{current}" d="{" ".join(d)}"/>')
    return out


def rounded(path_d):
    return re.sub(r"-?\d+\.\d+", lambda m: f"{float(m.group()):.2f}".rstrip("0").rstrip("."), path_d)


def wordmark(ink, accent, muted, y_name, y_tag):
    bold, semi = f"{font_dir}/LatoWeb-Bold.ttf", f"{font_dir}/LatoWeb-Semibold.ttf"
    name = outline(bold, [("Euk", ink), ("ore", accent)], 80, 2, 300, y_name)
    tag = outline(semi, [("EUKARYOTE MARKER GENES", muted)], 15, 4.5, 300, y_tag)
    return ("  <!-- the name (Lato Bold and Semibold, SIL Open Font License, as outlines) -->\n  <g>\n    "
            + "\n    ".join(rounded(p) for p in name + tag) + "\n  </g>\n")


src = open(f"{logo_dir}/drafts/eukore_roots_text.svg").read()
head, _ = src.split("  <!-- the name -->")

light = head + wordmark("#1F2A44", "#1E8C8E", "#5B6478", 380, 416) + "</svg>\n"
open(f"{logo_dir}/eukore_logo.svg", "w", newline="\n").write(light)

# dark backgrounds: light branches and name, a brighter database
dark = head.replace("#1F2A44", "#E4EAF3").replace("#16324F", "#24577A").replace("#1D5C73", "#3A8FA8")
dark += wordmark("#F2F5FA", "#4FC9C4", "#A9B3C6", 380, 416) + "</svg>\n"
open(f"{logo_dir}/eukore_logo_dark.svg", "w", newline="\n").write(dark)
print("written: eukore_logo.svg, eukore_logo_dark.svg")
