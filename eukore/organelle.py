"""Eukore's organelle genes (since 2026-11-05; docs/FORMAT.md 3): each species' mitochondrial and plastid rRNA (SSU, LSU,
5S), its mitochondrial COX1 and (since 2026-11-07) CYTB, its plastid rbcL, matK and tufA, where NCBI has them
(GENES_OF).

- RefSeq's organelle genomes (its release's mitochondrion and plastid GenBank files, REFSEQ_RELEASE): their annotated
  rRNA and protein-coding genes' features (rrna_gene, cds_gene), parsed once into WORK/organelle (refseq_genes). A
  species takes, of each organelle, its records' genes (choose): those of the one with the most of GENES, then its own
  taxid's (not a strain's or subspecies'), then the longest, and a gene it lacks from the next (a mitogenome in
  several chromosomes). A gene annotated twice (a plastid's inverted repeat) is taken once, `copies` counting them:
  since 2026-11-10 a whole feature before a partial one, then the longest, then the first (until then the first in
  the file, a truncated ycf1 at the inverted repeat's border before the gene); a feature marked /pseudo or
  /pseudogene is none of them.
- Else the species' own assembly (assembly_pass): the molecules NCBI's FASTA titles name an organelle (TITLE:
  mitochondrion, chloroplast, plastid, apicoplast; since 2026-11-10 kinetoplast, maxicircle, minicircle, chromatophore,
  cyanelle). nhmmer finds each gene, queried with that gene of the NEAREST
  RefSeq species in NCBI's lineage (its exons joined); barrnap the rRNA nhmmer does not (--kingdom mito in animals'
  mitochondria, whose 12S and 16S are its models; bac elsewhere, as plastids' rRNA are bacterial in kind). barrnap
  alone finds no rRNA in yeast's mitochondrion, nhmmer with Kluyveromyces' all of it (2026-11-05).
A gene's record is its sequence from its first base to its last, introns and all (fungi's and plants' COX1 hold group
I and II introns), on its strand; one across a circular genome's origin is joined across it (its start after its end).
"""

import collections
import gzip
import hashlib
import os
import re
import subprocess
import tempfile
import time

from . import gencode, rrna, util

REFSEQ_RELEASE = "https://ftp.ncbi.nlm.nih.gov/refseq/release"
RELEASE_FILE = re.compile(r"\b((mitochondrion|plastid)\.\d+\.genomic\.gbff\.gz)\b")
KINDS = ("mitochondrion", "plastid")
GENES = ("SSU", "LSU", "5S", "COX1", "COX3", "CYTB", "nad1", "nad5", "atp6", "rbcL", "matK", "tufA", "psbA", "ycf1",
         "trnH-psbA")
# each organelle's genes, since 2026-11-07 also: a mitochondrion's CYTB (vertebrates' and insects' barcode), a plastid's
# rbcL and matK (land plants' core barcode; rbcL also algae's) and tufA (green algae's barcode); since 2026-11-07
# (format 5) a mitochondrion's COX3, nad1, nad5 and atp6, a plastid's psbA, ycf1 and the trnH-psbA spacer (plants')
GENES_OF = {"mitochondrion": ("SSU", "LSU", "5S", "COX1", "COX3", "CYTB", "nad1", "nad5", "atp6"),
            "plastid": ("SSU", "LSU", "5S", "rbcL", "matK", "tufA", "psbA", "ycf1", "trnH-psbA")}
OLD_GENES = ("SSU", "LSU", "5S", "COX1")  # what an assembly's search found before 2026-11-07 (no #searched line)
RRNA_GENES = ("SSU", "LSU", "5S")
SPACER = "trnH-psbA"
CODING = tuple(g for g in GENES if g not in RRNA_GENES and g != SPACER)  # translated (cds.fna.gz, proteins.faa.gz)
METAZOA = 33208
EMBRYOPHYTA = 3193  # land plants: their mitochondrion's rRNA are 18S and 26S
KINETOPLASTEA = 5653  # their maxicircle's rRNA are 9S (the small subunit's) and 12S (the large one's)
KINETOPLASTID = "kinetoplastid"  # rrna_gene's `animal_mitochondrion` of a kinetoplastid (rrna_naming)
NO_PLASTID = {METAZOA, 4751}  # animals, fungi
COMPARTMENT = {"mitochondrion": "mito", "plastid": "plastid"}  # ncbi.code_of's
TITLE = rrna.ORGANELLE_TITLE
MAX_MOLECULE = 15_000_000  # an organelle molecule this long at most (Silene's mitochondrial chromosomes: ~11 Mb)
NEAREST = 3                # RefSeq species whose gene is queried in an assembly's organelle
MAX_EVALUE = 1e-10         # nhmmer's hits taken
MAX_GENE = 40_000          # a gene's hits this close to its best are one gene (introns: a COX1's up to ~20 kb in fungi)
MIN_COVER = 0.8            # of its query's exons a gene covers at least, else partial
NHMMER_WHOLE = 0.9         # of its query nhmmer's single hit covers at least, else miniprot's model (with_miniprot)
BARRNAP = {"mito": {"12S_rRNA": "SSU", "16S_rRNA": "LSU"},
           "bac": {"16S_rRNA": "SSU", "23S_rRNA": "LSU", "5S_rRNA": "5S"}}
GENE_COLUMNS = ("record", "gene", "organelle", "taxid", "species_taxid", "source", "seqid", "accession", "start",
                "end", "strand", "length", "exons", "copies", "partial", "ambiguous",
                # since 2026-11-07 (format 5): the gene's common name, where its exons come from, its parts on its
                # source, a trans-spliced gene, the spacer's regions, its checks (check)
                "name", "model", "parts", "trans_spliced", "regions", "internal_stops", "own_species", "same_as",
                "length_check", "confidence", "flags",
                "cross")  # since 2026-11-09 (format 6): the other clade's species it is shared with (qc.verdict)
CACHE_COLUMNS = ("seqid", "organelle", "taxid", "record_length", "gene", "start", "end", "strand", "exons", "copies",
                 "partial", "sequence", "spliced", "parts", "trans_spliced", "regions",
                 "transl_table",  # since 2026-11-10: a CDS's own /transl_table, `-` where it names none (table_of)
                 "transl_except")  # since 2026-11-10 (cache 3): its /transl_except codons (transl_excepts), or `-`
# What RefSeq's parsed genes (refseq_genes) and an assembly's search file (assembly_pass) were made by, part of their
# keys since 2026-11-10: to be raised where record_genes', the searches' or the models' logic changes what they hold.
# Until then (1, never written) a cache outlived such a change. 2: /pseudo, the longest feature, /codon_start,
# /transl_table, the spacer across the origin. 3 (2026-11-10), on two branches at once: /transl_except on one; on
# the other a kinetoplastid's 12S its LSU (rrna_gene) and miniprot's model in place of nhmmer's single hit that is
# not the gene's coding sequence (with_miniprot). 4: the two together, so that a cache either branch wrote is made
# again.
CACHE_VERSION = 4
ASSEMBLY_COLUMNS = ("gene", "organelle", "seqid", "start", "end", "strand", "exons", "copies", "partial", "sequence",
                    "parts", "model")  # parts, model: since 2026-11-07


# --- RefSeq's release ---------------------------------------------------------------------------------------------

def folder_of(src_dir):
    return os.path.join(src_dir, "refseq_organelle")


def fetch(src_dir, download, log=print):
    """RefSeq's organelle release (its mitochondrion.N and plastid.N GenBank files, ~2.2 GB in 2026-09) into
    src_dir/refseq_organelle, kept when there; release.txt names them, written last."""
    folder = folder_of(src_dir)
    listing = os.path.join(folder, "release.txt")
    if os.path.isfile(listing):
        return
    os.makedirs(folder, exist_ok=True)
    names = []
    for kind in KINDS:
        index = os.path.join(folder, f"{kind}.html")
        download(f"{REFSEQ_RELEASE}/{kind}/", index)
        with open(index) as fh:
            names += sorted({m.group(1) for m in RELEASE_FILE.finditer(fh.read()) if m.group(2) == kind})
        os.remove(index)
    for name in names:
        dest = os.path.join(folder, name)
        if not os.path.isfile(dest):
            log(f"RefSeq: fetching {name}")
            download(f"{REFSEQ_RELEASE}/{name.split('.')[0]}/{name}", dest)
    with open(util.part_of(listing), "w", newline="\n") as fh:
        fh.writelines(n + "\n" for n in names)
    os.replace(util.part_of(listing), listing)


def release_files(src_dir):
    """The release's files (fetch), or [] where it was not fetched."""
    listing = os.path.join(folder_of(src_dir), "release.txt")
    if not os.path.isfile(listing):
        return []
    with open(listing) as fh:
        return [os.path.join(folder_of(src_dir), n.strip()) for n in fh if n.strip()]


def gbff_records(path):
    """A GenBank flat file -> dict(seqid (its VERSION), taxid, length, circular, sequence, features [(key, [location's
    lines], [(qualifier, value)])]) for each record (its organelle: the release file's, refseq_genes)."""
    rec = feat = qual = None
    section = None
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                f = line.split()
                rec = dict(seqid=f[1], taxid=None, length=int(f[2]) if f[2].isdigit() else 0,
                           circular=" circular " in line, features=[], seq=[])
                section = feat = qual = None
            elif rec is None:
                continue
            elif line.startswith("//"):
                rec["sequence"] = "".join(rec.pop("seq")).upper()
                rec["length"] = rec["length"] or len(rec["sequence"])
                for key, _loc, quals in rec["features"]:
                    if key == "source":
                        for name, value in quals:
                            if name == "db_xref" and value.startswith("taxon:") and value[6:].isdigit():
                                rec["taxid"] = int(value[6:])
                        break
                yield rec
                rec = None
            elif line.startswith("VERSION"):
                rec["seqid"] = (line.split() + [rec["seqid"]])[1]
            elif line.startswith("FEATURES"):
                section = "features"
            elif line.startswith("ORIGIN"):
                section = "origin"
            elif section == "origin":
                rec["seq"].append("".join(line.split()[1:]))
            elif section == "features" and line.startswith("     "):
                key, body = line[5:21].strip(), line[21:].rstrip("\n").strip()
                if key:
                    feat, qual = (key, [body], []), None
                    rec["features"].append(feat)
                elif body.startswith("/"):
                    name, _eq, value = body[1:].partition("=")
                    feat[2].append((name, value.strip('"')))
                    qual = len(feat[2]) - 1
                elif feat is not None and qual is None:
                    feat[1].append(body)
                elif feat is not None:
                    name, value = feat[2][qual]
                    feat[2][qual] = (name, f"{value} {body.strip(chr(34))}")
            elif section == "features":
                section = None


def parts_of(location):
    """A GenBank location -> (parts [(start, end, strand)] 1-based in the order transcribed, partial), or None where a
    part is another record's (a gene spliced across a plant mitochondrion's chromosomes) or a site between bases.
    complement() reverses its parts' order and strands; join() and order() keep theirs (since 2026-11-07: a part on
    either strand, a trans-spliced gene's; until then such a gene was left out)."""
    if ":" in location or "^" in location:
        return None
    text, pos = location.replace("<", "").replace(">", ""), 0

    def expr():
        nonlocal pos
        for op in ("complement(", "join(", "order("):
            if text.startswith(op, pos):
                pos += len(op)
                items = [expr()]
                while text[pos] == ",":
                    pos += 1
                    items.append(expr())
                if text[pos] != ")":
                    raise ValueError(location)
                pos += 1
                parts = [p for item in items for p in item]
                return [(a, b, "-" if s == "+" else "+") for a, b, s in reversed(parts)] if op == "complement(" \
                    else parts
        m = re.match(r"(\d+)(?:\.\.(\d+))?", text[pos:])
        if not m:
            raise ValueError(location)
        pos += m.end()
        return [(int(m.group(1)), int(m.group(2) or m.group(1)), "+")]
    try:
        parts = expr()
    except (ValueError, IndexError):
        return None
    return (parts, "<" in location or ">" in location) if parts and pos == len(text) else None


def trans_spliced(parts, circular=True):
    """Whether a gene's parts (parts_of's, in the order transcribed) are not one stretch of its molecule read along
    one strand (once across a circular molecule's origin at most): a trans-spliced gene, its record its exons joined."""
    if len({s for _a, _b, s in parts}) > 1:
        return True
    sign = 1 if parts[0][2] == "+" else -1
    turns = sum(1 for p, q in zip(parts, parts[1:]) if sign * (q[0] - p[0]) < 0)
    return turns > (1 if circular else 0)


def span_of(spans, length, circular=True):
    """A feature's parts -> (start, end), 1-based: its first base to its last, introns between; the start after the end
    where it crosses a circular genome's origin (its largest gap is not the one across the origin)."""
    spans = sorted(spans)
    first, last = spans[0][0], max(e for _s, e in spans)
    if len(spans) == 1 or not circular or not length:
        return first, last
    gap, i = max((spans[j + 1][0] - spans[j][1], j) for j in range(len(spans) - 1))
    if length - last + first >= gap:
        return first, last
    return spans[i + 1][0], spans[i][1]


def subsequence(seq, start, end, strand):
    """seq's 1-based start..end (across the origin where start > end), on the strand."""
    s = seq[start - 1:end] if start <= end else seq[start - 1:] + seq[:end]
    return rrna.revcomp(s) if strand == "-" else s


def spliced(seq, parts):
    """A feature's parts (parts_of's) joined in the order transcribed."""
    return "".join(subsequence(seq, a, b, s) for a, b, s in parts)


def parts_text(parts):
    """Parts [(start, end, strand)] -> `start-end:strand,...` (genes.tsv's `parts`)."""
    return ",".join(f"{a}-{b}:{s}" for a, b, s in parts)


def parts_read(text):
    """genes.tsv's `parts` -> [(start, end, strand)], [] for `-`."""
    out = []
    for p in (text or "-").split(","):
        if p and p != "-":
            span, _c, strand = p.partition(":")
            a, _d, b = span.partition("-")
            out.append((int(a), int(b), strand or "+"))
    return out


NUMBER = re.compile(r"\brrn(\d+(?:\.\d+)?)\b|\b(\d+(?:\.\d+)?)\s?s\b")
SMALL = re.compile(r"small subunit|\bs-rrna\b|\brrns\b|\brns\b|\bssu\b")
LARGE = re.compile(r"large subunit|\bl-rrna\b|\brrnl\b|\brnl\b|\blsu\b")
COX1_PRODUCT = re.compile(r"cytochrome (c )?oxidase subunit (1|i)\b")
COX1_GENES = {"cox1", "coi", "co1", "coxi", "cox-1", "mt-co1", "mtco1", "cox1a"}


def rrna_gene(text, animal_mitochondrion=False):
    """An rRNA feature's /gene and /product, lower case -> SSU, LSU, 5S or None (4.5S, unknown). "16S" is an animal
    mitochondrion's LSU (its SSU is the 12S), elsewhere the SSU. animal_mitochondrion is how the mitochondrion's rRNA
    are numbered (rrna_naming): True in an animal and, since 2026-11-10, KINETOPLASTID in a kinetoplastid, whose
    maxicircle's 12S is the LARGE subunit's rRNA and whose 9S is the small one's. Until then both were read as the
    SSU: the species' "SSU" was its 12S (the longer feature of the two, `copies` 2) and it held no LSU (RefSeq's
    NC_000894, Leishmania tarentolae: its "SSU" the 12S of 1,173 bases)."""
    if SMALL.search(text):
        return "SSU"
    if LARGE.search(text):
        return "LSU"
    m = NUMBER.search(text)
    if not m:
        return None
    s = float(m.group(1) or m.group(2))
    if s == 5:
        return "5S"
    if s < 9:
        return None
    kinetoplastid = animal_mitochondrion == KINETOPLASTID
    if (s == 12 and kinetoplastid) or (s == 16 and animal_mitochondrion and not kinetoplastid):
        return "LSU"
    return "SSU" if s < 20 else "LSU"


RBCL_PRODUCT = re.compile(r"ribulose[- ]1,5[- ]bis-?phosphate carboxylase(/oxygenase)? large|"
                          r"ribulose bisphosphate carboxylase large|rubisco large")
MATK_PRODUCT = re.compile(r"\bmaturase k\b")
CYTB_GENES = {"cytb", "cob", "cyb", "mt-cyb", "mtcyb"}
CYTB_PRODUCT = re.compile(r"\b(apo)?cytochrome b\b(?![-\d])")  # not b5, b6 (petB), b559, b-245
TUFA_PRODUCT = re.compile(r"\belongation factor tu\b")
# the protein-coding genes (since 2026-11-07 one table, CDS_RULES): gene, organelle, its /gene names (lower case), its
# /product, and whether the product names it beside another /gene (COX1's, as until then)
CDS_RULES = (
    ("COX1", "mitochondrion", COX1_GENES, COX1_PRODUCT, True),
    ("COX3", "mitochondrion", {"cox3", "coiii", "co3", "coxiii", "cox-3", "mt-co3", "mtco3"},
     re.compile(r"cytochrome (c )?oxidase subunit (3|iii)\b"), False),
    ("CYTB", "mitochondrion", CYTB_GENES, CYTB_PRODUCT, False),
    ("nad1", "mitochondrion", {"nad1", "nd1", "mt-nd1", "nadh1", "nadh-1"},
     re.compile(r"\bnadh (dehydrogenase|ubiquinone oxidoreductase)( chain| subunit) (1|i)\b(?!-)"), False),
    ("nad5", "mitochondrion", {"nad5", "nd5", "mt-nd5", "nadh5", "nadh-5"},
     re.compile(r"\bnadh (dehydrogenase|ubiquinone oxidoreductase)( chain| subunit) (5|v)\b(?!-)"), False),
    ("atp6", "mitochondrion", {"atp6", "mt-atp6", "atpase6", "atp-6", "atpase 6"},
     re.compile(r"\batp synthase (f0 |fo )?(subunit|chain) (6|a)\b"), False),
    ("rbcL", "plastid", {"rbcl"}, RBCL_PRODUCT, False),
    ("matK", "plastid", {"matk"}, MATK_PRODUCT, False),
    ("tufA", "plastid", {"tufa"}, TUFA_PRODUCT, False),
    ("psbA", "plastid", {"psba"}, re.compile(r"photosystem ii (protein d1|d1 protein|reaction center protein d1)\b"),
     False),
    ("ycf1", "plastid", {"ycf1", "tic214"}, re.compile(r"\bycf1\b|\btic214\b"), False),
)
NAMES = {("mitochondrion", "SSU", True): "12S", ("mitochondrion", "LSU", True): "16S",
         ("mitochondrion", "SSU", False): "SSU", ("mitochondrion", "LSU", False): "LSU",
         ("plastid", "SSU", False): "16S", ("plastid", "LSU", False): "23S"}


def name_of(organelle, gene, animal=False, plant=False):
    """A gene's common name (genes.tsv's `name`): an animal mitochondrion's 12S and 16S, a land plant mitochondrion's
    18S and 26S, a plastid's 16S and 23S; a fungus' or a protist's mitochondrial rRNA stays SSU or LSU (rns, rnl)."""
    if organelle == "mitochondrion" and plant and gene in ("SSU", "LSU"):
        return "18S" if gene == "SSU" else "26S"
    return NAMES.get((organelle, gene, animal and organelle == "mitochondrion"), gene)


def cds_gene(quals, organelle):
    """A CDS feature's gene of GENES_OF[organelle] (CDS_RULES): by its /gene, else (none named) its /product; or None."""
    gene = {v.lower() for n, v in quals if n == "gene"}
    product = " ".join(v.lower() for n, v in quals if n == "product")
    for name, kind, names, products, by_product in CDS_RULES:
        if kind == organelle and (gene & names or ((by_product or not gene) and products.search(product))):
            return name
    return None


def is_his_trna(quals):
    """A tRNA feature of histidine's (trnH-GUG: the trnH-psbA spacer's end)."""
    text = " ".join(v.lower() for n, v in quals if n in ("gene", "product"))
    return bool(re.search(r"\btrnh\b|trnh-|trna-his\b", text))


SPACER_FLANK = 200    # of psbA's 3' end the trnH-psbA record holds (the barcode's psbA primer sits there)
SPACER_MAX = 3000     # trnH and psbA this far apart at most: neighbours


def in_frame(parts, quals):
    """A CDS's parts (parts_of's) from its first whole codon on: its /codon_start 2 or 3 (a CDS whose 5' end is not
    known begins inside a codon) takes that many bases less one off its first part. Since 2026-11-10: until then
    /codon_start was not read, and such a gene's coding sequence and protein were read out of frame."""
    shift = next((int(v) - 1 for n, v in quals if n == "codon_start" and v in ("2", "3")), 0)
    a, b, strand = parts[0]
    if not shift or b - a + 1 <= shift:
        return parts
    return [(a + shift, b, strand) if strand == "+" else (a, b - shift, strand)] + list(parts[1:])


TRANSL_EXCEPT = re.compile(r"pos:(complement\()?(\d+)(?:\.\.(\d+))?\)?,\s*aa:(\w+)")
EXCEPT_AA = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
             "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
             "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U", "Pyl": "O", "Other": "X", "Xaa": "X"}


def transl_excepts(parts, quals):
    """A CDS's /transl_except codons read otherwise than by its table (`(pos:a..b,aa:Trp)`, `(pos:complement(a..b),
    aa:Sec)`; since 2026-11-10, until then not read and such a codon a stop inside the protein) -> `codon:aa,...` (the
    codon's 0-based number in its coding sequence (parts: in_frame's), its amino acid's letter), `-` for none. A TERM
    (an animal mitochondrion's stop completed by the mRNA's poly-A: a codon cut short at the gene's end) is no amino
    acid, and a codon not on the coding sequence's frame is left out."""
    out = []
    for value in (v for n, v in quals if n == "transl_except"):
        for minus, a, b, aa in TRANSL_EXCEPT.findall(value):
            if aa not in EXCEPT_AA:
                continue
            first = int(b or a) if minus else int(a)  # the codon's first base, in the order transcribed
            before = 0
            for x, y, s in parts:
                if x <= first <= y and (s == "-") == bool(minus):
                    offset = before + (y - first if s == "-" else first - x)
                    if offset % 3 == 0:
                        out.append(f"{offset // 3}:{EXCEPT_AA[aa]}")
                    break
                before += y - x + 1
    return ",".join(out) or "-"


def with_excepts(protein, stops, excepts):
    """A protein (gencode.protein_of's, its internal stops) with its /transl_except codons' amino acids (transl_excepts'
    text) in place of the stops their codons read as -> (protein, internal stops)."""
    if not excepts or excepts == "-":
        return protein, stops
    p = list(protein)
    for item in excepts.split(","):
        i, _colon, aa = item.partition(":")
        if i.isdigit() and int(i) < len(p) and p[int(i)] == "*":
            p[int(i)] = aa
    p = "".join(p)
    return p, p.count("*")


def record_genes(rec, organelle, animal=False):
    """A RefSeq organelle record's genes -> {gene: dict(start, end, strand, exons, copies, partial, sequence, spliced
    (its exons joined), parts, trans_spliced, regions, transl_table)}. A trans-spliced gene's sequence is its exons
    joined; the trnH-psbA spacer (since 2026-11-07) runs from trnH-GUG's first base to SPACER_FLANK bases into psbA's
    3' end. Of a gene's features (a plastid's inverted repeat holds it twice) a whole one is taken before a partial
    one, then the longest (its exons joined), then the first in the file; a CDS begins at its first whole codon
    (in_frame)."""
    found, counts, his = {}, collections.Counter(), []
    seq = rec["sequence"]
    for key, location, quals in rec["features"]:
        location = "".join(location)
        if key == "rRNA":
            text = " ".join(v.lower() for n, v in quals if n in ("gene", "product"))
            gene = rrna_gene(text, organelle == "mitochondrion" and animal) or \
                rrna_gene(" ".join(v.lower() for n, v in quals if n == "note"), organelle == "mitochondrion" and animal)
        elif key == "CDS":
            gene = cds_gene(quals, organelle)
        elif key == "tRNA" and organelle == "plastid" and is_his_trna(quals):
            parsed = parts_of(location)
            if parsed:
                his.append(parsed[0])
            continue
        else:
            continue
        parsed = parts_of(location) if gene else None
        # since 2026-11-10: a feature RefSeq marks /pseudo or /pseudogene is not a copy of the gene (a truncated ycf1
        # at the inverted repeat's border, a fragment of a broken gene); until then it was counted, and where it
        # came first in the file it was the species' record and the assemblies' query
        if not parsed or any(n in ("pseudo", "pseudogene") for n, _v in quals):
            continue
        counts[gene] += 1
        parts, partial = parsed
        if max(b for _a, b, _s in parts) > len(seq):
            continue
        if key == "CDS":
            parts = in_frame(parts, quals)
        trans = trans_spliced(parts, rec["circular"])
        joined = spliced(seq, parts)
        if trans:  # its exons on different strands or out of order: the record is its exons joined
            start, end, strand, whole = min(a for a, _b, _s in parts), max(b for _a, b, _s in parts), parts[0][2], joined
        else:
            strand = parts[0][2]
            start, end = span_of([(a, b) for a, b, _s in parts], rec["length"], rec["circular"])
            whole = subsequence(seq, start, end, strand)
        # the feature's own genetic code (a CDS's /transl_table; table_of reads it before NCBI's of the taxid)
        table = next((v for n, v in quals if n == "transl_table" and v.isdigit()), "-") if key == "CDS" else "-"
        row = dict(start=start, end=end, strand=strand, exons=len(parts), partial=partial, sequence=whole,
                   spliced=joined, parts=parts_text(parts), trans_spliced=trans, regions="-", transl_table=table,
                   transl_except=transl_excepts(parts, quals) if key == "CDS" else "-")
        # since 2026-11-10 the gene's best feature, not its first: whole before partial, then the longest; the first
        # of equals (an inverted repeat's two copies), as until then
        if gene not in found or (not partial, len(joined)) > (not found[gene]["partial"], len(found[gene]["spliced"])):
            found[gene] = row
    for gene, row in found.items():
        row["copies"] = counts[gene]
    if "psbA" in found and his and not found["psbA"]["trans_spliced"]:
        spacer = spacer_of(seq, found["psbA"], his, rec["circular"])
        if spacer:
            found[SPACER] = spacer
    return found


def spacer_of(seq, psba, his, circular=False):
    """The trnH-psbA record of a plastid record (its psbA's row, its trnH-GUG parts) -> a gene row, or None where
    trnH is not psbA's neighbour on its strand. psbA's 3' end is its end on its strand; trnH lies beyond it. On a
    circular record the origin may lie anywhere between psbA's 3' end and trnH (or inside psbA): the record is
    joined across it, its start after its end, as a gene's is (since 2026-11-10; until then such a record, opened
    between trnH and psbA, gave no spacer)."""
    strand, length = psba["strand"], len(seq)
    if psba["start"] > psba["end"] and not circular:
        return None
    three = psba["start"] if strand == "-" else psba["end"]  # psbA's last base, on the molecule
    psba_length = (psba["end"] - psba["start"]) % length + 1
    for parts in his:
        if len(parts) != 1 or parts[0][2] != strand:
            continue
        a, b, _s = parts[0]
        near = three - b if strand == "-" else a - three  # trnH's near end beyond psbA's last base, along its strand
        if circular and not 0 < near <= SPACER_MAX:
            # the way round the origin; on a small circle a trnH before psbA is also "beyond" it that way: taken only
            # where it is nearer to psbA's 3' end than, the other way round, to its 5' end
            near %= length
            if length - psba_length - (near - 1) - (b - a + 1) < near - 1:
                continue
        if not 0 < near <= SPACER_MAX:
            continue
        size = SPACER_FLANK + near - 1 + b - a + 1  # psbA's 3' end, the spacer, trnH
        start, end = (a, a + size - 1) if strand == "-" else (b - size + 1, b)
        if circular and size < length:  # an end beyond the record's: back on it, the start after the end
            start, end = (start - 1) % length + 1, (end - 1) % length + 1
        elif start < 1 or end > length:
            continue
        regions = [("psbA3'", 0, SPACER_FLANK), ("spacer", SPACER_FLANK, SPACER_FLANK + near - 1),
                   ("trnH", SPACER_FLANK + near - 1, size)]
        whole = subsequence(seq, start, end, strand)
        pieces = [(start, end, strand)] if start <= end else \
            [(start, length, strand), (1, end, strand)][::1 if strand == "+" else -1]  # in the order transcribed
        return dict(start=start, end=end, strand=strand, exons=1, partial=False, sequence=whole, spliced=whole,
                    parts=parts_text(pieces), trans_spliced=False, copies=1, transl_table="-", transl_except="-",
                    regions=",".join(f"{n}:{x}-{y}" for n, x, y in regions))
    return None


def is_animal(ncbi, taxid):
    return taxid is not None and METAZOA in ncbi.ancestors(ncbi.merged.get(taxid, taxid))


def rrna_naming(ncbi, taxid):
    """How a taxon's mitochondrial rRNA are numbered, as rrna_gene reads it (record_genes' `animal`, which hands it on
    for a mitochondrion's features): KINETOPLASTID under Kinetoplastea (since 2026-11-10: their 12S is the LSU), else
    whether it is an animal (its 16S is the LSU), as until then."""
    if taxid is not None and KINETOPLASTEA in ncbi.ancestors(ncbi.merged.get(taxid, taxid)):
        return KINETOPLASTID
    return is_animal(ncbi, taxid)


def table_of(row, ncbi):
    """The genetic code a gene's coding sequence is read by: its RefSeq CDS's own /transl_table where it names one of
    gencode's tables (since 2026-11-10: the code its annotation was made with; until then never read), else NCBI's
    of its taxid's organelle, as before (an assembly's gene; a CDS that names none)."""
    own = str(row.get("transl_table") or "-")
    if own.isdigit() and int(own) in gencode.TABLES:
        return int(own)
    return ncbi.code_of(int(row["taxid"]), COMPARTMENT[row["organelle"]])


def release_sizes(files):
    """The release's files' names and sizes, as text (refseq_genes' key: their content's stand-in; an md5 of the
    ~2.2 GB at each build's start would cost more than a release of the same sizes is likely)."""
    return ",".join(f"{os.path.basename(p)}:{os.path.getsize(p) if os.path.isfile(p) else '-'}" for p in files)


def rules_text():
    """What a RefSeq record's genes are named and cut by, as text (refseq_genes' key): the protein-coding genes' rules,
    the rRNA's patterns, the spacer's lengths. A rule changed parses the release again without CACHE_VERSION raised."""
    rules = [(name, kind, sorted(names), products.pattern, by_product)
             for name, kind, names, products, by_product in CDS_RULES]
    return repr((rules, NUMBER.pattern, SMALL.pattern, LARGE.pattern, SPACER_FLANK, SPACER_MAX))


MADE_POLL = 30  # seconds between a waiting build's looks for the file another build makes (made_once)


def made_once(path, make, log=print, poll=MADE_POLL):
    """path made by make() where it is not there, by one build at a time of all builds on its folder (since 2026-11-11:
    a SLURM array's tasks run the organelle pass, and each parsed RefSeq's ~2.2 GB release at once): the build that
    takes the lock (util.single_run: folder/.<name>.lock) makes it, the others wait for it; a waiting build whose maker
    ended without it (killed) makes it in its place (a dead build's lock is taken over)."""
    folder, name = os.path.split(path)
    os.makedirs(folder, exist_ok=True)
    said = False
    while not os.path.isfile(path):
        try:
            with util.single_run(folder, name):
                if not os.path.isfile(path):  # made by another meanwhile
                    make()
            return
        except util.RunHeld as e:
            if not said:
                log(f"{name}: being made by another build ({e.host or '?'}, process {e.pid or '?'}): waiting for it")
                said = True
            time.sleep(poll)


def refseq_genes(src_dir, work_folder, ncbi, log=print):
    """Every gene of RefSeq's organelle release (record_genes), parsed once into WORK/organelle/refseq_<release>.tsv.gz
    (CACHE_COLUMNS) -> [row dicts], [] where the release was not fetched."""
    files = release_files(src_dir)
    if not files:
        return []
    # the release's and the genes': parsed again where either changes (CYTB, rbcL, matK, tufA since 2026-11-07; the
    # format-5 genes and the cache's parts, trans-spliced genes and regions since 2026-11-07). Since 2026-11-10 the
    # key also holds the release's files' sizes (release.txt names them alone, and RefSeq's names are the same in
    # every release: a release fetched anew was read from the old cache), the rules the genes are named by
    # (rules_text) and CACHE_VERSION. The caches of before are parsed again once (~2.2 GB of GenBank files), which is
    # accepted; their files (refseq_<old key>.tsv.gz) stay in WORK/organelle and can be deleted.
    key = hashlib.md5((util.md5(os.path.join(folder_of(src_dir), "release.txt")) + ",".join(GENES + CACHE_COLUMNS)
                       + release_sizes(files) + rules_text() + f"version:{CACHE_VERSION}").encode()).hexdigest()[:10]
    cache = os.path.join(work_folder, f"refseq_{key}.tsv.gz")

    def parse():
        log(f"RefSeq organelles: parsing {len(files)} GenBank files")
        with util.bgzf_open(util.part_of(cache), level=util.WORK_LEVEL) as fh:
            fh.write("\t".join(CACHE_COLUMNS) + "\n")
            for path in files:
                organelle = "mitochondrion" if os.path.basename(path).startswith("mitochondrion") else "plastid"
                for rec in gbff_records(path):
                    genes = record_genes(rec, organelle, rrna_naming(ncbi, rec["taxid"]))
                    for gene, g in genes.items():
                        row = (rec["seqid"], organelle, rec["taxid"] or "-", rec["length"], gene, g["start"], g["end"],
                               g["strand"], g["exons"], g["copies"], "t" if g["partial"] else "f", g["sequence"],
                               g["spliced"], g["parts"], "t" if g["trans_spliced"] else "f", g["regions"],
                               g["transl_table"], g["transl_except"])
                        fh.write("\t".join(map(str, row)) + "\n")
        os.replace(util.part_of(cache), cache)
    made_once(cache, parse, log)
    rows = []
    with gzip.open(cache, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["taxid"].isdigit():
                rows.append(r)
    return rows


def choose(rows, ncbi):
    """RefSeq's gene rows -> {(species taxid, organelle): [its rows]}: of each species' records of an organelle, the one
    with the most genes, then its own taxid's, then the longest, then the first by accession; a gene it lacks from the
    next records in that order that hold it (a mitogenome in several chromosomes: lice's minicircles, Silene's; until
    2026-11-06 one record's genes alone)."""
    records = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        sp = ncbi.species_of(int(r["taxid"]))
        records[(sp, r["organelle"])][r["seqid"]].append(r)
    out = {}
    for (sp, organelle), by_seqid in records.items():
        ranked = sorted(by_seqid.items(), key=lambda kv: (-len(kv[1]), int(kv[1][0]["taxid"]) != sp,
                                                          -int(kv[1][0]["record_length"]), kv[0]))
        taken, held = [], set()
        for _seqid, genes in ranked:
            new = [r for r in genes if r["gene"] not in held]
            taken += new
            held |= {r["gene"] for r in new}
        out[(sp, organelle)] = taken
    return out


# --- a species' own assembly ----------------------------------------------------------------------------------------

def saved_molecules(path):
    """The organelle molecules metadata.describe kept of a genome (WORK/done/<accession>/organelle.fna.gz: >name
    organelle) -> [(name, organelle, sequence)], as organelle_molecules finds them in the genome."""
    out = []
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].split()
                out.append([f[0] if f else "", f[1] if len(f) > 1 else "mitochondrion", []])
            elif out:
                out[-1][2].append(line.strip())
    return [(n, k, "".join(p).upper()) for n, k, p in out]


def organelle_molecules(path):
    """A genome's organelle molecules (TITLE in their FASTA titles) -> [(name, organelle, sequence)]."""
    out, cur = [], None
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                kind = rrna.organelle_of(line, TITLE)  # since 2026-11-10 one reading of the title's word
                cur = None
                if kind:
                    cur = (line[1:].split() or [""])[0], kind, []
                    out.append(cur)
            elif cur is not None:
                cur[2].append(line.strip())
    return [(n, k, "".join(p).upper()) for n, k, p in out if sum(map(len, p)) <= MAX_MOLECULE]


class Nearest:
    """RefSeq's genes (whole, of the chosen records; their exons joined), by organelle, gene and NCBI lineage."""

    def __init__(self, chosen, ncbi):
        self.ncbi, self.seq, self.protein = ncbi, {}, {}
        self.under = collections.defaultdict(lambda: collections.defaultdict(list))
        for (sp, organelle), rows in sorted(chosen.items()):
            for r in rows:
                if r["partial"] == "f" and r.get("trans_spliced", "f") != "t":
                    self.seq[(r["seqid"], r["gene"])] = r["spliced"]
                    if r["gene"] in CODING:  # its protein, by its own code (table_of; miniprot_genes' query)
                        p, stops = with_excepts(*gencode.protein_of(r["spliced"], table_of(
                            dict(r, organelle=organelle), ncbi)), r.get("transl_except"))  # since 2026-11-10
                        if not stops:
                            self.protein[(r["seqid"], r["gene"])] = p
                    for a in ncbi.ancestors(sp):
                        self.under[(organelle, r["gene"])][a].append(r["seqid"])

    def queries(self, taxid, organelle, n=NEAREST, exclude=()):
        """Each gene of the n RefSeq species nearest to taxid with it -> [(gene|seqid, its exons)]; exclude: seqids
        not taken (tools/organelle_models.py: a record's own, searched as if it were an assembly's)."""
        out = []
        for gene in GENES_OF[organelle]:
            under, got = self.under.get((organelle, gene), {}), []
            for a in self.ncbi.ancestors(self.ncbi.merged.get(taxid, taxid)):
                for s in under.get(a, ()):
                    if len(got) >= n:
                        break
                    if s not in got and s not in exclude:
                        got.append(s)
                if len(got) >= n:
                    break
            out += [(f"{gene}|{s}", self.seq[(s, gene)]) for s in got]
        return out

    def proteins(self, taxid, organelle, n=NEAREST, exclude=()):
        """Each protein-coding gene's protein of the n RefSeq species nearest to taxid with it -> {gene~seqid:
        protein} (miniprot's queries: loci.parse_alignments reads the gene before the ~)."""
        return {f"{q.split('|')[0]}~{q.split('|')[1]}": self.protein[(q.split("|")[1], q.split("|")[0])]
                for q, _s in self.queries(taxid, organelle, n, exclude)
                if (q.split("|")[1], q.split("|")[0]) in self.protein}


def barrnap_genes(molecules, barrnap, kingdom, threads=1):
    """barrnap over molecules [(name, organelle, sequence)] -> {(organelle, gene): (rows, copies)}, the longest whole
    hit of each gene first."""
    if not molecules:
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "m.fna")
        with open(fasta, "w", newline="\n") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, _k, s in molecules)
        feats = rrna.parse_gff(rrna.barrnap_gff(barrnap, fasta, threads, kingdom))
    seqs = {n: (k, s) for n, k, s in molecules}
    hits = collections.defaultdict(list)
    for name, by_model in feats.items():
        organelle, seq = seqs[name]
        for model, found in by_model.items():
            gene = BARRNAP[kingdom].get(model)
            for start, end, strand, partial, _share in found if gene else ():
                hits[(organelle, gene)].append(dict(gene=gene, organelle=organelle, seqid=name, start=start, end=end,
                                                    strand=strand, exons=1, partial=partial,
                                                    sequence=subsequence(seq, start, end, strand),
                                                    parts=f"{start}-{end}:{strand}", model="barrnap"))
    return {k: (min(v, key=lambda h: (h["partial"], h["start"] - h["end"], h["seqid"], h["start"])), len(v))
            for k, v in hits.items()}


def collinear(best, hits, gap=MAX_GENE):
    """The best hit's chain among hits of its molecule, query and strand: further along the query is further along the
    gene's strand, each within gap of the last. A hit out of that order is left out (until 2026-11-06 every hit within
    MAX_GENE was joined: a mitochondrion opened inside COX1 shipped whole as its "COX1"; a concatemer's next copy)."""
    sign = 1 if best["strand"] == "+" else -1
    chain = [best]
    for step in (1, -1):  # toward the query's end, then toward its start
        cur = best
        for h in sorted(hits, key=lambda h: step * sign * h["start"]):
            if any(h is c for c in chain):
                continue
            along = step * sign * (h["start"] - cur["start"]) > 0 and step * sign * (h["end"] - cur["end"]) > 0
            later = step * (h["q1"] - cur["q1"]) > 0 and step * (h["q2"] - cur["q2"]) > 0
            if along and later and h["start"] <= cur["end"] + gap and h["end"] >= cur["start"] - gap:
                chain.append(h)
                cur = h
    return chain


def nhmmer_genes(molecules, queries, nhmmer, threads=1):
    """nhmmer of RefSeq's genes (queries [(gene|seqid, its exons)], Nearest.queries) over one organelle's molecules
    [(name, organelle, sequence)] -> {gene: row}: of each gene, its best hit's query's hits on its molecule and strand
    within MAX_GENE of it, joined, its ends reaching the query's (introns between)."""
    if not molecules or not queries:
        return {}
    organelle = molecules[0][1]
    with tempfile.TemporaryDirectory() as tmp:
        q, t, tbl = (os.path.join(tmp, f) for f in ("q.fna", "t.fna", "hits.tbl"))
        for path, records in ((q, queries), (t, [(n, s) for n, _k, s in molecules])):
            with open(path, "w", newline="\n") as fh:
                fh.writelines(f">{n}\n{s}\n" for n, s in records)
        run = subprocess.run([nhmmer, "--dna", "--noali", "--cpu", str(max(1, threads)), "--tblout", tbl, "-o",
                              os.devnull, q, t], capture_output=True, text=True, check=False)
        if run.returncode != 0:
            raise RuntimeError(f"nhmmer failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        hits = []
        with open(tbl) as fh:
            for line in fh:
                f = line.split()
                if line.startswith("#") or len(f) < 14 or float(f[12]) > MAX_EVALUE:
                    continue
                a, b = sorted((int(f[6]), int(f[7])))
                hits.append(dict(target=f[0], query=f[2], q1=int(f[4]), q2=int(f[5]), start=a, end=b, strand=f[11],
                                 score=float(f[13])))
    seqs, q_lens = {n: s for n, _k, s in molecules}, {n: len(s) for n, s in queries}
    by_gene = collections.defaultdict(list)
    for h in hits:
        by_gene[h["query"].split("|", 1)[0]].append(h)
    out = {}
    for gene, gene_hits in by_gene.items():
        best = max(gene_hits, key=lambda h: (h["score"], h["target"], -h["start"]))
        seq = seqs[best["target"]]
        mine = [h for h in gene_hits if (h["target"], h["query"], h["strand"]) ==
                (best["target"], best["query"], best["strand"])]
        same = collinear(best, mine + across_origin(mine, len(seq)))
        q_len = q_lens[best["query"]]
        covered = set()
        for h in same:
            covered.update(range(h["q1"], h["q2"] + 1))
        start, end = gene_span(same, best["strand"], q_len, len(seq))
        out[gene] = dict(gene=gene, organelle=organelle, seqid=best["target"], start=start, end=end,
                         strand=best["strand"], exons=len(same), copies=1, partial=len(covered) < MIN_COVER * q_len,
                         sequence=subsequence(seq, start, end, best["strand"]),
                         parts=f"{start}-{end}:{best['strand']}" if len(same) == 1 else "-", model="nhmmer",
                         # not a column of the search file: the share of its query the hits cover (with_miniprot)
                         cover=len(covered) / max(1, q_len))
    return out


def miniprot_genes(molecules, proteins, miniprot, table=1, threads=1):
    """miniprot of RefSeq's proteins (Nearest.proteins) over one organelle's molecules by the organelle's genetic code
    (since 2026-11-07) -> {gene: row}: each protein-coding gene's best alignment, its exons miniprot's (exact, where
    nhmmer's hits' ends are not), partial where it covers under MIN_COVER of its query. A gene across a circular
    molecule's origin, or one miniprot does not align, is left to nhmmer's."""
    if not molecules or not proteins or not miniprot:
        return {}
    from . import loci  # noqa: PLC0415
    organelle, seqs = molecules[0][1], {n: s for n, _k, s in molecules}
    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "m.fna")
        with open(fasta, "w", newline="\n") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in seqs.items())
        alignments, why = loci.align(fasta, proteins, miniprot, threads, table, tmp)
    if why:
        raise RuntimeError(why)
    out = {}
    for a in sorted(alignments, key=lambda a: (-a["score"], a["contig"], a["start"])):
        gene = a["locus"]
        if gene in out or a["contig"] not in seqs or not a["exons"]:
            continue
        exons = sorted(a["exons"], reverse=a["strand"] == "-")  # in the order transcribed
        parts = [(x, y, a["strand"]) for x, y in exons]
        cover = (a["qend"] - a["qstart"]) / max(1, len(proteins[a["query"]]))
        start, end = min(x for x, _y in exons), max(y for _x, y in exons)
        out[gene] = dict(gene=gene, organelle=organelle, seqid=a["contig"], start=start, end=end, strand=a["strand"],
                         exons=len(exons), copies=1, partial=cover < MIN_COVER,
                         sequence=subsequence(seqs[a["contig"]], start, end, a["strand"]), parts=parts_text(parts),
                         model="miniprot")
    return out


def across_origin(hits, length):
    """Hits of a molecule of `length` bases as they lie one turn on and one back (its coordinates + and - length): a
    circular molecule opened inside a gene holds the gene's other end at its start (collinear chains them; until
    2026-11-06 such a gene was its half by the origin, partial)."""
    return [dict(h, start=h["start"] + k, end=h["end"] + k) for k in (length, -length) for h in hits]


def gene_span(chain, strand, q_len, length):
    """A gene's (start, end) on its molecule, 1-based, from its chained hits (collinear's; coordinates beyond the
    molecule where the chain crosses its origin), reaching the query's ends the hits leave out: start after end where
    it crosses the origin (subsequence joins it), else clipped to the molecule. The ends are reached however far
    (the audit's B2 proposed a few codons at most; measured 2026-11-10 on 190 RefSeq organelle genomes searched as
    assemblies, with_miniprot: of 642 genes without introns 426 came out exact and 86 over 6 codons off; with the
    ends reached by 9 bases at most, 401 and 117, by none where over 9, 400 and 129 - nhmmer's local hits stop
    short of a diverged end, and the query's length is the better guess of where the gene ends)."""
    q_first, q_last = min(h["q1"] for h in chain), max(h["q2"] for h in chain)
    start, end = min(h["start"] for h in chain), max(h["end"] for h in chain)
    lead, tail = q_first - 1, q_len - q_last  # the query's ends the hits leave out: the gene reaches them
    start, end = (start - lead, end + tail) if strand == "+" else (start - tail, end + lead)
    crossing = any(h["start"] < 1 or h["end"] > length for h in chain)
    if crossing and end - start + 1 < length:  # across the origin: each end back on the molecule
        return (start + length if start < 1 else start), (end - length if end > length else end)
    return max(1, start), min(length, end)


def assembly_file(work_folder, accession):
    return os.path.join(work_folder, "assembly", f"{accession}.tsv")


def search_key(organelle, animal, queries, proteins, code, models):
    """One organelle's search's inputs as a short digest (the #searched line's `inputs`, since 2026-11-10): the
    queries (the nearest RefSeq species' genes: other ones after a new release, or after their parsing changed),
    barrnap's kingdom, the models, with miniprot its proteins and the genetic code, the thresholds and CACHE_VERSION.
    Until then a file was kept by its genes' names and its models alone, whatever it had been searched with."""
    kingdom = "mito" if animal and organelle == "mitochondrion" else "bac"
    by_miniprot = (sorted(proteins.items()), code) if models == "miniprot" else ()
    text = repr((CACHE_VERSION, organelle, GENES_OF[organelle], kingdom, models, list(queries), by_miniprot,
                 MAX_EVALUE, MAX_GENE, MIN_COVER, MAX_MOLECULE, NHMMER_WHOLE))
    return hashlib.md5(text.encode()).hexdigest()[:10]


def searched_of(path):
    """An assembly's search file's genes searched -> {organelle: set of genes}: its #searched line's (and its other
    entries: `models`, since 2026-11-10 `inputs` {organelle=search_key}), else (a file of before 2026-11-07) OLD_GENES
    of each organelle."""
    with open(path) as fh:
        line = fh.readline()
    if not line.startswith("#searched"):
        return {k: set(OLD_GENES) for k in KINDS}
    out = {}
    for part in line.rstrip("\n").split("\t")[1:]:
        kind, _colon, genes = part.partition(":")
        out[kind] = {g for g in genes.split(",") if g}
    return out


def assembly_rows(path):
    """An assembly's search file -> its gene rows (ASSEMBLY_COLUMNS dicts)."""
    with open(path) as fh:
        lines = [line for line in fh if not line.startswith("#")]
    head = lines[0].rstrip("\n").split("\t") if lines else []
    return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in lines[1:] if line.strip()]


MODELS = ("nhmmer", "miniprot")  # the assembly pass' protein-coding gene models: nhmmer's hits alone, or with miniprot's
SAVED = "organelle.fna.gz"  # metadata.ORGANELLE: a genome's organelle molecules its scan kept


def stops_of(row, table):
    """The internal stops of a searched gene's coding sequence as the package will read it (coding_of: miniprot's
    exons joined, nhmmer's single hit whole) by the organelle's genetic code, or None where it ships none (an rRNA,
    nhmmer's hits joined, a gene across the origin, a row without its sequence)."""
    cds = coding_of(row, row["sequence"]) if row.get("sequence") else None
    return gencode.protein_of(cds, table)[1] if cds else None


def single_hit_off(mine, model, table):
    """Whether nhmmer's single hit of a gene (mine) is not the gene's coding sequence, by what miniprot's whole model
    of it (model) shows: the hit covers under NHMMER_WHOLE of its query (a gene with introns whose other exons stay
    below nhmmer's threshold: gene_span reaches for the query's ends into the intron), or, read as the coding
    sequence it would ship as, it holds more internal stops by the organelle's code than miniprot's model (an end
    reached out of frame or beyond the gene's stop, an intron inside the hit)."""
    if mine.get("cover", 1.0) < NHMMER_WHOLE:
        return True
    hit, exons = stops_of(mine, table), stops_of(model, table)
    return hit is not None and exons is not None and hit > exons


def with_miniprot(found, models, table=1):
    """nhmmer's (and barrnap's) genes of an organelle and miniprot's models of its protein-coding genes -> the genes
    taken: miniprot's where nhmmer joined several hits (a gene with introns: its exons unknown) and miniprot aligned the
    whole protein, or where nhmmer found none; nhmmer's single hit kept. tools/organelle_models.py on 100 RefSeq
    records of each organelle (2026-11-07), each searched as an assembly against RefSeq's own exons: with introns,
    miniprot's exons reach a median F1 of 0.96-0.99, nhmmer's joined hits 0.27-0.60; without, both 1.00, nhmmer's ends
    exact more often (COX3 67 of 86 genes, miniprot's 46).
    Since 2026-11-10 (the audit's B2) miniprot's whole model also replaces nhmmer's single hit where that hit is not
    the gene's coding sequence (single_hit_off; table: the organelle's genetic code). Until then every single hit was
    kept: also one exon of a gene with introns, reached to the query's length into the intron and shipped as the
    coding sequence, and a hit of a gene without introns whose end was reached out of frame or past its stop.
    Measured 2026-11-10 on 190 RefSeq organelle genomes (95 fungal, plant and protist mitochondria rich in introns,
    40 animal ones, 55 plastids), each run through run_assembly with its own species' records left out of the
    queries (in brackets: its genus'), 973 protein-coding genes against RefSeq's exons. Of 331 genes with introns 24
    were such a single hit (21 holding 30 to 599 intron bases, 14 without a flag), now 14 (13: miniprot aligns
    them in part, or RefSeq's exons are in doubt). Of 642 genes without introns the coding sequence was RefSeq's
    exactly in 426 (361), now 431 (366); over 6 codons off in 86 (104), now 76 (95); with internal stops RefSeq's
    has not in 59 (87), now 9 (11) - with introns in 52 (54), now 43 (43). 70 (97) rows changed: 22 (25) came
    nearer to RefSeq's, 3 (6) went from within 6 codons to beyond. The share of the query alone (the audit's rule)
    left 51 (74) of the 59 (87) with stops; miniprot's model always lost exact ends (426 to 418; animals 156 to
    141), as tools/organelle_models.py had found."""
    out = dict(found)
    for gene, row in models.items():
        mine = out.get(gene)
        if mine is None or (not row["partial"] and (str(mine.get("exons")) != "1"
                                                   or single_hit_off(mine, row, table))):
            out[gene] = row
    return out


def run_assembly(args):
    """One genome's organelle genes (a pool's task: accession, genome path, out file, animal, queries {organelle:
    Nearest.queries}, barrnap, nhmmer, threads, and since 2026-11-07 proteins {organelle: Nearest.proteins}, codes
    {organelle: genetic code}, miniprot or None) into its out file (ASSEMBLY_COLUMNS; written, even empty, once done)
    -> (accession, genes, why). nhmmer's genes first; barrnap's rRNA where nhmmer found none; with miniprot, its
    models of the protein-coding genes nhmmer found in several hits, or (since 2026-11-10) in a single hit that is not
    the gene's coding sequence (with_miniprot)."""
    acc, path, out, animal, queries, barrnap, nhmmer, threads, *more = args
    proteins, codes, miniprot = more if more else ({}, {}, None)
    try:  # the molecules metadata.describe kept (since 2026-11-08), else the genome read again
        molecules = saved_molecules(path) if os.path.basename(path) == SAVED else organelle_molecules(path)
        rows = []
        for organelle, kingdom in (("mitochondrion", "mito" if animal else "bac"), ("plastid", "bac")):
            if organelle not in queries:  # RefSeq's record taken, or none to hold
                continue
            group = [m for m in molecules if m[1] == organelle]
            found = nhmmer_genes(group, queries.get(organelle, []), nhmmer, threads)
            for (_k, gene), (row, copies) in barrnap_genes(group, barrnap, kingdom, threads).items():
                found.setdefault(gene, dict(row, copies=copies))
            if miniprot:
                found = with_miniprot(found, miniprot_genes(group, proteins.get(organelle, {}), miniprot,
                                                            codes.get(organelle, 1), threads),
                                      codes.get(organelle, 1))
            rows += found.values()
    except (OSError, EOFError, RuntimeError) as e:
        return acc, None, f"{type(e).__name__}: {e}"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    models = "miniprot" if miniprot else "nhmmer"
    inputs = [f"{k}={search_key(k, animal, queries[k], proteins.get(k, {}), codes.get(k, 1), models)}"
              for k in KINDS if k in queries]  # what it was searched with: assembly_pass searches again where other
    with open(util.part_of(out), "w", newline="\n") as fh:
        fh.write("\t".join(["#searched"] + [f"{k}:{','.join(GENES_OF[k])}" for k in KINDS if k in queries]
                           + [f"models:{models}", f"inputs:{','.join(inputs)}"]) + "\n")
        fh.write("\t".join(ASSEMBLY_COLUMNS) + "\n")
        for r in sorted(rows, key=lambda r: (KINDS.index(r["organelle"]), GENES.index(r["gene"]))):
            fh.write("\t".join(str("t" if r.get(c) is True else "f" if r.get(c) is False else r.get(c, "-"))
                               for c in ASSEMBLY_COLUMNS) + "\n")
    os.replace(util.part_of(out), out)
    return acc, len(rows), None


def organelles_of(ncbi, taxid):
    """The organelles a species can hold: no plastid in an animal or a fungus."""
    lineage = set(ncbi.ancestors(ncbi.merged.get(taxid, taxid)))
    return ("mitochondrion",) if lineage & NO_PLASTID else KINDS


def assembly_pass(work_folder, job_map, genome_of, chosen, ncbi, barrnap, nhmmer, workers=1, threads=1, log=print,
                  miniprot=None, every=False, saved_of=None, progress=None, claim=None):
    """The organelle genes of each species' assembly where RefSeq has not every organelle's record it can hold
    (chosen; organelles_of: until 2026-11-06 every animal and fungus was searched for a plastid, ~10,000 genomes), its
    organelle molecules there (saved_of(accession) -> the molecules metadata.describe kept, WORK/done/<accession>/
    organelle.fna.gz, or None: since 2026-11-08, so that a genome is not read again, nor needed on disk; else
    genome_of(accession) -> its genome) and not done yet -> how many were run, failed. Only the organelles RefSeq lacks
    are searched; with miniprot (since 2026-11-07) its protein-coding genes' exons are miniprot's, by the species'
    organelle code (an assembly searched with the other models is searched again; since 2026-11-10 also one searched
    with other inputs: search_key). saved_of and genome_of are asked by threads (util.threaded). claim: genomes.claimed
    of the WORK (picklable: a partial), each assembly searched under its genome's claim (run_claimed; since 2026-11-11:
    a SLURM array's tasks run the pass) -> those run here."""
    nearest = Nearest(chosen, ncbi)
    models = "miniprot" if miniprot else "nhmmer"
    asked = []
    for acc, job in sorted(job_map.items()):
        sp = int(job["species_taxid"])
        need = [k for k in organelles_of(ncbi, sp) if every or (sp, k) not in chosen]  # every: the panel's extras
        if job.get("path") or not need:
            continue
        animal = is_animal(ncbi, sp)
        queries, proteins = {k: nearest.queries(sp, k) for k in need}, {k: nearest.proteins(sp, k) for k in need}
        codes = {k: ncbi.code_of(sp, COMPARTMENT[k]) for k in need}
        # since 2026-11-10 searched again also where it was searched with other inputs (search_key: the queries of
        # another release, another genetic code, other thresholds, a CACHE_VERSION raised)
        same = {f"{k}={search_key(k, animal, queries[k], proteins[k], codes[k], models)}" for k in need}
        asked.append((acc, assembly_file(work_folder, acc), animal, queries, proteins, codes, (need, same, models)))

    def look(item):
        """An assembly's search file and its molecules, or None where it is searched already or not there: asked of
        READ_THREADS assemblies at once (since 2026-11-11: one after another, each a file of WORK/organelle and its
        packed genome's molecules read from the shared file system)."""
        acc, out, *_x, check = item
        if searched_already(out, *check):
            return None
        # asked of the assemblies to search alone (since 2026-11-10: of every assembly; a packed genome's saved
        # molecules are written out for the search, store.py)
        path = (saved_of(acc) if saved_of else None) or genome_of(acc)
        return path if os.path.isfile(path) else None
    todo = [(acc, path, out, animal, queries, barrnap, nhmmer, threads, proteins, codes, miniprot)
            for (acc, out, animal, queries, proteins, codes, _c), path in zip(asked, util.threaded(look, asked))
            if path]
    if not todo:
        return 0, 0
    log(f"organelle genes: {len(todo)} assemblies to search" + (", each under its genome's claim" if claim else ""))
    failed, elsewhere = 0, 0
    # under its genome's claim (since 2026-11-11), where several builds on a WORK (a SLURM array's tasks) share them
    checks = {item[0]: item[-1] for item in asked}
    tasks = [(claim, checks[t[0]], t) for t in todo] if claim else todo
    # progress(what, weights) -> its util.Progress (the build's: every --report_minutes, with the resources used;
    # until 2026-11-08 at each 5% of the assemblies)
    progress = (progress or (lambda what, w: util.Progress(what, w, log, unit="", scale=1)))(
        "assemblies searched", [1] * len(todo))
    with util.process_pool(workers) as pool:  # a worker that dies stops the build, not a wait forever
        for acc, n, why in util.unordered(pool, run_claimed if claim else run_assembly, tasks, progress,
                                          lambda: f", {min(max(1, workers), progress.n - progress.done)} running"):
            if why:
                failed += 1
                log(f"  {acc}: organelle genes failed: {why}")
            elsewhere += n == ELSEWHERE
            progress.step(1, failed=bool(why))
    if elsewhere:
        log(f"  organelle genes: {elsewhere} assemblies searched by another build, or under its claim")
    return len(todo) - elsewhere, failed


ELSEWHERE = "elsewhere"  # run_claimed's genes where another build holds the genome or searched it meanwhile


def searched_already(out, need, same, models):
    """Whether an assembly's search file (out) holds the organelles of need searched with these inputs (same:
    search_key's) and models: searched again for genes added since (CYTB, rbcL, matK, tufA: 2026-11-07), other models,
    or other inputs (since 2026-11-10; a file of before holds no `inputs` and is searched again, once: accepted, as
    nothing says what it was searched with)."""
    if not os.path.isfile(out):
        return False
    got = searched_of(out)
    return all(set(GENES_OF[k]) <= got.get(k, set()) for k in need) and got.get("models", {"nhmmer"}) == {models} \
        and set(same) <= got.get("inputs", set())


def run_claimed(args):
    """run_assembly under its genome's claim (a pool's task: claim, a function of an accession -> a context that
    yields whether this build holds it (genomes.claimed's), searched_already's need, same and models, run_assembly's
    args), since 2026-11-11: several builds on a WORK (a SLURM array's tasks) search each assembly once -> run_assembly's
    (accession, genes, why), genes ELSEWHERE where another build holds the genome or searched it meanwhile."""
    claim, check, item = args
    with claim(item[0]) as held:
        if not held or searched_already(item[2], *check):
            return item[0], ELSEWHERE, None
        return run_assembly(item)


# --- the package ------------------------------------------------------------------------------------------------------

def collect(ncbi, job_map, species_set, chosen, work_folder):
    """The package's organelle genes -> records [(name, sequence, row of GENE_COLUMNS)]: of each species in
    species_set and each organelle, RefSeq's record's genes (chosen), else its assembly's."""
    by_species = {int(j["species_taxid"]): (acc, j) for acc, j in job_map.items() if not j.get("path")}
    out = []
    for sp in sorted(species_set):
        for organelle in KINDS:
            rows = []
            if (sp, organelle) in chosen:
                for r in chosen[(sp, organelle)]:
                    rows.append(dict(r, taxid=int(r["taxid"]), source="refseq", accession="-"))
            elif sp in by_species:
                acc, job = by_species[sp]
                path = assembly_file(work_folder, acc)
                if os.path.isfile(path):
                    rows += [dict(r, taxid=int(job.get("taxid") or sp), source="assembly", accession=acc)
                             for r in assembly_rows(path) if r["organelle"] == organelle]
            animal, plant = is_animal(ncbi, sp), EMBRYOPHYTA in ncbi.ancestors(ncbi.merged.get(sp, sp))
            for r in sorted(rows, key=lambda r: GENES.index(r["gene"])):
                seq = r["sequence"]
                name = f"{r['seqid']}_{r['gene']}"
                model = "refseq" if r["source"] == "refseq" else r.get("model") or "-"
                out.append((name, seq, dict(record=name, gene=r["gene"], organelle=organelle, taxid=r["taxid"],
                                            species_taxid=sp, source=r["source"], seqid=r["seqid"],
                                            accession=r["accession"], start=r["start"], end=r["end"],
                                            strand=r["strand"], length=len(seq), exons=r["exons"], copies=r["copies"],
                                            partial=r["partial"], ambiguous=f"{rrna.ambiguous(seq):.4f}",
                                            name=name_of(organelle, r["gene"], animal, plant), model=model,
                                            parts=r.get("parts") or "-", trans_spliced=r.get("trans_spliced") or "f",
                                            regions=r.get("regions") or "-", cds=coding_of(r, seq),
                                            # not a column of genes.tsv: check reads it (table_of)
                                            transl_table=r.get("transl_table") or "-",
                                            transl_except=r.get("transl_except") or "-")))
    return out


def coding_of(r, seq):
    """A protein-coding gene's coding sequence (its exons joined): RefSeq's spliced, or its record cut by its parts
    (miniprot's exons, one exon's record); None for an rRNA, the spacer, or nhmmer's hits joined (their ends are not
    the exons')."""
    if r["gene"] not in CODING:
        return None
    if r.get("spliced"):
        return r["spliced"]
    parts = parts_read(r.get("parts"))
    if not parts:
        return seq if str(r.get("exons")) == "1" else None
    start, end = int(r["start"]), int(r["end"])
    if start > end:
        return None
    return "".join(seq[x - start:y - start + 1] if s == "+" else seq[end - y:end - x + 1] for x, y, s in parts)


MAX_AMBIGUOUS = 0.005  # a record with this share of other bases than ACGT: low


def check(records, ncbi, vsearch=None, threads=1):
    """The organelle genes' checks across species (since 2026-11-07; qc.check, each gene of each organelle a locus;
    a protein-coding gene's coding sequence compared where it has one, else its record) and their proteins by the
    species' organelle code (since 2026-11-10 a RefSeq CDS's own /transl_table before it: table_of) -> records, their
    rows completed with own_species, same_as, length_check, internal_stops, confidence and flags (`protein` kept in
    the row for write_genes)."""
    from . import qc  # noqa: PLC0415
    seqs, species, lengths = {}, {}, {}
    for name, seq, row in records:
        cds = row.get("cds")
        if row["gene"] in CODING and cds:
            # a /transl_except codon its amino acid, not a stop inside the protein (since 2026-11-10)
            p, stops = with_excepts(*gencode.protein_of(cds, table_of(row, ncbi)), row.get("transl_except"))
            row.update(protein=p, internal_stops=stops)
            lengths[(name, row["organelle"] + ":" + row["gene"])] = len(p)
        key = (name, row["organelle"] + ":" + row["gene"])
        seqs[key], species[name] = (cds or seq), int(row["species_taxid"])
        lengths.setdefault(key, len(cds or seq) // 3)
    checks = qc.check(seqs, species, ncbi, vsearch, threads, lengths)
    for name, _seq, row in records:
        c = checks[(name, row["organelle"] + ":" + row["gene"])]
        row.update(own_species=c["own_species"], same_as=c["same_as"], length_check=c["length"],
                   cross=c["cross"] or "-")
        row["confidence"], row["flags"] = confidence(row)
    return records


def confidence(row):
    """An organelle gene's tier and why (package.tier's four since 2026-11-09): low where a check calls it another
    organism's or not a working gene (own_species f; ambiguous bases over MAX_AMBIGUOUS); medium where something may
    be wrong (partial, its length off its relatives', trans-spliced, nhmmer's hits joined: their exons not exact, a
    sequence shared with another clade's, its owner undecided; a stop inside its protein by the organelle's code, low
    until 2026-11-09: RNA editing (a hornwort's, a fern's plastid genes) or a code NCBI does not give leaves stops in a
    working gene, and an organelle molecule is no NUMT); unverified where nothing is wrong but its own-species check
    could not decide; else high."""
    from .package import tier  # noqa: PLC0415
    low, medium, unverified = [], [], []
    if row.get("own_species") == "f":
        low.append("own_species")
    if float(row.get("ambiguous") or 0) > MAX_AMBIGUOUS:
        low.append("ambiguous")
    if int(row.get("internal_stops") or 0):
        medium.append("internal_stops")
    if row.get("partial") in ("t", True):
        medium.append("partial")
    if row.get("length_check", "ok") != "ok":
        medium.append(f"length_{row['length_check']}")
    if row.get("trans_spliced") in ("t", True):
        medium.append("trans_spliced")
    if row["gene"] in CODING and row.get("model") == "nhmmer" and str(row.get("exons")) != "1":
        medium.append("exons_approximate")
    if row.get("own_species") in ("-", None):
        if row.get("cross") not in (None, "-", ""):
            medium.append("cross_shared")
        else:
            unverified.append("own_species_undecided")
    return tier(low, medium, unverified)


def write_genes(folder, records):
    """records [(name, sequence, row)] -> folder/genes.fna.gz and genes.tsv (GENE_COLUMNS); since 2026-11-07 the
    protein-coding genes' coding sequences (cds.fna.gz) and proteins (proteins.faa.gz, by the organelle's code: check)
    where they have them."""
    os.makedirs(folder, exist_ok=True)
    paths = {n: os.path.join(folder, n) for n in ("genes.fna.gz", "genes.tsv", "cds.fna.gz", "proteins.faa.gz")}
    with util.bgzf_open(util.part_of(paths["genes.fna.gz"])) as fa, \
            open(util.part_of(paths["genes.tsv"]), "w", newline="\n") as tb, \
            util.bgzf_open(util.part_of(paths["cds.fna.gz"])) as fc, \
            util.bgzf_open(util.part_of(paths["proteins.faa.gz"])) as fp:
        tb.write("\t".join(GENE_COLUMNS) + "\n")
        for name, seq, row in records:
            fa.write(f">{name}\n{seq}\n")
            tb.write("\t".join(str("t" if row.get(c) is True else "f" if row.get(c) is False else row.get(c, "-"))
                               for c in GENE_COLUMNS) + "\n")
            if row.get("cds"):
                fc.write(f">{name}\n{row['cds'].upper()}\n")
            if row.get("protein"):
                fp.write(f">{name}\n{row['protein']}\n")
    for path in paths.values():
        os.replace(util.part_of(path), path)
