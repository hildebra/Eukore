"""NCBI's genetic codes (since 2026-11-07): every translation table of NCBI's taxdump (gencode.dmp), so that a copy is
translated by its genome's own table: a nuclear gene by the taxon's code, a mitochondrial or plastid gene by its
organelle's (ncbi.Ncbi.code_of). The tables are NCBI's amino acids of the 64 codons in TCAG order (TTT, TTC, TTA, TTG,
TCT, ...); `*` is a stop. Tables 27, 28 and 31 read a stop codon as an amino acid (their stops end a gene by context):
the amino acid is taken here.
"""

import functools

BASES = "TCAG"
# NCBI's taxdump gencode.dmp, 2026-11 (the `cde` column); table 34 (archaeal pyrrolysine) is a prokaryote's
TABLES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # standard
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",   # vertebrate mitochondrial
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # yeast mitochondrial
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # mold, protozoan, coelenterate mito.
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",   # invertebrate mitochondrial
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",   # ciliate, dasycladacean, Hexamita nuclear
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",   # echinoderm, flatworm mitochondrial
    10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # euplotid nuclear
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # bacterial, archaeal, plant plastid
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # alternative yeast nuclear (CTG clade)
    13: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",  # ascidian mitochondrial
    14: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # alternative flatworm mitochondrial
    15: "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Blepharisma macronuclear
    16: "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # chlorophycean mitochondrial
    21: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",  # trematode mitochondrial
    22: "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Scenedesmus obliquus mitochondrial
    23: "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Thraustochytrium mitochondrial
    24: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",  # Rhabdopleuridae mitochondrial
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # SR1 and Gracilibacteria
    26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Pachysolen tannophilus nuclear
    27: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # karyorelict nuclear
    28: "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Condylostoma nuclear
    29: "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Mesodinium nuclear
    30: "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # peritrich nuclear
    31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Blastocrithidia nuclear
    32: "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",  # Balanophoraceae plastid
    33: "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",  # Cephalodiscidae mitochondrial
}
CODONS = tuple(a + b + c for a in BASES for b in BASES for c in BASES)


@functools.lru_cache(maxsize=None)
def code(table=1):
    """{codon: amino acid} of an NCBI table; an unknown table reads as the standard one."""
    return dict(zip(CODONS, TABLES.get(table, TABLES[1])))


def translate(cds, table=1):
    """A coding sequence's protein (X for a codon with another base than ACGT; a trailing partial codon dropped)."""
    c, cds = code(table), cds.upper()
    return "".join(c.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 2, 3))


def protein_of(cds, table=1):
    """A copy's protein as the package writes it: its translation without the stop codon that ends it (an internal
    stop stays `*`) -> (protein, internal stops)."""
    p = translate(cds, table)
    if p.endswith("*"):
        p = p[:-1]
    return p, p.count("*")
