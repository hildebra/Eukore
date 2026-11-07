"""NCBI's taxonomy, Eukore's backbone (since 2026-11-05): each taxid's parent, rank, scientific name and synonyms, and the
taxids merged into others, read from NCBI's taxdump.tar.gz or from a package's taxonomy/ (the subset of every taxon on
its lineages, write_subset), so that a reader of the package has the same taxonomy the package was built on.
"""

import collections
import os
import re
import tarfile

from . import util

# names.dmp classes that name the same taxon (anamorph and teleomorph: a fungus's asexual and sexual names)
SYNONYM_CLASSES = ("synonym", "equivalent name", "basionym", "genbank synonym", "anamorph", "teleomorph",
                   "genbank anamorph", "includes")
# NCBI's placeholders for organisms it cannot name (uncultured eukaryote): one taxid for many unrelated genomes
PLACEHOLDER = re.compile(r"^(uncultured|unidentified|unclassified|environmental)\b|metagenome", re.IGNORECASE)


class Ncbi:
    """NCBI's taxonomy from a taxdump.tar.gz (path) or a package's taxonomy folder (nodes.tsv, names.tsv, merged.tsv)."""

    def __init__(self, source):
        self.parent, self.rank, self.name, self.merged, self.gencode = {}, {}, {}, {}, {}
        self.mgcode, self.pgcode = {}, {}  # the mitochondrial and plastid codes (since 2026-11-07)
        self.synonyms = collections.defaultdict(list)
        if os.path.isdir(source):
            self._read_subset(source)
        else:
            self._read_taxdump(source)

    def _read_taxdump(self, path):
        with tarfile.open(path, "r:gz") as tar:
            for member in tar:
                base = os.path.basename(member.name)
                if base not in ("nodes.dmp", "names.dmp", "merged.dmp"):
                    continue
                for raw in tar.extractfile(member):
                    fields = raw.decode().rstrip("\t|\n").split("\t|\t")
                    taxid = int(fields[0])
                    if base == "merged.dmp":
                        self.merged[taxid] = int(fields[1])
                    elif base == "nodes.dmp":
                        self.parent[taxid], self.rank[taxid] = int(fields[1]), fields[2]
                        for column, codes in ((6, self.gencode), (8, self.mgcode), (13, self.pgcode)):
                            if len(fields) > column and fields[column].strip().isdigit():
                                codes[taxid] = int(fields[column])
                    elif fields[3] == "scientific name":
                        self.name[taxid] = fields[1]
                    elif fields[3] in SYNONYM_CLASSES:
                        self.synonyms[taxid].append(fields[1])

    def _read_subset(self, folder):
        for name, handle in (("nodes.tsv", self._node), ("names.tsv", self._name), ("merged.tsv", self._merge)):
            path = os.path.join(folder, name)
            if os.path.isfile(path):
                with util.open_text(path) as fh:
                    next(fh)
                    for line in fh:
                        handle(line.rstrip("\n").split("\t"))

    def _node(self, f):
        self.parent[int(f[0])], self.rank[int(f[0])] = int(f[1]), f[2]
        for column, codes in ((3, self.gencode), (4, self.mgcode), (5, self.pgcode)):
            if len(f) > column and f[column].isdigit():
                codes[int(f[0])] = int(f[column])

    def _name(self, f):
        if f[2] == "scientific name":
            self.name[int(f[0])] = f[1]
        else:
            self.synonyms[int(f[0])].append(f[1])

    def _merge(self, f):
        self.merged[int(f[0])] = int(f[1])

    def ancestors(self, taxid):
        seen = 0
        while taxid in self.parent and seen < 100:
            yield taxid
            if self.parent[taxid] == taxid:
                return
            taxid, seen = self.parent[taxid], seen + 1

    def species_of(self, taxid):
        """The species a taxid is (a strain's, a merged taxid's); the taxid itself above species."""
        taxid = self.merged.get(taxid, taxid)
        for t in self.ancestors(taxid):
            if self.rank.get(t) == "species":
                return t
        return taxid

    def code_of(self, taxid, compartment="nuclear"):
        """A taxon's genetic code (NCBI's translation table) of its nucleus, mitochondrion ("mito") or plastid
        ("plastid"); where the taxonomy gives none: 1, 1, 11 (a taxdump before NCBI's plastid column gives none)."""
        codes, default = {"nuclear": (self.gencode, 1), "mito": (self.mgcode, 1),
                          "plastid": (self.pgcode, 11)}[compartment]
        return codes.get(self.merged.get(taxid, taxid), default)

    def ranks(self, taxid):
        """{rank: name} of the taxid's lineage."""
        at = {}
        for t in self.ancestors(self.merged.get(taxid, taxid)):
            at.setdefault(self.rank.get(t), self.name.get(t))
        return at

    def rank_ids(self, taxid):
        """{rank: taxid} of the taxid's lineage (names can be homonyms: Morus the plant and the bird)."""
        at = {}
        for t in self.ancestors(self.merged.get(taxid, taxid)):
            at.setdefault(self.rank.get(t), t)
        return at

    def is_placeholder(self, taxid):
        return bool(PLACEHOLDER.search(self.name.get(taxid, "")))

    def write_subset(self, folder, taxids):
        """taxonomy/nodes.tsv, names.tsv and merged.tsv of every taxon on the lineages of `taxids`, with their synonyms,
        and the merged taxids that point into them -> the taxa written."""
        keep = set()
        for taxid in taxids:
            keep.update(self.ancestors(self.merged.get(taxid, taxid)))
        os.makedirs(folder, exist_ok=True)
        rows = {"nodes.tsv": ("taxid\tparent\trank\tgencode\tmgcode\tpgcode\n",
                              [f"{t}\t{self.parent[t]}\t{self.rank.get(t, '')}\t{self.code_of(t)}\t"
                               f"{self.code_of(t, 'mito')}\t{self.code_of(t, 'plastid')}\n" for t in sorted(keep)]),
                "names.tsv": ("taxid\tname\tclass\n",
                              [f"{t}\t{self.name[t]}\tscientific name\n" for t in sorted(keep) if t in self.name] +
                              [f"{t}\t{s}\tsynonym\n" for t in sorted(keep) for s in self.synonyms.get(t, ())]),
                "merged.tsv": ("old\tnew\n", [f"{o}\t{n}\n" for o, n in sorted(self.merged.items()) if n in keep])}
        for name, (header, lines) in rows.items():
            path = os.path.join(folder, name)
            with open(util.part_of(path), "w", newline="\n") as fh:
                fh.write(header + "".join(lines))
            os.replace(util.part_of(path), path)
        return len(keep)
