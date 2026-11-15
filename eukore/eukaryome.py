"""EUKARYOME's taxonomy as one tree, and every taxon of the package placed in it (since 2026-11-14, format 7), so that a
user (protal) names every eukaryote by EUKARYOME's lineage: its species (species.tsv), its rRNA records and organelle
genes, those labelled above species too.

The tree is EUKARYOME's classification table (superkingdom to genus; its version in the manifest), made one tree:
- an empty rank is filled as EUKARYOME fills its own: the nearest name above with the rank's abbreviation
  (`Chordata.subphy.incertae.sedis`, `Pezizomycotina.ord.incertae.sedis`; a superkingdom of none
  `Eukaryota.subdomain.incertae.sedis`, EUKARYOME's own);
- a name the table files under several parents (Ascomycota under Dikarya and Fungi.subreg.incertae.sedis) takes the
  parent of most of its rows: every name under one lineage. A name two kingdoms share (a homonym: Bursaria the plant
  and the ciliate) is two taxa, the one of fewer rows written with its kingdom, as EUKARYOME writes its own
  (`Lankesterella(Viridiplantae)`, which the table's kingdom-specific column gives where it has one); a name used in
  two lineages of one kingdom is written with its parent. A genus the table names as its family (EUKARYOME's
  `Solanaceae;Solanaceae`) is no genus;
- EUKARYOME's long reads' labels (k__ .. s__) add their species, and the families and genera the table lacks
  (`Densosporales.fam08`), under the table's taxa;
- an NCBI taxon of the package is placed by name: `species` (EUKARYOME names the species by NCBI's scientific name),
  `synonym` (by one of NCBI's synonyms), `genus` (EUKARYOME names its genus: the species added under it by its NCBI
  name), `upward` (EUKARYOME names neither, but one of its NCBI clades above: its family, order, class, ..., by
  NCBI's name or the name data/eukaryome_names.tsv gives EUKARYOME's of a deep lineage; the ranks below filled with
  NCBI's names, an empty one with EUKARYOME's placeholder), `none` (none of its clades: under
  `Eukaryota.subdomain.incertae.sedis`). A taxon above species is placed by its own name where EUKARYOME names it
  (`name`), else upward. A name is matched without case, `_` as a space and a bracket after it; of several taxa of a
  name, the one whose kingdom holds the NCBI taxon's group (data/eukaryome_names.tsv's groups: a homonym of another
  kingdom is never taken), then the one sharing the most names with NCBI's lineage. A name the tree holds in another
  lineage is written with its kingdom (or parent) when NCBI's adds it.

A lineage is EUKARYOME's ten ranks, superkingdom;kingdom;subkingdom;phylum;subphylum;class;order;family;genus;species
(the table's `subreg` and `subphy` are the subkingdom and subphylum), filled down to the taxon's rank, the ranks below
it empty; a species as EUKARYOME writes it, `Genus_epithet`.
"""

import collections
import os
import re

from . import labels, util

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
NAMES = os.path.join(DATA, "eukaryome_names.tsv")
RANKS = ("superkingdom", "kingdom", "subkingdom", "phylum", "subphylum", "class", "order", "family", "genus", "species")
GENUS, SPECIES = RANKS.index("genus"), RANKS.index("species")
# the ranks at which a node of a taxon's name is checked against its NCBI lineage (Tree.conflicts; since 2026-11-15)
CHECKED = ((RANKS.index("phylum"), "phylum"), (RANKS.index("class"), "class"))
# the NCBI ranks whose EUKARYOME taxon is learned from their species placed by name (Tree.learn; since 2026-11-15)
LEARNED = tuple((r, RANKS.index(r)) for r in ("genus", "family", "order", "class", "phylum"))
# EUKARYOME's way of telling two taxa of one name apart: the name and its parent's (Chilodontidae.Seguenziida)
DOTTED = re.compile(r"^([A-Z][a-z]+)\.[A-Z][a-z]+$")
# EUKARYOME's abbreviations of its ranks in the names of its placeholders (Chordata.subphy.incertae.sedis)
ABBR = {"superkingdom": "subdomain", "kingdom": "reg", "subkingdom": "subreg", "phylum": "phy", "subphylum": "subphy",
        "class": "cl", "order": "ord", "family": "fam", "genus": "gen"}
LABEL_RANKS = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
DOMAIN = "Eukaryota"
TOP = "Eukaryota.subdomain.incertae.sedis"  # EUKARYOME's superkingdom of the lineages of none
UNNAMED = {"", "unclassified", "unidentified", "unknown", "?", "-", "none"}
MATCHES = ("species", "synonym", "genus", "name", "upward", "none", "label")
TABLE_COLUMNS = ("name", "rank", "parent", "parent_rank", "source", "ncbi_taxid")
TABLE = "taxonomy/eukaryome.tsv"
INCERTAE = ".incertae.sedis"


def clean(name):
    name = (name or "").strip()
    return "" if name.lower() in UNNAMED else name


def placeholder(above, rank):
    """EUKARYOME's name of an empty rank below `above` (its nearest real name)."""
    return f"{above}.{ABBR[rank]}{INCERTAE}"


def is_placeholder(name):
    return INCERTAE in name


def read_names(path=NAMES):
    """data/eukaryome_names.tsv -> ({NCBI name key: EUKARYOME name}, {EUKARYOME name: set of groups})."""
    to_eukaryome, groups = {}, {}
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    for f in rows[1:]:
        if f[2] != "-":
            to_eukaryome[labels.key(f[2])] = f[0]
        if f[3] != "-":
            groups.setdefault(f[0], set()).update(g.strip() for g in f[3].split(","))
    return to_eukaryome, groups


def species_text(name):
    """A species' name as EUKARYOME writes it: Genus_epithet (a genus's bracket left out)."""
    return re.sub(r"\(.*?\)", "", name).strip().replace(" ", "_")


class Tree:
    """EUKARYOME's taxonomy, one tree: a node is (rank index, name), with one parent (None for a superkingdom)."""

    def __init__(self, version="2.1", path=None, names=NAMES):
        self.version, self.path = version, path
        self.md5 = util.md5(path) if path and os.path.isfile(path) else None
        self.parent, self.source, self.ncbi_taxid = {}, {}, {}
        self.by_name = collections.defaultdict(set)  # exact name -> nodes (at any rank)
        self.by_key = collections.defaultdict(set)   # labels.key(name) -> nodes
        self.weight = collections.Counter()          # a node's table rows (which of two homonyms keeps the name)
        self.to_eukaryome, self.groups = read_names(names) if names else ({}, {})
        self._placed, self._labels, self._groups, self._kingdom = {}, {}, {}, {}
        self.matches = collections.Counter()
        self.learned = {}  # NCBI taxid -> the one EUKARYOME node of its rank its species stand under (learn)

    # --- the tree ------------------------------------------------------------------------------------------------

    def add(self, node, parent, source):
        self.parent[node] = parent
        self.source[node] = source
        self.by_name[node[1]].add(node)
        self.by_key[labels.key(node[1])].add(node)
        dotted = DOTTED.match(node[1])
        if dotted:  # EUKARYOME's homonym Chilodontidae.Seguenziida: found by Chilodontidae too (since 2026-11-15)
            self.by_key[labels.key(dotted.group(1))].add(node)
        self._kingdom.pop(node, None)

    def lineage(self, node):
        """The nodes from its superkingdom down to `node`."""
        out = []
        while node is not None:
            out.append(node)
            node = self.parent.get(node)
        return out[::-1]

    def kingdom(self, node):
        if node not in self._kingdom:
            line = self.lineage(node)
            self._kingdom[node] = line[1][1] if len(line) > 1 else None
        return self._kingdom[node]

    def text(self, node):
        """A node's lineage: the ten ranks' names, those below its rank empty ("" for none)."""
        if node is None:
            return ""
        names = [""] * len(RANKS)
        for i, name in self.lineage(node):
            names[i] = name
        return ";".join(names)

    def make(self, i, name, parent, source, taxid=None):
        """The node of `name` at rank i under parent: an existing one, else added; a name the tree holds in another
        lineage written with its kingdom (EUKARYOME's way: Morus(Viridiplantae)), else with its parent."""
        node = (i, name)
        if node in self.parent and self.parent[node] == parent:
            return node
        line = self.lineage(parent) if parent is not None else []
        clash = [n for n in self.by_name.get(name, ()) if n not in line and not (n in self.parent and
                                                                                  self.parent[n] == parent)]
        if clash or node in self.parent:
            mine = line[1][1] if len(line) > 1 else (name if i == 1 else None)
            theirs = {self.kingdom(n) for n in clash} | ({self.kingdom(node)} if node in self.parent else set())
            tags = ([mine] if mine and mine not in theirs else []) + ([parent[1]] if parent else []) + ["2", "3", "4"]
            for tag in tags:
                node = (i, f"{name}({tag})")
                if node not in self.parent or self.parent[node] == parent:
                    break
        if node not in self.parent:
            self.add(node, parent, source)
            if taxid is not None:
                self.ncbi_taxid[node] = taxid
        return node

    @classmethod
    def from_classification(cls, rows, version="2.1", path=None, names=NAMES):
        """EUKARYOME's classification table (rows: its header, then superkingdom .. genus and the kingdom-specific
        genus; labels.read_xlsx_rows) -> Tree (module docstring)."""
        tree = cls(version, path, names)
        header = [str(h).strip().lower() for h in rows[0]] if rows else []
        paths, rows_names = [], []
        for r in rows[1:]:
            r = [str(x or "").strip() for x in r] + [""] * (10 - len(r))
            if any(r[i] and r[i].lower() == header[i] for i in range(min(9, len(header)))):
                continue  # the header again (a row of the table names its columns)
            names = [clean(x) for x in r[:9]]
            if r[9]:  # a homonym's kingdom-specific name (Lankesterella(Viridiplantae))
                names[8] = r[9]
            if names[8] and names[7] and labels.key(names[8]) == labels.key(names[7]):
                names[8] = ""  # a genus named as its family (Solanaceae;Solanaceae): no genus
            if any(names):
                rows_names.append(names)
        # a name's rank: the column most of its rows give it; a row whose names stand a column off (EUKARYOME 2.1's
        # Metazoa;Arthropoda;Insecta;... with Arthropoda as its subkingdom) put back where they stand in most rows,
        # where that keeps their order
        votes = collections.defaultdict(collections.Counter)
        for names in rows_names:
            for i, n in enumerate(names):
                if n:
                    votes[n][i] += 1
        column = {n: max(c.items(), key=lambda kv: (kv[1], -kv[0]))[0] for n, c in votes.items()}
        for names in rows_names:
            at = [(i, n) for i, n in enumerate(names) if n]
            moved = []
            for i, n in at:  # a name of no other row (the row's own genus) goes along with those before it
                moved.append(max(column[n] if votes[n][column[n]] > votes[n][i] else i, moved[-1] + 1 if moved else 0))
            if moved != [i for i, _n in at] and all(a < b for a, b in zip(moved, moved[1:])) and moved[-1] < 9:
                names[:] = [""] * 9
                for i, (_old, n) in zip(moved, at):
                    names[i] = n
            deepest = max((i for i, n in enumerate(names) if n), default=-1)
            path, real = [], DOMAIN
            for i in range(deepest + 1):
                name = names[i] or (TOP if i == 0 else placeholder(real, RANKS[i]))
                path.append(name)
                if not is_placeholder(name) and i:
                    real = name
            paths.append(path)
        # a node of the rows: (rank, name, kingdom) below the kingdom (a name two kingdoms share: two taxa); below the
        # subphylum a name whose rows stand under two real phyla of its kingdom, each borne out by the name's parent
        # (a taxon of that phylum by most of its own rows), is two taxa too (since 2026-11-15: EUKARYOME 2.1 files the
        # golden moles as Mammalia;Alcyonacea;Chrysochloridae, and the one Alcyonacea of most rows, the soft corals',
        # took them into Anthozoa). A row whose phylum alone is another (Coniferophyta;Tracheophytina;Pinopsida, a
        # column off) is no taxon of its own: its names take the parent of most rows, as every name's rows do
        phyla = collections.defaultdict(collections.Counter)
        for path in paths:
            phylum = path[3] if len(path) > 3 and not is_placeholder(path[3]) else None
            for i, n in enumerate(path):
                if i > 3 and phylum:
                    phyla[(i, n, path[1])][phylum] += 1
        borne = collections.defaultdict(set)
        for path in paths:
            phylum = path[3] if len(path) > 3 and not is_placeholder(path[3]) else None
            for i in range(5, len(path)):
                up = phyla.get((i - 1, path[i - 1], path[1]))
                if phylum and not is_placeholder(path[i]) and up and up.most_common(1)[0][0] == phylum:
                    borne[(i, path[i], path[1])].add(phylum)
        split = {k: max(borne[k], key=lambda p: (c[p], p)) for k, c in phyla.items() if len(borne.get(k, ())) > 1}

        def key_of(i, n, kingdom, phylum):
            k = (i, n, kingdom if i >= 2 else None)
            if k in split:
                return k + (phylum if phylum in borne[k] else split[k],)
            return k + (None,)
        parents, weight = collections.defaultdict(collections.Counter), collections.Counter()
        for path in paths:
            kingdom = path[1] if len(path) > 1 else None
            phylum = path[3] if len(path) > 3 and not is_placeholder(path[3]) else None
            keys = [key_of(i, n, kingdom, phylum) for i, n in enumerate(path)]
            for i, k in enumerate(keys):
                weight[k] += 1
                parents[k][keys[i - 1] if i else None] += 1
        parent_of = {k: max(c.items(), key=lambda kv: (kv[1], str(kv[0])))[0] for k, c in parents.items()}

        def chain(k):
            out = []
            while k is not None:
                out.append(k)
                k = parent_of.get(k)
            return out
        name_of = {k: k[1] for k in parent_of}
        by_name = collections.defaultdict(list)
        for k in parent_of:
            by_name[k[1]].append(k)
        used = set(by_name)

        def kingdom_of(k):
            c = chain(k)
            return c[-2][1] if len(c) > 1 else None
        for name, keys in list(by_name.items()):
            if len(keys) < 2:
                continue
            keys.sort(key=lambda k: (-weight[k], str(k)))
            kept = [keys[0]]
            for k in keys[1:]:
                if any(k in chain(x) or x in chain(k) for x in kept):
                    kept.append(k)
                    continue
                king = kingdom_of(k)
                tag = king if king and king != kingdom_of(keys[0]) else (parent_of[k][1] if parent_of[k] else "2")
                base = re.sub(r"\(.*\)$", "", name)  # beside Triadinium(Alveolata): Triadinium(Goniodomataceae)
                new, n = f"{base}({tag})", 2
                while new in used:
                    new, n = f"{base}({tag}){n}", n + 1
                used.add(new)
                name_of[k] = new
        node_of = {k: (k[0], name_of[k]) for k in parent_of}
        for k in sorted(parent_of, key=lambda k: k[0]):
            tree.add(node_of[k], node_of[parent_of[k]] if parent_of[k] else None, "eukaryome")
            tree.weight[node_of[k]] += weight[k]
        return tree

    # --- EUKARYOME's labels -------------------------------------------------------------------------------------

    def kingdoms_apart(self, a, b):
        """Whether two EUKARYOME kingdoms hold none of the same groups (data/eukaryome_names.tsv; one of no groups
        there parts from none)."""
        if not a or not b or a == b:
            return False
        ga, gb = self.groups.get(a), self.groups.get(b)
        return bool(ga and gb and not ga & gb)

    def candidates(self, i, name, kingdom=None):
        """The nodes at rank i of a name (its exact name first, else by key); of a label's kingdom where given: a
        homonym of a kingdom apart from it is not taken, of several those of its own kingdom first."""
        exact = [n for n in self.by_name.get(name, ()) if n[0] == i]
        found = exact or [n for n in self.by_key.get(labels.key(name), ()) if n[0] == i]
        if kingdom and i > 1:
            found = [n for n in found if not self.kingdoms_apart(self.kingdom(n), kingdom)]
            found.sort(key=lambda n: (self.kingdom(n) != kingdom, -self.weight[n], n))
            return found
        return sorted(found, key=lambda n: (-self.weight[n], n))

    def place_label(self, label):
        """A EUKARYOME label (k__Fungi;p__...;g__Genus;s__epithet, unnamed ranks `unclassified`) -> its node (the
        deepest rank it names; the species, family or genus the tree lacks added under the deepest it holds), or None
        for a label naming nothing."""
        if label in self._labels:
            return self._labels[label]
        parts = {LABEL_RANKS[p[0]]: clean(p.split("__", 1)[1]) for p in label.split(";")
                 if "__" in p and p[:1] in LABEL_RANKS}
        if parts.get("species"):
            from .rrna import species_name  # noqa: PLC0415
            sp = species_name(parts.get("genus", ""), parts["species"].replace("_", " "))
            parts["species"] = species_text(sp) if sp and labels.is_species_name(sp) else ""
        named = [RANKS.index(r) for r, n in parts.items() if n]
        if not named:
            self._labels[label] = None
            return None
        deepest = max(named)
        kingdom = parts.get("kingdom") or None
        anchor = None
        for i in range(min(deepest, GENUS), 0, -1):
            name = parts.get(RANKS[i])
            if name:
                found = self.candidates(i, name, kingdom)
                if found:
                    # of several taxa of the name, the one sharing the most names with the label's lineage (since
                    # 2026-11-15: EUKARYOME's Triadinium(Dinoflagellata) is its dinoflagellates', not the ciliates')
                    above = {labels.key(n) for r, n in parts.items() if n and RANKS.index(r) < i}
                    anchor = max(found, key=lambda n: len({labels.key(re.sub(r"\(.*\)$", "", x))
                                                           for _j, x in self.lineage(n)} & above))
                    break
        if anchor is None:
            anchor = self.make(0, TOP, None, "eukaryome")
        node = anchor
        real = next((n for _i, n in reversed(self.lineage(anchor)) if not is_placeholder(n)), DOMAIN)
        for i in range(anchor[0] + 1, deepest + 1):
            name = parts.get(RANKS[i]) or placeholder(real, RANKS[i])
            node = self.make(i, name, node, "eukaryome_label")
            if not is_placeholder(name):
                real = name
        self._labels[label] = node
        return node

    # --- NCBI's taxa --------------------------------------------------------------------------------------------

    def ncbi_group(self, ncbi, taxid):
        """An NCBI taxon's (supergroup, group) (groups.classify)."""
        if taxid not in self._groups:
            from . import groups  # noqa: PLC0415
            self._groups[taxid] = groups.classify(ncbi, taxid)
        return self._groups[taxid]

    def fits(self, node, ncbi, taxid):
        """Whether a node may be an NCBI taxon's: its kingdom's groups (data/eukaryome_names.tsv) hold the taxon's
        group, where the table gives them and NCBI places the taxon in a supergroup (an NCBI lineage of incertae
        sedis, Microheliella's, fits any). Not its superkingdom's: EUKARYOME's superkingdoms hold other lineages than
        the deep lineages of Eukore's groups (Malawimonada under Discoba)."""
        supergroup, group = self.ncbi_group(ncbi, taxid)
        from .groups import INCERTAE  # noqa: PLC0415
        if supergroup == INCERTAE:
            return True
        line = self.lineage(node)
        allowed = self.groups.get(line[1][1]) if len(line) > 1 else None
        return not allowed or group in allowed

    def conflicts(self, node, ncbi, taxid, names):
        """Whether a node of the taxon's name is another lineage's, a homonym or a name EUKARYOME misfiles (since
        2026-11-15): at its phylum or class (CHECKED) it names a taxon not of the taxon's NCBI lineage (names: their
        keys) while EUKARYOME holds the taxon's own NCBI phylum or class at that rank, NCBI's taxon stands not above
        it (EUKARYOME's subkingdom Rhodophyta holds its Cyanidiophyta), a class is in NCBI's phylum but of another
        subphylum than every class of NCBI's name, and below that rank its lineage shares no name with NCBI's. So a
        genus name EUKARYOME holds in another phylum, or a class of another subphylum of NCBI's phylum, alone (a
        mammal's genus EUKARYOME knows as an anthozoan's, an insect's as a barnacle's: MIDORI2's species placed by their
        genus) is not the taxon's, while a genus EUKARYOME classifies otherwise within NCBI's taxon (Neocystis in
        Trebouxiophyceae, NCBI's Chlorophyceae; Galdieria in Cyanidiophyta, NCBI's Rhodophyta's Bangiophyceae) is:
        EUKARYOME's classification is Eukore's."""
        line = self.text(node).split(";")
        ids = ncbi.rank_ids(taxid)
        bare = [labels.key(re.sub(r"\(.*\)$", "", x)) for x in line]
        for i, rank in CHECKED:
            name = re.sub(r"\(.*\)$", "", line[i]) if i < len(line) else ""
            if not name or is_placeholder(name) or labels.key(name) in names:
                continue
            own = ids.get(rank)
            key = labels.key(ncbi.name.get(own, "")) if own else ""
            keys = {key} | ({labels.key(self.to_eukaryome[key])} if key in self.to_eukaryome else set())
            held = [n for k in keys for n in self.by_key.get(k, ()) if n[0] == i]
            if not key or not held:
                continue  # EUKARYOME holds not the taxon's own phylum (class) at that rank: nothing to tell
            if keys & set(bare[:i]):
                continue  # inside NCBI's taxon, which EUKARYOME ranks higher
            if i > CHECKED[0][0] and (bare[CHECKED[0][0]] not in names or
                                      any(self.text(n).split(";")[i - 1] == line[i - 1] for n in held)):
                continue  # a class outside NCBI's phylum, or of its class' subphylum: EUKARYOME's own classification
            below = {labels.key(re.sub(r"\(.*\)$", "", x)) for x in line[i + 1:GENUS] if x and not is_placeholder(x)}
            if not below & names:
                return True
        return False

    def best(self, nodes, ncbi, taxid, names):
        """Of a name's nodes, those that may be the taxon's (fits; since 2026-11-15 not another lineage's:
        conflicts), the one sharing the most names with its NCBI lineage (names: their keys), then of the most
        rows."""
        ok = [n for n in nodes if self.fits(n, ncbi, taxid) and not self.conflicts(n, ncbi, taxid, names)]
        if not ok:
            return None
        return max(ok, key=lambda n: (len({labels.key(x) for _i, x in self.lineage(n)} & names), self.weight[n],
                                      n[1]))

    def ncbi_names(self, ncbi, taxid):
        """The keys of the names of an NCBI taxon's lineage, with EUKARYOME's names of its deep clades."""
        out = set()
        for t in ncbi.ancestors(taxid):
            k = labels.key(ncbi.name.get(t, ""))
            out.add(k)
            if k in self.to_eukaryome:
                out.add(labels.key(self.to_eukaryome[k]))
        return out

    def place(self, ncbi, taxid):
        """An NCBI taxon -> (its node, how it was found: MATCHES): a species or one below it as its species; a taxon
        above species to its own rank (one of no rank to the deepest the tree places it at)."""
        taxid = ncbi.merged.get(int(taxid), int(taxid))
        if taxid in self._placed:
            return self._placed[taxid]
        rank = ncbi.rank.get(taxid)
        sp = ncbi.species_of(taxid)
        at_species = ncbi.rank.get(sp) == "species"
        target = sp if at_species else taxid
        names = self.ncbi_names(ncbi, target)
        node, how = self._by_name(ncbi, target, sp, at_species, names)
        if node is None:
            node, how = self._upward(ncbi, target, names, at_species, rank)
        self._placed[taxid] = (node, how)
        self.matches[how] += 1
        return node, how

    def _by_name(self, ncbi, target, sp, at_species, names, add=True):
        """A taxon placed by its own name (place): a species by its name or a synonym, else by its genus' (the
        species added under it, unless not `add`: learn's, the genus node then); a genus by its name -> (node, how) or
        (None, None)."""
        node, how = None, None
        if at_species:
            sci = ncbi.name.get(sp, "")
            for how_, candidates in (("species", [sci]), ("synonym", ncbi.synonyms.get(sp, []))):
                for name in candidates:
                    found = self.best([n for n in self.by_key.get(labels.key(species_text(name)), ())
                                       if n[0] == SPECIES], ncbi, sp, names)
                    if found:
                        node, how = found, how_
                        break
                if node:
                    break
        ids = ncbi.rank_ids(target)
        if node is None:
            genus = ids.get("genus")
            gname = ncbi.name.get(genus) if genus else ((ncbi.name.get(sp, "").split() or [""])[0]
                                                        if at_species else "")
            if gname and (at_species or genus == target):
                found = self.best([n for n in self.by_key.get(labels.key(gname), ()) if n[0] == GENUS], ncbi,
                                  target, names)
                if found:
                    node, how = found, ("genus" if at_species else "name")
                    if at_species and add:
                        node = self.make(SPECIES, species_text(ncbi.name.get(sp, str(sp))), node, "ncbi", sp)
        return node, how

    def learn(self, ncbi, taxids, share=0.01):
        """The EUKARYOME taxon of each NCBI genus, family, order, class and phylum of the taxa's species (LEARNED;
        since 2026-11-15): the one taxon of its rank its species placed by their own name (species, synonym, genus)
        stand under, a taxon of under `share` of them a stray (a misfiled or homonym genus) -> self.learned {NCBI
        taxid: node}. NCBI's Copepoda is so EUKARYOME's Hexanauplia. A taxon placed upward then stands under its
        deepest NCBI clade's learned taxon before one EUKARYOME names alike (_upward), so that a clade EUKARYOME names
        otherwise is one taxon, not two (a second class Copepoda beside Hexanauplia). A clade EUKARYOME splits (NCBI's
        Pinopsida, its Pinopsida and Cupressopsida) is not learned: its taxon above the parts would leave the ranks
        between them empty, so it is found by its name as before. Call it with every taxon to be placed, before
        placing them."""
        self.forget_placed()  # since 2026-11-15: what earlier placements added under NCBI's names votes for nothing
        votes = collections.defaultdict(collections.Counter)
        for t in sorted({ncbi.merged.get(int(x), int(x)) for x in taxids}):
            sp = ncbi.species_of(t)
            if ncbi.rank.get(sp) != "species":
                continue
            node, _how = self._by_name(ncbi, sp, sp, True, self.ncbi_names(ncbi, sp), add=False)
            if node is None:
                continue
            at = {i: n for i, n in self.lineage(node)}
            for rank, i in LEARNED:
                clade = ncbi.rank_ids(sp).get(rank)
                n = at.get(i)
                if clade and n and not is_placeholder(n):
                    votes[clade][(i, n)] += 1
        self.learned = {}
        for clade, c in votes.items():
            total = sum(c.values())
            held = [node for node, n in c.items() if n >= share * total]
            if len(held) == 1:
                self.learned[clade] = held[0]
        self._placed = {}  # placed again with what was learned
        return self.learned

    def forget_placed(self):
        """The nodes placements added under NCBI's names (source ncbi: a species, a genus EUKARYOME lacks, NCBI's
        clades below a taxon placed upward) taken out of the tree, but those holding another source's node, and the
        placements made forgotten (since 2026-11-15; the audit's T1). The audit placed the published reads before
        the package learned (learn): a copepod genus EUKARYOME lacks stood under a class `Copepoda` of NCBI's name,
        its species then voted for it beside EUKARYOME's Hexanauplia, and Copepoda was learned as nothing."""
        children = collections.defaultdict(list)
        for node, parent in self.parent.items():
            if parent is not None:
                children[parent].append(node)

        def only_ncbi(node):
            return self.source.get(node) == "ncbi" and all(only_ncbi(c) for c in children[node])
        for node in [n for n in self.parent if only_ncbi(n)]:
            del self.parent[node], self.source[node]
            self.ncbi_taxid.pop(node, None)
            self._kingdom.pop(node, None)
            for index, key in ((self.by_name, node[1]), (self.by_key, labels.key(node[1]))):
                index[key].discard(node)
            dotted = DOTTED.match(node[1])
            if dotted:
                self.by_key[labels.key(dotted.group(1))].discard(node)
        self._placed, self._labels = {}, {}

    def _upward(self, ncbi, target, names, at_species, rank):
        """A taxon whose species and genus EUKARYOME does not name: under the deepest of its NCBI clades EUKARYOME
        names (by NCBI's name or eukaryome_names.tsv's), NCBI's names below it, EUKARYOME's placeholder where NCBI
        names none -> (node, upward, name or none)."""
        ids = ncbi.rank_ids(target)
        genus = ids.get("genus")
        start = genus if (at_species and genus) else target
        found, matched = None, None
        for a in ncbi.ancestors(start):
            if at_species and a in (genus, target):
                continue  # its species and genus: not named (place)
            if ncbi.rank.get(a) in ("superkingdom", "domain") or ncbi.name.get(a) in labels.GENERIC:
                break
            key = labels.key(ncbi.name.get(a, ""))
            keys = [key] + ([labels.key(self.to_eukaryome[key])] if key in self.to_eukaryome else [])
            nodes = [n for k in keys for n in self.by_key.get(k, ()) if n[0] < GENUS or (n[0] == GENUS and a == target)]
            found = self.best(nodes, ncbi, target, names)
            if a in self.learned:  # the EUKARYOME taxon its species stand under (learn), or a taxon of the clade's
                # name below it (NCBI's phylum Microsporidia: EUKARYOME's Rozellomycota, its subphylum Microsporidia)
                learned = self.learned[a]
                if not (found and found[0] > learned[0] and learned in self.lineage(found)):
                    found = learned
            if found:
                matched = a
                break
        how = "upward" if found else "none"
        if found is None:
            found = self.make(0, TOP, None, "eukaryome")
            matched = None
        elif matched == target:
            return found, "name"
        last = SPECIES if at_species else RANKS.index(rank) if rank in RANKS[1:] else found[0]
        above = set(ncbi.ancestors(matched)) if matched is not None else set()
        node = found
        real = next((n for _i, n in reversed(self.lineage(found)) if not is_placeholder(n)), DOMAIN)
        for i in range(found[0] + 1, last + 1):
            if i == SPECIES:
                node = self.make(i, species_text(ncbi.name.get(target, str(target))), node, "ncbi", target)
                break
            t = ids.get(RANKS[i])
            if t is not None and t not in above:
                node = self.make(i, ncbi.name[t], node, "ncbi", t)
                real = node[1]
            else:
                node = self.make(i, placeholder(real, RANKS[i]), node, "ncbi")
        return node, how

    # --- output -------------------------------------------------------------------------------------------------

    def rows(self, nodes):
        """taxonomy/eukaryome.tsv's rows: every node on the lineages of `nodes`, by rank and name; the domain first."""
        keep = set()
        for n in nodes:
            if n is not None:
                keep.update(self.lineage(n))
        out = [dict(name=DOMAIN, rank="domain", parent="-", parent_rank="-", source="eukaryome", ncbi_taxid="2759")]
        for n in sorted(keep, key=lambda n: (n[0], n[1])):
            p = self.parent.get(n)
            out.append(dict(name=n[1], rank=RANKS[n[0]], parent=p[1] if p else DOMAIN,
                            parent_rank=RANKS[p[0]] if p else "domain", source=self.source.get(n, "eukaryome"),
                            ncbi_taxid=self.ncbi_taxid.get(n, "-")))
        return out

    def stamp(self, nodes=()):
        """The manifest's record of it: version, file, MD5, its nodes by source, the NCBI taxa by how placed."""
        by_source = collections.Counter(self.source.values())
        return {"version": self.version, "file": os.path.basename(self.path or ""), "md5": self.md5,
                "nodes": dict(by_source), "placed": dict(self.matches)}


def to_ncbi(path=NAMES):
    """{EUKARYOME name's key: NCBI's name} of data/eukaryome_names.tsv (rrna.label_taxon's step upward)."""
    out = {}
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    for f in rows[1:]:
        if f[2] != "-":
            out[labels.key(f[0])] = f[2]
    return out


def record_node(tree, ncbi, taxid, label=None):
    """A record's node and how it was found: its NCBI species' (place, as species.tsv's) where its taxon is a
    species; else its own label's (place_label: a EUKARYOME record's) where that reaches a rank deeper
    than its NCBI taxon's node (`label`: a read labelled by a genus or family NCBI does not name); else its taxon's."""
    node, how = tree.place(ncbi, taxid)
    if label and node is not None and node[0] < SPECIES:
        own = tree.place_label(label)
        if own is not None and own[0] > node[0]:
            return own, "label"
    return node, how


def load(eukaryome_folder, version="2.1", names=NAMES, log=print):
    """The tree of EUKARYOME's classification table in its source folder (sources.ready's), with its long reads'
    labels placed (their species, and the families and genera the table lacks: the species an NCBI species is found
    by) -> Tree, or None where there is none (or no openpyxl to read it)."""
    if not eukaryome_folder:
        return None
    xlsx = os.path.join(eukaryome_folder, "EUKARYOME2_1_classification-2.xlsx")
    if not os.path.isfile(xlsx):
        log(f"WARNING: no {xlsx}: no EUKARYOME lineages")
        return None
    try:
        rows = labels.read_xlsx_rows(xlsx)
    except ImportError:
        log("WARNING: no openpyxl: no EUKARYOME lineages")
        return None
    tree = Tree.from_classification(rows, version, xlsx, names)
    from . import sources  # noqa: PLC0415
    if os.path.isfile(os.path.join(eukaryome_folder, "General_EUK_longread_v2.1.fasta.gz")):
        for _acc, label, _seq in sources.eukaryome_long(eukaryome_folder):
            tree.place_label(label)
    else:
        log("WARNING: no EUKARYOME long reads: EUKARYOME's taxonomy without its species")
    return tree
