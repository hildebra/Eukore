#!/usr/bin/env python3
"""Tests of the reference screen (eukore/refscreen.py; since 2026-11-10): a record that is another order's gene under
its name is misnamed and the true owner of its sequence is not; a batch of contaminated records, copies of one
sequence, does not make their source's record another clade's; a record with another order's piece is cut to its own
part; a conserved stretch several orders hold tells no clades apart; another gene's record under a gene's name;
MIDORI2's two header layouts; the lists and the screened FASTA.

    python3 -m unittest discover -s tests -t tests
"""

import argparse
import gzip
import io
import os
import random
import sys
import tarfile
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import labels, ncbi as ncbi_mod, refscreen  # noqa: E402
from eukore.refscreen import Record, Verdict  # noqa: E402
import fixtures as fx  # noqa: E402

# NCBI's taxonomy of the test: two insect orders (four moths, three beetles), a mite order of another class, a fish
NODES = [
    (1, 1, "no rank", "root"), (2759, 1, "superkingdom", "Eukaryota"), (33208, 2759, "kingdom", "Metazoa"),
    (6656, 33208, "phylum", "Arthropoda"), (50557, 6656, "class", "Insecta"),
    (7088, 50557, "order", "Lepidoptera"), (7089, 7088, "family", "Notodontidae"), (7090, 7089, "genus", "Nystalea"),
    (101, 7090, "species", "Nystalea prima"), (102, 7090, "species", "Nystalea secunda"),
    (103, 7090, "species", "Nystalea tertia"), (104, 7090, "species", "Nystalea quarta"),
    (7041, 50557, "order", "Coleoptera"), (7042, 7041, "family", "Carabidae"), (7043, 7042, "genus", "Carabus"),
    (201, 7043, "species", "Carabus primus"), (202, 7043, "species", "Carabus secundus"),
    (203, 7043, "species", "Carabus tertius"),
    (6854, 6656, "class", "Arachnida"), (6946, 6854, "order", "Sarcoptiformes"), (6947, 6946, "genus", "Galumna"),
    (301, 6947, "species", "Galumna prima"), (302, 6947, "species", "Galumna secunda"),
    (303, 6947, "species", "Galumna tertia"),
    (7711, 33208, "phylum", "Chordata"), (7898, 7711, "class", "Actinopteri"), (8000, 7898, "order", "Gadiformes"),
    (8001, 8000, "genus", "Gadus"), (401, 8001, "species", "Gadus primus"),
    (402, 8001, "species", "Gadus secundus"),
]

TOOL = '''#!/usr/bin/env python3
"""%s's stand-in: a query placed on a target along the diagonals their shared 12-mers name, without gaps; its
identity is the share of equal bases over their overlap."""
import sys
args = sys.argv[1:]


def read(path):
    seqs, name = {}, None
    for line in (sys.stdin if path == "-" else open(path)):
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs


def place(q, t, index):
    best = None
    for d in {p - i for i in range(0, len(q) - 11, 4) for p in index.get(q[i:i + 12], ())}:
        a, b = max(0, -d), min(len(q), len(t) - d)
        same = sum(q[i] == t[i + d] for i in range(a, b))
        if b > a and (best is None or same > best[0]):
            best = (same, a, b, d)
    return best


def index_of(seqs):
    out = {}
    for name, t in seqs.items():
        ix = out.setdefault(name, {})
        for p in range(len(t) - 11):
            ix.setdefault(t[p:p + 12], []).append(p)
    return out
'''

MINIMAP2 = TOOL % "minimap2 -c" + '''
targets = read(args[-2])
ix = index_of(targets)
for name, q in read(args[-1]).items():
    for tname, t in targets.items():
        got = place(q, t, ix[tname])
        if got and got[0] / (got[2] - got[1]) >= 0.9:
            same, a, b, d = got
            print("\\t".join(map(str, [name, len(q), a, b, "+", tname, len(t), a + d, b + d, same, b - a, 60])))
'''

VSEARCH = TOOL % "vsearch --usearch_global" + '''
q, db, out = (args[args.index(o) + 1] for o in ("--usearch_global", "--db", "--userout"))
low = 100 * float(args[args.index("--id") + 1])
high = 100 * float(args[args.index("--maxid") + 1]) if "--maxid" in args else 101.0
cover = float(args[args.index("--query_cov") + 1])
targets = read(db)
ix = index_of(targets)
with open(out, "w") as fh:
    for name, s in read(q).items():
        for tname, t in targets.items():
            got = place(s, t, ix[tname])
            if got and (got[2] - got[1]) / len(s) >= cover:
                ident = 100.0 * got[0] / (got[2] - got[1])
                if low <= ident <= high:
                    fh.write(f"{name}\\t{tname}\\t{ident:.1f}\\n")
'''


def taxdump(path):
    with tarfile.open(path, "w:gz") as tar:
        for name, lines in (("nodes.dmp", [f"{t}\t|\t{p}\t|\t{r}\t|\t\t|\t0\t|\t1\t|\t1\t|\t1\t|\t5\t|\t1\t|\t0\t|\t0\t|"
                                           f"\t\t|\n" for t, p, r, _ in NODES]),
                            ("names.dmp", [f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, _, _, n in NODES])):
            data = "".join(lines).encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


def mutate(seq, rate, rng):
    return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in seq)


def world(rng, length=900):
    """Each order's gene: an ancestor, its species 4% from it (8% from one another: relatives); the species'
    records and RefSeq's (trusted) of the first moth and the first beetle."""
    ancestor = {o: "".join(rng.choice("ACGT") for _ in range(length)) for o in ("moth", "beetle", "mite", "fish")}
    species = {"moth": [101, 102, 103, 104], "beetle": [201, 202, 203], "mite": [301, 302, 303], "fish": [401]}
    seqs = {t: mutate(ancestor[o], 0.04, rng) for o, taxa in species.items() for t in taxa}
    records = [Record(f"R{t}", t, "COX1", seqs[t]) for t in sorted(seqs)]
    trusted = [Record("NC_101_COX1", 101, "COX1", mutate(seqs[101], 0.002, rng)),
               Record("NC_201_COX1", 201, "COX1", mutate(seqs[201], 0.002, rng))]
    return seqs, records, trusted


class JudgeTest(unittest.TestCase):
    def setUp(self):
        self.ref = refscreen.References()
        self.moth, self.beetle, self.mite = (1, 2, 10), (1, 2, 11), (1, 3, 12)

    def add(self, species, lineage, trusted=False, gene="COX1", twice=False):
        return self.ref.add(f"r{len(self.ref.name)}", species, gene, lineage, trusted, twice)

    def test_windows(self):
        self.assertEqual(refscreen.window_starts(200), [0])
        self.assertEqual(refscreen.window_starts(658), [0, 150, 300, 358])
        self.assertEqual(refscreen.window_starts(600), [0, 150, 300])
        self.assertEqual(refscreen.parting((1, 2, 10), (1, 2, 11)), 2)
        self.assertEqual(refscreen.parting((1, None, 10), (1, 3, 10)), None)  # a rank one names none of: not compared
        self.assertEqual(refscreen.parting((None, 2, 10), (5, 3, 11)), 1)
        self.assertEqual(refscreen.kept_part([(600, 900)], 900), (0, 600))
        self.assertEqual(refscreen.kept_part([(0, 100), (700, 900)], 900), (100, 700))

    def test_a_clade_holds_a_window_by_its_relatives(self):
        own = [self.add(s, self.beetle) for s in (201, 202, 203)]
        moths = [self.add(s, self.moth) for s in (101, 102, 103)]
        tested = self.ref.tested(201, "COX1", self.beetle)
        self.assertEqual(tested, 2)

        def judge(hits):
            return refscreen.judge(201, "COX1", self.beetle, tested, hits, self.ref)
        # a "beetle" record that is a moth's: a copy of one moth, the other moths its relatives (the second best
        # species'), no beetle near
        v = judge([(100.0, moths[0]), (92.0, moths[1]), (91.0, moths[2])])
        self.assertEqual((v.kind, v.clade, v.identity, v.target, v.held, v.own), ("claim", (2, 10), 100.0, moths[0],
                                                                                 91.0, 0.0))
        # one record alone makes no window another clade's (it may be the misnamed one), nor do copies under the
        # other clade's names (a batch of contaminated specimens): no relatives of theirs, the beetles' at 93%
        v = judge([(100.0, moths[0]), (99.5, moths[1]), (93.0, own[1]), (92.0, own[2])])
        self.assertEqual((v.kind, v.held, v.own), ("claim", 0.0, 93.0))
        # a record of its own clade about as near as the copy: the window does not tell the two apart
        self.assertEqual(judge([(100.0, moths[0]), (92.0, moths[1]), (91.0, moths[2]), (96.0, own[1])]).kind, "own")
        # a trusted record of its own species holding it settles it; a trusted relative counts as two species
        mine = self.add(201, self.beetle, trusted=True)
        self.assertEqual(judge([(100.0, moths[0]), (92.0, moths[1]), (91.0, moths[2]), (99.8, mine)]).kind, "own")
        curated = self.add(104, self.moth, trusted=True)
        self.assertEqual(judge([(100.0, moths[0]), (93.0, curated)]).held, 93.0)
        # its own species' other record in the database vouches for nothing
        twin = self.add(201, self.beetle)
        self.assertEqual(judge([(100.0, moths[0]), (92.0, moths[1]), (91.0, moths[2]), (100.0, twin)]).kind, "claim")
        # no copy of another clade's: its own where its own clade's records are the nearest
        self.assertEqual(judge([(93.0, own[1]), (90.0, moths[0])]).kind, "own")
        self.assertEqual(judge([(84.0, own[1]), (96.0, moths[0])]).kind, "none")
        self.assertEqual(judge([]).kind, "none")

    def test_a_conserved_stretch_tells_no_clades_apart(self):
        [self.add(s, self.beetle) for s in (201, 202)]
        moths = [self.add(s, self.moth) for s in (101, 102, 103)]
        mites = [self.add(s, self.mite) for s in (301, 302, 303)]
        tested = self.ref.tested(201, "COX1", self.beetle)
        hits = [(100.0, moths[0]), (96.0, moths[1]), (95.5, moths[2]), (99.0, mites[0]), (96.0, mites[1]),
                (95.0, mites[2])]
        v = refscreen.judge(201, "COX1", self.beetle, tested, hits, self.ref)
        self.assertEqual((v.kind, v.third), ("claim", True))
        self.assertFalse(refscreen.judge(201, "COX1", self.beetle, tested, hits[:3] + [(99.0, mites[0])],
                                         self.ref).third)  # one record of a third clade is no hold on it
        # a lone order's record is tested at its class: another order of it is not another clade
        lone = (1, 2, 99)
        self.assertEqual(self.ref.tested(900, "COX1", lone), 1)
        v = refscreen.judge(900, "COX1", lone, 1, [(100.0, moths[0]), (92.0, moths[1]), (91.0, moths[2])], self.ref)
        self.assertEqual(v.kind, "own")
        v = refscreen.judge(900, "COX1", lone, 1, [(100.0, mites[0]), (92.0, mites[1]), (91.0, mites[2])], self.ref)
        self.assertEqual((v.kind, v.clade), ("claim", (1, 3)))

    def test_another_gene(self):
        self.add(201, self.beetle)
        cytb = self.add(202, self.beetle, trusted=True, gene="CYTB")
        cox = self.add(202, self.beetle, trusted=True)

        def judge(hits):
            return refscreen.judge(201, "COX1", self.beetle, 2, hits, self.ref)
        # the trusted records tell the gene: another gene's nearer than any of its own gene's by GENE_MARGIN
        v = judge([(95.0, cytb)])
        self.assertEqual((v.kind, v.clade), ("gene", "CYTB"))
        self.assertEqual(judge([(95.0, cytb), (94.0, cox)]).kind, "own")
        # a database's record under the window's gene's name does not outweigh them (it may be misfiled alike)
        self.assertEqual(judge([(95.0, cytb), (99.0, self.add(203, self.beetle))]).kind, "gene")
        # no trusted record near: the database's records of the other gene, of two species, each filed once
        once = [self.add(s, self.beetle, gene="nad1") for s in (202, 203)]
        twice = [self.add(s, self.beetle, gene="atp6", twice=True) for s in (202, 203)]
        self.assertEqual(judge([(96.0, once[0]), (90.0, once[1])]).clade, "nad1")
        self.assertEqual(judge([(96.0, once[0])]).kind, "none")
        self.assertEqual(judge([(96.0, twice[0]), (90.0, twice[1])]).kind, "none")

    def test_calls(self):
        moth = (2, 10)
        own, none = Verdict("own", None, 95.0), Verdict("none")

        def claim(identity=98.0, held=92.0, own=0.0, third=False):
            return Verdict("claim", moth, identity, 7, held, own, third)
        # every window another clade's, its relatives nearer than the record's own clade's: misnamed; a single
        # window of a short record too
        c = refscreen.call_of(600, {s: claim() for s in range(0, 301, 50)})
        self.assertEqual((c["call"], c["action"], c["clade"], c["relatives"], c["own_clade"]),
                         ("misnamed", "drop", moth, 92.0, 80.0))
        self.assertEqual(refscreen.call_of(250, {0: claim()})["call"], "misnamed")
        # a stretch is decided by all its windows: one window's relatives a point ahead by chance decide nothing
        v = {s: claim(held=84.0, own=83.0) for s in range(0, 301, 50)}
        v[150] = claim(held=88.0, own=0.0)
        c = refscreen.call_of(600, v)
        self.assertEqual((c["call"], c["action"], c["windows"]), ("weak", "review", {"undecided": 7}))
        # relatives under the search's limit on both sides tell nothing; a third clade holding it neither
        self.assertEqual(refscreen.call_of(600, {s: claim(held=81.0, own=0.0) for s in range(0, 301, 50)})["call"],
                         "weak")
        self.assertEqual(refscreen.call_of(600, {s: claim(third=True) for s in range(0, 301, 50)})["call"], "weak")
        # its own clade holds the relatives: the copy is the other clade's misnamed record
        self.assertEqual(refscreen.call_of(600, {s: claim(held=0.0, own=93.0) for s in range(0, 301, 50)})["call"],
                         "cleared")
        # its own windows and a certain piece: a chimera, cut to its own part where enough of it stays
        v = {s: own for s in range(0, 601, 50)}
        v.update({s: claim() for s in (500, 550, 600)})
        c = refscreen.call_of(900, v)
        self.assertEqual((c["call"], c["action"], c["keep"]), ("chimera", "trim", (0, 500)))
        v = {0: own, 50: claim(), 100: claim(), 150: claim(), 200: claim()}
        self.assertEqual(refscreen.call_of(500, v)["action"], "drop")  # 50 bases of its own would stay
        # one foreign window under 99%: not certain; at 99% or more it is
        v = {s: own for s in range(0, 601, 50)}
        v[600] = claim(98.0)
        self.assertEqual((refscreen.call_of(900, v)["call"], refscreen.call_of(900, v)["action"]), ("weak", "review"))
        v[600] = claim(99.5)
        self.assertEqual(refscreen.call_of(900, v)["call"], "chimera")
        # windows nothing places, a certain piece of under half of them, none its own: for review
        v = {s: none for s in range(0, 601, 50)}
        v.update({0: claim(), 50: claim()})
        self.assertEqual(refscreen.call_of(900, v)["call"], "weak")
        self.assertEqual(refscreen.call_of(900, {s: own for s in range(0, 601, 50)})["call"], "cleared")
        # another gene's record whole, and a piece of one
        gene = Verdict("gene", "CYTB", 96.0, 3)
        self.assertEqual(refscreen.call_of(600, {s: gene for s in range(0, 301, 50)})["call"], "other_gene")
        v = {s: own for s in range(0, 601, 50)}
        v.update({550: gene, 600: gene})
        c = refscreen.call_of(900, v)
        self.assertEqual((c["call"], c["action"], c["keep"]), ("other_gene", "trim", (0, 550)))


class ScreenTest(unittest.TestCase):
    def run_screen(self, tmp, records, trusted):
        ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(tmp, "taxdump.tar.gz")))
        minimap2, vsearch = fx.tool(tmp, "minimap2", MINIMAP2), fx.tool(tmp, "vsearch", VSEARCH)
        return refscreen.screen([lambda: iter(records)], trusted, ncbi, minimap2, vsearch, tmp=tmp)

    def test_the_screen(self):
        rng = random.Random(11)
        seqs, records, trusted = world(rng)
        records += [
            # a "beetle" that is the first moth's gene (a contaminant); two "mites", copies of the second moth's
            Record("misnamed", 203, "COX1", mutate(seqs[101], 0.004, rng)),
            Record("batch1", 302, "COX1", mutate(seqs[102], 0.003, rng)),
            Record("batch2", 303, "COX1", mutate(seqs[102], 0.003, rng)),
            # a moth's record whose last 300 bases are a beetle's; a beetle's that ends in its own species' CYTB
            Record("chimera", 104, "COX1", mutate(seqs[104][:600], 0.003, rng) + mutate(seqs[202][600:], 0.003, rng)),
            Record("short", 102, "COX1", seqs[102][:100]), Record("unknown", 999999, "COX1", seqs[101]),
        ]
        cytb = "".join(rng.choice("ACGT") for _ in range(900))
        trusted += [Record("NC_201_CYTB", 201, "CYTB", cytb)]
        records.append(Record("fused", 201, "COX1",
                              mutate(seqs[201][:500], 0.003, rng) + mutate(cytb[:400], 0.01, rng)))
        with tempfile.TemporaryDirectory() as tmp:
            calls, counts = self.run_screen(tmp, records, trusted)
            got = {name: (c["call"], c["action"]) for (_g, name), c in calls.items() if c["call"] != "cleared"}
            calls_of = {name: c for (_gene, name), c in calls.items()}
            self.assertEqual(got, {"misnamed": ("misnamed", "drop"), "batch1": ("misnamed", "drop"),
                                   "batch2": ("misnamed", "drop"), "chimera": ("chimera", "trim"),
                                   "fused": ("other_gene", "trim")}, calls)
            # the sequences' true owners, the records the contaminants are copies of, are not called
            self.assertEqual(calls_of["misnamed"]["clade"], "order Lepidoptera")
            self.assertEqual(calls_of["misnamed"]["label_clade"], "order Coleoptera")
            self.assertEqual(calls_of["batch1"]["clade"], "class Insecta")
            self.assertEqual((calls_of["chimera"]["keep_start"], calls_of["chimera"]["keep_end"]), (1, 600))
            # cut at the first window the other gene covers COVER of (bases 450-750: 250 of CYTB)
            self.assertEqual((calls_of["fused"]["keep_start"], calls_of["fused"]["keep_end"], calls_of["fused"]["clade"]),
                             (1, 450, "CYTB"))
            self.assertEqual((counts["short"], counts["no_taxid"], counts["trusted"]), (1, 1, 3))
            self.assertEqual(refscreen.dropped(calls), {("COX1", n) for n in ("misnamed", "batch1", "batch2")})
            # the lists, as KSGP's; the screened FASTA: the dropped records out, the cut ones cut
            refscreen.write_lists(os.path.join(tmp, "lists"), calls)
            rows = {n: fx.read_tsv(os.path.join(tmp, "lists", f"{n}.tsv")) for n in ("drop", "trim", "review", "calls")}
            self.assertEqual(sorted(r["record"] for r in rows["drop"]), ["batch1", "batch2", "misnamed"])
            self.assertEqual([(r["record"], r["start"], r["end"]) for r in rows["trim"]],
                             [("chimera", "1", "600"), ("fused", "1", "450")])
            self.assertIn("order Lepidoptera", rows["drop"][2]["reason"])
            self.assertIn("dropped 3, cut 2", refscreen.summary(calls, counts))
            fasta = os.path.join(tmp, "db.fasta")
            with open(fasta, "w") as fh:
                fh.writelines(f">{r.name}\tlineage\n{r.sequence}\n" for r in records)
            n = refscreen.write_screened(fasta, os.path.join(tmp, "db.screened.fasta.gz"), calls, "COX1")
            with gzip.open(os.path.join(tmp, "db.screened.fasta.gz"), "rt") as fh:
                text = fh.read()
            self.assertEqual(n, len(records) - 3)
            self.assertNotIn(">misnamed", text)
            self.assertIn(">chimera\tlineage screen:1-600\n", text)
            cut = "".join(text.split(">chimera")[1].split(">")[0].split("\n")[1:])
            self.assertEqual(cut, records[[r.name for r in records].index("chimera")].sequence[:600])

    def test_a_batch_of_copies_does_not_take_their_source(self):
        # two mite species whose only records are copies of a moth's gene (a batch of contaminated specimens): the
        # moth's own record, nearer its relatives than any mite's, is not called another clade's
        rng = random.Random(3)
        seqs, records, trusted = world(rng)
        records = [r for r in records if r.taxid not in (302, 303)]
        records += [Record("batch1", 302, "COX1", mutate(seqs[102], 0.003, rng)),
                    Record("batch2", 303, "COX1", mutate(seqs[102], 0.003, rng))]
        with tempfile.TemporaryDirectory() as tmp:
            calls, _counts = self.run_screen(tmp, records, trusted)
            self.assertEqual(calls[("COX1", "R102")]["call"], "cleared")

    def test_a_record_filed_under_two_genes(self):
        # MIDORI2 files one record under lrRNA and srRNA (GB273: 9,189 names): each file's copy is its own record;
        # RefSeq's genes tell which it is, the database's other copy of it is no witness
        rng = random.Random(8)
        ssu, lsu = ("".join(rng.choice("ACGT") for _ in range(700)) for _ in range(2))
        trusted = [Record("NC_201_SSU", 201, "SSU", ssu), Record("NC_201_LSU", 201, "LSU", lsu)]
        twice = mutate(ssu, 0.03, rng)  # the second beetle's SSU
        files = {"SSU": [Record("twice", 202, "SSU", twice), Record("S203", 203, "SSU", mutate(ssu, 0.04, rng))],
                 "LSU": [Record("twice", 202, "LSU", twice), Record("L203", 203, "LSU", mutate(lsu, 0.04, rng))]}
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            minimap2, vsearch = fx.tool(tmp, "minimap2", MINIMAP2), fx.tool(tmp, "vsearch", VSEARCH)
            calls, counts = refscreen.screen([lambda g=g: iter(files[g]) for g in ("SSU", "LSU")], trusted, ncbi,
                                             minimap2, vsearch, tmp=tmp)
            self.assertEqual((counts["screened"], counts["filed_twice"], counts["duplicate_names"]), (4, 1, 0))
            got = calls[("LSU", "twice")]
            self.assertEqual((got["call"], got["action"], got["clade"]), ("other_gene", "drop", "SSU"))
            self.assertEqual(calls[("SSU", "twice")]["call"], "cleared")
            self.assertEqual(refscreen.dropped(calls), {("LSU", "twice")})
            self.assertEqual(set(refscreen.read_calls(self.lists(tmp, calls))), set(calls))

    @staticmethod
    def lists(tmp, calls):
        refscreen.write_lists(os.path.join(tmp, "lists2"), calls)
        return os.path.join(tmp, "lists2", "calls.tsv")

    def test_a_clean_database_has_no_calls(self):
        rng = random.Random(5)
        _seqs, records, trusted = world(rng)
        with tempfile.TemporaryDirectory() as tmp:
            calls, counts = self.run_screen(tmp, records, trusted)
            self.assertEqual({c["call"] for c in calls.values()} - {"cleared"}, set())
            self.assertEqual(counts["screened"], len(records))


class BuildTest(unittest.TestCase):
    def test_the_build_screens_midori2(self):
        # eukore_build.reference_screen: the screen kept in WORK/reference_screen (a rerun reads it); the record it
        # drops, a fish species' only one (a moth's gene), names no species in MIDORI2's lineages
        import eukore_build  # noqa: PLC0415
        rng = random.Random(21)
        seqs, records, _trusted = world(rng)
        records.append(Record("FISH.1.<1.>900", 402, "COX1", mutate(seqs[101], 0.004, rng)))
        name = {t: n for t, _p, _r, n in NODES}
        with tempfile.TemporaryDirectory() as tmp:
            src, work = os.path.join(tmp, "sources"), os.path.join(tmp, "work")
            os.makedirs(src)
            ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(src, "taxdump.tar.gz")))
            path = os.path.join(tmp, "MIDORI2_LONGEST_NUC_GB273_CO1_RAW.fasta")
            with open(path, "w") as fh:
                for r in records:
                    genus = next(t for t in ncbi.ancestors(r.taxid) if ncbi.rank[t] == "genus")
                    fh.write(f">{r.name}\troot_1;domain_Eukaryota_2759;genus_{name[genus]}_{genus};"
                             f"species_{name[r.taxid]}_{r.taxid}\n{r.sequence}\n")
            fx.tool(tmp, "minimap2", MINIMAP2)
            fx.tool(tmp, "vsearch", VSEARCH)
            opts = argparse.Namespace(midori2=[path], no_reference_screen=False, minimap2="minimap2", work=work,
                                      threads=1)
            said, before = [], os.environ["PATH"]
            eukore_build.say, say = said.append, eukore_build.say
            os.environ["PATH"] = tmp + os.pathsep + before
            try:
                calls, counts = eukore_build.reference_screen(opts, src, ncbi)
                again, counts_again = eukore_build.reference_screen(opts, src, ncbi)
                self.assertIsNone(eukore_build.reference_screen(
                    argparse.Namespace(**dict(vars(opts), no_reference_screen=True)), src, ncbi))
            finally:
                os.environ["PATH"], eukore_build.say = before, say
            self.assertEqual(refscreen.dropped(calls), {("COX1", "FISH.1.<1.>900")})
            self.assertEqual(refscreen.dropped(again), refscreen.dropped(calls))
            self.assertEqual(counts_again["action_drop"], counts["action_drop"])
            self.assertTrue(any("is there" in line for line in said), said)
            folder = os.path.join(work, "reference_screen", os.listdir(os.path.join(work, "reference_screen"))[0])
            self.assertEqual(sorted(os.listdir(folder)), ["calls.tsv", "counts.json", "drop.tsv", "review.tsv",
                                                          "trim.tsv"])
            rows = refscreen.table_rows("MIDORI2", again)
            self.assertEqual([(r["reference"], r["record"], r["call"], r["clade"]) for r in rows],
                             [("MIDORI2", "FISH.1.<1.>900", "misnamed", "phylum Arthropoda")])
            self.assertEqual(set(refscreen.TABLE_COLUMNS) - set(rows[0]), set())
            log = []
            unscreened = labels.references(src, None, None, log.append, [path])[-1]
            screened = labels.references(src, None, None, log.append, [path],
                                         (refscreen.dropped(again), refscreen.stamp(counts_again)))[-1]
            self.assertEqual((unscreened.label(ncbi, 402)[1], screened.label(ncbi, 402)[1]), ("taxid", "genus"))
            self.assertEqual(screened.label(ncbi, 401)[1], "taxid")
            stamp = screened.stamp()["screen"]
            self.assertEqual((stamp["dropped"], stamp["misnamed"], stamp["version"]), (1, 1, refscreen.VERSION))
            self.assertNotIn("screen", unscreened.stamp())


class Midori2Test(unittest.TestCase):
    OLD = (">MK1.1.1.658###root_1;Eukaryota_2759;Arthropoda_6656;Insecta_50557;Lepidoptera_7088;Notodontidae_7089;"
           "Nystalea_7090;Nystalea prima_101\n")
    NEW = (">MG559732.1.<1.>690\troot_1;cellularroot_cellular_organisms_131567;domain_Eukaryota_2759;"
           "clade_Opisthokonta_33154;kingdom_Metazoa_33208;phylum_Arthropoda_6656;class_Insecta_50557;"
           "order_Lepidoptera_7088;family_Notodontidae_7089;genus_Nystalea_7090;species_Nystalea secunda_102\n")

    def test_both_header_layouts(self):
        self.assertEqual(labels.midori2_header(self.OLD),
                         ("MK1.1.1.658", 101, ["root", "Eukaryota", "Arthropoda", "Insecta", "Lepidoptera",
                                               "Notodontidae", "Nystalea", "Nystalea prima"]))
        self.assertEqual(labels.midori2_header(self.NEW),
                         ("MG559732.1.<1.>690", 102, ["root", "Eukaryota", "Arthropoda", "Insecta", "Lepidoptera",
                                                      "Notodontidae", "Nystalea", "Nystalea secunda"]))
        # a rank the later layout does not give is empty: the eight ranks stay in their places
        names = labels.midori2_header(">X.1\troot_1;domain_Eukaryota_2759;genus_Mysterium_951;"
                                      "species_Mysterium obscurum_952\n")[2]
        self.assertEqual(names, ["root", "Eukaryota", "", "", "", "", "Mysterium", "Mysterium obscurum"])
        self.assertEqual(labels.midori2_gene("MIDORI2_UNIQ_NUC_GB273_CO1_RAW.fasta.gz"), "COX1")
        self.assertEqual(labels.midori2_gene("MIDORI2_LONGEST_NUC_GB273_srRNA_RAW.fasta.gz"), "SSU")
        self.assertEqual(labels.midori2_gene("my_coi.fasta"), "COX1")

    def test_dropped_records_name_no_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            path = os.path.join(tmp, "MIDORI2_LONGEST_NUC_GB273_CO1_RAW.fasta")
            with open(path, "w") as fh:
                fh.write(self.OLD + "ACGT\nACGT\n" + self.NEW + "ACGT\n")
            self.assertEqual([(a, t, s) for a, t, s in labels.midori2_records(path)],
                             [("MK1.1.1.658", 101, "ACGTACGT"), ("MG559732.1.<1.>690", 102, "ACGT")])
            ref = labels.read_midori2(path)
            self.assertEqual((ref.version, ref.label(ncbi, 101)[1], ref.label(ncbi, 102)[1]),
                             ("GB273", "taxid", "taxid"))
            self.assertEqual(ref.label(ncbi, 102)[0][-1], "Nystalea_secunda")
            ref = labels.read_midori2(path, dropped={("COX1", "MG559732.1.<1.>690")})
            self.assertEqual((ref.label(ncbi, 101)[1], ref.label(ncbi, 102)[1]), ("taxid", "genus"))


if __name__ == "__main__":
    unittest.main()
