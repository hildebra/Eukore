"""Eukore's organelle genes (since 2026-11-05; FORMAT.md 3): each species' mitochondrial and plastid rRNA (SSU, LSU, 5S)
and its mitochondrial COX1, where NCBI has them.

- RefSeq's organelle genomes (its release's mitochondrion and plastid GenBank files, REFSEQ_RELEASE): their annotated
  rRNA and COX1 features (rrna_gene, is_cox1), parsed once into WORK/organelle (refseq_genes). A species takes, of
  each organelle, one record (choose): the one with the most of GENES, then its own taxid's (not a strain's or
  subspecies'), then the longest. A gene annotated twice (a plastid's inverted repeat) is taken once, its first,
  `copies` counting them.
- Else the species' own assembly (assembly_pass): the molecules NCBI's FASTA titles name an organelle (TITLE:
  mitochondrion, chloroplast, plastid, apicoplast). nhmmer finds each gene, queried with that gene of the NEAREST
  RefSeq species in NCBI's lineage (its exons joined); barrnap the rRNA nhmmer does not (--kingdom mito in animals'
  mitochondria, whose 12S and 16S are its models; bac elsewhere, as plastids' rRNA are bacterial in kind). barrnap
  alone finds no rRNA in yeast's mitochondrion, nhmmer with Kluyveromyces' all of it (2026-11-05).
A gene's record is its sequence from its first base to its last, introns and all (fungi's and plants' COX1 hold group
I and II introns), on its strand; one across a circular genome's origin is joined across it (its start after its end).
"""

import collections
import gzip
import multiprocessing
import os
import re
import subprocess
import tempfile

from . import rrna, util

REFSEQ_RELEASE = "https://ftp.ncbi.nlm.nih.gov/refseq/release"
RELEASE_FILE = re.compile(r"\b((mitochondrion|plastid)\.\d+\.genomic\.gbff\.gz)\b")
KINDS = ("mitochondrion", "plastid")
GENES = ("SSU", "LSU", "5S", "COX1")
METAZOA = 33208
TITLE = re.compile(r"\b(mitochondrion|chloroplast|plastid|apicoplast)\b", re.IGNORECASE)
MAX_MOLECULE = 15_000_000  # an organelle molecule this long at most (Silene's mitochondrial chromosomes: ~11 Mb)
NEAREST = 3                # RefSeq species whose gene is queried in an assembly's organelle
MAX_EVALUE = 1e-10         # nhmmer's hits taken
MAX_GENE = 40_000          # a gene's hits this close to its best are one gene (introns: a COX1's up to ~20 kb in fungi)
MIN_COVER = 0.8            # of its query's exons a gene covers at least, else partial
BARRNAP = {"mito": {"12S_rRNA": "SSU", "16S_rRNA": "LSU"},
           "bac": {"16S_rRNA": "SSU", "23S_rRNA": "LSU", "5S_rRNA": "5S"}}
GENE_COLUMNS = ("record", "gene", "organelle", "taxid", "species_taxid", "source", "seqid", "accession", "start",
                "end", "strand", "length", "exons", "copies", "partial", "ambiguous")
CACHE_COLUMNS = ("seqid", "organelle", "taxid", "record_length", "gene", "start", "end", "strand", "exons", "copies",
                 "partial", "sequence", "spliced")
ASSEMBLY_COLUMNS = ("gene", "organelle", "seqid", "start", "end", "strand", "exons", "copies", "partial", "sequence")


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
    """A GenBank flat file -> dict(seqid (its VERSION), taxid, organelle (its source's /organelle), length, circular,
    sequence, features [(key, [location's lines], [(qualifier, value)])]) for each record."""
    rec = feat = qual = None
    section = None
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                f = line.split()
                rec = dict(seqid=f[1], taxid=None, organelle="", length=int(f[2]) if f[2].isdigit() else 0,
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
                            elif name == "organelle":
                                rec["organelle"] = value
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
    """A GenBank location -> (parts [(start, end)] 1-based in the order written, strand, partial), or None where it
    is another record's, or on both strands (a trans-spliced gene)."""
    if ":" in location or "^" in location:
        return None
    spans = [(int(a), int(b or a)) for a, b in re.findall(r"(\d+)(?:\.\.[<>]?(\d+))?", location.replace("<", ""))]
    if not spans:
        return None
    n_complement = location.count("complement(")
    if n_complement and not location.startswith("complement(") and n_complement != len(spans):
        return None
    return spans, "-" if n_complement else "+", "<" in location or ">" in location


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


def spliced(seq, spans, strand, location):
    """A feature's parts joined in the order transcribed."""
    order = list(reversed(spans)) if strand == "-" and location.startswith("complement(") else spans
    return "".join(subsequence(seq, a, b, strand) for a, b in order)


NUMBER = re.compile(r"\brrn(\d+(?:\.\d+)?)\b|\b(\d+(?:\.\d+)?)\s?s\b")
SMALL = re.compile(r"small subunit|\bs-rrna\b|\brrns\b|\brns\b|\bssu\b")
LARGE = re.compile(r"large subunit|\bl-rrna\b|\brrnl\b|\brnl\b|\blsu\b")
COX1_PRODUCT = re.compile(r"cytochrome (c )?oxidase subunit (1|i)\b")
COX1_GENES = {"cox1", "coi", "co1", "coxi", "cox-1", "mt-co1", "mtco1", "cox1a"}


def rrna_gene(text, animal_mitochondrion=False):
    """An rRNA feature's /gene and /product, lower case -> SSU, LSU, 5S or None (4.5S, unknown). "16S" is an animal
    mitochondrion's LSU (its SSU is the 12S), elsewhere the SSU."""
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
    if s == 16 and animal_mitochondrion:
        return "LSU"
    return "SSU" if s < 20 else "LSU"


def is_cox1(quals):
    gene = {v.lower() for n, v in quals if n == "gene"}
    product = " ".join(v.lower() for n, v in quals if n == "product")
    return bool(gene & COX1_GENES or COX1_PRODUCT.search(product))


def record_genes(rec, organelle, animal=False):
    """A RefSeq organelle record's genes -> {gene: dict(start, end, strand, exons, copies, partial, sequence, spliced
    (its exons joined))}."""
    found, counts = {}, collections.Counter()
    for key, location, quals in rec["features"]:
        location = "".join(location)
        if key == "rRNA":
            text = " ".join(v.lower() for n, v in quals if n in ("gene", "product"))
            gene = rrna_gene(text, organelle == "mitochondrion" and animal) or \
                rrna_gene(" ".join(v.lower() for n, v in quals if n == "note"), organelle == "mitochondrion" and animal)
        elif key == "CDS" and organelle == "mitochondrion" and is_cox1(quals):
            gene = "COX1"
        else:
            continue
        parsed = parts_of(location) if gene else None
        if not parsed:
            continue
        counts[gene] += 1
        if gene in found:
            continue
        spans, strand, partial = parsed
        if max(e for _s, e in spans) > len(rec["sequence"]):
            continue
        start, end = span_of(spans, rec["length"], rec["circular"])
        found[gene] = dict(start=start, end=end, strand=strand, exons=len(spans), partial=partial,
                           sequence=subsequence(rec["sequence"], start, end, strand),
                           spliced=spliced(rec["sequence"], spans, strand, location))
    for gene, row in found.items():
        row["copies"] = counts[gene]
    return found


def is_animal(ncbi, taxid):
    return taxid is not None and METAZOA in ncbi.ancestors(ncbi.merged.get(taxid, taxid))


def refseq_genes(src_dir, work_folder, ncbi, log=print):
    """Every gene of RefSeq's organelle release (record_genes), parsed once into WORK/organelle/refseq_<release>.tsv.gz
    (CACHE_COLUMNS) -> [row dicts], [] where the release was not fetched."""
    files = release_files(src_dir)
    if not files:
        return []
    key = util.md5(os.path.join(folder_of(src_dir), "release.txt"))[:10]
    cache = os.path.join(work_folder, f"refseq_{key}.tsv.gz")
    if not os.path.isfile(cache):
        log(f"RefSeq organelles: parsing {len(files)} GenBank files")
        os.makedirs(work_folder, exist_ok=True)
        with gzip.open(util.part_of(cache), "wt", compresslevel=6) as fh:
            fh.write("\t".join(CACHE_COLUMNS) + "\n")
            for path in files:
                organelle = "mitochondrion" if os.path.basename(path).startswith("mitochondrion") else "plastid"
                for rec in gbff_records(path):
                    genes = record_genes(rec, organelle, is_animal(ncbi, rec["taxid"]))
                    for gene, g in genes.items():
                        row = (rec["seqid"], organelle, rec["taxid"] or "-", rec["length"], gene, g["start"], g["end"],
                               g["strand"], g["exons"], g["copies"], "t" if g["partial"] else "f", g["sequence"],
                               g["spliced"])
                        fh.write("\t".join(map(str, row)) + "\n")
        os.replace(util.part_of(cache), cache)
    rows = []
    with gzip.open(cache, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["taxid"].isdigit():
                rows.append(r)
    return rows


def choose(rows, ncbi):
    """RefSeq's gene rows -> {(species taxid, organelle): [the rows of its record]}: of each species' records of an
    organelle, the one with the most genes, then its own taxid's, then the longest, then the first by accession."""
    records = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        sp = ncbi.species_of(int(r["taxid"]))
        records[(sp, r["organelle"])][r["seqid"]].append(r)
    out = {}
    for (sp, organelle), by_seqid in records.items():
        best = min(by_seqid.items(), key=lambda kv: (-len(kv[1]), int(kv[1][0]["taxid"]) != sp,
                                                     -int(kv[1][0]["record_length"]), kv[0]))
        out[(sp, organelle)] = best[1]
    return out


# --- a species' own assembly ----------------------------------------------------------------------------------------

def organelle_molecules(path):
    """A genome's organelle molecules (TITLE in their FASTA titles) -> [(name, organelle, sequence)]."""
    out, cur = [], None
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                m = TITLE.search(line)
                cur = None
                if m:
                    cur = (line[1:].split() or [""])[0], \
                        "mitochondrion" if m.group(1).lower() == "mitochondrion" else "plastid", []
                    out.append(cur)
            elif cur is not None:
                cur[2].append(line.strip())
    return [(n, k, "".join(p).upper()) for n, k, p in out if sum(map(len, p)) <= MAX_MOLECULE]


class Nearest:
    """RefSeq's genes (whole, of the chosen records; their exons joined), by organelle, gene and NCBI lineage."""

    def __init__(self, chosen, ncbi):
        self.ncbi, self.seq = ncbi, {}
        self.under = collections.defaultdict(lambda: collections.defaultdict(list))
        for (sp, organelle), rows in sorted(chosen.items()):
            for r in rows:
                if r["partial"] == "f":
                    self.seq[(r["seqid"], r["gene"])] = r["spliced"]
                    for a in ncbi.ancestors(sp):
                        self.under[(organelle, r["gene"])][a].append(r["seqid"])

    def queries(self, taxid, organelle, n=NEAREST):
        """Each gene of the n RefSeq species nearest to taxid with it -> [(gene|seqid, its exons)]."""
        out = []
        for gene in GENES:
            under, got = self.under.get((organelle, gene), {}), []
            for a in self.ncbi.ancestors(self.ncbi.merged.get(taxid, taxid)):
                for s in under.get(a, ()):
                    if len(got) >= n:
                        break
                    if s not in got:
                        got.append(s)
                if len(got) >= n:
                    break
            out += [(f"{gene}|{s}", self.seq[(s, gene)]) for s in got]
        return out


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
            for start, end, strand, partial in found if gene else ():
                hits[(organelle, gene)].append(dict(gene=gene, organelle=organelle, seqid=name, start=start, end=end,
                                                    strand=strand, exons=1, partial=partial,
                                                    sequence=subsequence(seq, start, end, strand)))
    return {k: (min(v, key=lambda h: (h["partial"], h["start"] - h["end"], h["seqid"], h["start"])), len(v))
            for k, v in hits.items()}


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
        same = [h for h in gene_hits if (h["target"], h["query"], h["strand"]) ==
                (best["target"], best["query"], best["strand"])
                and h["start"] <= best["end"] + MAX_GENE and h["end"] >= best["start"] - MAX_GENE]
        q_len = q_lens[best["query"]]
        covered = set()
        for h in same:
            covered.update(range(h["q1"], h["q2"] + 1))
        seq = seqs[best["target"]]
        q_first, q_last = min(h["q1"] for h in same), max(h["q2"] for h in same)
        start, end = min(h["start"] for h in same), max(h["end"] for h in same)
        lead, tail = q_first - 1, q_len - q_last  # the query's ends the hits leave out: the gene reaches them
        if best["strand"] == "+":
            start, end = max(1, start - lead), min(len(seq), end + tail)
        else:
            start, end = max(1, start - tail), min(len(seq), end + lead)
        out[gene] = dict(gene=gene, organelle=organelle, seqid=best["target"], start=start, end=end,
                         strand=best["strand"], exons=len(same), copies=1, partial=len(covered) < MIN_COVER * q_len,
                         sequence=subsequence(seq, start, end, best["strand"]))
    return out


def assembly_file(work_folder, accession):
    return os.path.join(work_folder, "assembly", f"{accession}.tsv")


def run_assembly(args):
    """One genome's organelle genes (a pool's task: accession, genome path, out file, animal, queries {organelle:
    Nearest.queries}, barrnap, nhmmer, threads) into its out file (ASSEMBLY_COLUMNS; written, even empty, once done)
    -> (accession, genes, why). nhmmer's genes first; barrnap's where nhmmer found none."""
    acc, path, out, animal, queries, barrnap, nhmmer, threads = args
    try:
        molecules = organelle_molecules(path)
        rows = []
        for organelle, kingdom in (("mitochondrion", "mito" if animal else "bac"), ("plastid", "bac")):
            group = [m for m in molecules if m[1] == organelle]
            found = nhmmer_genes(group, queries.get(organelle, []), nhmmer, threads)
            for (_k, gene), (row, copies) in barrnap_genes(group, barrnap, kingdom, threads).items():
                found.setdefault(gene, dict(row, copies=copies))
            rows += found.values()
    except (OSError, EOFError, RuntimeError) as e:
        return acc, None, f"{type(e).__name__}: {e}"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(util.part_of(out), "w", newline="\n") as fh:
        fh.write("\t".join(ASSEMBLY_COLUMNS) + "\n")
        for r in sorted(rows, key=lambda r: (KINDS.index(r["organelle"]), GENES.index(r["gene"]))):
            fh.write("\t".join(str("t" if r[c] is True else "f" if r[c] is False else r[c]) for c in ASSEMBLY_COLUMNS)
                     + "\n")
    os.replace(util.part_of(out), out)
    return acc, len(rows), None


def assembly_pass(work_folder, job_map, genome_of, chosen, ncbi, barrnap, nhmmer, workers=1, threads=1, log=print):
    """The organelle genes of each species' assembly where RefSeq has neither of its organelles' records (chosen), its
    genome on disk (genome_of(accession) -> path) and not done yet -> how many were run, failed."""
    nearest = Nearest(chosen, ncbi)
    todo = []
    for acc, job in sorted(job_map.items()):
        sp = int(job["species_taxid"])
        out = assembly_file(work_folder, acc)
        path = genome_of(acc)
        if job.get("path") or all((sp, k) in chosen for k in KINDS) or os.path.isfile(out) or not os.path.isfile(path):
            continue
        todo.append((acc, path, out, is_animal(ncbi, sp), {k: nearest.queries(sp, k) for k in KINDS}, barrnap,
                     nhmmer, threads))
    if not todo:
        return 0, 0
    log(f"organelle genes: {len(todo)} assemblies to search")
    failed = 0
    with multiprocessing.Pool(max(1, workers)) as pool:
        for i, (acc, _n, why) in enumerate(pool.imap_unordered(run_assembly, todo), 1):
            if why:
                failed += 1
                log(f"  {acc}: organelle genes failed: {why}")
            if i % max(1, len(todo) // 20) == 0 or i == len(todo):
                log(f"  {i}/{len(todo)} assemblies searched ({failed} failed)")
    return len(todo), failed


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
                    with open(path) as fh:
                        head = fh.readline().rstrip("\n").split("\t")
                        for line in fh:
                            r = dict(zip(head, line.rstrip("\n").split("\t")))
                            if r["organelle"] == organelle:
                                rows.append(dict(r, taxid=int(job.get("taxid") or sp), source="assembly",
                                                 accession=acc))
            for r in sorted(rows, key=lambda r: GENES.index(r["gene"])):
                seq = r["sequence"]
                name = f"{r['seqid']}_{r['gene']}"
                out.append((name, seq, dict(record=name, gene=r["gene"], organelle=organelle, taxid=r["taxid"],
                                            species_taxid=sp, source=r["source"], seqid=r["seqid"],
                                            accession=r["accession"], start=r["start"], end=r["end"],
                                            strand=r["strand"], length=len(seq), exons=r["exons"], copies=r["copies"],
                                            partial=r["partial"], ambiguous=f"{rrna.ambiguous(seq):.4f}")))
    return out


def write_genes(folder, records):
    """records [(name, sequence, row)] -> folder/genes.fna.gz and genes.tsv (GENE_COLUMNS)."""
    os.makedirs(folder, exist_ok=True)
    fasta, table = os.path.join(folder, "genes.fna.gz"), os.path.join(folder, "genes.tsv")
    with gzip.open(util.part_of(fasta), "wt", compresslevel=6) as fa, \
            open(util.part_of(table), "w", newline="\n") as tb:
        tb.write("\t".join(GENE_COLUMNS) + "\n")
        for name, seq, row in records:
            fa.write(f">{name}\n{seq}\n")
            tb.write("\t".join(str(row.get(c, "-")) for c in GENE_COLUMNS) + "\n")
    os.replace(util.part_of(fasta), fasta)
    os.replace(util.part_of(table), table)
