"""Eukore's own marker loci (since 2026-11-07): classic phylogenetic and barcode genes that are not BUSCOs, RPB1, RPB2,
EF-2, EF-1a and others (data/markers.tsv), as a dataset in BUSCO's layout that odb.Odb reads beside eukaryota_odb12,
so that loci.py calls them in genomes (introns and all) and transcriptomes as it calls the BUSCOs.

Their OrthoDB groups are single-copy in few species (docs/loci.md): OrthoDB groups each with its paralogs (RNA
polymerase II's largest subunit with polymerases I's and III's), and lineages duplicated it. So each marker is defined
by KEGG's KOfam profiles: its ortholog's (targets) and its paralogs' (decoys). The decoys' profiles join the HMMs a
genome's candidates are scored against (odb.Odb.hmm_db): a region whose best profile is a decoy's is not the marker
(loci.confirm's rule, the BUSCOs' own: the locus' HMM must score it best).

The build (SOURCES/eukore_markers): its inputs fetched by the fetch stage (fetch(), fetch_inputs()), the set built of
them by the process stage, which fetches nothing (prepare(), build(); since 2026-11-11: the fetch stage built it too,
its hmmsearch of ~230,000 candidates on the download node):
- KOfam's ko_list (each KO's threshold and score type) and the targets' and decoys' profiles, taken out of its
  profiles archive (~1.6 GB, fetched in ranges at once and deleted; a copy beside the folder is read and kept);
- the candidates: Swiss-Prot's eukaryote proteins (~90 MB from EBI's FTP; UniProt's REST export of the Pfam families
  of reference proteomes was tried and was too slow, 2026-11-07), and NCBI's proteins named as the marker (the
  table's Entrez query `ncbi`, within NCBI_SCOPE: RefSeq's of every eukaryote, GenBank's of the protists; of an
  organism the PER_TAXON of the length nearest the median: its isoforms; ncbi/<locus>.faa.gz): Swiss-Prot holds the
  markers of 14 to 82 species, mostly models; a name only gathers, the profiles decide;
- each candidate scored with every target and decoy profile (hmmsearch); one is the marker's where its best profile
  is one of the marker's targets and one of them reaches its KO's threshold (KOfam's rule: a name's error, a paralog,
  is caught here); a member is of 0.6 to 1.6 times the median length and complete (its profile's coverage, a genome
  call's rule: a truncated model is not), both medians of Swiss-Prot's and RefSeq's records (GenBank's protists' are
  many PCR fragments), one per species: the best checked record (Swiss-Prot's, RefSeq's curated NP_, RefSeq's
  predicted XP_, GenBank's), then one of typical length, then the best scoring; a GenBank member nearest to other
  phyla's alone at near identity is dropped (foreign(): a contaminant of its source, a mislabelled organism);
  members.tsv: every candidate assigned, its source, score, coverage and fate (member, foreign, decoy: scoring a decoy
  of data/marker_decoys.tsv over its target, since 2026-11-09, incomplete, other);
- hmms/<locus>.hmm (the target's profile, named the locus), decoys/<KO>.hmm (named decoy-<KO>), scores_cutoff (half
  the members' 5th percentile score: fragments pass, the decoys tell paralogs), hmm_spans (each target profile's
  median coverage of its members: KOfam's profiles are longer than the proteins, EF-1a's 772 positions for ~460
  residues, so a call's coverage is read over it, odb.Odb.by_locus), refseq_db.faa.gz (<locus>_<taxid>_<accession>:
  miniprot's queries, the nearest species' first), links_to_ODB12.txt (the names), markers.json (the inputs, the
  members per marker).
These files (SET) are written to a folder beside it and put in its place when whole (since 2026-11-10: build(),
swap()); the folder's other files are the build's caches, which a new set takes over.
"""

import collections
import functools
import gzip
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import tarfile
import tempfile
import time
import urllib.request

from . import util

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE = os.path.join(HERE, "data", "markers.tsv")
DECOYS = os.path.join(HERE, "data", "marker_decoys.tsv")  # decoys KEGG has no KO of (decoy_profiles; since 2026-11-09)
DECOY_SCOPE = "eukaryota[orgn]"  # their NCBI proteins': any eukaryote's, environmental samples' too
DECOY_UNLIKE = 0.5   # a decoy's seed shares under this share of its 5-mers (both ways) with every member of the marker
DECOY_MIN = 10       # seeds a decoy profile is built of at least
FOLDER = "eukore_markers"
STATE = "markers.json"
# the set itself in the folder (what odb.Odb, recheck() and the manifest read); the rest of it is the build's caches
SET = ("hmms", "decoys", "scores_cutoff", "links_to_ODB12.txt", "hmm_spans", "members.tsv", "refseq_db.faa.gz", STATE)
OLD = ".old"  # the set a build replaces, for the moment of the swap (swap())
RECHECK = "recheck."  # recheck()'s cache in the folder: recheck.<digest>.tsv.gz
RECIPE = 1  # the build's rules: raised with a change of them that none of recipe()'s constants shows
KOFAM = "https://www.genome.jp/ftp/db/kofam"
SWISSPROT = ("https://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/complete/"
             "uniprot_sprot.fasta.gz")  # its eukaryote proteins: the candidates (UniProt's REST export is too slow)
RELDATE = "uniprot_reldate.txt"  # Swiss-Prot's reldate.txt, fetched beside it (its release, the state's `swissprot`)
WAIT = 30  # seconds between looks at the set another run builds (prepare)
LENGTH_RANGE = (0.6, 1.6)
TYPICAL = 0.1  # of a species' records, those within this share of the median length first
CUTOFF_SHARE = 0.5
CUTOFF_QUANTILE = 0.05
MIN_MEMBERS = 10


def read_table(path=TABLE):
    """data/markers.tsv -> [dict(locus, name, targets, decoys, kind)] in its order."""
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    head = lines[0]
    out = []
    for f in lines[1:]:
        r = dict(zip(head, f))
        out.append(dict(r, targets=r["targets"].split(","), decoys=[k for k in r["decoys"].split(",") if k]))
    return out


def read_decoys(path=DECOYS):
    """data/marker_decoys.tsv -> [dict(locus, decoy, name, ncbi)] ([] without the file)."""
    if not path or not os.path.isfile(path):
        return []
    with open(path) as fh:
        lines = [line.rstrip("\n").split("\t") for line in fh if line.strip() and not line.startswith("#")]
    return [dict(zip(lines[0], f)) for f in lines[1:]]


def folder_of(src_dir):
    return os.path.join(src_dir, FOLDER)


def recipe():
    """A digest of what shapes the set beside its tables and inputs: RECIPE and the constants of the build's rules (the
    state's `recipe`, since 2026-11-10: a set built under another NCBI_SCOPE, PER_TAXON or threshold read as this
    code's)."""
    from . import loci as loci_mod  # noqa: PLC0415
    rules = (RECIPE, NCBI_SCOPE, PER_TAXON, DECOY_SCOPE, DECOY_UNLIKE, DECOY_MIN, LENGTH_RANGE, TYPICAL, CUTOFF_SHARE,
             CUTOFF_QUANTILE, MIN_MEMBERS, FOREIGN_SIM, FOREIGN_QUANTILE, NEAR_SIM, MIN_PAIRS, loci_mod.COMPLETE)
    return hashlib.md5(repr(rules).encode()).hexdigest()[:12]


def md5_or_none(path):
    return util.md5(path) if os.path.isfile(path) else None


def ready(src_dir, table=TABLE, decoys=None):
    """The marker dataset's folder where it is built from this table (and, given decoys, since 2026-11-09, from this
    decoy table: current(), the fetch stage's test; a set built before its decoys still serves the other stages, the
    build saying so), else None. Given decoys, since 2026-11-10 it is also built by this code's rules (recipe()) and
    from the inputs its folder holds now (the state's `inputs`: KOfam's ko_list and profiles, Swiss-Prot's file,
    NCBI's proteins' files, by their MD5s; one that is replaced or deleted makes the set stale, so a newer Swiss-Prot
    put there, or ncbi/ emptied, is built from by the next process stage, after a fetch stage where one is missing):
    a set of before has neither and is built again there, its caches read."""
    folder = folder_of(src_dir)
    try:
        with open(os.path.join(folder, STATE)) as fh:
            state = json.load(fh)
    except (OSError, ValueError):
        return None
    if state.get("table_md5") != util.md5(table):
        return None
    if decoys is not None:
        inputs = state.get("inputs")
        if state.get("decoys_md5") != md5_or_none(decoys) or state.get("recipe") != recipe() or \
                not isinstance(inputs, dict):
            return None
        if any(digest != md5_or_none(os.path.join(folder, name)) for name, digest in inputs.items()):
            return None
    return folder


def current(src_dir):
    """Whether the marker dataset is built from this table and this decoy table (data/marker_decoys.tsv), by this
    code's rules and from the inputs its folder holds (ready())."""
    return ready(src_dir, TABLE, DECOYS) is not None


def ko_list(path):
    """KOfam's ko_list -> {KO: (threshold, score type full or domain, definition)}."""
    out = {}
    with util.open_text(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r.get("threshold", "-") not in ("-", ""):
                out[r["knum"]] = (float(r["threshold"]), r["score_type"], r.get("definition", ""))
    return out


PARTS = 16  # KOfam's server serves one connection at ~0.2 MB/s (2026-11-07): its archive fetched in ranges at once


def fetch_ranges(url, dest, parts=PARTS, log=print):
    """A large file fetched in `parts` byte ranges at once into dest (one stream where the server takes no ranges)."""
    import concurrent.futures  # noqa: PLC0415
    head = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(head, timeout=120) as r:
        size, ranges = int(r.headers.get("Content-Length", 0)), r.headers.get("Accept-Ranges", "") == "bytes"
    if not ranges or not size:
        urllib.request.urlretrieve(url, util.part_of(dest))
        os.replace(util.part_of(dest), dest)
        return dest
    step = (size + parts - 1) // parts
    pieces = [(i, i * step, min(size, (i + 1) * step) - 1) for i in range(parts)]
    log(f"markers: fetching {url} ({size / 1e9:.1f} GB) in {parts} ranges")

    def get(piece):
        i, a, b = piece
        path = f"{dest}.part{i:02d}"
        if os.path.isfile(path) and os.path.getsize(path) == b - a + 1:
            return path
        req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b}"})
        with urllib.request.urlopen(req, timeout=600) as r, open(path, "wb") as fh:
            while True:
                block = r.read(1 << 20)
                if not block:
                    break
                fh.write(block)
        if os.path.getsize(path) != b - a + 1:
            raise OSError(f"{url}: range {a}-{b} cut short")
        return path
    with concurrent.futures.ThreadPoolExecutor(parts) as pool:
        paths = list(pool.map(get, pieces))
    with open(util.part_of(dest), "wb") as out:
        for path in paths:
            with open(path, "rb") as fh:
                while True:
                    block = fh.read(1 << 24)
                    if not block:
                        break
                    out.write(block)
            os.remove(path)
    os.replace(util.part_of(dest), dest)
    return dest


def stream_profiles(url, wanted, folder, log=print, archive=None):
    """KOfam's profiles archive (archive: a local copy; else fetched in ranges beside folder and deleted after), the
    wanted KOs' profiles written to folder/<KO>.hmm (the rest not kept)."""
    os.makedirs(folder, exist_ok=True)
    left = {k for k in wanted if not os.path.isfile(os.path.join(folder, f"{k}.hmm"))}
    if not left:
        return
    fetched = archive is None or not os.path.isfile(archive)
    if fetched:
        archive = fetch_ranges(url, os.path.join(os.path.dirname(folder), "profiles.tar.gz"), log=log)
    log(f"markers: taking KOfam's profiles of {len(left)} KOs out of {archive}")
    with tarfile.open(archive, mode="r|gz") as tar:
        for member in tar:
            ko = os.path.basename(member.name)[:-4] if member.name.endswith(".hmm") else None
            if ko in left:
                data = tar.extractfile(member).read()
                path = os.path.join(folder, f"{ko}.hmm")
                with open(util.part_of(path), "wb") as fh:
                    fh.write(data)
                os.replace(util.part_of(path), path)
                left.discard(ko)
                if not left:
                    break
    if fetched:
        os.remove(archive)  # 1.6 GB; fetched again only where the table names a KO not taken yet
    if left:
        raise OSError(f"KOfam's archive holds no profile of {sorted(left)}")


def swissprot(path, eukaryote):
    """Swiss-Prot's FASTA (>sp|accession|ID name OS=... OX=taxid ...) -> its eukaryotes' proteins (eukaryote(taxid))
    [dict(accession, organism_id, sequence)]."""
    out, cur = [], None
    with util.open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                cur = None
                f = line[1:].split("|")
                ox = line.split(" OX=", 1)[1].split()[0] if " OX=" in line else ""
                if len(f) >= 2 and ox.isdigit() and eukaryote(int(ox)):
                    cur = dict(accession=f[1], organism_id=ox, sequence=[])
                    out.append(cur)
            elif cur is not None:
                cur["sequence"].append(line.strip())
    for r in out:
        r["sequence"] = "".join(r["sequence"])
    return out


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_SCOPE = ("(refseq[filter] OR (eukaryota[orgn] NOT (fungi[orgn] OR metazoa[orgn] OR viridiplantae[orgn]))) "
              "AND eukaryota[orgn] NOT \"environmental samples\"[orgn]")
# RefSeq's of every eukaryote, GenBank's of the protists, none of an environmental sample (a row's query sets its
# length range: [slen]). Entrez's NOT joins two terms, it negates none alone: until 2026-11-09 the scope's
# "refseq[filter] OR NOT (fungi[orgn] OR ...)" took GenBank's fungi, animals and plants too (TEF1's query: 32,708
# records, every one of GenBank's 17,851 fungal; now 3,109, of fungi RefSeq's 104)
SUMMARIES = 5000  # esummary's per call (XML: JSON's is 500)
BATCH = 300  # efetch's ids per call
PER_TAXON = 3  # an organism's proteins fetched (its isoforms and paralogs named alike: the length nearest the median)


def eutils(endpoint, **params):
    """One NCBI E-utilities call (POST) -> its text; tried again on a network error (a connection cut mid-answer
    too); 0.4 s apart (NCBI's 3 a second)."""
    import http.client  # noqa: PLC0415
    import time  # noqa: PLC0415
    import urllib.parse  # noqa: PLC0415
    data = urllib.parse.urlencode(dict(params, tool="eukore")).encode()
    for attempt in range(5):
        try:
            time.sleep(0.4)
            with urllib.request.urlopen(f"{EUTILS}/{endpoint}", data=data, timeout=300) as r:
                return r.read().decode()
        except (OSError, http.client.HTTPException):
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))
    return ""


def chosen(summaries, per_taxon=PER_TAXON):
    """NCBI's summaries [(accession, taxid, length)] -> the accessions fetched: per organism its per_taxon of the
    length nearest the median of all."""
    if not summaries:
        return []
    median = statistics.median(n for _a, _t, n in summaries)
    by_taxon = collections.defaultdict(list)
    for a, t, n in summaries:
        by_taxon[t].append((abs(n - median), a))
    return sorted(a for got in by_taxon.values() for _d, a in sorted(got)[:per_taxon])


EFETCH_TRIES = 3  # a batch whose reply is not whole is asked for again this often


def fasta_records(text):
    """FASTA text -> [(header without its ">", sequence)]. A record starts at a line's ">" alone (until 2026-11-10 the
    text was split at every ">": a title holding one, "5'->3' exoribonuclease", was cached with an empty sequence)."""
    out = []
    for line in text.splitlines():
        if line.startswith(">"):
            out.append([line[1:], []])
        elif out:
            out[-1][1].append(line.strip())
    return [(head, "".join(seq)) for head, seq in out]


def efetch(ids, length, tries=EFETCH_TRIES):
    """One efetch call's proteins of ids (accessions) -> [(header, sequence)], every id answered: as many records as
    ids, none of a known accession shorter than its summary's length (length: {accession: residues, 0 not known});
    asked again where it is not, then an OSError. Since 2026-11-10: NCBI's reply cut short inside an answer of status
    200 (fewer records, or the last one's sequence cut) was written to the cache as it came, and the cache, once
    renamed, was read as complete for good: a marker's candidates then lacked whole batches of organisms."""
    why = ""
    for _attempt in range(tries):
        records = fasta_records(eutils("efetch.fcgi", db="protein", id=",".join(ids), rettype="fasta",
                                       retmode="text"))
        cut = [head.split()[0] for head, seq in records if head.split() and len(seq) < length.get(head.split()[0], 0)]
        if len(records) >= len(ids) and not cut:
            return records
        why = f"{len(records)} records for {len(ids)} ids" + (f", {cut[0]}'s sequence cut short" if cut else "")
    raise OSError(f"NCBI's efetch answered in part {tries} times ({why}; the first id {ids[0]})")


def ncbi_proteins(term, path, log=print, scope=NCBI_SCOPE):
    """NCBI's proteins of an Entrez query within scope (default NCBI_SCOPE), PER_TAXON per organism (chosen()), kept
    in path (>accession taxid) -> [dict(accession, organism_id, sequence)]. path is written under another name and
    renamed once every batch came whole (efetch()): a download that fails leaves no cache."""
    import xml.etree.ElementTree as ET  # noqa: PLC0415
    if not os.path.isfile(path):
        found = eutils("esearch.fcgi", db="protein", term=f"({term}) AND {scope}", usehistory="y", retmax=0)
        count = int(re.search(r"<Count>(\d+)</Count>", found).group(1))
        env, key = re.search(r"<WebEnv>(\S+)</WebEnv>", found).group(1), re.search(r"<QueryKey>(\d+)", found).group(1)
        summaries = []
        for start in range(0, count, SUMMARIES):
            for doc in ET.fromstring(eutils("esummary.fcgi", db="protein", WebEnv=env, query_key=key,
                                            retstart=start, retmax=SUMMARIES)).iter("DocSum"):
                item = {i.get("Name"): i.text for i in doc.iter("Item")}
                if item.get("AccessionVersion") and item.get("TaxId", "0") != "0":
                    summaries.append((item["AccessionVersion"], item["TaxId"], int(item.get("Length") or 0)))
        wanted = chosen(summaries)
        taxid = {a: t for a, t, _n in summaries}
        length = {a: n for a, _t, n in summaries}
        log(f"markers: fetching {len(wanted):,} of NCBI's {count:,} proteins of {term[:50]}...")
        try:
            with util.bgzf_open(util.part_of(path), level=util.WORK_LEVEL) as out:
                for i in range(0, len(wanted), BATCH):
                    for head, seq in efetch(wanted[i:i + BATCH], length):
                        acc = head.split()[0] if head.split() else ""
                        if acc in taxid:
                            out.write(f">{acc} {taxid[acc]}\n{seq}\n")
        except BaseException:  # no half file left under the writer's name either
            if os.path.isfile(util.part_of(path)):
                os.remove(util.part_of(path))
            raise
        os.replace(util.part_of(path), path)
    return read_proteins(path)


def read_proteins(path):
    """ncbi_proteins' cache -> [dict(accession, organism_id, sequence)]."""
    out = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                acc, taxid = line[1:].split()[:2]
                out.append(dict(accession=acc, organism_id=taxid, sequence=next(fh).strip()))
    return out


def release_id(text):
    """UniProt's release id (2026_04) of its reldate.txt's first line ("UniProt Knowledgebase Release 2026_04 consists
    of:"), the text itself where it names none."""
    found = re.search(r"Release (\S+)", text or "")
    return found.group(1) if found else text


def release_of(path):
    """UniProt's release (release_id of its reldate.txt, fetched into path: fetch_inputs), or "?" (since 2026-11-11 a
    file: the build asked UniProt's site, which the process stage does not)."""
    try:
        with open(path) as fh:
            return release_id(fh.readline().strip()) or "?"
    except OSError:
        return "?"


def score(proteins, hmm_db, hmmsearch="hmmsearch", threads=1, cache=None):
    """Each protein's full-sequence score, best domain's score and HMM coverage (its domains' union over the profile's
    length, loci.hmmsearch's) against each profile -> {protein: {profile: (full, domain, coverage)}}. cache: a file the
    scores are kept in, read where it is there (named by what it holds)."""
    if cache and os.path.isfile(cache):
        out = collections.defaultdict(dict)
        with gzip.open(cache, "rt") as fh:
            for line in fh:
                p, h, a, b, c = line.split("\t")
                out[p][h] = (float(a), float(b), float(c))
        return out
    best, spans, qlen = {}, collections.defaultdict(list), {}
    with tempfile.TemporaryDirectory() as tmp:
        faa, dom = os.path.join(tmp, "c.faa"), os.path.join(tmp, "hits.dom")
        with open(faa, "w") as fh:
            fh.writelines(f">{n}\n{s}\n" for n, s in proteins.items())
        subprocess.run([hmmsearch, "--cpu", str(max(1, threads)), "--noali", "-o", os.devnull, "--domtblout", dom,
                        hmm_db, faa], check=True, capture_output=True)
        with open(dom) as fh:
            for line in fh:
                if not line.startswith("#"):
                    f = line.split()
                    key = (f[0], f[3])
                    qlen[f[3]] = int(f[5])
                    best[key] = (float(f[7]), max(float(f[13]), best.get(key, (0.0, float("-inf")))[1]))
                    spans[key].append((int(f[15]), int(f[16])))
    out = collections.defaultdict(dict)
    for (p, h), (full, domain) in best.items():
        covered = {i for a, b in spans[(p, h)] for i in range(a, b + 1)}
        out[p][h] = (full, domain, len(covered) / max(qlen[h], 1))
    if cache:
        keep_scores(out, cache)
    return out


def keep_scores(scores, cache):
    """score()'s result written to its cache file (under another name, then renamed)."""
    with util.bgzf_open(util.part_of(cache), level=util.WORK_LEVEL) as fh:
        fh.writelines(f"{p}\t{h}\t{v[0]}\t{v[1]}\t{v[2]:.4f}\n" for p, hs in scores.items() for h, v in hs.items())
    os.replace(util.part_of(cache), cache)


def classify(scores, kos, targets):
    """A protein's profile scores {KO: (full, domain, ...)} (every marker's targets and decoys) -> the marker's target
    it belongs to: its best by full score, where that is one of the marker's targets and one of them reaches its KO's
    threshold (by its score type, KOfam's rule; a marker of two targets is their union: a plant's actin scores best
    against vertebrates' actin profile, under its threshold, and reaches the other eukaryotes' one), else None."""
    if not scores:
        return None
    best = max(scores, key=lambda k: scores[k][0])
    if best not in targets:
        return None

    def reaches(k):
        threshold, kind, _d = kos.get(k, (float("inf"), "full", ""))
        return k in scores and (scores[k][0] if kind == "full" else scores[k][1]) >= threshold
    return best if any(reaches(k) for k in targets) else None


UNIPROT = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}")
# UniProt's accession (its own pattern: P04050, A0A0B4J2F0), as swissprot() reads it from Swiss-Prot's headers


def source_rank(accession):
    """How far a protein's record is checked: 0 Swiss-Prot's (a UniProt accession), 1 RefSeq's curated (NP_), 2
    RefSeq's predicted (XP_), 3 another (GenBank's: a submitter's model, often of a transcriptome). Since 2026-11-10 by
    UniProt's accession pattern: until then any accession without a version was Swiss-Prot's, and NCBI's protein
    database names some so (a PDB chain 4CBU_A, PIR's and PRF's records): the best checked rank, the basis of a
    marker's median length and never foreign(), for a record nobody curated as the marker."""
    if UNIPROT.fullmatch(accession):
        return 0
    return 1 if accession.startswith("NP_") else 2 if accession.startswith("XP_") else 3


def members(rows, assigned, own_scores, species_of=lambda taxid: taxid, complete=None, basis=None):
    """A marker's members (rows assigned to its targets) -> one per species (species_of(organism_id)): the best
    checked record (source_rank), then the best scoring; of 0.6 to 1.6 times the median length (of basis, the
    accessions of whole records; None: of all), and complete (complete: the accessions whose profile coverage is a
    complete call's, None: all)."""
    mine = [r for r in rows if r["accession"] in assigned]
    if not mine:
        return []
    median = statistics.median(len(r["sequence"]) for r in mine if basis is None or r["accession"] in basis)
    best = {}

    def rank(r):  # a typical length before the score: Arabidopsis's plastid GAPCp (a transit peptide, 422 residues)
        typical = abs(len(r["sequence"]) - median) <= TYPICAL * median  # scores 1.7 bits over its cytosolic GAPC (338)
        return source_rank(r["accession"]), not typical, -own_scores[r["accession"]]
    for r in mine:
        if not LENGTH_RANGE[0] * median <= len(r["sequence"]) <= LENGTH_RANGE[1] * median:
            continue
        if complete is not None and r["accession"] not in complete:
            continue
        sp = species_of(r["organism_id"])
        if sp not in best or rank(r) < rank(best[sp]):
            best[sp] = r
    return sorted(best.values(), key=lambda r: (int(r["organism_id"]), r["accession"]))


FOREIGN_SIM = 0.9       # a member this alike (5-mers, both ways) to its nearest, of another phylum: another organism's
FOREIGN_QUANTILE = 0.99  # ... and more alike than the marker's members of two phyla are but for this share: in the
#                          same kingdom, that kingdom's pairs (animals' β-tubulins are alike across phyla, fungi's are
#                          not), else all pairs of two kingdoms (H3 is alike across them)
NEAR_SIM = 0.15          # a member of its own phylum this near the best keeps it: the true owner of a contaminant's copy
MIN_PAIRS = 100          # a kingdom's pairs to read its own quantile from


def foreign(proteins, taxa, depth_of, block=500, trusted=()):
    """A marker's members (proteins, their taxa) -> the indices of those another organism's protein: its nearest (by
    shared 5-mers, both ways) of another phylum at FOREIGN_SIM or more and more than the marker's members of two phyla
    (of that kingdom, or of two kingdoms) are but for 1 - FOREIGN_QUANTILE, and no member of its own phylum within
    NEAR_SIM of it, where its phylum has others in the marker when the nearest is of its kingdom (depth_of: qc.depth's,
    0 none .. 5 genus; .phylum and .kingdom, -1 for none): a contaminant of the record's source, a mislabelled organism
    (qc.py's own_species rule, read on proteins). trusted: indices never called foreign (Swiss-Prot's and RefSeq's
    records: RefSeq annotates assemblies NCBI screened; a contaminant's copy of them elsewhere is not theirs)."""
    import numpy as np  # noqa: PLC0415
    from scipy import sparse  # noqa: PLC0415
    from .benchmarks import kmers  # noqa: PLC0415
    n = len(proteins)
    if n < 3:
        return set()
    vocab, rows, cols, sizes = {}, [], [], np.zeros(n)
    for i, p in enumerate(proteins):
        ks = kmers(p)
        sizes[i] = max(len(ks), 1)
        for k in ks:
            rows.append(i)
            cols.append(vocab.setdefault(k, len(vocab)))
    m = sparse.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(n, len(vocab)))
    mt = m.T.tocsc()
    phylum = np.array([depth_of.phylum(t) for t in taxa])  # -1 where NCBI names none
    kingdom = np.array([depth_of.kingdom(t) for t in taxa])
    hist, nearest = collections.defaultdict(lambda: np.zeros(101)), []  # kingdom (-2: two kingdoms) -> sims' counts
    for a in range(0, n, block):
        sim = (m[a:a + block] @ mt).toarray() / np.maximum(sizes[a:a + block, None], sizes[None, :])
        for i in range(sim.shape[0]):
            sim[i, a + i] = -1.0
        ph, kg = phylum[a:a + block, None], kingdom[a:a + block, None]
        apart = (ph != phylum[None, :]) & (ph >= 0) & (phylum[None, :] >= 0)
        same = (kg == kingdom[None, :]) & (kg >= 0)
        for k in set(kingdom[a:a + block]) - {-1}:
            hist[k] += np.bincount((np.clip(sim[apart & same & (kg == k)], 0, 1) * 100).astype(int), minlength=101)
        hist[-2] += np.bincount((np.clip(sim[apart & ~same], 0, 1) * 100).astype(int), minlength=101)
        for i in range(sim.shape[0]):
            j = int(sim[i].argmax())
            nearest.append((sim[i, j], j, np.flatnonzero(sim[i] >= sim[i, j] - NEAR_SIM)))

    def threshold(key):
        h = hist[key] if hist[key].sum() >= MIN_PAIRS else hist[-2]
        cum = np.cumsum(h) / max(h.sum(), 1)
        return max(FOREIGN_SIM, (np.searchsorted(cum, FOREIGN_QUANTILE) + 1) / 100)
    cut = {k: threshold(k) for k in list(hist)}
    own = collections.Counter(phylum[phylum >= 0].tolist())  # members per phylum
    trusted = set(trusted)
    out = set()
    for i, (top, j, near) in enumerate(nearest):
        key = kingdom[i] if kingdom[i] >= 0 and kingdom[i] == kingdom[j] else -2
        if i in trusted or depth_of(taxa[i], taxa[j]) != 0 or top < cut.get(key, cut[-2]):
            continue
        if key != -2 and own[phylum[i]] < 2:  # in its kingdom, a lone phylum's conserved protein (an animal's
            continue                           # β-tubulin, H3) is not told from another's: no evidence
        if all(depth_of(taxa[i], taxa[x]) == 0 for x in near):
            out.add(i)
    return out


SOURCES = ("swissprot", "refseq_curated", "refseq_predicted", "genbank")  # source_rank's names
MEMBER_COLUMNS = ("locus", "accession", "taxid", "species_taxid", "source", "score", "coverage", "length", "status")


def profile_name(locus, i):
    """A marker's i-th target profile's name: the locus, then <locus>.2 ..."""
    return locus if i == 0 else f"{locus}.{i + 1}"


class Depth:
    """qc.depth over NCBI's taxonomy, its lineages kept; .phylum(taxid): its phylum's taxid, -1 for none."""

    def __init__(self, ncbi):
        self.ncbi, self.cache = ncbi, {}

    def __call__(self, a, b):
        """qc.depth's, but 0 only where the two are known to part (qc.apart; since 2026-11-09: two protists of
        sister genera, NCBI ranking neither's phylum or class, were another phylum's to foreign()), else 1 at least."""
        from .qc import apart, depth  # noqa: PLC0415
        d = depth(self.ncbi, int(a), int(b), self.cache)
        return d if d or apart(self.ncbi, int(a), int(b), self.cache) else 1

    def phylum(self, taxid):
        taxid = int(taxid)
        return self.cache.setdefault(taxid, self.ncbi.rank_ids(taxid)).get("phylum", -1)

    def kingdom(self, taxid):
        taxid = int(taxid)
        return self.cache.setdefault(taxid, self.ncbi.rank_ids(taxid)).get("kingdom", -1)


def renamed(profile_text, name):
    """An HMMER3 profile with its NAME set."""
    return "".join(f"NAME  {name}\n" if line.startswith("NAME ") else line
                   for line in profile_text.splitlines(keepends=True))


def salvage(path, folder):
    """The build's caches (all but the set's own, SET) in path moved into folder, where it has none of the name."""
    for name in os.listdir(path):
        if name not in SET and not os.path.lexists(os.path.join(folder, name)):
            os.rename(os.path.join(path, name), os.path.join(folder, name))


def recover(src_dir, log=print):
    """What builds that ended half-way left beside the folder (util.dead_partials: dead runs' alone, another build's
    going on is not touched): a whole new set whose swap was cut between its two renames is put in place, a replaced
    set and a half-built one are removed, the caches a cut swap had moved into them taken back."""
    folder = folder_of(src_dir)
    for path in util.dead_partials(src_dir, FOLDER) + util.dead_partials(src_dir, FOLDER + OLD):  # new sets first
        if not os.path.isdir(path):
            continue
        if not os.path.isdir(folder) and os.path.isfile(os.path.join(path, STATE)):
            os.rename(path, folder)
            log(f"markers: {folder} was missing (a build cut while it put its set in place): {path} taken")
            continue
        os.makedirs(folder, exist_ok=True)
        salvage(path, folder)
        shutil.rmtree(path, ignore_errors=True)


def swap(folder, new):
    """A whole new set (new: SET's files) put in the folder's place: the folder's caches moved into it, the folder
    renamed away and new renamed to it, the old set removed. A reader sees the old set or the new one, whole; between
    the two renames (no file is written there) there is no folder, as there is none before the first build."""
    old = util.part_of(folder + OLD)
    shutil.rmtree(old, ignore_errors=True)
    for name in os.listdir(folder):
        if name not in SET:
            os.rename(os.path.join(folder, name), os.path.join(new, name))
    os.rename(folder, old)
    try:
        os.rename(new, folder)
    except OSError:
        os.rename(old, folder)  # the old set back in its place
        raise
    shutil.rmtree(old, ignore_errors=True)


def build(src_dir, table=TABLE, hmmsearch="hmmsearch", threads=1, log=print, decoys=DECOYS):
    """The marker dataset (module docstring) of the inputs fetched (fetch_inputs; nothing is fetched here) -> its
    folder. decoys: the table of decoys KEGG has no KO of (decoy_profile; since 2026-11-09). Since 2026-11-10 the set is
    written to a folder of its own beside the final one (util.part_of's name) and put in its place when it is whole
    (swap()), the caches staying in the final one until then: until then hmms/ and decoys/ were emptied and written
    anew in the folder itself, under the old markers.json, so that a build that failed half-way, or one going on, left
    a half set ready() took."""
    folder, new = folder_of(src_dir), util.part_of(folder_of(src_dir))
    recover(src_dir, log)
    os.makedirs(folder, exist_ok=True)
    shutil.rmtree(new, ignore_errors=True)
    os.mkdir(new)
    try:
        write_set(src_dir, new, table, hmmsearch, threads, log, decoys)
        swap(folder, new)
    except BaseException:
        if os.path.isdir(folder) and os.path.isdir(new):  # else a swap cut between its renames: recover()'s
            try:
                salvage(new, folder)  # the caches a failed swap had moved
                shutil.rmtree(new, ignore_errors=True)
            except OSError:
                pass  # left to the next build's recover()
        raise
    return folder


def query_tag(term, scope):
    """An Entrez query's tag in its cache's name: an edited query, or another scope, is fetched anew."""
    return hashlib.md5(f"{term} {scope}".encode()).hexdigest()[:8]


def kofam_wanted(markers):
    """The KOs whose profiles a marker table names (targets and decoys), sorted."""
    return sorted({k for m in markers for k in m["targets"] + m["decoys"]})


def marker_proteins(marker):
    """A marker's NCBI proteins' file in the marker folder (its query's: an edited query is fetched anew)."""
    return f"ncbi/{marker['locus']}.{query_tag(marker['ncbi'], NCBI_SCOPE)}.faa.gz"


def decoy_proteins(row):
    """A decoy's (data/marker_decoys.tsv's row) NCBI proteins' file in the marker folder."""
    return f"ncbi/decoy-{row['decoy']}.{query_tag(row['ncbi'], DECOY_SCOPE)}.faa.gz"


def inputs_wanted(table=TABLE, decoys=DECOYS):
    """The files in the marker folder the build reads (fetch_inputs fetches them) -> [name]: KOfam's ko_list and the
    table's KOs' profiles, Swiss-Prot's proteins, NCBI's proteins of each marker's query and of each decoy's."""
    markers = read_table(table)
    loci = {m["locus"] for m in markers}
    return (["ko_list.gz", "uniprot_sprot.fasta.gz"] + [f"kofam/{k}.hmm" for k in kofam_wanted(markers)]
            + [marker_proteins(m) for m in markers if m.get("ncbi")]
            + [decoy_proteins(d) for d in read_decoys(decoys) if d["locus"] in loci])


def inputs_missing(src_dir, table=TABLE, decoys=DECOYS):
    """inputs_wanted's files not in the marker folder -> [name]."""
    folder = folder_of(src_dir)
    return [n for n in inputs_wanted(table, decoys) if not os.path.isfile(os.path.join(folder, n))]


def fetch_inputs(src_dir, table=TABLE, download=None, log=print, decoys=DECOYS):
    """The files the build reads (inputs_wanted) and Swiss-Prot's release date fetched into the marker folder, each
    kept where it is there; nothing is computed (the fetch stage's, since 2026-11-11: build() fetched them)."""
    folder = folder_of(src_dir)
    os.makedirs(os.path.join(folder, "ncbi"), exist_ok=True)
    get = download or (lambda u, d: urllib.request.urlretrieve(u, d))
    reldate = os.path.join(folder, RELDATE)
    if not os.path.isfile(reldate):
        try:
            get(SWISSPROT.rsplit("/", 1)[0] + "/reldate.txt", reldate)
        except OSError as e:  # the state's release then "?"
            log(f"WARNING: markers: Swiss-Prot's release date not fetched ({e})")
    if not inputs_missing(src_dir, table, decoys):
        log(f"markers: the set's inputs are there ({folder})")
        return folder
    markers = read_table(table)
    listing = os.path.join(folder, "ko_list.gz")
    if not os.path.isfile(listing):
        get(f"{KOFAM}/ko_list.gz", listing)
    kos = ko_list(listing)
    wanted = kofam_wanted(markers)
    missing = [k for k in wanted if k not in kos]
    if missing:
        raise ValueError(f"data/markers.tsv: KOs KOfam does not hold: {missing}")
    stream_profiles(f"{KOFAM}/profiles.tar.gz", wanted, os.path.join(folder, "kofam"), log,
                    os.path.join(folder, "profiles.tar.gz"))  # a copy there is read (and kept)
    sprot = os.path.join(folder, "uniprot_sprot.fasta.gz")
    if not os.path.isfile(sprot):
        log(f"markers: fetching Swiss-Prot ({SWISSPROT})")
        get(SWISSPROT, sprot)
    loci = {m["locus"] for m in markers}
    for term, name, scope in ([(m["ncbi"], marker_proteins(m), NCBI_SCOPE) for m in markers if m.get("ncbi")]
                              + [(d["ncbi"], decoy_proteins(d), DECOY_SCOPE) for d in read_decoys(decoys)
                                 if d["locus"] in loci]):
        if not os.path.isfile(os.path.join(folder, name)):
            ncbi_proteins(term, os.path.join(folder, name), log, scope)
    return folder


def write_set(src_dir, new, table=TABLE, hmmsearch="hmmsearch", threads=1, log=print, decoys=DECOYS):
    """build()'s work: the set (SET) written into new of the inputs in the marker folder (an input missing is an
    OSError: nothing is fetched here), the caches into the marker folder, markers.json last."""
    folder = folder_of(src_dir)
    missing = inputs_missing(src_dir, table, decoys)
    if missing:
        raise OSError(f"{len(missing)} of its inputs not fetched into {folder} (e.g. {missing[:3]}): --stage fetch "
                      "fetches them")
    markers = read_table(table)
    decoy_rows = read_decoys(decoys)
    listing = os.path.join(folder, "ko_list.gz")
    kos = ko_list(listing)
    wanted = kofam_wanted(markers)
    missing = [k for k in wanted if k not in kos]
    if missing:
        raise ValueError(f"data/markers.tsv: KOs KOfam does not hold: {missing}")
    sprot = os.path.join(folder, "uniprot_sprot.fasta.gz")
    from .ncbi import Ncbi  # noqa: PLC0415
    ncbi = Ncbi(os.path.join(src_dir, "taxdump.tar.gz"))
    eukaryote = functools.lru_cache(maxsize=None)(lambda t: 2759 in ncbi.ancestors(ncbi.merged.get(t, t)))
    rows = {r["accession"]: r for r in swissprot(sprot, eukaryote)}
    release = release_of(os.path.join(folder, RELDATE))
    from_ncbi = {}
    inputs = ["ko_list.gz", "kofam.hmm", "uniprot_sprot.fasta.gz"]  # the state's `inputs` (ready()), with ncbi/'s
    for m in markers:  # NCBI's proteins gathered by name (Swiss-Prot holds few species): the profiles decide
        if m.get("ncbi"):
            name = marker_proteins(m)
            got = read_proteins(os.path.join(folder, name))
            inputs.append(name)
            from_ncbi[m["locus"]] = len(got)
            for r in got:
                if r["organism_id"].isdigit() and eukaryote(int(r["organism_id"])):
                    rows.setdefault(r["accession"], r)
    db = os.path.join(folder, "kofam.hmm")
    with open(util.part_of(db), "w", newline="\n") as out:
        for k in wanted:
            out.write(open(os.path.join(folder, "kofam", f"{k}.hmm")).read())
    os.replace(util.part_of(db), db)
    log(f"markers: scoring {len(rows):,} candidates with {len(wanted)} profiles")
    with open(db) as fh:
        key = hashlib.md5((fh.read() + "\n".join(sorted(rows))).encode()).hexdigest()[:10]
    scores = score({a: r["sequence"] for a, r in rows.items()}, db, hmmsearch, threads,
                   os.path.join(folder, f"scores.{key}.tsv.gz"))
    for sub in ("hmms", "decoys"):  # this table's alone (since 2026-11-08: a decoy dropped from it kept competing)
        os.makedirs(os.path.join(new, sub))
    from . import loci as loci_mod  # noqa: PLC0415
    species_of = functools.lru_cache(maxsize=None)(lambda t: ncbi.species_of(int(t)))
    depth = Depth(ncbi)
    counts, cut_lines, names, records, spans, table_rows = {}, [], [], [], [], []
    for m in markers:
        targets = m["targets"]
        assigned = {a for a, s in scores.items() if classify(s, kos, set(targets))}
        best = {a: max(targets, key=lambda k: scores[a].get(k, (0.0, 0.0, 0.0))[0]) for a in assigned}
        own = {a: scores[a][best[a]][0] for a in assigned}
        # what a whole protein is, from whole records' (Swiss-Prot's, RefSeq's): GenBank's protists are many PCR
        # fragments (EF-1a's: 351 residues on average, 0.44 of its profile, against RefSeq's 460 and 0.58)
        basis = {a for a in assigned if source_rank(a) <= 2}
        span = {}
        for k in set(best.values()):
            got = [scores[a][k][2] for a in basis if best[a] == k]
            span[k] = statistics.median(got if len(got) >= MIN_MEMBERS else
                                        [scores[a][k][2] for a in assigned if best[a] == k])
        complete = {a for a in assigned if scores[a][best[a]][2] >= loci_mod.COMPLETE * span[best[a]]}
        kept = members(list(rows.values()), assigned, own, species_of, complete,
                       basis if len(basis) >= MIN_MEMBERS else None)
        bad = foreign([r["sequence"] for r in kept], [int(r["organism_id"]) for r in kept], depth,
                      trusted=[i for i, r in enumerate(kept) if source_rank(r["accession"]) <= 2])
        dropped = {kept[i]["accession"] for i in bad}
        kept = [r for r in kept if r["accession"] not in dropped]
        # the decoys KEGG has no KO of (decoy_profile; since 2026-11-09: EFL for EF-1a), built before the cutoff: a
        # member that scores one over its target is that paralog's (EFL filed as EF-1a), not the marker's
        built, paralogs = {}, set()
        for d in [d for d in decoy_rows if d["locus"] == m["locus"]]:
            built[d["decoy"]], over = decoy_profile(d, m, kept, sorted(own[r["accession"]] for r in kept), folder, db,
                                                    hmmsearch, threads, log, new)
            inputs.append(decoy_proteins(d))
            paralogs |= over
        kept = [r for r in kept if r["accession"] not in paralogs]
        if len(kept) < MIN_MEMBERS:
            raise ValueError(f"marker {m['locus']}: {len(kept)} members (fewer than {MIN_MEMBERS})")
        is_member = {r["accession"] for r in kept}
        for a in sorted(assigned):
            r = rows[a]
            table_rows.append((m["locus"], a, r["organism_id"], species_of(r["organism_id"]), SOURCES[source_rank(a)],
                               f"{own[a]:.1f}", f"{scores[a][best[a]][2]:.3f}", len(r["sequence"]),
                               "member" if a in is_member else "foreign" if a in dropped else
                               "decoy" if a in paralogs else "incomplete" if a not in complete else "other"))
        own_scores = sorted(own[r["accession"]] for r in kept)
        cutoff = CUTOFF_SHARE * own_scores[int(CUTOFF_QUANTILE * (len(own_scores) - 1))]
        hmm = os.path.join(new, "hmms", f"{m['locus']}.hmm")
        with open(hmm, "w", newline="\n") as fh:
            for i, k in enumerate(targets):  # a second target profile (actin's) named <locus>.2 ...
                fh.write(renamed(open(os.path.join(folder, "kofam", f"{k}.hmm")).read(), profile_name(m["locus"], i)))
        for i, k in enumerate(targets):  # each target profile's median coverage of the members it scores best
            got = [scores[r["accession"]][k][2] for r in kept if best[r["accession"]] == k]
            if got:
                spans.append(f"{profile_name(m['locus'], i)}\t{statistics.median(got):.3f}\n")
        for k in m["decoys"]:
            with open(os.path.join(new, "decoys", f"{k}.hmm"), "w", newline="\n") as fh:
                fh.write(renamed(open(os.path.join(folder, "kofam", f"{k}.hmm")).read(), f"decoy-{k}"))
        cut_lines.append(f"{m['locus']}\t{cutoff:.1f}\n")
        names.append(f"{m['locus']}\t{m['name']}\n")
        records += [(f"{m['locus']}_{r['organism_id']}_{r['accession']}", r["sequence"]) for r in kept]
        sources = collections.Counter(SOURCES[source_rank(r["accession"])] for r in kept)
        counts[m["locus"]] = dict(members=len(kept), assigned=len(assigned), complete=len(complete),
                                  foreign=sorted(dropped), sources=dict(sources), cutoff=round(cutoff, 1),
                                  targets=m["targets"], decoys=m["decoys"], kind=m["kind"],
                                  built_decoys=built, decoy_members=sorted(paralogs))
        log(f"markers: {m['locus']}: {len(kept)} species' proteins (of {len(assigned)} assigned, {len(complete)} "
            f"complete, {len(dropped)} foreign, {len(paralogs)} a decoy's; {dict(sources)}), cutoff {cutoff:.0f}")
    for name, lines in (("scores_cutoff", cut_lines), ("links_to_ODB12.txt", names), ("hmm_spans", spans)):
        with open(os.path.join(new, name), "w", newline="\n") as fh:
            fh.writelines(lines)
    with open(os.path.join(new, "members.tsv"), "w", newline="\n") as fh:  # every candidate assigned, its fate
        fh.write("\t".join(MEMBER_COLUMNS) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in table_rows)
    with util.bgzf_open(os.path.join(new, "refseq_db.faa.gz"), level=util.WORK_LEVEL) as fh:
        fh.writelines(f">{n}\n{s}\n" for n, s in records)
    # recipe, inputs (since 2026-11-10): the rules and the files it is built of, which ready() holds it against
    state = dict(table_md5=util.md5(table), kofam_ko_list_md5=util.md5(listing), swissprot=release,
                 decoys_md5=util.md5(decoys) if decoys and os.path.isfile(decoys) else None,
                 ncbi=dict(scope=NCBI_SCOPE, per_taxon=PER_TAXON, fetched=from_ncbi), candidates=len(rows),
                 markers=counts, recipe=recipe(),
                 inputs={name: md5_or_none(os.path.join(folder, name)) for name in sorted(set(inputs))},
                 built=__import__("datetime").datetime.now().isoformat(timespec="seconds"))
    with open(os.path.join(new, STATE), "w", newline="\n") as fh:  # the last file: a folder holding it is a whole set
        json.dump(state, fh, indent=1)
    return new


def sibling(binary, name):
    """A HMMER tool beside binary (hmmsearch's: its hmmalign, hmmbuild), else name on the PATH."""
    path = shutil.which(binary) if binary else None
    near = os.path.join(os.path.dirname(path), name) if path else None
    return near if near and os.path.isfile(near) else shutil.which(name) or name


def decoy_profile(row, marker, members_, member_scores, folder, kofam_db, hmmsearch="hmmsearch", threads=1,
                  log=print, out_dir=None):
    """A decoy KEGG has no KO of (data/marker_decoys.tsv's row; since 2026-11-09) -> (<out_dir>/decoys/<decoy>.hmm
    (out_dir: the set being built, since 2026-11-10; None: folder, which holds the caches), its
    profile named decoy-<decoy>, and dict(seeds, members_over), the members scoring it over their own target: a
    paralog's protein filed as the marker's, the build drops them: 2026-11-09, 8 green algae's and Salpingoeca's
    "EF-1a", EFL); (None, set()) where under DECOY_MIN seeds are found. Its seeds: NCBI's proteins of its query
    (DECOY_SCOPE) whose best profile of every marker's KOs (kofam_db) is one of the marker's targets, scoring under the
    members' 5th percentile (member_scores, sorted) and unlike every member (under DECOY_UNLIKE of shared 5-mers, both
    ways: not the marker's fragment), aligned to the marker's first target profile (hmmalign --trim) and built
    (hmmbuild)."""
    from .benchmarks import kmers  # noqa: PLC0415
    tag = query_tag(row["ncbi"], DECOY_SCOPE)
    got = read_proteins(os.path.join(folder, decoy_proteins(row)))  # fetch_inputs fetched them
    seqs = {r["accession"]: r["sequence"] for r in got if r["sequence"]}
    scores = score(seqs, kofam_db, hmmsearch, threads,
                   os.path.join(folder, f"scores.decoy-{row['decoy']}.{tag}.tsv.gz"))
    targets, p5 = set(marker["targets"]), member_scores[int(CUTOFF_QUANTILE * (len(member_scores) - 1))]
    theirs = [kmers(r["sequence"]) for r in members_]
    seeds = {}
    for a, s in sorted(scores.items()):
        best = max(s, key=lambda k: s[k][0]) if s else None
        if best not in targets or s[best][0] >= p5:
            continue
        mine = kmers(seqs[a])
        if mine and all(len(mine & t) < DECOY_UNLIKE * min(len(mine), len(t)) for t in theirs if t):
            seeds[a] = seqs[a]
    log(f"markers: {marker['locus']}'s decoy {row['decoy']}: {len(seeds)} seeds of NCBI's {len(seqs)} proteins")
    if len(seeds) < DECOY_MIN:
        log(f"WARNING: markers: {marker['locus']}'s decoy {row['decoy']} not built (under {DECOY_MIN} seeds)")
        return None, set()
    out = os.path.join(out_dir or folder, "decoys", f"{row['decoy']}.hmm")
    with tempfile.TemporaryDirectory(dir=out_dir or folder) as tmp:
        faa, sto, pair = (os.path.join(tmp, n) for n in ("seeds.faa", "seeds.sto", "pair.hmm"))
        with open(faa, "w") as fh:
            fh.writelines(f">{a}\n{s}\n" for a, s in seeds.items())
        target = os.path.join(folder, "kofam", f"{marker['targets'][0]}.hmm")
        for cmd in ([sibling(hmmsearch, "hmmalign"), "--amino", "--trim", "-o", sto, target, faa],
                    [sibling(hmmsearch, "hmmbuild"), "--amino", "-n", f"decoy-{row['decoy']}", "-o", os.devnull,
                     util.part_of(out), sto]):
            run = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if run.returncode != 0:
                raise RuntimeError(f"{os.path.basename(cmd[0])} failed (exit {run.returncode}): "
                                   f"{run.stderr.strip()[-200:]}")
        os.replace(util.part_of(out), out)
        # the marker's members against it and its target: none should score it over the target (it took the marker)
        with open(pair, "w") as fh:
            for path in (target, out):
                with open(path) as src:
                    fh.write(src.read())
        mine = score({r["accession"]: r["sequence"] for r in members_}, pair, hmmsearch, threads)
    name = f"decoy-{row['decoy']}"
    own = {a: max((s for k, s in scores_of.items() if k != name), default=(0.0,))[0] for a, scores_of in mine.items()}
    over = {a for a, scores_of in mine.items() if scores_of.get(name, (0.0,))[0] > own.get(a, 0.0)}
    if over:  # a paralog's protein filed as the marker's (EFL as EF-1a): the build drops it
        log(f"markers: {len(over)} of {marker['locus']}'s members score its decoy {row['decoy']} over the target: "
            f"the decoy's, not members (e.g. {sorted(over)[:3]})")
    return dict(seeds=len(seeds), members_over=len(over)), over


def recheck(folder, hmmsearch="hmmsearch", threads=1):
    """The marker copies' check by the marker set's profiles, decoys too (since 2026-11-09; a genome run before a decoy
    joined its marker, EFL EF-1a's, kept its calls) -> a function of ({(accession, locus): coding sequence}, genomes
    {accession: row with gencode}) -> the keys whose best profile is a decoy's (package.collect's `paralog`), or None
    without a marker set. The search's scores are kept in the folder (since 2026-11-10: recheck.<digest>.tsv.gz, the
    digest of the profiles and of the proteins searched; until then every package run searched every marker copy
    again, a second run over the same genomes too); the last search's alone, and none where the folder is not to be
    written (a shared sources folder): searched as before then."""
    if not folder:
        return None

    def paralogs(cds, genomes):
        from .gencode import protein_of  # noqa: PLC0415
        if not cds:
            return set()
        names = {f"p{i}": key for i, key in enumerate(sorted(cds))}
        proteins = {n: protein_of(cds[k], genomes.get(k[0], {}).get("gencode", 1))[0].replace("*", "X")
                    for n, k in names.items()}
        wanted = {n: p for n, p in proteins.items() if p}
        profiles = []
        for sub in ("hmms", "decoys"):
            for f in sorted(os.listdir(os.path.join(folder, sub))):
                if f.endswith(".hmm"):
                    with open(os.path.join(folder, sub, f)) as fh:
                        profiles.append(fh.read())
        digest = hashlib.md5(("".join(profiles) + "\n" + "".join(f">{n}\n{p}\n" for n, p in wanted.items()))
                             .encode()).hexdigest()[:16]
        cache, got = os.path.join(folder, f"{RECHECK}{digest}.tsv.gz"), None
        if os.path.isfile(cache):
            try:
                got = score({}, "", cache=cache)
            except (OSError, ValueError, EOFError):  # not readable: searched again
                got = None
        if got is None:
            with tempfile.TemporaryDirectory() as tmp:
                db = os.path.join(tmp, "markers.hmm")
                with open(db, "w") as out:
                    out.write("".join(profiles))
                got = score(wanted, db, hmmsearch, threads)
            try:
                keep_scores(got, cache)
                for f in os.listdir(folder):
                    if f.startswith(RECHECK) and f.endswith(".tsv.gz") and f != os.path.basename(cache):
                        os.remove(os.path.join(folder, f))
            except OSError:
                pass
        out = set()
        for n, s in got.items():
            if s and max(s, key=lambda k: s[k][0]).startswith("decoy-"):
                out.add(names[n])
        return out
    return paralogs


def fetch(src_dir, download=None, log=print):
    """The marker dataset's inputs fetched (fetch_inputs; the fetch stage's) -> the marker folder. The set is built of
    them by the process stage (prepare; since 2026-11-11: here, its hmmsearch on the download node). A fetch that fails
    is an OSError (the fetch stage goes on without it), whatever failed (since 2026-11-08: an NCBI reply without its
    count, or of bad XML, stopped the fetch stage and every source after this one)."""
    try:
        return fetch_inputs(src_dir, TABLE, download, log, DECOYS)
    except Exception as e:  # noqa: BLE001
        raise OSError(f"marker set's inputs not fetched: {type(e).__name__}: {e}") from e


def prepare(src_dir, hmmsearch="hmmsearch", threads=1, log=print, table=TABLE, decoys=DECOYS, wait=WAIT):
    """The marker dataset built where it is not built from this table, decoy table, rules and inputs yet (ready()), of
    the inputs the fetch stage fetched (the process stage's, since 2026-11-11) -> its folder. One run builds it at a
    time (a lock in the marker folder: a SLURM array's tasks start together); the others wait for it and take its set,
    or build it themselves where its build failed. A build that fails is an OSError (the stage goes on without it, or
    with the set of before), whatever failed."""
    said = False
    while True:
        folder = ready(src_dir, table, decoys)
        if folder:
            return folder
        try:
            # the lock beside the marker folder, not in it: swap() moves the folder's files that are not the set's
            with util.single_run(src_dir, "markers"):
                folder = ready(src_dir, table, decoys)
                if folder:  # built by the run that held the lock
                    return folder
                missing = inputs_missing(src_dir, table, decoys)
                if missing:
                    raise OSError(f"marker set not built: {len(missing)} of its inputs not fetched into "
                                  f"{folder_of(src_dir)} (e.g. {missing[:3]}): --stage fetch fetches them")
                try:
                    return build(src_dir, table, hmmsearch, threads, log, decoys)
                except Exception as e:  # noqa: BLE001
                    raise OSError(f"marker set not built: {type(e).__name__}: {e}") from e
        except util.RunHeld as e:
            if not said:
                log(f"markers: another run builds the set ({e}): waiting for it")
                said = True
            time.sleep(wait)


def state_of(folder):
    """What a built marker dataset was made of, for the package manifest's inputs.markers (None for none): its table's,
    KOfam's ko_list's and members.tsv's MD5s, Swiss-Prot's release, NCBI's query scope and proteins fetched per
    marker, when it was built, and per marker its members and cutoff."""
    if not folder:
        return None
    with open(os.path.join(folder, STATE)) as fh:
        state = json.load(fh)
    members_tsv = os.path.join(folder, "members.tsv")
    return dict(table_md5=state.get("table_md5"), kofam_ko_list_md5=state.get("kofam_ko_list_md5"),
                swissprot=release_id(state.get("swissprot")), ncbi=state.get("ncbi"), built=state.get("built"),
                members_md5=util.md5(members_tsv) if os.path.isfile(members_tsv) else None,
                markers={locus: dict(kind=m.get("kind"), members=m.get("members"), cutoff=m.get("cutoff"),
                                     foreign=len(m.get("foreign", ())))
                         for locus, m in state.get("markers", {}).items()})


def kinds(folder):
    """{marker locus: clean or family} of a built marker dataset ({} for none)."""
    if not folder:
        return {}
    with open(os.path.join(folder, STATE)) as fh:
        return {locus: m.get("kind", "clean") for locus, m in json.load(fh).get("markers", {}).items()}
