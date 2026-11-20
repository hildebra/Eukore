"""protal's genome sketch of the package's species (since 2026-11-19; eukore_build.py --protal BIN): made by the package
stage with protal's own `BIN sketch` while the genomes are still on disk (before the prune deletes them), and shipped
as the package's metadata/protal_sketch.bin. protal (the metagenome profiler) merges this fixed eukaryote sketch into
the sketch of its own database builds, so that the eukaryote part is made once per Eukore release; Eukore never reads
it.

The sketch's species (table_rows): each species' representative genome (genomes.tsv's representative t) whose file is
in the genomes' folder; not a transcriptome (source tsa) or EukProt's records (no_assembly), nor a genome whose copies
the package leaves out (tiers poor, partial and contaminated, unless --keep_poor_genomes: package.collect's dropped).
Its key is the representative's accession as genomes.tsv has it, which protal maps to its own species names. The
species' further assemblies still on disk (source assembly or intraspecific, and since 2026-11-20 a genome merged into
the species: an isolate-named species' genome, same_species.py) follow it under the same key: protal keeps their
windows out of the other species' sketches, not out of their own species'. The windows on the BUSCO and marker copies
are left out (protal's alignment covers those genes): each copy's span on its contig, first to last codon, introns and
all, from its genome's loci.tsv and alternates.tsv in WORK/done (positions; a copy joined across contigs gives its best
piece's span alone, as loci.tsv does). A species' entries scale with its representative's size (entries()).

Its rerun cache (build), WORK/protal_sketch/: genomes.tsv and gene_positions.tsv (protal's --genomes and
--gene_positions), key.json (KEY: the binary's real path, size and mtime, the entries' options, the two files' MD5s;
and what made the sketch), sketch.log (protal's output) and protal_sketch.bin. A rerun with the same key reuses the
sketch; with another it is made again, unless genomes of the kept sketch's table are gone since (a prune): it cannot be
made again from them, and the kept one serves, said stale and why.
"""

import collections
import datetime
import functools
import json
import os
import shutil
import subprocess

from . import store, util

FOLDER = "protal_sketch"           # in WORK
SKETCH = "protal_sketch.bin"
IN_PACKAGE = "metadata/protal_sketch.bin"
TABLE, POSITIONS, KEY, LOG = "genomes.tsv", "gene_positions.tsv", "key.json", "sketch.log"
NEXT = ".next"                     # a table being written, before the key says whether it is taken
PARTIAL = ".partial"               # the sketch being written by protal
# Entries per species grow with its representative genome's size (since 2026-11-19; protal counts the bases as it
# reads the genome): KMERS x (its size / MB Mb)^EXPONENT, at least KMERS, with no upper bound (entries(); since
# 2026-11-20, protal 2.3.10's --kmers_exponent; KMERS at 500 Mb in proportion, KMERS / 5 to 20 x KMERS, before): a
# prokaryote's 5,000 at a prokaryote's size, 6,600 for a yeast of 12 Mb, 16,700 for 500 Mb, 26,200 for 3 Gb. For all of
# Eukore (33,341 species, 2026-11-17) ~490 M entries, ~5.9 GB on disk, ~37 GB of memory while protal sketch runs (4
# candidates of 19 bytes an entry); protal's database build takes them all (up to 20 x its prokaryotes' entries).
KMERS = 5000
MB = 4
EXPONENT = 0.25
NOT_GENOMES = ("tsa", "no_assembly")  # a transcriptome, EukProt's records: no genome to sketch
# a species' further assemblies (genomes.tsv's representative f); since 2026-11-20 a genome merged into a named species
# (same_species.SOURCE): the named species' further genome, never a species keeping its windows out of that one's
FURTHER = ("assembly", "intraspecific", "merged")
POSITION_COLUMNS = ("contig", "start", "end", "genome", "locus")  # protal reads contig, start, end (1-based, inclusive)
SUMMARY = "protal sketch: "        # protal's summary line starts so and holds "in the sketch ("
VERSION_SECONDS = 60


def table_rows(genomes, genome_file, left_out=()):
    """The sketch's genomes -> (rows [(species key, accession, path)]: each species' representative first, then its
    further assemblies on disk, the species in their representatives' order; missing [the representatives taken
    without a genome file]). genomes: package.collect's {accession: row} (the intraspecific panel's set by then);
    genome_file: accession -> its path; left_out: the accessions whose copies the package leaves out (collect's
    dropped): their species are not taken."""
    left_out = set(left_out)
    reps, further = {}, collections.defaultdict(list)
    for acc, g in sorted(genomes.items()):
        sp = int(g["species_taxid"])
        if g.get("representative") == "t" and g.get("source") not in NOT_GENOMES and acc not in left_out \
                and sp not in reps:
            reps[sp] = acc
        elif g.get("representative") != "t" and g.get("source") in FURTHER:
            further[sp].append(acc)
    rows, missing = [], []
    for sp, acc in sorted(reps.items(), key=lambda kv: kv[1]):
        path = genome_file(acc)
        if not os.path.isfile(path):
            missing.append(acc)
            continue
        rows.append((acc, acc, path))
        rows += [(acc, other, genome_file(other)) for other in further.get(sp, ())
                 if os.path.isfile(genome_file(other))]
    return rows, missing


def spans_of(got):
    """A done genome's tables (store.members' {name: bytes}: loci.tsv, alternates.tsv) -> its copies' spans
    {(contig, start, end, locus)}, the Missing loci and those without a place left out."""
    out = set()
    for name, data in got.items():
        for r in store.table_of(store.decoded(name, data)):
            if r.get("status") not in (None, "Missing") and str(r.get("start", "")).isdigit() \
                    and str(r.get("end", "")).isdigit() and r.get("contig") not in (None, "", "-"):
                out.add((r["contig"], int(r["start"]), int(r["end"]), r["locus"]))
    return out


def positions(work, accessions, path, threads=util.READ_THREADS):
    """The copies' spans on the genomes (accessions, in WORK/done) into path (POSITION_COLUMNS, a header first) -> (the
    spans written, the accessions without a done genome's loci.tsv: none of their spans known)."""
    from .genomes import ALTERNATES, done_dir  # noqa: PLC0415
    read = functools.partial(store.members, wanted=("loci.tsv", ALTERNATES[0]))
    n, without = 0, []
    with open(path, "w", newline="\n") as fh:
        fh.write("\t".join(POSITION_COLUMNS) + "\n")
        for acc, got in zip(accessions, util.threaded(read, [done_dir(work, a) for a in accessions], threads)):
            if not got or got.get("loci.tsv") is None:
                without.append(acc)
                continue
            spans = sorted(spans_of(got))
            fh.writelines(f"{contig}\t{start}\t{end}\t{acc}\t{locus}\n" for contig, start, end, locus in spans)
            n += len(spans)
    return n, without


def entries(kmers=KMERS, mb=MB, exponent=EXPONENT):
    """protal sketch's options of the entries per species: kmers x (the genome's size / mb Mb)^exponent, at least kmers,
    no upper bound (kmers_max 0) -> {kmers, kmers_mb, kmers_exponent, kmers_min, kmers_max} (the key's and the
    manifest's)."""
    kmers, mb = int(kmers), int(mb)
    return dict(kmers=kmers, kmers_mb=mb, kmers_exponent=float(exponent), kmers_min=kmers, kmers_max=0)


def binary_state(path):
    """The protal binary as the key takes it: its real path, size and mtime."""
    real = os.path.realpath(path)
    st = os.stat(real)
    return dict(path=real, bytes=st.st_size, mtime=int(st.st_mtime))


def version_of(binary):
    """`BIN --version`'s first line, or None where it does not run."""
    try:
        done = subprocess.run([binary, "--version"], capture_output=True, text=True, timeout=VERSION_SECONDS)
    except (OSError, subprocess.SubprocessError):
        return None
    lines = [s.strip() for s in (done.stdout or "").splitlines() + (done.stderr or "").splitlines() if s.strip()]
    return lines[0] if done.returncode == 0 and lines else None


def read_key(folder):
    """The kept sketch's key.json ({key, about}), or None."""
    try:
        with open(os.path.join(folder, KEY)) as fh:
            got = json.load(fh)
        return got if isinstance(got, dict) and "key" in got else None
    except (OSError, ValueError):
        return None


def differences(old, new):
    """What changed between two keys, in words."""
    out = []
    if old.get("protal") != new.get("protal"):
        o, n = old.get("protal") or {}, new.get("protal") or {}
        out.append(f"another protal ({n.get('path')})" if o.get("path") != n.get("path") else
                   f"protal {n.get('path')} changed since")
    if old.get("entries") != new.get("entries"):
        o, n = old.get("entries") or {}, new.get("entries") or {}
        out.append(", ".join(f"--{k} {o.get(k)} -> {n.get(k)}" for k in dict.fromkeys([*n, *o])
                             if o.get(k) != n.get(k)))
    if old.get("table_md5") != new.get("table_md5"):
        out.append("another genome table")
    if old.get("positions_md5") != new.get("positions_md5"):
        out.append("other gene positions")
    return ", ".join(out) or "its key"


def table_paths(path):
    """The genome files a table names."""
    try:
        with open(path) as fh:
            return [line.rstrip("\n").split("\t")[1] for line in fh if "\t" in line]
    except OSError:
        return []


def summary_of(log_path):
    """protal's summary line in its log (the last starting with SUMMARY and holding `in the sketch (`), or None."""
    try:
        with open(log_path, errors="replace") as fh:
            lines = [s.strip() for s in fh if s.startswith(SUMMARY) and "in the sketch (" in s]
    except OSError:
        return None
    return lines[-1] if lines else None


def tail_of(path, n=5):
    try:
        with open(path, errors="replace") as fh:
            return " | ".join(s.strip() for s in fh.readlines()[-n:] if s.strip())
    except OSError:
        return ""


def about_of(kept, state, stale=None):
    """The manifest's inputs.protal_sketch of a kept sketch (key.json's) -> dict."""
    out = dict(file=IN_PACKAGE, protal=kept["key"]["protal"]["path"], version=kept["about"].get("version"),
               **(kept["key"].get("entries") or {}), species=kept["about"].get("species"),
               genomes=kept["about"].get("genomes"), gene_spans=kept["about"].get("gene_spans"),
               made=kept["about"].get("made"), state=state)
    if stale:
        out["stale"] = stale
    return out


def build(binary, genomes, genome_file, work, kmers=KMERS, mb=MB, threads=1, left_out=(), log=print, exponent=EXPONENT):
    """protal's sketch of the package's species (module docstring), made by binary in WORK/protal_sketch/ or kept from
    an earlier run -> (its path, or None where there is none: a failure is said, never raised past here; the
    manifest's inputs.protal_sketch, or None). genomes, genome_file, left_out: table_rows'; kmers, mb, exponent:
    entries'."""
    folder = os.path.join(work, FOLDER)
    os.makedirs(folder, exist_ok=True)
    sketch, partial = os.path.join(folder, SKETCH), os.path.join(folder, SKETCH + PARTIAL)
    rows, missing = table_rows(genomes, genome_file, left_out)
    keys = list(dict.fromkeys(k for k, _a, _p in rows))
    log(f"protal sketch: {len(keys):,} species' genomes on disk ({len(rows) - len(keys):,} further assemblies of "
        f"them)" + (f"; {len(missing):,} species of the package without their genome file (pruned?), not in it "
                    f"(e.g. {missing[:3]})" if missing else ""))
    nxt = {name: os.path.join(folder, name + NEXT) for name in (TABLE, POSITIONS)}
    with open(nxt[TABLE], "w", newline="\n") as fh:
        fh.writelines(f"{k}\t{p}\n" for k, _a, p in rows)
    n_spans, without = positions(work, [a for _k, a, _p in rows], nxt[POSITIONS])
    reps_without = [a for a in without if a in keys]
    log(f"protal sketch: {n_spans:,} BUSCO and marker copies' spans left out of the windows"
        + (f"; no copies known (no done loci.tsv) of {len(reps_without):,} representatives (e.g. {reps_without[:3]}), "
           "the sketch made without their spans" if reps_without else "")
        + (f"; {len(without) - len(reps_without):,} further assemblies without (none called)"
           if len(without) > len(reps_without) else ""))
    options = entries(kmers, mb, exponent)
    key = dict(protal=binary_state(binary), entries=options, table_md5=util.md5(nxt[TABLE]),
               positions_md5=util.md5(nxt[POSITIONS]))
    kept = read_key(folder)
    if kept and os.path.isfile(sketch):
        if kept["key"] == key:
            for path in nxt.values():
                os.remove(path)
            log(f"protal sketch: {sketch} kept (its protal, entries, genomes and gene positions the same; made "
                f"{kept['about'].get('made')})")
            return sketch, about_of(kept, "reused")
        gone = [p for p in table_paths(os.path.join(folder, TABLE)) if not os.path.isfile(p)]
        if gone:  # its genomes pruned since: it cannot be made again from them
            why = differences(kept["key"], key)
            for path in nxt.values():
                os.remove(path)
            log(f"WARNING: protal sketch: {sketch} kept though stale ({why}): {len(gone):,} genomes of its table are "
                f"gone since it was made ({kept['about'].get('made')}; e.g. {gone[:2]}), so it cannot be made again")
            return sketch, about_of(kept, "stale", f"{why}; {len(gone):,} genomes of its table gone")
    for path in (sketch, partial, os.path.join(folder, KEY)):  # made anew: the kept one is not this key's
        if os.path.lexists(path):
            os.remove(path)
    for name, path in nxt.items():
        os.replace(path, os.path.join(folder, name))
    if not rows:
        log("WARNING: protal sketch: no genome of the package's species on disk: none made")
        return None, None
    cmd = [binary, "sketch", "--genomes", os.path.join(folder, TABLE), "--out", partial,
           "--gene_positions", os.path.join(folder, POSITIONS)] \
        + [a for k, v in options.items() for a in (f"--{k}", str(v))] + ["-t", str(max(1, threads))]
    log("protal sketch: " + " ".join(cmd) + f" (its output: {os.path.join(folder, LOG)})")
    version = version_of(binary)
    try:
        with open(os.path.join(folder, LOG), "w") as out:
            code = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL).returncode
        why = None if code == 0 else f"exit {code}"
    except OSError as e:
        why = f"{type(e).__name__}: {e}"
    if why is None and not os.path.isfile(partial):
        why = f"no {os.path.basename(partial)} written"
    if why:
        if os.path.lexists(partial):
            os.remove(partial)
        log(f"WARNING: protal sketch failed ({why}): no protal sketch in the package; protal's output "
            f"{os.path.join(folder, LOG)}: {tail_of(os.path.join(folder, LOG))[-500:]}")
        return None, None
    os.replace(partial, sketch)
    summary = summary_of(os.path.join(folder, LOG))
    if summary:
        log("  " + summary)
    kept = dict(key=key, about=dict(version=version, species=len(keys), genomes=len(rows), gene_spans=n_spans,
                                    made=datetime.datetime.now().isoformat(timespec="seconds"), summary=summary))
    with open(os.path.join(folder, KEY + NEXT), "w", newline="\n") as fh:
        json.dump(kept, fh, indent=1)
    os.replace(os.path.join(folder, KEY + NEXT), os.path.join(folder, KEY))
    log(f"protal sketch: {sketch} made ({os.path.getsize(sketch):,} bytes)")
    return sketch, about_of(kept, "made")


def into_package(sketch, part):
    """A package writer (package.write's writers): the sketch as the unfinished package's metadata/protal_sketch.bin,
    a hard link to it (a copy where the file system makes none), before the manifest lists its size and MD5."""
    dest = os.path.join(part, IN_PACKAGE)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.lexists(dest):
        os.remove(dest)
    try:
        os.link(sketch, dest)
    except OSError:
        shutil.copyfile(sketch, dest)
