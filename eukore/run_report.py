"""How a Eukore run went, and what to improve (since 2026-11-11; until 2026-11-12 tools/run_report.py): one report
from a build's work folder and, where given, the package it wrote. It only reads: nothing in WORK or the package is
changed (genomes.failed_table, which the build runs, removes the notes of genomes done since; this report leaves
them). The package stage writes it into its evaluation bundle (eukore_build.evaluation_bundle); any run can be
reported on again:

    python3 -m eukore.run_report --work WORK [--package OUT/Eukore_<v>] [--genomes DIR] [-t 16] [--out report.md]

The run (WORK):
- the plan's genomes (WORK/genome_set.tsv) done, and of those not done: failed (a note in WORK/failed/), running
  (a live claim), killed (a claim gone stale, or a partial done folder), not fetched (--genomes: no genome file there),
  else never tried; each by its accession's kind (an NCBI assembly, a transcriptome, records alone);
- the failures' reasons, counted (paths and quoted names left out, as eukore_build.failure_kinds does), by host and
  SLURM job: one node failing them all is the node's;
- of the done genomes: their done.json's rrna_failed, metadata_failed and tiara "failed" (the backfill tries again),
  an older genome record (run again, or --reversion), their BUSCOs' completeness (the poorest listed) and run time;
- the WORK's file system's free space.
The package (--package): benchmarks/summary.json (what was left out and why; completeness, by supergroup; the
OrthoDB verdicts, held out too; barcode gap), composition.tsv (supergroups mostly without a genome), flags.tsv (the
copies' confidence by supergroup) and completeness.tsv (the poorest genomes, the most duplicated); with it the run's
not-done genomes are also split by supergroup (species.tsv's group of genome_set.tsv's species).

Findings come first, each with its severity, what was seen and what to do; the tables after. Thresholds are the
constants below.
"""

import argparse
import collections
import glob
import json
import os
import re
import shutil
import statistics
import sys

from . import genomes, store

FREE_SHARE = 0.05        # a file system with less free is full for a build (the 2026-11-11 run: scratch 100%)
FREE_BYTES = 200e9       # or less than this free
FAILED_SHARE = 0.01      # not-done genomes over this share of the plan: a finding
HOST_SHARE = 0.5         # one host with over half the failures (of HOST_MIN or more): that node's
HOST_MIN = 20
POOR_COMPLETE = 0.5      # a done genome with under half its BUSCOs complete
POOR_SHARE = 0.05        # poor genomes over this share of the done: a finding
SUPERGROUP_COMPLETE = 0.7    # a supergroup's median completeness under this
CONCORDANT = 0.9         # OrthoDB's verdicts: concordant of called under this
HELD_OUT_CONCORDANT = 0.8
CALLED = 0.9
GAP_RESOLVED = 0.5       # a marker's barcode gap: under half the panel's species resolved
GENOME_SHARE = 0.5       # a supergroup (of SUPERGROUP_MIN species or more) with under half its species with a genome
SUPERGROUP_MIN = 20
LOW_SHARE = 0.1          # a supergroup's copies low or unverified over this share
DUPLICATED = 0.2         # a genome with this share of its BUSCOs duplicated: two genomes in one, or a polyploid
SHOWN = 15               # rows of each list

TSA = re.compile(r"^[A-Z]{4}\d{2}0{6}$")
PARTIAL = re.compile(r"^\..+?-\d+\.(.+)$")


def kind_of(accession):
    if accession.startswith(("GCA_", "GCF_")):
        return "assembly"
    if accession.startswith("EUK_"):
        return "records alone"
    if TSA.match(accession):
        return "transcriptome"
    return "other"


def reason_of(why):
    """A failure's reason without its paths and quoted names (eukore_build.failure_kinds' normalisation)."""
    return re.sub(r"'[^']*'|\"[^\"]*\"|/\S+", "...", " ".join(why.split()))[:120]


def rows_of(path):
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "-"


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


# --- the run --------------------------------------------------------------------------------------------------------

def failure_notes(work):
    """WORK/failed/<accession>.txt's rows -> {accession: dict(genomes.FAILED_COLUMNS)}."""
    notes = {}
    for path in glob.glob(os.path.join(work, genomes.FAILED, "*.txt")):
        try:
            with open(path) as fh:
                row = fh.readline().rstrip("\n").split("\t")
        except OSError:
            continue
        if len(row) == len(genomes.FAILED_COLUMNS):
            notes[row[0]] = dict(zip(genomes.FAILED_COLUMNS, row))
    return notes


def claims(work):
    """{accession: live?} of WORK/done's claims, and the accessions of its partial done folders."""
    done = os.path.join(work, "done")
    got = {os.path.basename(c)[1:-len(".claim")]: genomes.claim_live(c)
           for c in glob.glob(os.path.join(done, ".*.claim"))}
    partials = set()
    for p in glob.glob(os.path.join(done, ".*.partial")):
        name = os.path.basename(p)[:-len(".partial")]
        name = name[:-len(store.EXT)] if name.endswith(store.EXT) else name
        m = PARTIAL.match(name)  # .<host>-<pid>.<accession> (util.part_of)
        if m:
            partials.add(m.group(1))
    return got, partials


def run_facts(work, genomes_folder=None, threads=16):
    """The run's state: the plan's genomes, done ones' done.json, the not-done ones' fate."""
    planned = {r["accession"]: r["species_taxid"]
               for r in rows_of(os.path.join(work, genomes.GENOME_SET)) if r.get("accession")}
    with genomes.cached_summaries():
        genomes.prefetch(work, planned, threads=threads)
        summaries = {a: genomes.done_summary(work, a) for a in planned}
    notes = failure_notes(work)
    held, partials = claims(work)
    fate = {}
    for acc, summary in summaries.items():
        if summary is not None:
            fate[acc] = "done"
        elif held.get(genomes.safe(acc)):
            fate[acc] = "running"
        elif acc in notes:
            fate[acc] = "failed"
        elif genomes.safe(acc) in held or genomes.safe(acc) in partials:
            fate[acc] = "killed"
        elif genomes_folder and not os.path.isfile(genomes.genome_file(genomes_folder, acc)):
            fate[acc] = "not fetched"
        else:
            fate[acc] = "never tried"
    stale = sum(1 for live in held.values() if not live)
    return dict(planned=planned, summaries={a: s for a, s in summaries.items() if s is not None}, fate=fate,
                notes={a: n for a, n in notes.items() if fate.get(a) != "done"}, stale_claims=stale,
                partials=len(partials), stale_notes=sum(1 for a in notes if fate.get(a) == "done"))


def complete_of(summary):
    """A done genome's share of its searched BUSCOs Single or Duplicated, and of them Duplicated."""
    status = summary.get("status") or {}
    searched = len(summary.get("searched") or ()) or sum(status.values())
    if not searched:
        return None, None
    return ((status.get("Single", 0) + status.get("Duplicated", 0)) / searched,
            status.get("Duplicated", 0) / searched)


def free_space(path):
    try:
        u = shutil.disk_usage(path)
    except OSError:
        return None
    return u.free, u.total


# --- the report -----------------------------------------------------------------------------------------------------

class Report:
    def __init__(self):
        self.findings, self.sections = [], []

    def find(self, severity, what, todo):
        self.findings.append((("high", "medium", "low").index(severity), severity, what, todo))

    def section(self, title, *parts):
        self.sections.append(f"## {title}\n\n" + "\n\n".join(p for p in parts if p))

    def text(self, title):
        out = [f"# {title}", "", "## Findings", ""]
        if not self.findings:
            out.append("Nothing over the thresholds.")
        for i, (_r, sev, what, todo) in enumerate(sorted(self.findings, key=lambda f: f[0]), 1):
            out.append(f"{i}. **{sev}**: {what}  \n   *To do:* {todo}")
        return "\n".join(out + [""] + [s + "\n" for s in self.sections])


def report_run(rep, work, facts, supergroup_of=None):
    planned, fate, notes = facts["planned"], facts["fate"], facts["notes"]
    n = len(planned)
    if not n:
        rep.find("high", f"no plan in {work} ({genomes.GENOME_SET} missing or empty)",
                 "give the build's work folder (--build_dir), after a process stage has planned")
        return
    space = free_space(work)
    if space and (space[0] < FREE_BYTES or space[0] < FREE_SHARE * space[1]):
        rep.find("high", f"the WORK's file system has {space[0] / 1e9:,.0f} GB free of {space[1] / 1e9:,.0f} GB",
                 "free space (or set --scratch / the settings' SCRATCH elsewhere) before running again: a full disk "
                 "fails genomes and keeps their notes from being written")
    by_fate = collections.Counter(fate.values())
    by_kind = collections.defaultdict(collections.Counter)
    for acc, f in fate.items():
        by_kind[kind_of(acc)][f] += 1
    fates = ("done", "failed", "killed", "running", "not fetched", "never tried")
    rows = [[k] + [c[f] for f in fates] + [sum(c.values())] for k, c in sorted(by_kind.items())]
    rows.append(["all"] + [by_fate[f] for f in fates] + [n])
    parts = [table(["kind"] + list(fates) + ["planned"], rows)]
    not_done = n - by_fate["done"]
    if not_done > FAILED_SHARE * n:
        worst = max((k for k in by_kind), key=lambda k: sum(c for f, c in by_kind[k].items() if f != "done"))
        rep.find("high" if not_done > 5 * FAILED_SHARE * n else "medium",
                 f"{not_done:,} of {n:,} planned genomes not done ({pct(not_done, n)}): "
                 + ", ".join(f"{by_fate[f]:,} {f}" for f in fates[1:] if by_fate[f])
                 + f"; most of them of kind '{worst}'",
                 "run the process stage again (it resumes); look at the failures' reasons below first")
    if by_fate["killed"] or facts["stale_claims"]:
        rep.find("medium", f"{by_fate['killed']:,} genomes killed mid-run ({facts['stale_claims']:,} stale claims, "
                           f"{facts['partials']:,} partial done folders)",
                 "the next build removes the stale ones and runs them again; many at once point to a time or memory "
                 "limit (sacct's State TIMEOUT / OUT_OF_MEMORY)")
    if supergroup_of:
        sg = collections.defaultdict(collections.Counter)
        for acc, f in fate.items():
            sg[supergroup_of.get(planned[acc], "-")][f == "done"] += 1
        rows = sorted(([g, c[True] + c[False], c[False], pct(c[False], c[True] + c[False])] for g, c in sg.items()),
                      key=lambda r: -r[2])
        parts.append("Not done by supergroup:\n\n" + table(["supergroup", "planned", "not done", "share"], rows))
    rep.section("The run: the plan's genomes", *parts)

    if notes:
        reasons = collections.Counter(reason_of(r["why"]) for r in notes.values())
        hosts = collections.Counter(r["host"] for r in notes.values())
        jobs = collections.Counter(r["slurm_job"] for r in notes.values())
        kinds = collections.defaultdict(collections.Counter)
        for acc, r in notes.items():
            kinds[reason_of(r["why"])][kind_of(acc)] += 1
        top = [[c, why.replace("|", "\\|"), ", ".join(f"{k} {m}" for k, m in kinds[why].most_common())]
               for why, c in reasons.most_common(SHOWN)]
        when = sorted(r["when"] for r in notes.values())
        rep.section("The failures (WORK/failed)",
                    f"{len(notes):,} genomes' notes, {when[0]} .. {when[-1]}; {len(reasons):,} reasons.",
                    table(["genomes", "reason", "kinds"], top),
                    table(["host", "failures"], hosts.most_common(SHOWN)),
                    table(["SLURM job", "failures"], jobs.most_common(SHOWN)),
                    "Examples: " + ", ".join(sorted(notes)[:SHOWN]))
        why, c = reasons.most_common(1)[0]
        rep.find("high" if c > 0.5 * len(notes) else "medium",
                 f"the commonest failure ({c:,} of {len(notes):,}): `{why}`",
                 "fix its cause (a full disk: free space; a tool's error: run one of its genomes by hand with "
                 "--species/--accessions and keep its tmp)")
        host, hc = hosts.most_common(1)[0]
        if hc >= HOST_MIN and hc > HOST_SHARE * len(notes) and len(hosts) > 1:
            rep.find("medium", f"{hc:,} of {len(notes):,} failures on one node, {host}",
                     "exclude the node (sbatch --exclude) and look at its local disk and memory")
    if facts["stale_notes"]:
        rep.section("Notes of genomes done since", f"{facts['stale_notes']:,} failure notes are of genomes done "
                    "since (another build): the next build's failed_table removes them.")

    done = facts["summaries"]
    if not done:
        return
    flags = collections.Counter()
    for s in done.values():
        flags["rrna_failed"] += bool(s.get("rrna_failed"))
        flags["metadata_failed"] += bool(s.get("metadata_failed"))
        flags["tiara failed"] += s.get("tiara") == "failed"
        flags["older genome record"] += int(s.get("genome_record", 1)) < genomes.GENOME_RECORD
        flags["no rescue pass"] += "rescue" in s and not s.get("rescue")
    shares = {a: complete_of(s) for a, s in done.items()}
    complete = {a: c for a, (c, _d) in shares.items() if c is not None}
    poor = sorted((c, a) for a, c in complete.items() if c < POOR_COMPLETE)
    seconds = sorted(((s.get("seconds") or 0, a) for a, s in done.items()), reverse=True)
    total = sum(t for t, _a in seconds)
    by_kind = collections.defaultdict(list)
    for a, c in complete.items():
        by_kind[kind_of(a)].append(c)
    rep.section(
        "The done genomes",
        table(["check", "genomes"], [[k, f"{v:,}"] for k, v in flags.items()]),
        table(["kind", "genomes", "median complete", f"under {POOR_COMPLETE:.0%}"],
              [[k, len(v), f"{statistics.median(v):.3f}", sum(x < POOR_COMPLETE for x in v)]
               for k, v in sorted(by_kind.items())]),
        "Poorest:\n\n" + table(["accession", "complete"], [[a, f"{c:.3f}"] for c, a in poor[:SHOWN]]),
        f"Run time: {total / 3600:,.0f} hours of genome runs (wall seconds; times threads for CPU); slowest:\n\n"
        + table(["accession", "hours", "threads", "Gb"],
                [[a, f"{t / 3600:.1f}", done[a].get("threads", "-"), f"{(done[a].get('bytes') or 0) / 1e9:.2f}"]
                 for t, a in seconds[:SHOWN]]))
    for what, todo in (("rrna_failed", "eukore_build.py --backfill tries their rRNA again"),
                       ("metadata_failed", "the backfill tries their metadata again"),
                       ("tiara failed", "the backfill classifies their windows again"),
                       ("older genome record", "run them again, or --reversion where the change does not touch them")):
        if flags[what]:
            rep.find("low" if flags[what] < 0.01 * len(done) else "medium",
                     f"{flags[what]:,} done genomes with {what}", todo)
    assemblies = [a for a in complete if kind_of(a) == "assembly"]
    poor_assemblies = [a for c, a in poor if kind_of(a) == "assembly"]
    if assemblies and len(poor_assemblies) > POOR_SHARE * len(assemblies):
        rep.find("medium", f"{len(poor_assemblies):,} of {len(assemblies):,} done assemblies with under "
                           f"{POOR_COMPLETE:.0%} of their BUSCOs complete",
                 "look at them by supergroup (the package's completeness.tsv): a clade's divergent genes (the "
                 "rescue pass, nearer query proteins) or poor assemblies (a better one, genomes.rank)")


# --- the package ----------------------------------------------------------------------------------------------------

def report_package(rep, folder):
    bench = os.path.join(folder, "benchmarks")
    path = os.path.join(bench, "summary.json")
    if not os.path.isfile(path):
        rep.find("high", f"no benchmarks in {folder} (benchmarks/summary.json)",
                 f"python3 -m eukore.benchmarks {folder} --build_dir WORK -t 8")
        return
    with open(path) as fh:
        s = json.load(fh)
    parts = [f"{s.get('species', 0):,} species with a genome, {s.get('busco_loci', 0):,} BUSCOs."]
    left = s.get("left_out") or {}
    if left:
        parts.append("Left out:\n\n" + table(["benchmark", "why"], sorted(left.items())))
        rep.find("medium", f"{len(left)} benchmarks left out: {', '.join(sorted(left))}",
                 "give their inputs (--build_dir for BUSCO's dataset and the run times, vsearch on PATH) and measure "
                 "again: python3 -m eukore.benchmarks PACKAGE --build_dir WORK")
    comp = s.get("completeness") or {}
    weak = {g: m for g, m in (comp.get("by_supergroup") or {}).items() if m is not None and m < SUPERGROUP_COMPLETE}
    parts.append(f"Completeness: median {comp.get('median')}, {comp.get('at_least_90')} of genomes 90% or more.")
    if weak:
        parts.append(table(["supergroup", "median complete"], sorted(weak.items(), key=lambda kv: kv[1])))
        rep.find("medium", "supergroups with median completeness under "
                           f"{SUPERGROUP_COMPLETE:.0%}: " + ", ".join(f"{g} {m:.2f}" for g, m in
                                                                   sorted(weak.items(), key=lambda kv: kv[1])),
                 "their genes are far from BUSCO's query proteins: the rescue pass (HMM) and the held-out benchmark "
                 "say whether calls are missed or the genomes poor")
    odb = s.get("orthodb") or {}
    for label, d, want in (("the build's calls", odb, CONCORDANT), ("the held-out calls", odb.get("held_out") or {},
                                                                    HELD_OUT_CONCORDANT)):
        if not d:
            continue
        parts.append(f"OrthoDB, {label}: called {d.get('called')}, concordant of called "
                     f"{d.get('concordant_of_called')}, median recall {d.get('median_recall')}; {d.get('verdicts')}")
        if d.get("concordant_of_called") is not None and d["concordant_of_called"] < want:
            rep.find("medium", f"OrthoDB: {label} concordant {d['concordant_of_called']} of called (under {want})",
                     "orthodb.tsv's verdicts by supergroup and locus: a wrong paralog or a partial model")
        if d.get("called") is not None and d["called"] < CALLED:
            rep.find("medium", f"OrthoDB: {label} called {d['called']} of OrthoDB's genes (under {CALLED})",
                     "orthodb.tsv's missed rows: the seeding (miniprot's queries) or the cut-offs")
    gap = s.get("barcode_gap") or {}
    poor_gap = {m: v["resolved"] for m, v in gap.items() if v.get("resolved") is not None
                and v["resolved"] < GAP_RESOLVED}
    if poor_gap:
        rep.find("low", "markers resolving under half the panel's species: "
                        + ", ".join(f"{m} {v:.2f}" for m, v in sorted(poor_gap.items(), key=lambda kv: kv[1])),
                 "expected for the 18S; for a gene, its copies (resolution.tsv)")
    rep.section(f"The package: {os.path.basename(os.path.normpath(folder))}", *parts)

    parts = []
    comp_rows = rows_of(os.path.join(bench, "composition.tsv"))
    thin = [r for r in comp_rows if r["supergroup"] != "all" and int(r["species"]) >= SUPERGROUP_MIN
            and int(r["with_genome"]) < GENOME_SHARE * int(r["species"])]
    if comp_rows:
        parts.append(table(["supergroup", "species", "with genome", "rRNA only", "median complete", "operon",
                            "organelle"],
                           [[r["supergroup"], r["species"], r["with_genome"], r["rrna_only"], r["complete_median"],
                             r["with_operon"], r["with_organelle"]] for r in comp_rows]))
    if thin:
        rep.find("low", "supergroups mostly of rRNA alone: " + ", ".join(
            f"{r['supergroup']} ({r['with_genome']} of {r['species']} with a genome)" for r in thin),
                 "transcriptomes (TSA) or no-assembly records where NCBI has no genome; the gene markers miss these "
                 "species")
    flag_rows = rows_of(os.path.join(bench, "flags.tsv"))
    if flag_rows:
        by_group = collections.defaultdict(collections.Counter)
        reasons = collections.Counter()
        for r in flag_rows:
            n = int(r["copies"])
            by_group[r["supergroup"]][r["confidence"]] += n
            if r["confidence"] in ("low", "unverified"):
                for f in r["flags"].split(","):
                    if f != "-":
                        reasons[f] += n
        rows = []
        for g, c in sorted(by_group.items()):
            total = sum(c.values())
            bad = c["low"] + c["unverified"]
            rows.append([g, total, c["high"], c["medium"], c["unverified"], c["low"], pct(bad, total)])
            if total >= 100 and bad > LOW_SHARE * total:
                rep.find("low", f"{g}: {pct(bad, total)} of its BUSCO copies low or unverified",
                         "busco/copies.tsv's flags of its genomes: own_species f is another clade's gene (a "
                         "contaminated or misnamed genome), undecided a clade the others cannot tell")
        parts.append("The BUSCO copies' confidence:\n\n"
                     + table(["supergroup", "copies", "high", "medium", "unverified", "low", "low+unverified"], rows))
        if reasons:
            parts.append("Why low or unverified:\n\n" + table(["flag", "copies"], reasons.most_common(SHOWN)))
    genome_rows = rows_of(os.path.join(bench, "completeness.tsv"))
    if genome_rows:
        def f(r, k):
            try:
                return float(r[k])
            except (KeyError, ValueError):
                return 0.0
        dup = sorted((r for r in genome_rows if f(r, "loci") and f(r, "duplicated") / f(r, "loci") >= DUPLICATED),
                     key=lambda r: -f(r, "duplicated") / f(r, "loci"))
        poorest = sorted(genome_rows, key=lambda r: f(r, "complete"))[:SHOWN]
        parts.append("Poorest genomes:\n\n" + table(
            ["accession", "name", "supergroup", "source", "level", "complete", "missing"],
            [[r["accession"], r["name"], r["supergroup"], r["source"], r["level"], r["complete"], r["missing"]]
             for r in poorest]))
        if dup:
            parts.append(f"Duplicated {DUPLICATED:.0%} or more of their BUSCOs ({len(dup):,}):\n\n" + table(
                ["accession", "name", "supergroup", "duplicated", "loci"],
                [[r["accession"], r["name"], r["supergroup"], r["duplicated"], r["loci"]] for r in dup[:SHOWN]]))
            rep.find("low", f"{len(dup):,} genomes with {DUPLICATED:.0%} or more of their BUSCOs duplicated",
                     "a hybrid, a polyploid or two organisms in one assembly: their copies' own_species and cross "
                     "flags say which")
    rep.section("The package's species and copies", *parts)


def supergroups(folder):
    return {r["taxid"]: r.get("group") or r.get("supergroup", "-")
            for r in rows_of(os.path.join(folder, "species.tsv"))}


def report(work, package=None, genomes_folder=None, threads=16):
    """The report of a run (WORK) and the package it wrote, as Markdown."""
    rep = Report()
    facts = run_facts(work, genomes_folder, threads)
    report_run(rep, work, facts, supergroups(package) if package else None)
    if package:
        report_package(rep, package)
    return rep.text(f"Eukore run report: {work}" + (f", {package}" if package else ""))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--work", required=True, help="the build's work folder (--build_dir)")
    p.add_argument("--package", help="the package it wrote (OUT/Eukore_<version>)")
    p.add_argument("--genomes", help="the genomes' folder (OUT/genomes): the not-done genomes never fetched told apart")
    p.add_argument("-t", "--threads", type=int, default=16, help="done.json read at once (shared file system)")
    p.add_argument("--out", help="the report (Markdown); else printed")
    a = p.parse_args(argv)
    text = report(a.work, a.package, a.genomes, a.threads)
    if a.out:
        with open(a.out, "w", newline="\n") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
