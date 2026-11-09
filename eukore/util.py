#!/usr/bin/env python3
"""Eukore's small helpers (since 2026-11-05; copied here so that Eukore imports nothing from protal's scripts/): text
files plain, gzipped or zstd-compressed (open_text; zstd by Python 3.14's compression.zstd where there is one, else the
zstd command through a pipe), a writer's temporary name, MD5s, downloads, a long loop's progress.
"""

import gzip
import shutil
import subprocess

ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"


def is_zstd(head):
    """Whether bytes start a zstd frame or a skippable frame."""
    return head[:4] == ZSTD_MAGIC or (len(head) >= 4 and head[0] & 0xf0 == 0x50 and head[1:4] == b"\x2a\x4d\x18")


def _stdlib_zstd():
    try:
        from compression import zstd  # Python 3.14
        return zstd
    except ImportError:
        return None


def _zstd_command():
    command = shutil.which("zstd")
    if not command:
        raise OSError("zstd files need Python 3.14 (compression.zstd) or the zstd command on PATH")
    return command


def read_zstd(path):
    """A zstd file's decompressed content, as text (compression.zstd, else the zstd command: raises OSError where it
    failed)."""
    zstd = _stdlib_zstd()
    if zstd is not None:
        with zstd.open(path, "rb") as fh:
            return fh.read().decode()
    run = subprocess.run([_zstd_command(), "-dcq", path], capture_output=True)
    if run.returncode:
        raise OSError(f"zstd -dcq {path} failed ({run.returncode}): {run.stderr.decode(errors='replace')[-300:]}")
    return run.stdout.decode()


# --- Eukore's own helpers (since 2026-11-05; self-contained: nothing imported from protal's scripts/) ---------------

import collections  # noqa: E402
import concurrent.futures  # noqa: E402
import contextlib  # noqa: E402
import datetime  # noqa: E402
import glob  # noqa: E402
import hashlib  # noqa: E402
import itertools  # noqa: E402
import os  # noqa: E402
import socket  # noqa: E402
import time  # noqa: E402
import urllib.error  # noqa: E402
import urllib.parse  # noqa: E402
import urllib.request  # noqa: E402
import http.client  # noqa: E402
import zlib  # noqa: E402

import threading  # noqa: E402

USER_AGENT = "Eukore-build (https://protal.earlham.ac.uk)"
ATTEMPTS = 5       # tries of a download before it is given up on
RETRY_WAIT = 2.0   # seconds before the second try, doubled for each further one (or the server's Retry-After)
# set on a stop (Ctrl-C: eukore_build.at_a_time): a download not tried again (until 2026-11-06 curl's SIGINT was a
# transient failure, and a 15 GB genome began again while the stop waited for it)
STOP = threading.Event()


def part_of(path):
    """A writer's own temporary name for `path` (renamed into place once written): .<host>-<pid>.<name>.partial beside
    it."""
    folder, name = os.path.split(path)
    return os.path.join(folder, f".{socket.gethostname()}-{os.getpid()}.{name}.partial")


PARTIAL_STALE = 2 * 86400  # another host's partial untouched this long is a dead run's (dead_partials)


def dead_partials(folder, name):
    """The partials (part_of's .<host>-<pid>.<name>.partial) in folder of runs that ended: this host's whose process is
    gone, another host's untouched for PARTIAL_STALE (since 2026-11-08: a run deleted every one, another build's going
    on too) -> their paths."""
    out, host = [], socket.gethostname()
    for path in glob.glob(os.path.join(glob.escape(folder), f".*.{glob.escape(name)}.partial")):
        their_host, _dash, pid = os.path.basename(path)[1:-len(f".{name}.partial")].rpartition("-")
        if their_host == host and pid.isdigit():
            if not pid_alive(int(pid)):
                out.append(path)
            continue
        try:
            if time.time() - os.path.getmtime(path) >= PARTIAL_STALE:
                out.append(path)
        except OSError:  # gone meanwhile
            pass
    return out


def pid_alive(pid):
    """Whether a process of this host is alive (another user's too)."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        return True
    return True


def ends_in_newline(path):
    """Whether a file's last byte is a line end (a cache a killed run appended to: its last line cut where not)."""
    with open(path, "rb") as fh:
        fh.seek(-1, os.SEEK_END)
        return fh.read(1) == b"\n"


def md5(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path):
    """A text file, plain, gzipped or zstd-compressed (by its first bytes; a zstd file read whole)."""
    import io  # noqa: PLC0415
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"\x1f\x8b":
        return gzip.open(path, "rt")
    if is_zstd(head):
        return io.StringIO(read_zstd(path))
    return open(path)


# fetch_file's transfers where there: ~3x urllib's rate (2026-11-06, NBI's HPC), else urllib
CURL = shutil.which("curl")


def curl_get(url, path):
    """url into path by curl (a stall of a minute is a failure; a transfer cut short of its announced length curl's own
    error): raises urllib.error.HTTPError (a status of 400 or more) or OSError (curl failed)."""
    run = subprocess.run([CURL, "-sS", "-L", "-A", USER_AGENT, "--connect-timeout", "60", "-y", "60", "-Y", "1024",
                          "-o", path, "-w", "%{http_code}", url], capture_output=True, text=True)
    if run.returncode:
        raise OSError(f"curl exit {run.returncode}: {run.stderr.strip()}")
    code = int(run.stdout.strip() or 0)
    if code >= 400:
        raise urllib.error.HTTPError(url, code, f"HTTP {code}", None, None)


# fetch_file's gzip check where there: ISA-L's igzip (conda-forge isa-l, apt isal; ~3.5x Python's zlib, 2026-11-06,
# WSL), else pigz (~1.1x), else Python's gzip
GZIP_TEST = shutil.which("igzip") or shutil.which("pigz")


def gzip_whole(path, tool=None):
    """Reads a gzip file to its end, by tool's -t (GZIP_TEST; from stdin: no suffix to check) or Python's gzip: raises
    (EOFError, OSError, zlib.error) where it is cut short or damaged."""
    if tool:
        with open(path, "rb") as fh:
            run = subprocess.run([tool, "-t"], stdin=fh, capture_output=True, text=True)
        if run.returncode:
            raise OSError(f"{os.path.basename(tool)} -t: {run.stderr.strip() or f'exit {run.returncode}'}")
        return
    with gzip.open(path, "rb") as gz:
        while gz.read(1 << 22):
            pass


class KeepAuthorization(urllib.request.HTTPRedirectHandler):
    """urllib's redirects, a request's Authorization (an API key) dropped where one leaves its host or https (since
    2026-11-08: urllib passes every header on, to whatever host a server redirects to)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        new = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new is not None:
            old, to = urllib.parse.urlsplit(req.full_url), urllib.parse.urlsplit(new.full_url)
            if to.scheme != "https" or (to.hostname or "").lower() != (old.hostname or "").lower():
                new.remove_header("Authorization")
        return new


OPENER = urllib.request.build_opener(KeepAuthorization)


def fetch_file(url, dest, attempts=None, gzip_check=True, scratch=None, checker=None, times=None, headers=None):
    """Downloads url to dest (by curl where there, else urllib): "" if it arrived whole (the length it announced, and, a
    gzip file, gzip reads it to its end), else why not. A missing file ("not found") or another refusal (HTTP 4xx but
    429) is not asked for again; other failures are (ATTEMPTS times), after growing waits, or as long as the server
    asks (Retry-After, a minute at most; urllib's alone).
    scratch: the folder of the download's .part file (eukore_build.py --scratch: dest's storage then gets the whole
    file once, written in one go), else beside dest; checker: a process pool the gzip check runs in (a core of its own),
    else this thread; times: a dict (a Counter) its seconds of transfer, check (with its wait for the checker) and
    copy (from scratch) are added to, and the bytes of the file it placed; headers: the request's others (an API
    key's: by urllib, never on curl's command line, where other users of the node could read it)."""
    why, asked = "", 0.0
    times = collections.Counter() if times is None else times
    part = os.path.join(scratch, os.path.basename(dest) + ".part") if scratch else dest + ".part"
    for attempt in range(attempts or ATTEMPTS):
        if attempt:
            time.sleep(max(RETRY_WAIT * 2 ** (attempt - 1), asked))
        if STOP.is_set():
            why = why or "stopped"
            break
        asked = 0.0
        try:
            expected, started = None, time.monotonic()
            if CURL and not headers:
                curl_get(url, part)
            else:
                request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, **(headers or {})})
                with OPENER.open(request, timeout=60) as response, open(part, "wb") as fh:
                    expected = response.headers.get("Content-Length")
                    shutil.copyfileobj(response, fh, 1 << 20)
            times["transfer"] += time.monotonic() - started
            if expected is not None and os.path.getsize(part) != int(expected):
                why = f"{os.path.getsize(part)} of {expected} bytes arrived"
                continue
            started = time.monotonic()
            if gzip_check and dest.endswith(".gz"):
                if checker:
                    checker.submit(gzip_whole, part, GZIP_TEST).result()
                else:
                    gzip_whole(part, GZIP_TEST)
            times["check"] += time.monotonic() - started
            started, size = time.monotonic(), os.path.getsize(part)
            if scratch:  # another file system: copied beside dest, then renamed into place
                shutil.copyfile(part, dest + ".part")
                os.replace(dest + ".part", dest)
                os.remove(part)
            else:
                os.replace(part, dest)
            times["copy"] += time.monotonic() - started
            times["bytes"] += size
            return ""
        except urllib.error.HTTPError as e:
            e.close()  # its response (urllib's; a ResourceWarning where left open)
            if e.code == 404:
                why = "not found"  # its .part (curl's: the error page) removed below
                break
            why = f"HTTP {e.code}"
            if 400 <= e.code < 500 and e.code != 429:
                break
            retry = e.headers.get("Retry-After", "") if e.headers else ""
            if retry.isdigit():
                asked = min(60.0, float(retry))
        except (OSError, EOFError, http.client.HTTPException, zlib.error) as e:
            why = f"{type(e).__name__}: {e}"
    for left in {part, dest + ".part"}:
        if os.path.exists(left):
            os.remove(left)
    return why


# --- progress of a long loop (since 2026-11-06: the genomes' downloads and runs) ----------------------------------

REPORT_EVERY = 600  # seconds between progress reports (eukore_build.py --report_minutes)
# a pass's first reports, seconds in: every minute of its first 10 (is the rate right?), then every REPORT_EVERY
# (since 2026-11-09; until then at 1, 2, 3, 5 and 8 minutes)
REPORT_EARLY = tuple(range(60, 601, 60))


def early_reports(every):
    """REPORT_EARLY's reports up to `every` seconds in, the later ones `every` apart (Progress's early)."""
    return [s for s in REPORT_EARLY if s <= every]


def duration(seconds):
    """Seconds as a person reads them: 45s, 12m 05s, 3h 07m, 2d 04h."""
    s = max(0, int(seconds))
    d, s = divmod(s, 86400)
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return f"{d}d {h:02d}h" if d else f"{h}h {m:02d}m" if h else f"{m}m {s:02d}s" if m else f"{s}s"


class Progress:
    """A long loop's progress, said every `every` seconds (after the item that passes it) and at its end: items done,
    failed and left, the work left, the time so far, the time left and when it should end. weights: each item's share
    of the work (a genome's size), so that the time left follows the work left, not the count of items; an item
    without one weighs the median. early: seconds in of the first reports (REPORT_EARLY), `every` seconds apart after
    the last of them. clock: time.time (a test's own). usage: a function whose text ends each report (Usage.line: the
    cores, memory and disk the build uses; since 2026-11-08). unit "": no work left said, the items left alone."""

    def __init__(self, what, weights, say, every=REPORT_EVERY, unit="Gb", scale=1e9, clock=time.time, early=(),
                 usage=None):
        known = sorted(w for w in weights if w)
        self.typical = typical = known[len(known) // 2] if known else 1
        self.weights = [w or typical for w in weights]
        self.what, self.say, self.every, self.unit, self.scale, self.clock = what, say, every, unit, scale, clock
        self.usage = usage
        self.total = sum(self.weights)
        self.n, self.done, self.failed, self.skipped, self.worked = len(self.weights), 0, 0, 0, 0.0
        self.started = self.last = clock()
        self.early = sorted(early)

    def _due(self, now):
        return now - self.started >= self.early[0] if self.early else now - self.last >= self.every

    def _report(self, now, note):
        self.last = now
        while self.early and now - self.started >= self.early[0]:
            self.early.pop(0)
        self.say(self.line(now) + note + (self.usage() if self.usage else ""))

    def tick(self, note=""):
        """No item finished, but a report is due (since 2026-11-08: none came until one finished, hours in where the
        largest genomes run first): said with note (e.g. how many are running)."""
        now = self.clock()
        if self._due(now) and self.done < self.n:
            self._report(now, note)

    def step(self, weight, failed=False, skipped=False, note=""):
        """One item finished (weight: its own; skipped: no work done, e.g. a genome not fetched yet, its work out of
        the estimate); said when `every` seconds have passed since the last report, or at the end."""
        weight = weight or self.typical  # as the item weighed in the total (until 2026-11-06 the mean: drifting)
        self.done += 1
        self.failed += bool(failed)
        if skipped:
            self.skipped += 1
            self.total -= weight
        else:
            self.worked += weight
        now = self.clock()
        if self._due(now) or self.done == self.n:
            self._report(now, note)

    def line(self, now=None):
        now = self.clock() if now is None else now
        elapsed, left_n = now - self.started, self.n - self.done
        text = f"  {self.done:,}/{self.n:,} {self.what} ({self.failed:,} failed), {left_n:,} left"
        if not left_n:
            return text + f"; took {duration(elapsed)}"
        if self.unit:
            text += f" ({max(0.0, self.total - self.worked) / self.scale:,.1f} {self.unit})"
        if self.worked and elapsed > 0:
            left = max(0.0, self.total - self.worked) / (self.worked / elapsed)
            end = datetime.datetime.fromtimestamp(now + left)
            text += f"; {duration(elapsed)} so far, ~{duration(left)} left (until ~{end:%Y-%m-%d %H:%M})"
        else:
            text += f"; {duration(elapsed)} so far"
        return text


def vsearch_global(vsearch, queries, db, ident, threads=1, accepts=32, rejects=256, options=(), tmp=None):
    """vsearch --usearch_global (--iddef 2) of queries {name: sequence} against db {name: sequence} (the same dict:
    one file for both), sequences of ident or more (a share), at most `accepts` hits a query after `rejects` misses (0:
    no limit), with further options (e.g. --self), its files in a folder in tmp (default TMPDIR) -> [(query, target,
    identity %)]; raises OSError where vsearch failed. The one vsearch search of qc.py, screen.py, intraspecific.py
    and tools/ (since 2026-11-08)."""
    import tempfile  # noqa: PLC0415
    with tempfile.TemporaryDirectory(dir=tmp) as work:
        q, d, out = (os.path.join(work, n) for n in ("q.fna", "db.fna", "hits.tsv"))
        for path, seqs in ((q, queries),) + (() if db is queries else ((d, db),)):
            with open(path, "w") as fh:
                fh.writelines(f">{n}\n{s}\n" for n, s in seqs.items() if s)
        run = subprocess.run([vsearch, "--usearch_global", q, "--db", q if db is queries else d, "--id", str(ident),
                              "--iddef", "2", "--maxaccepts", str(accepts), "--maxrejects", str(rejects), *options,
                              "--userout", out, "--userfields", "query+target+id", "--threads", str(max(1, threads)),
                              "--quiet"], capture_output=True, text=True, check=False)
        if run.returncode != 0 or not os.path.isfile(out):
            raise OSError(f"vsearch failed (exit {run.returncode}): {run.stderr.strip()[-200:]}")
        with open(out) as fh:
            return [(f[0], f[1], float(f[2])) for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) >= 3]


@contextlib.contextmanager
def process_pool(workers):
    """A pool of `workers` processes (concurrent.futures'; since 2026-11-08). A worker that dies (killed by the OOM
    killer or SLURM, a crash in a C library) breaks the pool: the build stops (SystemExit), a rerun goes on; until then
    multiprocessing.Pool lost its task and the build waited for it forever. On any other error or a stop (Ctrl-C) the
    workers are ended, not waited for."""
    pool = concurrent.futures.ProcessPoolExecutor(max(1, workers))
    pool.size = max(1, workers)
    try:
        yield pool
    except concurrent.futures.BrokenExecutor as e:
        _end_workers(pool)
        raise SystemExit(f"a worker process died ({e}): killed for memory? Rerun to go on; work done is kept") from e
    except BaseException:
        _end_workers(pool)
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)


def _end_workers(pool):
    for process in list((getattr(pool, "_processes", None) or {}).values()):
        process.terminate()


def unordered(pool, task, items, progress=None, note=lambda: "", poll=60):
    """task(item) of each item in a process_pool -> the results as they finish, at most twice the pool's workers
    submitted at a time; progress said on time between them too (Progress.tick, its note from note()), every poll
    seconds where none finished."""
    items, running = iter(items), set()
    while True:
        running |= {pool.submit(task, item) for item in itertools.islice(items, max(0, 2 * pool.size - len(running)))}
        if not running:
            return
        done, running = concurrent.futures.wait(running, timeout=poll,
                                                return_when=concurrent.futures.FIRST_COMPLETED)
        if not done and progress is not None:
            progress.tick(note())
        for future in done:
            yield future.result()


def memory_budget(environ=None, proc="/proc", cgroup="/sys/fs/cgroup"):
    """The memory this build may use (bytes) and where that came from -> (bytes, what) or (None, None): SLURM's
    allocation (SLURM_MEM_PER_NODE, or SLURM_MEM_PER_CPU times the cores; MB as SLURM counts them, MiB), else the
    cgroup's limit (v2's memory.max, v1's memory.limit_in_bytes), else the node's memory available; never more than
    the node has available."""
    env = os.environ if environ is None else environ
    got, what = None, None
    try:
        if env.get("SLURM_MEM_PER_NODE", "").isdigit() and int(env["SLURM_MEM_PER_NODE"]):
            got, what = int(env["SLURM_MEM_PER_NODE"]) * 2 ** 20, "SLURM_MEM_PER_NODE"
        elif env.get("SLURM_MEM_PER_CPU", "").isdigit():
            cpus = int(env.get("SLURM_CPUS_ON_NODE") or env.get("SLURM_CPUS_PER_TASK") or 1)
            got, what = int(env["SLURM_MEM_PER_CPU"]) * 2 ** 20 * cpus, "SLURM_MEM_PER_CPU"
    except ValueError:
        pass
    if got is None:
        for name in ("memory.max", "memory/memory.limit_in_bytes"):
            try:
                with open(os.path.join(cgroup, name)) as fh:
                    text = fh.read().strip()
            except OSError:
                continue
            if text.isdigit() and int(text) < 2 ** 60:  # v1's "none" is a huge number
                got, what = int(text), "the cgroup's limit"
                break
    try:
        with open(os.path.join(proc, "meminfo")) as fh:
            info = {line.split(":")[0]: int(line.split()[1]) * 1024 for line in fh if line.strip()}
        available = info["MemAvailable"]
    except (OSError, KeyError, ValueError, IndexError):
        available = None
    if available and (got is None or available < got):
        got, what = available, "the node's memory available"
    return got, what


def size_text(n):
    """Bytes as a person reads them: 512 MB, 84.1 GB, 2.1 TB."""
    for unit, scale in (("TB", 1e12), ("GB", 1e9)):
        if n >= scale:
            return f"{n / scale:,.1f} {unit}"
    return f"{n / 1e6:,.0f} MB"


class Usage:
    """The node's resources the build uses, for the progress reports (since 2026-11-08): the cores its processes kept
    busy since the last report (the CPU seconds of this process and every process under it, those ended included, over
    the wall's), their memory (resident, summed: pages they share counted in each) and its peak over the reports, the
    node's memory available and its load, and the free space of folders [(label, path)] (one line for those on one
    file system). From Linux's /proc; elsewhere the disk alone. cores: the build's own (-t)."""

    def __init__(self, folders=(), cores=None, clock=time.time, proc="/proc"):
        self.folders, self.cores, self.clock, self.proc = list(folders), cores, clock, proc
        self.peak = 0
        self.last = (clock(), self.tree()[0])

    def tree(self):
        """(CPU seconds, resident bytes) of this process and its descendants, or (None, None) without /proc."""
        try:
            pids = [int(e) for e in os.listdir(self.proc) if e.isdigit()]
        except OSError:
            return None, None
        children, stat = collections.defaultdict(list), {}
        for pid in pids:
            try:
                with open(os.path.join(self.proc, str(pid), "stat")) as fh:
                    text = fh.read()
            except OSError:  # ended meanwhile
                continue
            f = text[text.rfind(")") + 2:].split()  # field 3 on (the name may hold spaces and brackets)
            children[int(f[1])].append(pid)
            stat[pid] = (sum(int(x) for x in f[11:15]), int(f[21]))  # utime stime cutime cstime; rss (pages)
        cpu = rss = 0
        todo = [os.getpid()]
        while todo:
            pid = todo.pop()
            if pid in stat:
                cpu, rss = cpu + stat[pid][0], rss + stat[pid][1]
            todo += children.get(pid, [])
        return cpu / os.sysconf("SC_CLK_TCK"), rss * os.sysconf("SC_PAGE_SIZE")

    def node(self):
        """(memory available, memory in all) of the node, bytes, or None."""
        try:
            with open(os.path.join(self.proc, "meminfo")) as fh:
                info = {line.split(":")[0]: int(line.split()[1]) * 1024 for line in fh if line.strip()}
            return info["MemAvailable"], info["MemTotal"]
        except (OSError, KeyError, ValueError, IndexError):
            return None

    def disk(self):
        seen = {}
        for label, path in self.folders:
            try:
                dev, free = os.stat(path).st_dev, shutil.disk_usage(path).free
            except OSError:
                continue
            seen.setdefault(dev, [[], free])[0].append(label)
        return [(", ".join(labels), free) for labels, free in seen.values()]

    def line(self):
        """'; 30.2 cores busy (of 32), 84.1 GB resident (peak 120.3 GB); node: 300.0 GB of 512.0 GB free, load 31.5;
        disk free: WORK 2.1 TB, OUT 10.0 TB'."""
        now = self.clock()
        cpu, rss = self.tree()
        parts = []
        if cpu is not None:
            then, cpu_then = self.last
            self.last = (now, cpu)
            self.peak = max(self.peak, rss)
            busy = max(0.0, cpu - (cpu_then or 0)) / max(1e-9, now - then)
            parts.append(f"{busy:,.1f} cores busy" + (f" (of {self.cores})" if self.cores else "")
                         + f", {size_text(rss)} resident (peak {size_text(self.peak)})")
        node = self.node()
        if node:
            load = os.getloadavg()[0] if hasattr(os, "getloadavg") else None
            parts.append(f"node: {size_text(node[0])} of {size_text(node[1])} free"
                         + (f", load {load:,.1f}" if load is not None else ""))
        disk = self.disk()
        if disk:
            parts.append("disk free: " + ", ".join(f"{label} {size_text(free)}" for label, free in disk))
        return "".join("; " + p for p in parts)
