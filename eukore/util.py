#!/usr/bin/env python3
"""Eukore's small helpers (since 2026-11-05; copied here so that Eukore imports nothing from protal's scripts/):
compressed files (zstd or gzip, as protal's scripts/compressed.py), a writer's temporary name, MD5s, downloads.

open_write(path) writes zstd if the name ends in .zst, gzip if in .gz (a trailing .partial does not count), plain
otherwise; open_read(path) reads a file by its first bytes (zstd, gzip or plain), as protal reads read files. zstd is
Python 3.14's compression.zstd where there is one, else the zstd command through a pipe (in the protal-db-build
environment and on Ubuntu: apt install zstd).

usage: import compressed; with compressed.open_write("reads.fq.zst") as fh: fh.write(b"...")
"""

import gzip
import shutil
import subprocess

ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
SUFFIXES = {"zstd": ".zst", "gzip": ".gz"}
LEVELS = {"zstd": 3, "gzip": 1}  # fast: the files are read once, by protal


def compression_of(path):
    """zstd, gzip or None by a file name (without a trailing .partial)."""
    name = path[:-len(".partial")] if path.endswith(".partial") else path
    return "zstd" if name.endswith(".zst") else "gzip" if name.endswith(".gz") else None


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


class _Pipe:
    """A binary file object over a zstd process: written to (compressing into a file) or read from (decompressing
    one). close() waits for the process and raises OSError if it failed."""

    def __init__(self, command, writing):
        self.writing = writing
        self.process = subprocess.Popen(command, stdin=subprocess.PIPE if writing else subprocess.DEVNULL,
                                        stdout=None if writing else subprocess.PIPE, stderr=subprocess.PIPE)
        self.stream = self.process.stdin if writing else self.process.stdout
        self.command = command

    def write(self, data):
        return self.stream.write(data)

    def read(self, size=-1):
        return self.stream.read(size)

    def readline(self, size=-1):
        return self.stream.readline(size)

    def __iter__(self):
        return iter(self.stream)

    def close(self):
        if self.process is None:
            return
        process, self.process = self.process, None
        if self.writing:
            process.stdin.close()
        else:
            process.stdout.read()  # to the end, so that the process does not wait on a full pipe
            process.stdout.close()
        error = process.stderr.read().decode(errors="replace").strip()
        process.stderr.close()
        if process.wait() != 0:
            raise OSError(f"{' '.join(self.command)} failed ({process.returncode}): {error[-300:]}")

    def __enter__(self):
        return self

    def __exit__(self, kind, value, traceback):
        if kind is None:
            self.close()
        else:  # an exception already: the process ends without hiding it
            try:
                self.close()
            except OSError:
                pass


def open_write(path, level=None):
    """A binary file object writing `path`: zstd (.zst), gzip (.gz) or plain, by its name; level: the codec's
    (default LEVELS)."""
    kind = compression_of(path)
    if kind is None:
        return open(path, "wb")
    level = LEVELS[kind] if level is None else level
    if kind == "gzip":
        return gzip.open(path, "wb", compresslevel=level)
    zstd = _stdlib_zstd()
    if zstd is not None:
        return zstd.open(path, "wb", level=level)
    return _Pipe([_zstd_command(), "-q", "-f", f"-{level}", "-o", path], writing=True)


def open_read(path):
    """A binary file object reading `path` decompressed: zstd, gzip or plain, by its first bytes."""
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"\x1f\x8b":
        return gzip.open(path, "rb")
    if is_zstd(head):
        zstd = _stdlib_zstd()
        if zstd is not None:
            return zstd.open(path, "rb")
        return _Pipe([_zstd_command(), "-dcq", path], writing=False)
    return open(path, "rb")


def read_text(path):
    """The decompressed content of a file, as text."""
    with open_read(path) as fh:
        return fh.read().decode()


# --- Eukore's own helpers (since 2026-11-05; self-contained: nothing imported from protal's scripts/) ---------------

import datetime  # noqa: E402
import hashlib  # noqa: E402
import os  # noqa: E402
import socket  # noqa: E402
import time  # noqa: E402
import urllib.error  # noqa: E402
import urllib.request  # noqa: E402
import http.client  # noqa: E402
import zlib  # noqa: E402

USER_AGENT = "Eukore-build (https://protal.earlham.ac.uk)"
ATTEMPTS = 5       # tries of a download before it is given up on
RETRY_WAIT = 2.0   # seconds before the second try, doubled for each further one (or the server's Retry-After)


def part_of(path):
    """A writer's own temporary name for `path` (renamed into place once written): .<host>-<pid>.<name>.partial beside
    it, whose compression open_write still tells by the name."""
    folder, name = os.path.split(path)
    return os.path.join(folder, f".{socket.gethostname()}-{os.getpid()}.{name}.partial")


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
        return io.StringIO(read_text(path))
    return open(path)


CURL = shutil.which("curl")  # fetch_file's transfers where there: ~3x urllib's rate (2026-11-06, NBI's HPC), else urllib


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


# fetch_file's gzip check where there: ISA-L's igzip (bioconda/conda-forge isa-l, apt isal; ~3.5x Python's zlib, 2026-11-06,
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


def fetch_file(url, dest, attempts=None, gzip_check=True, scratch=None, checker=None):
    """Downloads url to dest (by curl where there, else urllib): "" if it arrived whole (the length it announced, and, a
    gzip file, gzip reads it to its end), else why not. A missing file ("not found") or another refusal (HTTP 4xx but 429) is not asked for again; other
    failures are (ATTEMPTS times), after growing waits, or as long as the server asks (Retry-After, a minute at most).
    scratch: the folder of the download's .part file (eukore_build.py --scratch: dest's storage then gets the whole
    file once, written in one go), else beside dest; checker: a process pool the gzip check runs in (a core of its own),
    else this thread."""
    why, asked = "", 0.0
    part = os.path.join(scratch, os.path.basename(dest) + ".part") if scratch else dest + ".part"
    for attempt in range(attempts or ATTEMPTS):
        if attempt:
            time.sleep(max(RETRY_WAIT * 2 ** (attempt - 1), asked))
        asked = 0.0
        try:
            expected = None
            if CURL:
                curl_get(url, part)
            else:
                request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(request, timeout=60) as response, open(part, "wb") as fh:
                    expected = response.headers.get("Content-Length")
                    shutil.copyfileobj(response, fh, 1 << 20)
            if expected is not None and os.path.getsize(part) != int(expected):
                why = f"{os.path.getsize(part)} of {expected} bytes arrived"
                continue
            if gzip_check and dest.endswith(".gz"):
                if checker:
                    checker.submit(gzip_whole, part, GZIP_TEST).result()
                else:
                    gzip_whole(part, GZIP_TEST)
            if scratch:  # another file system: copied beside dest, then renamed into place
                shutil.copyfile(part, dest + ".part")
                os.replace(dest + ".part", dest)
                os.remove(part)
            else:
                os.replace(part, dest)
            return ""
        except urllib.error.HTTPError as e:
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
REPORT_EARLY = (60, 120, 180, 300, 480)  # the genomes' downloads' first reports, seconds in: is the rate right?


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
    the last of them. clock: time.time (a test's own)."""

    def __init__(self, what, weights, say, every=REPORT_EVERY, unit="Gb", scale=1e9, clock=time.time, early=()):
        known = sorted(w for w in weights if w)
        typical = known[len(known) // 2] if known else 1
        self.weights = [w or typical for w in weights]
        self.what, self.say, self.every, self.unit, self.scale, self.clock = what, say, every, unit, scale, clock
        self.total = sum(self.weights)
        self.n, self.done, self.failed, self.skipped, self.worked = len(self.weights), 0, 0, 0, 0.0
        self.started = self.last = clock()
        self.early = sorted(early)

    def step(self, weight, failed=False, skipped=False, note=""):
        """One item finished (weight: its own; skipped: no work done, e.g. a genome not fetched yet, its work out of
        the estimate); said when `every` seconds have passed since the last report, or at the end."""
        weight = weight or (self.total / self.n if self.n else 0)
        self.done += 1
        self.failed += bool(failed)
        if skipped:
            self.skipped += 1
            self.total -= weight
        else:
            self.worked += weight
        now = self.clock()
        due = now - self.started >= self.early[0] if self.early else now - self.last >= self.every
        if due or self.done == self.n:
            self.last = now
            while self.early and now - self.started >= self.early[0]:
                self.early.pop(0)
            self.say(self.line(now) + note)

    def line(self, now=None):
        now = self.clock() if now is None else now
        elapsed, left_n = now - self.started, self.n - self.done
        text = f"  {self.done:,}/{self.n:,} {self.what} ({self.failed:,} failed), {left_n:,} left"
        if not left_n:
            return text + f"; took {duration(elapsed)}"
        text += f" ({max(0.0, self.total - self.worked) / self.scale:,.1f} {self.unit})"
        if self.worked and elapsed > 0:
            left = max(0.0, self.total - self.worked) / (self.worked / elapsed)
            end = datetime.datetime.fromtimestamp(now + left)
            text += f"; {duration(elapsed)} so far, ~{duration(left)} left (until ~{end:%Y-%m-%d %H:%M})"
        return text
