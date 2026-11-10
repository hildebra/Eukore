#!/usr/bin/env python3
"""Eukore's small helpers (since 2026-11-05; copied here so that Eukore imports nothing from protal's scripts/): text
files plain, gzipped or zstd-compressed (open_text; zstd by Python 3.14's compression.zstd where there is one, else the
zstd command through a pipe), a writer's temporary name, MD5s, downloads, a long loop's progress; since 2026-11-09
every compressed file Eukore writes as BGZF (bgzf_open), on many threads where asked (compression_threads).
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
import signal  # noqa: E402
import socket  # noqa: E402
import time  # noqa: E402
import urllib.error  # noqa: E402
import urllib.parse  # noqa: E402
import urllib.request  # noqa: E402
import http.client  # noqa: E402
import io  # noqa: E402
import struct  # noqa: E402
import zlib  # noqa: E402

import threading  # noqa: E402
import weakref  # noqa: E402

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


def dead_partials(folder, name, end=".partial"):
    """The partials (part_of's .<host>-<pid>.<name>.partial) in folder of runs that ended: this host's whose process is
    gone, another host's untouched for PARTIAL_STALE (since 2026-11-08: a run deleted every one, another build's going
    on too) -> their paths. end: a partial's own ending where a tool reads its file's format off the name
    (metadata.describe's .<host>-<pid>.sketch.partial.sig.gz)."""
    out, host = [], socket.gethostname()
    for path in glob.glob(os.path.join(glob.escape(folder), f".*.{glob.escape(name)}{glob.escape(end)}")):
        their_host, _dash, pid = os.path.basename(path)[1:-len(f".{name}{end}")].rpartition("-")
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


def publish(part, dest):
    """A file written under a temporary name (part) renamed into place (dest), its bytes on the disk first and the
    folder's entry after (fsync; since 2026-11-10, for the files a release publishes: a node that went down after the
    rename could leave a file of the final name empty or cut on a network file system, which a rerun then kept)."""
    fd = os.open(part, os.O_RDWR)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)
    os.replace(part, dest)
    with contextlib.suppress(OSError):  # a folder is not opened or synced everywhere (Windows, some network systems)
        fd = os.open(os.path.dirname(os.path.abspath(dest)), os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
    return dest


# --- one run at a time (since 2026-11-10) ----------------------------------------------------------------------------
# A stage that two runs must not go through at once on one folder (the fetch stage: its downloads have no claims of
# their own) holds a lock file there, as a genome's run holds its claim (genomes.claim, genomes.claimed).

RUN_BEAT = 60    # seconds between a held lock's touches (single_run)
RUN_STALE = 300  # another host's lock untouched this long is a dead run's
RUN_YOUNG = 60   # a lock that cannot be read (written this instant) is held while younger; a takeover's lock too
_RUNS = set()    # the locks this process holds (one whose file names this process but is not here: a dead run's)
_RUNS_LOCK = threading.Lock()


class RunHeld(RuntimeError):
    """Another run holds a single_run lock: its path, and the host and process it names (None where unreadable)."""

    def __init__(self, path, host=None, pid=None):
        who = f"process {pid} on {host}" if host else "a run that is starting"
        super().__init__(f"{path}: held by {who}; one run at a time (a dead run's lock is taken over: this host's "
                         f"when its process is gone, another host's {RUN_STALE} s after its last touch)")
        self.path, self.host, self.pid = path, host, pid


def run_lock(folder, name):
    """The path of single_run's lock file: folder/.<name>.lock."""
    return os.path.join(folder, f".{name}.lock")


def lock_holder(path):
    """(host, pid) a lock file names, or None (not there, or not readable)."""
    try:
        with open(path) as fh:
            host, pid = fh.read().split()
        return host, int(pid)
    except (OSError, ValueError):
        return None


def _lock_mine():
    return socket.gethostname(), os.getpid()


def lock_live(path, stale=RUN_STALE):
    """Whether a lock file is held: by a process alive on this host, or by another host's that touched it within
    `stale` seconds (single_run touches its lock every RUN_BEAT seconds)."""
    try:
        age = time.time() - os.path.getmtime(path)
    except OSError:  # not there
        return False
    holder = lock_holder(path)
    if holder is None:  # being written this instant, or damaged: held while young
        return age < RUN_YOUNG
    if holder == _lock_mine():  # this process' own; a file that names it but is not held here: an earlier process'
        return path in _RUNS
    if holder[0] == socket.gethostname():
        return pid_alive(holder[1])
    return age < stale


def _lock_new(path):
    """This process' lock file made where none is there (O_EXCL: one process alone can) -> whether it was."""
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as fh:
        fh.write("{}\t{}\n".format(*_lock_mine()))
    return True


def _lock_taken(path, stale):
    """The lock file made for this process -> whether it was. A dead run's (lock_live) is taken over under a lock of
    the takeover (<path>.take, O_EXCL), so that two runs that both found it dead do not each remove what is there and
    make their own, the second removing the first's new lock (as genomes.claim)."""
    if _lock_new(path):
        return True
    if lock_live(path, stale):
        return False
    take = path + ".take"
    try:
        os.close(os.open(take, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644))
    except FileExistsError:  # another run takes it over this instant; a killed takeover's lock removed for the next
        with contextlib.suppress(OSError):
            if time.time() - os.path.getmtime(take) > RUN_YOUNG:
                os.remove(take)
        return False
    try:
        if lock_live(path, stale):  # taken over since it was found dead
            return False
        with contextlib.suppress(FileNotFoundError):
            os.remove(path)
        return _lock_new(path)
    finally:
        with contextlib.suppress(OSError):
            os.remove(take)


@contextlib.contextmanager
def single_run(folder, name, stale=RUN_STALE, beat=RUN_BEAT):
    """One run at a time through the block, of all runs on folder, this host's and other hosts' (since 2026-11-10): a
    lock file folder/.<name>.lock (run_lock) naming this host and process, made with O_CREAT|O_EXCL, touched every
    `beat` seconds by a thread while the block runs and removed after it, while it is this process' still -> yields
    the lock's path. Raises RunHeld where another run holds it: a process alive on this host, or another host's whose
    lock was touched within `stale` seconds; a dead run's lock is taken over.
        try:
            with util.single_run(opts.genomes, "fetch"):
                ...  # the stage
        except util.RunHeld as e:
            raise SystemExit(f"fetch: {e}") from None
    """
    os.makedirs(folder, exist_ok=True)
    path = run_lock(folder, name)
    with _RUNS_LOCK:  # this process' threads one at a time: the file names the process, _RUNS says whether it holds it
        if not _lock_taken(path, stale):
            raise RunHeld(path, *(lock_holder(path) or ()))
        _RUNS.add(path)
    stop = threading.Event()

    def touch():
        while not stop.wait(beat):
            if lock_holder(path) == _lock_mine():  # not another run's, which took over a run stalled past `stale`
                with contextlib.suppress(OSError):
                    os.utime(path)
    heart = threading.Thread(target=touch, daemon=True)
    heart.start()
    try:
        yield path
    finally:
        stop.set()
        heart.join()
        _RUNS.discard(path)
        if lock_holder(path) == _lock_mine():
            with contextlib.suppress(OSError):
                os.remove(path)


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


# --- BGZF: every file Eukore compresses (since 2026-11-09; gzip's own stream until then) -----------------------------
# BGZF (htslib's blocked gzip, SAM/BAM specification 4.1) is gzip, one member per 64 KiB block, each block's size in
# its header: gzip, zcat, Python's gzip and every zlib reader read it as they read gzip; bgzip -@, samtools, seqkit and
# rapidgzip decompress it on many threads, and samtools faidx and tabix index it. Its blocks are compressed apart,
# so a writer compresses them on many threads too (compression_threads).

BGZF_LEVEL = 9  # the package's and the release's files, those a user downloads: deflate's best (since 2026-11-09)
WORK_LEVEL = 6  # the work folder's own files, read by the build alone: faster
BGZF_BLOCK = 0xff00  # uncompressed bytes a block (htslib's BGZF_BLOCK_SIZE: its deflate always fits the 64 KiB a
#                      block may take, header and trailer included)
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")  # the empty last block
BGZF_PENDING = 32  # a threaded writer's blocks compressing at once at most (32 x 64 KiB a writer)
_BGZF_POOL = {"threads": 1, "pid": None, "executor": None}
_BGZF_LOCK = [threading.Lock()]  # the executor made once where several threads open their first writers at once
_BGZF_OPEN = weakref.WeakSet()  # this process' writers on its compression threads, open (bgzf_quiet)


def _bgzf_lock_anew():
    _BGZF_LOCK[0] = threading.Lock()  # a forked child's own: the parent's may be held by a thread the child has not


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_bgzf_lock_anew)


def compression_threads(n):
    """The threads every BGZF writer of this process compresses its blocks on from now (bgzf_open; 1: in the writer's
    own thread). zlib lets go of Python's lock while it compresses, so threads compress in parallel; the blocks are
    written in their order whatever thread made them, the file the same."""
    old = _BGZF_POOL["executor"]
    if old is not None and _BGZF_POOL["pid"] == os.getpid():
        old.shutdown(wait=True)
    _BGZF_POOL.update(threads=max(1, int(n or 1)), pid=None, executor=None)


def _bgzf_executor():
    """The process' compression threads (compression_threads), None for none; made at first use, again in a child
    process (a forked child has the parent's executor without its threads)."""
    if _BGZF_POOL["threads"] <= 1:
        return None
    with _BGZF_LOCK[0]:  # since 2026-11-10: the release's threads (release.downloads) each made one at their first file
        if _BGZF_POOL["pid"] != os.getpid():
            _BGZF_POOL.update(pid=os.getpid(), executor=concurrent.futures.ThreadPoolExecutor(
                _BGZF_POOL["threads"], thread_name_prefix="bgzf"))
        return _BGZF_POOL["executor"]


def bgzf_quiet():
    """This process' compression threads ended where no BGZF writer of it is open (since 2026-11-10; made again at the
    next writer's first block) -> whether they are: a process forked while threads run may hang (a lock a thread
    held, never let go of in the child; Python 3.12 warns of it), and they stayed, idle, once made."""
    if any(not w.closed for w in _BGZF_OPEN):
        return False
    old = _BGZF_POOL["executor"]
    if old is not None and _BGZF_POOL["pid"] == os.getpid():
        old.shutdown(wait=True)
    _BGZF_POOL.update(pid=None, executor=None)
    return True


def bgzf_block(data, level=BGZF_LEVEL):
    """One BGZF block of data (at most BGZF_BLOCK bytes): a gzip member with a BC extra field giving its size."""
    deflate = zlib.compressobj(level, zlib.DEFLATED, -15)
    body = deflate.compress(data) + deflate.flush()
    # ID1 ID2 CM FLG(FEXTRA) MTIME XFL OS(unknown) XLEN, subfield B C SLEN BSIZE (the block's size less 1)
    return (struct.pack("<4BIBBH2BHH", 31, 139, 8, 4, 0, 0, 255, 6, 66, 67, 2, len(body) + 25) + body
            + struct.pack("<II", zlib.crc32(data), len(data)))


class BgzfWriter(io.RawIOBase):
    """A BGZF file being written (bytes): its data cut into BGZF_BLOCK blocks, compressed at level on the process'
    compression threads (compression_threads) or in this thread, closed by the empty EOF block. target: a path or a
    binary file object (left open)."""

    def __init__(self, target, level=BGZF_LEVEL):
        super().__init__()
        self._own = isinstance(target, (str, bytes, os.PathLike))
        self._fh = open(target, "wb") if self._own else target
        self._level = level
        self._buf = bytearray()
        self._pending = collections.deque()
        self._pool = _bgzf_executor()
        if self._pool is not None:
            _BGZF_OPEN.add(self)

    def writable(self):
        return True

    def write(self, data):
        self._buf += data
        if len(self._buf) >= BGZF_BLOCK:
            view, at = memoryview(self._buf), 0
            while len(self._buf) - at >= BGZF_BLOCK:
                self._block(bytes(view[at:at + BGZF_BLOCK]))
                at += BGZF_BLOCK
            view.release()
            del self._buf[:at]
        return len(data)

    def _block(self, data):
        if self._pool is None:
            self._fh.write(bgzf_block(data, self._level))
            return
        self._pending.append(self._pool.submit(bgzf_block, data, self._level))
        while self._pending and (len(self._pending) > BGZF_PENDING or self._pending[0].done()):
            self._fh.write(self._pending.popleft().result())

    def close(self):
        if self.closed:
            return
        try:
            if self._buf:
                self._block(bytes(self._buf))
                self._buf.clear()
            while self._pending:
                self._fh.write(self._pending.popleft().result())
            self._fh.write(BGZF_EOF)
        finally:
            if self._own:
                self._fh.close()
            super().close()


def bgzf_open(target, mode="wt", level=BGZF_LEVEL, encoding="utf-8", newline="\n"):
    """A BGZF file opened for writing (BgzfWriter): mode "wb" (bytes) or "wt" (text, newline "\\n" by default); level:
    BGZF_LEVEL for a file of the package or the release, WORK_LEVEL for the work folder's own. Read it back with
    gzip.open or open_text."""
    if mode not in ("wb", "wt"):
        raise ValueError(f"bgzf_open: mode {mode!r} (wb or wt)")
    binary = io.BufferedWriter(BgzfWriter(target, level), buffer_size=1 << 20)
    return binary if mode == "wb" else io.TextIOWrapper(binary, encoding=encoding, newline=newline)


def is_bgzf(path):
    """Whether a file starts with a BGZF block (a gzip header whose first extra subfield is BC)."""
    with open(path, "rb") as fh:
        head = fh.read(18)
    return len(head) == 18 and head[:4] == b"\x1f\x8b\x08\x04" and head[12:14] == b"BC"


def bgzf_whole(path):
    """Whether a file ends with BGZF's empty end block (BGZF_EOF), as every file bgzf_open wrote to its end does: one
    cut short (a killed writer, a node that went down, a copy that broke off) does not, nor does a file that is not
    there (since 2026-11-10; nothing looked at a kept file's end)."""
    try:
        with open(path, "rb") as fh:
            if fh.seek(0, os.SEEK_END) < len(BGZF_EOF):
                return False
            fh.seek(-len(BGZF_EOF), os.SEEK_END)
            return fh.read(len(BGZF_EOF)) == BGZF_EOF
    except OSError:
        return False


def bgzf_copy(src, dest, level=BGZF_LEVEL):
    """src (plain, gzipped or zstd-compressed) as a BGZF file dest, written beside it and renamed into place."""
    part = part_of(dest)
    with open_binary(src) as fh, bgzf_open(part, "wb", level) as out:
        shutil.copyfileobj(fh, out, 1 << 22)
    os.replace(part, dest)
    return dest


def open_binary(path):
    """A file's bytes, plain or decompressed (gzip and BGZF; zstd read whole), by its first bytes."""
    with open(path, "rb") as fh:
        head = fh.read(4)
    if head[:2] == b"\x1f\x8b":
        return gzip.open(path, "rb")
    if is_zstd(head):
        return io.BytesIO(read_zstd(path).encode())
    return open(path, "rb")


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


def fetch_part(path):
    """A download's own temporary name, <path>.<host>-<pid>.part (since 2026-11-10; until then <path>.part, which two
    runs fetching one file both wrote and renamed: one's bytes in the other's file). No download goes on from a part
    an earlier one left (curl -o and urllib write it anew), so nothing needs one name. It ends in .part still: the
    fetch stage removes the old *.part files of killed runs."""
    return f"{path}.{socket.gethostname()}-{os.getpid()}.part"


def fetch_file(url, dest, attempts=None, gzip_check=True, scratch=None, checker=None, times=None, headers=None):
    """Downloads url to dest (by curl where there, else urllib): "" if it arrived whole (the length it announced, and, a
    gzip file, gzip reads it to its end), else why not. A missing file ("not found") or another refusal (HTTP 4xx but
    429) is not asked for again; other failures are (ATTEMPTS times), after growing waits, or as long as the server
    asks (Retry-After, a minute at most; urllib's alone).
    scratch: the folder of the download's .part file (fetch_part: this process' own; eukore_build.py --local_tmp:
    dest's storage then gets the whole file once, written in one go), else beside dest; checker: a process pool the
    gzip check runs in (a core of its own), else this thread; times: a dict (a Counter) its seconds of transfer, check
    (with its wait for the checker) and copy (from scratch) are added to, and the bytes of the file it placed;
    headers: the request's others (an API key's: by urllib, never on curl's command line, where other users of the
    node could read it)."""
    why, asked = "", 0.0
    times = collections.Counter() if times is None else times
    beside = fetch_part(dest)
    part = fetch_part(os.path.join(scratch, os.path.basename(dest))) if scratch else beside
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
                shutil.copyfile(part, beside)
                os.replace(beside, dest)
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
    for left in {part, beside}:
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
    workers are ended, not waited for, and the tools they started with them (_end_workers). This process'
    idle compression threads ended first (bgzf_quiet; since 2026-11-10)."""
    bgzf_quiet()
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


@contextlib.contextmanager
def forked_pool(workers):
    """A process_pool of `workers` whose workers are all forked on entering (since 2026-11-10), None for one worker:
    the compression threads ended first (process_pool), so that a caller opening BGZF writers inside it forks
    no process while their threads run. The parent's objects (the package's tens of GB) out of the collector's reach
    while workers fork from it: a collection in a worker would touch every one of them, and copy their pages into
    each worker (Python's gc docs)."""
    if workers <= 1:
        yield None
        return
    import gc  # noqa: PLC0415
    gc.freeze()
    try:
        with process_pool(workers) as pool:
            pool.submit(os.getpid).result()  # forks every worker where processes fork (Linux until Python 3.14)
            yield pool
    finally:
        gc.unfreeze()


def ordered(task, items, workers=1, ahead=None, pool=None):
    """task(item) of each item -> the results in the items' order (since 2026-11-10: the package's reading of every
    done genome, until then one after another): in this process where workers is 1, else in a forked_pool of
    `workers` (or `pool`, one entered before: its size), at most `ahead` (default 4 x workers) submitted beyond the
    result awaited, so that the parent holds a few items' results, not all. task: a module's function (a worker
    imports it); an item's error is raised here."""
    if pool is None and workers <= 1:
        for item in items:
            yield task(item)
        return
    with contextlib.nullcontext(pool) if pool is not None else forked_pool(workers) as pool:
        ahead = max(1, ahead or 4 * pool.size)
        pending = collections.deque()
        for item in items:
            pending.append(pool.submit(task, item))
            if len(pending) > ahead:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()


def descendants(pids, proc="/proc"):
    """The processes under pids (their children, and theirs), by /proc's parents -> their pids; [] without /proc."""
    try:
        entries = [int(e) for e in os.listdir(proc) if e.isdigit()]
    except OSError:
        return []
    children = collections.defaultdict(list)
    for pid in entries:
        try:
            with open(os.path.join(proc, str(pid), "stat")) as fh:
                text = fh.read()
            children[int(text[text.rfind(")") + 2:].split()[1])].append(pid)  # field 4: the parent
        except (OSError, ValueError, IndexError):  # ended meanwhile
            continue
    out, todo = [], list(pids)
    while todo:
        below = children.get(todo.pop(), [])
        out += below
        todo += below
    return out


def _end_workers(pool):
    """The pool's workers ended (SIGTERM), and the processes under them: the tools they run (miniprot, vsearch;
    since 2026-11-10: the workers alone were ended, and on a workstation their tools went on to their end, writing
    into folders the next run uses; SLURM ends a job's with its cgroup). The tools are listed before the workers end:
    after it they are init's children. Where there is no /proc (macOS) the workers alone, as before."""
    workers = list((getattr(pool, "_processes", None) or {}).values())
    below = descendants([p.pid for p in workers if p.pid])
    for process in workers:
        process.terminate()
    for pid in below:
        with contextlib.suppress(OSError):  # ended meanwhile
            os.kill(pid, signal.SIGTERM)


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
