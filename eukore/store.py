"""A done genome's results as one file (since 2026-11-10).

A genome's run leaves 8-14 small files (done.json, loci.tsv, copies.fna, features.tsv, windows.fna.gz, rrna.tsv, its
held-out tables, annotation.tsv, stats.json, organelle.fna.gz, sketch.sig.gz, fcs.tsv): for all eukaryotes some
300,000 files in 30,000 folders, on a shared file system that is slow with many small files, and the package stage
opens each of them. Packed, they are the members of one zip file beside where the folder was:

    WORK/done/<accession>/        the folder (every run until 2026-11-10; --done_folders)
    WORK/done/<accession>.zip     the same files as one (text deflated, a .gz member stored as it is)

Every function here takes the folder's path, whichever form the genome has (the zip where both are there: a pack
the kill of which left the folder behind), so that a WORK of folders is read as before and packed as the build meets
its genomes. A change (put, pack) writes the whole zip anew beside it and renames it into place: a reader sees the
old file or the new one, never a half-written one. Two changes of one genome at once (two builds on a WORK) wait
for each other (a lock file beside the zip), so that neither is lost.
"""
import collections
import contextlib
import gzip
import io
import os
import shutil
import threading
import time
import zipfile
import zlib

from . import util

EXT = ".zip"
KEEP_OPEN = 4    # zips kept open by a process: a genome's tables are read one after another (package.collect)
LOCK_STALE = 120  # a change's lock older than this is a killed process': taken over (a change takes a second)
_LOCK = threading.RLock()
_OPEN = collections.OrderedDict()  # {zip path: ((mtime ns, bytes, inode), ZipFile)}


def _after_fork():
    """A pool's worker (a fork) opens its own zips: its parent's open files share their place in the file, and a
    lock another of its parent's threads held at the fork would never be given back."""
    global _LOCK
    _LOCK = threading.RLock()
    _OPEN.clear()


if hasattr(os, "register_at_fork"):
    os.register_at_fork(after_in_child=_after_fork)


def packed(folder):
    """The zip a done folder is packed into."""
    return folder.rstrip("/\\") + EXT


def is_packed(folder):
    return os.path.isfile(packed(folder))


def exists(folder):
    """Whether the genome has results in either form."""
    return is_packed(folder) or os.path.isdir(folder)


def _forget(path):
    old = _OPEN.pop(path, None)
    if old is not None:
        old[1].close()


def _zip(folder):
    """The folder's zip, open for reading (kept open while the file is the same: its time, size and inode), or None
    where the genome is not packed. Called under _LOCK."""
    path = packed(folder)
    try:
        st = os.stat(path)
    except OSError:
        _forget(path)
        return None
    signature = (st.st_mtime_ns, st.st_size, st.st_ino)
    got = _OPEN.get(path)
    if got is not None and got[0] == signature:
        _OPEN.move_to_end(path)
        return got[1]
    _forget(path)
    try:
        opened = zipfile.ZipFile(path)
    except (OSError, zipfile.BadZipFile):
        return None
    _OPEN[path] = (signature, opened)
    while len(_OPEN) > KEEP_OPEN:
        _forget(next(iter(_OPEN)))
    return opened


def close():
    """This process' open zips closed (a test's end; a long pass' end)."""
    with _LOCK:
        for path in list(_OPEN):
            _forget(path)


def _of(folder, from_zip, from_folder, none):
    """A read of the genome in the form it has: from_zip(its open zip), else from_folder() (OSError: none). A folder
    packed between the two looks (another build's pack) is read from its zip."""
    for _again in range(2):
        with _LOCK:
            z = _zip(folder)
            if z is not None:
                return from_zip(z)
        try:
            return from_folder()
        except OSError:
            if not is_packed(folder):
                break
    return none


def names(folder):
    """The names of the genome's files, set() where it has none."""
    return _of(folder, lambda z: set(z.namelist()),
               lambda: {n for n in os.listdir(folder) if os.path.isfile(os.path.join(folder, n))}, set())


def _stat(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    return True


def has(folder, name):
    return _of(folder, lambda z: name in z.NameToInfo, lambda: _stat(os.path.join(folder, name)), False)


def _member(z, name):
    try:
        return z.read(name)
    except KeyError:
        return None


def _file(path):
    with open(path, "rb") as fh:
        return fh.read()


def read(folder, name):
    """A file's bytes, or None where the genome has none of the name."""
    return _of(folder, lambda z: _member(z, name), lambda: _file(os.path.join(folder, name)), None)


def text(folder, name):
    """A file's text (a .gz one's unpacked), or None."""
    return decoded(name, read(folder, name))


def decoded(name, data):
    """A file's bytes (read's, members') as its text (a .gz one's unpacked), or None for None."""
    if data is None:
        return None
    return (gzip.decompress(data) if name.endswith(".gz") else data).decode()


def open_text(folder, name):
    """A file's text as a file to read lines from, or None."""
    got = text(folder, name)
    return None if got is None else io.StringIO(got)


def rows(folder, name):
    """A table of the genome (header first) -> its rows, split; [] where it has none."""
    fh = open_text(folder, name)
    if fh is None:
        return []
    fh.readline()
    return [line.rstrip("\n").split("\t") for line in fh if line.strip()]


def table(folder, name):
    """A table of the genome (header first) -> [row dicts]; [] where it has none."""
    return table_of(text(folder, name))


def table_of(got):
    """A table's text (header first; text's, decoded's) -> [row dicts]; [] for None."""
    if got is None:
        return []
    fh = io.StringIO(got)
    head = fh.readline().rstrip("\n").split("\t")
    return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def extract(folder, name, dest):
    """A file of the genome written as dest (a tool that reads a path) -> whether it has one."""
    data = read(folder, name)
    if data is None:
        return False
    with open(dest, "wb") as fh:
        fh.write(data)
    return True


def stamp(folder, name="done.json"):
    """What tells a file of the genome from one written since: its checksum (CRC-32, the same in a folder and packed:
    a pack is no new result), or None where it has none."""
    def from_zip(z):
        info = z.NameToInfo.get(name)
        return None if info is None else info.CRC
    return _of(folder, from_zip, lambda: zlib.crc32(_file(os.path.join(folder, name))), None)


def snapshot(folder, name="done.json"):
    """The genome as one look sees it -> (its files' names, name's bytes, their stamp), or None where it has no
    results. Read without the process' lock and its open zips (since 2026-11-11), so that many threads read many
    genomes at once: a plan's reads of some 30,000 genomes wait on the shared file system, one at a time half an hour
    with 32 builds reading together."""
    for _again in range(2):
        try:
            with zipfile.ZipFile(packed(folder)) as z:
                info = z.NameToInfo.get(name)
                return set(z.namelist()), None if info is None else z.read(info), None if info is None else info.CRC
        except FileNotFoundError:
            pass
        except (OSError, zipfile.BadZipFile):
            return None
        try:
            files = {n for n in os.listdir(folder) if os.path.isfile(os.path.join(folder, n))}
            data = _file(os.path.join(folder, name)) if name in files else None
            return files, data, None if data is None else zlib.crc32(data)
        except OSError:
            if not is_packed(folder):  # neither form; else packed between the two looks: its zip read
                return None
    return None


def members(folder, wanted):
    """Several files of the genome in one look, read as snapshot reads (without the process' lock and its open zips:
    threads read many genomes at once; since 2026-11-11 the package stage's passes over every genome, until then one
    genome after another) -> {name: its bytes, or None where the genome has none of the name}, or None where it has no
    results."""
    for _again in range(2):
        try:
            with zipfile.ZipFile(packed(folder)) as z:
                return {n: z.read(z.NameToInfo[n]) if n in z.NameToInfo else None for n in wanted}
        except FileNotFoundError:
            pass
        except (OSError, zipfile.BadZipFile):
            return None
        try:
            files = {n for n in os.listdir(folder) if os.path.isfile(os.path.join(folder, n))}
            return {n: _file(os.path.join(folder, n)) if n in files else None for n in wanted}
        except OSError:
            if not is_packed(folder):  # neither form; else packed between the two looks: its zip read
                return None
    return None


def _method(name):
    return zipfile.ZIP_STORED if name.endswith((".gz", ".zst", ".zip")) else zipfile.ZIP_DEFLATED


@contextlib.contextmanager
def _changing(folder):
    """The genome's results held for a change by one thread of one process: a lock file beside its zip (made with
    O_EXCL; one older than LOCK_STALE is a killed process', taken over). Other genomes' changes go on beside it."""
    lock = packed(folder) + ".lock"
    os.makedirs(os.path.dirname(lock) or ".", exist_ok=True)
    while True:
        try:
            os.close(os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644))
            break
        except FileExistsError:
            with contextlib.suppress(OSError):
                if time.time() - os.path.getmtime(lock) > LOCK_STALE:
                    os.remove(lock)
                    continue
            time.sleep(0.05)
    try:
        yield
    finally:
        with contextlib.suppress(OSError):
            os.remove(lock)


def _write(path, members, old=None, drop=()):
    """A zip written whole as path (its partial name, then renamed): members {name: bytes, or a file's path}, after
    those of old (an open zip) that are neither in members nor in drop, kept as they are."""
    part = util.part_of(path)
    try:
        with zipfile.ZipFile(part, "w", compresslevel=6) as out:
            if old is not None:
                for info in old.infolist():
                    if info.filename not in members and info.filename not in drop:
                        out.writestr(info, old.read(info.filename), compress_type=info.compress_type)
            for name in sorted(members):
                if isinstance(members[name], bytes):
                    out.writestr(name, members[name], compress_type=_method(name))
                else:
                    out.write(members[name], name, compress_type=_method(name))
        os.replace(part, path)
    except BaseException:
        if os.path.isfile(part):
            os.remove(part)
        raise


def pack(source, folder=None):
    """A folder of results (source: a run's partial folder, or the done folder itself) packed as the zip of folder
    (default: source), the done folder of the same genome and source removed after -> the zip's path. Nothing is
    removed before the zip is in place: a kill leaves the genome as it was, or packed with its folder still there
    (the zip is read, the folder removed by the next pack or put). The done folder of a genome packed already (that
    left-behind folder, or another build's pack meanwhile) is removed, its zip kept: the zip is the one changed
    since. The folder is read under the change's lock: a put meanwhile is in the zip whole or not at all."""
    folder = folder or source
    with _changing(folder):
        if source == folder and is_packed(folder):
            shutil.rmtree(folder, ignore_errors=True)
            return packed(folder)
        members = {}
        for top, _dirs, files in os.walk(source):
            for f in files:
                path = os.path.join(top, f)
                if not f.endswith(".partial"):  # a killed change's half-written file
                    members[os.path.relpath(path, source).replace(os.sep, "/")] = path
        if not members:  # no such folder (another build's pack took it)
            return packed(folder) if is_packed(folder) else None
        _write(packed(folder), members)
        for old in {source, folder}:
            shutil.rmtree(old, ignore_errors=True)
    return packed(folder)


def unpacked(done):
    """The done folders under done (WORK/done) not packed yet, or left behind by a killed pack: those with a
    done.json, no partial's or claim's name."""
    try:
        with os.scandir(done) as entries:
            found = [e.path for e in entries if not e.name.startswith(".") and e.is_dir()]
    except OSError:
        return []
    return sorted(p for p in found if os.path.isfile(os.path.join(p, "done.json")))


def put(folder, files):
    """The genome's files changed: files {name: its bytes, a file's path (read; left where it is), or None to remove
    it}, in the form the genome has (a zip written anew and renamed into place; in a folder each file written under
    its partial name and renamed). A genome with nothing yet gets a folder."""
    with _changing(folder):
        with _LOCK:  # the open zip is this process' threads' shared one
            z = _zip(folder)
            if z is not None:
                _write(packed(folder), {n: v for n, v in files.items() if v is not None}, z,
                       {n for n, v in files.items() if v is None})
        if z is not None:
            if os.path.isdir(folder):  # a killed pack's
                shutil.rmtree(folder, ignore_errors=True)
            return
        os.makedirs(folder, exist_ok=True)
        for name, value in files.items():
            path = os.path.join(folder, name)
            if value is None:
                if os.path.isfile(path):
                    os.remove(path)
                continue
            if isinstance(value, bytes):
                with open(util.part_of(path), "wb") as fh:
                    fh.write(value)
            else:
                shutil.copyfile(value, util.part_of(path))
            os.replace(util.part_of(path), path)
