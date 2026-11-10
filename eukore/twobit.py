"""UCSC's .2bit genomes (TOGA2's assemblies outside NCBI, since 2026-11-06) as gzipped FASTA, in Python alone.

The format (kent's twoBit.h): a header (signature 0x1A412743 in the file's byte order, version 0 or 1, the sequences'
count, a reserved word), an index (per sequence its name's length, its name and its record's offset: 4 bytes, 8 in
version 1), and per sequence its length, its N blocks (starts, sizes), its lower-case (soft-masked) blocks, a reserved
word and its bases, four a byte (T, C, A, G as 0-3, the first base in the high bits).

Memory: a sequence is held once (its bytearray), unpacked and written out in pieces, ~1.1 bytes a base at the peak; until
2026-11-08 a bytes copy, a list of its lines and the joined record made it ~5.8, and the fetch stage converts
--downloads genomes at once: 12 mammalian chromosomes' ~17 GB killed a build on a 10 GB node.
"""

import os
import struct

from . import util

SIGNATURE = 0x1A412743
BASES = b"TCAG"
LINE = 80
PIECE = 1 << 22  # bases unpacked at a time (a multiple of 4)
# a byte's four bases, one table per place in the byte
PLANES = [bytes(BASES[(b >> shift) & 3] for b in range(256)) for shift in (6, 4, 2, 0)]


def records(path):
    """A .2bit file -> [(name, sequence bytearray)], one at a time: N blocks as N, masked blocks in lower case. The
    bytearray is the one copy of the sequence (no bytes copy), let go of before the next record's is made."""
    with open(path, "rb") as fh:
        head = fh.read(16)
        for order in "<>":
            sig, version, count, _r = struct.unpack(order + "IIII", head)
            if sig == SIGNATURE:
                break
        else:
            raise ValueError(f"{path}: not a .2bit file")
        index = []
        for _i in range(count):
            n = fh.read(1)[0]
            name = fh.read(n).decode()
            offset = struct.unpack(order + ("Q" if version == 1 else "I"), fh.read(8 if version == 1 else 4))[0]
            index.append((name, offset))
        for name, offset in index:
            fh.seek(offset)
            size, blocks = struct.unpack(order + "II", fh.read(8))
            n_starts = struct.unpack(f"{order}{blocks}I", fh.read(4 * blocks))
            n_sizes = struct.unpack(f"{order}{blocks}I", fh.read(4 * blocks))
            masks = struct.unpack(order + "I", fh.read(4))[0]
            m_starts = struct.unpack(f"{order}{masks}I", fh.read(4 * masks))
            m_sizes = struct.unpack(f"{order}{masks}I", fh.read(4 * masks))
            fh.read(4)
            seq = bytearray(4 * ((size + 3) // 4))
            # since 2026-11-10 named before the loop: a record of length 0 never enters it, and `del packed` below
            # failed on the unbound name (UnboundLocalError), the whole genome's conversion with it
            packed = b""
            for at in range(0, len(seq), PIECE):  # the packed bases read and unpacked a piece at a time: no copy
                packed = fh.read(min(PIECE, len(seq) - at) // 4)
                for i, plane in enumerate(PLANES):
                    seq[at + i:at + 4 * len(packed):4] = packed.translate(plane)
            del packed, seq[size:]
            for s, n in zip(n_starts, n_sizes):
                seq[s:s + n] = b"N" * n
            for s, n in zip(m_starts, m_sizes):
                seq[s:s + n] = seq[s:s + n].lower()
            yield name, seq
            del seq


def to_fasta(path, dest, chunk=LINE * 16384):
    """A .2bit file as dest (gzipped FASTA, LINE bases a line), written whole or not at all -> its bases. chunk: the
    bases lined and written at a time (a multiple of LINE), the sequence itself never copied."""
    total = 0
    # this process' own partial name (since 2026-11-10; until then dest.part, which two runs converting one genome
    # both wrote); it ends in .part still, so the fetch stage removes a killed conversion's
    part = util.fetch_part(dest)
    with util.bgzf_open(part, "wb", util.WORK_LEVEL) as out:
        for name, seq in records(path):
            out.write(b">" + name.encode() + b"\n")
            view = memoryview(seq)
            for i in range(0, len(seq), chunk):
                block = view[i:i + chunk]
                out.write(b"\n".join(block[j:j + LINE] for j in range(0, len(block), LINE)) + b"\n")
                block.release()
            view.release()
            total += len(seq)
            del seq  # before the next record's sequence is made: one in memory at a time
    os.replace(part, dest)
    return total
