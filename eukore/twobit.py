"""UCSC's .2bit genomes (TOGA2's assemblies outside NCBI, since 2026-11-06) as gzipped FASTA, in Python alone.

The format (kent's twoBit.h): a header (signature 0x1A412743 in the file's byte order, version 0 or 1, the sequences'
count, a reserved word), an index (per sequence its name's length, its name and its record's offset: 4 bytes, 8 in
version 1), and per sequence its length, its N blocks (starts, sizes), its lower-case (soft-masked) blocks, a reserved
word and its bases, four a byte (T, C, A, G as 0-3, the first base in the high bits).
"""

import gzip
import os
import struct

SIGNATURE = 0x1A412743
BASES = b"TCAG"
LINE = 80
# a byte's four bases, one table per place in the byte
PLANES = [bytes(BASES[(b >> shift) & 3] for b in range(256)) for shift in (6, 4, 2, 0)]


def records(path):
    """A .2bit file -> [(name, sequence bytes)], one at a time: N blocks as N, masked blocks in lower case."""
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
            packed = fh.read((size + 3) // 4)
            seq = bytearray(4 * len(packed))
            for i, plane in enumerate(PLANES):
                seq[i::4] = packed.translate(plane)
            del seq[size:]
            for s, n in zip(n_starts, n_sizes):
                seq[s:s + n] = b"N" * n
            for s, n in zip(m_starts, m_sizes):
                seq[s:s + n] = seq[s:s + n].lower()
            yield name, bytes(seq)


def to_fasta(path, dest):
    """A .2bit file as dest (gzipped FASTA, LINE bases a line), written whole or not at all -> its bases."""
    total = 0
    part = dest + ".part"
    with gzip.open(part, "wb", compresslevel=6) as out:
        for name, seq in records(path):
            out.write(b">" + name.encode() + b"\n")
            out.write(b"\n".join(seq[i:i + LINE] for i in range(0, len(seq), LINE)) + b"\n")
            total += len(seq)
    os.replace(part, dest)
    return total
