import gzip
import struct
import logging
import numpy as np


QS2_MAGIC = b"\x0b\x0e\x0a\xc1"
QDATA_MAGIC = b"\x0b\x0e\x0a\xcd"
QS_LEGACY_MAGIC = b"\x0b\x0e\x0a\x0c"
QS2_MAX_BLOCKSIZE = 1048576
QS2_SHUFFLE_MASK = 1 << 31
QS2_SHUFFLE_ELEMSIZE = 8


class ByteStream:
    """Wraps a binary file object and keeps the byte order of the serialization stream ('>' for XDR, '<'/'>' for native binary)."""

    def __init__(self, f):
        self.f = f
        self.read = f.read
        self.tell = f.tell
        self.bo = ">"

    def close(self):
        self.f.close()


def _zstd_decompressor():
    try:
        from compression import zstd  # python >= 3.14

        return zstd.decompress
    except ImportError:
        pass
    try:
        import zstandard
    except ImportError:
        raise ImportError(
            "reading qs2 files requires the 'zstandard' package (pip install zstandard)"
        )
    dctx = zstandard.ZstdDecompressor()
    return lambda data: dctx.decompress(data, max_output_size=QS2_MAX_BLOCKSIZE)


class Qs2Reader:
    """File-like reader of qs2 files (qs2::qs_save), see https://github.com/qsbase/qs2.
    File is 24 bytes header followed by zstd compressed blocks (uint32 compressed size, highest bit marks byte-shuffled blocks)
    that together contain a standard R serialization stream (native binary format).
    """

    def __init__(self, file_path):
        self.f = open(file_path, "rb")
        header = self.f.read(24)
        if len(header) < 24 or header[:4] != QS2_MAGIC:
            self.f.close()
            raise ValueError(f"{file_path} is not a qs2 file")
        format_version, compression, endian, shuffle = header[4:8]
        logging.debug(
            f"qs2 format version {format_version}; compression {compression}; endian {endian}; shuffle {shuffle}"
        )
        if format_version > 1:
            logging.warning(
                f"qs2 format version {format_version} is newer than supported (1)"
            )
        if compression != 1:
            self.f.close()
            raise NotImplementedError(
                f"Unknown qs2 compression algorithm '{compression}'"
            )
        # blocks sizes are written in native order of the machine that saved the file
        self.file_bo = ">" if endian == 1 else "<"
        self.bo = ">"
        self._decompress = _zstd_decompressor()
        self._buf = b""
        self._pos = 0
        self._offset = 0  # decompressed bytes preceding current block

    def _next_block(self):
        zsize_bytes = self.f.read(4)
        if len(zsize_bytes) < 4:
            return False
        zsize = struct.unpack(self.file_bo + "I", zsize_bytes)[0]
        shuffled = bool(zsize & QS2_SHUFFLE_MASK)
        zsize &= ~QS2_SHUFFLE_MASK
        max_zsize = QS2_MAX_BLOCKSIZE + QS2_MAX_BLOCKSIZE // 256
        if not 0 < zsize <= max_zsize:
            raise ValueError(f"Invalid qs2 compressed block size: {zsize}")
        zblock = self.f.read(zsize)
        if len(zblock) != zsize:
            raise EOFError("Unexpected end of qs2 file while reading block")
        block = self._decompress(zblock)
        if shuffled:
            # undo blosc byte-shuffle (element size 8), trailing bytes are not shuffled
            n = len(block) - len(block) % QS2_SHUFFLE_ELEMSIZE
            unshuffled = (
                np.frombuffer(block, dtype=np.uint8, count=n)
                .reshape(QS2_SHUFFLE_ELEMSIZE, -1)
                .T.tobytes()
            )
            block = unshuffled + block[n:]
        self._offset += len(self._buf)
        self._buf = block
        self._pos = 0
        return True

    def read(self, n):
        end = self._pos + n
        if end <= len(self._buf):
            out = self._buf[self._pos : end]
            self._pos = end
            return out
        parts = [self._buf[self._pos :]]
        need = n - len(parts[0])
        self._pos = len(self._buf)
        while need > 0 and self._next_block():
            part = self._buf[:need]
            parts.append(part)
            self._pos = len(part)
            need -= len(part)
        return b"".join(parts)

    def tell(self):
        return self._offset + self._pos

    def close(self):
        self.f.close()


def _open_stream(file_path):
    with open(file_path, "rb") as f:
        magic_number = f.read(4)

    if magic_number[:2] == b"\x1f\x8b":
        logging.debug("RDS is compressed")
        return ByteStream(gzip.open(file_path, "rb"))
    if magic_number == QS2_MAGIC:
        logging.debug("qs2 format detected")
        return Qs2Reader(file_path)
    if magic_number == QDATA_MAGIC:
        raise NotImplementedError(
            "qdata format (qs2::qd_save) is not supported, use qs2::qs_save"
        )
    if magic_number == QS_LEGACY_MAGIC:
        raise NotImplementedError(
            "legacy qs format (qs::qsave) is not supported, use qs2::qs_save"
        )
    return ByteStream(open(file_path, "rb"))
