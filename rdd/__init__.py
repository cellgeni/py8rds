import gzip
import struct
import logging

__version__ = "0.0.1"

logging.basicConfig(level="DEBUG", format="[%(asctime)s][%(levelname)s] %(message)s")

# REFERENCES:
# https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats
# https://github.com/LTLA/rds2cpp/blob/master/include/rds2cpp/parse_object.hpp#L30
# https://blog.djnavarro.net/posts/2021-11-15_serialisation-with-rds/


SEXPTYPEs = {
    0: {"SEXPTYPE": "NILSXP", "description": "NULL"},
    1: {"SEXPTYPE": "SYMSXP", "description": "symbols"},
    2: {"SEXPTYPE": "LISTSXP", "description": "pairlists"},
    3: {"SEXPTYPE": "CLOSXP", "description": "closures"},
    4: {"SEXPTYPE": "ENVSXP", "description": "environments"},
    5: {"SEXPTYPE": "PROMSXP", "description": "promises"},
    6: {"SEXPTYPE": "LANGSXP", "description": "languageobjects"},
    7: {"SEXPTYPE": "SPECIALSXP", "description": "special functions"},
    8: {"SEXPTYPE": "BUILTINSXP", "description": "builtinfunctions"},
    9: {"SEXPTYPE": "CHARSXP", "description": "internalcharacterstrings"},
    10: {"SEXPTYPE": "LGLSXP", "description": "logicalvectors"},
    13: {"SEXPTYPE": "INTSXP", "description": "integervectors"},
    14: {"SEXPTYPE": "REALSXP", "description": "numericvectors"},
    15: {"SEXPTYPE": "CPLXSXP", "description": "complexvectors"},
    16: {"SEXPTYPE": "STRSXP", "description": "charactervectors"},
    17: {"SEXPTYPE": "DOTSXP", "description": "dot-dot-dotobject"},
    18: {"SEXPTYPE": "ANYSXP", "description": "make'any'argswork"},
    19: {"SEXPTYPE": "VECSXP", "description": "list(genericvector)"},
    20: {"SEXPTYPE": "EXPRSXP", "description": "expressionvector"},
    21: {"SEXPTYPE": "BCODESXP", "description": "bytecode"},
    22: {"SEXPTYPE": "EXTPTRSXP", "description": "externalpointer"},
    23: {"SEXPTYPE": "WEAKREFSXP", "description": "weak"},
    24: {"SEXPTYPE": "RAWSXP", "description": "reference raw vector"},
    25: {"SEXPTYPE": "S4SXP", "description": "S4 classes"},
    254: {"SEXPTYPE": "NILSXP", "description": "NULL"},
}


class RdsFile:
    def __init__(self):
        self.format_version = 0
        self.writer_version = [0, 0, 0]
        self.reader_version = [0, 0, 0]
        self.encoding = ""
        self.object = None
        self.environments = []
        self.symbols = []
        self.external_pointers = []


def parse_object_header(reader):
    logging.debug(f"---> start of object header (0x{reader.tell():02X})")

    object_header = [reader.read(1) for _ in range(4)]
    header_bin = [f"{struct.unpack('>B',b)[0]:08b}" for b in object_header]
    logging.debug(f"header {header_bin}")

    sexp_type = struct.unpack("<b", object_header[3])[0]  # bit 0-7
    if sexp_type in SEXPTYPEs:
        logging.debug(f"SEXPTYPE   {SEXPTYPEs[sexp_type]['SEXPTYPE']} ({sexp_type})")
    else:
        logging.warning(f"SEXPTYPE   UNKNOWN! ({sexp_type})")

    other_bits = struct.unpack(">B", object_header[2])[0]
    object_bit = bool(other_bits & 0b00000001)  # bit 9
    attributes_bit = bool((other_bits & 0b00000010) >> 1)  # bit 10
    tag_bit = bool((other_bits & 0b00000100) >> 2)  # bit 11 (pairlist types)
    gp_bit = bool(struct.unpack(">B", object_header[1])[0] & 0b00000001)
    logging.debug(f"is object  {object_bit}")
    logging.debug(f"has attr   {attributes_bit}")
    logging.debug(f"has tag    {tag_bit}")
    logging.debug(f"gp_bit     {gp_bit}")

    logging.debug("<--- end of object header")

    return {
        "bytes": object_header,
        "sexptype": SEXPTYPEs[sexp_type]["SEXPTYPE"] if sexp_type in SEXPTYPEs else -1,
        "flags": {
            "object": object_bit,
            "attributes": attributes_bit,
            "tag": tag_bit,
            "gp": gp_bit,
        },
    }


def parse_rds(file_path: str) -> RdsFile:
    logging.info(f"Reading {file_path}")
    rds = RdsFile()
    with open(file_path, "rb") as f:
        magic_number = f.read(2)

    if magic_number == b"\x1f\x8b":
        logging.debug("RDS is compressed")
        reader = gzip.open(file_path, "rb")
    else:
        reader = open(file_path, "rb")

    logging.debug("---> start of file header")
    format = reader.read(2)
    if format == b"X\n":
        logging.debug("XDR binary format detected")
    elif format == b"A\n":
        raise NotImplementedError("ASCII format detected")
    elif format == b"B\n":
        raise NotImplementedError("Native word-order binary format detected")
    else:
        raise NotImplementedError(f"Unknown Format '{format}'")

    rds.format_version = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"version: {rds.format_version}")

    rds.writer_version = struct.unpack("BBB", reader.read(4)[1:])
    logging.debug(
        f"version of R which wrote the file: {'.'.join(map(str,rds.writer_version))}"
    )

    rds.reader_version = struct.unpack("BBB", reader.read(4)[1:])
    logging.debug(
        f"minimal version of R needed to read the format: {'.'.join(map(str,rds.reader_version))}"
    )

    if rds.format_version == 3:
        encoding_length = struct.unpack(">I", reader.read(4))[0]
        rds.encoding = reader.read(encoding_length).decode()
        logging.debug(f"encoding: {rds.encoding}")

    logging.debug("<--- end of file header")

    x = parse_object(reader)


def parse_object(reader):
    header = parse_object_header(reader)

    if header["flags"]["attributes"]:
        attributes = []
        while True:
            attr = parse_object(reader)
            if attr == None:
                # NILSXP? re-read the hader befor breaking
                # so i can be processed and we're not stuck in a loop
                header = parse_object_header(reader)
                break
            attributes.append(attr)

    if header["sexptype"] == "S4SXP":
        logging.debug(f">> S4SXP")
        return parse_S4SXP(reader)

    if header["sexptype"] == "LISTSXP":
        logging.debug(f">> LISTSXP")
        return parse_LISTSXP(reader)

    if header["sexptype"] == "LANGSXP":
        logging.debug(f">> LANGSXP")
        return parse_LANGSXP(reader)

    if header["sexptype"] == "SYMSXP":
        logging.debug(f">> SYMSXP expecting CHARSXP after")
        return parse_SYMSXP(reader)

    if header["sexptype"] == "CHARSXP":
        logging.debug(f">> CHARSXP")
        return parse_CHARSXP(reader)

    if header["sexptype"] == "LGLSXP":
        logging.debug(f">> LGLSXP")
        return parse_LGLSXP(reader)

    if header["sexptype"] == "STRSXP":
        logging.debug(f">> STRSXP")
        return parse_STRSXP(reader)

    if header["sexptype"] == "REALSXP":
        logging.debug(f">> REALSXP")
        return parse_REALSXP(reader)

    if header["sexptype"] == "INTSXP":
        logging.debug(f">> INTSXP")
        return parse_INTSXP(reader)

    if header["sexptype"] == "BUILTINSXP":
        logging.debug(f">> BUILTINSXP")
        return parse_BUILTINSXP(reader)

    if header["sexptype"] == "VECSXP":
        logging.debug(f">> VECSXP")
        return parse_VECSXP(reader)

    if header["sexptype"] == "CLOSXP":
        logging.debug(f">> CLOSXP")
        return parse_CLOSXP(reader)

    if header["sexptype"] == "NILSXP":
        return None

    logging.debug(header)


def parse_S4SXP(reader):
    return parse_object(reader)


def parse_LISTSXP(reader):
    return parse_object(reader)


def parse_SYMSXP(reader):
    return parse_object(reader)


def parse_CHARSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = reader.read(size).decode()
    logging.debug(f"value      '{value}'")
    return value


def parse_STRSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        value.append(parse_object(reader))
    logging.debug(f"STRSXP value      {value}")
    # value = reader.read(size)
    # try:
    #     value = value.decode()
    # except:
    #     pass
    # logging.debug(f"value      '{value}'")
    return value


def parse_REALSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        value.append(struct.unpack(">d", reader.read(8))[0])
    if len(value) > 10:
        logging.debug(f"value      '{value[:10]}...[more]'")
    else:
        logging.debug(f"value      '{value}'")
    return value


def parse_INTSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        value.append(struct.unpack(">i", reader.read(4))[0])
    if len(value) > 10:
        logging.debug(f"value      '{value[:10]}...[more]'")
    else:
        logging.debug(f"value      '{value}'")
    return value


def parse_BUILTINSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = reader.read(size).decode()
    logging.debug(f"value      '{value}'")
    return value


def parse_LGLSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = "fake"
    reader.seek(size * 4, 1)  # seek fakes reading into memory
    # for _ in range(size):
    #    value.append(struct.unpack('<I', reader.read(4))[0])
    # logging.debug(f"value      '{value[:10]}'")
    return value


def parse_VECSXP(reader):
    size = struct.unpack(">Q", reader.read(8))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        value.append(parse_object(reader))
    logging.debug(f"value      '{value}'")
    return value


def parse_LANGSXP(reader):
    return parse_object(reader)


def parse_CLOSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = "fake"
    reader.seek(size, 1)  # seek fakes reading into memory
    # value = reader.read(size)
    # logging.debug(f"value      '{value}'")
    return value


if __name__ == "__main__":
    ## TESTS
    rds_file = parse_rds("test_rds/pbmc_tutorial.rds")
    # print("============")
    # rds_file = parse_rds("test_rds/211014KTOrganoid.rds")
    # print("============")
    # rds_file = parse_rds("test_rds/pbmc3k_final.rds")
    # print("============")
    # rds_file = parse_rds("test_rds/data.frame.rds")
    # print("============")
    # rds_file = parse_rds("test_rds/a.rds")
    # print("============")
    # rds_file = parse_rds("test_rds/list.rds")
