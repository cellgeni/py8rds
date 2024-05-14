import gzip
import struct
import logging

__version__ = "0.0.1"

logging.basicConfig(level="DEBUG", format="[%(asctime)s][%(levelname)s] %(message)s")
logging.getLogger().setLevel('INFO')

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
    255: {"SEXPTYPE": "SEQPOOL", "description": "NULL"}, # it is not official sextype, but in it will make parsing seqpooled symbols smoother
}

# pairlist names seems to be stored into pool and be referenced on next occurences
# it should be probably placed into RdsFile or Robj, but they are not passed to parse functions, so I'll keep it on the top level
SEQPOOL = []


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

class Robj:
    def __init__(self):
        self.value = None
        # attributes has unique names and then can be stored in dict, but they are serialized as pairlist which names are not necessary unique. 
        # So we will store attributes as list of [string,Robj] pairs
        self.attributes = []
    def __str__(self):
      return self.toString()
    
    def __repr__(self):
      return self.toString()
    
    # it is to print Robj in human readable manner
    def isSimple(self,o):
      if isinstance(o,list):
        types = list({type(x) for x in o if x != None})
      else:
        return False
      return (len(types) == 1) and any([types[0] == type(t) for t in ['a',1,1.1,True]]) # ugly but working
    
    def toString(self,indent='',withAttr=True,maxItems=5,header='',printBrackets = True,sep='.'):
      result = ''
      if printBrackets:
        result += indent + header +"(\n"
      result += indent + sep+"value: "
      # atomic lists, should be always stored as list of same type
      if self.isSimple(self.value):
        v = self.value[:min(maxItems,len(self.value))]
        v = '[' + ','.join([str(x) for x in v]) + (',...]' if len(self.value) > maxItems else ']')
        result += v
      else:
        # list, stored as list of Robj or pairlist stored as list of [name,value]
        if isinstance(self.value,list):
          if len(self.value)==0:
            result += '[]'
          # list
          elif isinstance(self.value[0],Robj):
            result += "[\n"
            for v in self.value:
              if v != None:
                result += v.toString(indent=indent+sep+sep,withAttr=withAttr,maxItems=maxItems,sep=sep) + ",\n"
              else:
                result += indent+sep+sep +"(None),\n"
            result = result[:(len(result)-2)]
            result += "\n" + indent + sep+"]"
          # pairlist (including S4)
          else:
            result += "{\n"
            for k in self.value:
              result += k[1].toString(indent=indent+sep+sep,withAttr=withAttr,maxItems=maxItems,header=k[0]+":",sep=sep) + ",\n"
            result = result[:(len(result)-2)]
            result += "\n" + indent + sep+"}"
        # nested Robj
        if isinstance(self.value,Robj):
          result += "(\n"
          result += self.value.toString(indent=indent+sep,withAttr=withAttr,maxItems=maxItems,printBrackets=False,sep=sep)
          result += "\n" + indent + sep+")"
        
      if withAttr:
        result += "\n" + indent + sep+'attrs: {'
        for k in self.attributes:
            result += "\n" 
            if isinstance(k[1],Robj):
              result += k[1].toString(indent=indent+sep+sep,withAttr=withAttr,maxItems=maxItems,header=k[0]+":",sep=sep) + ","
            else:
              result += indent + sep +sep + k[0] +": " + str(k[1]) + ","
        if len(self.attributes)>0:
            result = result[:(len(result)-1)]
        result += "}"
      if printBrackets:
        result += "\n" + indent + ")"
      return result
      

def parse_object_header(reader):
    logging.debug(f"---> start of object header (0x{reader.tell():02X})")

    object_header = [reader.read(1) for _ in range(4)]
    header_bin = [f"{struct.unpack('>B',b)[0]:08b}" for b in object_header]
    logging.debug(f"header {header_bin}")

    sexp_type = struct.unpack("<B", object_header[3])[0]  # bit 0-7
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
        "sexp_type_code": sexp_type,
        "sexptype": SEXPTYPEs[sexp_type]["SEXPTYPE"] if sexp_type in SEXPTYPEs else -1,
        "flags": {
            "object": object_bit,
            "attributes": attributes_bit & (sexp_type < 255),  # sexp_type==255 is symbol, symbol has no attributes (I hope), but have sepool index in preceeding bytes
            "tag": tag_bit,
            "gp": gp_bit,
        },
    }


def parse_rds(file_path: str) -> RdsFile:
    logging.info(f"Reading {file_path}")
    # clean seqpool before parsing new file
    global SEQPOOL
    SEQPOOL = [] 
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

    rds.object = parse_object(reader)
    logging.info(f"Done reading {file_path}")
    return rds

def parse_object(reader):
    header = parse_object_header(reader)
    robj = Robj()
    if header["sexptype"] == "S4SXP":
        logging.debug(f">> S4SXP")
        robj.value = parse_S4SXP(reader)

    elif header["sexptype"] == "LISTSXP":
        logging.debug(f">> LISTSXP")
        robj.value = parse_LISTSXP(reader)

    elif header["sexptype"] == "LANGSXP":
        logging.debug(f">> LANGSXP")
        robj.value = parse_LANGSXP(reader)

    elif header["sexptype"] == "SYMSXP":
        logging.debug(f">> SYMSXP expecting CHARSXP after")
        robj.value = parse_SYMSXP(reader)

    elif header["sexptype"] == "CHARSXP":
        logging.debug(f">> CHARSXP")
        robj.value = parse_CHARSXP(reader)

    elif header["sexptype"] == "LGLSXP":
        logging.debug(f">> LGLSXP")
        robj.value = parse_LGLSXP(reader)

    elif header["sexptype"] == "STRSXP":
        logging.debug(f">> STRSXP")
        robj.value = parse_STRSXP(reader)

    elif header["sexptype"] == "REALSXP":
        logging.debug(f">> REALSXP")
        robj.value = parse_REALSXP(reader)

    elif header["sexptype"] == "INTSXP":
        logging.debug(f">> INTSXP")
        robj.value = parse_INTSXP(reader)

    elif header["sexptype"] == "BUILTINSXP":
        logging.debug(f">> BUILTINSXP")
        robj.value = parse_BUILTINSXP(reader)

    elif header["sexptype"] == "VECSXP":
        logging.debug(f">> VECSXP")
        robj.value = parse_VECSXP(reader)

    elif header["sexptype"] == "CLOSXP":
        logging.debug(f">> CLOSXP")
        robj.value = parse_CLOSXP(reader)
    
    elif header["sexptype"] == "ENVSXP":
        logging.debug(f">> ENVSXP")
        robj.value = parse_ENVSXP(reader)

    elif header["sexptype"] == "NILSXP":
        robj.value = None
    # assume it is symbol from seqpool
    elif header['sexptype'] == 'SEQPOOL':
        size = struct.unpack(">I",b''.join(header['bytes']))[0]
        inx = (size+1)//(2**8)-2
        logging.debug(f"seqpool      header='{size}'; inx='{inx}'")
        robj.value = SEQPOOL[inx]
    else:
      raise NotImplementedError(f"unimplemented SEXPTYPE '{header['sexptype']}'")

    logging.debug(header)
    # S4 fields are stored as attribures (pairlist), so attr flag is true, but we will store them as data.
    # so at this point for S4 we've already read attributes
    if header["flags"]["attributes"] & (header["sexptype"] != 'S4SXP'):
        logging.debug("---> start of object attributes")
        attributes = parse_object(reader).value
        if attributes != None:
          robj.attributes = attributes
        logging.debug("<--- end of object attributes")
    
    robj.attributes.append(["sexptype", header["sexptype"]])
    return robj

def parse_S4SXP(reader):
    return parse_object(reader) # S4 is stored as LISTSXP, so just continue

def parse_LISTSXP(reader):
    # pairlist consist of sence of pairs [symbol,any R object], ended with NULL (254)
    value = []
    # assume that attribute name doesn't have own attributes
    pair_name = parse_object(reader).value
    pair_value = parse_object(reader)
    value.append([pair_name,pair_value])
    logging.debug(f"pair      '{pair_name}': '{pair_value}'")
    # parse rest of pairs, assume that individaul pair doesn't have its own attributes
    rest = parse_object(reader).value
    if rest != None:
        value.extend(rest)       # I assume attribute names are unique
    return value

def parse_SYMSXP(reader):
    # should be text representation of literals
    value = parse_object(reader).value
    SEQPOOL.append(value)
    return value


def parse_CHARSXP(reader):
    size = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"size       {size}")
    if size != -1:
      value = reader.read(size).decode()
    else:
      value = None
    logging.debug(f"value      '{value}'")
    return value


def parse_STRSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        value.append(parse_object(reader).value)
    logging.debug(f"STRSXP value      {value[:min(10,len(value))]}")
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
    value = [x if x != -2147483648 else None for x in value] # -2147483648 is NA
    return value


def parse_BUILTINSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = reader.read(size).decode()
    logging.debug(f"value      '{value}'")
    return value


def parse_LGLSXP(reader):
    # logical are stored as integers 
    value = parse_INTSXP(reader)
    return [v == 1 if v!= None else v for v in value]


def parse_VECSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        value.append(parse_object(reader))
    logging.debug(f"value      '{value[:min(10,len(value))]}'")
    return value


def parse_LANGSXP(reader):
    raise NotImplementedError("parse_LANGSXP is not implemented")

def parse_ENVSXP(reader):
    raise NotImplementedError("parse_ENVSXP is not implemented")


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
    rds_file = parse_rds("test_rds/seurat.rds") # should work on 10xPBMC
    print(rds_file.object) # human readable representation
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
    # print("============")
    # rds_file = parse_rds("test_rds/list.rds")
    # rds below are not working
    print('!!!!!!!!!!!!!!!!!!!!!!!!')
    rds_file = parse_rds("test_rds/env.rds") # environments are not implemented
    # these two are not working for yet unkown reasons
    rds_file = parse_rds('/nfs/users/nfs_p/pm19/nfs/projects/2303.bcc.skin/data.nfs/tic-2123/bcc.epi-3-CG.rds')
    rds_file = parse_rds('/nfs/users/nfs_p/pm19/nfs/projects/2303.bcc.skin/data.nfs/nmf.frq.r5.n100.rds')

    
