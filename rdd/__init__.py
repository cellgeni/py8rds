import gzip
import struct
import logging

__version__ = "0.0.1"

logging.basicConfig(level="DEBUG", format="[%(asctime)s][%(levelname)s] %(message)s")
logging.getLogger().setLevel('DEBUG')

# REFERENCES:
# https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats
# https://github.com/LTLA/rds2cpp/blob/master/include/rds2cpp/parse_object.hpp#L30
# https://blog.djnavarro.net/posts/2021-11-15_serialisation-with-rds/
# https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/main/serialize.c#L1855

SEXPTYPEs = {
   #-1: {"SEXPTYPE": "MINUS_ONE", "description": "NULL"}, # no idea what it is
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
    238: {"SEXPTYPE": "ALTREP", "description": "ALTREP_SXP / custom ALTREP"}, 
    242: {"SEXPTYPE": "EMPTYENV", "description": "NULL"},
    249: {"SEXPTYPE": "NAMESPACESXP", "description": "NULL"},
    251: {"SEXPTYPE": "MISSINGARG", "description": "NULL"},
    253: {"SEXPTYPE": "GLOBALENV", "description": "NULL"},
    254: {"SEXPTYPE": "NILSXP", "description": "NULL"},
    255: {"SEXPTYPE": "REF", "description": "reference to previous object"}, 
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
        
    def __str__(self):
      return self.object.toString()
    
    def __repr__(self):
      return self.object.toString()

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
              # not sure I understand where this None come from
              if k == None:
                result += indent+sep+sep + 'None'
              elif not isinstance(k,list):
                result += indent+sep+sep + str(k)
              elif isinstance(k[1],Robj):
                result += k[1].toString(indent=indent+sep+sep,withAttr=withAttr,maxItems=maxItems,header=str(k[0])+":",sep=sep) + ",\n"
              else:
                result += indent+sep+sep + str(k[1])
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
              result += k[1].toString(indent=indent+sep+sep,withAttr=withAttr,maxItems=maxItems,header=str(k[0])+":",sep=sep) + ","
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
    if any(len(b) == 0 for b in object_header):
        raise EOFError("Unexpected end of file while reading object header")

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
            "attributes": attributes_bit & (sexp_type < 255),  # sexp_type==255 is symbol, symbol has no attributes (I hope), but have seqpool index in preceeding bytes
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
        logging.debug(">> LISTSXP")
        # For pairlists, attributes are serialized BEFORE TAG/CAR/CDR.
        attrs = None
        if header["flags"]["attributes"]:
            logging.debug("---> start of LISTSXP attributes")
            attrs = parse_object(reader).value
            logging.debug("<--- end of LISTSXP attributes")
        robj.value = parse_LISTSXP(reader, header["flags"]["tag"])
        if attrs is not None:
            robj.attributes = attrs

    elif header["sexptype"] == "LANGSXP":
        logging.debug(">> LANGSXP")
        # Language objects have the same pairlist layout as LISTSXP.
        attrs = None
        if header["flags"]["attributes"]:
            logging.debug("---> start of LANGSXP attributes")
            attrs = parse_object(reader).value
            logging.debug("<--- end of LANGSXP attributes")
        robj.value = parse_LISTSXP(reader, header["flags"]["tag"])
        if attrs is not None:
            robj.attributes = attrs

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
        robj.value = parse_CLOSXP(reader,header["flags"]["attributes"])
    
    elif header["sexptype"] == "ENVSXP":
        logging.debug(f">> ENVSXP")
        robj.value = parse_ENVSXP(reader)
    
    elif header["sexptype"] == "BCODESXP":
        logging.debug(">> BCODESXP (bytecode)")
        robj.value = parse_BCODESXP(reader)

    elif header["sexptype"] in ["NILSXP",'MISSINGARG','GLOBALENV','EMPTYENV','MINUS_ONE']:
        robj.value = None
        
    elif header["sexptype"] == "NAMESPACESXP":
        logging.debug(">> NAMESPACESXP")
        robj.value = parse_NAMESPACESXP(reader)
    # REFSXP / pooled symbol reference
    elif header["sexptype"] == "REF":
        # In the original hack this was treated purely as a SEQPOOL
        # entry for symbols. That works for simple RDS files but fails
        # once REF is also used for environments, functions, etc.
        #
        # The reference index is encoded in the same 4 bytes as the flag,
        # and the original code derived a 'pool index' like this:
        #     inx = (size + 1) // 256 - 2
        size = struct.unpack(">I", b"".join(header["bytes"]))[0]
        inx = (size + 1) // 256 - 2
        logging.debug(f"REF header='{size}'; derived index='{inx}'")

        if 0 <= inx < len(SEQPOOL):
            # Looks like this REF is actually pointing into the symbol pool.
            robj.value = SEQPOOL[inx]
            logging.debug(f"REF resolved from SEQPOOL[{inx}] = {robj.value!r}")
        else:
            # This REF is not a pooled symbol (e.g. environment, other object),
            # and we don’t currently maintain a full reference table. To keep
            # parsing robust and avoid crashes, we return a placeholder.
            logging.warning(
                f"REF index {inx} out of SEQPOOL range (len={len(SEQPOOL)}); "
                "treating as opaque reference placeholder"
            )
            robj.value = f"<REF_{inx}>"
    elif header["sexptype"] == "ALTREP":
        logging.debug(">> ALTREP")
        robj.value = parse_ALTREP(reader)

    else:
      raise NotImplementedError(f"unimplemented SEXPTYPE '{header['sexptype']}'")

    logging.debug(header)
    # S4 fields are stored as attribures (pairlist), so attr flag is true, but we will store them as data.
    # so at this point for S4 we've already read attributes
    if header["flags"]["attributes"] and header["sexptype"] not in ('S4SXP','ENVSXP','CLOSXP', 'LISTSXP', 'LANGSXP','ALTREP','NAMESPACESXP'):
        logging.debug("---> start of object attributes")
        attributes = parse_object(reader).value
        if attributes is not None:
          robj.attributes = attributes
        logging.debug("<--- end of object attributes")
    
    robj.attributes.append(["sexptype", header["sexptype"]])
    return robj

def parse_S4SXP(reader):
    return parse_object(reader) # S4 is stored as LISTSXP, so just continue

def parse_LISTSXP(reader, has_tag):
    """
    Parse a dotted pair / pairlist object (LISTSXP, also used for LANGSXP).

    Serialization layout for each cons cell (serialize.c):

        [ATTRIB] (handled outside this function)
        [TAG]  (only if has_tag == True)
        [CAR]
        [CDR]  (another LISTSXP cell, or NILSXP terminator)

    'has_tag' is constant for the entire list – either every cell has a
    tag or none of them do.
    """
    value = []

    # TAG or None
    if has_tag:
        pair_name = parse_object(reader).value
    else:
        pair_name = None

    # CAR (always present)
    pair_value = parse_object(reader)
    value.append([pair_name, pair_value])
    logging.debug(f"pair      '{pair_name}': '{pair_value}'")

    # CDR: either another pairlist cell or NULL terminator
    rest = parse_object(reader)

    if rest.value is not None:
        # A proper pairlist from our parser should be represented as a list
        # of [name, Robj] pairs (and not as a list of Robj instances).
        if isinstance(rest.value, list) and (
            len(rest.value) == 0 or not isinstance(rest.value[0], Robj)
        ):
            value.extend(rest.value)
        else:
            # Unexpected structure in CDR; preserve it as a single item so
            # we don't lose information.
            value.append([None, rest])

    return value

def parse_NAMESPACESXP(reader):
    """
    Parse a NAMESPACESXP (code 249).

    R's ReadItem does:

        SEXP s = InStringVec(stream, ref_table);
        s = R_FindNamespace1(s);

    InStringVec reads:
        - placeholder int (must be 0),
        - length (int),
        - then 'length' strings via ReadItem (usually CHARSXP).

    Here we:
      * Read the placeholder and length.
      * Read 'length' string objects with parse_object.
      * Return a Python list of strings (namespace spec), similar to STRSXP.

    The enclosing Robj will still get an attribute ['sexptype', 'NAMESPACESXP'],
    so you can distinguish it from a normal STRSXP if needed.
    """

    # Placeholder (should be 0).
    placeholder = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"NAMESPACESXP placeholder: {placeholder}")
    if placeholder != 0:
        logging.warning(
            f"NAMESPACESXP: expected placeholder 0 but found {placeholder}"
        )

    # Length of the string vector.
    length = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"NAMESPACESXP string vector length: {length}")

    strings = []
    for i in range(length):
        elt = parse_object(reader)   # typically CHARSXP Robj
        strings.append(elt.value)
        logging.debug(f"NAMESPACESXP element {i}: {elt.value!r}")

    # Just return the list of strings, STRSXP-style.
    return strings


def parse_SYMSXP(reader):
    # should be text representation of literals
    value = parse_object(reader).value
    SEQPOOL.append(value)
    return value

def parse_ALTREP(reader):
    """
    Parse an ALTREP object (ALTREP_SXP, code 238).

    According to R's serialize.c (ReadItem):

        case ALTREP_SXP:
        {
            R_ReadItemDepth++;
            SEXP info  = ReadItem(ref_table, stream);
            SEXP state = ReadItem(ref_table, stream);
            SEXP attr  = ReadItem(ref_table, stream);
            s = ALTREP_UNSERIALIZE_EX(info, state, attr, objf, levs);
            ...
        }

    We cannot call ALTREP_UNSERIALIZE_EX from Python, so we expose the
    three components directly:

        value = [
            ["info",  <Robj>],
            ["state", <Robj>],
            ["attr",  <Robj>],
        ]

    This keeps the stream aligned and lets Robj.toString() print something
    sensible for ALTREP objects.
    """

    logging.debug("ALTREP: parsing info")
    info = parse_object(reader)   # usually something indicating ALTREP class

    logging.debug("ALTREP: parsing state")
    state = parse_object(reader)  # ALTREP-specific state (e.g., length/start/step)

    logging.debug("ALTREP: parsing attr")
    attr = parse_object(reader)   # attributes that R would attach to the final object

    value = [
        ["info", info],
        ["state", state],
        ["attr", attr],
    ]

    logging.debug("ALTREP parsed into info/state/attr triple")
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
        logging.debug(f"VECSXP: {_}/{size}")
        value.append(parse_object(reader))
    logging.debug(f"value      '{value[:min(10,len(value))]}'")
    return value


def parse_ENVSXP(reader):
    """
    Parse an R environment (ENVSXP) according to R's serialize.c ReadItem:

        int locked = InInteger(stream);
        ENCLOS  = ReadItem(...)
        FRAME   = ReadItem(...)
        HASHTAB = ReadItem(...)
        ATTRIB  = ReadItem(...)

    We expose a simple structure as a pairlist-like list of [name, Robj]
    pairs, and include the locked flag and environment-specific attributes
    inside the value.
    """

    # locked flag
    locked = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"ENVSXP locked flag: {locked}")

    # ENCLOS: enclosing environment
    logging.debug("ENVSXP: parsing enclosure (ENCLOS)")
    enclos = parse_object(reader)

    # FRAME: symbol/value bindings pairlist
    logging.debug("ENVSXP: parsing frame (FRAME)")
    frame = parse_object(reader)

    # HASHTAB: hash table, often NULL or VECSXP
    logging.debug("ENVSXP: parsing hash table (HASHTAB)")
    hashtab = parse_object(reader)

    # ATTRIB: environment-specific attributes
    logging.debug("ENVSXP: parsing attributes (ATTRIB)")
    env_attrs = parse_object(reader)

    env_value = [
        ["locked", locked],
        ["enclos", enclos],
        ["frame", frame],
        ["hashtab", hashtab],
        ["env_attrs", env_attrs],
    ]

    return env_value

def parse_CLOSXP(reader, has_attributes=False):
    """
    Parse an R closure (CLOSXP).

    Layout:

        [CLOSXP header]
        [ATTRIB] ?  if has_attributes
        [CLOENV]
        [FORMALS]
        [BODY]
    """

    if has_attributes:
        logging.debug("CLOSXP has attributes -> consuming attribute object")
        _attr_robj = parse_object(reader)  # ignored, just advance stream
        logging.debug("CLOSXP attributes consumed (ignored)")

    logging.debug("CLOSXP: parsing environment (CLOENV)")
    env_robj = parse_object(reader)

    logging.debug("CLOSXP: parsing formals (FORMALS)")
    formals_robj = parse_object(reader)

    logging.debug("CLOSXP: parsing body (BODY)")
    body_robj = parse_object(reader)

    value = [
        ["env", env_robj],
        ["formals", formals_robj],
        ["body", body_robj],
    ]

    logging.debug("CLOSXP parsed into env/formals/body triple")
    return value


# start of BCODESXP ####
# Bytecode helper constants from serialize.c
BCODESXP_CODE    = 21   # standard R SEXP type
LANGSXP_CODE     = 6
LISTSXP_CODE     = 2

BCREPDEF_CODE    = 244
BCREPREF_CODE    = 243
ATTRLANGSXP_CODE = 240
ATTRLISTSXP_CODE = 239


def parse_BCODESXP(reader):
    """
    Parse a BCODESXP (bytecode object).

    R's ReadItem for BCODESXP calls ReadBC, which does:

        reps_len = InInteger(stream);
        reps     = allocVector(VECSXP, reps_len);
        s        = ReadBC1(ref_table, reps, stream);

    ReadBC1 then reads:
        code   = ReadItem(...)
        consts = ReadBCConsts(...)

    We mirror that structure but return a pairlist-style value:

        [
          ["code",   <Robj>],
          ["consts", <Robj with value = list of const Robj>]
        ]
    """

    # ---- reps length (for shared subtrees in bytecode language) ----
    reps_len = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"BCODESXP: reps_len = {reps_len}")
    reps = [None] * max(reps_len, 0)

    value = _read_bc1(reader, reps)
    return value


def _read_bc1(reader, reps):
    """
    Python version of ReadBC1(ref_table, reps, stream).

    Layout:

        code   = ReadItem(...)
        consts = ReadBCConsts(...)
    """

    logging.debug("BCODESXP: reading bytecode 'code' object")
    code_robj = parse_object(reader)

    logging.debug("BCODESXP: reading bytecode const pool")
    consts_list = _read_bc_consts(reader, reps)  # list of Robj / nested BC

    # Wrap consts list into an Robj so that it prints as a list.
    consts_robj = Robj()
    consts_robj.value = consts_list
    consts_robj.attributes = []

    # BCODESXP value is represented as a pairlist-style list of [name, Robj]
    # so that Robj.toString() prints it nicely.
    value = [
        ["code", code_robj],
        ["consts", consts_robj],
    ]

    logging.debug("BCODESXP: finished ReadBC1")
    return value


def _read_bc_consts(reader, reps):
    """
    Python version of ReadBCConsts(ref_table, reps, stream).

    Layout:

        n = InInteger(stream)
        for i in 0..n-1:
            type_code = InInteger(stream)
            ... then type-specific payload ...
    """
    n = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"BCODESXP consts: n = {n}")

    consts = []
    for i in range(n):
        type_code = struct.unpack(">i", reader.read(4))[0]
        logging.debug(f"BCODESXP const[{i}] type_code = {type_code}")

        if type_code == BCODESXP_CODE:
            # Nested bytecode object; use BC1 (no new reps_len).
            logging.debug("BCODESXP const: nested BCODESXP")
            const_val = _read_bc1(reader, reps)

            # Wrap nested BC in an Robj for consistency.
            nested = Robj()
            nested.value = const_val
            nested.attributes = [["sexptype", "BCODESXP"]]
            consts.append(nested)

        elif type_code in (
            LANGSXP_CODE,
            LISTSXP_CODE,
            BCREPDEF_CODE,
            BCREPREF_CODE,
            ATTRLANGSXP_CODE,
            ATTRLISTSXP_CODE,
        ):
            logging.debug("BCODESXP const: using ReadBCLang-style parser")
            const_val = _read_bc_lang(reader, type_code, reps)
            consts.append(const_val)

        else:
            # Default case in ReadBCConsts: ReadItem(...)
            logging.debug("BCODESXP const: delegating to parse_object()")
            const_val = parse_object(reader)
            consts.append(const_val)

    return consts


def _read_bc_lang(reader, type_code, reps):
    """
    Python version of ReadBCLang(int type, ...).

    This handles the compressed representation of language / pairlist
    objects inside the bytecode constant pool.

    It returns an Robj whose .value is a pairlist-like list of
    [name, Robj/primitive] pairs:

        [
          ["bctype", "LANGSXP" | "LISTSXP" | ...],
          ["tag",    <Robj>],
          ["car",    <Robj or bc-lang>],
          ["cdr",    <Robj or bc-lang>],
        ]
    """

    # BCREPREF: reference to a previously defined bc-lang node.
    if type_code == BCREPREF_CODE:
        idx = struct.unpack(">i", reader.read(4))[0]
        logging.debug(f"ReadBCLang: BCREPREF idx={idx}")
        if 0 <= idx < len(reps):
            return reps[idx]
        logging.warning(f"ReadBCLang: BCREPREF idx={idx} out of range")
        return None

    # BCREPDEF / LANGSXP / LISTSXP / ATTRLANGSXP / ATTRLISTSXP
    pos = -1
    hasattr = False
    t = type_code

    if t == BCREPDEF_CODE:
        # BCREPDEF: position + actual type
        pos = struct.unpack(">i", reader.read(4))[0]
        t = struct.unpack(">i", reader.read(4))[0]
        logging.debug(f"ReadBCLang: BCREPDEF pos={pos}, underlying type={t}")

    # Handle attribute-wrapped types.
    if t == ATTRLANGSXP_CODE:
        t = LANGSXP_CODE
        hasattr = True
    elif t == ATTRLISTSXP_CODE:
        t = LISTSXP_CODE
        hasattr = True

    # Determine a human-readable type name.
    type_name = SEXPTYPEs.get(t, {}).get("SEXPTYPE", f"TYPE_{t}")

    node = Robj()
    node.attributes = []

    # Optional attributes on the bc-lang node.
    attr_robj = None
    if hasattr:
        logging.debug("ReadBCLang: reading attributes")
        attr_robj = parse_object(reader)

    logging.debug("ReadBCLang: reading TAG")
    tag_robj = parse_object(reader)

    logging.debug("ReadBCLang: reading CAR type")
    car_type = struct.unpack(">i", reader.read(4))[0]
    car_val = _read_bc_lang(reader, car_type, reps)

    logging.debug("ReadBCLang: reading CDR type")
    cdr_type = struct.unpack(">i", reader.read(4))[0]
    cdr_val = _read_bc_lang(reader, cdr_type, reps)

    # Store in reps[] if this was a BCREPDEF node.
    if 0 <= pos < len(reps):
        reps[pos] = node

    # Represent the bc-lang node as a pairlist-style value.
    value = [
        ["bctype", type_name],
        ["tag", tag_robj],
        ["car", car_val],
        ["cdr", cdr_val],
    ]
    node.value = value

    # If there were attributes, attach them as normal Robj attributes.
    if attr_robj is not None:
        # If the attributes object is itself a pairlist, we can just copy its
        # .value; otherwise wrap it in a simple pair.
        if isinstance(attr_robj.value, list):
            node.attributes = attr_robj.value
        else:
            node.attributes = [["attr", attr_robj]]

    return node

# end of BCODESXP ####

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

    
