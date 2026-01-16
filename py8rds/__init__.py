import gzip
import struct
import logging
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import numpy as np
from numpy.dtypes import StringDType

__version__ = "0.0.1"

logging.basicConfig(level="DEBUG", format="[%(asctime)s][%(levelname)s] %(message)s")
logging.getLogger().setLevel('INFO')

# REFERENCES:
# https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats
# https://github.com/LTLA/rds2cpp/blob/master/include/rds2cpp/parse_object.hpp#L30
# https://blog.djnavarro.net/posts/2021-11-15_serialisation-with-rds/
# https://github.com/wch/r-source/blob/79298c499218846d14500255efd622b5021c10ec/src/main/serialize.c#L1855

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
    238: {"SEXPTYPE": "ALTREP", "description": "ALTREP_SXP / custom ALTREP"}, 
    242: {"SEXPTYPE": "EMPTYENV", "description": "NULL"},
    249: {"SEXPTYPE": "NAMESPACESXP", "description": "NULL"},
    251: {"SEXPTYPE": "MISSINGARG", "description": "NULL"},
    253: {"SEXPTYPE": "GLOBALENV", "description": "NULL"},
    254: {"SEXPTYPE": "NILSXP", "description": "NULL"},
    255: {"SEXPTYPE": "REFSXP", "description": "reference to previous object"}, 
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
        

class Robj:
    def __init__(self):
        # can be: Robj or list of scalars or Robjs or list of [scalar,Robj] pairs/Nones
        self.value = None
        # attributes has unique names and then can be stored in dict, but they are serialized as pairlist which names are not necessary unique. 
        # So we will store attributes as list of [scalar,Robj] pairs
        self.attributes = []
    
    
    def get(self,inxs):
      """
        Recursively navigate into an Robj using indices and/or names.
        Indeces are used to only search in value, names are used to search both in attributes (first) 
        and then in values if they are paired list (should be the case for S4 slots)

        Parameters
        ----------
        inxs : int | str | list[int | str]
            A single index/name or a list describing a path withint R obj

        Returns
        -------
        Any
            - Another `Robj`,
            - a primitive Python value (int, float, bool, str, etc.),
            - or `None` if any step along the path cannot be resolved.

        Examples
        --------
        Assuming an R list like `list(a = 1:3, b = 4:5)`:

        >>> robj.get("a")
        <Robj for 1:3>
        >>> robj.get(["a", 1])
        2
        """
      if not isinstance(inxs,list):
        inxs = [inxs]
      r = self._get(inxs[0])
      if (r is None) or (len(inxs) == 1):
        return r
      return r.get(inxs[1:])
    
    def _get(self,inx):
      # look by index
      if isinstance(inx,int):
        if isinstance(self.value,list):
          return self.value[inx]
        elif inx == 0:
          return self.value
        else:
          return None
      # look by name
      elif isinstance(inx,str):
        # in attributes
        for a in self.attributes:
          if inx == a[0]:
            return a[1]
        # in values
        if not isinstance(self.value,list):
          return None
        for v in self.value:
          if isinstance(v,list):
            if inx == v[0]:
              return v[1]
      return None
    
    def getClass(self):
      r = self._get('class')
      if r is None:
        r = self._get('sexptype')
      else:
        r = ','.join(r.value)
      return r
      
    def is_primitive(self,x):
      return isinstance(x, (int, float, bool, str, bytes, type(None)))
    
    def show(self,level=1):
      print(self.toString(level=level))
    
    def toString(self,name='Robj',indent='',maxItems=5,level=1e3):
      if level < 0:
        return ""
      level -= 1
      result = indent + name + '('+str(self.getClass())+"): "
      indent = indent.replace('+','|').replace('*','|').replace('&','|')
      # value
      values = self.value
      if (not isinstance(values,list)) and (not isinstance(values, np.ndarray)) :
        values = [values]
      if len(values) == 0:
        result += '[]\n'
      elif self.is_primitive(values[0]):
        result += '['
        result += ','.join([str(v) for v in values[:min(3,len(values))]])
        if len(values) > 3:
          result += ',...'
        result += ']\n'
      elif isinstance(values[0],Robj):
        result += "\n"
        for i in range(len(values)):
          result += values[i].toString(name=str(i),indent=indent+"+",maxItems=maxItems,level=level)
      else:
        result += "\n"
        for i in range(len(values)):
          if self.is_primitive(values[i]):
            result += indent+"&" + str(values[i]) +"\n"
          elif isinstance(values[i][1],Robj):
            result += values[i][1].toString(name=str(values[i][0]),indent=indent+"&",maxItems=maxItems,level=level)
          else:
            result += indent+"&"+str(values[i][0])+":"+str(values[i][1])+"\n"
      # attributes
      for a in self.attributes:
        if len(a) != 2:
          raise RuntimeError("Atributes are not in pair:"+str(a))
        if a[0] != 'sexptype':
          result += a[1].toString(name=str(a[0]),indent=indent+"*",maxItems=maxItems,level=level)
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


def parse_rds(file_path: str) -> Robj:
    logging.info(f"Reading {file_path}")
    # clean seqpool before parsing new file
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

    rds.object = parse_object(reader,rds)
    logging.info(f"Done reading {file_path}")
    return rds.object

def parse_object(reader,rds: RdsFile):
    header = parse_object_header(reader)
    robj = Robj()
    if header["sexptype"] == "S4SXP":
        logging.debug(f">> S4SXP")
        robj = parse_S4SXP(reader,rds)

    elif header["sexptype"] == "LISTSXP":
        logging.debug(">> LISTSXP")
        # For pairlists, attributes are serialized BEFORE TAG/CAR/CDR.
        attrs = None
        if header["flags"]["attributes"]:
            logging.debug("---> start of LISTSXP attributes")
            attrs = parse_object(reader).value
            logging.debug("<--- end of LISTSXP attributes")
        robj.value = parse_LISTSXP(reader,rds, header["flags"]["tag"])
        if attrs is not None:
            robj.attributes = attrs

    elif header["sexptype"] == "LANGSXP":
        logging.debug(">> LANGSXP")
        # Language objects have the same pairlist layout as LISTSXP.
        attrs = None
        if header["flags"]["attributes"]:
            logging.debug("---> start of LANGSXP attributes")
            attrs = parse_object(reader,rds).value
            logging.debug("<--- end of LANGSXP attributes")
        robj.value = parse_LISTSXP(reader, rds, header["flags"]["tag"])
        if attrs is not None:
            robj.attributes = attrs

    elif header["sexptype"] == "SYMSXP":
        logging.debug(f">> SYMSXP expecting CHARSXP after")
        robj.value = parse_SYMSXP(reader,rds)

    elif header["sexptype"] == "CHARSXP":
        logging.debug(f">> CHARSXP")
        robj.value = parse_CHARSXP(reader)

    elif header["sexptype"] == "LGLSXP":
        logging.debug(f">> LGLSXP")
        robj.value = parse_LGLSXP(reader)

    elif header["sexptype"] == "STRSXP":
        logging.debug(f">> STRSXP")
        robj.value = parse_STRSXP(reader,rds)

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
        robj.value = parse_VECSXP(reader,rds)

    elif header["sexptype"] == "CLOSXP":
        logging.debug(f">> CLOSXP")
        robj.value = parse_CLOSXP(reader,rds,header["flags"]["attributes"])
    
    elif header["sexptype"] == "ENVSXP":
        logging.debug(f">> ENVSXP")
        robj.value = parse_ENVSXP(reader,rds)
        rds.symbols.append(robj.value)
    
    elif header["sexptype"] == "BCODESXP":
        logging.debug(">> BCODESXP (bytecode)")
        robj.value = parse_BCODESXP(reader,rds)

    elif header["sexptype"] in ["NILSXP",'MISSINGARG','GLOBALENV','EMPTYENV','MINUS_ONE']:
        robj.value = None
        
    elif header["sexptype"] == "NAMESPACESXP":
        logging.debug(">> NAMESPACESXP")
        robj.value = parse_NAMESPACESXP(reader,rds)
        rds.symbols.append(robj.value)
    # REFSXP / pooled symbol reference
    elif header["sexptype"] == "REFSXP":
        size = struct.unpack(">I", b"".join(header["bytes"]))[0]
        inx = (size + 1) // 256 - 2
        logging.debug(f"REFSXP header='{size}'; derived index='{inx}'")

        if 0 <= inx < len(rds.symbols):
            # Looks like this REFSXP is actually pointing into the symbol pool.
            robj.value = rds.symbols[inx]
            logging.debug(f"REFSXP resolved from rds.symbols[{inx}] = {robj.value!r}")
        else:
            # This REFSXP is not a pooled symbol (e.g. environment, other object),
            # and we don’t currently maintain a full reference table. To keep
            # parsing robust and avoid crashes, we return a placeholder.
            logging.warning(
                f"REFSXP index {inx} out of rds.symbols range (len={len(rds.symbols)}); "
                "treating as opaque reference placeholder"
            )
            robj.value = f"<REFSXP_{inx}>"
    elif header["sexptype"] == "ALTREP":
        logging.debug(">> ALTREP")
        robj.value = parse_ALTREP(reader,rds)
    elif header["sexptype"] == "EXTPTRSXP":
      logging.debug(">> EXTPTRSXP")
      val, attrs = parse_EXTPTRSXP(reader, rds)
      robj.value = val
      if attrs is not None:
          robj.attributes = attrs
          
    else:
      raise NotImplementedError(f"unimplemented SEXPTYPE '{header['sexptype']}'")

    logging.debug(header)
    # S4 fields are stored as attribures (pairlist), so attr flag is true, but we will store them as data.
    # so at this point for S4 we've already read attributes
    if header["flags"]["attributes"] and header["sexptype"] not in ('S4SXP','ENVSXP','CLOSXP', 'LISTSXP', 'LANGSXP','ALTREP','NAMESPACESXP','EXTPTRSXP'):
        logging.debug("---> start of object attributes")
        attributes = parse_object(reader,rds).value
        if attributes is not None:
          robj.attributes = attributes
        logging.debug("<--- end of object attributes")
    
    # shit happens
    if not isinstance(robj.attributes, list):
      robj.attributes = [["__raw_attributes__",robj.attributes]]
            
    robj.attributes.append(["sexptype", header["sexptype"]])
    
    return robj

def parse_S4SXP(reader,rds):
    return parse_object(reader,rds) # S4 is stored as LISTSXP, so just continue

def parse_LISTSXP(reader, rds, has_tag):
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
        pair_name = parse_object(reader,rds).value
    else:
        pair_name = None

    # CAR (always present)
    pair_value = parse_object(reader,rds)
    value.append([pair_name, pair_value])

    # CDR: either another pairlist cell or NULL terminator
    rest = parse_object(reader,rds)

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

def parse_NAMESPACESXP(reader,rds):
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
        elt = parse_object(reader,rds)   # typically CHARSXP Robj
        strings.append(elt.value)
        logging.debug(f"NAMESPACESXP element {i}: {elt.value!r}")

    # Just return the list of strings, STRSXP-style.
    return strings


def parse_SYMSXP(reader,rds):
    # should be text representation of literals
    value = parse_object(reader,rds).value
    rds.symbols.append(value)
    return value

def parse_ALTREP(reader,rds):
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
    info = parse_object(reader,rds)   # usually something indicating ALTREP class

    logging.debug("ALTREP: parsing state")
    state = parse_object(reader,rds)  # ALTREP-specific state (e.g., length/start/step)

    logging.debug("ALTREP: parsing attr")
    attr = parse_object(reader,rds)   # attributes that R would attach to the final object

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

    if size == -1:
        value = None
    else:
        raw = reader.read(size)
        # Try UTF-8 first (what most modern R sessions use)
        try:
            value = raw.decode("utf-8")
        except UnicodeDecodeError:
            # Fall back to Latin-1; this matches how R can keep 8-bit data
            # around in non-UTF-8 encodings without blowing up.
            try:
                value = raw.decode("latin-1")
                logging.debug("CHARSXP: decoded using latin-1 fallback")
            except UnicodeDecodeError:
                # As an absolute fallback, replace invalid bytes so we never crash.
                value = raw.decode("utf-8", errors="replace")
                logging.warning(
                    "CHARSXP: could not decode as utf-8 or latin-1 cleanly; "
                    "using utf-8 with replacement characters"
                )
    return value



def parse_STRSXP(reader,rds):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = np.empty(size, dtype=StringDType())

    for i in range(size):
        value[i] = parse_object(reader,rds).value
    return value


def parse_REALSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = np.zeros(size, dtype=np.float64)
    for i in range(size):
        value[i] = struct.unpack(">d", reader.read(8))[0]
    return value


def parse_INTSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = np.zeros(size, dtype=np.int32)
    for i in range(size):
        # -2147483648 is NA but there is no nans in int32 so we just ignore it...
        value[i] = struct.unpack(">i", reader.read(4))[0]
    return value


def parse_BUILTINSXP(reader):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = reader.read(size).decode()
    return value


def parse_LGLSXP(reader):
    # logical are stored as integers 
    value = parse_INTSXP(reader)
    return [v == 1 if v!= None else v for v in value]


def parse_VECSXP(reader,rds):
    size = struct.unpack(">I", reader.read(4))[0]
    logging.debug(f"size       {size}")
    value = []
    for _ in range(size):
        logging.debug(f"VECSXP: {_}/{size}")
        value.append(parse_object(reader,rds))
    return value


def parse_ENVSXP(reader,rds):
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
    enclos = parse_object(reader,rds)

    # FRAME: symbol/value bindings pairlist
    logging.debug("ENVSXP: parsing frame (FRAME)")
    frame = parse_object(reader,rds)

    # HASHTAB: hash table, often NULL or VECSXP
    logging.debug("ENVSXP: parsing hash table (HASHTAB)")
    hashtab = parse_object(reader,rds)

    # ATTRIB: environment-specific attributes
    logging.debug("ENVSXP: parsing attributes (ATTRIB)")
    env_attrs = parse_object(reader,rds)

    env_value = [
        ["locked", locked],
        ["enclos", enclos],
        ["frame", frame],
        ["hashtab", hashtab],
        ["env_attrs", env_attrs],
    ]

    return env_value

def parse_CLOSXP(reader,rds, has_attributes=False):
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
        _attr_robj = parse_object(reader,rds)  # ignored, just advance stream
        logging.debug("CLOSXP attributes consumed (ignored)")

    logging.debug("CLOSXP: parsing environment (CLOENV)")
    env_robj = parse_object(reader,rds)

    logging.debug("CLOSXP: parsing formals (FORMALS)")
    formals_robj = parse_object(reader,rds)

    logging.debug("CLOSXP: parsing body (BODY)")
    body_robj = parse_object(reader,rds)

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


def parse_BCODESXP(reader,rds):
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

    value = _read_bc1(reader, reps,rds)
    return value


def _read_bc1(reader, reps,rds):
    """
    Python version of ReadBC1(ref_table, reps, stream).

    Layout:

        code   = ReadItem(...)
        consts = ReadBCConsts(...)
    """

    logging.debug("BCODESXP: reading bytecode 'code' object")
    code_robj = parse_object(reader,rds)

    logging.debug("BCODESXP: reading bytecode const pool")
    consts_list = _read_bc_consts(reader, reps,rds)  # list of Robj / nested BC

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


def _read_bc_consts(reader, reps,rds):
    n = struct.unpack(">i", reader.read(4))[0]
    logging.debug(f"BCODESXP consts: n = {n}")

    consts = []
    for i in range(n):
        type_code = struct.unpack(">i", reader.read(4))[0]
        logging.debug(f"BCODESXP const[{i}] type_code = {type_code}")

        if type_code == BCODESXP_CODE:
            logging.debug("BCODESXP const: nested BCODESXP")
            const_val = _read_bc1(reader, reps)
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
            const_val = _read_bc_lang(reader, type_code, reps,rds)
            consts.append(const_val)

        else:
            logging.debug("BCODESXP const: delegating to parse_object()")
            const_val = parse_object(reader,rds)
            consts.append(const_val)

    return consts


def _read_bc_lang(reader, type_code, reps,rds):
    """
    Python port of ReadBCLang(int type, SEXP ref_table, SEXP reps, R_inpstream_t stream).

    type_code is the integer 'type' that R read with InInteger().

    Semantics:

      - If type == BCREPREF:
            return reps[ InInteger(stream) ]

      - If type in {BCREPDEF, LANGSXP, LISTSXP, ATTRLANGSXP, ATTRLISTSXP}:
            build a LANG/LIST node with optional attributes and
            recursively-read CAR/CDR (which themselves are given by
            type codes and handled again via ReadBCLang).

      - Otherwise:
            default: just ReadItem(...) and return that SEXP
            (i.e. in Python, call parse_object(reader)).
    """

    # 1. BCREPREF: reference to an existing node in 'reps'
    if type_code == BCREPREF_CODE:
        idx = struct.unpack(">i", reader.read(4))[0]
        logging.debug(f"ReadBCLang: BCREPREF idx={idx}")
        if 0 <= idx < len(reps):
            return reps[idx]
        logging.warning(f"ReadBCLang: BCREPREF idx={idx} out of range")
        return None

    # 2. Default branch: NOT a BCLang-compressed node
    if type_code not in (
        BCREPDEF_CODE,
        LANGSXP_CODE,
        LISTSXP_CODE,
        ATTRLANGSXP_CODE,
        ATTRLISTSXP_CODE,
    ):
        logging.debug(f"ReadBCLang: default -> parse_object (type_code={type_code})")
        # This mirrors the 'default' case in C code (ReadItem).
        return parse_object(reader,rds)

    # 3. Now we're in the "real" BCLang path: BCREPDEF/LANG/LIST/ATTRLANG/ATTRLIST
    pos = -1
    hasattr = False
    t = type_code

    if t == BCREPDEF_CODE:
        # BCREPDEF: first int is position, second is underlying type
        pos = struct.unpack(">i", reader.read(4))[0]
        t = struct.unpack(">i", reader.read(4))[0]
        logging.debug(f"ReadBCLang: BCREPDEF pos={pos}, inner_type={t}")

    # Handle "attribute-wrapped" types.
    if t == ATTRLANGSXP_CODE:
        t = LANGSXP_CODE
        hasattr = True
    elif t == ATTRLISTSXP_CODE:
        t = LISTSXP_CODE
        hasattr = True

    # Human-readable type name for debugging / printing
    type_name = SEXPTYPEs.get(t, {}).get("SEXPTYPE", f"TYPE_{t}")

    node = Robj()
    node.attributes = []

    # Optional attributes for the bc-lang node
    attr_robj = None
    if hasattr:
        logging.debug("ReadBCLang: reading attributes")
        attr_robj = parse_object(reader,rds)

    logging.debug("ReadBCLang: reading TAG")
    tag_robj = parse_object(reader,rds)

    logging.debug("ReadBCLang: reading CAR type")
    car_type = struct.unpack(">i", reader.read(4))[0]
    car_val = _read_bc_lang(reader, car_type, reps,rds)

    logging.debug("ReadBCLang: reading CDR type")
    cdr_type = struct.unpack(">i", reader.read(4))[0]
    cdr_val = _read_bc_lang(reader, cdr_type, reps,rds)

    # Store this node in 'reps' if it was a BCREPDEF
    if 0 <= pos < len(reps):
        reps[pos] = node

    # Represent this bc-lang node as a pairlist-like structure so your
    # existing toString() can print it nicely.
    node.value = [
        ["bctype", type_name],
        ["tag", tag_robj],
        ["car", car_val],
        ["cdr", cdr_val],
    ]

    # Attach attributes if present.
    if attr_robj is not None:
        if isinstance(attr_robj.value, list):
            node.attributes = attr_robj.value
        else:
            node.attributes = [["attr", attr_robj]]

    return node
  
def parse_EXTPTRSXP(reader, rds):
    """
    Parse EXTPTRSXP (external pointer).

    According to R serialize.c:
      - pointer value itself is NOT serialized (restored as NULL)
      - prot, tag, and attributes ARE serialized

    Layout on disk:
        [EXTPTRSXP header]
        PROT  = ReadItem(...)
        TAG   = ReadItem(...)
        ATTR  = ReadItem(...)

    We represent this as:
        value = [
            ["ptr",  None],
            ["prot", <Robj>],
            ["tag",  <Robj>],
        ]
    """

    logging.debug("EXTPTRSXP: parsing prot")
    prot = parse_object(reader, rds)

    logging.debug("EXTPTRSXP: parsing tag")
    tag = parse_object(reader, rds)

    logging.debug("EXTPTRSXP: parsing attributes")
    attrs = parse_object(reader, rds)

    value = [
        ["ptr", None],   # pointer is always NULL after unserialization
        ["prot", prot],
        ["tag", tag],
    ]

    return value, attrs


# end of BCODESXP ####

# helper functions to convert to python types ###############
def as_data_frame(robj):
  """
  Converts Robj to pandas DataFrame

  Parameters
  ----------
  robj : Robj

  Returns
  -------
  DataFrame

  Examples
  --------
  as_data_frame(robj)
  """
  cols = {}
  names = robj.get('names').value
  for i in range(len(names)):
    val = robj.value[i]
    cl = val.get('class')
    if (cl is not None) and 'factor' in cl.value:
      levels = val.get('levels').value
      val = [levels[i-1] for i in val.value]
    else:
      val = val.value
    cols[names[i]] = val
    
  r = pd.DataFrame(cols)
  r.index = robj.get('row.names').value
  return r

def as_numpy(robj):
  """
  Converts Robj with array/Matrix into dense/sparse numpy array
  
  Parameters
  ----------
  robj : Robj

  Returns
  -------
  np.array

  Examples
  --------
  as_numpy(robj)
  """
  cl = robj.getClass()
  if cl == 'dgCMatrix':
    return _dgCMatrix2numpy(robj)
  if cl == 'REALSXP':
    return _array2numpy(robj)

def as_anndata(robj):
  """
  Converts Robj with array/Matrix into dense/sparse anndata keeping dimnames if any
  
  Parameters
  ----------
  robj : Robj

  Returns
  -------
  ad.AnnData

  Examples
  --------
  as_anndata(robj)
  """
  X = as_numpy(robj).T
  adata = ad.AnnData(X=X)
  dimnames = robj.get('dimnames')
  if( dimnames is None):
    dimnames = robj.get('Dimnames') # for sparse Matrices 
  
  if dimnames is not None:
    dimnames = dimnames.value
    if dimnames[0].value is not None:
      adata.var_names = dimnames[0].value
    if dimnames[1].value is not None:
      adata.obs_names = dimnames[1].value
  return adata
    
def seurat2adata(robj,assay=0,layer='counts'):
  """
  Converts Seurat Robj into anndata
  it loads:
  1. specified assay/layer as data
  2. cell metadata
  3. var metadata if any
  4. all reduced dimensions
  
  Parameters
  ----------
  robj : Robj or str (path 2 rds file)
  assay : int - assay index
  layer : str - named of layer to use

  Returns
  -------
  AnnData

  Examples
  --------
  seurat2adata(robj)
  """
  if isinstance(robj,str):
    robj = parse_rds(robj)
  cnts = robj.get(['assays',assay,layer])
  obs = as_data_frame(robj.get('meta.data'))
  
  if cnts is not None:
    adata = as_anndata(cnts)
    adata.var = as_data_frame(robj.get(['assays',assay,'meta.features']))
    adata.obs = obs.loc[adata.obs_names,:]
  # try Assay5
  else:
    names = robj.get(['assays',assay,'layers','names']).value
    layer = np.where(names == layer)[0]
    cnts  = robj.get(['assays',assay,'layers',layer])
    adata = as_anndata(cnts)
    adata.var_names = robj.get(['assays',0,'features','dimnames',0]).value
    adata.obs = obs
  # load reduced dims
  rdims = robj.get(['reductions'])
  rdims_names = rdims.get('names')
  if rdims_names is not None:
    rdims_names = rdims_names.value
    for i in range(len(rdims_names)):
      adata.obsm[rdims_names[i]] = _array2numpy(rdims.get([i,'cell.embeddings']))
  
  return adata

def _array2numpy(robj):
  """
  Converts Robj with array into dense numpy array
  Expect R array (dense)
  
  Parameters
  ----------
  robj : Robj

  Returns
  -------
  np.array

  Examples
  --------
  _array2numpy(robj)
  """
  dim = robj.get('dim')
  if dim is None:
    dim = len(robj.value)
  else:
    dim = dim.value
  X = np.array(robj.value).reshape(dim,order='F')
  return X
  
def _dgCMatrix2numpy(robj):
  i = robj.get('i').value
  p = robj.get('p').value
  dim = robj.get('Dim').value
  x = robj.get('x').value
  dimnames = robj.get('Dimnames').value
  X = sparse.csc_matrix((x, i, p), dim)
  return X



    
