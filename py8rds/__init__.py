import logging
from importlib.metadata import version

__version__ = version("py8rds")

logging.basicConfig(level="DEBUG", format="[%(asctime)s][%(levelname)s] %(message)s")
logging.getLogger().setLevel("INFO")

from .robj import INT_NA, Robj, RdsFile
from .parser import parse_rds, parse_object
from .convert import (
    factor2array,
    as_data_frame,
    as_numpy,
    as_anndata,
    as_dict,
    seurat2adata,
    seurat2adata_spatial,
    is_default_index,
)
