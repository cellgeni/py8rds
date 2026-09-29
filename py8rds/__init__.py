import logging
from importlib import import_module
from importlib.metadata import version

__version__ = version("py8rds")

logging.basicConfig(level="DEBUG", format="[%(asctime)s][%(levelname)s] %(message)s")
logging.getLogger().setLevel("INFO")

from .robj import INT_NA, Robj, RdsFile
from .parser import parse_rds, parse_object

# converters depend on pandas/anndata/scipy which are slow to import, so they are loaded on first use
_CONVERTERS = {
    "factor2array",
    "as_data_frame",
    "as_numpy",
    "as_anndata",
    "as_dict",
    "seurat2adata",
    "seurat2adata_spatial",
    "is_default_index",
}

__all__ = ["INT_NA", "Robj", "RdsFile", "parse_rds", "parse_object", *sorted(_CONVERTERS)]


def __getattr__(name):
    if name in _CONVERTERS:
        value = getattr(import_module(".convert", __name__), name)
        globals()[name] = value  # cache, so __getattr__ is not called again
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | _CONVERTERS)
