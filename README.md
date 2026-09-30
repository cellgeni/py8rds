# 🐍 🍽️ RDS

`py8rds` *(python ate RDS)* provides pure-Python deserialization of R `.rds` and `.qs2` files, allowing you to load R data directly into Python without requiring an R installation.


## Prerequisites

- Python 3.8+

## Installation

Use pip to install directly from the repository.
```bash
pip install git+https://github.com/cellgeni/py8rds.git
```

## Usage

```python
import py8rds

# to load data.frames save
df = py8rds.as_data_frame("data_frame.rds")

# to load dense/sparse named matrices; anndata allows to keep dimnames.
mtx = py8rds.as_anndata("matrix.rds")

# to read Seurat. Please note that it reads only one assay:layer (first assay and `counts` layer by default, see assay and layer parameters)
adata = py8rds.seurat2adata("seurat.rds")

# general function to read rds into python
robj = py8rds.parse_rds("data.rds")
robj.show(level=2)
# you can use converters on robj once you know what is in:
metadata = py8rds.as_data_frame(robj.get(["meta.data"]))
```
All functions also accept [qs2](https://github.com/qsbase/qs2) files saved by `qs2::qs_save` (file format is detected automatically). Files saved by `qs2::qd_save` (qdata format) or by the legacy `qs::qsave` are not supported.
Please see the [tutorial](tutorials/tutorial.ipynb).


## Details

`parse_rds` is the base function that reads an RDS/qs2 file into a Python `Robj`, which has a tree-like structure. `Robj` has two main methods:
1. `show(level=1)` shows the object structure to the specified level.
2. `get([inx1,key2])` recursively subsets the object by the provided keys/indices. Returns Robj or None.
Each `Robj` has values that are indexed by integers (shown as `+N` by `show` function) and slots/attributes that are indexed by keys (shown as `&/*<key>` by `show` function).

All convertor functions (such as `as_data_frame`, `as_numpy`,`as_anndata`, `seurat2adata`, `seurat2adata_spatial`) can take as input both, file name and Robj, so if you are unsure about rds/qs2 file content you may first load it with `parse_rds`, browse and subset by `show` and `get` and then convert. This approach can save time on file reading.

Convertor functions named `as_` are designed to keep all content of rds/qs2 in python representation. Convertors like `seurat2adata` only keeps some data, as Seurat object cannot be mapped completely into anndata.


## Acknowledgements
[Amazing blog](https://blog.djnavarro.net/posts/2021-11-15_serialisation-with-rds/) by Danielle Navarro helped us kick-start the project, other projects such as [rds2cpp](https://github.com/LTLA/rds2cpp/blob/master/include/rds2cpp/parse_object.hpp) helped us move forward, [R source code](https://github.com/wch/r-source/blob/trunk/src/main/serialize.c) became the last resort after meeting with the [R documentation](https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats), and at the very end, ChatGPT came to our rescue.

## Similar projects
1. [rds2py](https://github.com/BiocPy/rds2py), based on rds2cpp, cannot read functions so fails on complex objects such as Seurat
2. [pyreadr](https://github.com/ofajardo/pyreadr), focused on simple data types such as data.frames, cannot read complex objects such as Seurat.
3. [rdata](https://github.com/vnmabus/rdata) in addition to reading rds files it can also save python objects into rds, but it fails to read Seurat objects.

So, none of the alternatives seem able (at least at the moment) to read Seurat objects (see this [notebook](tutorials/alternatives.ipynb)).
