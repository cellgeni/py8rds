# 🐍 🍽️ RDS

`py8rds` *(python ate RDS)* provides pure-Python deserialization of R `.rds` files, allowing you to load R data directly into Python without requiring an R installation.


## Prerequisites

- Python 3.8+

## Installation

Use pip to install directly from the repo
```bash
pip install git+https://github.com/cellgeni/py8rds.git
```

## Usage

```python
import py8rds
df = py8rds.as_data_frame('data_frame.rds')
adata = py8rds.as_anndata('seurat.rds')
robj = py8rds.parse_rds('data.rds')
robj.show()
```
Please check the [tutorial](tutorials/tutorial.ipynb). 


## Details

`parse_rds` is base function that reads rds into python Robj object that has tree-like structure. Robj has two main functions:
1. `show(level=1)` shows object structure to specified level. 
2. `get([inx1,key2])` recursively subsets object by provided keys/indeces.
Each Robj has values that are indexed by integers (shown as `+N` by `show` function) and slotes/attributes that are indexed by keys (shown as `&/*<key>` by `show` function)


## Acknowledgements
[Amazing blog](https://blog.djnavarro.net/posts/2021-11-15_serialisation-with-rds/) by Danielle Navarro allowed us to kick start the project, alternative projects such as [rds2cpp](https://github.com/LTLA/rds2cpp/blob/master/include/rds2cpp/parse_object.hpp) helped to move forward, [r source code](https://github.com/wch/r-source/blob/trunk/src/main/serialize.c) became the last resort after meeting with [R documentation](https://cran.r-project.org/doc/manuals/r-release/R-ints.html#Serialization-Formats) and at the very end chatGPT came to our rescue.

## Similar projects
1. [rds2py](https://github.com/BiocPy/rds2py), based on rds2cpp, cannot read functions so fails on complex objects such as Seurat
2. [pyreadr](https://github.com/ofajardo/pyreadr), focused on simple data types such as data.frames, cannot read complex objects such as Seurat.
3. [rdata](https://github.com/vnmabus/rdata) in addition to reading rds files is can also save python objects into rds, but it fails to read Seurat objects.
So, neither of alternative seems to be able (at least at given moment) to read Seurat object (see this [notebook](tutorials/alternatives.ipynb)).