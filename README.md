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
2. `get([inx1,key2])` subsets object by keys.
Each Robj has values that are indexed by integers (shown as `+N` by `show` function) and slotes/attributes that are indexed by keys (shown as `&/*<key>` by `show` function)
