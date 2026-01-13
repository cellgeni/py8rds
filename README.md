# 🐍 🍽️ RDS

`py8rds` *(python ate RDS)* provides pure-Python deserialization of R `.rds` files, allowing you to load R data directly into Python without requiring an R installation.

## Getting Started

To get a local copy up and running follow these steps.

### Prerequisites

- Python 3.8+

### Installation

Use pip to install directly from the repo
```bash
pip install git+hhttps://github.com/cellgeni/py8rds.git
```

### Development Installation

1. Clone the repo
   ```bash
   git clone hhttps://github.com/cellgeni/py8rds.git
   ```
2. Usse pip to install in [editable mode](https://setuptools.pypa.io/en/latest/userguide/development_mode.html)
   ```bash
   pip install --editable ./py8rds
   ```

## Usage

```python
import rdd

result = rdd.parse_rds("myfile.rds")

print(result.values)
```


### Testing

From the package folder, run:
```bash
python -m unittest tests
```

Test data is created using `tests/create_test_data.R`  which requires a working R installation.
Seurat object creation requires [Seurat installed](https://satijalab.org/seurat/articles/install.html).
