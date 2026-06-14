![logo](assets/logo.png)
## a Python package to parse and manipulate the BRENDA database

[![tests](https://github.com/Robaina/BRENDApyrser/actions/workflows/tests.yml/badge.svg)](https://github.com/Robaina/BRENDApyrser/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/Robaina/BRENDApyrser/graph/badge.svg?token=214SPFXRTG)](https://codecov.io/gh/Robaina/BRENDApyrser)
![PyPI](https://img.shields.io/pypi/v/brendapyrser)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/Robaina/Brendapyrser)
[![GitHub license](https://img.shields.io/github/license/Robaina/BRENDApyrser)](https://github.com/Robaina/BRENDApyrser/blob/master/LICENSE)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)
[![DOI](https://zenodo.org/badge/299416438.svg)](https://zenodo.org/badge/latestdoi/299416438)

## What is Brendapyrser?
This project provides python classes and functions to parse the JSON file containing the entire BRENDA enzyme database (https://www.brenda-enzymes.org)

> **Note**: Since the 2024.x releases, BRENDA is distributed as a structured JSON document
> (schema [2.0.0](https://www.brenda-enzymes.org/schemas/docs/2.0.0/brenda.schema.html))
> rather than the legacy flat text file. BRENDApyrser `>= 0.1.0` parses this JSON format. For
> the old text format, use BRENDApyrser `0.0.4`.

This is an ongoing project!

## Installation
1. ```pip install brendapyrser```

or

2. Git clone project to local directory.

   In terminal navigate to directory and enter: ```python setup.py install```

## Usage

Due to BRENDA's license, BRENDA's database cannot be downloaded directly by the parser, instead, the user is asked to download the database as a JSON file after accepting usage conditions [here](https://www.brenda-enzymes.org/download.php).

The parser accepts the JSON file directly, as well as `.json.gz` and the `.json.tar.gz`
archive in which BRENDA ships the database — no manual extraction required:

```python
from brendapyrser import BRENDA

brenda = BRENDA("brenda_2026_1.json.tar.gz")

reaction = brenda.reactions.get_by_id("1.1.1.1")
print(reaction.name)            # alcohol dehydrogenase
print(reaction.organisms[:3])   # source organisms
print(reaction.KMvalues)        # KM values keyed by substrate
print(reaction.bibliography)    # structured references (PMID, authors, year, ...)
```

You can find a jupyter notebook with usage examples [here](docs/examples.ipynb).

## Contribute

Contributions are always more than welcome! Feel free to fork, make changes and pull requests. If you don't know where to start, you can always browse open issues. Thanks!

## Citation

If you use this software, please cite it as below:

Robaina-Estévez, S. (2026). BRENDApyrser: a Python package to parse and manipulate the BRENDA database (Version 0.1.0)[Computer software]. https://doi.org/10.5281/zenodo.7026555
