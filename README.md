# Molecular selector

Identify functional groups and interaction sites in molecular SDF records using
RDKit. Selection returns the matched atom indices and highlighted site indices;
axis generation and molecular alignment belong to molecular-aligner.

## Installation

Install Git and Python 3.10 or later, then install from this repository:

```sh
python -m pip install .
```

The installation fetches `molecular-utilities` from GitHub as specified in
`pyproject.toml`. No sibling checkout is required. The distribution and CLI
command are named `molecular-selector`; the old `molecule-selector` command
is no longer provided.

For development, use `python -m pip install -e .`, then run
`python -m unittest discover -s tests` for the offline tests.

## Usage

```sh
molecular-selector --list-selectors
molecular-selector molecule.sdf --name thiol -o matches.csv
molecular-selector --csv molecules.csv -o matches.csv
```

The batch CSV contains one SDF path per row. Relative SDF paths are resolved
relative to the CSV file. Output SDF paths are relative to the working directory.
The output includes zero-based record, match and atom indices; atom index arrays
are encoded as JSON. Existing output files are not overwritten.

The `molecular_selector` Python package exposes `MolecularSelector`,
`SelectorRegistry`, `BatchSelector`, `CsvResultWriter`, and `carboxylate_sites`.
For example:

```python
from molecular_selector import MolecularSelector
```

`MolecularSelector.read_sdf(path)`
provides a public iterator over SDF records with explicit hydrogens retained.
Selection CSVs can be consumed by `molecular-axes` from molecular-aligner.

The CLI configures its shared utilities logger under the `molecular_selector`
namespace and writes log messages to standard error; CSV output remains on
standard output. Importing the library does not configure logging.
