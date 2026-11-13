# Molecular selector

Identify functional groups and interaction sites in molecular SDF records using
RDKit. Selection returns the matched atom indices and highlighted site indices;
axis generation and molecular alignment belong to molecular-aligner.

## Installation

Install Git and Python 3.10 or later, then install from this repository:

```sh
python -m pip install .
```

The installation fetches `molecular-utilities` and a pinned revision of
`molecular-functional-groups` from GitHub as specified in `pyproject.toml`.
No sibling checkout is required. The distribution and CLI
command are named `molecular-selector`; the old `molecule-selector` command
is no longer provided.

For development, use `python -m pip install -e .`, then run
`python -m unittest discover -s tests` for the offline tests.

## Examples

The [examples](examples/) directory contains sample SDFs, a one-column CSV
that runs all selectors, and a CSV that chooses selectors per molecule. After
installation, run the examples or the tests from any working directory:

```sh
./examples/run-example.sh
./examples/run-example.sh /path/to/new-output-directory
./examples/run-tests.sh
```

The example runner prints its output directory and writes CSVs for a named
selector, a custom SMARTS query, selected batch selectors, and all selectors.
Without an argument it creates a temporary output directory. With an argument,
existing output files are not overwritten. Set `PYTHON` to select the
interpreter where the package and its dependencies are installed.

## Usage

```sh
molecular-selector --list-selectors
molecular-selector molecule.sdf --name thiol -o matches.csv
molecular-selector molecule.sdf --name "primary amide" -o amides.csv
molecular-selector molecule.sdf --name "interaction halogen bond sigma hole donor" -o halogens.csv
molecular-selector --csv molecules.csv -o matches.csv
```

The batch CSV contains one SDF path per row. Relative SDF paths are resolved
relative to the CSV file. Output SDF paths are relative to the working directory.
The output includes zero-based record, match and atom indices; atom index arrays
are encoded as JSON. Existing output files are not overwritten.

## Shared definitions and interaction sites

The named functional-group selectors are built from the atom-mapped patterns
in `molecular-functional-groups`. Names accept spaces, underscores or hyphens.
Every shared interaction/role combination is available under
`interaction <interaction type> <role>`, as shown by `--list-selectors`;
existing names such as `hydrogen bond donor` are also available. One-column
batch input selects **all** named selectors, including the newly added ones.

Site indices for shared interaction selectors come from the role-specific map
labels, not all labeled atoms. A role can have several group queries. For a
multi-query selector, the Python result's `selector.smarts` is a JSON array
encoded as a string, and each matching entry has a `smarts` field identifying
the query that produced it. The CSV `smarts` column contains that individual
query for a match, or the JSON array if no query matched. Results with the same
set of site atoms are deduplicated in deterministic query order.

This switches applicable existing selectors to the shared definitions and
therefore changes matches and atom order. For example, shared ether excludes
esters, the shared carboxylate pattern can match carbonate, and a halogen
donor can include fluorine. The shared aldehyde pattern excludes formaldehyde;
the shared primary-alcohol pattern excludes methanol;
the shared chalcogen donor patterns cover sulfur groups rather than the
previous selenium/tellurium query. `atom_indices` follow each shared query's
atom order. The standalone `carboxylate_sites` helper retains its narrower
carboxylate definition for compatibility with alignment.

Matches describe **potential** sites, not confirmed interactions. The shared
patterns and role mappings have limitations: protonation, resonance, local
environment and geometry matter, and not every mapped atom is necessarily a
valid interaction partner. In particular, atom-level acceptor checks can be
needed for phosphate and amine groups. Selector-specific `carbonyl`, charge,
metal-center and aromatic-ring selectors remain available because the shared
package does not provide equivalent definitions.

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
