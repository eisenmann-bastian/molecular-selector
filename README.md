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
`pyproject.toml`. Functional-group definitions are included in this package;
no separate functional-groups installation or sibling checkout is required.
The distribution and CLI command are named `molecular-selector`; the old
`molecule-selector` command is no longer provided.

For development, use `python -m pip install -e .`, then run
`python -m unittest discover -s tests` for the offline tests.

## Examples

The [examples](examples/) directory contains sample SDFs, a one-column CSV
that runs all selectors, and a CSV that chooses selectors per molecule. Run
the examples from the project directory:

```sh
./examples/run-example.sh
./examples/run-example.sh /path/to/new-output-directory
```

The example runner creates the project's `.venv` and installs this checkout
in editable mode if the environment does not exist. It activates the
environment within the script; it does not change the calling shell. If
`.venv` already exists, it is used as-is rather than modified; an incomplete
or uninstalled environment must be repaired or recreated manually. Python
3.10 or later with `venv`, Git and network access are required for initial
setup. Set `PYTHON` to select the interpreter used for that creation.

To set up an environment explicitly and run the tests:

```sh
./examples/setup-venv.sh                  # creates .venv in the project
# Or choose a different, new path:
./examples/setup-venv.sh /path/to/selector-venv
. .venv/bin/activate                     # or activate the chosen path
./examples/run-tests.sh
```

The setup script installs dependencies from `pyproject.toml`, including the
GitHub dependencies. It refuses to modify an existing path; if installation
fails after creation, inspect or remove that environment before retrying.
The test runner uses the active interpreter or the interpreter specified by
`PYTHON`. The example runner always uses the project `.venv`, even if another
environment is active.

The example runner prints its output directory and writes CSVs for a named
selector, a custom SMARTS query, selected batch selectors, and all selectors.
Without an argument it creates a temporary output directory. With an argument,
existing output files are not overwritten. Set `PYTHON` to select the
interpreter used to create the project environment when it is missing.

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

## Functional groups and interaction sites

The `molecular_selector` package exposes functional-group SMARTS constants,
`FunctionalGroup`, `InteractionType`, `InteractionRole`, and their mappings.
Its public functions include
`get_smarts_pattern`, `get_groups_for_role` and `get_interaction_site`:

```python
from molecular_selector import (
    FunctionalGroup, InteractionRole, InteractionType, MolecularSelector,
    get_groups_for_role, get_interaction_site, get_smarts_pattern,
)
from rdkit import Chem

pattern = get_smarts_pattern(FunctionalGroup.CARBOXYLATE)
donors = get_groups_for_role(InteractionType.HYDROGEN_BOND, InteractionRole.DONOR)
site = get_interaction_site("hydrogen_bond", "donor", "phenol")
print(pattern, donors, site.smarts, site.atom_map_numbers)
result = MolecularSelector().select_molecules(
    [Chem.MolFromSmiles("Oc1ccccc1")], name=FunctionalGroup.PHENOL
)
```

`FUNCTIONAL_GROUP_SMARTS` maps groups to patterns, and
`MAPPED_FUNCTIONAL_GROUP_SMARTS` supplies atom-mapped versions.
`INTERACTION_ROLES` (also available as `BOND_TYPE_ROLES`) maps interactions
and roles to groups; `INTERACTION_SITES` adds each group's role-specific
query and site labels. Atom-map numbers are local query labels, not molecule
atom indices. A role-specific query may be more restrictive than the group
query, for example an ammonium hydrogen-bond donor must have a hydrogen.
Unknown enum values raise `ValueError`; `get_interaction_site` raises
`KeyError` for a known but unmapped combination. `get_groups_for_role`
returns an empty set for an unmapped role.

Code previously importing the standalone package or the interim
`molecular_selector.functional_groups` subpackage must import these names from
`molecular_selector` instead. The standalone package is not installed or
replaced by this distribution.

Named selectors are derived from these mappings; compatibility names and
composite selectors are handled by the selector. Names accept spaces,
underscores or hyphens. Every interaction/role combination is available under
`interaction <interaction type> <role>`. Where an existing short name such
as `hydrogen bond donor` exists, the longer name is its alias and results use
the short canonical name. `--list-selectors` lists canonical names.
One-column batch input selects each canonical selector once.

For a named functional group, highlighted sites are the union of its
hydrogen-bond donor and acceptor atom labels. If neither role is defined, the
entire matched group is highlighted. Group selection uses the general group
pattern, not any more restrictive role-specific query; to select sites for a
specific role, use an interaction selector instead. This derivation can
change site indices for groups that previously highlighted the whole match.

Site indices for interaction selectors come from the role-specific map
labels, not all labeled atoms. A role can have several group queries. For a
multi-query selector, the Python result's `selector.smarts` is a JSON array
encoded as a string, and each matching entry has a `smarts` field identifying
the query that produced it. The CSV `smarts` column contains that individual
query for a match, or the JSON array if no query matched. Results with the same
set of site atoms are deduplicated in deterministic query order.

These patterns can differ from older selector releases. For example, ether
excludes esters, carboxylate can match carbonate, and a halogen donor can
include fluorine. Aldehyde excludes formaldehyde; primary alcohol excludes
methanol; chalcogen donor patterns cover sulfur groups rather than the
previous selenium/tellurium query. `atom_indices` follow each query's
atom order. The standalone `carboxylate_sites` helper retains its narrower
carboxylate definition for compatibility with alignment.

Matches describe **potential** sites, not confirmed interactions. The
patterns and role mappings have limitations: protonation, resonance, local
environment and geometry matter, and not every mapped atom is necessarily a
valid interaction partner. In particular, atom-level acceptor checks can be
needed for phosphate and amine groups. Carbonyl, charge and metal-center
selectors use general site patterns from this package. Aromatic-ring
selection remains selector-specific to support aromatic heterocycles and
fused rings.

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
