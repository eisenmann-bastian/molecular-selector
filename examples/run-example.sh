#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
if [ "$#" -gt 1 ]; then
    echo "Usage: $0 [output-directory]" >&2
    exit 2
fi

venv="$root/.venv"
if [ ! -e "$venv" ] && [ ! -L "$venv" ]; then
    "$root/examples/setup-venv.sh" "$venv"
fi
if [ ! -f "$venv/bin/activate" ] || [ ! -x "$venv/bin/python" ]; then
    echo "Invalid virtual environment at $venv; inspect it before retrying." >&2
    exit 1
fi
. "$venv/bin/activate"
python="$venv/bin/python"

if [ "$#" -eq 1 ]; then
    mkdir -p -- "$1"
    output_dir=$(CDPATH= cd -- "$1" && pwd)
else
    output_dir=$(mktemp -d "${TMPDIR:-/tmp}/molecular-selector-example.XXXXXXXX")
fi

cd "$root"
"$python" -m molecular_selector examples/molecules/ethanethiol.sdf \
    --name thiol -o "$output_dir/thiol.csv"
"$python" -m molecular_selector examples/molecules/acetic-acid.sdf \
    --smarts '[CX3](=O)[OX2H1]' -o "$output_dir/custom-smarts.csv"
"$python" -m molecular_selector --csv examples/selected-molecules.csv \
    -o "$output_dir/selected.csv"
"$python" -m molecular_selector --csv examples/molecules.csv \
    -o "$output_dir/all-selectors.csv"
printf 'Example output: %s\n' "$output_dir"
