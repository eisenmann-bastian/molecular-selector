# Project instructions

- Keep functional-group patterns, matching, interaction-site definitions,
  batch handling, and output in their existing `molecular_selector/` modules.
  Preserve public selector names and documented match semantics.
- Configure `molecular_utilities.logging` only at the CLI entry point using
  the `molecular_selector` namespace; avoid import-time configuration.
- Add tests and update examples and the [README](README.md) when behavior
  changes. Run `./examples/run-tests.sh` in an installed environment.
- Preserve existing output files and virtual environments.
