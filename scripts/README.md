## AFLOW Dev Scripts

### Linting

- `format-changed.sh`: Runs clang-format on changed lines since the `staging` branch.
- `tidy-changed.sh`: Runs clang-tidy using `.clang-tidy-lint` on changed lines since 
the `staging` branch.

Both of these scripts should be run from the repo root as the working directory (i.e. run 
as `scripts/format-changed.sh` and `scripts/tidy-changed.sh`).

You can install the clang-format and clang-tidy tools easily using the [./dev/requirements.txt](./dev/requirements.txt)
file. Setup a venv and just run `pip install -r scripts/dev/requirements.txt`.
