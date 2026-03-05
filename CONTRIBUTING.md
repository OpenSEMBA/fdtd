# Contributing to OpenSEMBA / fdtd

Thank you for your interest in contributing to the fdtd solver!
This document describes the basic workflow and a few conventions to
follow when contributing.

## Code of conduct

Be respectful and constructive when interacting in issues, pull
requests, and reviews. Harassment or discrimination is not tolerated.

## How to get started

- Fork the repository on GitHub.
- Clone your fork and add the main repository as an `upstream` remote.
- Make sure submodules are initialized:
  - `git submodule init`
  - `git submodule update`
- Create a feature branch from the latest `main` or `dev` branch
  (see the repository README for current workflow).

Example:

```bash
git clone https://github.com/<your-username>/fdtd.git
cd fdtd
git remote add upstream https://github.com/OpenSEMBA/fdtd.git
git checkout -b my-feature-branch upstream/dev
```

Once you have prepared your branch, create a pull request (PR) to `dev`.
Before being merged PRs must:

- Be reviewed and approved by at least one of the top 3 contributors.
- Pass all existing tests.

## Development environment

The project uses CMake and Fortran (with optional MPI, HDF5 and MTLN
support). For details about compilation, pre‑compiled libraries and
platform‑specific notes, refer to:

- `doc/development.md`

In short:

- Always update submodules before configuring CMake.
- Use a separate `build/` directory (as in the examples in
  `doc/development.md`).
- Prefer reproducible build configurations by passing the same CMake
  options you expect CI to use.

### Python tools and tests

If you work with the Python wrapper or Python tests, it is recommended
to use a virtual environment instead of the system Python:

```bash
cd PATH_TO_FOLDER/fdtd  # or your fdtd clone
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Then run tests with:

```bash
pytest test/
```

(You can also use `python -m pytest test/`.)

## Making changes

- Keep changes focused and as small as reasonably possible.
- Follow the existing code style of the surrounding Fortran and Python
  code instead of introducing new styles.
- Update or add documentation in `doc/` when behaviour or usage
  changes.
- When modifying public interfaces (input formats, command line
  options, etc.), please describe the change clearly in your pull
  request.

## Testing your changes

Before opening a pull request:

- Build the project using CMake (see `doc/development.md`).
- Run unit tests, if they apply to your changes. For example:
  - `build/bin/fdtd_tests` (depending on your setup).
  - `pytest test/` for Python tests.
- If you add new functionality, add or update tests when possible.

## Commit and pull request guidelines

- Write clear, descriptive commit messages in English.
- Squash small fixup commits when possible, especially before merging.
- Open a pull request against the `dev` or `main` branch, following
  the project’s current workflow.
- In the pull request description, include:
  - A short summary of the change.
  - Motivation / context (what problem it solves or what feature it
    adds).
  - How it was tested (commands run, platforms used).
- Be responsive to review comments and keep the discussion technical
  and respectful.

## Reporting bugs and requesting features

If you are not yet ready to contribute code, you can still help by:

- Opening detailed bug reports with:
  - Steps to reproduce.
  - Expected vs. actual behaviour.
  - Platform, compiler, and relevant CMake options.
- Proposing enhancements or new features, including motivation and
  potential use cases.

Thank you again for contributing to fdtd!