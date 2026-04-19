# Contributing

This page summarizes the repository contribution workflow. The authoritative source is the top-level `CONTRIBUTING.md`, but the key expectations are collected here for easier navigation inside the docs site.

## Basic workflow

1. Fork the repository.
2. Clone your fork.
3. Initialize and update submodules.
4. Create a feature branch from the active development branch.
5. Open a pull request once the change is documented and tested.

## Expectations for changes

- Keep changes focused.
- Follow the surrounding Fortran and Python style.
- Update documentation when behavior or public interfaces change.
- Add or update tests when new functionality is introduced.

## Before opening a pull request

- Build the project with CMake.
- Run the relevant compiled tests and Python tests.
- Explain the motivation, scope, and validation of the change in the PR description.

## Review policy

Pull requests are expected to be reviewed and approved before merge. The repository also requires extra human review for significant AI-generated contributions, as described in the root contribution guide.
