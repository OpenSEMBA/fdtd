---
name: code-reviewer
description: PR review checklist for semba-fdtd: CI verification across Ubuntu/Windows matrix, Fortran code review patterns, conditional compilation guards, test coverage expectations, MPI/OpenMP safety checks, and AI-assisted contribution policies
---

## When to use

- Reviewing a pull request before merging
- Evaluating code quality and correctness of Fortran/Python changes
- Checking that tests cover new functionality
- Verifying conditional compilation guards are correct
- Ensuring documentation is updated for user-facing changes
- Reviewing AI-assisted contributions for correctness

## Key Files to Read

### Project Standards
- `CONTRIBUTING.md` — Contribution guidelines, review requirements
- `AGENTS.md` — AI coding agent guidance (this file)
- `README.md` — Project overview and quality expectations

### CI/CD Pipelines
- `.github/workflows/ubuntu.yml` — Ubuntu CI: GCC + Intel, MPI ON/OFF, MTLN ON/OFF
- `.github/workflows/windows.yml` — Windows CI: IntelLLVM + Ninja
- `.github/workflows/check-submodules-default-branch.yml` — Submodule integrity checks

### Code Style References
- `CMakeLists.txt` — Compiler flags, build configuration
- `doc/fdtdjson.md` — Input format spec (check for input format changes)
- `doc/development.md` — Development setup and conventions

## What CI Checks Must Pass

Every PR must pass **all** of the following before merging:

### Ubuntu CI (`ubuntu.yml`)
14 build configurations in matrix:
- Compilers: GCC 11 + Intel 2025.1
- MPI: ON and OFF
- MTLN: ON and OFF
- HDF: ON (always)
- Double precision: OFF (default) + 1 Intel config with ON

Each configuration must:
- Build successfully
- Pass `build/bin/fdtd_tests` (unit tests)
- Pass `pytest test/ --durations=20` (integration tests)

### Windows CI (`windows.yml`)
6 build configurations:
- Compiler: IntelLLVM 2025.2.1 + Ninja
- MTLN: ON and OFF
- HDF: ON
- Double precision: ON and OFF

Each configuration must:
- Build successfully
- Pass `build/bin/fdtd_tests.exe`
- Pass `pytest -m 'not codemodel' test/ --durations=20` (codemodel excluded on Windows)

### Submodule Check
- All submodule commits must be ancestors of their `main`/`master` branch history

## Review Checklist

### 1. General PR Quality

- [ ] PR description includes a short summary of the change
- [ ] PR description includes motivation/context (what problem it solves)
- [ ] PR description explains how it was tested (commands run, platforms used)
- [ ] PR is opened against `dev` (or `main` per current workflow)
- [ ] Commit messages are clear, descriptive, and in English
- [ ] Commit messages reference issue numbers where applicable (e.g., `#392`)
- [ ] Changes are focused and small (not mixing unrelated changes)
- [ ] Small fixup commits are squashed

### 2. Fortran Code Review

- [ ] **Code style matches surrounding code** — no new indentation, naming, or structural conventions introduced
- [ ] **Conditional compilation guards are correct** — check that `#ifdef` guards match the CMake flags:
  - `SEMBA_FDTD_ENABLE_MPI` — MPI support
  - `SEMBA_FDTD_ENABLE_HDF` — HDF5 output
  - `SEMBA_FDTD_ENABLE_MTLN` — MTLN solver + ngspice
  - `SEMBA_FDTD_ENABLE_SMBJSON` — JSON input parser
  - `SEMBA_FDTD_ENABLE_DOUBLE_PRECISION` — 8-byte reals
- [ ] **Module dependencies are respected** — code doesn't depend on modules below it in the library chain:
  ```
  semba-types -> semba-reports -> smbjson/conformal/semba-components
    -> mtlnsolver/semba-outputs -> semba-main -> semba-fdtd
  ```
- [ ] **No new unused variables or dead code**
- [ ] **Array bounds are safe** — FDTD uses extensive pointer aliasing; verify no out-of-bounds access
- [ ] **Memory management is correct** — no leaks in derived type allocation/deallocation
- [ ] **OpenMP directives are correct** — if adding parallel loops, verify `COLLAPSE`, `DEFAULT(SHARED)`, and private variables
- [ ] **MPI communication is safe** — if modifying MPI code, verify non-blocking patterns and buffer sizes

### 3. Build System Review

- [ ] **CMakeLists.txt changes are minimal** — only add new targets, don't restructure existing ones
- [ ] **New optional features have CMake options** — follow the pattern of existing flags
- [ ] **Compiler flags are correct** — match existing patterns for the target compiler
- [ ] **Windows compatibility** — if changes affect Windows, verify IntelLLVM + Ninja builds
- [ ] **Submodule changes are intentional** — verify submodule commits are on accepted branch history

### 4. Test Coverage

- [ ] **New functionality has corresponding tests** — both unit tests and integration tests where applicable
- [ ] **Bug fixes have regression tests** — a test that would have caught the original bug
- [ ] **Test data is appropriate** — minimal but sufficient to exercise the code path
- [ ] **Tests use correct pytest markers** — `mtln`, `hdf`, `mpi`, `codemodel`
- [ ] **Tests are isolated** — Python tests use `tmp_path`, don't write to repo directory
- [ ] **Reference outputs are updated** — if expected outputs change, `testData/outputs/` is updated
- [ ] **All existing tests still pass** — no regressions

### 5. Documentation

- [ ] **`doc/` is updated** for user-facing behavior changes
- [ ] **`.fdtd.json` changes reference `doc/fdtdjson.md`** — verify input format is documented
- [ ] **New CMake options are documented** in `AGENTS.md` and `CLAUDE.md`
- [ ] **Tutorial examples are updated** if the workflow changes

### 6. AI-Assisted Contributions

If the PR has significant AI-generated content:
- [ ] PR description clearly states AI assistance was used
- [ ] Extent of AI assistance is described
- [ ] **At least two human reviewers** have approved (per CONTRIBUTING.md)
- [ ] Human author has taken responsibility for correctness
- [ ] Code quality is indistinguishable from human-written code in this repo

## Common Issues to Watch For

### Conditional Compilation Gotchas

    ! BAD: Missing guard for optional feature
    call some_mtl_function()   ! Only available with SEMBA_FDTD_ENABLE_MTLN

    ! GOOD: Properly guarded
    #ifdef SEMBA_FDTD_ENABLE_MTLN
    call some_mtl_function()
    #endif

### MPI Safety

- [ ] Non-blocking operations (`MPI_ISEND`/`MPI_IRECV`) are followed by `MPI_WAITALL`
- [ ] Buffer sizes account for ghost cell exchange correctly
- [ ] Wire MPI communication is updated if wire models are modified
- [ ] `MPI_Barrier` usage is appropriate (not causing deadlocks)

### OpenMP Safety

- [ ] `COLLAPSE(n)` dimensions are valid for the loop nest
- [ ] Private/shared variable clauses are correct
- [ ] No data races on shared arrays
- [ ] Reduction clauses used where accumulators are involved

### Memory Safety

- [ ] No pointer aliasing that could cause false dependencies in parallel regions
- [ ] Derived types are properly initialized before use
- [ ] Arrays are allocated before being accessed
- [ ] `contiguous` attribute used on field pointers for cache-friendly access

### Input Format Changes

If `.fdtd.json` format is modified:
- [ ] `doc/fdtdjson.md` is updated with the new format
- [ ] Backward compatibility is maintained if possible
- [ ] Migration guide is provided for existing input files
- [ ] Parser tests (`test/smbjson/`) are updated

## Review Process

### Step-by-Step

1. **Read the PR description** — understand the intent and context
2. **Check CI status** — verify all CI jobs pass (or are expected to fail)
3. **Review the diff** — line by line for logic errors, style issues
4. **Verify test coverage** — ensure new code has tests
5. **Check conditional compilation** — verify all `#ifdef` guards are correct
6. **Verify documentation** — ensure user-facing changes are documented
7. **Run local tests** (if needed) — build and run tests for the specific change
8. **Leave constructive feedback** — be specific, suggest fixes, explain why

### Feedback Style

- Be specific about what needs changing and why
- Suggest concrete fixes when possible
- Distinguish between blocking issues (must fix) and suggestions (nice to have)
- Acknowledge good changes and improvements
- Keep discussion technical and respectful

## Quick Reference: What to Check by Change Type

### Adding a new feature
- New code follows existing conventions
- New CMake option (if optional feature)
- New tests (unit + integration)
- Documentation updated
- Test data added

### Fixing a bug
- Root cause identified and explained
- Regression test added
- Minimal changes (no over-engineering)
- All existing tests still pass

### Modifying optional feature
- `#ifdef` guards are correct for all build configurations
- Tested with feature ON and OFF
- CI matrix includes the relevant configurations

### Modifying MPI code
- Non-blocking I/O pattern preserved
- Buffer sizes correct for all domain sizes
- Wire MPI communication updated if needed
- No deadlocks or race conditions

### Modifying OpenMP code
- Parallel correctness verified
- Performance impact assessed
- Thread scaling tested
- No false sharing or cache contention

### Modifying input format
- `doc/fdtdjson.md` updated
- Parser tests updated
- Backward compatibility maintained
- Migration guide provided
