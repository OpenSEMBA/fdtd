---
name: fortran-refactor-cleanliness
description: Use whenever cleaning up, simplifying, renaming, reorganizing, or behaviour-preservingly refactoring Fortran FDTD code in semba-fdtd. Trigger for .F90 modules, solver loops, parser code, output code, CMake Fortran target changes, tests, dead-code cleanup, import cleanup, naming cleanup, and readability-only edits. Use this even when the user says "just tidy this" or "make it easier to read"; combine with fortran-performance-awareness for hot loops and with fortran-module-readability for module API or boundary changes.
---

# Fortran Refactor Cleanliness

Use this skill when refactoring this repository's Fortran code, especially when the intended outcome is easier reading or maintenance rather than new behaviour.
The goal is cleaner, more readable code without changing numerical behaviour unless the user explicitly asks for a behavioural fix.
For performance-sensitive refactors in hot loops, MPI/OpenMP paths, allocation-heavy code, or numerical kernels, also use `fortran-performance-awareness`.

## First Pass

Before editing, classify the refactor so the diff stays reviewable:

1. Name the single readability problem being fixed, such as duplicated local logic, unclear branching, misleading names, broad imports, or tangled setup steps.
2. Check whether the code is a hot numerical path, public module API, parser/output compatibility path, or optional-feature branch.
3. Search for callers before renaming procedures, types, components, public exports, files, or CMake source entries.
4. Decide the smallest behaviour-preserving change that solves the readability problem.
5. Keep feature work, bug fixes, formatting churn, and broad style normalization out of the same patch unless the user explicitly asked for them.

## Project Context

This project is `semba-fdtd`, a Finite-Difference Time-Domain electromagnetic solver written primarily in Fortran.
It contains legacy numerical kernels, newer typed modules, conditional compilation, MPI/OpenMP paths, optional HDF5 output, SMBJSON parsing, and MTLN/SPICE coupling.

Key areas:

- `src_main_pub/`: core solver, preprocessing, postprocessing, time stepping, geometry, main types.
- `src_json_parser/`: `.fdtd.json` parser and conversion helpers.
- `src_output/`: probe, VTK, XDMF, HDF5, and output utility code.
- `src_mtln/`: multiconductor transmission-line solver and ngspice integration.
- `src_conformal/`: conformal mapping and staircase reduction.
- `src_wires_pub/`: wire and thin-wire models.
- `test/`: Fortran/C++ unit tests and Python integration tests.

The build is layered through CMake static libraries.
Respect the existing dependency direction and avoid introducing upward dependencies between lower-level libraries and higher-level solver/output code.

## Refactoring Priorities

Prefer small, behaviour-preserving changes:

- Improve names when the domain meaning is clear from nearby code.
- Reduce duplicated local logic when a helper makes the numerical intent easier to read.
- Clarify long conditionals, `select case` blocks, and repeated string/coordinate formatting.
- Add `only:` to `use` statements when touching imports and when it does not create excessive churn.
- Tighten visibility with `private` and explicit `public` lists when working in already-structured modules.
- Add or improve `intent` declarations where missing and obvious.
- Replace magic literals only when their meaning is certain and the new name is local or already established.
- Remove dead local variables only after checking preprocessor branches and nearby compile options.
- Improve CMake source organization only when it directly follows from a file move, new module, or target-boundary cleanup.

Avoid broad rewrites:

- Do not redesign modules, data ownership, or library boundaries as part of a readability task.
- Do not change physics formulas, update ordering, time-step sequencing, boundary semantics, MPI exchange order, or output formats unless explicitly requested.
- Do not rename domain terms with historical or input-file significance without asking.
- Do not split every long routine mechanically; extract only natural concepts with clear names.
- Do not introduce compatibility shims or abstraction layers unless there is a concrete caller or persisted-data need.
- Do not normalize formatting across an entire file when the functional refactor touches only a small region; formatting-only churn hides semantic review.

## Fortran-Specific Guardrails

Treat numerical and memory semantics as part of behaviour:

- Preserve `implicit none`.
- Preserve `kind` choices such as `RKIND`, `RKIND_tiempo`, `SINGLE`, and existing integer kinds unless the task is specifically about precision or portability.
- Preserve array rank, shape, lower bounds, allocation ownership, pointer association, and `contiguous` assumptions.
- Be careful when changing `pointer` to `allocatable` or vice versa; this can alter aliasing and ownership.
- Do not reorder loops over field arrays unless there is a measured performance or correctness reason.
- Preserve OpenMP and MPI assumptions around shared data, halo exchanges, reductions, and output ordering.
- Keep preprocessor branches such as `CompileWithMPI`, `CompileWithMTLN`, `CompileWithSMBJSON`, and `CompileWithReal8` valid even if the local build uses only one configuration.
- Preserve file formats, exact labels, and serialized names used by JSON, VTK, XDMF, HDF5, probe `.dat` files, or legacy `.fdtd` inputs.
- New API and implementation must target strict Fortran 2018 or newer; preserve its CMake standard requirement and disabled compiler extensions.

## Workflow

Before editing:

1. Read the target module and enough neighbouring code to understand ownership and callers.
2. Search for procedure/type/module references before changing names, signatures, public exports, or file-level interfaces.
3. Identify relevant compile flags around the code, especially `#ifdef` blocks.
4. Choose one coherent refactor with a small reviewable diff.

During editing:

1. Keep the public interface stable unless changing it is the point of the task.
2. Prefer local changes over new global helpers.
3. Preserve surrounding formatting style unless the formatting itself harms readability.
4. Add comments only to explain non-obvious domain constraints, numerical assumptions, or compiler/preprocessor constraints.
5. Do not combine cleanup with unrelated feature work.
6. Preserve exact strings and ordering in user-visible output unless the task is explicitly about changing output.

After editing:

1. Build or at least compile the affected target when feasible.
2. Run the most targeted tests available.
3. State that no behavioural change is intended, or describe the exact intended behaviour change if there is one.
4. Mention any build/test coverage gaps caused by unavailable dependencies or disabled options.
5. Re-read the diff from a reviewer perspective and remove incidental churn that does not support the stated refactor.

## Readability Checklist

Use this checklist before finishing a refactor:

- The code reads in the same order as the operation it performs.
- Names reflect domain meaning, not just type or storage.
- Conditionals have clear cases and meaningful default/error handling.
- Public module surface is no larger than necessary.
- Imports are understandable and preferably constrained with `only:` where practical.
- Local variables are declared near the routine where they are used and are not misleadingly reused.
- Comments explain why the code exists or why a surprising choice is necessary.
- The diff is small enough to review for numerical equivalence.
- The final response names the intended behaviour-preserving nature of the change.

## Testing And Verification

Useful build commands:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
```

Unit tests:

```bash
./build/bin/fdtd_tests
```

Python integration tests:

```bash
pytest test/ --durations=20
```

Marker-specific tests:

```bash
pytest test/ -m mtln
pytest test/ -m hdf
pytest test/ -m mpi
```

Run the narrowest relevant command first.
For parser changes, prefer `test/smbjson` tests.
For output changes, prefer `test/output` and VTK/XDMF tests.
For MTLN changes, prefer `test/mtln` tests.
For solver loop changes, build and run the broad unit suite if feasible.

## Good Refactor Examples

Good changes:

- Replace repeated coordinate string assembly in one module with a local helper that preserves exact output text.
- Extract a named local predicate from a long conditional when it is used multiple times in the same routine.
- Convert an unclear `if/elseif` chain to `select case` without changing defaults.
- Add `intent(in)`, `intent(out)`, or `intent(inout)` to arguments where usage is unambiguous.
- Split a long parser routine into parse, validate, and convert steps when those steps already exist conceptually.

Risky changes that need explicit justification:

- Changing loop order over `Ex`, `Ey`, `Ez`, `Hx`, `Hy`, or `Hz` arrays.
- Replacing pointer arrays with allocatables.
- Changing `real` or `integer` kinds.
- Renaming JSON labels, probe labels, field prefixes, or legacy NFDE terms.
- Removing preprocessor branches that are not active in the current build.
- Introducing generic abstractions across solver, parser, output, and MTLN code without a concrete repeated pattern.

## Final Response Expectations

When reporting a completed refactor, include:

- The files changed.
- The readability improvement made.
- Whether behaviour was intended to change.
- The build or tests run, including failures or skipped verification.

Keep the response concise and factual.
