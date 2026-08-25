---
name: fortran-module-readability
description: Use whenever creating, splitting, moving, or substantially restructuring Fortran modules in semba-fdtd. This is the right skill for .F90 module boundaries, public/private APIs, use/import lists, derived-type ownership, procedure ordering, CMake source placement, and naming. Prefer this skill for module-level design questions even if the user only says "clean up this module" or "add a new type"; combine it with fortran-performance-awareness when the module contains hot FDTD kernels, MPI/OpenMP paths, or large array operations.
---

# Fortran Module Readability

Use this skill when writing a new Fortran module, moving code across modules, adding exported types or procedures, or substantially restructuring an existing module in this repository.
The goal is a module that has one clear purpose, a small public surface, readable internal organization, and dependencies that fit the existing `semba-fdtd` library layering.
For modules containing hot kernels, large array operations, MPI/OpenMP code, or repeated allocation paths, also use `fortran-performance-awareness`.

## First Pass

Before changing a module boundary, build enough context to avoid accidental API or dependency damage:

1. Identify the module's current responsibility, public exports, and callers.
2. Check nearby modules in the same source directory for established naming, import, and visibility style.
3. Search for public type, component, procedure, and module-name references before renaming or moving anything.
4. Note optional preprocessor branches that may expose different dependencies from the local build.
5. Decide whether the task is module design, ordinary refactoring, performance work, or a mix; use the narrowest applicable change.

## Core Principles

Write modules as explicit boundaries:

- Give each module one primary responsibility that can be stated in one sentence.
- Keep public APIs small, intentional, and stable.
- Prefer clear domain names over generic implementation names.
- Separate data definitions, parsing/conversion, numerical operations, and output concerns unless they are naturally part of the same abstraction.
- Avoid introducing new global state; pass data through arguments or owned derived-type components when practical.
- Keep changes local unless a broader module boundary change is explicitly requested.
- New API and implementations must target strict Fortran 2018 or newer; keep its CMake targets standard-required with extensions disabled.

## Module Skeleton

Prefer this structure for new or cleaned-up modules:

```fortran
module example_module
  use dependency_module, only: dependency_type, dependency_routine

  implicit none

  private


  public :: example_type
  public :: init_example


  type :: example_type
    ! Components here.
  end type example_type

contains

  subroutine init_example(example)
    type(example_type), intent(out) :: example


  end subroutine init_example

end module example_module
```

Use existing local formatting when editing established files, but keep the same conceptual order when it does not create churn:

1. `module` statement.
2. `use` statements.
3. `implicit none`.
4. default visibility, usually `private`.
5. explicit `public` declarations.
6. parameters, interfaces, and derived types.
7. `contains` procedures.
8. `end module module_name`.

## Public API Design

Default to a private module and explicit exports:

- Use `private` near the top of the module unless the existing module style makes that too disruptive.
- Export only the types, constants, interfaces, and procedures intended for callers.
- Keep helper procedures private by default.
- Avoid exporting implementation details just because tests or nearby code can reach them; prefer testing through the real public behaviour.
- Avoid renaming public procedures or derived-type components that are already used broadly unless the task explicitly includes API cleanup.
- Preserve names that map to input formats, output labels, physics terms, or legacy NFDE/JSON concepts unless the user approves a compatibility-impacting change.

Good public APIs read like a short capability list. If the public list is long, look for mixed responsibilities before adding more exports.

When changing an existing public API, distinguish between three cases:

- Internal-only exports with few callers can often be cleaned up in one small change after checking references.
- Input/output-facing names, serialized labels, and physics terms are compatibility-sensitive; preserve them unless the user asked for that change.
- Widely used interfaces should be changed only when the replacement makes callers clearer enough to justify the churn.

## Dependency Hygiene

Keep imports narrow and dependency direction clear:

- Prefer `use module_name, only: symbol_a, symbol_b` for new modules and touched imports.
- Do not introduce circular dependencies.
- Respect the existing CMake library layering: lower-level type/report/parser/component code should not depend on higher-level solver, launcher, or output orchestration code.
- Keep optional-feature boundaries intact for `CompileWithMPI`, `CompileWithMTLN`, `CompileWithSMBJSON`, `CompileWithReal8`, and similar preprocessor paths.
- Avoid a shared utility module unless there is a concrete repeated pattern across multiple callers.
- Prefer passing dependencies as arguments over reaching into unrelated modules for mutable state.

When a new module needs symbols from several distant layers, pause and check whether the responsibility belongs somewhere else.

## Derived Types

Use derived types to express ownership and domain concepts clearly:

- Name types after the domain object they represent, not just storage shape.
- Keep components cohesive; avoid derived types that become bags of unrelated solver, parser, and output state.
- Make allocation ownership obvious from initialization and finalization routines.
- Prefer type-bound procedures only when they clarify ownership or behaviour; do not add them mechanically.
- Keep pointer components only when aliasing or association is required. Prefer existing project patterns over changing ownership semantics.
- Document non-obvious invariants such as array bounds, coordinate ordering, field staggering, or MPI-local versus global indexing.

If a derived type needs many setup steps, provide one clear initializer instead of requiring callers to know internal ordering.

## Procedure Organization

Write procedures that expose intent without hiding numerical behaviour:

- Put the public, high-level procedures first when that matches the reading flow.
- Keep private helpers close to the public procedure that uses them when they are local to that concept.
- Use explicit `intent(in)`, `intent(out)`, or `intent(inout)` for arguments.
- Avoid long argument lists when a cohesive derived type already exists, but do not introduce a derived type just to reduce line length.
- Extract helpers for named concepts, not for every small block.
- Preserve loop order, update ordering, and array shape semantics in numerical code.
- Keep error handling and validation near the boundary where invalid data enters the module.

A good procedure name should let a caller understand what happens without reading the body. A good body should still make the numerical or data-flow steps visible.

## Naming And Comments

Prefer names that match the repository's domain language:

- Use established field, material, geometry, parser, output, MTLN, and boundary-condition terminology.
- Avoid overly broad names such as `manager`, `handler`, `processor`, or `data` unless that is already the local convention.
- Use consistent prefixes only when they help group related procedures or avoid ambiguity.
- Name boolean procedures and variables as predicates when practical, for example `is_valid`, `has_source`, or `uses_mpi`.
- Add comments for surprising constraints, physical assumptions, indexing conventions, preprocessor requirements, or file-format compatibility.
- Do not add comments that restate the code line-by-line.

## File And Build Integration

When adding a module:

1. Place the file in the directory matching its layer and responsibility.
2. Add it to the correct CMake target or source list without changing unrelated target boundaries.
3. Confirm optional-feature guards still compile for relevant configurations.
4. Search for similarly named modules before creating a new one.
5. Prefer a small module over modifying a broad catch-all module, but do not fragment one concept across many files.

When moving code between existing modules, update CMake only if the file list or target ownership changes. Avoid moving a source file into a higher-level library just to access a dependency; that often signals the responsibility belongs elsewhere.

## Modularity Guardrails

Avoid these patterns unless there is an explicit reason:

- A module that mixes parsing, solver state updates, and output formatting.
- A new dependency from `src_conformal`, `src_json_parser`, `src_mtln`, or lower-level type code into `src_main_pub` orchestration code.
- A public helper module created for one caller.
- A module-level mutable variable used as hidden communication between procedures.
- A generic abstraction that erases important electromagnetic, indexing, or file-format meaning.
- A rewrite from procedural module style to object-oriented style just for style consistency.

## Review Checklist

Before finishing a new or restructured module, check:

- The module purpose is clear from its name and public API.
- `implicit none` is present.
- Visibility is explicit, preferably `private` by default.
- Public exports are minimal and intentional.
- Imports use `only:` where practical.
- The module sits in the correct source directory and CMake layer.
- Derived types have clear ownership and initialization rules.
- Procedures have explicit argument intents.
- Numerical order, array bounds, precision kinds, and preprocessor branches are preserved.
- Comments explain constraints or intent, not obvious syntax.
- The module can be tested through public behaviour.

## Verification

For module-writing changes, verify with the narrowest practical command first:

```bash
cmake --build build -j
```

When the change affects parser, output, MTLN, or solver behaviour, run the corresponding focused tests when available. If a build or test cannot be run because the local environment is not configured, state that clearly.

## Final Response Expectations

When reporting module changes, include:

- The module files changed or added.
- The public API shape introduced or modified.
- The dependency or CMake integration point touched.
- Whether behaviour was intended to change.
- The build or tests run, including any skipped verification.
