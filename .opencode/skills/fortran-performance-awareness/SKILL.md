---
name: fortran-performance-awareness
description: Use whenever writing, reviewing, diagnosing, or refactoring performance-sensitive Fortran in semba-fdtd. Trigger for FDTD time-step loops, field/material/boundary/wire/MTLN updates, array access patterns, allocation or temporary-array concerns, MPI/OpenMP communication, I/O cadence, profiling results, or any user request about speed, scaling, memory use, or "making this faster." Use this even for small-looking cleanups inside numerical kernels because preserving update order and array semantics matters more than cosmetic refactoring.
---

# Fortran Performance Awareness

Use this skill when working on performance-sensitive Fortran code in this repository, including review and diagnosis tasks where no code may be edited.
The goal is efficient code that preserves numerical correctness, keeps the physics readable, and avoids premature or unmeasured rewrites.

## First Pass

Start by locating the performance risk before proposing changes:

1. Determine whether the code is in a repeated time-step path, setup path, output path, parser path, or test-only path.
2. Identify loop bounds, array ranks, lower bounds, halo regions, and any OpenMP or MPI ownership assumptions.
3. Look for obvious repeated costs such as allocation, I/O, string formatting, polymorphic dispatch, map lookups, or large temporary arrays.
4. Separate proven bottlenecks from plausible risks. A structural cleanup can be useful, but do not report it as a measured speedup without measurement.
5. If the task is mostly module/API restructuring, also use `fortran-module-readability`; if it is a broad cleanup, also use `fortran-refactor-cleanliness`.

## Core Principles

Optimize deliberately:

- Preserve numerical behaviour unless the task explicitly asks for a behaviour change.
- Prefer simple, predictable hot loops over clever abstractions.
- Keep performance-sensitive code readable enough to audit for physics and indexing correctness.
- Optimize the actual bottleneck when measurements or code structure make it clear.
- Avoid broad rewrites when a local change removes the cost.
- Treat MPI/OpenMP synchronization, allocation, and I/O as explicit performance costs.

## Hot Path Priorities

Pay special attention to code inside time-step loops, field-update kernels, material updates, boundary-condition application, MPI exchange paths, and output sampling.

Prefer code that:

- Avoids allocation and deallocation inside repeated update loops.
- Avoids repeated string operations, file operations, type conversions, and lookups in numerical kernels.
- Keeps loop bodies small and branch structure understandable.
- Reuses precomputed constants when they are truly invariant over the loop.
- Avoids unnecessary temporary arrays, especially large field-sized temporaries.
- Keeps I/O cadence intentional and outside kernels when possible.
- Keeps synchronization points minimal and tied to actual data dependencies.

When code is not on a hot path, avoid adding complexity for hypothetical speed. A parser or one-time setup routine can often favor clarity unless it allocates field-sized data repeatedly or dominates large-case startup.

## FDTD-Specific Guardrails

Treat update order as part of correctness:

- Do not reorder electric and magnetic field updates without understanding the time-stepping scheme.
- Do not reorder boundary-condition, material, wire, MTLN, source, or observation operations casually.
- Preserve Yee-grid staggering, coordinate indexing, lower bounds, and halo assumptions.
- Preserve dispersive and anisotropic material update sequencing.
- Preserve CPML, Mur, MPI halo exchange, and far-field sampling ordering.
- Treat output labels, sampling cadence, and serialized array layout as observable behaviour.

When an optimization touches these areas, verify with tests or a representative case whenever feasible.

## Fortran Efficiency Guidelines

Avoid accidental costs common in Fortran:

- Prefer contiguous memory access patterns in hot loops.
- Be careful with array slices passed to procedures; non-contiguous slices can create temporaries.
- Avoid whole-array expressions on large arrays when they obscure allocation or temporary creation in hot paths.
- Use explicit interfaces and argument `intent` so compilers can reason about calls.
- Avoid unnecessary `pointer` aliasing in numerical kernels.
- Prefer `allocatable` for owned storage in new code, but do not change existing pointer ownership semantics just for style.
- Preserve established kind choices such as `RKIND`, `RKIND_tiempo`, `SINGLE`, and integer kinds unless precision or portability is the task.
- Avoid converting scalar helper functions into calls inside deeply nested loops unless the compiler can inline them or the readability benefit is worth the cost.
- Keep frequently used scalar values local when that avoids repeated component dereferences in hot loops.
- Be cautious with assumed-shape dummy arguments and array expressions in helper procedures called from kernels; check whether they can introduce copying or inhibit compiler optimization.
- Keep data layout and loop nesting aligned with Fortran column-major storage when this does not conflict with stencil dependencies or established code style.

Do not make code less obviously correct for a theoretical speedup. If an optimization relies on a non-obvious compiler or memory-layout assumption, document it briefly.

## Loop And Array Changes

Before changing loops over field arrays:

1. Identify the array dimensions, lower bounds, and memory layout.
2. Check whether loop order is chosen for cache locality, stencil dependencies, MPI halos, or readability.
3. Confirm whether OpenMP directives, reductions, or private variables depend on the current structure.
4. Preserve boundary ranges and off-by-one behaviour exactly unless fixing a known bug.
5. Compare results after the change when feasible.

Avoid combining loop-order changes with unrelated cleanup. They should be reviewable as performance-sensitive changes.

When changing a loop for locality, vectorization, or OpenMP scheduling, state the intended performance property in the code review summary. This helps reviewers distinguish deliberate numerical-kernel work from incidental formatting.

## Parallel Code

For MPI and OpenMP paths:

- Keep data-sharing attributes explicit and correct.
- Preserve reductions and their numerical implications.
- Avoid moving MPI communication across computation phases unless dependencies are fully understood.
- Avoid adding barriers, critical sections, or atomics unless they are required for correctness.
- Watch for false sharing when introducing per-thread scratch data.
- Keep thread-local scratch allocation outside repeated parallel regions when practical.
- Preserve deterministic output ordering where the code currently guarantees it.

Parallel performance changes must not weaken correctness under configurations that are not active in the local build.

## Allocation And Ownership

Allocation changes are performance and correctness changes:

- Allocate once at setup time when the size is known and reused across time steps.
- Deallocate at clear ownership boundaries.
- Avoid hidden reallocation from assignment to allocatable arrays in hot paths.
- Keep scratch arrays local only when their size is small or the routine is not hot.
- Avoid module-level scratch state unless there is a clear ownership and thread-safety story.
- Preserve pointer association and aliasing semantics when editing legacy structures.

## Readability Balance

Efficient code should still be maintainable:

- Keep domain names visible in formulas and update steps.
- Prefer a clear local scalar or named coefficient over repeated dense expressions.
- Avoid abstractions that hide stencil shape, update ordering, or boundary handling.
- Add short comments for non-obvious performance constraints such as loop order, contiguous assumptions, or synchronization placement.
- Do not add comments that simply say code is faster.

## Measurement And Verification

Use the narrowest practical verification first:

```bash
cmake --build build -j
```

Then run targeted tests for the touched area when available. For performance-focused changes, also compare a representative case runtime when practical, ideally with the same input, build type, MPI rank count, thread count, and output cadence.

When reporting results, distinguish between:

- Measured speedup or reduced runtime.
- Structural improvement likely to reduce overhead, such as removing allocation from a loop.
- Readability-neutral cleanup that only prepares for future optimization.

Do not claim a performance improvement as measured unless it was actually measured.

If measurement is not practical, report the change as a structural improvement, for example "moves allocation out of the time-step loop" or "preserves loop order while reducing repeated component lookups."

## Review Checklist

Before finishing performance-sensitive Fortran work, check:

- Numerical behaviour is intended to remain unchanged, or the intended change is explicit.
- No unnecessary allocation, deallocation, I/O, string work, or large temporaries were added to hot paths.
- Array access and loop order are deliberate and preserve existing dependencies.
- MPI/OpenMP synchronization and data-sharing remain correct.
- Optional preprocessor branches remain valid.
- Public APIs were not expanded for performance helpers that only have one caller.
- Any non-obvious performance constraint is documented briefly.
- Verification or measurement is reported accurately.

## Final Response Expectations

When reporting performance-related changes, include:

- The files and hot paths touched.
- The specific overhead avoided or performance property preserved.
- Whether behaviour was intended to change.
- The build, tests, or measurements run.
- Any performance claims that are unmeasured and should be treated as expectations rather than results.
