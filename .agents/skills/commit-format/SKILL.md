---
name: commit-format
description: Provides this repository's required commit-message format, categories, and examples. Use when creating, reviewing, or proposing git commit messages for FDTD.
license: MIT
metadata:
  author: Elemwave
  version: "1.0"
---

When creating or reviewing commit messages, follow this workflow:

1. **Identify change layer**: Choose category for innermost layer changed. Use `Domain` for pure business logic, `ElemData` for data models and utilities, `Command` for FreeCAD commands and data containers, `Mapper` for domain-to-FreeCAD mappings, `Service` for export/validation/preferences, `Util` for FreeCAD utilities, `View` for task panels and widgets, `Test` for tests, `Refactor` for non-functional reorganisations, `CI` for automation, or `APIWrapper` for the SEMBA CLI interface. Use actual module name, such as `QtWrapper`, `Extensions`, or `Init`, if no category fits.

2. **Split changes semantically**: Each commit must be atomic: one logical, cohesive change. Cross-layer changes should be split where practical. A Domain model plus Command change normally becomes one `Domain` commit and one `Command` commit.

3. **Write repository commit message**: Use exact format `FDTD | <type> | #<Issue> | <Description>`. `FDTD` is literal. Include GitHub issue number, such as `#866`.
Allowed types
   - `feat`: new feature
   - `fix`: bug fix
   - `refactor`: code restructuring without behavior change
   - `docs`: documentation only
   - `style`: formatting, whitespace, semicolons (no logic change)
   - `test`: adding or updating tests
   - `chore`: tooling, configs, dependencies
   - `perf`: performance improvement
   - `ci`: CI/CD changes
   - `build`: build system changes
   - `revert`: reverting a previous commit
Write description in imperative mood, concise, and ideally 72 characters or fewer.
   - Good: `FDTD | chore | #866 | Use relative path for STEP file in TULIP JSON export`
   - Good: `FDTD | feat | #880 | Add validation for solver parameter ranges`
   - Bad: `fix(service): use relative path for STEP file`

4. **Apply commit rules**:
   - Use one-line messages only. Never use a commit body.
   - Never include trailers, `Co-Authored-By`, AI/tool attribution, generated-by markers, bot signatures, or AI metadata.
   - Preserve capitalisation for `FDTD`, categories, proper nouns, and identifiers.
   - Use `git add` with specific files for each commit. Do not use `git add .` unless all changes belong to same commit.
   - After each commit, run `git log -1 --format=%B` and verify message has no prohibited attribution or metadata. Then run `git status` to verify success.

5. **Execution order**: Stage and commit one semantic group at a time. Do not skip ahead.
