# Fortran style conventions

These conventions apply to all Fortran code in this repository, including
`src_*/`, `include/`, and `test/`. They were established in
[issue #283](https://github.com/OpenSEMBA/fdtd/issues/283).

1. **English only.** All names and comments must be written in English.

2. **Lowercase keywords.** Fortran keywords are always lowercase.
   Write `type`, not `TYPE`.

3. **No space before `(` in declarations.** Write `type(edge_t)`, not
   `type (edge_t)`.

4. **Type names end in `_t`.** Prefer `type(cell_t)` over `type(cell)`.

5. **Module names end in `_m`.** Prefer `module mesh_m` over `module mesh`.

6. **Two-word endings.** Prefer `end if` over `endif`, `end do` over
   `enddo`, `end subroutine` over `endsubroutine`, and so on.

7. **Do not use Fortran keywords or intrinsic names as identifiers.**
   For example, do not use `size`, `index`, `count`, `end`, `data`, or
   `type` as variable names.

8. **No space between `(` and its interior** in function calls and
   declarations. Prefer `integer(kind=4)` over `integer ( kind=4)`, and
   `null()` over `null ( )`.

9. **`parameter` names are uppercase.** Prefer `RKIND` over `rkind`.

## Notes

- Fortran is case-insensitive, so rules 2 and 9 are pure case changes.
- These rules are enforced by review today. When a rule cannot be checked
  automatically, keep changes consistent with the surrounding code.
