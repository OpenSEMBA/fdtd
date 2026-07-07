# Rewrite Boundary

| Concern | MUST remain identical | MAY change internally | Notes |
| --- | --- | --- | --- |
| Output registration | request-to-output mapping, skipped cases, fatal cases | allocation strategy | preserve observation-visible outcomes |
| Domain normalisation | time and frequency clamping rules | arithmetic organisation | preserve exact edge outcomes |
| Path generation | per-probe folder, prefixes, suffixes, separators, MPI rotations | helper organisation | paths are compatibility contracts |
| Time data files | folder placement, append semantics, and row shapes | writer buffering | preserve existing row order |
| Frequency point files | folder placement, replace semantics, and row shapes | writer buffering | repeated flush behaviour matters |
| Wire current output | six-column row contract and formulas | wire data access design | preserve sign and voltage calculations |
| Wire charge output | two-column row contract | wire data access design | preserve unsupported-flavour rejection |
| Block output | integral formulas and row contract | loop organisation | preserve boundary indexing |
| Volumetric selection | in-bounds and material validity rules | storage layout | selected point count affects outputs |
| Movie output | binary record shape and HDF5/XDMF names | HDF5 library wrapper | preserve misspelled dataset names |
| Frequency-slice output | accumulation, binary replacement, close-time HDF5 writes | memory layout | preserve compatibility quirks unless migrated |
| Far-field gateway | input values passed to external far-field behaviour | internal orchestration | external file contract needs more extraction |
| Geometry map | selected topology, VTU cells, metadata legend | grid writer implementation | preserve tags and cell types |
| Structured export | VTK, HDF5, and XDMF observable files | implementation language | preserve file readers' expectations |
| Error behaviour | fatal messages and no-data messages | reporting mechanism | exact process-exit semantics need validation |
| Database contract | no database tables or writes | not applicable | do not introduce database dependency |

## Modernisation Boundary

| May change | Condition |
| --- | --- |
| Internal type names | only if all output files and behaviours stay identical |
| In-memory buffers | only if flush ordering and overflow behaviour are specified |
| HDF5 writer abstraction | only if dataset names, ranks, and append semantics stay identical |
| VTK writer abstraction | only if generated topology and data arrays stay compatible |
| Reporting abstraction | only if fatal and warning behaviours stay observable-equivalent |

## Must Not Change Without Migration

| Contract | Reason |
| --- | --- |
| `CurrenDensity` dataset spelling | external readers may depend on it |
| frequency binary third-component repetition | existing consumers may parse it |
| `_tm.dat` and `_fq.dat` suffixes | wrapper tools discover probe files by suffix |
| per-probe output folders | generated artifacts must not be mixed beside the case root |
| coordinate suffix order | output file identity depends on it |
| output register final `END!` line | downstream tooling may use it as sentinel |
