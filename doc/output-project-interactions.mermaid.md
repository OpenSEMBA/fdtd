# Output Module Project Interactions

The output module is the runtime boundary between the FDTD solver state and
published simulation artifacts.
The diagram below describes the active new-output path selected by
`CompileWithNewOutputModule`.

```mermaid
flowchart TD
  LAUNCH["launcher.F90\nsemba-fdtd"] --> APP["SEMBA_FDTD_m\nsemba_fdtd_t"]
  APP --> PARSE["Input and preprocessing\nSMBJSON / NFDE / Preprocess_m"]
  PARSE --> MODEL["Simulation model\nSGGFDTDINFO_t, media, bounds, control"]
  MODEL --> SOLVER["Solver_m\ntimestepping.F90"]

  subgraph LIFE["Output lifecycle"]
    INIT["init_outputs()\ncreate probe objects and partitions"]
    UPDATE["update_outputs()\nsample fields and currents"]
    FLUSH["flush_outputs()\nwrite buffered artifacts"]
    CLOSE["close_outputs()\nclose XDMF/frequency writers"]
    INIT --> UPDATE --> FLUSH --> CLOSE
  end

  SOLVER -->|initialise observations| INIT
  SOLVER -->|each time step| UPDATE
  SOLVER -->|scheduled flush| FLUSH
  SOLVER -->|final flush| FLUSH
  SOLVER -->|end of run| CLOSE

  FIELDS["FDTD field state\nE/H arrays and time array"] --> UPDATE
  CONTROL["sim_control_t\noutput root, rank, cadence, paths"] --> INIT
  CONTROL --> UPDATE
  CONTROL --> FLUSH
  GEOMETRY["Geometry and material state\nmedia, tags, grid steps, PML bounds"] --> INIT
  GEOMETRY --> UPDATE
  WIRES["Wire / MTLN state\nwire segments and circuit probes"] --> INIT
  WIRES --> UPDATE
  FARFIELD["farfield_m\nnear-to-far inputs"] --> INIT
  FARFIELD --> UPDATE

  INIT --> COORD["output_m\nmanifest, lifecycle, ownership"]
  UPDATE --> PROBES["Probe implementations\npoint / wire / bulk / movie / frequency / far field"]
  FLUSH --> PROBES
  COORD --> PROBES

  PROBES --> BINARY["outputBinary_m\nportable binary records"]
  PROBES --> META["outputMetadata_m\nprobe JSON descriptors"]
  PROBES --> VIS["outputVisualisation_m\nXDMF + HDF5"]
  PROBES --> VTK["mapVTKOutput_m + vtkAPI_m\nVTK geometry"]
  PROBES --> TEXT["Human-readable probe files\n.dat and related files"]

  DECOMP["outputDecomposition_m\nrank partitions and offsets"] --> COORD
  COLLECTIVE["outputCollective_m\nowner and publication mode"] --> COORD
  TRANSPORT["outputTransport_m\nMPI eligibility and gathers"] --> COORD
  MPI["MPI runtime\nrank coordination"] --> TRANSPORT
  MPI --> COLLECTIVE
  COORD --> DECOMP
  COORD --> COLLECTIVE

  UTILS["fdtd-utils\ndirectoryUtils, allocation, paths"] --> BINARY
  UTILS --> META
  UTILS --> VIS
  UTILS --> TEXT
  XDMF["XDMF::HDF5 external library"] --> VIS

  BINARY --> ARTIFACTS["Run artifacts\n.bin, .dat, .xdmf, .h5, .vtk/.vtu"]
  MANIFEST["Root-owned run manifest\n*_output_manifest.json"]
  COORD --> MANIFEST
  ARTIFACTS --> CONSUMERS["Post-processing and visualisation\nParaview / scripts / users"]
  MANIFEST --> CONSUMERS

  CLEANUP["semba_end()\noptional intermediate cleanup"] --> DELETE["delete_run_output_manifest()"]
  APP --> CLEANUP
  DELETE --> MANIFEST

  subgraph TESTS["Verification boundary"]
    UNIT["test/unit/output\nunit and contract tests"]
    MPITEST["test/mpi/output\ncollective and root aggregation"]
  end
  UNIT -. "links and tests" .-> OUTLIB["fdtd-output / vtkAPI"]
  MPITEST -. "links and tests" .-> OUTLIB

  LEGACY["Legacy observation dispatch\nFlushObservationFiles / UpdateObservation"] -. "preprocessor fallback" .-> SOLVER

  classDef runtime fill:#dbeafe,stroke:#2563eb,stroke-width:2px,color:#111827;
  classDef output fill:#dcfce7,stroke:#16a34a,color:#111827;
  classDef artifact fill:#fef3c7,stroke:#d97706,color:#111827;
  classDef external fill:#f3e8ff,stroke:#9333ea,color:#111827;
  classDef legacy fill:#fee2e2,stroke:#dc2626,stroke-dasharray:5 5,color:#111827;
  class LAUNCH,APP,PARSE,MODEL,SOLVER,FIELDS,CONTROL,GEOMETRY,WIRES,FARFIELD,CLEANUP runtime;
  class INIT,UPDATE,FLUSH,CLOSE,COORD,PROBES,DECOMP,COLLECTIVE,TRANSPORT,UNIT,MPITEST,OUTLIB output;
  class BINARY,META,VIS,VTK,TEXT,ARTIFACTS,MANIFEST,CONSUMERS,DELETE artifact;
  class MPI,UTILS,XDMF external;
  class LEGACY legacy;
```

## Lifecycle Contract

Initialisation creates one output object for each supported observation request
and prepares volumetric partitions when needed.
Updates consume the current solver fields and advance probe buffers.
Flushes persist scalar, wire, bulk, far-field, binary, and visualisation data.
The final flush requests far-field completion and then closes long-lived
writers.

The coordinator also owns output lifecycle metadata, canonical scalar-writer
selection, volumetric publication mode, and the root-owned run manifest.
The manifest is removed by `semba_end()` when intermediate-data cleanup is
requested.

The legacy branch remains in `timestepping.F90` behind the inverse preprocessor
condition and is explicitly scheduled for retirement in task `T19` of the
robust-output-layer plan.
