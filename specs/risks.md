# Ambiguity and Risk Register

| Item | Type | Risk level | Evidence |
| --- | --- | --- | --- |
| Frequency-slice binary writes x as the third vector value | Compatibility quirk | Critical | VERIFIED by binary record behaviour |
| Frequency field accumulation assigns beyond one addressed value | Compatibility quirk | Critical | VERIFIED by accumulation behaviour |
| Far-field external accumulator contract is only partially specified here | Missing external contract | High | VERIFIED wrapper inputs; external writer not fully extracted |
| Wire flavour variants depend on compile-time options | Configuration ambiguity | High | VERIFIED conditional behaviour |
| MPI coordinate rotation depends on compile-time option | Configuration ambiguity | High | VERIFIED conditional behaviour |
| Register filename formatting may include whitespace before adjustment | Compatibility ambiguity | Medium | VERIFIED path construction pattern |
| Some fatal warnings come from shared reporting behaviour | Error contract ambiguity | Medium | VERIFIED message prefixes only |
| Concurrency and simultaneous flush behaviour are not defined | Missing behaviour | Medium | no locking behaviour observed |
| File creation failure handling is inconsistent | Failure ambiguity | Medium | some paths stop, some only print |
| Output buffer overflow beyond fixed buffer length is not specified | Missing edge case | High | fixed buffers observed, no overflow scenario found |
| Far-field spec says time domain while frequency fields are passed | Domain ambiguity | High | VERIFIED current wrapper behaviour |
| Existing scalar, wire, block, and map outputs may not yet use per-probe folders | Implementation drift | High | VERIFIED by previous extraction and refined product requirement |
| Geometry map metadata is written as one legend string | Compatibility quirk | Low | VERIFIED output side effect |
| Direct output tests are gated by a compile-time option | Coverage risk | Medium | VERIFIED test gating |

## Explicit Unknowns

| Unknown | Why it matters |
| --- | --- |
| Exact external far-field file formats | Required for a complete far-field rewrite |
| Exact shared fatal-report side effects | Required to reproduce process exit and logging |
| Exact HDF5 failure behaviour for every failed call | Required for failure-compatible exports |
| Exact behaviour when pending sample count exceeds buffer capacity | Required for long-running simulations |
| Exact concurrent writer behaviour | Required if parallel output writes share paths |

## Non-Applicable Risk Areas

| Area | Reason |
| --- | --- |
| User authentication risk | no user actor participates in this subsystem |
| Tenant or ownership risk | no tenant or owner concept participates here |
| Database migration risk | no database schema participates here |
| Notification risk | no email, SMS, push, or in-app notification exists here |
| Background queue risk | no queued background work exists here |
