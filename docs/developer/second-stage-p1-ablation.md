# P1 retained final-polish changes

The retained changes build on `6e7ac3c5f452812da8aeb1b1ae2497b5cd0143ec`:

1. Final dependency polish runs only for `converged` stops when enabled.
   Non-converged stops persist the existing selected base state directly.
   Polished candidates still require strict operator recertification.
2. Final polish uses a deterministic radius `1.0` for every member and round,
   referenced to that round's endpoint rather than outer radius history.
   Ordinary boundary correction retains its outer-previous trust reference.

The experimental removal of cooperative rescue and production Grow was reverted
at the user's request. Both mechanisms and their existing test assertions are
restored to the baseline. The rescue-removal experiment had failed the existing
boundary intensity-scaling regression; its historical evidence remains in
`build/p1-ablation` and does not describe the current implementation.

After restoration, the full existing default CTest suite passes 16/16 groups,
including the previously failing intensity-scaling case. No test cases were
added and no numerical assertions or tolerances were weakened. The only retained
test edit adapts the existing final-polish case to the internal interface change.
Current validation logs are under `build/p1-restored`.

The [execution map](second-stage-p0-structure.md#ownership-and-execution) shows
rescue and transaction publication followed by converged-only final polish.
Public options, CLI and model formats are unchanged. Existing-test acceptance
does not establish complete branch coverage or real-data quality/performance
improvement; no 6Z6U, fold-168 or corpus run was performed.
