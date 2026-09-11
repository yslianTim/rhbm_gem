# P1 retained final-polish and radius changes

The retained changes build on `6e7ac3c5f452812da8aeb1b1ae2497b5cd0143ec`:

1. Final dependency polish runs only for `converged` stops when enabled.
   Non-converged stops persist the existing selected base state directly.
   Polished candidates still require strict operator recertification.
2. Final polish uses a deterministic radius `1.0` for every member and round,
   referenced to that round's endpoint rather than outer radius history.
   Ordinary boundary correction retains its outer-previous trust reference.

3. Production radius updates use Keep or Shrink only. Shrink remains `0.5`,
   with minimum `0.0625`; new topology keys start at `1.0`. Grow-key storage,
   update branches and growth cancellation are removed. Rho shadow Grow remains
   diagnostic-only and cannot update production state.

Cooperative rescue remains enabled. The initial combined rescue/Grow removal
failed the boundary intensity-scaling regression and was reverted. A subsequent
independent Grow ablation, retaining rescue and both final-polish changes, passes
the full existing default CTest suite: 16/16 groups, including intensity scaling.
Historical combined-ablation evidence remains under `build/p1-ablation`;
the current independent Grow validation is under `build/p1-grow-only`.

No test cases were added or numerical tolerances weakened. Existing production
Grow expectations now assert Keep, while Shrink, saturation, shadow Grow and
rescue assertions remain. The final-polish case retains its interface adaptation.

The [execution map](second-stage-p0-structure.md#ownership-and-execution) shows
rescue and transaction publication followed by converged-only final polish.
Public options, CLI and model formats are unchanged. Existing-test acceptance
does not establish complete branch coverage or real-data quality/performance
improvement; no 6Z6U, fold-168 or corpus run was performed.
