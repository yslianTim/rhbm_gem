# Per-cluster history dependency audit and observer extraction

Baseline: `a58713d1d37fe326e55c1ccc5bc29928a45aae96`.
This is a structural cleanup, not a change to numerical acceptance policy.

## Dependency findings

| Path | Former history use | Decision dependency and retained behavior |
| --- | --- | --- |
| Initialization / partition reset | Reconcile keys and save previous parameters, reset provenance | History only; objective domain, topology, radius reconciliation and global-best refresh remain production operations. |
| Background refresh | Reset affected cluster history using the changed sample background | History only; domain revision, previous/global-best reevaluation and quarantine retry semantics remain unchanged. |
| Local search restart / acceptance | Copy baseline history, reevaluate historical parameters after the previous gate passes | Cannot select a factor or accept/reject a candidate; each restart restores observer baseline history. |
| Local polish | Update history only after strict improvement passes | Tie-break selects the retained historical record, not the polished patch. |
| Boundary / rescue | Evaluate member history against iteration-baseline references; retain provisional updates | Member gates and component/global gates are independent of history. Only the selected component observation is published. |
| Salvage / all-rejected | Restore baseline history with rolled-back keys | Model rollback and key selection do not read history. Observer receives the resulting rejection notification. |
| Commit / later iterations | Publish retained history and best-source IDs | History feeds only later history updates and diagnostics, not radius, quarantine, certificate, stop or persistence. |
| Phase audit | Copy historical references for counterfactual member comparisons | Observer-only replay; it owns a read-only baseline snapshot and isolated counters. |

The old path still performed allocations, map lookups and historical objective
reconstruction before returning production decisions. Those operations could
propagate exceptions. Removing numerical predicates alone was insufficient to
make history observation-only.

`best_audit_state` is separate and remains in production: global-best acceptance,
audit patience, non-converged final-state selection and final certification are
unchanged. Previous objective references, sample domains and all tolerances remain.

## Ownership and publication

`ClusterHistoryObserver` is created once for a non-quiet Debug second-stage run.
Info/quiet runs do not allocate its maps or evaluate historical patches. It owns
per-attempt baseline and staged maps, component observations and history-work
counters. Local workers update preallocated independent key entries; source IDs
use the existing per-key sequence and event mutex. Boundary notifications follow
worker completion.

Local evaluation returns only its production decision and numerical diagnostics.
The caller then requests a history diagnostic payload. Boundary records carry an
observation ID for diagnostic publication; numerical evaluation results carry no
history or history-dependent decision. Endpoint IDs are saved before later
correction records can grow the diagnostic vector. This preserves endpoint
history when correction fails. Component history always starts from the iteration
baseline, not local staged history. Rejected member observations still generate
the same best-comparison diagnostics; provisional source events retain their
existing retained/unretained semantics.

The production transaction publishes state, provenance, quarantine and radii
before notifying the observer to publish history. All observer entry points are
`noexcept`; failure disables further history work for the run and reports
unavailable diagnostics. No observer failure becomes a rejection, rollback or stop.
History helper functions that can throw are confined to this module and isolated
phase-audit replay; they are not called by production decisions.

## Diagnostic contract and verification

Existing best-source, best-comparison, publication and phase-audit field names and
schemas are preserved. Historical reference availability is an observer payload.
Historical sample-evaluation counts are observer-owned, so production work counts
can decrease even when numerical trajectories are identical. This is a counter
ownership change, not a measured performance claim.

Validation uses existing CTest groups with both diagnostic switches OFF/ON in all
four combinations, before and after extraction. Existing history assertions are
adapted to observer APIs; no test cases or numerical tolerances are added or
relaxed. Build-local captures of existing fixtures compare history events and
production traces. Commands, logs and the completed matrix are recorded under
`build/cluster-history/`. No new numerical corpus or fold-168 baseline is generated.

Final validation passed all 16 existing CTest groups in each of the four flag
combinations, before and after the change. The existing-fixture comparison passed
23 captures containing 25,703 production/diagnostic records; only timing fields
and unordered parallel-worker log order were excluded. No numerical comparison
was relaxed. The final-source build matrix and exact comparison rules are in
`build/cluster-history/report.md`.
