# FixedNeighbor endpoint retry fixture

This deterministic model/map pair reproduces the production FixedNeighbor
endpoint disagreement observed at the 107th sweep. The model contains 168
non-hydrogen atoms and produces 30,108 observations in one structural
component, partitioned into 23 blocks with a 12-atom core. `initial_widths.txt`
contains the corresponding first-stage width seeds; it is not simulation truth.

The original inputs match the reproduction record from 2026-10-12. The CIF
fixture has trailing whitespace removed; its parsed atom records are unchanged.

| File | SHA-256 |
| --- | --- |
| `model.cif` (normalized fixture) | `72fb35ea88f628e93678e6d05f880f8d09519381a6e5d7215d0836c29546ec5e` |
| `model.cif` (original input) | `156d35aa326f0d4408d726a999329d2ffede775489aeaa5d99a2cc9b9f663cab` |
| `map.map` | `aa0173a4aa6990bac5c367dcebdd8db232228e2a15b925c56a96b8327a0b815c` |

The regression checks use only solver evidence. No true parameter values are
included or used to control the search.
