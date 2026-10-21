# Joint component construction contract

Structural components are formed from observation support, independently of
current amplitudes, charges, residual values, or numerical sensitivity. Every
row connects all contributors recorded for that observation. The partition
preserves stable atom identity, original and fitting row identities, and
parent-to-component mappings.

Rows with no fitting contributor remain constant objective rows and are
included exactly once in the assembled loss. Unsupported atoms remain explicit
unobserved entries; they do not receive invented parameters or profile claims.
Components with failed numerical states remain isolated and do not fill missing
predictions with zero.

Compact component views map into the parent problem without copying its data.
Malformed or duplicate identities are reported explicitly. Parent observation
scale and coefficient directions are preserved across component solves.

The deterministic coverage in `Components_test.cpp` checks shared-row
connectivity, bridge and constant rows, stable IDs across storage permutations,
compact row mappings, parent scale, failure isolation, and unsupported atoms.
The current production and numerical policies are described in the
[Joint runtime guide](joint-component-runtime.md).
