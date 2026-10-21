# Joint compact SVD utilities

The compact SVD utilities remain part of the local numerical solver. They
support local active-set rank checks and bounded least-squares operations used
by Profile LM and local conditioning. They are not a global rank certification
pipeline.

`CompactSvd` dispatches to JacobiSVD for compact dimensions below 16 and Eigen
BDCSVD otherwise. Near-threshold cases and failed/nonfinite decompositions use
the existing bounded Jacobi retry. The caller's threshold policy and local
observation dimensions remain explicit. A failed retry remains unavailable;
a singular matrix is not made full rank with an arbitrary epsilon.

The local Profile LM keeps the complete profiled derivative, including the
nonzero-residual coefficient correction. The local derivative path retains its
existing tiled QR, cancellation handling, and compact least-squares solve.
Local conditioning analyzes one FixedNeighbor block at a time and reports its
conditional free A/C, projected-width, and profile derivative evidence.

Production does not request or retain global weak directions. Historical
right-vector and global projected-width studies are archived in the
[Joint evidence record](joint-component-evidence.md); their memory
estimates do not describe current runtime. The local SVD controls remain owned
by `Numerics_test.cpp`, `LocalConditioning_test.cpp`, and the Profile LM
derivative tests.
