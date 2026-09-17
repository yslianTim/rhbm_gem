#pragma once
#include "support/AtomCenteredVoxelUnion.hpp"

namespace second_stage_test::matched::atom_block {
boost::json::object Diagnostics(const joint_ac::Evidence &, const Eigen::VectorXd & residual,
    const Eigen::VectorXd & variances, const joint_ac::Blocks &, const unique_grid::Grid &,
    const std::vector<Atom> &);
} // namespace second_stage_test::matched::atom_block
