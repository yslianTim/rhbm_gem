#pragma once
#include "support/UniqueStencilGrid.hpp"

namespace second_stage_test::matched::atom_union {
inline constexpr double radius=2.5;
unique_grid::Grid BuildGrid(const std::vector<Atom> &, const std::vector<Stencil> &,
    const rhbm_gem::MapObject & generation, const rhbm_gem::MapObject & observation, double cutoff=radius);
Eigen::SparseMatrix<double> BuildDesign(const unique_grid::Grid &, const std::vector<Atom> &,
    const rhbm_gem::MapObject & generation, double cutoff=radius);
joint_ac::Blocks BuildBlocks(const unique_grid::Grid &, const std::vector<Atom> &,
    const rhbm_gem::MapObject & generation, const std::vector<double> & alphas, double cutoff=radius);
} // namespace second_stage_test::matched::atom_union
