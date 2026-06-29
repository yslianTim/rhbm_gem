#pragma once

#include <utility>
#include <vector>

#include <Eigen/Sparse>

namespace rhbm_gem::algorithm {

struct SparseRegressionRow
{
    std::vector<std::pair<Eigen::Index, double>> basis_entries;
    double response{ 0.0 };
};

} // namespace rhbm_gem::algorithm
