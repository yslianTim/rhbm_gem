#pragma once

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include <Eigen/Sparse>

#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>

#include "core/detail/LocalFittingStateView.hpp"

namespace rhbm_gem::core::detail {

class ReusableWeightedRidgeSolver
{
public:
    bool Solve(
        const algorithm::WeightedRidgeSystem & system,
        const Eigen::VectorXd & weight,
        Eigen::VectorXd & parameter)
    {
        const auto pattern{ BuildPattern(system.design_matrix) };
        if (m_row_count != system.design_matrix.rows() ||
            m_column_count != system.design_matrix.cols() ||
            m_pattern != pattern)
        {
            if (!m_solver.AnalyzePattern(system)) return false;
            m_row_count = system.design_matrix.rows();
            m_column_count = system.design_matrix.cols();
            m_pattern = pattern;
        }
        return m_solver.SolveNumeric(system, weight, parameter);
    }

    std::size_t GetSymbolicAnalysisCount() const
    {
        return m_solver.GetSymbolicAnalysisCount();
    }

private:
    static std::vector<std::pair<Eigen::Index, Eigen::Index>> BuildPattern(
        const Eigen::SparseMatrix<double> & matrix)
    {
        std::vector<std::pair<Eigen::Index, Eigen::Index>> pattern;
        pattern.reserve(static_cast<std::size_t>(matrix.nonZeros()));
        for (Eigen::Index outer = 0; outer < matrix.outerSize(); outer++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator iter(matrix, outer);
                iter;
                ++iter)
            {
                pattern.emplace_back(iter.row(), iter.col());
            }
        }
        return pattern;
    }

    algorithm::WeightedRidgeSolver m_solver{};
    Eigen::Index m_row_count{ -1 };
    Eigen::Index m_column_count{ -1 };
    std::vector<std::pair<Eigen::Index, Eigen::Index>> m_pattern{};
};

struct LocalFittingClusterSolverWorkspace
{
    ReusableWeightedRidgeSolver joint_offset{};
    ReusableWeightedRidgeSolver joint_polish{};
};

using LocalFittingClusterSolverWorkspaceMap =
    std::map<LocalFittingClusterKey, LocalFittingClusterSolverWorkspace>;

} // namespace rhbm_gem::core::detail
