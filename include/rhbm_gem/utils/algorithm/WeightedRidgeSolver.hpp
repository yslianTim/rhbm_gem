#pragma once

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace rhbm_gem::algorithm {

struct WeightedRidgeSystem
{
    Eigen::SparseMatrix<double> design_matrix;
    Eigen::VectorXd response;
    Eigen::VectorXd previous_parameter;
    Eigen::VectorXd ridge_diagonal;
};

class WeightedRidgeSolver
{
    struct NormalEquation
    {
        Eigen::SparseMatrix<double> normal_matrix;
        Eigen::VectorXd right_hand_side;
    };

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solver;
    bool m_analysis_success{ false };
    std::size_t m_symbolic_analysis_count{ 0 };
    Eigen::Index m_row_count{ -1 };
    Eigen::Index m_column_count{ -1 };
    std::vector<std::pair<Eigen::Index, Eigen::Index>> m_pattern{};

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

    static NormalEquation BuildNormalEquation(
        const WeightedRidgeSystem & system,
        const Eigen::VectorXd & weight)
    {
        auto weighted_design{ system.design_matrix };
        for (Eigen::Index column = 0; column < weighted_design.outerSize(); column++)
        {
            for (Eigen::SparseMatrix<double>::InnerIterator iter(weighted_design, column);
                 iter;
                 ++iter)
            {
                iter.valueRef() *= std::sqrt(weight(iter.row()));
            }
        }
        const auto weighted_response{
            weight.array().sqrt().matrix().cwiseProduct(system.response)
        };
        Eigen::SparseMatrix<double> normal_matrix{
            weighted_design.transpose() * weighted_design
        };
        const auto column_count{ normal_matrix.cols() };
        for (Eigen::Index i = 0; i < column_count; i++)
        {
            normal_matrix.coeffRef(i, i) += system.ridge_diagonal(i);
        }
        normal_matrix.makeCompressed();
        Eigen::VectorXd right_hand_side{
            weighted_design.transpose() * weighted_response +
            system.ridge_diagonal.cwiseProduct(system.previous_parameter)
        };

        return { std::move(normal_matrix), std::move(right_hand_side) };
    }

    bool AnalyzePatternAndCache(
        const WeightedRidgeSystem & system,
        std::vector<std::pair<Eigen::Index, Eigen::Index>> pattern)
    {
        const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
        auto equation{ BuildNormalEquation(system, weight) };
        m_solver.analyzePattern(equation.normal_matrix);
        m_analysis_success = m_solver.info() == Eigen::Success;
        m_symbolic_analysis_count++;
        if (m_analysis_success)
        {
            m_row_count = system.design_matrix.rows();
            m_column_count = system.design_matrix.cols();
            m_pattern = std::move(pattern);
        }
        return m_analysis_success;
    }

public:
    WeightedRidgeSolver() = default;

    explicit WeightedRidgeSolver(const WeightedRidgeSystem & system)
    {
        AnalyzePattern(system);
    }

    bool AnalyzePattern(const WeightedRidgeSystem & system)
    {
        return AnalyzePatternAndCache(system, BuildPattern(system.design_matrix));
    }

    bool SolveNumeric(
        const WeightedRidgeSystem & system,
        const Eigen::VectorXd & weight,
        Eigen::VectorXd & parameter)
    {
        if (!m_analysis_success) return false;

        auto equation{ BuildNormalEquation(system, weight) };
        m_solver.factorize(equation.normal_matrix);
        if (m_solver.info() != Eigen::Success) return false;
        parameter = m_solver.solve(equation.right_hand_side);
        return m_solver.info() == Eigen::Success && parameter.allFinite();
    }

    bool Solve(
        const WeightedRidgeSystem & system,
        const Eigen::VectorXd & weight,
        Eigen::VectorXd & parameter)
    {
        auto pattern{ BuildPattern(system.design_matrix) };
        if (!m_analysis_success ||
            m_row_count != system.design_matrix.rows() ||
            m_column_count != system.design_matrix.cols() ||
            m_pattern != pattern)
        {
            if (!AnalyzePatternAndCache(system, std::move(pattern))) return false;
        }
        return SolveNumeric(system, weight, parameter);
    }

    std::size_t GetSymbolicAnalysisCount() const
    {
        return m_symbolic_analysis_count;
    }
};

} // namespace rhbm_gem::algorithm
