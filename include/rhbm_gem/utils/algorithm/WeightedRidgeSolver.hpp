#pragma once

#include <cmath>
#include <utility>

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

public:
    explicit WeightedRidgeSolver(const WeightedRidgeSystem & system)
    {
        const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
        auto equation{ BuildNormalEquation(system, weight) };
        m_solver.analyzePattern(equation.normal_matrix);
        m_analysis_success = m_solver.info() == Eigen::Success;
    }

    bool Solve(
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
};

} // namespace rhbm_gem::algorithm
