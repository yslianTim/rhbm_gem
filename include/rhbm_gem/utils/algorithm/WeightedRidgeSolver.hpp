#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <rhbm_gem/utils/algorithm/WeightedRidgeNormalEquation.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSystem.hpp>

namespace rhbm_gem::algorithm {

class WeightedRidgeSolver
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_solver;
    bool m_analysis_success{ false };

public:
    explicit WeightedRidgeSolver(const WeightedRidgeSystem & system)
    {
        const Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
        auto equation{ BuildWeightedRidgeNormalEquation(system, weight) };
        m_solver.analyzePattern(equation.normal_matrix);
        m_analysis_success = m_solver.info() == Eigen::Success;
    }

    bool Solve(
        const WeightedRidgeSystem & system,
        const Eigen::VectorXd & weight,
        Eigen::VectorXd & parameter)
    {
        if (!m_analysis_success) return false;

        auto equation{ BuildWeightedRidgeNormalEquation(system, weight) };
        m_solver.factorize(equation.normal_matrix);
        if (m_solver.info() != Eigen::Success) return false;
        parameter = m_solver.solve(equation.right_hand_side);
        return m_solver.info() == Eigen::Success && parameter.allFinite();
    }
    
};

} // namespace rhbm_gem::algorithm
