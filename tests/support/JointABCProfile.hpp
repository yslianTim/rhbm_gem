#pragma once
#include "support/UniqueStencilGrid.hpp"
#include "support/JointABCContext.hpp"

namespace second_stage_test::matched::joint_abc {
using Sparse=Eigen::SparseMatrix<double>;
struct Support {Eigen::Index row; double square;};
// Fixed observation domain and contributor support; no parameter truth.
struct Domain
{
    Eigen::Index rows;
    std::vector<std::vector<Support>> atoms;
    Domain(const unique_grid::Grid &, const std::vector<Atom> &);
    Domain(Eigen::Index count, std::vector<std::vector<Support>> support)
        :rows(count),atoms(std::move(support)) {}
};
struct Evaluation
{
    Eigen::VectorXd eta, beta, residual, gradient;
    Sparse x, derivative;
    boost::json::object certificate;
    bool valid{};
    std::string reason;
};
struct Differential
{
    Eigen::MatrixXd projected, jacobian;
    bool valid{};
    std::string reason;
};
Evaluation Evaluate(const Domain &, const Eigen::VectorXd & y,
    const Eigen::VectorXd & eta, bool reference=false, const EvaluationContext * = nullptr,
    const std::vector<joint_ac::LinearBlock> * = nullptr);
Evaluation AtState(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & eta,
    const Eigen::VectorXd & beta, const EvaluationContext &);
boost::json::object MatrixSpectrum(const Sparse &, const RankPolicy &, Eigen::Index columns, bool normalize);
boost::json::object MatrixSpectrum(const Eigen::MatrixXd &, const RankPolicy &, Eigen::Index columns, bool normalize);
Eigen::VectorXd LocalCorrection(const Evaluation &, const Differential &, const EvaluationContext &,
    double absolute_threshold = -1);
Differential Differentiate(const Evaluation &, double scale, const EvaluationContext * = nullptr,
    double free_design_threshold = -1);
boost::json::object Trust(const Domain &, const Eigen::VectorXd &, const Evaluation &,
    const EvaluationContext * = nullptr);
boost::json::object Assess(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & eta,
    const EvaluationContext &, const Eigen::VectorXd * supplied_beta=nullptr);
boost::json::object Fit(const Domain &, const Eigen::VectorXd & y, const Eigen::VectorXd & initial_b,
    boost::json::object * resources = nullptr, const std::string & variant = "",
    const EvaluationContext * = nullptr);
void Run(const std::string & manifest, const std::string & map,
    const std::string & checkpoint, const std::string & output);
} // namespace second_stage_test::matched::joint_abc
