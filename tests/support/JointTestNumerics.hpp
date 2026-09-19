#pragma once
#include "support/JointTestGeometry.hpp"
#include "support/JointTestContext.hpp"

namespace second_stage_test::matched::joint_abc {
using Sparse=Eigen::SparseMatrix<double>;
using Support=runtime::Support;
struct Domain : runtime::Domain
{
    using runtime::Domain::Domain;
    Domain(const runtime::Domain & d):runtime::Domain(d) {}
    Domain(const unique_grid::Grid &,const std::vector<Atom> &);
};
struct Evaluation : runtime::Evaluation
{
    boost::json::object certificate;
    Evaluation()=default;
    explicit Evaluation(runtime::Evaluation e);
};
using Differential=runtime::Differential;
Evaluation Evaluate(const Domain &, const Eigen::VectorXd & y,
    const Eigen::VectorXd & eta, bool reference=false, const EvaluationContext * = nullptr,
    const std::vector<runtime::LinearBlock> * = nullptr);
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
    const EvaluationContext * = nullptr);
} // namespace second_stage_test::matched::joint_abc
