#pragma once
#include "support/UniqueStencilGrid.hpp"

namespace second_stage_test::matched::joint_abc {
using Sparse=Eigen::SparseMatrix<double>;
struct Support {Eigen::Index row; double square;};
// Fixed observation domain and contributor support; no parameter truth.
struct Domain
{
    Eigen::Index rows;
    std::vector<std::vector<Support>> atoms;
    Domain(const unique_grid::Grid &, const std::vector<Atom> &);
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
    const Eigen::VectorXd & eta, bool reference=false);
Differential Differentiate(const Evaluation &, double scale);
boost::json::object Fit(const Domain &, const Eigen::VectorXd & y, const Eigen::VectorXd & initial_b,
    boost::json::object * resources = nullptr);
void Run(const std::string & manifest, const std::string & map,
    const std::string & checkpoint, const std::string & output);
} // namespace second_stage_test::matched::joint_abc
