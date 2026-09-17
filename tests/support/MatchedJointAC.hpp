#pragma once
#include "support/ObservationMatchedExperiment.hpp"
#include <Eigen/SparseCore>

namespace second_stage_test::matched::joint_ac {
struct LinearResult
{
    Eigen::VectorXd beta;
    bool valid{};
    std::string reason;
    int rank{}, solves{}, releases{};
};
struct Block
{
    std::size_t owner{};
    double alpha{};
    std::vector<Eigen::Index> rows;
};
using Blocks = std::vector<Block>;
struct Evidence
{
    Eigen::VectorXd weights, linear_weights, prefactors, denominators, scaled;
    double objective{}, stationarity{};
    bool valid{};
    std::string reason;
    int failure_owner{-1};
};
Eigen::VectorXd BlockVariances(const Eigen::VectorXd & residual, const Blocks &);
double Objective(const Eigen::VectorXd & residual, const Eigen::VectorXd & variances, const Blocks &);
LinearResult WeightedSolve(const Eigen::MatrixXd &, const Eigen::VectorXd &,
    const Eigen::VectorXd & weights, bool svd = false, bool blocked_svd = false, const Eigen::SparseMatrix<double> * sparse_design = nullptr);
Evidence Evaluate(const Eigen::MatrixXd &, const Eigen::VectorXd &,
    const Eigen::VectorXd & beta, const Eigen::VectorXd & variances, const Blocks &);
boost::json::object Fit(const Eigen::MatrixXd &, const Eigen::VectorXd &,
    const Eigen::VectorXd & initial, const Blocks &, int budget = 2000, int reference_budget = 4000,
    bool blocked_svd = false, const Eigen::SparseMatrix<double> * sparse_design = nullptr);
struct Components
{
    std::vector<std::vector<std::size_t>> atoms, rows, contributors;
};
Components BuildComponents(const std::vector<Stencil> &, const std::vector<Atom> &, double cutoff,
    const std::vector<std::size_t> & owners);
Eigen::MatrixXd BuildDesign(const std::vector<Stencil> &, const std::vector<Atom> &,
    const std::vector<std::size_t> & rows, const std::vector<std::size_t> & members, double cutoff);
void Run(const std::string & manifest, const std::string & map,
    const std::string & state_index, const std::string & output);
} // namespace second_stage_test::matched::joint_ac
