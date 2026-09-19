#pragma once
#include "support/ObservationMatchedExperiment.hpp"
#include <Eigen/SparseCore>
#include "support/JointABCContext.hpp"

namespace second_stage_test::matched::joint_ac {
struct LinearResult
{
    Eigen::VectorXd beta;
    bool valid{};
    std::string reason;
    int rank{}, solves{}, releases{};
    int block_factorizations{};
};
struct LinearBlock {std::vector<Eigen::Index> rows, columns;};
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
    // Composite memberships are ordered by block, then by its row list.
    Eigen::VectorXd membership_weights, block_prefactors, log_prefactors;
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
// Sparse-only experiment path. Dense production-test entry points stay unchanged.
LinearResult WeightedSolve(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &, bool svd = false, bool blocked_svd = true,
    const Eigen::SparseMatrix<double> * sparse_design = nullptr,
    const joint_abc::LinearPolicy * = nullptr, const std::vector<LinearBlock> * = nullptr);
Evidence Evaluate(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &, const Eigen::VectorXd &, const Blocks &);
boost::json::object Fit(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &, const Blocks &, int budget = 100, int reference_budget = 100,
    bool blocked_svd = true, const Eigen::SparseMatrix<double> * sparse_design = nullptr);
boost::json::object SparseSpectrum(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd & weights);
Evidence EvaluateComposite(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &, const Eigen::VectorXd &, const Blocks &);
boost::json::object FitComposite(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &, const Blocks &, int budget = 100, int reference_budget = 100);
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
