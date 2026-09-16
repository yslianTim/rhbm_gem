#pragma once

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::mdpde_detail {

inline constexpr int kEndpointEquationBudget{ 128 };

struct MDPDEEquationEvidence
{
    Eigen::VectorXd raw, scaled, weights, singular_values;
    double denominator{}, condition{};
    int rank{}, floor_count{};
    bool valid{};
    std::string reason;
};

MDPDEEquationEvidence EvaluateMDPDEEquations(const RHBMMemberDataset &,
    double alpha, const Eigen::VectorXd & beta, double variance, double floor);
Eigen::VectorXd CalculateMDPDEBeta(const RHBMMemberDataset &, const Eigen::VectorXd & weights);
double CalculateMDPDEVariance(const RHBMMemberDataset &, double alpha,
    const Eigen::VectorXd & weights, const Eigen::VectorXd & beta);
RHBMDiagonalMatrix CalculateMDPDECovariance(double variance, const Eigen::VectorXd & weights);

struct RootEvaluation
{
    Eigen::VectorXd u;
    std::optional<double> residual_inf{};
    double denominator{};
    int floor_count{};
    std::string reason;
};

struct MDPDERootResult
{
    Eigen::VectorXd initial_u, u, beta;
    double variance{};
    MDPDEEquationEvidence equations;
    int native_status{}, iterations{}, equation_evaluations{}, verification_equation_evaluations{}, jacobian_evaluations{};
    std::string stop;
    std::optional<double> jacobian_condition{}, estimated_remaining_u_error{};
};

MDPDERootResult SolveMDPDERoot(const RHBMMemberDataset &, double alpha, double floor,
    const Eigen::VectorXd & initial, int total_budget, std::vector<RootEvaluation> * trace = nullptr);

struct MDPDEBranchComparison
{
    bool pass{};
    std::optional<std::array<double, 3>> relative_coordinate_difference{};
    std::optional<double> weight_max_difference{};
    std::optional<bool> floor_masks_equal{};
};

MDPDEBranchComparison CompareMDPDEBranches(const Eigen::VectorXd & candidate_beta, double candidate_variance,
    const MDPDEEquationEvidence & candidate, const Eigen::VectorXd & reference_beta, double reference_variance,
    const MDPDEEquationEvidence & reference, double floor);

// Detailed evidence is optional; production retains only the bounded scalar diagnostics.
struct EndpointRefinementEvidence
{
    MDPDEEquationEvidence original;
    std::optional<MDPDERootResult> root{};
    std::optional<MDPDEEquationEvidence> reference{};
    Eigen::VectorXd reference_beta;
    double reference_variance{};
    std::string reference_stop;
    std::optional<MDPDEBranchComparison> branch{};
    int candidate_linear_solves{};
};

RHBMBetaEstimateResult RefineMDPDEEndpoint(const RHBMMemberDataset &, double alpha,
    const RHBMExecutionOptions &, const RHBMBetaEstimateResult & endpoint,
    int equation_budget = kEndpointEquationBudget,
    const MDPDEEquationEvidence * initial = nullptr,
    EndpointRefinementEvidence * evidence = nullptr, std::vector<RootEvaluation> * trace = nullptr);

// Policy entry point: successful native solves do no additional equation work.
RHBMBetaEstimateResult ApplyFailedOnlyRefinement(const RHBMMemberDataset &, double alpha,
    const RHBMExecutionOptions &, const RHBMBetaEstimateResult & endpoint);

} // namespace rhbm_gem::mdpde_detail
