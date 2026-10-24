#pragma once
#include "support/JointTestNumerics.hpp"
#include <cstddef>
#include <optional>
#include <vector>

namespace rhbm_gem::core {class JointProblem; struct JointFitResult;}

namespace second_stage_test::matched::joint_abc {
rhbm_gem::core::JointFitResult FitBoundedJointComponentsForTesting(
    const rhbm_gem::core::JointProblem &,const std::vector<double> & initial_widths,double maximum_width);
struct FixedNeighborFitDiagnostic
{
    runtime::SearchResult search;
    runtime::Vector candidate_eta,candidate_beta;
    double candidate_objective{};
    bool candidate_available{};
    runtime::Assessment assessment;
    std::optional<runtime::Assessment> trusted_assessment;
    std::optional<runtime::Endpoint> trusted_state;
    std::optional<std::size_t> trusted_trial;
    double assessment_seconds{};
    std::optional<runtime::TrustEvidence> endpoint_trust;
    bool search_success{};
};
using ComponentView=runtime::ComponentView;
using ComponentPartition=runtime::ComponentPartition;
ComponentPartition BuildPartition(const Domain &, const std::vector<std::string> & atom_ids);
Eigen::VectorXd Select(const Eigen::VectorXd &, const std::vector<Eigen::Index> &);
EvaluationContext ComponentContext(const EvaluationContext &, const ComponentView &, bool independent);
boost::json::object Census(const Domain &, const ComponentPartition &, const EvaluationContext &);
Evaluation EvaluatePartitioned(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & eta,
    const ComponentPartition &, const EvaluationContext &);
Evaluation ComponentEvaluation(const Evaluation &, const ComponentView &);
Differential DifferentiatePartitioned(const Evaluation &, const ComponentPartition &, const EvaluationContext &);
boost::json::object SameState(const Domain &, const Eigen::VectorXd &, const Eigen::VectorXd & eta,
    const Eigen::VectorXd & beta, const ComponentPartition &, const EvaluationContext &);
runtime::FixedNeighborComponentResult SolveFixedNeighborComponentForTesting(const ComponentView &,
    const Eigen::VectorXd & parent_initial_b, const EvaluationContext &);
FixedNeighborFitDiagnostic SolveFixedNeighborComponentForDiagnostics(const ComponentView &,
    const Eigen::VectorXd & parent_initial_b, const EvaluationContext &);
boost::json::object FitComponent(const ComponentView &, const Eigen::VectorXd & parent_initial_b,
    const EvaluationContext &);
boost::json::object Assemble(const Domain &, const Eigen::VectorXd &, const ComponentPartition &,
    const EvaluationContext &, const boost::json::array & component_fits);
boost::json::object FitComponents(std::shared_ptr<const rhbm_gem::core::JointProblemInput>,
    const Eigen::VectorXd & initial_b, const ComponentPartition &, const EvaluationContext &);
} // namespace second_stage_test::matched::joint_abc
