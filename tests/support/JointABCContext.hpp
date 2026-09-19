#pragma once
#include <Eigen/Core>
#include <boost/json.hpp>
#include <memory>
#include <string>
#include <vector>

namespace second_stage_test::matched::joint_abc {
// Dimensions belong to the observation problem, not to a QR reduction tile.
// Absolute SVD thresholds are resolved afresh from each matrix's spectrum.
struct RankPolicy
{
    Eigen::Index rows{}, design_columns{}, width_columns{};
    double Relative(Eigen::Index columns) const;
    double Absolute(Eigen::Index columns, double maximum_singular) const;
};
struct AuditPlan
{
    bool trial_details{}, expanded_if_unverified{}, precision{}, boundary{};
    bool block_precision{}, cache_precision{};
    std::vector<Eigen::Index> boundary_atoms;
    Eigen::MatrixXd directions;
};
struct LinearPolicy
{
    double rank_relative{};
    int active_set_iteration_factor{20};
    double release_factor{128};
    double release_response_norm{-1};
};
struct EvaluationContext
{
    std::string snapshot_hash;
    std::shared_ptr<const Eigen::VectorXd> observations;
    std::vector<std::string> atom_ids, row_ids;
    double scale{1};
    RankPolicy rank;
    LinearPolicy linear;
    AuditPlan audit;
    bool independent_search{};
    int profile_budget{200}, update_budget{100};
};
// This is the experiment registry. Numerical kernels never inspect case names.
AuditPlan RegisteredAudit(Eigen::Index atoms, const std::string & dataset = "",
    const std::string & case_name = "");
EvaluationContext MakeContext(const Eigen::VectorXd &, Eigen::Index atoms,
    const std::string & snapshot_hash = "", const AuditPlan * = nullptr);
boost::json::object ContextEvidence(const EvaluationContext &);
} // namespace second_stage_test::matched::joint_abc
