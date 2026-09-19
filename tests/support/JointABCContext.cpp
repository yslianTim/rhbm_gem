#include "support/JointABCContext.hpp"
#include <algorithm>
#include <limits>

namespace second_stage_test::matched::joint_abc {
double RankPolicy::Relative(Eigen::Index columns) const
{return std::numeric_limits<double>::epsilon()*static_cast<double>(std::max(rows,columns));}
double RankPolicy::Absolute(Eigen::Index columns,double maximum) const
{return Relative(columns)*maximum;}
AuditPlan RegisteredAudit(Eigen::Index atoms,const std::string & dataset,const std::string & name)
{
    AuditPlan plan; plan.trial_details=atoms<=12; plan.expanded_if_unverified=atoms<=12;
    plan.precision=atoms<=12 && (dataset=="weak-1e-4" || dataset=="near-0.02" || dataset=="active-a" ||
        (dataset=="baseline" && name.starts_with("first-stage")));
    plan.boundary=dataset=="active-a";
    if(plan.boundary) plan.boundary_atoms={1,5,9};
    return plan;
}
EvaluationContext MakeContext(const Eigen::VectorXd & y,Eigen::Index atoms,const std::string & hash,const AuditPlan * plan)
{
    EvaluationContext c; c.snapshot_hash=hash; c.observations=std::make_shared<const Eigen::VectorXd>(y);
    c.scale=std::max(1.0,y.norm()); c.rank={y.size(),2*atoms,atoms};
    c.linear.rank_relative=c.rank.Relative(2*atoms);
    c.linear.release_response_norm=y.norm();
    c.audit=plan ? *plan : RegisteredAudit(atoms);
    for(Eigen::Index a=0;a<atoms;++a) c.atom_ids.push_back(std::to_string(a));
    for(Eigen::Index r=0;r<y.size();++r) c.row_ids.push_back(std::to_string(r));
    return c;
}
boost::json::object ContextEvidence(const EvaluationContext & c)
{
    namespace j=boost::json;
    j::array atoms,rows,boundary,directions;
    for(const auto & a:c.atom_ids) atoms.emplace_back(a);
    for(const auto & r:c.row_ids) rows.emplace_back(r);
    for(auto a:c.audit.boundary_atoms) boundary.push_back(a);
    for(Eigen::Index k=0;k<c.audit.directions.cols();++k)
    {j::array values; for(double x:c.audit.directions.col(k)) values.push_back(x); directions.push_back(values);}
    return {{"schema_version",1},{"parent_snapshot_sha256",c.snapshot_hash},{"observation_scale",c.scale},
        {"global_rows",c.observations ? c.observations->size() : c.rank.rows},
        {"rank_rows",c.rank.rows},{"design_columns",c.rank.design_columns},{"width_columns",c.rank.width_columns},
        {"atom_ids",atoms},{"row_ids",rows},{"search_policy",c.independent_search ? "component-local-linear-policy-parent-observation-scale" : "global-linear-policy"},
        {"svd_threshold","epsilon * max(rank_rows, family_columns) * current_global_sigma_max"},
        {"sparse_qr_threshold","epsilon * max(rank_rows, design_columns) * maximum_normalized_column_norm"},
        {"active_set_release_factor",c.linear.release_factor},{"active_set_iteration_factor",c.linear.active_set_iteration_factor},
        {"active_set_response_norm",c.linear.release_response_norm},
        {"active_set_release_rule","factor * epsilon * max(1, norm(Z) * (parent_response_norm + norm(Z) * norm(scaled_beta)))"},
        {"profile_budget",c.profile_budget},{"accepted_update_budget",c.update_budget},
        {"audit",j::object{{"trial_details",c.audit.trial_details},{"expanded_if_unverified",c.audit.expanded_if_unverified},
            {"precision",c.audit.precision},{"block_precision_reference",c.audit.block_precision},
            {"precision_cache",c.audit.cache_precision ? "exact-input, current-audit-case-only" : "none"},
            {"directions",directions},{"boundary",c.audit.boundary},{"boundary_atoms",boundary}}}};
}
} // namespace second_stage_test::matched::joint_abc
