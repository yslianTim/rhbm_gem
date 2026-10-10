#include "support/JointTestContext.hpp"
#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>

namespace second_stage_test::matched::joint_abc {
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
{return runtime::CreateContext(y,atoms,hash,plan ? *plan : RegisteredAudit(atoms));}
std::shared_ptr<rhbm_gem::core::JointProblemInput> MakeSyntheticJointProblemInput(
    const runtime::Domain & domain,const Eigen::VectorXd & observations,
    std::vector<std::string> atom_ids,std::vector<std::string> row_ids)
{
    if(domain.rows!=observations.size() || !observations.allFinite())
        throw std::invalid_argument("Invalid synthetic joint observations.");
    if(atom_ids.empty()) for(std::size_t atom=0;atom<domain.atoms.size();++atom)
        atom_ids.push_back("atom-"+std::to_string(atom));
    if(row_ids.empty()) for(Eigen::Index row=0;row<domain.rows;++row)
        row_ids.push_back("row-"+std::to_string(row));
    if(atom_ids.size()!=domain.atoms.size() || row_ids.size()!=static_cast<std::size_t>(domain.rows) ||
        std::set<std::string>(atom_ids.begin(),atom_ids.end()).size()!=atom_ids.size() ||
        std::set<std::string>(row_ids.begin(),row_ids.end()).size()!=row_ids.size())
        throw std::invalid_argument("Invalid synthetic joint identities.");
    auto input=std::make_shared<rhbm_gem::core::JointProblemInput>();
    input->observations.assign(observations.data(),observations.data()+observations.size());
    input->atom_ids.reserve(domain.atoms.size()); input->support.resize(domain.atoms.size());
    for(std::size_t atom=0;atom<domain.atoms.size();++atom)
    {
        input->atom_ids.push_back(std::move(atom_ids[atom]));
        for(const auto & support:domain.atoms[atom])
        {
            if(support.row<0 || support.row>=domain.rows || !std::isfinite(support.square) ||
                support.square<0 || support.square>6.25)
                throw std::invalid_argument("Invalid synthetic joint support.");
            input->support[atom].push_back({static_cast<std::size_t>(support.row),support.square});
        }
    }
    input->row_ids=std::move(row_ids);
    return input;
}
EvaluationContext MakeContext(std::shared_ptr<const rhbm_gem::core::JointProblemInput> input,
    const std::string & hash,const AuditPlan * plan)
{
    const auto atoms=static_cast<Eigen::Index>(input->atom_ids.size());
    return runtime::CreateContext(std::move(input),hash,plan ? *plan : RegisteredAudit(atoms));
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
