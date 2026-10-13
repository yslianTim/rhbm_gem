#pragma once
#include "core/detail/joint_component/Numerics.hpp"
#include <boost/json.hpp>
#include <cmath>

namespace second_stage_test::matched::runtime_json {
namespace r=rhbm_gem::core::joint_component;
namespace j=boost::json;
inline j::value Number(double x) {return std::isfinite(x) ? j::value(x) : j::value(nullptr);}
inline j::array Values(const Eigen::VectorXd & v) {j::array out; for(double x:v) out.push_back(Number(x)); return out;}
inline j::array Indices(const std::vector<Eigen::Index> & v) {j::array out; for(auto x:v) out.push_back(x); return out;}
inline j::object Certificate(const r::Certificate & c)
{
    if(!c.evaluated) return {};
    j::object out{{"feasible",c.feasible},{"kkt_passed",c.kkt_passed},{"projected_kkt",Number(c.projected_kkt)}};
    if(c.available)
    {
        out["rss"]=Number(c.rss); out["objective"]=Number(c.objective); out["residual_scale"]=Number(c.residual_scale);
        out["residual_rmse"]=Number(c.residual_rmse); out["residual_max"]=Number(c.residual_max);
        out["relative_residual"]=Number(c.relative_residual); out["active_atoms"]=Indices(c.active_atoms);
    }
    if(c.linear_solves) out["linear_solves"]=*c.linear_solves;
    if(c.free_rank) out["free_rank"]=*c.free_rank;
    if(c.block_factorizations) out["block_factorizations"]=*c.block_factorizations;
    return out;
}
inline j::object Endpoint(const r::Endpoint & e)
{
    auto out=Certificate(e.certificate);
    out["valid"]=e.valid; out["reason"]=e.reason; out["beta"]=Values(e.beta);
    out["eta"]=Values(e.eta); out["b"]=Values(e.eta.array().exp()); out["b_gradient"]=Values(e.gradient);
    out["b_gradient_inf"]=e.gradient.size() ? Number(e.gradient.lpNorm<Eigen::Infinity>()) : j::value(nullptr);
    return out;
}
inline j::object DesignSpectrum(const r::Spectrum & s)
{return {{"rank",s.rank},{"minimum_singular",Number(s.minimum)},{"condition",Number(s.condition)}};}
inline j::object Spectrum(const r::Spectrum & s)
{
    auto out=DesignSpectrum(s); out["singular_values"]=Values(s.singular_values); out["rank_threshold"]=Number(s.threshold); return out;
}
inline j::object StructuralRank(const r::StructuralRankEvidence & s)
{
    return {{"available",s.available},{"matching_available",s.matching_available},
        {"proven_deficient",s.proven_deficient},{"structurally_full_rank",s.structurally_full_rank},
        {"rows",s.rows},{"columns",s.columns},{"nonzero_entries",s.nonzero_entries},
        {"supported_rows",s.supported_rows},{"supported_columns",s.supported_columns},
        {"matching_size",s.matching_size},{"rank_upper_bound",s.rank_upper_bound},
        {"identical_support_pattern_groups",s.identical_support_pattern_groups},
        {"memory_estimate_bytes",s.memory_estimate_bytes},{"reason",s.reason}};
}
inline j::object MatrixSpectrum(const r::Spectrum & s)
{
    if(!s.available) return {{"available",false},{"reason",s.reason}};
    if(s.reason=="unobserved-columns") return {{"available",true},{"singular_values",j::array{}},{"rank",0},
        {"rank_threshold",0},{"global_rows",s.rows},{"global_columns",s.columns},{"reason",s.reason}};
    auto out=Spectrum(s); out["available"]=true; out["global_rows"]=s.rows; out["global_columns"]=s.columns;
    out["column_norms"]=Values(s.column_norms); return out;
}
inline j::object Trust(const r::TrustEvidence & t)
{
    j::object out{{"reference",t.reference ? j::value(Endpoint(*t.reference)) : j::value(nullptr)},
        {"passed",t.passed},{"reason",t.reason}};
    if(t.independent)
    {
        const auto & replay=*t.independent;
        out["reference_solve_calls"]=0;
        out["independent_replay"]=j::object{{"available",replay.available},{"valid",replay.valid},
            {"passed",replay.passed},{"feasible",replay.feasible},
            {"prediction_difference",Number(replay.prediction_difference)},
            {"residual_difference",Number(replay.residual_difference)},
            {"objective_difference",Number(replay.objective_difference)},
            {"primary_ac_kkt",Number(replay.primary_ac_kkt)},
            {"independent_ac_kkt",Number(replay.independent_ac_kkt)},
            {"kkt_difference",Number(replay.kkt_difference)},
            {"primary_width_gradient",Number(replay.primary_width_gradient)},
            {"independent_width_gradient",Number(replay.independent_width_gradient)},
            {"gradient_difference",Number(replay.gradient_difference)},
            {"cancellation_ratio",Number(replay.cancellation_ratio)},{"reason",replay.reason}};
    }
    out["scaled_coefficient_difference"]=Number(t.coefficient_difference);
    if(t.design) out["design_spectrum"]=DesignSpectrum(*t.design);
    if(!t.primary_valid) return out;
    out["prediction_passed"]=t.prediction_passed;
    out["gradient_passed"]=t.gradient_passed; out["kkt_replay_difference"]=Number(t.kkt_difference);
    out["maximum_prediction_difference"]=Number(t.prediction_difference);
    out["maximum_gradient_difference"]=Number(t.gradient_difference); out["maximum_cancellation_ratio"]=Number(t.cancellation_ratio); return out;
}
inline j::object Trial(const r::Trial & t)
{
    auto out=Endpoint(t.endpoint); out["evaluation"]=t.evaluation; out["accepted"]=t.accepted; out["seconds"]=t.seconds;
    if(t.accepted_update) out["accepted_update"]=*t.accepted_update;
    if(t.trust) out["trust"]=Trust(*t.trust);
    if(t.lm)
    {
        const auto & v=*t.lm;
        out["lm"]=j::object{{"accepted_eta",Values(v.accepted_eta)},{"step",Values(v.step)},{"diagonal",Values(v.diagonal)},
            {"radius",Number(v.radius)},{"damping",Number(v.damping)},{"actual_decrease",Number(v.actual_decrease)},
            {"predicted_decrease",Number(v.predicted_decrease)},{"ratio",Number(v.ratio)},{"proposed_acceptance",v.proposed_acceptance}};
    }
    return out;
}
inline const char * Status(rhbm_gem::core::JointCheckStatus status)
{
    switch(status)
    {
    case rhbm_gem::core::JointCheckStatus::Passed: return "passed";
    case rhbm_gem::core::JointCheckStatus::Failed: return "failed";
    case rhbm_gem::core::JointCheckStatus::Unavailable: return "unavailable";
    case rhbm_gem::core::JointCheckStatus::NotRun: return "not-run";
    }
    return "unavailable";
}
inline rhbm_gem::core::JointCheckStatus ParseStatus(const j::value & value)
{
    using S=rhbm_gem::core::JointCheckStatus;
    for(const auto status:{S::Passed,S::Failed,S::Unavailable,S::NotRun})
        if(value==Status(status)) return status;
    throw std::invalid_argument("Unknown runtime convergence status.");
}
inline j::object Assessment(const r::Assessment & a)
{
    const auto evidence=r::AssessmentEvidence(a,rhbm_gem::core::JointEvidenceScope::ComponentLocal);
    j::object out{{"runtime_convergence",Status(a.primary.valid ? r::ConvergenceStatus(evidence,rhbm_gem::core::JointEvidenceScope::ComponentLocal) : rhbm_gem::core::JointCheckStatus::Unavailable)},
        {"primary",Endpoint(a.primary)},{"directional_evaluations",0},{"runtime_failure",a.failure}};
    if(a.reference_free)
    {
        const auto & replay=a.independent_replay;
        out["reference"]=nullptr;
        out["endpoint_evaluations"]=1;
        out["reference_solve_calls"]=0;
        out["independent_replay"]=j::object{{"available",replay.available},{"valid",replay.valid},
            {"passed",replay.passed},{"feasible",replay.feasible},
            {"prediction_difference",Number(replay.prediction_difference)},
            {"residual_difference",Number(replay.residual_difference)},
            {"objective_difference",Number(replay.objective_difference)},
            {"primary_ac_kkt",Number(replay.primary_ac_kkt)},
            {"independent_ac_kkt",Number(replay.independent_ac_kkt)},
            {"kkt_difference",Number(replay.kkt_difference)},
            {"primary_width_gradient",Number(replay.primary_width_gradient)},
            {"independent_width_gradient",Number(replay.independent_width_gradient)},
            {"gradient_difference",Number(replay.gradient_difference)},
            {"cancellation_ratio",Number(replay.cancellation_ratio)},{"reason",replay.reason}};
        out["runtime_checks"]=j::object{{"ac_feasibility",replay.available && replay.feasible},
            {"conditional_ac_optimality",replay.primary_kkt_passed && replay.feasible},
            {"independent_ac_optimality",replay.independent_kkt_passed && replay.feasible},
            {"independent_kkt_consistency",replay.kkt_agrees},
            {"independent_replay",replay.valid && replay.passed},
            {"width_gradient_consistency",replay.gradient_agrees},
            {"width_stationarity",replay.primary_width_gradient<=1e-12},
            {"design_identifiability",a.design && a.design->rank==a.design->columns},
            {"numerical_identifiability",a.identified}};
    }
    else
    {
        out["reference"]=Endpoint(a.reference);
        out["endpoint_evaluations"]=2;
    }
    if(a.design_structure) out["design_structure"]=StructuralRank(*a.design_structure);
    if(a.design)
    {
        if(!a.reference_free) out["scaled_reference_difference"]=Number(a.coefficient_difference);
        out["design_spectrum"]=DesignSpectrum(*a.design);
    }
    if(a.widths)
    {
        auto widths=Spectrum(*a.widths); widths["column_norms"]=Values(a.widths->column_norms);
        widths["active_face_only"]=!a.primary.certificate.active_atoms.empty();
        j::array weak; for(Eigen::Index k=0;k<a.weak_directions.cols();++k) weak.push_back(Values(a.weak_directions.col(k)));
        widths["weak_directions"]=weak; out["width_spectrum"]=widths;
    }
    if(!a.reference_free) out["runtime_checks"]=j::object{{"inner",a.inner},{"b_gradient",a.gradient},
        {"identified",a.identified}};
    return out;
}
inline j::object Search(const r::SearchResult & s,const r::EvaluationContext & c,Eigen::Index rows)
{
    j::array trials; for(const auto & t:s.trials) trials.push_back(Trial(t));
    j::object out{{"schema_version",2},{"experiment","joint-abc-profile"},{"alpha",0},{"execution_complete",true},{"runtime_convergence","unavailable"},
        {"initial",Endpoint(s.initial)},{"lm_status",s.lm_status},{"stop_reason",s.stop_reason},{"profile_evaluations",s.evaluations},
        {"jacobian_evaluations",s.derivatives},{"accepted_updates",s.accepted},{"trials",trials},{"row_count",rows},
        {"settings",j::object{{"factor",.1},{"ftol",1e-14},{"xtol",1e-12},{"gtol",1e-12},{"profile_budget",c.profile_budget},{"accepted_update_budget",c.update_budget}}},
        {"linear_solver","sparse-qr-householder-1024"},{"reference_solver","independent-tsqr-8192-svd"},
        {"residual_scale",c.scale},{"variance_semantics","descriptive RSS/N; zero permitted"}};
    out["variant"]="guarded"; out["search_reference_evaluations"]=s.references;
    out["initial_accepted"]=s.initial_accepted; out["search_reference_seconds"]=s.reference_seconds;
    out["search_stopped_without_convergence"]=s.stopped;
    return out;
}
}
