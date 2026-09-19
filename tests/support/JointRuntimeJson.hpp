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
    j::object out{{"reference",Endpoint(t.reference)},{"passed",t.passed},{"reason",t.reason}};
    if(t.design) out["design_spectrum"]=DesignSpectrum(*t.design);
    if(!t.primary_valid) return out;
    out["scaled_coefficient_difference"]=Number(t.coefficient_difference); out["prediction_passed"]=t.prediction_passed;
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
inline j::object Assessment(const r::Assessment & a)
{
    j::object out{{"joint_qualified",a.qualified},{"primary",Endpoint(a.primary)},{"reference",Endpoint(a.reference)},
        {"endpoint_evaluations",2},{"directional_evaluations",2*a.derivatives.size()},{"qualification_failure",a.failure}};
    if(!a.design) return out;
    out["scaled_reference_difference"]=Number(a.coefficient_difference); out["design_spectrum"]=DesignSpectrum(*a.design);
    if(!a.widths) return out;
    auto widths=Spectrum(*a.widths); widths["column_norms"]=Values(a.widths->column_norms);
    widths["column_normalized"]=Spectrum(*a.normalized_widths); widths["active_face_only"]=!a.primary.certificate.active_atoms.empty();
    j::array weak; for(Eigen::Index k=0;k<a.weak_directions.cols();++k) weak.push_back(Values(a.weak_directions.col(k)));
    widths["weak_directions"]=weak; out["width_spectrum"]=widths;
    out["local_correction"]=Values(a.correction); out["local_correction_inf"]=Number(a.correction.lpNorm<Eigen::Infinity>());
    out["profile_jacobian_spectrum"]=Spectrum(*a.jacobian);
    j::array checks;
    for(const auto & d:a.derivatives) checks.push_back(j::object{{"direction",d.direction},{"h",d.h},{"relative_l2_difference",Number(d.error)},
        {"same_active_face",d.same_face},{"passed",d.passed},{"plus_valid",d.plus_valid},{"minus_valid",d.minus_valid}});
    out["derivative_checks"]=checks; out["derivative_verified"]=a.derivative_verified;
    out["qualification_checks"]=j::object{{"inner",a.inner},{"b_gradient",a.gradient},{"local_correction",a.local},
        {"identified",a.identified},{"derivative",a.derivative_verified}};
    return out;
}
inline j::object Search(const r::SearchResult & s,const r::EvaluationContext & c,Eigen::Index rows)
{
    j::array trials; for(const auto & t:s.trials) trials.push_back(Trial(t));
    j::object out{{"schema_version",1},{"experiment","joint-abc-profile"},{"alpha",0},{"execution_complete",true},{"joint_qualified",false},
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
