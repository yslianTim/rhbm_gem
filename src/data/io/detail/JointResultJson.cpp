#include "JointResultJson.hpp"
#include "core/detail/joint_component/ConvergenceCriteria.hpp"
#include <boost/json.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace rhbm_gem::joint_result_io {
namespace j = boost::json;
std::string_view StatusText(JointCheckStatus status)
{
    switch (status)
    {
    case JointCheckStatus::Passed: return "passed";
    case JointCheckStatus::Failed: return "failed";
    case JointCheckStatus::Unavailable: return "unavailable";
    case JointCheckStatus::NotRun: return "not-run";
    }
    throw std::invalid_argument("Invalid joint check status.");
}
namespace {
using Object = j::object;
template<class T> T Read(const Object & o, const char * key) {return j::value_to<T>(o.at(key));}
JointCheckStatus Status(const j::value & v)
{
    const auto text=j::value_to<std::string>(v);
    for (auto status:{JointCheckStatus::Passed,JointCheckStatus::Failed,JointCheckStatus::Unavailable,JointCheckStatus::NotRun})
        if (StatusText(status)==text) return status;
    throw std::invalid_argument("Invalid saved joint check status: "+text);
}
const char * ScopeText(JointEvidenceScope scope)
{
    switch (scope)
    {
    case JointEvidenceScope::ComponentLocal: return "component-local";
    case JointEvidenceScope::AssembledGlobal: return "assembled-global";
    }
    throw std::invalid_argument("Invalid joint evidence scope.");
}
JointEvidenceScope Scope(const Object & o)
{
    const auto s=Read<std::string>(o,"scope");
    if(s=="component-local") return JointEvidenceScope::ComponentLocal;
    if(s=="assembled-global") return JointEvidenceScope::AssembledGlobal;
    throw std::invalid_argument("Invalid saved joint evidence scope: "+s);
}
const char * HaloMedianSourceText(HaloMedianSource source)
{
    switch(source)
    {
    case HaloMedianSource::GroupKeyMedian: return "group-key-median";
    case HaloMedianSource::GlobalTargetMedian: return "global-target-median";
    }
    throw std::invalid_argument("Invalid Halo median source.");
}
HaloMedianSource HaloMedianSourceValue(const j::value & value)
{
    const auto source=j::value_to<std::string>(value);
    if(source=="group-key-median") return HaloMedianSource::GroupKeyMedian;
    if(source=="global-target-median") return HaloMedianSource::GlobalTargetMedian;
    throw std::invalid_argument("Invalid saved Halo median source: "+source);
}
j::value Number(std::optional<double> x)
{
    if(x && !std::isfinite(*x)) throw std::invalid_argument("Nonfinite joint result number.");
    return x ? j::value(*x) : j::value(nullptr);
}
std::optional<double> OptionalNumber(const j::value & v)
{return v.is_null() ? std::nullopt : std::optional<double>(j::value_to<double>(v));}
j::value OptionalBoolean(std::optional<bool> value)
{return value ? j::value(*value) : j::value(nullptr);}
std::optional<bool> ReadOptionalBoolean(const j::value & value)
{return value.is_null() ? std::nullopt : std::optional<bool>(j::value_to<bool>(value));}
const char * QualificationText(RHBMSolveQualification qualification)
{
    switch(qualification)
    {
    case RHBMSolveQualification::NativeSuccess: return "NativeSuccess";
    case RHBMSolveQualification::RefinedSuccess: return "RefinedSuccess";
    case RHBMSolveQualification::Unqualified: return "Unqualified";
    }
    throw std::invalid_argument("Invalid FirstStage solve qualification.");
}
RHBMSolveQualification Qualification(const j::value & value)
{
    const auto text=j::value_to<std::string>(value);
    if(text=="NativeSuccess") return RHBMSolveQualification::NativeSuccess;
    if(text=="RefinedSuccess") return RHBMSolveQualification::RefinedSuccess;
    if(text=="Unqualified") return RHBMSolveQualification::Unqualified;
    throw std::invalid_argument("Invalid saved FirstStage solve qualification.");
}
j::value EndpointRefinement(const RHBMEndpointRefinementDiagnostics & r)
{
    j::value coordinates=nullptr;
    if(r.relative_coordinate_difference)
        coordinates=j::array{(*r.relative_coordinate_difference)[0],(*r.relative_coordinate_difference)[1],
            (*r.relative_coordinate_difference)[2]};
    return Object{{"accepted",r.accepted},{"reason",r.reason},
        {"candidate_equation_evaluations",r.candidate_equation_evaluations},
        {"reference_updates",r.reference_updates},{"reference_stop",r.reference_stop},
        {"original_residual",Number(r.original_residual)},{"candidate_residual",Number(r.candidate_residual)},
        {"reference_residual",Number(r.reference_residual)},{"relative_coordinate_difference",std::move(coordinates)},
        {"weight_max_difference",Number(r.weight_max_difference)},
        {"floor_masks_equal",OptionalBoolean(r.floor_masks_equal)}};
}
RHBMEndpointRefinementDiagnostics ReadEndpointRefinement(const j::value & value)
{
    const auto & o=value.as_object(); RHBMEndpointRefinementDiagnostics r;
    r.accepted=Read<bool>(o,"accepted"); r.reason=Read<std::string>(o,"reason");
    r.candidate_equation_evaluations=Read<int>(o,"candidate_equation_evaluations");
    r.reference_updates=Read<int>(o,"reference_updates"); r.reference_stop=Read<std::string>(o,"reference_stop");
    r.original_residual=OptionalNumber(o.at("original_residual"));
    r.candidate_residual=OptionalNumber(o.at("candidate_residual"));
    r.reference_residual=OptionalNumber(o.at("reference_residual"));
    if(!o.at("relative_coordinate_difference").is_null())
    {
        const auto & a=o.at("relative_coordinate_difference").as_array();
        if(a.size()!=3) throw std::invalid_argument("Invalid endpoint refinement coordinate shape.");
        r.relative_coordinate_difference=std::array<double,3>{j::value_to<double>(a[0]),
            j::value_to<double>(a[1]),j::value_to<double>(a[2])};
    }
    r.weight_max_difference=OptionalNumber(o.at("weight_max_difference"));
    r.floor_masks_equal=ReadOptionalBoolean(o.at("floor_masks_equal"));
    return r;
}
j::value FailedOnlyRecovery(const RHBMFailedOnlyRecoveryDiagnostics & r)
{
    j::value endpoint=nullptr, retry=nullptr;
    if(r.endpoint_refinement) endpoint=EndpointRefinement(*r.endpoint_refinement);
    if(r.bounded_native_retry)
    {
        const auto & value=*r.bounded_native_retry;
        retry=Object{{"max_iterations",value.max_iterations},{"native_status",value.native_status ? j::value(*value.native_status) : j::value(nullptr)},
            {"iterations",value.iterations},{"squared_beta_change",Number(value.squared_beta_change)},
            {"relative_variance_change",Number(value.relative_variance_change)},
            {"equation_residual",Number(value.equation_residual)},{"seconds",value.seconds}};
    }
    return Object{{"accepted",r.accepted},{"reason",r.reason},
        {"endpoint_refinement_seconds",r.endpoint_refinement_seconds},
        {"endpoint_refinement",std::move(endpoint)},{"bounded_native_retry",std::move(retry)}};
}
RHBMFailedOnlyRecoveryDiagnostics ReadFailedOnlyRecovery(const j::value & value)
{
    const auto & o=value.as_object(); RHBMFailedOnlyRecoveryDiagnostics r;
    r.accepted=Read<bool>(o,"accepted"); r.reason=Read<std::string>(o,"reason");
    r.endpoint_refinement_seconds=Read<double>(o,"endpoint_refinement_seconds");
    if(!o.at("endpoint_refinement").is_null()) r.endpoint_refinement=ReadEndpointRefinement(o.at("endpoint_refinement"));
    if(!o.at("bounded_native_retry").is_null())
    {
        const auto & retry=o.at("bounded_native_retry").as_object(); RHBMNativeRetryDiagnostics x;
        x.max_iterations=Read<int>(retry,"max_iterations");
        if(!retry.at("native_status").is_null()) x.native_status=Read<int>(retry,"native_status");
        x.iterations=Read<int>(retry,"iterations"); x.squared_beta_change=OptionalNumber(retry.at("squared_beta_change"));
        x.relative_variance_change=OptionalNumber(retry.at("relative_variance_change"));
        x.equation_residual=OptionalNumber(retry.at("equation_residual")); x.seconds=Read<double>(retry,"seconds");
        r.bounded_native_retry=std::move(x);
    }
    return r;
}
// Initialization diagnostics may contain failed/nonfinite estimates. These are
// explicitly null on disk; a decoded null remains unavailable (NaN), never zero.
template<class Range> j::array Diagnostics(const Range & values)
{
    j::array out;
    for(double x:values) out.push_back(std::isfinite(x) ? j::value(x) : j::value(nullptr));
    return out;
}
std::vector<double> DiagnosticValues(const j::value & v)
{
    std::vector<double> out;
    for(const auto & x:v.as_array()) out.push_back(x.is_null() ? std::numeric_limits<double>::quiet_NaN() : j::value_to<double>(x));
    return out;
}
Object State(const JointState & x)
{
    return {{"ac",j::value_from(x.ac)},{"b",j::value_from(x.b)},{"log_b",j::value_from(x.log_b)},
        {"width_gradient",j::value_from(x.width_gradient)},{"objective",x.objective},
        {"nuisance_amplitudes",j::value_from(x.nuisance_amplitudes)}};
}
JointState ReadState(const j::value & v)
{
    const auto & o=v.as_object();
    return {Read<std::vector<double>>(o,"ac"),Read<std::vector<double>>(o,"b"),Read<std::vector<double>>(o,"log_b"),
        Read<std::vector<double>>(o,"width_gradient"),Read<double>(o,"objective"),
        o.contains("nuisance_amplitudes") ? Read<std::vector<double>>(o,"nuisance_amplitudes") : std::vector<double>{}};
}
j::value OptionalState(const std::optional<JointState> & x) {return x ? j::value(State(*x)) : j::value(nullptr);}
const char * CandidateStatusText(const JointCandidateDiagnostic & x)
{return x.certification_status==JointCheckStatus::Passed ? "certified" : "uncertified";}
Object CandidateDiagnostic(const JointCandidateDiagnostic & x)
{
    const auto & strict=core::joint_component::StrictJointConvergenceCriteria;
    return {{"candidate_status",CandidateStatusText(x)},
        {"candidate_available",x.candidate_available},
        {"candidate_eta",j::value_from(x.candidate_eta)},
        {"candidate_ac",j::value_from(x.candidate_ac)},
        {"candidate_objective",Number(x.candidate_objective)},
        {"candidate_ac_kkt",Number(x.candidate_ac_kkt)},
        {"candidate_width_gradient_inf_norm",Number(x.candidate_width_gradient_inf_norm)},
        {"convergence_criteria_id",x.convergence_criteria_id.empty() ? "strict-v1" : x.convergence_criteria_id},
        {"global_ac_kkt_threshold",x.global_ac_kkt_threshold.value_or(strict.global_ac_kkt)},
        {"global_width_gradient_threshold",x.global_width_gradient_threshold.value_or(strict.global_width_gradient)},
        {"confirmed_log_width_change_threshold",x.confirmed_log_width_change_threshold.value_or(
            strict.confirmed_log_width_change)},
        {"candidate_log_width_change_inf_norm",Number(x.candidate_log_width_change_inf_norm)},
        {"candidate_feasible",x.candidate_feasible ? j::value(*x.candidate_feasible) : j::value(nullptr)},
        {"certification_status",StatusText(x.certification_status)},
        {"endpoint_verification_status",StatusText(x.endpoint_verification_status)},
        {"numerical_qualification_status",StatusText(x.numerical_qualification_status)},
        {"verification_reason",x.verification_reason},
        {"resource_limit_reason",x.resource_limit_reason}};
}
std::optional<JointCandidateDiagnostic> ReadCandidateDiagnostic(const j::value & value,bool reference_free_schema,
    bool convergence_criteria_schema)
{
    if(value.is_null()) return {};
    const auto & o=value.as_object(); JointCandidateDiagnostic x;
    x.candidate_available=Read<bool>(o,"candidate_available");
    x.candidate_eta=Read<std::vector<double>>(o,"candidate_eta");
    x.candidate_ac=Read<std::vector<double>>(o,"candidate_ac");
    x.candidate_objective=OptionalNumber(o.at("candidate_objective"));
    x.candidate_ac_kkt=OptionalNumber(o.at("candidate_ac_kkt"));
    x.candidate_width_gradient_inf_norm=OptionalNumber(o.at("candidate_width_gradient_inf_norm"));
    if(convergence_criteria_schema)
    {
        x.convergence_criteria_id=Read<std::string>(o,"convergence_criteria_id");
        x.global_ac_kkt_threshold=Read<double>(o,"global_ac_kkt_threshold");
        x.global_width_gradient_threshold=Read<double>(o,"global_width_gradient_threshold");
        x.confirmed_log_width_change_threshold=Read<double>(o,"confirmed_log_width_change_threshold");
        x.candidate_log_width_change_inf_norm=OptionalNumber(o.at("candidate_log_width_change_inf_norm"));
    }
    else
    {
        const auto & strict=core::joint_component::StrictJointConvergenceCriteria;
        x.convergence_criteria_id="strict-v1";
        x.global_ac_kkt_threshold=strict.global_ac_kkt;
        x.global_width_gradient_threshold=strict.global_width_gradient;
        x.confirmed_log_width_change_threshold=strict.confirmed_log_width_change;
    }
    x.candidate_feasible=o.at("candidate_feasible").is_null() ? std::nullopt :
        std::optional<bool>(Read<bool>(o,"candidate_feasible"));
    x.certification_status=Status(o.at("certification_status"));
    if(reference_free_schema)
    {
        x.endpoint_verification_status=Status(o.at("endpoint_verification_status"));
        x.numerical_qualification_status=Status(o.at("numerical_qualification_status"));
        x.verification_reason=Read<std::string>(o,"verification_reason");
    }
    x.resource_limit_reason=Read<std::string>(o,"resource_limit_reason");
    const auto candidate_status=Read<std::string>(o,"candidate_status");
    if(candidate_status!=CandidateStatusText(x))
        throw std::invalid_argument("Invalid saved joint candidate status: "+candidate_status);
    return x;
}
std::optional<JointAnalysisComponent::LegacyCandidateEvidence> ReadLegacyCandidateEvidence(const j::object & component)
{
    const j::object * source{};
    if(component.contains("legacy_candidate_evidence") && !component.at("legacy_candidate_evidence").is_null())
        source=&component.at("legacy_candidate_evidence").as_object();
    else if(component.contains("candidate_diagnostic") && !component.at("candidate_diagnostic").is_null())
        source=&component.at("candidate_diagnostic").as_object();
    if(!source) return {};
    JointAnalysisComponent::LegacyCandidateEvidence out;
    if(source->contains("reference_status")) out.reference_status=Status(source->at("reference_status"));
    else if(source->contains("legacy_reference_status")) out.reference_status=Status(source->at("legacy_reference_status"));
    if(source->contains("rank_assessment_status")) out.rank_assessment_status=Status(source->at("rank_assessment_status"));
    else if(source->contains("legacy_rank_assessment_status")) out.rank_assessment_status=Status(source->at("legacy_rank_assessment_status"));
    if(!out.reference_status && !out.rank_assessment_status) return {};
    return out;
}
Object LegacyCandidateEvidence(const std::optional<JointAnalysisComponent::LegacyCandidateEvidence> & x)
{
    Object out;
    if(!x) return out;
    if(x->reference_status) out.emplace("reference_status",StatusText(*x->reference_status));
    if(x->rank_assessment_status) out.emplace("rank_assessment_status",StatusText(*x->rank_assessment_status));
    return out;
}
j::array Checks(const std::vector<JointCheck> & values)
{
    j::array out;
    for(const auto & x:values) out.emplace_back(Object{{"name",x.name},{"scope",ScopeText(x.scope)},
        {"status",StatusText(x.status)},{"value",Number(x.value)},{"threshold",Number(x.threshold)},{"reason",x.reason}});
    return out;
}
std::vector<JointCheck> ReadChecks(const j::value & v)
{
    std::vector<JointCheck> out;
    for(const auto & item:v.as_array())
    {
        const auto & o=item.as_object();
        out.push_back({Read<std::string>(o,"name"),Status(o.at("status")),Scope(o),OptionalNumber(o.at("value")),
            OptionalNumber(o.at("threshold")),Read<std::string>(o,"reason")});
    }
    return out;
}
j::array Ranks(const std::vector<JointRankEvidence> & values)
{
    j::array out;
    for(const auto & x:values)
    {
        for(double value:x.singular_values)
            if(!std::isfinite(value)) throw std::invalid_argument("Nonfinite joint rank evidence.");
        Object rank{{"name",x.name},{"scope",ScopeText(x.scope)},
            {"rank",x.rank},{"threshold",Number(x.threshold)},{"singular_values",j::value_from(x.singular_values)}};
        if(!x.context_identity.empty()) rank.emplace("context_identity",x.context_identity);
        out.emplace_back(std::move(rank));
    }
    return out;
}
std::vector<JointRankEvidence> ReadRanks(const j::value & v,bool has_context_identity)
{
    std::vector<JointRankEvidence> out;
    for(const auto & item:v.as_array())
    {
        const auto & o=item.as_object();
        out.push_back({Read<std::string>(o,"name"),Scope(o),Read<std::size_t>(o,"rank"),Read<double>(o,"threshold"),
            Read<std::vector<double>>(o,"singular_values"),has_context_identity && o.contains("context_identity") ?
                Read<std::string>(o,"context_identity") : std::string{}});
    }
    return out;
}
Object Initialization(const JointInitialization & x)
{
    j::array atoms;
    for(const auto & a:x.atoms) atoms.emplace_back(Object{{"id",a.id},{"ols",Diagnostics(a.ols)},{"mdpde",Diagnostics(a.mdpde)},
        {"alpha",std::isfinite(a.alpha) ? j::value(a.alpha) : j::value(nullptr)},{"sample_count",a.sample_count},
        { "native_status",a.native_status ? j::value(*a.native_status) : j::value(nullptr)},{"reason",a.reason},
        {"original_b",Number(a.original_b)},{"used_b",Number(a.used_b)},{"seed_source",a.seed_source},{"donor_count",a.donor_count},
        {"first_stage_fit_qualified",a.first_stage_fit_qualified},
        {"first_stage_qualification",a.first_stage_qualification ?
            j::value(QualificationText(*a.first_stage_qualification)) : j::value(nullptr)},
        {"first_stage_recovery",a.first_stage_recovery ?
            FailedOnlyRecovery(*a.first_stage_recovery) : j::value(nullptr)}});
    return {{"valid",x.valid},{"reason",x.reason},{"b",Diagnostics(x.b)},{"atoms",std::move(atoms)},{"data_scope",x.data_scope}};
}
JointInitialization ReadInitialization(const j::value & v)
{
    const auto & o=v.as_object(); JointInitialization x;
    x.valid=Read<bool>(o,"valid"); x.reason=Read<std::string>(o,"reason"); x.b=DiagnosticValues(o.at("b"));
    x.data_scope=Read<std::string>(o,"data_scope");
    for(const auto & item:o.at("atoms").as_array())
    {
        const auto & a=item.as_object(); JointInitializationAtom atom;
        atom.id=Read<std::string>(a,"id"); atom.reason=Read<std::string>(a,"reason");
        const auto ols=DiagnosticValues(a.at("ols")),mdpde=DiagnosticValues(a.at("mdpde"));
        if(ols.size()!=3 || mdpde.size()!=3) throw std::invalid_argument("Invalid joint initialization dimensions.");
        std::copy(ols.begin(),ols.end(),atom.ols.begin()); std::copy(mdpde.begin(),mdpde.end(),atom.mdpde.begin());
        atom.alpha=OptionalNumber(a.at("alpha")).value_or(std::numeric_limits<double>::quiet_NaN());
        atom.sample_count=Read<std::size_t>(a,"sample_count");
        if(a.contains("first_stage_fit_qualified")) atom.first_stage_fit_qualified=Read<bool>(a,"first_stage_fit_qualified");
        if(a.contains("first_stage_qualification") && !a.at("first_stage_qualification").is_null())
            atom.first_stage_qualification=Qualification(a.at("first_stage_qualification"));
        if(a.contains("first_stage_recovery") && !a.at("first_stage_recovery").is_null())
            atom.first_stage_recovery=ReadFailedOnlyRecovery(a.at("first_stage_recovery"));
        if(a.contains("seed_source"))
        {
            atom.original_b=OptionalNumber(a.at("original_b")); atom.used_b=OptionalNumber(a.at("used_b"));
            atom.seed_source=Read<std::string>(a,"seed_source"); atom.donor_count=Read<std::size_t>(a,"donor_count");
        }
        if(!a.at("native_status").is_null()) atom.native_status=Read<int>(a,"native_status");
        x.atoms.push_back(std::move(atom));
    }
    return x;
}
void Require(bool condition, const char * message)
{if(!condition) throw std::invalid_argument(std::string("Invalid joint result: ")+message);}
bool IsSha256(const std::string & value)
{
    return value.size()==64 && std::all_of(value.begin(),value.end(),[](char c) {
        return (c>='0' && c<='9') || (c>='a' && c<='f');
    });
}
j::value OptionalText(const std::optional<std::string> & value)
{return value ? j::value(*value) : j::value(nullptr);}
std::optional<std::string> ReadOptionalText(const j::value & value)
{return value.is_null() ? std::nullopt : std::optional<std::string>(j::value_to<std::string>(value));}
j::value OptionalSize(std::optional<std::size_t> value)
{return value ? j::value(*value) : j::value(nullptr);}
std::optional<std::size_t> ReadOptionalSize(const j::value & value)
{return value.is_null() ? std::nullopt : std::optional<std::size_t>(j::value_to<std::size_t>(value));}
j::array BlockConditioning(const std::vector<JointBlockConditioning> & values)
{
    j::array out;
    for(const auto & x:values) out.emplace_back(Object{{"block_id",x.block_id},{"status",x.status},{"reason",x.reason},
        {"conditioning_classification",x.conditioning_classification},
        {"core_atoms",j::value_from(x.core_atoms)},{"core_atom_count",x.core_atom_count},
        {"affected_row_count",x.affected_row_count},{"free_ac_columns",x.free_ac_columns},
        {"active_amplitude_count",x.active_amplitude_count},{"free_ac_rank",OptionalSize(x.free_ac_rank)},
        {"projected_width_rank",OptionalSize(x.projected_width_rank)},{"profile_rank",OptionalSize(x.profile_rank)},
        {"active_width_bound_targets",OptionalSize(x.active_width_bound_targets)},
        {"workspace_peak_memory_estimate_bytes",OptionalSize(x.workspace_peak_memory_estimate_bytes)},
        {"free_ac_sigma_min",Number(x.free_ac_sigma_min)},{"free_ac_sigma_max",Number(x.free_ac_sigma_max)},
        {"free_ac_rank_threshold",Number(x.free_ac_rank_threshold)},
        {"free_ac_rank_threshold_margin_relative",Number(x.free_ac_rank_threshold_margin_relative)},
        {"local_ac_projected_kkt",Number(x.local_ac_projected_kkt)},
        {"local_projected_width_kkt",Number(x.local_projected_width_kkt)},
        {"free_ac_threshold_sensitive",OptionalBoolean(x.free_ac_threshold_sensitive)},
        {"free_ac_condition",Number(x.free_ac_condition)},{"projected_width_condition",Number(x.projected_width_condition)},
        {"projected_width_sigma_min",Number(x.projected_width_sigma_min)},
        {"projected_width_sigma_max",Number(x.projected_width_sigma_max)},
        {"projected_width_rank_threshold",Number(x.projected_width_rank_threshold)},
        {"projected_width_rank_threshold_margin_relative",Number(x.projected_width_rank_threshold_margin_relative)},
        {"projected_width_threshold_sensitive",OptionalBoolean(x.projected_width_threshold_sensitive)},
        {"profile_sigma_min",Number(x.profile_sigma_min)},{"profile_sigma_max",Number(x.profile_sigma_max)},
        {"profile_rank_threshold",Number(x.profile_rank_threshold)},
        {"profile_rank_threshold_margin_relative",Number(x.profile_rank_threshold_margin_relative)},
        {"profile_threshold_sensitive",OptionalBoolean(x.profile_threshold_sensitive)},
        {"profile_condition",Number(x.profile_condition)},{"free_ac_status",x.free_ac_status},
        {"free_ac_reason",x.free_ac_reason},{"projected_width_status",x.projected_width_status},
        {"projected_width_reason",x.projected_width_reason},{"free_ac_conditioning",x.free_ac_conditioning},
        {"projected_width_conditioning",x.projected_width_conditioning},{"profile_conditioning",x.profile_conditioning},
        {"profile_status",x.profile_status},
        {"profile_reason",x.profile_reason}});
    return out;
}
std::vector<JointBlockConditioning> ReadBlockConditioning(const j::value & value)
{
    std::vector<JointBlockConditioning> out;
    for(const auto & item:value.as_array())
    {
        const auto & o=item.as_object(); JointBlockConditioning x;
        x.block_id=Read<std::string>(o,"block_id"); x.status=Read<std::string>(o,"status");
        x.reason=Read<std::string>(o,"reason"); x.core_atoms=Read<std::vector<std::size_t>>(o,"core_atoms");
        if(o.contains("conditioning_classification"))
            x.conditioning_classification=Read<std::string>(o,"conditioning_classification");
        x.core_atom_count=Read<std::size_t>(o,"core_atom_count");
        x.affected_row_count=Read<std::size_t>(o,"affected_row_count");
        x.free_ac_columns=Read<std::size_t>(o,"free_ac_columns");
        x.active_amplitude_count=Read<std::size_t>(o,"active_amplitude_count");
        x.free_ac_rank=ReadOptionalSize(o.at("free_ac_rank"));
        x.projected_width_rank=ReadOptionalSize(o.at("projected_width_rank"));
        x.profile_rank=ReadOptionalSize(o.at("profile_rank"));
        if(o.contains("active_width_bound_targets"))
            x.active_width_bound_targets=ReadOptionalSize(o.at("active_width_bound_targets"));
        x.workspace_peak_memory_estimate_bytes=ReadOptionalSize(o.at("workspace_peak_memory_estimate_bytes"));
        x.free_ac_sigma_min=OptionalNumber(o.at("free_ac_sigma_min"));
        x.free_ac_sigma_max=OptionalNumber(o.at("free_ac_sigma_max"));
        x.free_ac_condition=OptionalNumber(o.at("free_ac_condition"));
        if(o.contains("free_ac_rank_threshold")) x.free_ac_rank_threshold=OptionalNumber(o.at("free_ac_rank_threshold"));
        if(o.contains("free_ac_rank_threshold_margin_relative"))
            x.free_ac_rank_threshold_margin_relative=OptionalNumber(o.at("free_ac_rank_threshold_margin_relative"));
        if(o.contains("local_ac_projected_kkt"))
            x.local_ac_projected_kkt=OptionalNumber(o.at("local_ac_projected_kkt"));
        if(o.contains("local_projected_width_kkt"))
            x.local_projected_width_kkt=OptionalNumber(o.at("local_projected_width_kkt"));
        if(o.contains("free_ac_threshold_sensitive"))
            x.free_ac_threshold_sensitive=ReadOptionalBoolean(o.at("free_ac_threshold_sensitive"));
        x.projected_width_condition=OptionalNumber(o.at("projected_width_condition"));
        if(o.contains("projected_width_sigma_min")) x.projected_width_sigma_min=OptionalNumber(o.at("projected_width_sigma_min"));
        if(o.contains("projected_width_sigma_max")) x.projected_width_sigma_max=OptionalNumber(o.at("projected_width_sigma_max"));
        if(o.contains("projected_width_rank_threshold"))
            x.projected_width_rank_threshold=OptionalNumber(o.at("projected_width_rank_threshold"));
        if(o.contains("projected_width_rank_threshold_margin_relative"))
            x.projected_width_rank_threshold_margin_relative=OptionalNumber(o.at("projected_width_rank_threshold_margin_relative"));
        if(o.contains("projected_width_threshold_sensitive"))
            x.projected_width_threshold_sensitive=ReadOptionalBoolean(o.at("projected_width_threshold_sensitive"));
        x.profile_condition=OptionalNumber(o.at("profile_condition"));
        if(o.contains("profile_sigma_min")) x.profile_sigma_min=OptionalNumber(o.at("profile_sigma_min"));
        if(o.contains("profile_sigma_max")) x.profile_sigma_max=OptionalNumber(o.at("profile_sigma_max"));
        if(o.contains("profile_rank_threshold")) x.profile_rank_threshold=OptionalNumber(o.at("profile_rank_threshold"));
        if(o.contains("profile_rank_threshold_margin_relative"))
            x.profile_rank_threshold_margin_relative=OptionalNumber(o.at("profile_rank_threshold_margin_relative"));
        if(o.contains("profile_threshold_sensitive"))
            x.profile_threshold_sensitive=ReadOptionalBoolean(o.at("profile_threshold_sensitive"));
        x.free_ac_status=Read<std::string>(o,"free_ac_status"); x.free_ac_reason=Read<std::string>(o,"free_ac_reason");
        x.projected_width_status=Read<std::string>(o,"projected_width_status");
        x.projected_width_reason=Read<std::string>(o,"projected_width_reason");
        if(o.contains("free_ac_conditioning")) x.free_ac_conditioning=Read<std::string>(o,"free_ac_conditioning");
        if(o.contains("projected_width_conditioning"))
            x.projected_width_conditioning=Read<std::string>(o,"projected_width_conditioning");
        if(o.contains("profile_conditioning")) x.profile_conditioning=Read<std::string>(o,"profile_conditioning");
        x.profile_status=Read<std::string>(o,"profile_status"); x.profile_reason=Read<std::string>(o,"profile_reason");
        out.push_back(std::move(x));
    }
    return out;
}
j::object StabilityBlock(const JointBlockParameterStability & x)
{
    j::array probes;
    for(const auto & probe:x.probes)
        probes.emplace_back(Object{{"direction",probe.direction},{"source",probe.source},{"status",probe.status},
            {"minus_seed_status",probe.minus_seed_status},{"plus_seed_status",probe.plus_seed_status},
            {"minus_seed_reason",probe.minus_seed_reason},{"plus_seed_reason",probe.plus_seed_reason},
            {"unit_direction",j::value_from(probe.unit_direction)},{"evaluations",probe.evaluations},
            {"parameter_change",Number(probe.parameter_change)},{"prediction_change",Number(probe.prediction_change)},
            {"objective_change",Number(probe.objective_change)},
            {"minus_requested_displacement_inf",Number(probe.minus_requested_displacement_inf)},
            {"minus_effective_displacement_inf",Number(probe.minus_effective_displacement_inf)},
            {"plus_requested_displacement_inf",Number(probe.plus_requested_displacement_inf)},
            {"plus_effective_displacement_inf",Number(probe.plus_effective_displacement_inf)}});
    return {{"block_id",x.block_id},{"status",x.status},{"reason",x.reason},{"probe_selected",x.probe_selected},
        {"core_atoms",j::value_from(x.core_atoms)},{"prediction_status",x.prediction_status},
        {"objective_consistency_status",x.objective_consistency_status},
        {"log_width_parameter_change_inf",Number(x.log_width_parameter_change_inf)},
        {"scaled_ac_parameter_change_inf",Number(x.scaled_ac_parameter_change_inf)},
        {"local_prediction_change_relative",Number(x.local_prediction_change_relative)},
        {"objective_difference",Number(x.objective_difference)},{"cpu_seconds",x.cpu_seconds},
        {"baseline_local_ac_kkt",Number(x.baseline_local_ac_kkt)},
        {"baseline_projected_width_kkt",Number(x.baseline_projected_width_kkt)},
        {"baseline_width_bound_slack",Number(x.baseline_width_bound_slack)},
        {"baseline_local_objective",Number(x.baseline_local_objective)},
        {"baseline_local_residual_rmse_normalized",Number(x.baseline_local_residual_rmse_normalized)},
        {"baseline_active_ac_constraints",OptionalSize(x.baseline_active_ac_constraints)},
        {"baseline_active_width_bound_targets",OptionalSize(x.baseline_active_width_bound_targets)},
        {"baseline_numerical_validity_status",x.baseline_numerical_validity_status},
        {"baseline_conditional_optimality_status",x.baseline_conditional_optimality_status},
        {"baseline_eligibility_status",x.baseline_eligibility_status},
        {"baseline_eligibility_reason",x.baseline_eligibility_reason},
        {"baseline_independent_replay_status",x.baseline_independent_replay_status},
        {"baseline_independent_replay_reason",x.baseline_independent_replay_reason},
        {"conditional_reprofile_status",x.conditional_reprofile_status},
        {"conditional_reprofile_local_ac_kkt",Number(x.conditional_reprofile_local_ac_kkt)},
        {"conditional_reprofile_projected_width_kkt",Number(x.conditional_reprofile_projected_width_kkt)},
        {"conditional_reprofile_objective_difference",Number(x.conditional_reprofile_objective_difference)},
        {"conditional_reprofile_scaled_ac_parameter_difference",
            Number(x.conditional_reprofile_scaled_ac_parameter_difference)},
        {"conditional_reprofile_prediction_change_relative",
            Number(x.conditional_reprofile_prediction_change_relative)},
        {"conditional_reprofile_independent_replay_status",x.conditional_reprofile_independent_replay_status},
        {"conditional_reprofile_independent_replay_reason",x.conditional_reprofile_independent_replay_reason},
        {"probes",std::move(probes)},
        {"perturbation_directions",x.perturbation_directions},{"perturbation_seeds",x.perturbation_seeds},
        {"perturbation_solves",x.perturbation_solves},{"perturbation_evaluations",x.perturbation_evaluations}};
}
JointBlockParameterStability ReadStabilityBlock(const j::value & value)
{
    const auto & o=value.as_object(); JointBlockParameterStability x;
    x.block_id=Read<std::string>(o,"block_id"); x.status=Read<std::string>(o,"status");
    x.reason=Read<std::string>(o,"reason"); x.core_atoms=Read<std::vector<std::size_t>>(o,"core_atoms");
    if(o.contains("probe_selected")) x.probe_selected=Read<bool>(o,"probe_selected");
    x.prediction_status=Read<std::string>(o,"prediction_status");
    x.objective_consistency_status=Read<std::string>(o,"objective_consistency_status");
    x.log_width_parameter_change_inf=OptionalNumber(o.at("log_width_parameter_change_inf"));
    x.scaled_ac_parameter_change_inf=OptionalNumber(o.at("scaled_ac_parameter_change_inf"));
    x.local_prediction_change_relative=OptionalNumber(o.at("local_prediction_change_relative"));
    x.objective_difference=OptionalNumber(o.at("objective_difference"));
    x.cpu_seconds=Read<double>(o,"cpu_seconds");
    if(o.contains("baseline_local_ac_kkt")) x.baseline_local_ac_kkt=OptionalNumber(o.at("baseline_local_ac_kkt"));
    if(o.contains("baseline_projected_width_kkt"))
        x.baseline_projected_width_kkt=OptionalNumber(o.at("baseline_projected_width_kkt"));
    if(o.contains("baseline_width_bound_slack"))
        x.baseline_width_bound_slack=OptionalNumber(o.at("baseline_width_bound_slack"));
    if(o.contains("baseline_local_objective"))
        x.baseline_local_objective=OptionalNumber(o.at("baseline_local_objective"));
    if(o.contains("baseline_local_residual_rmse_normalized"))
        x.baseline_local_residual_rmse_normalized=OptionalNumber(o.at("baseline_local_residual_rmse_normalized"));
    if(o.contains("baseline_active_ac_constraints"))
        x.baseline_active_ac_constraints=ReadOptionalSize(o.at("baseline_active_ac_constraints"));
    if(o.contains("baseline_active_width_bound_targets"))
        x.baseline_active_width_bound_targets=ReadOptionalSize(o.at("baseline_active_width_bound_targets"));
    if(o.contains("baseline_numerical_validity_status"))
        x.baseline_numerical_validity_status=Read<std::string>(o,"baseline_numerical_validity_status");
    if(o.contains("baseline_conditional_optimality_status"))
        x.baseline_conditional_optimality_status=Read<std::string>(o,"baseline_conditional_optimality_status");
    if(o.contains("baseline_eligibility_status"))
        x.baseline_eligibility_status=Read<std::string>(o,"baseline_eligibility_status");
    if(o.contains("baseline_eligibility_reason"))
        x.baseline_eligibility_reason=Read<std::string>(o,"baseline_eligibility_reason");
    if(o.contains("baseline_independent_replay_status"))
        x.baseline_independent_replay_status=Read<std::string>(o,"baseline_independent_replay_status");
    if(o.contains("baseline_independent_replay_reason"))
        x.baseline_independent_replay_reason=Read<std::string>(o,"baseline_independent_replay_reason");
    if(o.contains("conditional_reprofile_status"))
        x.conditional_reprofile_status=Read<std::string>(o,"conditional_reprofile_status");
    if(o.contains("conditional_reprofile_local_ac_kkt"))
        x.conditional_reprofile_local_ac_kkt=OptionalNumber(o.at("conditional_reprofile_local_ac_kkt"));
    if(o.contains("conditional_reprofile_projected_width_kkt"))
        x.conditional_reprofile_projected_width_kkt=OptionalNumber(o.at("conditional_reprofile_projected_width_kkt"));
    if(o.contains("conditional_reprofile_objective_difference"))
        x.conditional_reprofile_objective_difference=OptionalNumber(o.at("conditional_reprofile_objective_difference"));
    if(o.contains("conditional_reprofile_scaled_ac_parameter_difference"))
        x.conditional_reprofile_scaled_ac_parameter_difference=
            OptionalNumber(o.at("conditional_reprofile_scaled_ac_parameter_difference"));
    if(o.contains("conditional_reprofile_prediction_change_relative"))
        x.conditional_reprofile_prediction_change_relative=
            OptionalNumber(o.at("conditional_reprofile_prediction_change_relative"));
    if(o.contains("conditional_reprofile_independent_replay_status"))
        x.conditional_reprofile_independent_replay_status=
            Read<std::string>(o,"conditional_reprofile_independent_replay_status");
    if(o.contains("conditional_reprofile_independent_replay_reason"))
        x.conditional_reprofile_independent_replay_reason=
            Read<std::string>(o,"conditional_reprofile_independent_replay_reason");
    if(o.contains("probes"))
        for(const auto & probe_value:o.at("probes").as_array())
        {
            const auto & p=probe_value.as_object(); JointStabilityProbe probe;
            probe.direction=Read<std::string>(p,"direction"); probe.source=Read<std::string>(p,"source");
            probe.status=Read<std::string>(p,"status");
            if(p.contains("minus_seed_status")) probe.minus_seed_status=Read<std::string>(p,"minus_seed_status");
            if(p.contains("plus_seed_status")) probe.plus_seed_status=Read<std::string>(p,"plus_seed_status");
            if(p.contains("minus_seed_reason")) probe.minus_seed_reason=Read<std::string>(p,"minus_seed_reason");
            if(p.contains("plus_seed_reason")) probe.plus_seed_reason=Read<std::string>(p,"plus_seed_reason");
            probe.unit_direction=Read<std::vector<double>>(p,"unit_direction");
            probe.evaluations=Read<std::size_t>(p,"evaluations");
            probe.parameter_change=OptionalNumber(p.at("parameter_change"));
            probe.prediction_change=OptionalNumber(p.at("prediction_change"));
            probe.objective_change=OptionalNumber(p.at("objective_change"));
            if(p.contains("minus_requested_displacement_inf"))
                probe.minus_requested_displacement_inf=OptionalNumber(p.at("minus_requested_displacement_inf"));
            if(p.contains("minus_effective_displacement_inf"))
                probe.minus_effective_displacement_inf=OptionalNumber(p.at("minus_effective_displacement_inf"));
            if(p.contains("plus_requested_displacement_inf"))
                probe.plus_requested_displacement_inf=OptionalNumber(p.at("plus_requested_displacement_inf"));
            if(p.contains("plus_effective_displacement_inf"))
                probe.plus_effective_displacement_inf=OptionalNumber(p.at("plus_effective_displacement_inf"));
            x.probes.push_back(std::move(probe));
        }
    if(o.contains("perturbation_directions")) x.perturbation_directions=Read<std::size_t>(o,"perturbation_directions");
    if(o.contains("perturbation_seeds")) x.perturbation_seeds=Read<std::size_t>(o,"perturbation_seeds");
    x.perturbation_solves=Read<std::size_t>(o,"perturbation_solves");
    x.perturbation_evaluations=Read<std::size_t>(o,"perturbation_evaluations");
    return x;
}
j::value ParameterStability(const std::optional<JointParameterStability> & value)
{
    if(!value) return nullptr;
    const auto & x=*value; j::array blocks;
    for(const auto & block:x.blocks) blocks.emplace_back(StabilityBlock(block));
    return Object{{"recent_sweep_count",x.recent_sweep_count},{"maximum_eta_change",Number(x.maximum_eta_change)},
        {"maximum_scaled_ac_change",Number(x.maximum_scaled_ac_change)},
        {"trajectory_stability_status",x.trajectory_stability_status},{"perturbation_status",x.perturbation_status},
        {"blocks",std::move(blocks)},{"seconds",x.seconds},{"cpu_seconds",x.cpu_seconds},
        {"perturbation_directions",x.perturbation_directions},{"perturbation_seeds",x.perturbation_seeds},
        {"perturbation_solves",x.perturbation_solves},{"perturbation_evaluations",x.perturbation_evaluations},
        {"workspace_peak_memory_estimate_bytes",OptionalSize(x.workspace_peak_memory_estimate_bytes)}};
}
std::optional<JointParameterStability> ReadParameterStability(const j::value & value)
{
    if(value.is_null()) return std::nullopt;
    const auto & o=value.as_object(); JointParameterStability x;
    x.recent_sweep_count=Read<std::size_t>(o,"recent_sweep_count");
    x.maximum_eta_change=OptionalNumber(o.at("maximum_eta_change"));
    x.maximum_scaled_ac_change=OptionalNumber(o.at("maximum_scaled_ac_change"));
    x.trajectory_stability_status=Read<std::string>(o,"trajectory_stability_status");
    x.perturbation_status=Read<std::string>(o,"perturbation_status");
    for(const auto & block:o.at("blocks").as_array()) x.blocks.push_back(ReadStabilityBlock(block));
    x.seconds=Read<double>(o,"seconds"); x.cpu_seconds=Read<double>(o,"cpu_seconds");
    if(o.contains("perturbation_directions")) x.perturbation_directions=Read<std::size_t>(o,"perturbation_directions");
    if(o.contains("perturbation_seeds")) x.perturbation_seeds=Read<std::size_t>(o,"perturbation_seeds");
    x.perturbation_solves=Read<std::size_t>(o,"perturbation_solves");
    x.perturbation_evaluations=Read<std::size_t>(o,"perturbation_evaluations");
    x.workspace_peak_memory_estimate_bytes=ReadOptionalSize(o.at("workspace_peak_memory_estimate_bytes"));
    return x;
}
j::value OptionalInt(std::optional<int> value)
{return value ? j::value(*value) : j::value(nullptr);}
std::optional<int> ReadOptionalInt(const j::value & value)
{return value.is_null() ? std::nullopt : std::optional<int>(j::value_to<int>(value));}
j::array FixedHalos(const std::vector<JointFixedHaloParameter> & values)
{
    j::array out;
    for(const auto & value:values)
        out.emplace_back(Object{{"contributor_index",value.contributor_index},{"atom_id",value.atom_id},
            {"amplitude",value.amplitude},{"width",value.width},{"charge",value.charge},
            {"first_stage_sample_count",value.first_stage_sample_count},
            {"first_stage_status",OptionalInt(value.first_stage_status)},
            {"source_contract",value.source_contract},{"first_stage_fit_qualified",value.first_stage_fit_qualified},
            {"first_stage_qualification",QualificationText(FirstStageQualification(value.first_stage_status,
                value.first_stage_fit_qualified,value.first_stage_qualification,value.first_stage_recovery))},
            {"first_stage_recovery",value.first_stage_recovery ?
                FailedOnlyRecovery(*value.first_stage_recovery) : j::value(nullptr)}});
    return out;
}
std::vector<JointFixedHaloParameter> ReadFixedHalos(const j::value & value)
{
    std::vector<JointFixedHaloParameter> out;
    for(const auto & item:value.as_array())
    {
        const auto & object=item.as_object();
        JointFixedHaloParameter fixed{Read<std::size_t>(object,"contributor_index"),Read<std::string>(object,"atom_id"),
            Read<double>(object,"amplitude"),Read<double>(object,"width"),Read<double>(object,"charge"),
            Read<std::size_t>(object,"first_stage_sample_count"),ReadOptionalInt(object.at("first_stage_status")),
            Read<std::string>(object,"source_contract"),Read<bool>(object,"first_stage_fit_qualified")};
        if(object.contains("first_stage_qualification") && !object.at("first_stage_qualification").is_null())
            fixed.first_stage_qualification=Qualification(object.at("first_stage_qualification"));
        if(object.contains("first_stage_recovery") && !object.at("first_stage_recovery").is_null())
            fixed.first_stage_recovery=ReadFailedOnlyRecovery(object.at("first_stage_recovery"));
        out.push_back(std::move(fixed));
    }
    return out;
}
j::value HaloMedianSnapshot(const std::optional<JointHaloMedianSnapshot> & value)
{
    if(!value) return nullptr;
    const auto & snapshot=*value; j::array halos;
    for(const auto & halo:snapshot.halos)
        halos.emplace_back(Object{{"halo_contributor_index",halo.halo_contributor_index},
            {"atom_id",halo.atom_id},{"group_key",halo.group_key},{"amplitude",halo.amplitude},
            {"width",halo.width},{"charge",halo.charge},{"source",HaloMedianSourceText(halo.source)},
            {"donor_count",halo.donor_count},{"joint_run_id",halo.joint_run_id},{"contract",halo.contract}});
    return Object{{"contract",snapshot.contract},{"joint_run_id",snapshot.joint_run_id},
        {"global_median",j::value_from(snapshot.global_median)},
        {"valid_target_count",snapshot.valid_target_count},{"halos",std::move(halos)}};
}
std::optional<JointHaloMedianSnapshot> ReadHaloMedianSnapshot(const j::value & value)
{
    if(value.is_null()) return std::nullopt;
    const auto & object=value.as_object(); JointHaloMedianSnapshot snapshot;
    snapshot.contract=Read<std::string>(object,"contract");
    snapshot.joint_run_id=Read<std::string>(object,"joint_run_id");
    snapshot.global_median=Read<std::array<double,3>>(object,"global_median");
    snapshot.valid_target_count=Read<std::size_t>(object,"valid_target_count");
    for(const auto & value_tmp:object.at("halos").as_array())
    {
        const auto & halo=value_tmp.as_object();
        snapshot.halos.push_back({Read<std::size_t>(halo,"halo_contributor_index"),Read<std::string>(halo,"atom_id"),
            Read<GroupKey>(halo,"group_key"),Read<double>(halo,"amplitude"),Read<double>(halo,"width"),
            Read<double>(halo,"charge"),HaloMedianSourceValue(halo.at("source")),Read<std::size_t>(halo,"donor_count"),
            Read<std::string>(halo,"joint_run_id"),Read<std::string>(halo,"contract")});
    }
    return snapshot;
}
Object Units()
{
    return {{"contract","joint-kernel-map-units-v1"},{"map_value","fit-map-unit"},
        {"geometry","angstrom"},{"A","fit-map-unit*angstrom^3"},{"B","angstrom"},
        {"C","fit-map-unit*angstrom"}};
}
Object Metadata(const JointAnalysisMetadata & m)
{
    j::value normalization=nullptr,software=nullptr,solver=nullptr;
    if(m.map_normalization)
    {
        const auto & n=*m.map_normalization;
        normalization=Object{{"requested",n.requested},{"applied",n.applied},{"divisor",n.divisor}};
    }
    if(m.software)
    {
        const auto & p=*m.software;
        software=Object{{"version",p.version},{"source_sha256",p.source_sha256},
            {"configuration_sha256",p.configuration_sha256},{"build_sha256",p.build_sha256}};
    }
    if(m.solver)
    {
        const auto & p=*m.solver;
        solver=Object{{"search_method",p.search_method},
            {"contract_version",OptionalText(p.contract_version)},
            {"sparse_backend",OptionalText(p.sparse_backend)},
            {"preconditioner",OptionalText(p.preconditioner)},
            {"operator_rank_mode",OptionalText(p.operator_rank_mode)},
            {"operator_rank_backend",OptionalText(p.operator_rank_backend)},
            {"operator_pcg_iterations",OptionalInt(p.operator_pcg_iterations)},
            {"operator_damping_trials",OptionalInt(p.operator_damping_trials)},
            {"operator_rank_budget_seconds",Number(p.operator_rank_budget_seconds)},
            {"operator_rank_budget_entries",OptionalSize(p.operator_rank_budget_entries)},
            {"operator_rank_budget_workspace_bytes",OptionalSize(p.operator_rank_budget_workspace_bytes)},
            {"schwarz_core_atoms",OptionalSize(p.schwarz_core_atoms)},
            {"schwarz_overlap_hops",OptionalSize(p.schwarz_overlap_hops)},
            {"schwarz_max_block_atoms",OptionalSize(p.schwarz_max_block_atoms)},
            {"schwarz_storage_bytes",OptionalSize(p.schwarz_storage_bytes)},
            {"schwarz_scratch_bytes",OptionalSize(p.schwarz_scratch_bytes)},
            {"fixed_neighbor_core_atoms",OptionalSize(p.fixed_neighbor_core_atoms)},
            {"fixed_neighbor_local_work",OptionalText(p.fixed_neighbor_local_work)},
            {"fixed_neighbor_policy_version",OptionalText(p.fixed_neighbor_policy_version)},
            {"fixed_neighbor_maximum_sweeps",OptionalSize(p.fixed_neighbor_maximum_sweeps)},
            {"fixed_neighbor_order",OptionalText(p.fixed_neighbor_order)},
            {"fixed_neighbor_local_search",OptionalText(p.fixed_neighbor_local_search)}};
    }
    return {{"model_path",m.model_path},{"map_path",m.map_path},
        {"grid_size",j::value_from(m.grid_size)},{"grid_spacing",j::value_from(m.grid_spacing)},
        {"origin",j::value_from(m.origin)},{"simulation",m.simulation},
        {"map_normalization",std::move(normalization)},{"model_sha256",OptionalText(m.model_sha256)},
        {"map_sha256",OptionalText(m.map_sha256)},{"software",std::move(software)},
        {"solver",std::move(solver)},{"units",Units()}};
}
JointAnalysisMetadata ReadMetadata(const Object & m)
{
    JointAnalysisMetadata out;
    out.model_path=Read<std::string>(m,"model_path"); out.map_path=Read<std::string>(m,"map_path");
    out.grid_size=Read<std::array<int,3>>(m,"grid_size");
    out.grid_spacing=Read<std::array<double,3>>(m,"grid_spacing");
    out.origin=Read<std::array<double,3>>(m,"origin"); out.simulation=Read<bool>(m,"simulation");
    out.model_sha256=ReadOptionalText(m.at("model_sha256")); out.map_sha256=ReadOptionalText(m.at("map_sha256"));
    Require(m.at("units")==Units(),"unsupported unit contract");
    if(!m.at("map_normalization").is_null())
    {
        const auto & n=m.at("map_normalization").as_object();
        out.map_normalization=JointMapNormalization{Read<bool>(n,"requested"),Read<bool>(n,"applied"),Read<double>(n,"divisor")};
    }
    if(!m.at("software").is_null())
    {
        const auto & p=m.at("software").as_object();
        out.software=JointSoftwareProvenance{Read<std::string>(p,"version"),Read<std::string>(p,"source_sha256"),
            Read<std::string>(p,"configuration_sha256"),Read<std::string>(p,"build_sha256")};
    }
    if(m.contains("solver") && !m.at("solver").is_null())
    {
        const auto & p=m.at("solver").as_object();
        JointSolverProvenance solver;
        solver.search_method=Read<std::string>(p,"search_method");
        solver.fixed_neighbor_core_atoms=p.contains("fixed_neighbor_core_atoms") ?
            ReadOptionalSize(p.at("fixed_neighbor_core_atoms")) : std::nullopt;
        solver.fixed_neighbor_local_work=p.contains("fixed_neighbor_local_work") ?
            ReadOptionalText(p.at("fixed_neighbor_local_work")) : std::nullopt;
        solver.contract_version=p.contains("contract_version") ? ReadOptionalText(p.at("contract_version")) : std::nullopt;
        solver.sparse_backend=p.contains("sparse_backend") ? ReadOptionalText(p.at("sparse_backend")) : std::nullopt;
        solver.preconditioner=p.contains("preconditioner") ? ReadOptionalText(p.at("preconditioner")) : std::nullopt;
        solver.operator_rank_mode=p.contains("operator_rank_mode") ? ReadOptionalText(p.at("operator_rank_mode")) : std::nullopt;
        solver.operator_rank_backend=p.contains("operator_rank_backend") ? ReadOptionalText(p.at("operator_rank_backend")) : std::nullopt;
        solver.operator_pcg_iterations=p.contains("operator_pcg_iterations") ? ReadOptionalInt(p.at("operator_pcg_iterations")) : std::nullopt;
        solver.operator_damping_trials=p.contains("operator_damping_trials") ? ReadOptionalInt(p.at("operator_damping_trials")) : std::nullopt;
        solver.operator_rank_budget_seconds=p.contains("operator_rank_budget_seconds") ? OptionalNumber(p.at("operator_rank_budget_seconds")) : std::nullopt;
        solver.operator_rank_budget_entries=p.contains("operator_rank_budget_entries") ? ReadOptionalSize(p.at("operator_rank_budget_entries")) : std::nullopt;
        solver.operator_rank_budget_workspace_bytes=p.contains("operator_rank_budget_workspace_bytes") ? ReadOptionalSize(p.at("operator_rank_budget_workspace_bytes")) : std::nullopt;
        solver.schwarz_core_atoms=p.contains("schwarz_core_atoms") ? ReadOptionalSize(p.at("schwarz_core_atoms")) : std::nullopt;
        solver.schwarz_overlap_hops=p.contains("schwarz_overlap_hops") ? ReadOptionalSize(p.at("schwarz_overlap_hops")) : std::nullopt;
        solver.schwarz_max_block_atoms=p.contains("schwarz_max_block_atoms") ? ReadOptionalSize(p.at("schwarz_max_block_atoms")) : std::nullopt;
        solver.schwarz_storage_bytes=p.contains("schwarz_storage_bytes") ? ReadOptionalSize(p.at("schwarz_storage_bytes")) : std::nullopt;
        solver.schwarz_scratch_bytes=p.contains("schwarz_scratch_bytes") ? ReadOptionalSize(p.at("schwarz_scratch_bytes")) : std::nullopt;
        solver.fixed_neighbor_policy_version=p.contains("fixed_neighbor_policy_version") ? ReadOptionalText(p.at("fixed_neighbor_policy_version")) : std::nullopt;
        solver.fixed_neighbor_maximum_sweeps=p.contains("fixed_neighbor_maximum_sweeps") ? ReadOptionalSize(p.at("fixed_neighbor_maximum_sweeps")) : std::nullopt;
        solver.fixed_neighbor_order=p.contains("fixed_neighbor_order") ? ReadOptionalText(p.at("fixed_neighbor_order")) : std::nullopt;
        solver.fixed_neighbor_local_search=p.contains("fixed_neighbor_local_search") ? ReadOptionalText(p.at("fixed_neighbor_local_search")) : std::nullopt;
        out.solver=std::move(solver);
    }
    return out;
}
void ValidateMetadata(const JointAnalysisMetadata & m)
{
    for(const auto * hash:{&m.model_sha256,&m.map_sha256})
        Require(!*hash || IsSha256(**hash),"invalid input SHA-256");
    if(m.software)
    {
        const auto & p=*m.software;
        Require(!p.version.empty() && IsSha256(p.source_sha256) && IsSha256(p.configuration_sha256) &&
            IsSha256(p.build_sha256),"invalid software provenance");
    }
    if(m.solver)
    {
        const auto & p=*m.solver;
        Require(p.search_method=="legacy-compact" || p.search_method=="operator-pcg" ||
            p.search_method=="fixed-neighbor","invalid solver search method");
        const bool current_contract=p.contract_version &&
            *p.contract_version==JointSolverProvenanceContractVersion;
        const bool previous_contract=p.contract_version &&
            *p.contract_version==JointSolverProvenancePreviousContractVersion;
        const bool reference_contract=p.contract_version &&
            *p.contract_version==JointSolverProvenanceReferenceContractVersion;
        const bool historical_contract=!p.contract_version || reference_contract ||
            *p.contract_version==JointSolverProvenanceHistoricalContractVersion;
        Require(current_contract || previous_contract || historical_contract,"unsupported solver provenance contract");
        if(current_contract || previous_contract || reference_contract)
        {
            Require(p.search_method=="fixed-neighbor",
                "current solver provenance must use FixedNeighbor");
            Require(p.sparse_backend && *p.sparse_backend=="EIGEN",
                "current solver provenance must use EIGEN");
            Require(!p.preconditioner && !p.operator_rank_mode && !p.operator_rank_backend &&
                !p.operator_pcg_iterations && !p.operator_damping_trials &&
                !p.operator_rank_budget_seconds && !p.operator_rank_budget_entries &&
                !p.operator_rank_budget_workspace_bytes && !p.schwarz_core_atoms &&
                !p.schwarz_overlap_hops && !p.schwarz_max_block_atoms &&
                !p.schwarz_storage_bytes && !p.schwarz_scratch_bytes,
                "current solver provenance contains retired operator controls");
            const auto policy=current_contract ? FixedNeighborPolicyContractVersion :
                previous_contract ? FixedNeighborPreviousPolicyContractVersion : FixedNeighborReferencePolicyContractVersion;
            Require(p.fixed_neighbor_policy_version &&
                *p.fixed_neighbor_policy_version==policy && !p.fixed_neighbor_local_work,
                "incomplete FixedNeighbor provenance");
        }
        else if(p.sparse_backend)
            Require(*p.sparse_backend=="EIGEN" || *p.sparse_backend=="SPQR","invalid historical solver sparse backend");
        if(p.preconditioner)
            Require(*p.preconditioner=="Identity" || *p.preconditioner=="Diagonal" || *p.preconditioner=="Schwarz",
                "invalid solver preconditioner");
        if(p.operator_rank_mode)
            Require(*p.operator_rank_mode=="Auto" || *p.operator_rank_mode=="Dense" ||
                *p.operator_rank_mode=="SpqrBounds","invalid operator rank mode");
        if(p.operator_rank_backend)
            Require(*p.operator_rank_backend=="Dense" || *p.operator_rank_backend=="SpqrBounds",
                "invalid operator rank backend");
        if(p.operator_pcg_iterations) Require(*p.operator_pcg_iterations>=-1,"invalid operator PCG iteration budget");
        if(p.operator_damping_trials) Require(*p.operator_damping_trials>=0,"invalid operator damping trial budget");
        if(p.operator_rank_budget_seconds)
            Require(std::isfinite(*p.operator_rank_budget_seconds) && *p.operator_rank_budget_seconds>0,
                "invalid operator rank time budget");
        const auto require_positive=[](const auto & value) {
            if(value) Require(*value>0,"invalid solver resource budget");
        };
        require_positive(p.operator_rank_budget_entries); require_positive(p.operator_rank_budget_workspace_bytes);
        require_positive(p.schwarz_core_atoms);
        // overlap_hops is size_t, so its type already enforces >= 0; zero is valid.
        require_positive(p.schwarz_max_block_atoms); require_positive(p.schwarz_storage_bytes);
        require_positive(p.schwarz_scratch_bytes);
        if(p.schwarz_core_atoms && p.schwarz_max_block_atoms)
            Require(*p.schwarz_max_block_atoms>=*p.schwarz_core_atoms,"invalid Schwarz block limits");
        if(p.search_method=="fixed-neighbor")
        {
            Require(p.fixed_neighbor_core_atoms && *p.fixed_neighbor_core_atoms>0,
                "incomplete fixed-neighbor provenance");
            const bool current_policy=current_contract && p.fixed_neighbor_policy_version &&
                *p.fixed_neighbor_policy_version==FixedNeighborPolicyContractVersion;
            const bool previous_policy=previous_contract && p.fixed_neighbor_policy_version &&
                *p.fixed_neighbor_policy_version==FixedNeighborPreviousPolicyContractVersion;
            const bool reference_policy=reference_contract && p.fixed_neighbor_policy_version &&
                *p.fixed_neighbor_policy_version==FixedNeighborReferencePolicyContractVersion;
            const bool legacy_policy=!p.fixed_neighbor_policy_version ||
                *p.fixed_neighbor_policy_version==FixedNeighborLegacyPolicyContractVersion;
            Require(current_policy || previous_policy || reference_policy || legacy_policy,"unsupported fixed-neighbor policy contract");
            if(current_policy || previous_policy || reference_policy)
            {
                if(p.fixed_neighbor_local_work)
                    Require(*p.fixed_neighbor_local_work=="one-accepted",
                        "current FixedNeighbor cannot claim retired local work");
            }
            else
            {
                Require(p.fixed_neighbor_local_work && !p.fixed_neighbor_local_work->empty(),
                    "incomplete legacy fixed-neighbor provenance");
                Require(*p.fixed_neighbor_local_work=="full" || *p.fixed_neighbor_local_work=="one-accepted" ||
                    *p.fixed_neighbor_local_work=="two-accepted","invalid legacy fixed-neighbor local work");
            }
            if(p.fixed_neighbor_maximum_sweeps) Require(*p.fixed_neighbor_maximum_sweeps>0,
                "invalid fixed-neighbor sweep budget");
            if(p.fixed_neighbor_order)
                Require(*p.fixed_neighbor_order=="forward" || *p.fixed_neighbor_order=="reverse",
                    "invalid fixed-neighbor block order");
            if(p.fixed_neighbor_local_search)
                Require(*p.fixed_neighbor_local_search==FixedNeighborLocalSearchName ||
                    *p.fixed_neighbor_local_search=="legacy-compact",
                    "invalid fixed-neighbor local search");
        }
    }
    if(m.map_normalization)
    {
        const auto & n=*m.map_normalization;
        Require(std::isfinite(n.divisor) && n.divisor>0,"invalid normalization divisor");
        Require(n.applied ? n.requested && !m.simulation : n.divisor==1,"inconsistent normalization");
    }
}
void ValidateState(const JointState & x,std::size_t atoms,std::size_t groups=0)
{
    Require(x.ac.size()==2*atoms && x.b.size()==atoms && x.log_b.size()==atoms && x.width_gradient.size()==atoms,"state dimensions");
    for(const auto * values:{&x.ac,&x.b,&x.log_b,&x.width_gradient})
        for(double value:*values) Require(std::isfinite(value),"nonfinite state");
    for(double b:x.b) Require(b>0,"nonpositive width");
    Require(std::isfinite(x.objective),"nonfinite objective");
    Require(x.nuisance_amplitudes.size()==groups,"nuisance amplitude dimensions");
    for(double value:x.nuisance_amplitudes) Require(std::isfinite(value),"nonfinite nuisance amplitude");
}
void ValidateCandidateDiagnostic(const JointCandidateDiagnostic & x,std::size_t atoms,bool trusted_state,
    bool reference_free_contract)
{
    const bool no_criteria=x.convergence_criteria_id.empty() && !x.global_ac_kkt_threshold &&
        !x.global_width_gradient_threshold && !x.confirmed_log_width_change_threshold;
    const bool has_criteria=!x.convergence_criteria_id.empty() && x.global_ac_kkt_threshold &&
        x.global_width_gradient_threshold && x.confirmed_log_width_change_threshold;
    Require(no_criteria || has_criteria,"incomplete convergence criteria provenance");
    if(has_criteria)
    {
        const auto & strict=core::joint_component::StrictJointConvergenceCriteria;
        const auto & practical=core::joint_component::PracticalJointConvergenceCriteria;
        Require(x.convergence_criteria_id=="strict-v1" || x.convergence_criteria_id=="practical-v1" ||
            x.convergence_criteria_id=="custom-v1","invalid convergence criteria identifier");
        Require(std::isfinite(*x.global_ac_kkt_threshold) && *x.global_ac_kkt_threshold>0 &&
            std::isfinite(*x.global_width_gradient_threshold) && *x.global_width_gradient_threshold>0 &&
            std::isfinite(*x.confirmed_log_width_change_threshold) && *x.confirmed_log_width_change_threshold>0,
            "invalid convergence criteria threshold");
        if(x.convergence_criteria_id=="strict-v1")
            Require(*x.global_ac_kkt_threshold==strict.global_ac_kkt &&
                *x.global_width_gradient_threshold==strict.global_width_gradient &&
                *x.confirmed_log_width_change_threshold==strict.confirmed_log_width_change,
                "strict convergence criteria do not match their identifier");
        if(x.convergence_criteria_id=="practical-v1")
            Require(*x.global_ac_kkt_threshold==practical.global_ac_kkt &&
                *x.global_width_gradient_threshold==practical.global_width_gradient &&
                *x.confirmed_log_width_change_threshold==practical.confirmed_log_width_change,
                "practical convergence criteria do not match their identifier");
    }
    Require(!x.candidate_log_width_change_inf_norm ||
        (std::isfinite(*x.candidate_log_width_change_inf_norm) && *x.candidate_log_width_change_inf_norm>=0),
        "invalid candidate log-width change");
    Require(atoms<=std::numeric_limits<std::size_t>::max()/2,"candidate dimensions");
    if(x.candidate_available)
    {
        Require(x.candidate_eta.size()==atoms && x.candidate_ac.size()==2*atoms && x.candidate_objective,
            "candidate dimensions");
        for(double eta:x.candidate_eta)
            Require(std::isfinite(eta) && std::isfinite(std::exp(eta)) && std::exp(eta)>0,"invalid candidate width");
        for(double value:x.candidate_ac) Require(std::isfinite(value),"nonfinite candidate A/C");
        Require(std::isfinite(*x.candidate_objective),"nonfinite candidate objective");
    }
    else Require(x.candidate_eta.empty() && x.candidate_ac.empty() && !x.candidate_objective,
        "unavailable candidate contains parameters");
    for(const auto * value:{&x.candidate_ac_kkt,&x.candidate_width_gradient_inf_norm})
        Require(!*value || (std::isfinite(**value) && **value>=0),"invalid candidate diagnostic value");
    Require((x.certification_status==JointCheckStatus::Passed)==trusted_state,
        "candidate certification disagrees with trusted state");
    Require(x.certification_status!=JointCheckStatus::Passed || x.candidate_available,
        "certified candidate is unavailable");
    if(reference_free_contract)
    {
        Require((x.numerical_qualification_status==JointCheckStatus::Passed)==trusted_state,
            "candidate numerical qualification disagrees with trusted state");
        Require(x.numerical_qualification_status!=JointCheckStatus::Passed || x.candidate_available,
            "qualified candidate is unavailable");
        Require(x.endpoint_verification_status!=JointCheckStatus::Passed || x.candidate_available,
            "verified candidate is unavailable");
    }
}
void ValidateBlockConditioning(const JointBlockConditioning & x)
{
    Require(!x.block_id.empty() && !x.status.empty() && !x.free_ac_status.empty() &&
        !x.projected_width_status.empty() && !x.profile_status.empty(),"local conditioning identity/status");
    const auto valid_conditioning_class=[](const std::string & value) {
        return value=="WellConditioned" || value=="IllConditioned" || value=="Singular" ||
            value=="ActiveFace" || value=="Unavailable" || value=="NotAssessed";
    };
    Require(valid_conditioning_class(x.conditioning_classification) &&
        valid_conditioning_class(x.free_ac_conditioning) &&
        valid_conditioning_class(x.projected_width_conditioning) &&
        valid_conditioning_class(x.profile_conditioning),"local conditioning classification");
    Require(!x.free_ac_rank || *x.free_ac_rank<=x.free_ac_columns,"local free-AC rank");
    Require(!x.projected_width_rank || *x.projected_width_rank<=x.core_atom_count,"local projected-width rank");
    Require(!x.profile_rank || *x.profile_rank<=x.core_atom_count,"local profile rank");
    Require(!x.active_width_bound_targets || *x.active_width_bound_targets<=x.core_atom_count,
        "local active width bound count");
    for(const auto * value:{&x.free_ac_sigma_min,&x.free_ac_sigma_max,&x.free_ac_condition,
        &x.free_ac_rank_threshold,&x.free_ac_rank_threshold_margin_relative,
        &x.local_ac_projected_kkt,&x.local_projected_width_kkt,
        &x.projected_width_sigma_min,&x.projected_width_sigma_max,&x.projected_width_condition,
        &x.projected_width_rank_threshold,&x.projected_width_rank_threshold_margin_relative,
        &x.profile_sigma_min,&x.profile_sigma_max,&x.profile_condition,
        &x.profile_rank_threshold,&x.profile_rank_threshold_margin_relative})
        Require(!*value || (std::isfinite(**value) && **value>=0),"local conditioning value");
}
void ValidateStabilityBlock(const JointBlockParameterStability & x)
{
    Require(!x.block_id.empty() && !x.status.empty() && !x.prediction_status.empty() &&
        !x.objective_consistency_status.empty(),"local stability identity/status");
    Require(std::isfinite(x.cpu_seconds) && x.cpu_seconds>=0,"local stability CPU time");
    for(const auto * value:{&x.log_width_parameter_change_inf,&x.scaled_ac_parameter_change_inf,
        &x.local_prediction_change_relative,&x.objective_difference,&x.baseline_local_ac_kkt,
        &x.baseline_projected_width_kkt,&x.baseline_local_objective,&x.baseline_local_residual_rmse_normalized,
        &x.conditional_reprofile_local_ac_kkt,&x.conditional_reprofile_projected_width_kkt,
        &x.conditional_reprofile_objective_difference,&x.conditional_reprofile_scaled_ac_parameter_difference,
        &x.conditional_reprofile_prediction_change_relative})
        Require(!*value || (std::isfinite(**value) && **value>=0),"local stability value");
    Require(!x.baseline_width_bound_slack || std::isfinite(*x.baseline_width_bound_slack),
        "local stability bound slack");
    for(const auto & probe:x.probes)
    {
        Require(!probe.direction.empty() && !probe.source.empty() &&
            (probe.status=="Assessed" || probe.status=="OneSidedAssessed" || probe.status=="SolverFailed" ||
                probe.status=="BudgetExhausted" || probe.status=="Unavailable"),"local stability probe status");
        double norm{};
        for(const double value:probe.unit_direction)
        {Require(std::isfinite(value),"local stability probe direction"); norm=std::hypot(norm,value);}
        Require(probe.unit_direction.empty() || std::abs(norm-1.)<=1e-10,"local stability probe unit direction");
        for(const auto * value:{&probe.parameter_change,&probe.prediction_change,&probe.objective_change})
            Require(!*value || (std::isfinite(**value) && **value>=0),"local stability probe value");
        for(const auto * value:{&probe.minus_requested_displacement_inf,&probe.minus_effective_displacement_inf,
            &probe.plus_requested_displacement_inf,&probe.plus_effective_displacement_inf})
            Require(!*value || (std::isfinite(**value) && **value>=0),"local stability probe displacement");
        Require(!probe.minus_seed_status.empty() && !probe.plus_seed_status.empty(),
            "local stability probe seed status");
    }
}
void ValidateParameterStability(const JointParameterStability & x)
{
    Require(!x.trajectory_stability_status.empty() && !x.perturbation_status.empty(),"parameter stability status");
    Require(std::isfinite(x.seconds) && x.seconds>=0 && std::isfinite(x.cpu_seconds) && x.cpu_seconds>=0,
        "parameter stability time");
    for(const auto * value:{&x.maximum_eta_change,&x.maximum_scaled_ac_change})
        Require(!*value || (std::isfinite(**value) && **value>=0),"parameter stability trajectory value");
    for(const auto & block:x.blocks) ValidateStabilityBlock(block);
}
j::value SelectionDomain(const std::optional<JointSelectionDomain> & domain)
{
    if(!domain) return nullptr;
    j::value halo_eligibility=nullptr,fixed_halo_selection=nullptr,pure_halo_support_exclusion=nullptr;
    if(domain->halo_eligibility)
    {
        const auto & provenance=*domain->halo_eligibility; j::array excluded;
        for(const auto & record:provenance.excluded_halos)
            excluded.emplace_back(Object{{"atom_id",record.atom_id},{"support_rows",record.support_rows},
                {"distinct_radial_distances",record.distinct_radial_distances},{"reason",record.reason}});
        halo_eligibility=Object{{"minimum_observations",provenance.minimum_observations},
            {"minimum_distinct_radial_information",provenance.minimum_distinct_radial_information},
            {"candidate_halos",provenance.candidate_halos},{"retained_halos",provenance.retained_halos},
            {"excluded_halos",std::move(excluded)}};
    }
    if(domain->fixed_halo_selection)
    {
        const auto & provenance=*domain->fixed_halo_selection;
        fixed_halo_selection=Object{{"minimum_observations",provenance.minimum_observations},
            {"candidate_halos",provenance.candidate_halos},{"retained_halos",provenance.retained_halos},
            {"zero_support_halos",provenance.zero_support_halos}};
    }
    if(domain->pure_halo_support_exclusion)
    {
        const auto & provenance=*domain->pure_halo_support_exclusion;
        pure_halo_support_exclusion=Object{{"original_rows",provenance.original_rows},
            {"fit_rows",provenance.fit_rows},{"excluded_rows",provenance.excluded_rows},
            {"candidate_halos",provenance.candidate_halos},{"retained_halos",provenance.retained_halos},
            {"zero_support_halos",provenance.zero_support_halos}};
    }
    Object out{{"contract",domain->contract},{"target_indices",j::value_from(domain->target_indices)},
        {"observation_radius",domain->observation_radius},{"support_radius",domain->support_radius},
        {"contributor_policy",domain->contributor_policy},{"halo_eligibility",std::move(halo_eligibility)},
        {"fixed_halo_selection",std::move(fixed_halo_selection)}};
    if(domain->pure_halo_support_exclusion)
        out.emplace("pure_halo_support_exclusion",std::move(pure_halo_support_exclusion));
    return out;
}
std::optional<JointSelectionDomain> ReadSelectionDomain(const j::value & value)
{
    if(value.is_null()) return std::nullopt;
    const auto & o=value.as_object();
    JointSelectionDomain out;
    out.contract=Read<std::string>(o,"contract");
    out.target_indices=Read<std::vector<std::size_t>>(o,"target_indices");
    out.observation_radius=Read<double>(o,"observation_radius");
    out.support_radius=Read<double>(o,"support_radius");
    out.contributor_policy=Read<std::string>(o,"contributor_policy");
    if(o.contains("halo_eligibility") && !o.at("halo_eligibility").is_null())
    {
        const auto & item=o.at("halo_eligibility").as_object(); JointHaloEligibilityProvenance provenance;
        provenance.minimum_observations=Read<std::size_t>(item,"minimum_observations");
        provenance.minimum_distinct_radial_information=Read<std::size_t>(item,"minimum_distinct_radial_information");
        provenance.candidate_halos=Read<std::size_t>(item,"candidate_halos");
        provenance.retained_halos=Read<std::size_t>(item,"retained_halos");
        for(const auto & excluded_value:item.at("excluded_halos").as_array())
        {
            const auto & excluded=excluded_value.as_object();
            provenance.excluded_halos.push_back({Read<std::string>(excluded,"atom_id"),
                Read<std::size_t>(excluded,"support_rows"),Read<std::size_t>(excluded,"distinct_radial_distances"),
                Read<std::string>(excluded,"reason")});
        }
        out.halo_eligibility=std::move(provenance);
    }
    if(o.contains("fixed_halo_selection") && !o.at("fixed_halo_selection").is_null())
    {
        const auto & item=o.at("fixed_halo_selection").as_object();
        out.fixed_halo_selection=JointFixedHaloSelectionProvenance{
            Read<std::size_t>(item,"minimum_observations"),Read<std::size_t>(item,"candidate_halos"),
            Read<std::size_t>(item,"retained_halos"),Read<std::size_t>(item,"zero_support_halos")};
    }
    if(o.contains("pure_halo_support_exclusion") && !o.at("pure_halo_support_exclusion").is_null())
    {
        const auto & item=o.at("pure_halo_support_exclusion").as_object();
        out.pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
            Read<std::size_t>(item,"original_rows"),Read<std::size_t>(item,"fit_rows"),
            Read<std::size_t>(item,"excluded_rows"),Read<std::size_t>(item,"candidate_halos"),
            Read<std::size_t>(item,"retained_halos"),Read<std::size_t>(item,"zero_support_halos")};
    }
    return out;
}
j::value Layout(const std::optional<JointParameterLayout> & layout)
{
    if(!layout) return nullptr;
    j::array groups;
    for(const auto & group:layout->groups) groups.emplace_back(Object{{"row",group.row},{"atoms",j::value_from(group.atoms)}});
    return Object{{"full_atoms",j::value_from(layout->full_atoms)},{"informative_rows",j::value_from(layout->informative_rows)},{"groups",std::move(groups)}};
}
std::optional<JointParameterLayout> ReadLayout(const Object & object)
{
    if(!object.contains("layout") || object.at("layout").is_null()) return std::nullopt;
    const auto & o=object.at("layout").as_object();
    JointParameterLayout out{Read<std::vector<std::size_t>>(o,"full_atoms"),Read<std::vector<std::size_t>>(o,"informative_rows"),{}};
    for(const auto & item:o.at("groups").as_array())
    {const auto & group=item.as_object(); out.groups.push_back({Read<std::size_t>(group,"row"),Read<std::vector<std::size_t>>(group,"atoms")});}
    return out;
}
void ValidateLayout(const JointParameterLayout & layout,const std::vector<std::size_t> & atoms,
    const std::vector<std::size_t> & rows,const JointSelectionDomain * selection)
{
    std::set<std::size_t> remaining_atoms(atoms.begin(),atoms.end()),remaining_rows(rows.begin(),rows.end());
    Require(std::is_sorted(layout.full_atoms.begin(),layout.full_atoms.end()),"unsorted full atom layout");
    Require(std::is_sorted(layout.informative_rows.begin(),layout.informative_rows.end()),"unsorted informative rows");
    for(auto a:layout.full_atoms) Require(remaining_atoms.erase(a)==1,"full atom mapping");
    for(auto r:layout.informative_rows) Require(remaining_rows.erase(r)==1,"informative row mapping");
    std::optional<std::size_t> previous;
    for(const auto & group:layout.groups)
    {
        Require(!previous || *previous<group.row,"unsorted or duplicate nuisance row"); previous=group.row;
        Require(!group.atoms.empty() && remaining_rows.erase(group.row)==1,"nuisance row mapping");
        Require(std::is_sorted(group.atoms.begin(),group.atoms.end()),"unsorted nuisance atoms");
        for(auto a:group.atoms)
        {
            Require(remaining_atoms.erase(a)==1,"nuisance atom mapping");
            Require(selection && !std::binary_search(selection->target_indices.begin(),selection->target_indices.end(),a),"nuisance atom must be a recorded halo");
        }
    }
    Require(remaining_atoms.empty() && remaining_rows.empty(),"incomplete parameter layout");
}
j::value Target(const std::optional<JointTargetEvidence> & value)
{
    if(!value) return nullptr;
    const auto & t=*value; j::array atoms,point_atoms;
    for(const auto & a:t.atoms) atoms.emplace_back(Object{{"atom",a.atom},{"status",StatusText(a.status)},
        {"null_space_leakage",Number(a.null_space_leakage)},{"reason",a.reason}});
    for(const auto & a:t.point_atoms)
        point_atoms.emplace_back(Object{{"atom",a.atom},{"component_convergence",StatusText(a.component_convergence)},
            {"parameter_validity",StatusText(a.parameter_validity)},{"local_conditioning",a.local_conditioning},
            {"parameter_stability",a.parameter_stability},{"point_quality",a.point_quality},
            {"local_conditioning_blocks",j::value_from(a.local_conditioning_blocks)},
            {"parameter_stability_blocks",j::value_from(a.parameter_stability_blocks)},{"reason",a.reason}});
    return Object{{"contract",t.contract},{"checks",Checks(t.checks)},{"ranks",Ranks(t.ranks)},
        {"atoms",atoms},{"point_atoms",point_atoms},{"column_scales",j::value_from(t.column_scales)},
        {"original_rows",t.original_rows}};
}
std::optional<JointTargetEvidence> ReadTarget(const j::value & value,int schema_version)
{
    if(value.is_null()) return {};
    const auto & o=value.as_object(); JointTargetEvidence t;
    t.contract=Read<std::string>(o,"contract"); t.checks=ReadChecks(o.at("checks"));
    t.ranks=ReadRanks(o.at("ranks"),schema_version>=9);
    t.column_scales=Read<std::vector<double>>(o,"column_scales"); t.original_rows=Read<std::size_t>(o,"original_rows");
    for(const auto & value_tmp:o.at("atoms").as_array())
    {
        const auto & a=value_tmp.as_object();
        t.atoms.push_back({Read<std::size_t>(a,"atom"),Status(a.at("status")),OptionalNumber(a.at("null_space_leakage")),Read<std::string>(a,"reason")});
    }
    if(o.contains("point_atoms")) for(const auto & value_tmp:o.at("point_atoms").as_array())
    {
        const auto & a=value_tmp.as_object(); JointAtomPointValidation point;
        point.atom=Read<std::size_t>(a,"atom");
        point.component_convergence=Status(a.at("component_convergence"));
        point.parameter_validity=Status(a.at("parameter_validity"));
        point.local_conditioning=Read<std::string>(a,"local_conditioning");
        point.parameter_stability=Read<std::string>(a,"parameter_stability");
        point.point_quality=Read<std::string>(a,"point_quality");
        point.local_conditioning_blocks=Read<std::vector<std::string>>(a,"local_conditioning_blocks");
        point.parameter_stability_blocks=Read<std::vector<std::string>>(a,"parameter_stability_blocks");
        point.reason=Read<std::string>(a,"reason"); t.point_atoms.push_back(std::move(point));
    }
    return t;
}
void ValidateTarget(const std::optional<JointTargetEvidence> & value,const std::vector<std::size_t> & atoms,
    const std::vector<std::size_t> & full_atoms,const std::optional<JointSelectionDomain> & selection,
    std::size_t full,std::size_t rows)
{
    if(!value) return;
    const auto & t=*value;
    const bool derived=t.contract=="target-estimability-v3" || t.contract=="target-estimability-v4";
    const bool point_contract=t.contract==JointTargetPointValidationContractVersion;
    Require(selection.has_value() && (t.contract=="target-estimability-v1" || t.contract=="target-estimability-v2" || derived || point_contract),
        "target evidence contract");
    Require(t.original_rows==rows,"target rank row context");
    if(point_contract)
    {
        Require(t.ranks.empty() && t.atoms.empty() && t.column_scales.empty(),"point evidence cannot contain rank evidence");
        std::set<std::size_t> expected(full_atoms.begin(),full_atoms.end()),seen;
        for(const auto & point:t.point_atoms)
        {
            Require(expected.contains(point.atom) && seen.insert(point.atom).second,"point atom mapping");
            Require(point.local_conditioning=="Adequate" || point.local_conditioning=="Warning" ||
                point.local_conditioning=="Unavailable" || point.local_conditioning=="NotAssessed","point local conditioning status");
            Require(point.parameter_stability=="Stable" || point.parameter_stability=="Sensitive" ||
                point.parameter_stability=="Unstable" || point.parameter_stability=="Unavailable" ||
                point.parameter_stability=="NotAssessed","point parameter stability status");
            Require(point.point_quality=="Available" || point.point_quality=="Warning" ||
                point.point_quality=="Unavailable" || point.point_quality=="NotAssessed","point quality status");
            Require(!point.reason.empty(),"point quality reason");
        }
        Require(seen==expected,"incomplete point atom evidence");
    }
    if(derived)
    {
        Require(t.column_scales.empty(),"derived target evidence cannot contain column scales");
        Require(std::none_of(t.atoms.begin(),t.atoms.end(),[](const auto & a){return a.null_space_leakage.has_value();}),
            "derived target evidence cannot contain null-space leakage");
        const bool derived_pass=std::any_of(t.atoms.begin(),t.atoms.end(),[](const auto & atom) {
            return atom.status==JointCheckStatus::Passed &&
                atom.reason=="derived-from-component-identifiability";
        });
        if(derived_pass)
            Require(t.ranks.size()==2,"derived target rank evidence dimensions");
        if(!t.ranks.empty())
        {
            Require(t.ranks.size()==2,"derived target rank evidence dimensions");
            const auto design=std::find_if(t.ranks.begin(),t.ranks.end(),[](const auto & rank){return rank.name=="design";});
            const auto widths=std::find_if(t.ranks.begin(),t.ranks.end(),[](const auto & rank){return rank.name=="projected-width";});
            Require(design!=t.ranks.end() && widths!=t.ranks.end(),"derived target rank evidence sources");
            if(derived_pass)
            {
                Require(design->scope==widths->scope && !design->context_identity.empty() &&
                    design->context_identity==widths->context_identity,
                    "derived target rank context");
                Require(full<=std::numeric_limits<std::size_t>::max()/2 && design->rank==2*full &&
                    widths->rank==full && design->singular_values.size()==2*full &&
                    widths->singular_values.size()==full && std::isfinite(design->threshold) &&
                    design->threshold>0 && std::isfinite(widths->threshold) && widths->threshold>0,
                    "derived target rank dimensions");
            }
        }
    }
    Require(t.column_scales.empty() || t.column_scales.size()==3*full,"target scale dimensions");
    for(double s:t.column_scales) Require(std::isfinite(s) && s>0,"target column scale");
    if(!point_contract)
    {
        std::set<std::size_t> expected(atoms.begin(),atoms.end()),seen;
        for(const auto & a:t.atoms)
        {
            Require(expected.contains(a.atom) && seen.insert(a.atom).second,"target atom mapping");
            if(a.null_space_leakage) Require(std::isfinite(*a.null_space_leakage) && *a.null_space_leakage>=0,"target null space leakage");
        }
        Require(seen==expected,"incomplete target atom evidence");
    }
}
void Validate(const JointAnalysisResult & x,int schema_version)
{
    ValidateMetadata(x.metadata);
    if(schema_version<9)
    {
        const auto old_target_contract=[](const auto & evidence) {
            return !evidence || (evidence->contract!="target-estimability-v3" &&
                evidence->contract!="target-estimability-v4");
        };
        Require(old_target_contract(x.target_evidence) &&
            std::all_of(x.components.begin(),x.components.end(),[&](const auto & component) {
                return old_target_contract(component.target_evidence);
            }),"derived target evidence requires schema 9");
    }
    if(schema_version<10)
    {
        const auto new_target_contract=[](const auto & evidence) {
            return !evidence || evidence->contract!="target-estimability-v4";
        };
        Require(new_target_contract(x.target_evidence) &&
            std::all_of(x.components.begin(),x.components.end(),[&](const auto & component) {
                return new_target_contract(component.target_evidence);
            }),"target-estimability-v4 requires schema 10");
    }
    if(schema_version<8)
    {
        Require(!x.metadata.solver || !x.metadata.solver->contract_version ||
            *x.metadata.solver->contract_version!=JointSolverProvenanceContractVersion,
            "correction-free solver contract requires schema 8");
        const auto old_target_contract=[](const auto & evidence) {
            return !evidence || evidence->contract=="target-estimability-v1";
        };
        Require(old_target_contract(x.target_evidence) &&
            std::all_of(x.components.begin(),x.components.end(),[&](const auto & component) {
                return old_target_contract(component.target_evidence);
            }),"target-estimability-v2 requires schema 8");
    }
    const bool reference_free_contract=x.metadata.solver && x.metadata.solver->contract_version &&
        (*x.metadata.solver->contract_version==JointSolverProvenanceContractVersion ||
         *x.metadata.solver->contract_version==JointSolverProvenancePreviousContractVersion);
    const auto atoms=x.atom_ids.size(), rows=x.row_ids.size();
    Require(!x.selection_domain || x.selection_domain->IsValid(atoms),"selection domain");
    const bool pruned_policy=x.selection_domain &&
        x.selection_domain->contributor_policy==FullABCHaloContributorPolicy;
    const bool fixed_halo_policy=x.selection_domain &&
        x.selection_domain->contributor_policy==FixedHaloContributorPolicy;
    const bool pure_exclusion_policy=x.selection_domain &&
        x.selection_domain->contributor_policy==PureHaloSupportExclusionContributorPolicy;
    const bool target_only_fixed_halo=x.parameterization_contract==JointTargetOnlyFixedHaloParameterizationContract;
    const bool target_only_pure_exclusion=x.parameterization_contract==
        JointPureHaloSupportExclusionParameterizationContract;
    Require(schema_version>=11 || !pruned_policy,"halo pruning policy requires schema 11");
    Require(schema_version>=13 || !fixed_halo_policy,"fixed-halo policy requires schema 13");
    Require(schema_version>=17 || !pure_exclusion_policy,"Pure Exclusion policy requires schema 17");
    Require(schema_version>=17 || !target_only_pure_exclusion,"Pure Exclusion parameterization requires schema 17");
    Require(schema_version>=14 || (!target_only_fixed_halo && x.fixed_halos.empty()),
        "target-only fixed-halo state requires schema 14");
    if(pruned_policy)
    {
        Require(x.selection_domain->halo_eligibility.has_value(),"halo pruning requires eligibility provenance");
        Require(!x.layout && x.parameterization_contract=="full-abc-v1",
            "halo pruning results require FullABC parameterization");
    }
    Require(x.initialization.data_scope=="caller-provided-widths" ||
        x.initialization.data_scope=="contributor-local-sampling-may-read-outside-target-domain","initialization data scope");
    std::set<std::string> initialized;
    for(const auto & atom:x.initialization.atoms)
        Require(std::find(x.atom_ids.begin(),x.atom_ids.end(),atom.id)!=x.atom_ids.end() && initialized.insert(atom.id).second,
            "initialization atom identity");
    const std::set<std::string> atom_id_set(x.atom_ids.begin(),x.atom_ids.end());
    Require(atoms>0 && atom_id_set.size()==atoms,"atom identities");
    if(pruned_policy) for(const auto & excluded:x.selection_domain->halo_eligibility->excluded_halos)
        Require(!atom_id_set.contains(excluded.atom_id),"excluded halo is a retained contributor");
    Require(std::set<std::string>(x.row_ids.begin(),x.row_ids.end()).size()==rows,"row identities");
    Require(x.available_row_mask.size()==rows,"row mask dimensions");
    Require(std::isfinite(x.observation_scale) && x.observation_scale>=1,"observation scale");
    Require(!x.objective || std::isfinite(*x.objective),"nonfinite objective");
    if(target_only_fixed_halo)
        Require(fixed_halo_policy && !x.layout,"target-only fixed-halo parameterization contract");
    else if(target_only_pure_exclusion)
        Require(pure_exclusion_policy && !x.layout && x.fixed_halos.empty(),
            "target-only Pure Exclusion parameterization contract");
    else
    {
        Require(x.parameterization_contract==(x.layout ? "singleton-halo-profile-v1" : "full-abc-v1"),"parameterization contract");
        Require(x.fixed_halos.empty(),"fixed-halo parameters require target-only parameterization");
    }
    if(target_only_fixed_halo)
    {
        const auto & targets=x.selection_domain->target_indices;
        std::vector<std::size_t> expected_halos;
        for(std::size_t atom=0;atom<atoms;++atom)
            if(!std::binary_search(targets.begin(),targets.end(),atom)) expected_halos.push_back(atom);
        Require(x.fixed_halos.empty() || x.fixed_halos.size()==expected_halos.size(),"incomplete fixed-halo state");
        if(!x.fixed_halos.empty())
            for(std::size_t k=0;k<x.fixed_halos.size();++k)
            {
                const auto & fixed=x.fixed_halos[k];
                Require(fixed.contributor_index==expected_halos[k] && fixed.atom_id==x.atom_ids[expected_halos[k]],
                    "fixed-halo identity mapping");
                Require(std::isfinite(fixed.amplitude) && fixed.amplitude>=0 &&
                    std::isfinite(fixed.width) && fixed.width>0 && fixed.charge==0.0,
                    "invalid fixed-halo parameters");
                Require(IsQualifiedFirstStageFixedHalo(fixed) &&
                    fixed.source_contract==FixedHaloParameterSourceContractVersion,
                    "invalid fixed-halo FirstStage provenance");
            }
        Require(x.fixed_halos.size()==expected_halos.size() ||
            (!x.assembled_state && !x.objective && x.components.empty() &&
                x.runtime_convergence!=JointCheckStatus::Passed),
            "fitted target-only result requires every fixed halo");
    }
    if(target_only_pure_exclusion)
    {
        const auto & targets=x.selection_domain->target_indices;
        const auto & provenance=*x.selection_domain->pure_halo_support_exclusion;
        Require(provenance.original_rows==rows,"Pure Exclusion original row provenance");
        Require(x.fitting_row_ids.size()==x.fitting_to_original_rows.size() &&
            x.fitting_row_ids.size()==provenance.fit_rows && !x.fitting_row_ids.empty(),
            "Pure Exclusion fitting identity dimensions");
        std::set<std::size_t> fitting_rows;
        for(std::size_t k=0;k<x.fitting_row_ids.size();++k)
        {
            const auto row=x.fitting_to_original_rows[k];
            Require(row<rows && fitting_rows.insert(row).second && x.fitting_row_ids[k]==x.row_ids[row],
                "Pure Exclusion fitting row identity");
        }
        Require(x.fitting_to_original_contributors.size()==targets.size(),
            "Pure Exclusion fitting contributor dimensions");
        for(std::size_t k=0;k<targets.size();++k)
            Require(x.fitting_to_original_contributors[k]<atoms &&
                x.fitting_to_original_contributors[k]==targets[k],
                "Pure Exclusion fitting contributor mapping");
        const bool has_halos=atoms>targets.size();
        if(has_halos && x.halo_median_reconstruction)
        {
            const auto & snapshot=*x.halo_median_reconstruction;
            Require(snapshot.contract==HaloMedianReconstructionContractVersion &&
                !snapshot.joint_run_id.empty() && snapshot.valid_target_count>0,
                "Pure Exclusion Halo median snapshot contract");
            for(const double value:snapshot.global_median)
                Require(std::isfinite(value),"Pure Exclusion global median finiteness");
            Require(snapshot.global_median[0]>=0 && snapshot.global_median[1]>0,
                "Pure Exclusion global median validity");
            Require(snapshot.halos.size()==atoms-targets.size(),"Pure Exclusion Halo median coverage");
            std::size_t previous{}; bool first=true;
            for(const auto & halo:snapshot.halos)
            {
                Require(halo.halo_contributor_index<atoms && (first || previous<halo.halo_contributor_index) &&
                    x.atom_ids[halo.halo_contributor_index]==halo.atom_id &&
                    !std::binary_search(targets.begin(),targets.end(),halo.halo_contributor_index),
                    "Pure Exclusion Halo median identity");
                Require(halo.contract==HaloMedianReconstructionContractVersion &&
                    halo.joint_run_id==snapshot.joint_run_id && halo.donor_count>0,
                    "Pure Exclusion Halo median provenance");
                Require(std::isfinite(halo.amplitude) && halo.amplitude>=0 &&
                    std::isfinite(halo.width) && halo.width>0 && std::isfinite(halo.charge),
                    "Pure Exclusion Halo median parameters");
                if(halo.source==HaloMedianSource::GlobalTargetMedian)
                    Require(halo.donor_count==snapshot.valid_target_count,
                        "Pure Exclusion global median donor count");
                previous=halo.halo_contributor_index; first=false;
            }
        }
        else if(has_halos)
        {
            Require(!IsQualifiedJointAnalysisResult(x) ||
                HasUnavailableJointHaloMedianReconstruction(x),
                "Pure Exclusion unavailable Halo median reconstruction");
        }
        else Require(!x.halo_median_reconstruction,"Pure Exclusion no-Halo median snapshot");
    }
    else Require(!x.halo_median_reconstruction,"unexpected Halo median snapshot");
    if(x.layout)
    {
        std::vector<std::size_t> all_atoms(atoms),all_rows(rows);
        std::iota(all_atoms.begin(),all_atoms.end(),0); std::iota(all_rows.begin(),all_rows.end(),0);
        Require(!x.layout->groups.empty(),"empty observable reduction");
        ValidateLayout(*x.layout,all_atoms,all_rows,x.selection_domain ? &*x.selection_domain : nullptr);
    }
    std::vector<std::size_t> all(atoms); std::iota(all.begin(),all.end(),0);
    const bool target_only_parameters=target_only_fixed_halo || target_only_pure_exclusion;
    const auto & full_atoms=target_only_parameters ? x.selection_domain->target_indices :
        x.layout ? x.layout->full_atoms : all;
    ValidateTarget(x.target_evidence,target_only_parameters ? full_atoms : all,full_atoms,x.selection_domain,
        target_only_parameters ? full_atoms.size() : x.layout ? x.layout->full_atoms.size() : atoms,rows);
    if(x.assembled_state) ValidateState(*x.assembled_state,
        target_only_parameters ? full_atoms.size() : x.layout ? x.layout->full_atoms.size() : atoms,
        x.layout ? x.layout->groups.size() : 0);
    if(x.initialization.valid)
    {
        Require(x.initialization.b.size()==(target_only_parameters ? full_atoms.size() : atoms),"initial width dimensions");
        if(target_only_parameters)
            for(const double width:x.initialization.b)
                Require(std::isfinite(width) && width>0,"invalid initial width");
        else for(std::size_t atom=0;atom<atoms;++atom)
            if(!x.layout || std::binary_search(x.layout->full_atoms.begin(),x.layout->full_atoms.end(),atom))
                Require(std::isfinite(x.initialization.b[atom]) && x.initialization.b[atom]>0,"invalid initial width");
    }
    Require(x.initialization.atoms.size()<=atoms,"initialization dimensions");
    std::set<std::size_t> used_atoms,used_rows;
    std::set<std::string> components;
    for(const auto & c:x.components)
    {
        Require(components.insert(c.id).second,"duplicate component");
        for(auto a:c.atoms) Require(a<atoms && used_atoms.insert(a).second,"component atom mapping");
        for(auto r:c.rows) Require(r<rows && used_rows.insert(r).second,"component row mapping");
        Require(c.layout.has_value()==x.layout.has_value(),"component layout contract");
        if(target_only_parameters)
            for(const auto atom:c.atoms)
                Require(std::binary_search(x.selection_domain->target_indices.begin(),
                    x.selection_domain->target_indices.end(),atom),"fixed halo cannot enter target component");
        if(c.layout)
        {
            ValidateLayout(*c.layout,c.atoms,c.rows,x.selection_domain ? &*x.selection_domain : nullptr);
            for(auto a:c.layout->full_atoms) Require(std::binary_search(x.layout->full_atoms.begin(),x.layout->full_atoms.end(),a),"component full atom disagreement");
            for(auto r:c.layout->informative_rows) Require(std::binary_search(x.layout->informative_rows.begin(),x.layout->informative_rows.end(),r),"component informative row disagreement");
            for(const auto & group:c.layout->groups)
            {
                const auto found=std::find_if(x.layout->groups.begin(),x.layout->groups.end(),[&](const auto & g){return g.row==group.row;});
                Require(found!=x.layout->groups.end() && found->atoms==group.atoms,"component nuisance group disagreement");
            }
        }
        const auto & component_full_atoms=c.layout ? c.layout->full_atoms : c.atoms;
        ValidateTarget(c.target_evidence,c.atoms,component_full_atoms,x.selection_domain,
            c.layout ? c.layout->full_atoms.size() : c.atoms.size(),c.rows.size());
        Require(std::isfinite(c.local_conditioning_seconds) && c.local_conditioning_seconds>=0,
            "local conditioning time");
        for(const auto & block:c.local_conditioning) ValidateBlockConditioning(block);
        if(c.parameter_stability) ValidateParameterStability(*c.parameter_stability);
        if(c.state) ValidateState(*c.state,c.layout ? c.layout->full_atoms.size() : c.atoms.size(),c.layout ? c.layout->groups.size() : 0);
        if(c.candidate_diagnostic) ValidateCandidateDiagnostic(*c.candidate_diagnostic,
            c.layout ? c.layout->full_atoms.size() : c.atoms.size(),c.state.has_value(),reference_free_contract);
    }
}
}
std::string Encode(const JointAnalysisResult & x)
{
    const bool pure_exclusion=x.parameterization_contract==JointPureHaloSupportExclusionParameterizationContract;
    Validate(x,pure_exclusion ? 17 : 16);
    const auto & m=x.metadata; const auto & t=x.costs;
    j::array components;
    for(const auto & c:x.components)
    {
        Object component{
        {"id",c.id},{"stop_reason",c.stop_reason},{"atoms",j::value_from(c.atoms)},{"rows",j::value_from(c.rows)},
        {"search_completed",c.search_completed},{"search_sweeps",c.search_sweeps},
        {"first_order_stationarity_sweep",c.first_order_stationarity_sweep},
        {"confirmed_stationarity_sweep",c.confirmed_stationarity_sweep},
        {"endpoint_certification_sweep",c.endpoint_certification_sweep},
        {"profile_evaluations",c.profile_evaluations},{"reference_evaluations",c.reference_evaluations},
        {"accepted_updates",c.accepted_updates},{"native_status",c.native_status},{"state",OptionalState(c.state)},
        {"candidate_diagnostic",c.candidate_diagnostic ? j::value(CandidateDiagnostic(*c.candidate_diagnostic)) : j::value(nullptr)},
        {"evidence",Checks(c.evidence)},{"legacy_global_rank_evidence",Ranks(c.legacy_global_rank_evidence)},
        {"local_conditioning",BlockConditioning(c.local_conditioning)},
        {"local_conditioning_seconds",Number(c.local_conditioning_seconds)},
        {"local_conditioning_peak_memory_estimate_bytes",OptionalSize(c.local_conditioning_peak_memory_estimate_bytes)},
        {"parameter_stability",ParameterStability(c.parameter_stability)},
        {"regular_certificate",StatusText(c.regular_certificate)},
        { "runtime_convergence",StatusText(c.runtime_convergence)},{"layout",Layout(c.layout)},
        {"target_evidence",Target(c.target_evidence)},{"target_runtime_convergence",StatusText(c.target_runtime_convergence)}};
        if(c.legacy_candidate_evidence)
            component.emplace("legacy_candidate_evidence",LegacyCandidateEvidence(c.legacy_candidate_evidence));
        components.emplace_back(std::move(component));
    }
    Object out{{"schema_version",pure_exclusion ? 17 : 16},{"estimator","joint-components"},{"estimator_contract","guarded-joint-ls-v1"},{"objective_contract","parent-normalized-half-rss-v1"},
        { "support_contract","sphere-fma-v1"},{"parameterization_contract",x.parameterization_contract},{"layout",Layout(x.layout)},{"metadata",Metadata(m)},{"selection_domain",SelectionDomain(x.selection_domain)},
        {"fixed_halos",FixedHalos(x.fixed_halos)},
        {"atom_ids",j::value_from(x.atom_ids)},{"row_ids",j::value_from(x.row_ids)},{"initialization",Initialization(x.initialization)},
        {"costs",Object{{"construction_seconds",t.construction_seconds},{"initialization_seconds",t.initialization_seconds},{"search_seconds",t.search_seconds},
            {"search_reference_seconds",t.search_reference_seconds},{"assessment_seconds",t.assessment_seconds},{"assembly_seconds",t.assembly_seconds}}},
        {"components",std::move(components)},{"assembled_state",OptionalState(x.assembled_state)},{"objective",Number(x.objective)},
        {"available_row_mask",j::value_from(x.available_row_mask)},{"evidence",Checks(x.evidence)},
        {"legacy_global_rank_evidence",Ranks(x.legacy_global_rank_evidence)},
        {"search_completed",x.search_completed},{"observation_scale",x.observation_scale},
        {"regular_certificate",StatusText(x.regular_certificate)},{"runtime_convergence",StatusText(x.runtime_convergence)},
        {"target_evidence",Target(x.target_evidence)},{"target_runtime_convergence",StatusText(x.target_runtime_convergence)}};
    if(pure_exclusion)
    {
        out.emplace("fitting_row_ids",j::value_from(x.fitting_row_ids));
        out.emplace("fitting_to_original_contributors",j::value_from(x.fitting_to_original_contributors));
        out.emplace("fitting_to_original_rows",j::value_from(x.fitting_to_original_rows));
        out.emplace("halo_median_reconstruction",HaloMedianSnapshot(x.halo_median_reconstruction));
    }
    return j::serialize(out);
}
JointAnalysisResult Decode(std::string_view text)
{
    j::parse_options options; options.numbers=j::number_precision::precise;
    const auto parsed=j::parse(text,{},options); const auto & o=parsed.as_object();
    const int version=Read<int>(o,"schema_version");
    Require(version==3 || version==4 || version==5 || version==6 || version==7 || version==8 || version==9 || version==10 || version==11 || version==12 || version==13 || version==14 || version==15 || version==16 || version==17,
        "unsupported result schema version (expected 3 through 17)");
    if(version>=15)
        for(const auto & fixed:o.at("fixed_halos").as_array())
            Require(fixed.as_object().contains("first_stage_qualification"),
                "schema 15 fixed halo requires qualification provenance");
    Require(Read<std::string>(o,"estimator")=="joint-components" &&
        Read<std::string>(o,"estimator_contract")=="guarded-joint-ls-v1" &&
        Read<std::string>(o,"objective_contract")=="parent-normalized-half-rss-v1" &&
        Read<std::string>(o,"support_contract")=="sphere-fma-v1","unsupported estimator contract");
    JointAnalysisResult x;
    if(version>=4) {x.layout=ReadLayout(o); x.parameterization_contract=Read<std::string>(o,"parameterization_contract");}
    x.metadata=ReadMetadata(o.at("metadata").as_object());
    x.selection_domain=ReadSelectionDomain(o.at("selection_domain"));
    if(version>=14) x.fixed_halos=ReadFixedHalos(o.at("fixed_halos"));
    x.atom_ids=Read<std::vector<std::string>>(o,"atom_ids"); x.row_ids=Read<std::vector<std::string>>(o,"row_ids");
    if(version>=17)
    {
        x.fitting_row_ids=Read<std::vector<std::string>>(o,"fitting_row_ids");
        x.fitting_to_original_contributors=Read<std::vector<std::size_t>>(o,"fitting_to_original_contributors");
        x.fitting_to_original_rows=Read<std::vector<std::size_t>>(o,"fitting_to_original_rows");
        x.halo_median_reconstruction=ReadHaloMedianSnapshot(o.at("halo_median_reconstruction"));
    }
    x.initialization=ReadInitialization(o.at("initialization"));
    const auto & t=o.at("costs").as_object();
    x.costs={Read<double>(t,"initialization_seconds"),Read<double>(t,"search_seconds"),Read<double>(t,"search_reference_seconds"),
        Read<double>(t,"assessment_seconds"),Read<double>(t,"assembly_seconds")};
    x.costs.construction_seconds=Read<double>(t,"construction_seconds");
    for(const auto & item:o.at("components").as_array())
    {
        const auto & c=item.as_object(); JointAnalysisComponent component;
        if(version>=4) component.layout=ReadLayout(c);
        component.id=Read<std::string>(c,"id"); component.stop_reason=Read<std::string>(c,"stop_reason");
        component.atoms=Read<std::vector<std::size_t>>(c,"atoms"); component.rows=Read<std::vector<std::size_t>>(c,"rows");
        component.search_completed=Read<bool>(c,"search_completed");
        component.search_sweeps=c.contains("search_sweeps") ? Read<std::size_t>(c,"search_sweeps") : 0;
        component.first_order_stationarity_sweep=c.contains("first_order_stationarity_sweep") ?
            Read<std::size_t>(c,"first_order_stationarity_sweep") : 0;
        component.confirmed_stationarity_sweep=c.contains("confirmed_stationarity_sweep") ?
            Read<std::size_t>(c,"confirmed_stationarity_sweep") : 0;
        component.endpoint_certification_sweep=c.contains("endpoint_certification_sweep") ?
            Read<std::size_t>(c,"endpoint_certification_sweep") : 0;
        component.profile_evaluations=Read<int>(c,"profile_evaluations");
        component.reference_evaluations=Read<int>(c,"reference_evaluations"); component.accepted_updates=Read<int>(c,"accepted_updates");
        component.native_status=Read<int>(c,"native_status");
        if(!c.at("state").is_null()) component.state=ReadState(c.at("state"));
        if(version>=6)
        {
            component.candidate_diagnostic=ReadCandidateDiagnostic(c.at("candidate_diagnostic"),version>=7,version>=16);
            component.legacy_candidate_evidence=ReadLegacyCandidateEvidence(c);
        }
        component.evidence=ReadChecks(c.at("evidence"));
        component.legacy_global_rank_evidence=ReadRanks(c.at(version>=12 ? "legacy_global_rank_evidence" : "ranks"),version>=9);
        if(version>=12)
        {
            component.local_conditioning=ReadBlockConditioning(c.at("local_conditioning"));
            component.local_conditioning_seconds=Read<double>(c,"local_conditioning_seconds");
            component.local_conditioning_peak_memory_estimate_bytes=
                ReadOptionalSize(c.at("local_conditioning_peak_memory_estimate_bytes"));
            component.parameter_stability=ReadParameterStability(c.at("parameter_stability"));
        }
        component.regular_certificate=Status(c.at("regular_certificate")); component.runtime_convergence=Status(c.at("runtime_convergence"));
        if(version>=5) {component.target_evidence=ReadTarget(c.at("target_evidence"),version); component.target_runtime_convergence=Status(c.at("target_runtime_convergence"));}
        x.components.push_back(std::move(component));
    }
    if(!o.at("assembled_state").is_null()) x.assembled_state=ReadState(o.at("assembled_state"));
    x.objective=OptionalNumber(o.at("objective")); x.available_row_mask=Read<std::vector<bool>>(o,"available_row_mask");
    x.evidence=ReadChecks(o.at("evidence"));
    x.legacy_global_rank_evidence=ReadRanks(o.at(version>=12 ? "legacy_global_rank_evidence" : "ranks"),version>=9);
    x.search_completed=Read<bool>(o,"search_completed"); x.observation_scale=Read<double>(o,"observation_scale");
    x.regular_certificate=Status(o.at("regular_certificate")); x.runtime_convergence=Status(o.at("runtime_convergence"));
    if(version>=5) {x.target_evidence=ReadTarget(o.at("target_evidence"),version); x.target_runtime_convergence=Status(o.at("target_runtime_convergence"));}
    Validate(x,version); return x;
}
}
