#include "JointResultJson.hpp"
#include "PersistenceSchema.hpp"
#include "core/detail/joint_component/ConvergenceCriteria.hpp"
#include <boost/json.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
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
void RequireOnlyKeys(const Object &,std::initializer_list<const char *>,const char *,
    std::initializer_list<const char *> optional_keys={});
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
    if(source==HaloMedianSource::GeometryLocalTargetMedian) return "geometry-local-target-median";
    throw std::invalid_argument("Invalid current Halo median source.");
}
HaloMedianSource HaloMedianSourceValue(const j::value & value)
{
    if(j::value_to<std::string>(value)=="geometry-local-target-median")
        return HaloMedianSource::GeometryLocalTargetMedian;
    throw std::invalid_argument("Invalid current Halo median source.");
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
        {"width_gradient",j::value_from(x.width_gradient)},{"objective",x.objective}};
}
JointState ReadState(const j::value & v)
{
    const auto & o=v.as_object();
    return {Read<std::vector<double>>(o,"ac"),Read<std::vector<double>>(o,"b"),
        Read<std::vector<double>>(o,"log_b"),Read<std::vector<double>>(o,"width_gradient"),
        Read<double>(o,"objective")};
}
j::value OptionalState(const std::optional<JointState> & x) {return x ? j::value(State(*x)) : j::value(nullptr);}
const char * CandidateStatusText(const JointCandidateDiagnostic & x)
{return x.certification_status==JointCheckStatus::Passed ? "certified" : "uncertified";}
Object CandidateDiagnostic(const JointCandidateDiagnostic & x)
{
    const auto & strict=core::joint_component::StrictJointConvergenceCriteria;
    Object out{{"candidate_status",CandidateStatusText(x)},
        {"candidate_available",x.candidate_available},
        {"candidate_eta",j::value_from(x.candidate_eta)},
        {"candidate_ac",j::value_from(x.candidate_ac)},
        {"candidate_objective",Number(x.candidate_objective)},
        {"candidate_ac_kkt",Number(x.candidate_ac_kkt)},
        {"candidate_width_gradient_inf_norm",Number(x.candidate_width_gradient_inf_norm)},
        {"convergence_criteria_id",x.convergence_criteria_id.empty() ? "strict-v2" : x.convergence_criteria_id},
        {"global_ac_kkt_threshold",x.global_ac_kkt_threshold.value_or(strict.global_ac_kkt)},
        {"global_width_gradient_threshold",x.global_width_gradient_threshold.value_or(strict.global_width_gradient)},
        {"candidate_feasible",x.candidate_feasible ? j::value(*x.candidate_feasible) : j::value(nullptr)},
        {"certification_status",StatusText(x.certification_status)},
        {"endpoint_verification_status",StatusText(x.endpoint_verification_status)},
        {"numerical_qualification_status",StatusText(x.numerical_qualification_status)},
        {"verification_reason",x.verification_reason},
        {"resource_limit_reason",x.resource_limit_reason}};
    return out;
}
std::optional<JointCandidateDiagnostic> ReadCandidateDiagnostic(const j::value & value)
{
    if(value.is_null()) return {};
    const auto & o=value.as_object(); JointCandidateDiagnostic x;
    x.candidate_available=Read<bool>(o,"candidate_available");
    x.candidate_eta=Read<std::vector<double>>(o,"candidate_eta");
    x.candidate_ac=Read<std::vector<double>>(o,"candidate_ac");
    x.candidate_objective=OptionalNumber(o.at("candidate_objective"));
    x.candidate_ac_kkt=OptionalNumber(o.at("candidate_ac_kkt"));
    x.candidate_width_gradient_inf_norm=OptionalNumber(o.at("candidate_width_gradient_inf_norm"));
    x.convergence_criteria_id=Read<std::string>(o,"convergence_criteria_id");
    x.global_ac_kkt_threshold=Read<double>(o,"global_ac_kkt_threshold");
    x.global_width_gradient_threshold=Read<double>(o,"global_width_gradient_threshold");
    x.candidate_feasible=o.at("candidate_feasible").is_null() ? std::nullopt :
        std::optional<bool>(Read<bool>(o,"candidate_feasible"));
    x.certification_status=Status(o.at("certification_status"));
    x.endpoint_verification_status=Status(o.at("endpoint_verification_status"));
    x.numerical_qualification_status=Status(o.at("numerical_qualification_status"));
    x.verification_reason=Read<std::string>(o,"verification_reason");
    x.resource_limit_reason=Read<std::string>(o,"resource_limit_reason");
    const auto candidate_status=Read<std::string>(o,"candidate_status");
    if(candidate_status!=CandidateStatusText(x))
        throw std::invalid_argument("Invalid saved joint candidate status: "+candidate_status);
    return x;
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
        RequireOnlyKeys(o,{"name","scope","status","value","threshold","reason"},
            "unknown current Joint evidence field");
        out.push_back({Read<std::string>(o,"name"),Status(o.at("status")),Scope(o),OptionalNumber(o.at("value")),
            OptionalNumber(o.at("threshold")),Read<std::string>(o,"reason")});
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
    RequireOnlyKeys(o,{"valid","reason","b","atoms","data_scope"},"unknown current Joint initialization field");
    x.valid=Read<bool>(o,"valid"); x.reason=Read<std::string>(o,"reason"); x.b=DiagnosticValues(o.at("b"));
    x.data_scope=Read<std::string>(o,"data_scope");
    for(const auto & item:o.at("atoms").as_array())
    {
        const auto & a=item.as_object(); JointInitializationAtom atom;
        RequireOnlyKeys(a,{"id","ols","mdpde","alpha","sample_count","native_status","reason","original_b",
            "used_b","seed_source","donor_count","first_stage_fit_qualified","first_stage_qualification",
            "first_stage_recovery"},"unknown current FirstStage initialization field");
        atom.id=Read<std::string>(a,"id"); atom.reason=Read<std::string>(a,"reason");
        const auto ols=DiagnosticValues(a.at("ols")),mdpde=DiagnosticValues(a.at("mdpde"));
        if(ols.size()!=3 || mdpde.size()!=3) throw std::invalid_argument("Invalid joint initialization dimensions.");
        std::copy(ols.begin(),ols.end(),atom.ols.begin()); std::copy(mdpde.begin(),mdpde.end(),atom.mdpde.begin());
        atom.alpha=OptionalNumber(a.at("alpha")).value_or(std::numeric_limits<double>::quiet_NaN());
        atom.sample_count=Read<std::size_t>(a,"sample_count");
        atom.first_stage_fit_qualified=Read<bool>(a,"first_stage_fit_qualified");
        if(!a.at("first_stage_qualification").is_null())
            atom.first_stage_qualification=Qualification(a.at("first_stage_qualification"));
        if(!a.at("first_stage_recovery").is_null())
            atom.first_stage_recovery=ReadFailedOnlyRecovery(a.at("first_stage_recovery"));
        atom.original_b=OptionalNumber(a.at("original_b")); atom.used_b=OptionalNumber(a.at("used_b"));
        atom.seed_source=Read<std::string>(a,"seed_source"); atom.donor_count=Read<std::size_t>(a,"donor_count");
        if(!a.at("native_status").is_null()) atom.native_status=Read<int>(a,"native_status");
        x.atoms.push_back(std::move(atom));
    }
    return x;
}
void Require(bool condition, const char * message)
{if(!condition) throw std::invalid_argument(std::string("Invalid joint result: ")+message);}
void RequireOnlyKeys(const Object & object,std::initializer_list<const char *> keys,const char * message,
    std::initializer_list<const char *> optional_keys)
{
    for(const auto key:keys) Require(object.contains(key),message);
    for(auto it=object.begin();it!=object.end();++it)
        Require(std::find(keys.begin(),keys.end(),it->key())!=keys.end() ||
            std::find(optional_keys.begin(),optional_keys.end(),it->key())!=optional_keys.end(),message);
}
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
        x.reason=Read<std::string>(o,"reason"); x.conditioning_classification=Read<std::string>(o,"conditioning_classification");
        x.core_atoms=Read<std::vector<std::size_t>>(o,"core_atoms");
        x.core_atom_count=Read<std::size_t>(o,"core_atom_count");
        x.affected_row_count=Read<std::size_t>(o,"affected_row_count");
        x.free_ac_columns=Read<std::size_t>(o,"free_ac_columns");
        x.active_amplitude_count=Read<std::size_t>(o,"active_amplitude_count");
        x.free_ac_rank=ReadOptionalSize(o.at("free_ac_rank"));
        x.projected_width_rank=ReadOptionalSize(o.at("projected_width_rank"));
        x.profile_rank=ReadOptionalSize(o.at("profile_rank"));
        x.active_width_bound_targets=ReadOptionalSize(o.at("active_width_bound_targets"));
        x.workspace_peak_memory_estimate_bytes=ReadOptionalSize(o.at("workspace_peak_memory_estimate_bytes"));
        x.free_ac_sigma_min=OptionalNumber(o.at("free_ac_sigma_min"));
        x.free_ac_sigma_max=OptionalNumber(o.at("free_ac_sigma_max"));
        x.free_ac_condition=OptionalNumber(o.at("free_ac_condition"));
        x.free_ac_rank_threshold=OptionalNumber(o.at("free_ac_rank_threshold"));
        x.free_ac_rank_threshold_margin_relative=OptionalNumber(o.at("free_ac_rank_threshold_margin_relative"));
        x.local_ac_projected_kkt=OptionalNumber(o.at("local_ac_projected_kkt"));
        x.local_projected_width_kkt=OptionalNumber(o.at("local_projected_width_kkt"));
        x.free_ac_threshold_sensitive=ReadOptionalBoolean(o.at("free_ac_threshold_sensitive"));
        x.projected_width_condition=OptionalNumber(o.at("projected_width_condition"));
        x.projected_width_sigma_min=OptionalNumber(o.at("projected_width_sigma_min"));
        x.projected_width_sigma_max=OptionalNumber(o.at("projected_width_sigma_max"));
        x.projected_width_rank_threshold=OptionalNumber(o.at("projected_width_rank_threshold"));
        x.projected_width_rank_threshold_margin_relative=OptionalNumber(o.at("projected_width_rank_threshold_margin_relative"));
        x.projected_width_threshold_sensitive=ReadOptionalBoolean(o.at("projected_width_threshold_sensitive"));
        x.profile_condition=OptionalNumber(o.at("profile_condition"));
        x.profile_sigma_min=OptionalNumber(o.at("profile_sigma_min"));
        x.profile_sigma_max=OptionalNumber(o.at("profile_sigma_max"));
        x.profile_rank_threshold=OptionalNumber(o.at("profile_rank_threshold"));
        x.profile_rank_threshold_margin_relative=OptionalNumber(o.at("profile_rank_threshold_margin_relative"));
        x.profile_threshold_sensitive=ReadOptionalBoolean(o.at("profile_threshold_sensitive"));
        x.free_ac_status=Read<std::string>(o,"free_ac_status"); x.free_ac_reason=Read<std::string>(o,"free_ac_reason");
        x.projected_width_status=Read<std::string>(o,"projected_width_status");
        x.projected_width_reason=Read<std::string>(o,"projected_width_reason");
        x.free_ac_conditioning=Read<std::string>(o,"free_ac_conditioning");
        x.projected_width_conditioning=Read<std::string>(o,"projected_width_conditioning");
        x.profile_conditioning=Read<std::string>(o,"profile_conditioning");
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
    x.probe_selected=Read<bool>(o,"probe_selected");
    x.prediction_status=Read<std::string>(o,"prediction_status");
    x.objective_consistency_status=Read<std::string>(o,"objective_consistency_status");
    x.log_width_parameter_change_inf=OptionalNumber(o.at("log_width_parameter_change_inf"));
    x.scaled_ac_parameter_change_inf=OptionalNumber(o.at("scaled_ac_parameter_change_inf"));
    x.local_prediction_change_relative=OptionalNumber(o.at("local_prediction_change_relative"));
    x.objective_difference=OptionalNumber(o.at("objective_difference")); x.cpu_seconds=Read<double>(o,"cpu_seconds");
    x.baseline_local_ac_kkt=OptionalNumber(o.at("baseline_local_ac_kkt"));
    x.baseline_projected_width_kkt=OptionalNumber(o.at("baseline_projected_width_kkt"));
    x.baseline_width_bound_slack=OptionalNumber(o.at("baseline_width_bound_slack"));
    x.baseline_local_objective=OptionalNumber(o.at("baseline_local_objective"));
    x.baseline_local_residual_rmse_normalized=OptionalNumber(o.at("baseline_local_residual_rmse_normalized"));
    x.baseline_active_ac_constraints=ReadOptionalSize(o.at("baseline_active_ac_constraints"));
    x.baseline_active_width_bound_targets=ReadOptionalSize(o.at("baseline_active_width_bound_targets"));
    x.baseline_numerical_validity_status=Read<std::string>(o,"baseline_numerical_validity_status");
    x.baseline_conditional_optimality_status=Read<std::string>(o,"baseline_conditional_optimality_status");
    x.baseline_eligibility_status=Read<std::string>(o,"baseline_eligibility_status");
    x.baseline_eligibility_reason=Read<std::string>(o,"baseline_eligibility_reason");
    x.baseline_independent_replay_status=Read<std::string>(o,"baseline_independent_replay_status");
    x.baseline_independent_replay_reason=Read<std::string>(o,"baseline_independent_replay_reason");
    x.conditional_reprofile_status=Read<std::string>(o,"conditional_reprofile_status");
    x.conditional_reprofile_local_ac_kkt=OptionalNumber(o.at("conditional_reprofile_local_ac_kkt"));
    x.conditional_reprofile_projected_width_kkt=OptionalNumber(o.at("conditional_reprofile_projected_width_kkt"));
    x.conditional_reprofile_objective_difference=OptionalNumber(o.at("conditional_reprofile_objective_difference"));
    x.conditional_reprofile_scaled_ac_parameter_difference=OptionalNumber(o.at("conditional_reprofile_scaled_ac_parameter_difference"));
    x.conditional_reprofile_prediction_change_relative=OptionalNumber(o.at("conditional_reprofile_prediction_change_relative"));
    x.conditional_reprofile_independent_replay_status=Read<std::string>(o,"conditional_reprofile_independent_replay_status");
    x.conditional_reprofile_independent_replay_reason=Read<std::string>(o,"conditional_reprofile_independent_replay_reason");
    for(const auto & probe_value:o.at("probes").as_array())
    {
        const auto & p=probe_value.as_object(); JointStabilityProbe probe;
        probe.direction=Read<std::string>(p,"direction"); probe.source=Read<std::string>(p,"source");
        probe.status=Read<std::string>(p,"status");
        probe.minus_seed_status=Read<std::string>(p,"minus_seed_status");
        probe.plus_seed_status=Read<std::string>(p,"plus_seed_status");
        probe.minus_seed_reason=Read<std::string>(p,"minus_seed_reason");
        probe.plus_seed_reason=Read<std::string>(p,"plus_seed_reason");
        probe.unit_direction=Read<std::vector<double>>(p,"unit_direction");
        probe.evaluations=Read<std::size_t>(p,"evaluations");
        probe.parameter_change=OptionalNumber(p.at("parameter_change"));
        probe.prediction_change=OptionalNumber(p.at("prediction_change"));
        probe.objective_change=OptionalNumber(p.at("objective_change"));
        probe.minus_requested_displacement_inf=OptionalNumber(p.at("minus_requested_displacement_inf"));
        probe.minus_effective_displacement_inf=OptionalNumber(p.at("minus_effective_displacement_inf"));
        probe.plus_requested_displacement_inf=OptionalNumber(p.at("plus_requested_displacement_inf"));
        probe.plus_effective_displacement_inf=OptionalNumber(p.at("plus_effective_displacement_inf"));
        x.probes.push_back(std::move(probe));
    }
    x.perturbation_directions=Read<std::size_t>(o,"perturbation_directions");
    x.perturbation_seeds=Read<std::size_t>(o,"perturbation_seeds");
    x.perturbation_solves=Read<std::size_t>(o,"perturbation_solves");
    x.perturbation_evaluations=Read<std::size_t>(o,"perturbation_evaluations");
    return x;
}
j::value ParameterStability(const std::optional<JointParameterStability> & value)
{
    if(!value) return nullptr;
    const auto & x=*value; j::array blocks;
    for(const auto & block:x.blocks) blocks.emplace_back(StabilityBlock(block));
    return Object{{"perturbation_status",x.perturbation_status},
        {"blocks",std::move(blocks)},{"seconds",x.seconds},{"cpu_seconds",x.cpu_seconds},
        {"perturbation_directions",x.perturbation_directions},{"perturbation_seeds",x.perturbation_seeds},
        {"perturbation_solves",x.perturbation_solves},{"perturbation_evaluations",x.perturbation_evaluations},
        {"workspace_peak_memory_estimate_bytes",OptionalSize(x.workspace_peak_memory_estimate_bytes)}};
}
std::optional<JointParameterStability> ReadParameterStability(const j::value & value)
{
    if(value.is_null()) return std::nullopt;
    const auto & o=value.as_object(); JointParameterStability x;
    RequireOnlyKeys(o,{"perturbation_status","blocks","seconds","cpu_seconds","perturbation_directions",
        "perturbation_seeds","perturbation_solves","perturbation_evaluations",
        "workspace_peak_memory_estimate_bytes"},"unknown current parameter stability field");
    x.perturbation_status=Read<std::string>(o,"perturbation_status");
    for(const auto & block:o.at("blocks").as_array()) x.blocks.push_back(ReadStabilityBlock(block));
    x.seconds=Read<double>(o,"seconds"); x.cpu_seconds=Read<double>(o,"cpu_seconds");
    x.perturbation_directions=Read<std::size_t>(o,"perturbation_directions");
    x.perturbation_seeds=Read<std::size_t>(o,"perturbation_seeds");
    x.perturbation_solves=Read<std::size_t>(o,"perturbation_solves");
    x.perturbation_evaluations=Read<std::size_t>(o,"perturbation_evaluations");
    x.workspace_peak_memory_estimate_bytes=ReadOptionalSize(o.at("workspace_peak_memory_estimate_bytes"));
    return x;
}


j::value HaloMedianSnapshot(const std::optional<JointHaloMedianSnapshot> & value)
{
    if(!value) return nullptr;
    const auto & snapshot=*value; j::array halos;
    for(const auto & halo:snapshot.halos)
    {
        Object encoded{{"halo_contributor_index",halo.halo_contributor_index},
            {"atom_id",halo.atom_id},{"amplitude",halo.amplitude},
            {"width",halo.width},{"charge",halo.charge},{"source",HaloMedianSourceText(halo.source)},
            {"donor_count",halo.donor_count},{"joint_run_id",halo.joint_run_id},
            {"contract",halo.contract},{"donor_atom_ids",j::value_from(halo.donor_atom_ids)}};
        halos.emplace_back(std::move(encoded));
    }
    Object encoded{{"contract",snapshot.contract},{"joint_run_id",snapshot.joint_run_id},
        {"global_median",j::value_from(snapshot.global_median)},
        {"valid_target_count",snapshot.valid_target_count},{"halos",std::move(halos)},
        {"neighbor_count",snapshot.neighbor_count},{"donor_policy",snapshot.donor_policy}};
    return encoded;
}
std::optional<JointHaloMedianSnapshot> ReadHaloMedianSnapshot(const j::value & value)
{
    if(value.is_null()) return std::nullopt;
    const auto & object=value.as_object(); JointHaloMedianSnapshot snapshot;
    RequireOnlyKeys(object,{"contract","joint_run_id","global_median","valid_target_count","halos",
        "neighbor_count","donor_policy"},"unknown current Halo median field");
    snapshot.contract=Read<std::string>(object,"contract");
    snapshot.joint_run_id=Read<std::string>(object,"joint_run_id");
    snapshot.global_median=Read<std::array<double,3>>(object,"global_median");
    snapshot.valid_target_count=Read<std::size_t>(object,"valid_target_count");
    snapshot.neighbor_count=Read<std::size_t>(object,"neighbor_count");
    snapshot.donor_policy=Read<std::string>(object,"donor_policy");
    for(const auto & item:object.at("halos").as_array())
    {
        const auto & halo=item.as_object();
        RequireOnlyKeys(halo,{"halo_contributor_index","atom_id","amplitude","width","charge","source",
            "donor_count","joint_run_id","contract","donor_atom_ids"},"unknown current Halo median point field");
        JointHaloMedianReconstruction decoded;
        decoded.halo_contributor_index=Read<std::size_t>(halo,"halo_contributor_index");
        decoded.atom_id=Read<std::string>(halo,"atom_id");
        decoded.amplitude=Read<double>(halo,"amplitude");
        decoded.width=Read<double>(halo,"width");
        decoded.charge=Read<double>(halo,"charge");
        decoded.source=HaloMedianSourceValue(halo.at("source"));
        decoded.donor_count=Read<std::size_t>(halo,"donor_count");
        decoded.joint_run_id=Read<std::string>(halo,"joint_run_id");
        decoded.contract=Read<std::string>(halo,"contract");
        decoded.donor_atom_ids=Read<std::vector<std::string>>(halo,"donor_atom_ids");
        snapshot.halos.push_back(std::move(decoded));
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
            {"fixed_neighbor_core_atoms",OptionalSize(p.fixed_neighbor_core_atoms)},
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
        RequireOnlyKeys(n,{"requested","applied","divisor"},"unknown current map normalization field");
        out.map_normalization=JointMapNormalization{Read<bool>(n,"requested"),Read<bool>(n,"applied"),Read<double>(n,"divisor")};
    }
    if(!m.at("software").is_null())
    {
        const auto & p=m.at("software").as_object();
        RequireOnlyKeys(p,{"version","source_sha256","configuration_sha256","build_sha256"},
            "unknown current software provenance field");
        out.software=JointSoftwareProvenance{Read<std::string>(p,"version"),Read<std::string>(p,"source_sha256"),
            Read<std::string>(p,"configuration_sha256"),Read<std::string>(p,"build_sha256")};
    }
    if(!m.at("solver").is_null())
    {
        const auto & p=m.at("solver").as_object();
        RequireOnlyKeys(p,{"search_method","contract_version","sparse_backend","fixed_neighbor_core_atoms",
            "fixed_neighbor_policy_version","fixed_neighbor_maximum_sweeps","fixed_neighbor_order",
            "fixed_neighbor_local_search"},"unknown current solver provenance field");
        JointSolverProvenance solver;
        solver.search_method=Read<std::string>(p,"search_method");
        solver.contract_version=ReadOptionalText(p.at("contract_version"));
        solver.sparse_backend=ReadOptionalText(p.at("sparse_backend"));
        solver.fixed_neighbor_core_atoms=ReadOptionalSize(p.at("fixed_neighbor_core_atoms"));
        solver.fixed_neighbor_policy_version=ReadOptionalText(p.at("fixed_neighbor_policy_version"));
        solver.fixed_neighbor_maximum_sweeps=ReadOptionalSize(p.at("fixed_neighbor_maximum_sweeps"));
        solver.fixed_neighbor_order=ReadOptionalText(p.at("fixed_neighbor_order"));
        solver.fixed_neighbor_local_search=ReadOptionalText(p.at("fixed_neighbor_local_search"));
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
        Require(p.search_method=="fixed-neighbor" &&
            p.contract_version==std::optional<std::string>{std::string(JointSolverProvenanceContractVersion)} &&
            p.sparse_backend==std::optional<std::string>{"EIGEN"} &&
            p.fixed_neighbor_core_atoms==std::optional<std::size_t>{12} &&
            p.fixed_neighbor_policy_version==std::optional<std::string>{std::string(FixedNeighborPolicyContractVersion)} &&
            p.fixed_neighbor_maximum_sweeps==std::optional<std::size_t>{200} &&
            p.fixed_neighbor_order==std::optional<std::string>{"forward"} &&
            p.fixed_neighbor_local_search==std::optional<std::string>{std::string(FixedNeighborBoundedLocalSearchName)},
            "invalid current FixedNeighbor provenance");
    }
    if(m.map_normalization)
    {
        const auto & n=*m.map_normalization;
        Require(std::isfinite(n.divisor) && n.divisor>0,"invalid normalization divisor");
        Require(n.applied ? n.requested && !m.simulation : n.divisor==1,"inconsistent normalization");
    }
}
void ValidateState(const JointState & x,std::size_t atoms)
{
    Require(x.ac.size()==2*atoms && x.b.size()==atoms && x.log_b.size()==atoms && x.width_gradient.size()==atoms,"state dimensions");
    for(const auto * values:{&x.ac,&x.b,&x.log_b,&x.width_gradient})
        for(double value:*values) Require(std::isfinite(value),"nonfinite state");
    for(std::size_t atom=0;atom<atoms;++atom)
        Require(x.ac[2*atom]>=0 && x.b[atom]>0,"infeasible target parameters");
    Require(std::isfinite(x.objective),"nonfinite objective");
}
void ValidateCandidateDiagnostic(const JointCandidateDiagnostic & x,std::size_t atoms,bool trusted_state)
{
    const auto & strict=core::joint_component::StrictJointConvergenceCriteria;
    const auto & practical=core::joint_component::PracticalJointConvergenceCriteria;
    const bool default_criteria=x.convergence_criteria_id.empty() && !x.global_ac_kkt_threshold &&
        !x.global_width_gradient_threshold;
    const bool complete_criteria=!x.convergence_criteria_id.empty() && x.global_ac_kkt_threshold &&
        x.global_width_gradient_threshold;
    Require(default_criteria || complete_criteria,"incomplete convergence criteria provenance");
    const auto & criteria_id=default_criteria ? std::string_view{"strict-v2"} :
        std::string_view{x.convergence_criteria_id};
    const auto ac_threshold=default_criteria ? strict.global_ac_kkt : *x.global_ac_kkt_threshold;
    const auto width_threshold=default_criteria ? strict.global_width_gradient : *x.global_width_gradient_threshold;
    const bool current_criteria=criteria_id=="strict-v2" || criteria_id=="practical-v2";
    Require(current_criteria,
        "invalid convergence criteria identifier");
    Require(std::isfinite(ac_threshold) && ac_threshold>0 && std::isfinite(width_threshold) && width_threshold>0,
        "invalid convergence criteria thresholds");
    if(criteria_id=="strict-v2")
        Require(ac_threshold==strict.global_ac_kkt && width_threshold==strict.global_width_gradient,
            "strict convergence criteria do not match their identifier");
    if(criteria_id=="practical-v2")
        Require(ac_threshold==practical.global_ac_kkt && width_threshold==practical.global_width_gradient,
            "practical convergence criteria do not match their identifier");
    Require(atoms<=std::numeric_limits<std::size_t>::max()/2,"candidate dimensions");
    if(x.candidate_available)
    {
        Require(x.candidate_eta.size()==atoms && x.candidate_ac.size()==2*atoms && x.candidate_objective,
            "candidate dimensions");
        for(double eta:x.candidate_eta)
            Require(std::isfinite(eta) && std::isfinite(std::exp(eta)) && std::exp(eta)>0,"invalid candidate width");
        for(std::size_t atom=0;atom<atoms;++atom)
            Require(std::isfinite(x.candidate_ac[2*atom]) && x.candidate_ac[2*atom]>=0 &&
                std::isfinite(x.candidate_ac[2*atom+1]),"infeasible candidate A/C");
        Require(std::isfinite(*x.candidate_objective),"nonfinite candidate objective");
    }
    else Require(x.candidate_eta.empty() && x.candidate_ac.empty() && !x.candidate_objective,
        "unavailable candidate contains parameters");
    for(const auto * value:{&x.candidate_ac_kkt,&x.candidate_width_gradient_inf_norm})
        Require(!*value || (std::isfinite(**value) && **value>=0),"invalid candidate diagnostic value");
    Require((x.certification_status==JointCheckStatus::Passed)==trusted_state,
        "candidate certification disagrees with trusted state");
    Require((x.numerical_qualification_status==JointCheckStatus::Passed)==trusted_state,
        "candidate numerical qualification disagrees with trusted state");
    Require(x.certification_status!=JointCheckStatus::Passed || x.candidate_available,
        "certified candidate is unavailable");
    Require(x.numerical_qualification_status!=JointCheckStatus::Passed || x.candidate_available,
        "qualified candidate is unavailable");
    Require(x.endpoint_verification_status!=JointCheckStatus::Passed || x.candidate_available,
        "verified candidate is unavailable");
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
    Require(!x.perturbation_status.empty(),"parameter stability status");
    Require(std::isfinite(x.seconds) && x.seconds>=0 && std::isfinite(x.cpu_seconds) && x.cpu_seconds>=0,
        "parameter stability time");
    for(const auto & block:x.blocks) ValidateStabilityBlock(block);
}
j::value SelectionDomain(const std::optional<JointSelectionDomain> & domain)
{
    if(!domain) return nullptr;
    const auto & provenance=*domain->pure_halo_support_exclusion;
    const auto pure_halo_support_exclusion=Object{{"original_rows",provenance.original_rows},
        {"fit_rows",provenance.fit_rows},{"excluded_rows",provenance.excluded_rows},
        {"candidate_halos",provenance.candidate_halos},{"retained_halos",provenance.retained_halos},
        {"zero_support_halos",provenance.zero_support_halos}};
    Object out{{"contract",domain->contract},{"target_indices",j::value_from(domain->target_indices)},
        {"observation_radius",domain->observation_radius},{"support_radius",domain->support_radius},
        {"contributor_policy",domain->contributor_policy},
        {"pure_halo_support_exclusion",pure_halo_support_exclusion}};
    return out;
}
std::optional<JointSelectionDomain> ReadSelectionDomain(const j::value & value)
{
    if(value.is_null()) return std::nullopt;
    const auto & o=value.as_object();
    RequireOnlyKeys(o,{"contract","target_indices","observation_radius","support_radius","contributor_policy",
        "pure_halo_support_exclusion"},"unknown current selection domain field");
    JointSelectionDomain out;
    out.contract=Read<std::string>(o,"contract");
    out.target_indices=Read<std::vector<std::size_t>>(o,"target_indices");
    out.observation_radius=Read<double>(o,"observation_radius");
    out.support_radius=Read<double>(o,"support_radius");
    out.contributor_policy=Read<std::string>(o,"contributor_policy");
    Require(!o.at("pure_halo_support_exclusion").is_null(),"missing Pure Exclusion provenance");
    const auto & item=o.at("pure_halo_support_exclusion").as_object();
    RequireOnlyKeys(item,{"original_rows","fit_rows","excluded_rows","candidate_halos","retained_halos",
        "zero_support_halos"},"unknown current Pure Exclusion provenance field");
    out.pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        Read<std::size_t>(item,"original_rows"),Read<std::size_t>(item,"fit_rows"),
        Read<std::size_t>(item,"excluded_rows"),Read<std::size_t>(item,"candidate_halos"),
        Read<std::size_t>(item,"retained_halos"),Read<std::size_t>(item,"zero_support_halos")};
    return out;
}


j::value Target(const std::optional<JointTargetEvidence> & value)
{
    if(!value) return nullptr;
    const auto & t=*value; j::array point_atoms;
    for(const auto & a:t.point_atoms)
        point_atoms.emplace_back(Object{{"atom",a.atom},{"component_convergence",StatusText(a.component_convergence)},
            {"parameter_validity",StatusText(a.parameter_validity)},{"local_conditioning",a.local_conditioning},
            {"parameter_stability",a.parameter_stability},{"point_quality",a.point_quality},
            {"local_conditioning_blocks",j::value_from(a.local_conditioning_blocks)},
            {"parameter_stability_blocks",j::value_from(a.parameter_stability_blocks)},{"reason",a.reason}});
    return Object{{"contract",t.contract},{"checks",Checks(t.checks)},
        {"point_atoms",std::move(point_atoms)},{"original_rows",t.original_rows}};
}
std::optional<JointTargetEvidence> ReadTarget(const j::value & value)
{
    if(value.is_null()) return {};
    const auto & o=value.as_object();
    RequireOnlyKeys(o,{"contract","checks","point_atoms","original_rows"},"unknown current target evidence field");
    JointTargetEvidence t;
    t.contract=Read<std::string>(o,"contract"); t.checks=ReadChecks(o.at("checks"));
    t.original_rows=Read<std::size_t>(o,"original_rows");
    for(const auto & item:o.at("point_atoms").as_array())
    {
        const auto & a=item.as_object(); JointAtomPointValidation point;
        RequireOnlyKeys(a,{"atom","component_convergence","parameter_validity","local_conditioning",
            "parameter_stability","point_quality","local_conditioning_blocks","parameter_stability_blocks","reason"},
            "unknown current target point field");
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
void ValidateTarget(const std::optional<JointTargetEvidence> & value,const std::vector<std::size_t> & expected_atoms,
    std::size_t rows)
{
    if(!value) return;
    const auto & t=*value;
    Require(t.contract==JointTargetPointValidationContractVersion,"target point validation contract");
    Require(t.original_rows==rows,"target validation row context");
    std::set<std::size_t> expected(expected_atoms.begin(),expected_atoms.end()),seen;
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
void Validate(const JointAnalysisResult & x)
{
    ValidateMetadata(x.metadata);
    const auto atoms=x.atom_ids.size(),rows=x.row_ids.size();
    Require(x.parameterization_contract==JointPureHaloSupportExclusionParameterizationContract,
        "current Joint result requires target-only Pure Exclusion");
    Require(x.selection_domain && x.selection_domain->IsValid(atoms) &&
        x.selection_domain->contributor_policy==PureHaloSupportExclusionContributorPolicy,
        "current Joint result requires Pure Exclusion selection domain");
    for(const double seconds:{x.costs.construction_seconds,x.costs.initialization_seconds,x.costs.search_seconds,
        x.costs.assessment_seconds,x.costs.assembly_seconds})
        Require(std::isfinite(seconds) && seconds>=0,"invalid Joint cost timing");
    const auto & selection=*x.selection_domain;
    const auto & provenance=*selection.pure_halo_support_exclusion;
    Require(atoms>0 && std::set<std::string>(x.atom_ids.begin(),x.atom_ids.end()).size()==atoms,
        "atom identities");
    Require(std::set<std::string>(x.row_ids.begin(),x.row_ids.end()).size()==rows,"row identities");
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
    const auto & targets=selection.target_indices;
    Require(x.fitting_to_original_contributors==targets,"Pure Exclusion fitting contributor mapping");
    Require(x.available_row_mask.size()==rows,"row mask dimensions");
    Require(std::isfinite(x.observation_scale) && x.observation_scale>=1,"observation scale");
    Require(!x.objective || std::isfinite(*x.objective),"nonfinite objective");
    Require(x.initialization.data_scope=="caller-provided-widths" ||
        x.initialization.data_scope=="contributor-local-sampling-may-read-outside-target-domain","initialization data scope");
    std::set<std::string> initialized;
    for(const auto & atom:x.initialization.atoms)
    {
        Require(std::find(x.atom_ids.begin(),x.atom_ids.end(),atom.id)!=x.atom_ids.end() &&
            initialized.insert(atom.id).second,"initialization atom identity");
        Require(!atom.original_b || std::isfinite(*atom.original_b),"nonfinite original initialization width provenance");
        Require(!atom.used_b || (std::isfinite(*atom.used_b) && *atom.used_b>0),
            "invalid used initialization width provenance");
    }
    if(x.initialization.valid)
    {
        Require(x.initialization.b.size()==targets.size(),"initial width dimensions");
        for(const double width:x.initialization.b)
            Require(std::isfinite(width) && width>0,"invalid initial width");
    }
    for(const auto & check:x.evidence)
    {
        Require(!check.name.empty(),"empty Joint evidence name");
        Require(!check.value || std::isfinite(*check.value),"nonfinite Joint evidence value");
        Require(!check.threshold || std::isfinite(*check.threshold),"nonfinite Joint evidence threshold");
    }
    const auto evidence=[&](const std::vector<JointCheck> & checks,const std::string & name) -> const JointCheck * {
        const auto found=std::find_if(checks.begin(),checks.end(),[&](const auto & check){return check.name==name;});
        return found==checks.end() ? nullptr : &*found;
    };
    std::optional<double> width_bound;
    if(x.metadata.solver)
    {
        const auto * calibration=evidence(x.evidence,"global-width-calibration");
        const auto * coverage=evidence(x.evidence,"qualified-first-stage-donor-coverage");
        Require(calibration && calibration->status==JointCheckStatus::Passed && calibration->value &&
            calibration->threshold && *calibration->value>0 &&
            *calibration->threshold==4.0*(*calibration->value) &&
            calibration->reason.find("global-first-stage-p95-k4")!=std::string::npos,
            "missing global width bound provenance");
        Require(coverage && coverage->status==JointCheckStatus::Passed && coverage->value &&
            coverage->threshold && *coverage->threshold==0.95 &&
            std::isfinite(*coverage->value) && *coverage->value>=0.95 && *coverage->value<=1.0,
            "missing qualified FirstStage donor coverage");
        width_bound=calibration->threshold;
    }
    else
        Require(!x.search_completed && x.components.empty() && !x.assembled_state && !x.objective,
            "solver result is missing current provenance");
    const auto validate_bound=[&](const JointState & state) {
        if(width_bound)
            for(double width:state.b)
                Require(std::log(width)<=std::log(*width_bound),"state exceeds global width bound");
    };
    if(x.assembled_state)
    {
        ValidateState(*x.assembled_state,targets.size());
        validate_bound(*x.assembled_state);
    }
    ValidateTarget(x.target_evidence,targets,rows);
    Require(!(x.runtime_convergence==JointCheckStatus::Passed) ||
        (x.initialization.valid && x.assembled_state && x.objective),
        "passed convergence evidence is missing its assembled state");
    const bool has_halos=atoms>targets.size();
    if(x.halo_median_reconstruction)
    {
        const auto & snapshot=*x.halo_median_reconstruction;
        Require(snapshot.valid_target_count<=targets.size(),"Halo median valid target count");
        Require(has_halos && snapshot.contract==HaloMedianReconstructionContractVersion &&
            snapshot.neighbor_count==GeometryLocalTargetMedianNeighborCount &&
            snapshot.donor_policy==GeometryLocalTargetMedianDonorPolicy,
            "invalid geometry-local Halo median contract");
        Require(!snapshot.joint_run_id.empty() && snapshot.valid_target_count>0,
            "invalid Halo median snapshot identity");
        for(double value:snapshot.global_median) Require(std::isfinite(value),"nonfinite Halo median global value");
        Require(snapshot.global_median[0]>=0 && snapshot.global_median[1]>0,
            "invalid Halo median global value");
        Require(snapshot.halos.size()==atoms-targets.size(),"Halo median coverage");
        std::size_t previous{}; bool first=true;
        for(const auto & halo:snapshot.halos)
        {
            Require(halo.halo_contributor_index<atoms && (first || previous<halo.halo_contributor_index) &&
                x.atom_ids[halo.halo_contributor_index]==halo.atom_id &&
                !std::binary_search(targets.begin(),targets.end(),halo.halo_contributor_index),
                "Halo median atom mapping");
            Require(halo.joint_run_id==snapshot.joint_run_id && halo.donor_count>0 &&
                halo.donor_count<=std::min(snapshot.neighbor_count,snapshot.valid_target_count) &&
                halo.donor_count==halo.donor_atom_ids.size(),"Halo median donor provenance");
            Require(halo.contract==HaloMedianReconstructionContractVersion &&
            halo.source==HaloMedianSource::GeometryLocalTargetMedian,
            "Halo median source contract");
            std::set<std::string> donors;
            for(const auto & donor_id:halo.donor_atom_ids)
                Require(!donor_id.empty() && donors.insert(donor_id).second &&
                    std::any_of(targets.begin(),targets.end(),[&](const auto index) {
                        return x.atom_ids[index]==donor_id;
                    }),"Halo median donor identity");
            Require(std::isfinite(halo.amplitude) && halo.amplitude>=0 &&
                std::isfinite(halo.width) && halo.width>0 && std::isfinite(halo.charge),
                "invalid Halo median parameters");
            previous=halo.halo_contributor_index; first=false;
        }
    }
    else if(has_halos)
        Require(!IsQualifiedJointAnalysisResult(x) || HasUnavailableJointHaloMedianReconstruction(x),
            "qualified Halo result lacks reconstruction or unavailable evidence");
    std::set<std::size_t> used_atoms,used_rows;
    std::set<std::string> components;
    for(const auto & c:x.components)
    {
        Require(!c.id.empty() && components.insert(c.id).second,"component identity");
        for(const auto atom:c.atoms)
            Require(atom<atoms && std::binary_search(targets.begin(),targets.end(),atom) && used_atoms.insert(atom).second,
                "component target mapping");
        for(const auto row:c.rows) Require(row<rows && used_rows.insert(row).second,"component row mapping");
        Require(std::isfinite(c.local_conditioning_seconds) && c.local_conditioning_seconds>=0,
            "local conditioning time");
        for(const auto & block:c.local_conditioning)
        {
            ValidateBlockConditioning(block);
            for(const auto atom:block.core_atoms)
                Require(std::binary_search(c.atoms.begin(),c.atoms.end(),atom),"local conditioning atom mapping");
        }
        if(c.parameter_stability) ValidateParameterStability(*c.parameter_stability);
        ValidateTarget(c.target_evidence,c.atoms,c.rows.size());
        if(c.state)
        {
            ValidateState(*c.state,c.atoms.size());
            validate_bound(*c.state);
        }
        if(c.candidate_diagnostic)
        {
            ValidateCandidateDiagnostic(*c.candidate_diagnostic,c.atoms.size(),c.state.has_value());
            if(width_bound && c.candidate_diagnostic->candidate_available &&
                c.candidate_diagnostic->candidate_feasible.value_or(false))
                for(const double eta:c.candidate_diagnostic->candidate_eta)
                    Require(eta<=std::log(*width_bound),"candidate exceeds global width bound");
        }
    }
}
}
std::string Encode(const JointAnalysisResult & x)
{
    Validate(x);
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
        {"profile_evaluations",c.profile_evaluations},
        {"accepted_updates",c.accepted_updates},{"native_status",c.native_status},{"state",OptionalState(c.state)},
        {"candidate_diagnostic",c.candidate_diagnostic ? j::value(CandidateDiagnostic(*c.candidate_diagnostic)) : j::value(nullptr)},
        {"evidence",Checks(c.evidence)},
        {"local_conditioning",BlockConditioning(c.local_conditioning)},
        {"local_conditioning_seconds",Number(c.local_conditioning_seconds)},
        {"local_conditioning_peak_memory_estimate_bytes",OptionalSize(c.local_conditioning_peak_memory_estimate_bytes)},
        {"parameter_stability",ParameterStability(c.parameter_stability)},
        {"runtime_convergence",StatusText(c.runtime_convergence)},
        {"target_evidence",Target(c.target_evidence)},
        {"target_runtime_convergence",StatusText(c.target_runtime_convergence)}};
        components.emplace_back(std::move(component));
    }
    Object out{{"schema_version",persistence_schema::JointResultJson},
        {"estimator","joint-components"},{"estimator_contract","guarded-joint-ls-v1"},
        {"objective_contract","parent-normalized-half-rss-v1"},{"support_contract","sphere-fma-v1"},
        {"parameterization_contract",x.parameterization_contract},{"metadata",Metadata(m)},
        {"selection_domain",SelectionDomain(x.selection_domain)},
        {"atom_ids",j::value_from(x.atom_ids)},{"row_ids",j::value_from(x.row_ids)},
        {"initialization",Initialization(x.initialization)},
        {"costs",Object{{"construction_seconds",t.construction_seconds},
            {"initialization_seconds",t.initialization_seconds},{"search_seconds",t.search_seconds},
            {"assessment_seconds",t.assessment_seconds},{"assembly_seconds",t.assembly_seconds}}},
        {"components",std::move(components)},{"assembled_state",OptionalState(x.assembled_state)},
        {"objective",Number(x.objective)},{"available_row_mask",j::value_from(x.available_row_mask)},
        {"evidence",Checks(x.evidence)},{"search_completed",x.search_completed},
        {"observation_scale",x.observation_scale},
        {"runtime_convergence",StatusText(x.runtime_convergence)},
        {"target_evidence",Target(x.target_evidence)},
        {"target_runtime_convergence",StatusText(x.target_runtime_convergence)},
        {"fitting_row_ids",j::value_from(x.fitting_row_ids)},
        {"fitting_to_original_contributors",j::value_from(x.fitting_to_original_contributors)},
        {"fitting_to_original_rows",j::value_from(x.fitting_to_original_rows)},
        {"halo_median_reconstruction",HaloMedianSnapshot(x.halo_median_reconstruction)}};
    return j::serialize(out);
}
JointAnalysisResult DecodeCurrent(const Object & o)
{
    RequireOnlyKeys(o,{"schema_version","estimator","estimator_contract","objective_contract","support_contract",
        "parameterization_contract","metadata","selection_domain","atom_ids","row_ids","fitting_row_ids",
        "fitting_to_original_contributors","fitting_to_original_rows","initialization","costs","components",
        "assembled_state","objective","available_row_mask","evidence","search_completed","observation_scale",
        "runtime_convergence","target_evidence","target_runtime_convergence",
        "halo_median_reconstruction"},"unknown current Joint result field");
    const auto & metadata=o.at("metadata").as_object();
    if(!metadata.at("solver").is_null())
    {
        RequireOnlyKeys(metadata.at("solver").as_object(),{"search_method","contract_version","sparse_backend",
            "fixed_neighbor_core_atoms","fixed_neighbor_policy_version","fixed_neighbor_maximum_sweeps",
            "fixed_neighbor_order","fixed_neighbor_local_search"},"unknown current solver provenance field");
    }
    RequireOnlyKeys(metadata,{"model_path","map_path","grid_size","grid_spacing","origin","simulation",
        "map_normalization","model_sha256","map_sha256","software","solver","units"},
        "unknown current Joint metadata field");
    RequireOnlyKeys(o.at("selection_domain").as_object(),{"contract","target_indices","observation_radius",
        "support_radius","contributor_policy","pure_halo_support_exclusion"},"unknown current selection field");
    RequireOnlyKeys(o.at("costs").as_object(),{"construction_seconds","initialization_seconds","search_seconds",
        "assessment_seconds","assembly_seconds"},"unknown current Joint cost field");
    JointAnalysisResult x;
    Require(Read<int>(o,"schema_version")==persistence_schema::JointResultJson,
        "unsupported current Joint result schema");
    Require(Read<std::string>(o,"estimator")=="joint-components" &&
        Read<std::string>(o,"estimator_contract")=="guarded-joint-ls-v1" &&
        Read<std::string>(o,"objective_contract")=="parent-normalized-half-rss-v1" &&
        Read<std::string>(o,"support_contract")=="sphere-fma-v1","unsupported estimator contract");
    x.parameterization_contract=Read<std::string>(o,"parameterization_contract");
    x.metadata=ReadMetadata(metadata);
    x.selection_domain=ReadSelectionDomain(o.at("selection_domain"));
    x.atom_ids=Read<std::vector<std::string>>(o,"atom_ids");
    x.row_ids=Read<std::vector<std::string>>(o,"row_ids");
    x.fitting_row_ids=Read<std::vector<std::string>>(o,"fitting_row_ids");
    x.fitting_to_original_contributors=Read<std::vector<std::size_t>>(o,"fitting_to_original_contributors");
    x.fitting_to_original_rows=Read<std::vector<std::size_t>>(o,"fitting_to_original_rows");
    x.initialization=ReadInitialization(o.at("initialization"));
    const auto & costs=o.at("costs").as_object();
    x.costs.construction_seconds=Read<double>(costs,"construction_seconds");
    x.costs.initialization_seconds=Read<double>(costs,"initialization_seconds");
    x.costs.search_seconds=Read<double>(costs,"search_seconds");
    x.costs.assessment_seconds=Read<double>(costs,"assessment_seconds");
    x.costs.assembly_seconds=Read<double>(costs,"assembly_seconds");
    for(const auto & item:o.at("components").as_array())
    {
        const auto & c=item.as_object();
        RequireOnlyKeys(c,{"id","stop_reason","atoms","rows","search_completed","search_sweeps",
            "first_order_stationarity_sweep","confirmed_stationarity_sweep","endpoint_certification_sweep",
            "profile_evaluations","accepted_updates","native_status","state","candidate_diagnostic","evidence",
            "local_conditioning","local_conditioning_seconds","local_conditioning_peak_memory_estimate_bytes",
            "parameter_stability","runtime_convergence","target_evidence",
            "target_runtime_convergence"},"unknown current Joint component field");
        JointAnalysisComponent component;
        component.id=Read<std::string>(c,"id");
        component.stop_reason=Read<std::string>(c,"stop_reason");
        component.atoms=Read<std::vector<std::size_t>>(c,"atoms");
        component.rows=Read<std::vector<std::size_t>>(c,"rows");
        component.search_completed=Read<bool>(c,"search_completed");
        component.search_sweeps=Read<std::size_t>(c,"search_sweeps");
        component.first_order_stationarity_sweep=Read<std::size_t>(c,"first_order_stationarity_sweep");
        component.confirmed_stationarity_sweep=Read<std::size_t>(c,"confirmed_stationarity_sweep");
        component.endpoint_certification_sweep=Read<std::size_t>(c,"endpoint_certification_sweep");
        component.profile_evaluations=Read<int>(c,"profile_evaluations");
        component.accepted_updates=Read<int>(c,"accepted_updates");
        component.native_status=Read<int>(c,"native_status");
        if(!c.at("state").is_null())
        {
            RequireOnlyKeys(c.at("state").as_object(),{"ac","b","log_b","width_gradient","objective"},
                "unknown current Joint state field");
                component.state=ReadState(c.at("state"));
        }
        if(!c.at("candidate_diagnostic").is_null())
        {
            RequireOnlyKeys(c.at("candidate_diagnostic").as_object(),{"candidate_status","candidate_available",
                "candidate_eta","candidate_ac","candidate_objective","candidate_ac_kkt",
                "candidate_width_gradient_inf_norm","convergence_criteria_id","global_ac_kkt_threshold",
                "global_width_gradient_threshold",
                "candidate_feasible","certification_status",
                "endpoint_verification_status","numerical_qualification_status","verification_reason",
                "resource_limit_reason"},"unknown current candidate diagnostic field");
            component.candidate_diagnostic=ReadCandidateDiagnostic(c.at("candidate_diagnostic"));
        }
        component.evidence=ReadChecks(c.at("evidence"));
        component.local_conditioning=ReadBlockConditioning(c.at("local_conditioning"));
        component.local_conditioning_seconds=Read<double>(c,"local_conditioning_seconds");
        component.local_conditioning_peak_memory_estimate_bytes=
            ReadOptionalSize(c.at("local_conditioning_peak_memory_estimate_bytes"));
        component.parameter_stability=ReadParameterStability(c.at("parameter_stability"));
        component.runtime_convergence=Status(c.at("runtime_convergence"));
        component.target_evidence=ReadTarget(c.at("target_evidence"));
        component.target_runtime_convergence=Status(c.at("target_runtime_convergence"));
        x.components.push_back(std::move(component));
    }
    if(!o.at("assembled_state").is_null())
    {
        RequireOnlyKeys(o.at("assembled_state").as_object(),{"ac","b","log_b","width_gradient","objective"},
            "unknown current assembled state field");
        x.assembled_state=ReadState(o.at("assembled_state"));
    }
    x.objective=OptionalNumber(o.at("objective"));
    x.available_row_mask=Read<std::vector<bool>>(o,"available_row_mask");
    x.evidence=ReadChecks(o.at("evidence"));
    x.search_completed=Read<bool>(o,"search_completed");
    x.observation_scale=Read<double>(o,"observation_scale");
    x.runtime_convergence=Status(o.at("runtime_convergence"));
    x.target_evidence=ReadTarget(o.at("target_evidence"));
    x.target_runtime_convergence=Status(o.at("target_runtime_convergence"));
    x.halo_median_reconstruction=ReadHaloMedianSnapshot(o.at("halo_median_reconstruction"));
    Validate(x);
    return x;
}
JointAnalysisResult Decode(std::string_view text)
{
    j::parse_options options; options.numbers=j::number_precision::precise;
    const auto parsed=j::parse(text,{},options); const auto & o=parsed.as_object();
    const int version=Read<int>(o,"schema_version");
    Require(version==persistence_schema::JointResultJson,"unsupported Joint result schema version");
    return DecodeCurrent(o);
}
}
