#include "StageResultJson.hpp"
#include "PersistenceSchema.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/detail/PostJointLocalValidation.hpp"
#include "data/detail/PostJointQScoreValidation.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <boost/json.hpp>
#include <cmath>
#include <algorithm>
#include <limits>
#include <optional>
#include <set>

namespace rhbm_gem::stage_result_io {
namespace {
namespace j = boost::json;
using O = j::object;
using V = j::value;
using A = j::array;
template<class T> T Read(const O & o, const char * key) { return j::value_to<T>(o.at(key)); }
double Number(const V & v) { const double x = j::value_to<double>(v); if (!std::isfinite(x)) throw std::invalid_argument("Nonfinite stage value."); return x; }
template<class E> E Enum(const V & v, int maximum) { const int x = j::value_to<int>(v); if (x < 0 || x > maximum) throw std::invalid_argument("Unknown stage enum."); return static_cast<E>(x); }
// Failed native fits can contain nonfinite diagnostics; final stage points cannot.
V NativeNumber(double x) { return std::isfinite(x) ? V(x) : V(std::isnan(x) ? "nan" : x>0 ? "+inf":"-inf"); }
double NativeNumber(const V & v)
{
    if(!v.is_string()) return Number(v);
    if(v.as_string()=="nan") return std::numeric_limits<double>::quiet_NaN();
    if(v.as_string()=="+inf") return std::numeric_limits<double>::infinity();
    if(v.as_string()=="-inf") return -std::numeric_limits<double>::infinity();
    throw std::invalid_argument("Invalid native diagnostic number.");
}
V Pack(const GaussianModel3D & p) { return A{NativeNumber(p.GetAmplitude()),NativeNumber(p.GetWidth()),NativeNumber(p.GetOffset())}; }
GaussianModel3D Point(const V & v, bool positive_width=true)
{
    const auto & a=v.as_array(); if(a.size()!=3) throw std::invalid_argument("Invalid stage point.");
    const auto read=positive_width ? Number : static_cast<double(*)(const V &)>(NativeNumber);
    GaussianModel3D p(read(a[0]),read(a[1]),read(a[2]));
    if(positive_width) GaussianModel3D::RequireFinitePositiveWidthModel(p);
    return p;
}
template<class T> V Matrix(const T & m) { A a; for (Eigen::Index r=0;r<m.rows();++r) for (Eigen::Index c=0;c<m.cols();++c) a.push_back(m(r,c)); return a; }
template<int R,int C> Eigen::Matrix<double,R,C> Matrix(const V & v) { const auto & a=v.as_array(); if(a.size()!=R*C) throw std::invalid_argument("Stage matrix shape mismatch."); Eigen::Matrix<double,R,C> m; std::size_t i=0; for(int r=0;r<R;++r) for(int c=0;c<C;++c) m(r,c)=Number(a[i++]); return m; }
V Pack(const EstimateSource & s)
{
    O out{{"method",static_cast<int>(s.method)},{"atom",s.atom_id},{"component",s.component_id},
        {"run",s.run_id},{"role",static_cast<int>(s.role)}};
    if(s.fixed_from_initial_stage) out.emplace("fixed_from_initial_stage",true);
    return out;
}
EstimateSource Source(const V & v)
{
    const auto & o=v.as_object();
    const bool has_fixed_source=o.contains("fixed_from_initial_stage");
    if(o.size()!=(has_fixed_source ? 6u : 5u)) throw std::invalid_argument("Invalid stage estimate source fields.");
    const auto method=Enum<EstimateMethod>(o.at("method"),static_cast<int>(EstimateMethod::JointComponents));
    if (method!=EstimateMethod::Unspecified && method!=EstimateMethod::LocalMDPDE &&
        method!=EstimateMethod::JointComponents)
        throw std::invalid_argument("Unknown estimate method.");
    return {method,Read<std::string>(o,"atom"),Read<std::string>(o,"component"),Read<std::string>(o,"run"),
        Enum<FittingRole>(o.at("role"),2),has_fixed_source && Read<bool>(o,"fixed_from_initial_stage")};
}
V Pack(const PostJointQScoreResult & q)
{
    A halos;
    for (const auto & halo:q.extra_halos)
        halos.push_back(O{{"atom",halo.atom_id},{"donors",j::value_from(halo.donor_atom_ids)}});
    return O{{"source",Pack(q.source)},{"contract",q.contract},{"status",static_cast<int>(q.status)},
        {"scores",q.scores ? V(A{q.scores->standard,q.scores->joint}):V{}},{"reason",q.reason},
        {"reference",A{q.reference_height,q.reference_offset}},{"samples",q.sample_count},
        {"contributors",q.contributor_count},{"missing_contributors",q.missing_contributor_count},
        {"contributions",q.evaluated_contributions},{"halos",q.halo_count},{"extra_halos",std::move(halos)}};
}
PostJointQScoreResult QScore(const V & value)
{
    const auto & o=value.as_object();
    if (o.size()!=12) throw std::invalid_argument("Invalid QScore result fields.");
    PostJointQScoreResult q;
    q.source=Source(o.at("source")); q.contract=Read<std::string>(o,"contract");
    q.status=Enum<EvidenceStatus>(o.at("status"),3); q.reason=Read<std::string>(o,"reason");
    if (!o.at("scores").is_null())
    {
        const auto & scores=o.at("scores").as_array();
        if (scores.size()!=2) throw std::invalid_argument("Invalid paired QScore shape.");
        q.scores=PostJointQScoreValues{Number(scores[0]),Number(scores[1])};
    }
    const auto & reference=o.at("reference").as_array();
    if (reference.size()!=2) throw std::invalid_argument("Invalid QScore reference shape.");
    q.reference_height=Number(reference[0]); q.reference_offset=Number(reference[1]);
    q.sample_count=Read<std::size_t>(o,"samples"); q.contributor_count=Read<std::size_t>(o,"contributors");
    q.missing_contributor_count=Read<std::size_t>(o,"missing_contributors");
    q.evaluated_contributions=Read<std::size_t>(o,"contributions"); q.halo_count=Read<std::size_t>(o,"halos");
    for (const auto & halo:o.at("extra_halos").as_array())
    {
        const auto & h=halo.as_object();
        if (h.size()!=2) throw std::invalid_argument("Invalid QScore Halo fields.");
        q.extra_halos.push_back({Read<int>(h,"atom"),Read<std::vector<std::string>>(h,"donors")});
    }
    return q;
}
V Pack(const LocalStageEstimate & s)
{
    const auto & u=s.uncertainty;
    O uncertainty{{"status",static_cast<int>(u.status)},{"method",u.method},{"reason",u.reason},{"rank",u.rank},{"df",u.degrees_of_freedom},{"threshold",u.rank_threshold},
        {"covariance",u.covariance ? Matrix(*u.covariance) : V{}},{"variance",u.residual_variance ? V(*u.residual_variance) : V{}}};
    return O{{"point",s.point ? Pack(*s.point) : V{}},{"source",Pack(s.source)},{"reason",s.reason},{"convergence",static_cast<int>(s.convergence)},{"uncertainty",uncertainty}};
}
LocalStageEstimate Stage(const V & v)
{
    const auto & o=v.as_object(); LocalStageEstimate s;
    if(!o.at("point").is_null()) s.point=Point(o.at("point"));
    s.source=Source(o.at("source")); s.reason=Read<std::string>(o,"reason"); s.convergence=Enum<JointCheckStatus>(o.at("convergence"),3);
    const auto & u=o.at("uncertainty").as_object(); auto & t=s.uncertainty;
    t.status=Enum<EvidenceStatus>(u.at("status"),3); t.method=Read<std::string>(u,"method"); t.reason=Read<std::string>(u,"reason");
    t.rank=Read<std::size_t>(u,"rank"); t.degrees_of_freedom=Read<std::size_t>(u,"df"); t.rank_threshold=Number(u.at("threshold"));
    if(!u.at("covariance").is_null()) t.covariance=Matrix<3,3>(u.at("covariance"));
    if(!u.at("variance").is_null()) t.residual_variance=Number(u.at("variance"));
    return s;
}
V Pack(const GaussianModel3DWithUncertainty & p)
{
    const auto & s=p.GetStandardDeviationModel();
    return O{{"point",Pack(p.GetModel())},{"sd",A{NativeNumber(s.GetAmplitude()),NativeNumber(s.GetWidth()),NativeNumber(s.GetOffset())}}};
}
V Pack(const JointFixedChargeProfile & p) { return A{NativeNumber(p.coefficient), NativeNumber(p.width)}; }
JointFixedChargeProfile FixedChargeProfile(const V & v)
{
    const auto & a = v.as_array();
    if (a.size() != 2) throw std::invalid_argument("Invalid fixed charge profile.");
    return {Number(a[0]), Number(a[1])};
}
GaussianModel3DWithUncertainty Uncertain(const V & v)
{
    const auto & o=v.as_object(); const auto & a=o.at("sd").as_array();
    if(a.size()!=3) throw std::invalid_argument("Invalid uncertainty shape.");
    return {Point(o.at("point"),false),GaussianModel3DUncertainty(NativeNumber(a[0]),NativeNumber(a[1]),a[2].is_null() ? std::numeric_limits<double>::quiet_NaN():NativeNumber(a[2]))};
}
V OptionalNumber(const std::optional<double> & v) { return v ? NativeNumber(*v) : V{}; }
std::optional<double> OptionalNumber(const V & v) { return v.is_null() ? std::nullopt : std::optional<double>{NativeNumber(v)}; }
V Pack(const RHBMEndpointRefinementDiagnostics & r)
{
    V coordinates;
    if(r.relative_coordinate_difference)
        coordinates=A{NativeNumber((*r.relative_coordinate_difference)[0]),NativeNumber((*r.relative_coordinate_difference)[1]),
            NativeNumber((*r.relative_coordinate_difference)[2])};
    return O{{"accepted",r.accepted},{"reason",r.reason},{"evaluations",r.candidate_equation_evaluations},
        {"updates",r.reference_updates},{"stop",r.reference_stop},{"original",OptionalNumber(r.original_residual)},
        {"candidate",OptionalNumber(r.candidate_residual)},{"reference",OptionalNumber(r.reference_residual)},
        {"coordinate_difference",std::move(coordinates)},{"weight_difference",OptionalNumber(r.weight_max_difference)},
        {"floor_masks_equal",r.floor_masks_equal ? V(*r.floor_masks_equal):V{}}};
}
RHBMEndpointRefinementDiagnostics ReadEndpointRefinement(const V & value)
{
    const auto & r=value.as_object(); RHBMEndpointRefinementDiagnostics x;
    x.accepted=Read<bool>(r,"accepted"); x.reason=Read<std::string>(r,"reason");
    x.candidate_equation_evaluations=Read<int>(r,"evaluations"); x.reference_updates=Read<int>(r,"updates");
    x.reference_stop=Read<std::string>(r,"stop"); x.original_residual=OptionalNumber(r.at("original"));
    x.candidate_residual=OptionalNumber(r.at("candidate")); x.reference_residual=OptionalNumber(r.at("reference"));
    if(!r.at("coordinate_difference").is_null())
    {
        const auto & a=r.at("coordinate_difference").as_array();
        if(a.size()!=3) throw std::invalid_argument("Invalid refinement coordinate shape.");
        x.relative_coordinate_difference=std::array<double,3>{NativeNumber(a[0]),NativeNumber(a[1]),NativeNumber(a[2])};
    }
    x.weight_max_difference=OptionalNumber(r.at("weight_difference"));
    if(!r.at("floor_masks_equal").is_null()) x.floor_masks_equal=Read<bool>(r,"floor_masks_equal");
    return x;
}
V Pack(const RHBMFailedOnlyRecoveryDiagnostics & r)
{
    V endpoint, retry;
    if(r.endpoint_refinement) endpoint=Pack(*r.endpoint_refinement);
    if(r.bounded_native_retry)
    {
        const auto & x=*r.bounded_native_retry;
        retry=O{{"max_iterations",x.max_iterations},{"native_status",x.native_status ? V(*x.native_status):V{}},
            {"iterations",x.iterations},{"squared_beta_change",OptionalNumber(x.squared_beta_change)},
            {"relative_variance_change",OptionalNumber(x.relative_variance_change)},
            {"equation_residual",OptionalNumber(x.equation_residual)},{"seconds",NativeNumber(x.seconds)}};
    }
    return O{{"accepted",r.accepted},{"reason",r.reason},{"endpoint_refinement_seconds",NativeNumber(r.endpoint_refinement_seconds)},
        {"endpoint_refinement",std::move(endpoint)},{"bounded_native_retry",std::move(retry)}};
}
RHBMFailedOnlyRecoveryDiagnostics ReadFailedOnlyRecovery(const V & value)
{
    const auto & r=value.as_object(); RHBMFailedOnlyRecoveryDiagnostics x;
    x.accepted=Read<bool>(r,"accepted"); x.reason=Read<std::string>(r,"reason");
    x.endpoint_refinement_seconds=NativeNumber(r.at("endpoint_refinement_seconds"));
    if(!r.at("endpoint_refinement").is_null()) x.endpoint_refinement=ReadEndpointRefinement(r.at("endpoint_refinement"));
    if(!r.at("bounded_native_retry").is_null())
    {
        const auto & item=r.at("bounded_native_retry").as_object(); RHBMNativeRetryDiagnostics retry;
        retry.max_iterations=Read<int>(item,"max_iterations");
        if(!item.at("native_status").is_null()) retry.native_status=Read<int>(item,"native_status");
        retry.iterations=Read<int>(item,"iterations"); retry.squared_beta_change=OptionalNumber(item.at("squared_beta_change"));
        retry.relative_variance_change=OptionalNumber(item.at("relative_variance_change"));
        retry.equation_residual=OptionalNumber(item.at("equation_residual")); retry.seconds=NativeNumber(item.at("seconds"));
        x.bounded_native_retry=std::move(retry);
    }
    return x;
}
V Pack(const LocalFitDiagnostics & d)
{
    V refinement;
    if (d.refinement)
    {
        const auto & r=*d.refinement;
        refinement=O{{"accepted",r.accepted},{"reason",r.reason},{"evaluations",r.candidate_equation_evaluations},
            {"updates",r.reference_updates},{"stop",r.reference_stop},{"original",OptionalNumber(r.original_residual)},
            {"candidate",OptionalNumber(r.candidate_residual)},{"reference",OptionalNumber(r.reference_residual)},
            {"coordinate_difference",r.relative_coordinate_difference ? V(A{NativeNumber((*r.relative_coordinate_difference)[0]),NativeNumber((*r.relative_coordinate_difference)[1]),NativeNumber((*r.relative_coordinate_difference)[2])}):V{}},
            {"weight_difference",OptionalNumber(r.weight_max_difference)},{"floor_masks_equal",r.floor_masks_equal ? V(*r.floor_masks_equal):V{}}};
    }
    V failed_only_recovery;
    if(d.failed_only_recovery) failed_only_recovery=Pack(*d.failed_only_recovery);
    return O{{"status",static_cast<int>(d.status)},{"variance",NativeNumber(d.sigma_square)},{"iterations",d.iterations.iterations},
        {"beta_change",OptionalNumber(d.iterations.squared_beta_change)},{"variance_change",OptionalNumber(d.iterations.relative_variance_change)},
        {"refinement",refinement},{"failed_only_recovery",std::move(failed_only_recovery)}};
}
LocalFitDiagnostics Diagnostics(const V & v)
{
    const auto & o=v.as_object(); LocalFitDiagnostics d;
    d.status=Enum<RHBMEstimationStatus>(o.at("status"),4); d.sigma_square=NativeNumber(o.at("variance"));
    d.iterations={Read<int>(o,"iterations"),OptionalNumber(o.at("beta_change")),OptionalNumber(o.at("variance_change"))};
    if (!o.at("refinement").is_null())
    {
        d.refinement=ReadEndpointRefinement(o.at("refinement"));
    }
    if(o.contains("failed_only_recovery") && !o.at("failed_only_recovery").is_null())
        d.failed_only_recovery=ReadFailedOnlyRecovery(o.at("failed_only_recovery"));
    return d;
}
V Pack(const PostJointLocalResult & result)
{
    V information;
    if (result.information)
        information=O{{"precision",Matrix(result.information->precision)},
            {"moment",Matrix(result.information->moment)}};
    A failures;
    for (const auto & failure:result.sample_failures)
        failures.push_back(O{{"original_row",failure.original_row_index},{"grid",failure.global_grid_index},{"distance",failure.distance},{"reason",failure.reason}});
    return O{{"source",Pack(result.source)}, {"status",static_cast<int>(result.status)}, {"reason",result.reason},
        {"alpha",result.fit.alpha_r}, {"ols",Pack(result.fit.ols)}, {"mdpde",Pack(result.fit.mdpde)},
        {"uncertainty_recorded",result.fit.uncertainty_recorded},
        {"diagnostics",result.fit.diagnostics ? Pack(*result.fit.diagnostics):V{}},
        {"fixed_charge",A{result.fixed_charge.coefficient,result.fixed_charge.width}},
        {"information",information}, {"samples",result.sample_count},
        {"effective_samples",result.effective_sample_count}, {"corrected_samples",result.corrected_sample_count},
        {"positive_samples",result.positive_sample_count}, {"fitting_radius",result.fitting_radius}, {"sample_failures",failures}};
}
PostJointLocalResult PostJointLocal(const V & value)
{
    const auto & o=value.as_object();
    PostJointLocalResult result;
    result.source=Source(o.at("source"));
    result.status=Enum<EvidenceStatus>(o.at("status"),3);
    result.reason=Read<std::string>(o,"reason");
    result.fit.alpha_r=Number(o.at("alpha"));
    result.fit.ols=Uncertain(o.at("ols"));
    result.fit.mdpde=Uncertain(o.at("mdpde"));
    result.fit.uncertainty_recorded=Read<bool>(o,"uncertainty_recorded");
    if(!o.at("diagnostics").is_null()) result.fit.diagnostics=Diagnostics(o.at("diagnostics"));
    const auto & charge=o.at("fixed_charge").as_array();
    if(charge.size()!=2) throw std::invalid_argument("Invalid fixed charge profile.");
    result.fixed_charge={Number(charge[0]),Number(charge[1])};
    if(!o.at("information").is_null())
    {
        const auto & information=o.at("information").as_object();
        RHBMInformation parsed;
        parsed.precision=Matrix<2,2>(information.at("precision"));
        parsed.moment=Matrix<2,1>(information.at("moment"));
        result.information=std::move(parsed);
    }
    result.sample_count=Read<std::size_t>(o,"samples");
    result.effective_sample_count=Read<std::size_t>(o,"effective_samples");
    result.corrected_sample_count=Read<std::size_t>(o,"corrected_samples");
    result.positive_sample_count=Read<std::size_t>(o,"positive_samples");
    result.fitting_radius=Number(o.at("fitting_radius"));
    for(const auto & item:o.at("sample_failures").as_array())
    {
        const auto & failure=item.as_object();
        result.sample_failures.push_back({Read<std::size_t>(failure,"original_row"),Number(failure.at("distance")),
            Read<std::string>(failure,"reason"),Read<std::size_t>(failure,"grid")});
    }
    return result;
}
V InitialNative(const LocalPotentialEntry & e)
{
    const auto r=e.InitialStageGaussianResult(false); const auto & sd=r.mdpde.GetStandardDeviationModel();
    return O{{"alpha",r.alpha_r},{"ols",Pack(r.ols.GetModel())},{"ols_sd",r.uncertainty_recorded ? Pack(r.ols).as_object().at("sd"):V{}},
        {"mdpde_sd",r.uncertainty_recorded ? V(A{NativeNumber(sd.GetAmplitude()),NativeNumber(sd.GetWidth()),NativeNumber(sd.GetOffset())}):V{}},
        {"unfitted",e.InitialStageEstimate().point ? V{} : Pack(r.mdpde.GetModel())},
        {"diagnostics",r.diagnostics ? Pack(*r.diagnostics):V{}}};
}
void RestoreInitialNative(LocalPotentialEntry & e, const V & native, LocalStageEstimate estimate)
{
    const auto & o=native.as_object(); LocalGaussianResult r; r.alpha_r=Number(o.at("alpha"));
    r.uncertainty_recorded=!o.at("mdpde_sd").is_null();
    if(r.uncertainty_recorded!=!o.at("ols_sd").is_null()) throw std::invalid_argument("Incomplete native uncertainty.");
    r.ols={Point(o.at("ols"),false),{}};
    r.mdpde={estimate.point ? *estimate.point : Point(o.at("unfitted"),false),{}};
    if(r.uncertainty_recorded)
    {
        r.ols=Uncertain(O{{"point",o.at("ols")},{"sd",o.at("ols_sd")}});
        r.mdpde=Uncertain(O{{"point",Pack(r.mdpde.GetModel())},{"sd",o.at("mdpde_sd")}});
    }
    if(!o.at("diagnostics").is_null()) r.diagnostics=Diagnostics(o.at("diagnostics"));
    e.SetInitialStageGaussianResult(std::move(r));
    e.SetInitialStageEstimate(std::move(estimate));
}
V Pack(const GroupGaussianMemberResult & m) { return O{{"posterior",Pack(m.posterior)},{"outlier",m.is_outlier},{"distance",m.statistical_distance},{"source",m.evidence_source_id},{"charge_inferred",m.charge_inferred},{"covariance",m.parameter_covariance ? Matrix(*m.parameter_covariance):V{}},{"fixed_charge",m.fixed_charge ? Pack(*m.fixed_charge):V{}}}; }
GroupGaussianMemberResult Member(const V & v)
{
    const auto & o=v.as_object(); GroupGaussianMemberResult m;
    m.posterior=Uncertain(o.at("posterior")); m.is_outlier=Read<bool>(o,"outlier");
    m.statistical_distance=Number(o.at("distance")); m.evidence_source_id=Read<std::string>(o,"source");
    m.charge_inferred=Read<bool>(o,"charge_inferred");
    if(!o.at("covariance").is_null()) m.parameter_covariance=Matrix<2,2>(o.at("covariance"));
    if(!o.at("fixed_charge").is_null()) m.fixed_charge=FixedChargeProfile(o.at("fixed_charge"));
    return m;
}
V Pack(const GroupParameterSummary & s)
{
    V inference;
    if(s.inference) { const auto & i=*s.inference; inference=O{{"alpha",i.alpha_g},{"mean",Pack(i.mean)},{"mdpde",Pack(i.mdpde)},{"prior",Pack(i.prior)}}; }
    return O{{"status",static_cast<int>(s.status)},{"reason",s.reason},{"source",s.source_id},{"correlation",s.correlation_policy},{"points",s.point_count},{"eligible",s.eligible_count},{"excluded",s.excluded_count},{"mean",s.descriptive_mean ? Pack(*s.descriptive_mean):V{}},{"charge_median",OptionalNumber(s.descriptive_charge_median)},{"inference",inference},{"ids",j::value_from(s.member_ids)}};
}
GroupParameterSummary Summary(const V & v)
{
    const auto & o=v.as_object(); GroupParameterSummary s;
    s.status=Enum<EvidenceStatus>(o.at("status"),3);
    s.reason=Read<std::string>(o,"reason");
    s.source_id=Read<std::string>(o,"source");
    s.correlation_policy=Read<std::string>(o,"correlation");
    s.point_count=Read<std::size_t>(o,"points");
    s.eligible_count=Read<std::size_t>(o,"eligible");
    s.excluded_count=Read<std::size_t>(o,"excluded");
    s.member_ids=Read<std::vector<int>>(o,"ids");
    s.descriptive_charge_median=OptionalNumber(o.at("charge_median"));
    if(!o.at("mean").is_null()) s.descriptive_mean=Point(o.at("mean"));
    if(!o.at("inference").is_null())
    {
        const auto & i=o.at("inference").as_object(); GroupGaussianResult r;
        r.alpha_g=Number(i.at("alpha")); r.mean=Point(i.at("mean"));
        r.mdpde=Point(i.at("mdpde")); r.prior=Uncertain(i.at("prior"));
        s.inference=std::move(r);
    }
    return s;
}
V Samples(const LocalPotentialSampleList & samples, bool geometry)
{
    A out; for(const auto & s:samples) out.push_back(O{{"y",s.response},{"distance",s.point.distance},{"selected",s.point.is_selected},{"position",geometry ? j::value_from(s.point.position):V{}}}); return out;
}
LocalPotentialSampleList Samples(const V & value, bool geometry)
{
    LocalPotentialSampleList out; for(const auto & v:value.as_array()) { const auto & o=v.as_object(); LocalPotentialSample s; s.response=Number(o.at("y")); s.point.distance=Number(o.at("distance")); s.point.is_selected=Read<bool>(o,"selected"); if(geometry) {s.point.position=Read<std::array<double,3>>(o,"position"); for(double x:s.point.position) if(!std::isfinite(x)) throw std::invalid_argument("Invalid sample position.");} else if(!o.at("position").is_null()) throw std::invalid_argument("Unexpected sample geometry."); out.push_back(s); } return out;
}
}
std::string Encode(const ModelObject & model)
{
    Validate(model);
    const auto & data=ModelAnalysisData::Of(model); A atoms,groups;
    for(const auto & atom:model.GetAtomList())
    {
        const auto * e=data.FindAtomLocalEntry(*atom); if(!e) continue;
        V peeling;
        if(e->PostFitPeeling())
        {
            const auto & p=*e->PostFitPeeling();
            A samples;
            for(const auto & sample:p.samples)
            {
                const auto & observation=sample.observation;
                samples.push_back(O{{"response",sample.response ? V(*sample.response):V{}},{"reason",sample.reason},
                    {"original_row",observation.original_row_index},{"grid",observation.global_grid_index},
                    {"squared_distance",observation.squared_distance},{"raw_response",observation.raw_response},
                    {"position",j::value_from(observation.position)}});
            }
            peeling=O{{"source",Pack(p.source)},{"mode",p.mode},{"neighbors",p.neighbor_count},{"samples",samples}};
        }
        atoms.push_back(O{{"id",atom->GetSerialID()},{"initial",Pack(e->InitialStageEstimate())},{"joint",Pack(e->JointStageEstimate())},
            {"initial_native",InitialNative(*e)},
            {"geometry",e->SampleGeometryAvailable()},{"raw",Samples(e->RawSamplingEntries(),e->SampleGeometryAvailable())},
            {"post_joint_qscore",e->PostJointQScore() ? Pack(*e->PostJointQScore()):V{}},
            {"peeling",peeling},{"post_joint_local",e->PostJointLocalResult() ? Pack(*e->PostJointLocalResult()):V{}},
            {"posterior",e->GroupMemberResult() ? Pack(*e->GroupMemberResult()):V{}}});
    }
    const auto & g=data.AtomGroupEntry();
    for(const auto key:g.CollectGroupKeys()) { A ids; for(const auto * a:g.GetMembers(key)) ids.push_back(a->GetSerialID()); groups.push_back(O{{"key",static_cast<std::uint64_t>(key)},{"members",ids},{"alpha",g.GetParameterSummary(key) && !g.GetParameterSummary(key)->inference ? V(g.GetAlphaG(key)):V{}},{"native",g.GetParameterSummary(key) ? V{}:V(O{{"alpha",g.GetAlphaG(key)},{"mean",Pack(g.GetMean(key))},{"mdpde",Pack(g.GetMDPDE(key))},{"prior",Pack(g.GetPriorWithUncertainty(key))}})},{"summary",g.GetParameterSummary(key) ? Pack(*g.GetParameterSummary(key)):V{}}}); }
    return j::serialize(O{{"version",persistence_schema::StageResultJson},{"atoms",atoms},{"groups",groups}});
}
void Decode(ModelObject & model, const std::string & text)
{
    j::parse_options parse_options; parse_options.numbers=j::number_precision::precise;
    const auto document=j::parse(text,{},parse_options); const auto & root=document.as_object();
    if(Read<int>(root,"version")!=persistence_schema::StageResultJson)
        throw std::invalid_argument("Unsupported Stage result schema version.");
    auto & data=ModelAnalysisData::Of(model); std::set<int> seen;
    for(const auto & value:root.at("atoms").as_array())
    {
        const auto & o=value.as_object(); const int id=Read<int>(o,"id"); auto * atom=model.FindAtomPtr(id);
        if(!atom || !seen.insert(id).second) throw std::invalid_argument("Invalid stage atom identity.");
        if(o.size()!=10 || !o.contains("initial") || !o.contains("joint") || !o.contains("initial_native"))
            throw std::invalid_argument("Invalid Stage atom fields.");
        auto & entry=data.EnsureAtomLocalEntry(*atom);
        RestoreInitialNative(entry,o.at("initial_native"),Stage(o.at("initial")));
        entry.SetJointStageEstimate(Stage(o.at("joint")));
        const bool geometry=Read<bool>(o,"geometry"); entry.SetSampleGeometryAvailable(geometry);
        entry.SetRawSamplingEntries(Samples(o.at("raw"),geometry));
        if(!o.at("peeling").is_null())
        {
            const auto & p=o.at("peeling").as_object(); PostFitPeelingResult result;
            result.source=Source(p.at("source")); result.mode=Read<std::string>(p,"mode");
            result.neighbor_count=Read<std::size_t>(p,"neighbors");
            for(const auto & sample:p.at("samples").as_array())
            {
                const auto & item=sample.as_object(); PeelingSampleEstimate estimate;
                estimate.reason=Read<std::string>(item,"reason");
                estimate.observation={Read<std::size_t>(item,"original_row"),Read<std::size_t>(item,"grid"),
                    Number(item.at("squared_distance")),Number(item.at("raw_response")),
                    Read<std::array<double,3>>(item,"position")};
                if(!item.at("response").is_null()) estimate.response=Number(item.at("response"));
                result.samples.push_back(std::move(estimate));
            }
            entry.SetPostFitPeeling(std::move(result));
        }
        if(!o.at("post_joint_local").is_null())
            entry.SetPostJointLocalResult(PostJointLocal(o.at("post_joint_local")));
        entry.ClearPostJointQScore();
        if(!o.at("post_joint_qscore").is_null()) entry.SetPostJointQScore(QScore(o.at("post_joint_qscore")));
        entry.ClearGroupMemberResult();
        if(!o.at("posterior").is_null()) entry.SetGroupMemberResult(Member(o.at("posterior")));
    }
    if(seen.size()!=data.AtomLocalEntries().size())
        throw std::invalid_argument("Incomplete Stage local records.");

    auto & groups=data.AtomGroupEntry(); groups.ClearMembers(); std::set<GroupKey> keys;
    for(const auto & value:root.at("groups").as_array())
    {
        const auto & o=value.as_object(); const auto key=Read<GroupKey>(o,"key");
        if(!keys.insert(key).second) throw std::invalid_argument("Duplicate Stage group.");
        std::set<int> members;
        for(const int id:Read<std::vector<int>>(o,"members"))
        {
            auto * atom=model.FindAtomPtr(id);
            if(!atom || !members.insert(id).second)
                throw std::invalid_argument("Invalid Stage group member identity.");
            groups.AddMember(key,*atom);
        }
        if(!o.at("alpha").is_null()) groups.SetAlphaG(key,Number(o.at("alpha")));
        if(!o.at("native").is_null())
        {
            const auto & n=o.at("native").as_object(); GroupGaussianResult result;
            result.alpha_g=Number(n.at("alpha")); result.mean=Point(n.at("mean"),false);
            result.mdpde=Point(n.at("mdpde"),false); result.prior=Uncertain(n.at("prior"));
            groups.SetGaussianResult(key,result);
        }
        if(!o.at("summary").is_null())
        {
            auto summary=Summary(o.at("summary"));
            for(const int id:summary.member_ids)
                if(!members.contains(id)) throw std::invalid_argument("Posterior outside fitted group.");
            groups.SetParameterSummary(key,std::move(summary));
        }
    }
    Validate(model);
}
void Validate(const ModelObject & model)
{
    const auto & data=ModelAnalysisData::Of(model);
    const auto * joint_snapshot=data.joint_result ? &*data.joint_result : nullptr;
    const bool unqualified_joint=joint_snapshot && !IsQualifiedJointAnalysisResult(*joint_snapshot);
    std::optional<std::string> joint_run_id;
    std::optional<data_internal::PostJointQScoreValidation> qscore_validation;
    if (std::any_of(data.AtomLocalEntries().begin(),data.AtomLocalEntries().end(),
            [](const auto & item){return item.second->PostJointQScore().has_value();}))
        qscore_validation.emplace(model);
    for(const auto & [id,e]:data.AtomLocalEntries())
    {
        const auto & initial=e->InitialStageEstimate();
        const auto & s=e->JointStageEstimate();
        if ((initial.source.method!=EstimateMethod::Unspecified && initial.source.method!=EstimateMethod::LocalMDPDE) ||
            (s.source.method!=EstimateMethod::JointComponents &&
                (s.source.method!=EstimateMethod::Unspecified || s.point)))
            throw std::invalid_argument("Estimate method does not match the fitting stage.");
        if (s.source.method == EstimateMethod::JointComponents && !s.source.run_id.empty())
        {
            if (!joint_run_id) joint_run_id=s.source.run_id;
            else if (*joint_run_id != s.source.run_id)
                throw std::invalid_argument("Mixed Joint run identities.");
            if (joint_snapshot && joint_snapshot->halo_median_reconstruction &&
                s.source.run_id != joint_snapshot->halo_median_reconstruction->joint_run_id)
                throw std::invalid_argument("Joint stage and Halo median run identities differ.");
        }
        const bool formal_joint_point = s.point && s.convergence == JointCheckStatus::Passed;
        if (unqualified_joint && s.source.method == EstimateMethod::JointComponents &&
            (formal_joint_point || e->PostFitPeeling() || e->PostJointLocalResult() || e->GroupMemberResult()))
            throw std::invalid_argument("Unqualified Joint snapshot conflicts with current Joint stage state.");
        if(s.uncertainty.status==EvidenceStatus::Available && (!s.point || !s.uncertainty.covariance || !s.uncertainty.covariance->allFinite())) throw std::invalid_argument("Available uncertainty lacks finite point/covariance.");
        if(!s.source.atom_id.empty() && s.source.atom_id!=std::to_string(id)) throw std::invalid_argument("Stage atom identity mismatch.");
        if(e->PostFitPeeling())
        {
            if(!joint_snapshot) throw std::invalid_argument("Original-grid peeling requires its Joint snapshot.");
            if(joint_snapshot->selection_domain &&
                joint_snapshot->atom_ids.size()>joint_snapshot->selection_domain->target_indices.size() &&
                !joint_snapshot->halo_median_reconstruction &&
                std::any_of(e->PostFitPeeling()->samples.begin(),e->PostFitPeeling()->samples.end(),
                    [](const auto & sample){return sample.response.has_value();}))
                throw std::invalid_argument("Available grid peeling requires Halo reconstruction provenance.");
            const auto position=model.FindAtomPtr(id)->GetPosition();
            const auto reason=data_internal::ValidatePostFitPeelingResult(*e->PostFitPeeling(),s,joint_snapshot,&position);
            if(!reason.empty()) throw std::invalid_argument(reason);
        }
        if(e->PostJointQScore())
        {
            const auto reason=qscore_validation->Validate(*e->PostJointQScore(),*model.FindAtomPtr(id),s);
            if(!reason.empty()) throw std::invalid_argument(reason);
        }
        if(e->PostJointLocalResult())
        {
            const auto & local=*e->PostJointLocalResult();
            const auto & peeling=e->PostFitPeeling();
            const auto assessment=data_internal::AssessPostJointLocalResult(
                local,id,s,peeling ? &*peeling:nullptr);
            if(!assessment.valid) throw std::invalid_argument(assessment.reason);
        }
        if(e->GroupMemberResult() && s.source.method==EstimateMethod::JointComponents && e->GroupMemberResult()->evidence_source_id!=s.source.run_id) throw std::invalid_argument("Posterior source mismatch.");
    }
    const auto & groups=data.AtomGroupEntry();
    for(const auto key:groups.CollectGroupKeys())
    {
        const auto & summary=groups.GetParameterSummary(key); if(!summary) continue;
        if(summary->eligible_count!=summary->member_ids.size()) throw std::invalid_argument("Group eligibility count mismatch.");
        std::set<int> group_ids, selected_ids;
        std::vector<int> eligible_ids;
        std::vector<double> coefficients;
        std::size_t target_count{};
        for(const auto * atom:model.GetSelectedAtoms())
            if(data_internal::GetGroupKey(atom)==key) selected_ids.insert(atom->GetSerialID());
        for(const auto * atom:groups.GetMembers(key))
        {
            if(!group_ids.insert(atom->GetSerialID()).second || data_internal::GetGroupKey(atom)!=key)
                throw std::invalid_argument("Joint group member identity mismatch.");
            const auto * entry=data.FindAtomLocalEntry(*atom);
            if(!entry) continue;
            const auto & stage=entry->JointStageEstimate();
            if(stage.source.method!=EstimateMethod::JointComponents || stage.source.role!=FittingRole::Target) continue;
            ++target_count;
            if(stage.source.run_id.empty() || stage.source.run_id!=summary->source_id)
                throw std::invalid_argument("Joint group source mismatch.");
            const auto & local=entry->PostJointLocalResult();
            if(!local) continue;
            const auto & peeling=entry->PostFitPeeling();
            const auto assessment=data_internal::AssessPostJointLocalResult(
                *local,atom->GetSerialID(),stage,peeling ? &*peeling:nullptr);
            if(!assessment.valid) throw std::invalid_argument(assessment.reason);
            if(!assessment.group_eligible) continue;
            eligible_ids.push_back(atom->GetSerialID());
            coefficients.push_back(local->fixed_charge.coefficient);
        }
        if(group_ids!=selected_ids || eligible_ids!=summary->member_ids ||
            summary->point_count!=eligible_ids.size() ||
            summary->excluded_count!=target_count-eligible_ids.size())
            throw std::invalid_argument("Joint group eligibility/member identity mismatch.");
        if((summary->status==EvidenceStatus::Available)!=summary->inference.has_value() ||
            summary->status==EvidenceStatus::NotRun ||
            (summary->inference ? !summary->reason.empty() : summary->reason.empty()) ||
            (summary->inference && eligible_ids.size()<2))
            throw std::invalid_argument("Joint group inference status mismatch.");
        if(coefficients.empty())
        {
            if(summary->descriptive_charge_median || summary->descriptive_mean)
                throw std::invalid_argument("Empty Joint group has descriptive parameters.");
        }
        else
        {
            // Precise JSON double round-trips preserve the endpoint values and median exactly.
            const double median=array_helper::ComputeMedian(coefficients);
            if(!std::isfinite(median) || !summary->descriptive_charge_median ||
                *summary->descriptive_charge_median!=median || !summary->descriptive_mean)
                throw std::invalid_argument("Joint group charge median does not match eligible evidence.");
            std::vector<const GaussianModel3D *> models{&*summary->descriptive_mean};
            if(summary->inference)
            {
                const auto & inference=*summary->inference;
                models.insert(models.end(),{&inference.mean,&inference.mdpde,&inference.prior.GetModel()});
                if(!std::isfinite(inference.alpha_g) || inference.alpha_g<0.0 ||
                    !std::isnan(inference.prior.GetStandardDeviationModel().GetOffset()))
                    throw std::invalid_argument("Joint group has invalid alpha or inferred charge uncertainty.");
            }
            for(const auto * point:models)
            {
                GaussianModel3D::RequireFinitePositiveWidthModel(*point);
                if(point->GetAmplitude()<=0.0 || point->GetOffset()!=median)
                    throw std::invalid_argument("Joint group model does not match its charge median.");
            }
        }
        if(!summary->inference)
            for(const auto * atom:groups.GetMembers(key))
            {
                const auto * entry=data.FindAtomLocalEntry(*atom);
                if(entry && entry->GroupMemberResult())
                    throw std::invalid_argument("Unavailable Joint group has member posteriors.");
            }
        if (unqualified_joint && summary->inference)
        {
            const bool has_joint_member=std::any_of(summary->member_ids.begin(),summary->member_ids.end(),
                [&](const int id) {
                    const auto * atom=model.FindAtomPtr(id);
                    const auto * entry=atom ? data.FindAtomLocalEntry(*atom) : nullptr;
                    return entry && entry->JointStageEstimate().source.method == EstimateMethod::JointComponents;
                });
            if (has_joint_member)
                throw std::invalid_argument("Unqualified Joint snapshot conflicts with current group inference.");
        }
        if(!summary->inference) continue;
        const auto members=model.GetAnalysisView().GetGroupParameterSummary(key)->inference->member_results;
        if(members.size()!=summary->member_ids.size()) throw std::invalid_argument("Group posterior count mismatch.");
        for(std::size_t i=0;i<members.size();++i)
        {
            const auto * atom=model.FindAtomPtr(summary->member_ids[i]);
            const auto * entry=atom ? data.FindAtomLocalEntry(*atom):nullptr;
            if(!entry || !entry->GroupMemberResult() || members[i].evidence_source_id!=summary->source_id ||
                entry->GroupMemberResult()->posterior.GetModel().ToVector()!=members[i].posterior.GetModel().ToVector())
                throw std::invalid_argument("Group/member posterior source mismatch.");
            const auto & stage=entry->JointStageEstimate();
            if(stage.source.method==EstimateMethod::JointComponents)
            {
                const auto & local=entry->PostJointLocalResult();
                const auto & member=*entry->GroupMemberResult();
                if(!local || stage.source.run_id!=summary->source_id)
                    throw std::invalid_argument("Group posterior lacks current post-Joint local measurement information.");
                const auto & peeling=entry->PostFitPeeling();
                const auto assessment=data_internal::AssessPostJointLocalResult(
                    *local,atom->GetSerialID(),stage,peeling ? &*peeling:nullptr);
                if(!assessment.valid || !assessment.group_eligible ||
                    member.charge_inferred || !member.fixed_charge ||
                    member.fixed_charge->coefficient!=local->fixed_charge.coefficient ||
                    member.fixed_charge->width!=local->fixed_charge.width ||
                    !std::isfinite(member.fixed_charge->coefficient) ||
                    !std::isfinite(member.fixed_charge->width) || member.fixed_charge->width<=0.0 ||
                    member.posterior.GetModel().GetOffset()!=0.0 ||
                    !std::isnan(member.posterior.GetStandardDeviationModel().GetOffset()))
                    throw std::invalid_argument("Group posterior lacks current post-Joint local measurement information.");
            }
        }
    }
    if(!data.joint_result) return;
    if(!IsQualifiedJointAnalysisResult(*data.joint_result)) return;
    const auto mapped=data_internal::BuildJointStageEstimates(*data.joint_result,"");
    for(const auto & [id,expected]:mapped)
    {
        const auto * e=data.FindAtomLocalEntry(*model.FindAtomPtr(id));
        if(!e || e->JointStageEstimate().source.method!=EstimateMethod::JointComponents) throw std::invalid_argument("Joint snapshot lacks a matching JointStage record.");
        const auto & actual=e->JointStageEstimate();
        if(actual.convergence!=expected.convergence || actual.source.component_id!=expected.source.component_id || actual.source.role!=expected.source.role ||
            actual.source.fixed_from_initial_stage!=expected.source.fixed_from_initial_stage || actual.point.has_value()!=expected.point.has_value() ||
            (actual.point && actual.point->ToVector()!=expected.point->ToVector())) throw std::invalid_argument("Joint snapshot and JointStage record mismatch.");
    }
}
}
