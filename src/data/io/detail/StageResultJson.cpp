#include "StageResultJson.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <boost/json.hpp>
#include <cmath>
#include <algorithm>
#include <limits>
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
V Pack(const EstimateSource & s) { return O{{"method",static_cast<int>(s.method)},{"atom",s.atom_id},{"component",s.component_id},{"run",s.run_id},{"role",static_cast<int>(s.role)}}; }
EstimateSource Source(const V & v) { const auto & o=v.as_object(); return {Enum<EstimateMethod>(o.at("method"),3),Read<std::string>(o,"atom"),Read<std::string>(o,"component"),Read<std::string>(o,"run"),Enum<FittingRole>(o.at("role"),2)}; }
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
V Pack(const GroupParameterEvidence & e) { return O{{"atom",e.atom_id},{"component",e.component_id},{"source",e.source_id},{"status",static_cast<int>(e.status)},{"reason",e.reason},{"method",e.uncertainty_method},{"correlation",e.correlation_policy},{"estimate",e.estimate ? Matrix(*e.estimate):V{}},{"covariance",e.covariance ? Matrix(*e.covariance):V{}}}; }
GroupParameterEvidence Evidence(const V & v)
{
    const auto & o=v.as_object(); GroupParameterEvidence e;
    e.atom_id=Read<std::string>(o,"atom"); e.component_id=Read<std::string>(o,"component"); e.source_id=Read<std::string>(o,"source");
    e.status=Enum<EvidenceStatus>(o.at("status"),3); e.reason=Read<std::string>(o,"reason"); e.uncertainty_method=Read<std::string>(o,"method"); e.correlation_policy=Read<std::string>(o,"correlation");
    if(!o.at("estimate").is_null()) e.estimate=Matrix<2,1>(o.at("estimate"));
    if(!o.at("covariance").is_null()) e.covariance=Matrix<2,2>(o.at("covariance"));
    return e;
}
V Pack(const GaussianModel3DWithUncertainty & p)
{
    const auto & s=p.GetStandardDeviationModel();
    return O{{"point",Pack(p.GetModel())},{"sd",A{NativeNumber(s.GetAmplitude()),NativeNumber(s.GetWidth()),NativeNumber(s.GetOffset())}}};
}
GaussianModel3DWithUncertainty Uncertain(const V & v)
{
    const auto & o=v.as_object(); const auto & a=o.at("sd").as_array();
    if(a.size()!=3) throw std::invalid_argument("Invalid uncertainty shape.");
    return {Point(o.at("point"),false),GaussianModel3DUncertainty(NativeNumber(a[0]),NativeNumber(a[1]),a[2].is_null() ? std::numeric_limits<double>::quiet_NaN():NativeNumber(a[2]))};
}
V OptionalNumber(const std::optional<double> & v) { return v ? NativeNumber(*v) : V{}; }
std::optional<double> OptionalNumber(const V & v) { return v.is_null() ? std::nullopt : std::optional<double>{NativeNumber(v)}; }
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
    return O{{"status",static_cast<int>(d.status)},{"variance",NativeNumber(d.sigma_square)},{"iterations",d.iterations.iterations},
        {"beta_change",OptionalNumber(d.iterations.squared_beta_change)},{"variance_change",OptionalNumber(d.iterations.relative_variance_change)},
        {"refinement",refinement}};
}
LocalFitDiagnostics Diagnostics(const V & v)
{
    const auto & o=v.as_object(); LocalFitDiagnostics d;
    d.status=Enum<RHBMEstimationStatus>(o.at("status"),4); d.sigma_square=NativeNumber(o.at("variance"));
    d.iterations={Read<int>(o,"iterations"),OptionalNumber(o.at("beta_change")),OptionalNumber(o.at("variance_change"))};
    if (!o.at("refinement").is_null())
    {
        const auto & r=o.at("refinement").as_object(); RHBMEndpointRefinementDiagnostics x;
        x.accepted=Read<bool>(r,"accepted"); x.reason=Read<std::string>(r,"reason"); x.candidate_equation_evaluations=Read<int>(r,"evaluations");
        x.reference_updates=Read<int>(r,"updates"); x.reference_stop=Read<std::string>(r,"stop");
        x.original_residual=OptionalNumber(r.at("original")); x.candidate_residual=OptionalNumber(r.at("candidate")); x.reference_residual=OptionalNumber(r.at("reference"));
        if (!r.at("coordinate_difference").is_null())
        {
            const auto & a=r.at("coordinate_difference").as_array();
            if(a.size()!=3) throw std::invalid_argument("Invalid refinement coordinate shape.");
            x.relative_coordinate_difference=std::array<double,3>{NativeNumber(a[0]),NativeNumber(a[1]),NativeNumber(a[2])};
        }
        x.weight_max_difference=OptionalNumber(r.at("weight_difference"));
        if (!r.at("floor_masks_equal").is_null()) x.floor_masks_equal=Read<bool>(r,"floor_masks_equal");
        d.refinement=std::move(x);
    }
    return d;
}
V Native(const LocalPotentialEntry & e, FittingStage stage)
{
    if(e.StageEstimate(stage).source.method==EstimateMethod::JointComponents) return {};
    const auto r=e.GaussianResult(stage,false); const auto & sd=r.mdpde.GetStandardDeviationModel();
    return O{{"alpha",r.alpha_r},{"ols",Pack(r.ols.GetModel())},{"ols_sd",r.uncertainty_recorded ? Pack(r.ols).as_object().at("sd"):V{}},
        {"mdpde_sd",r.uncertainty_recorded ? V(A{NativeNumber(sd.GetAmplitude()),NativeNumber(sd.GetWidth()),NativeNumber(sd.GetOffset())}):V{}},
        {"unfitted",e.StageEstimate(stage).point ? V{} : Pack(r.mdpde.GetModel())},
        {"diagnostics",r.diagnostics ? Pack(*r.diagnostics):V{}}};
}
void RestoreNative(LocalPotentialEntry & e, FittingStage stage, const V & native, LocalStageEstimate estimate)
{
    if(!native.is_null())
    {
        if(estimate.source.method==EstimateMethod::JointComponents) throw std::invalid_argument("Joint stage contains native diagnostics.");
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
        e.SetGaussianResult(stage,std::move(r));
    }
    e.SetStageEstimate(stage,std::move(estimate));
}
V Pack(const GroupGaussianMemberResult & m) { return O{{"posterior",Pack(m.posterior)},{"outlier",m.is_outlier},{"distance",m.statistical_distance},{"source",m.evidence_source_id},{"charge_inferred",m.charge_inferred},{"covariance",m.parameter_covariance ? Matrix(*m.parameter_covariance):V{}}}; }
GroupGaussianMemberResult Member(const V & v) { const auto & o=v.as_object(); GroupGaussianMemberResult m; m.posterior=Uncertain(o.at("posterior")); m.is_outlier=Read<bool>(o,"outlier"); m.statistical_distance=Number(o.at("distance")); m.evidence_source_id=Read<std::string>(o,"source"); m.charge_inferred=Read<bool>(o,"charge_inferred"); if(!o.at("covariance").is_null()) m.parameter_covariance=Matrix<2,2>(o.at("covariance")); return m; }
V Pack(const GroupParameterSummary & s)
{
    V inference;
    if(s.inference) { const auto & i=*s.inference; inference=O{{"alpha",i.alpha_g},{"mean",Pack(i.mean)},{"mdpde",Pack(i.mdpde)},{"prior",Pack(i.prior)}}; }
    return O{{"status",static_cast<int>(s.status)},{"reason",s.reason},{"source",s.source_id},{"correlation",s.correlation_policy},{"points",s.point_count},{"eligible",s.eligible_count},{"excluded",s.excluded_count},{"mean",s.descriptive_mean ? Pack(*s.descriptive_mean):V{}},{"inference",inference},{"ids",j::value_from(s.member_ids)}};
}
GroupParameterSummary Summary(const V & v, bool legacy)
{
    const auto & o=v.as_object(); GroupParameterSummary s; s.status=Enum<EvidenceStatus>(o.at("status"),3); s.reason=Read<std::string>(o,"reason"); s.source_id=Read<std::string>(o,"source"); s.correlation_policy=Read<std::string>(o,"correlation"); s.point_count=Read<std::size_t>(o,"points"); s.eligible_count=Read<std::size_t>(o,"eligible"); s.excluded_count=Read<std::size_t>(o,"excluded"); s.member_ids=Read<std::vector<int>>(o,"ids");
    if(!o.at("mean").is_null()) s.descriptive_mean=Point(o.at("mean"));
    if(!o.at("inference").is_null()) { const auto & i=o.at("inference").as_object(); GroupGaussianResult r; r.alpha_g=Number(i.at("alpha")); r.mean=Point(i.at("mean")); r.mdpde=Point(i.at("mdpde")); r.prior=Uncertain(i.at("prior")); if(legacy) { for(const auto & m:i.at("members").as_array()) r.member_results.push_back(Member(m)); if(r.member_results.size()!=s.member_ids.size()) throw std::invalid_argument("Posterior member count mismatch."); } s.inference=std::move(r); }
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
        if(e->PostFitPeeling()) { const auto & p=*e->PostFitPeeling(); A samples; for(const auto & s:p.samples) samples.push_back(O{{"response",s.response ? V(*s.response):V{}},{"reason",s.reason}}); peeling=O{{"source",Pack(p.source)},{"mode",p.mode},{"neighbors",p.neighbor_count},{"samples",samples}}; }
        atoms.push_back(O{{"id",atom->GetSerialID()},{"first",Pack(e->StageEstimate(FittingStage::First))},{"second",Pack(e->StageEstimate(FittingStage::Second))},
            {"first_native",Native(*e,FittingStage::First)},{"second_native",Native(*e,FittingStage::Second)},{"neighbors",(e->PostFitPeeling() || e->StageEstimate(FittingStage::Second).source.method==EstimateMethod::JointComponents) ? V{}:V(e->NeighborCountForPeeling())},
            {"geometry",e->SampleGeometryAvailable()},{"raw",Samples(e->RawSamplingEntries(),e->SampleGeometryAvailable())},{"peeled",(e->PostFitPeeling() || e->StageEstimate(FittingStage::Second).source.method==EstimateMethod::JointComponents) ? V{}:Samples(e->PeelingSamplingEntries(),e->SampleGeometryAvailable())},
            {"peeling",peeling},{"evidence",e->GroupEvidence() ? Pack(*e->GroupEvidence()):V{}},{"posterior",e->GroupMemberResult() ? Pack(*e->GroupMemberResult()):V{}}});
    }
    const auto & g=data.AtomGroupEntry();
    for(const auto key:g.CollectGroupKeys()) { A ids; for(const auto * a:g.GetMembers(key)) ids.push_back(a->GetSerialID()); groups.push_back(O{{"key",static_cast<std::uint64_t>(key)},{"members",ids},{"alpha",g.GetParameterSummary(key) && !g.GetParameterSummary(key)->inference ? V(g.GetAlphaG(key)):V{}},{"native",g.GetParameterSummary(key) ? V{}:V(O{{"alpha",g.GetAlphaG(key)},{"mean",Pack(g.GetMean(key))},{"mdpde",Pack(g.GetMDPDE(key))},{"prior",Pack(g.GetPriorWithUncertainty(key))}})},{"summary",g.GetParameterSummary(key) ? Pack(*g.GetParameterSummary(key)):V{}}}); }
    return j::serialize(O{{"version",2},{"atoms",atoms},{"groups",groups}});
}
void Decode(ModelObject & model, const std::string & text, int expected_version, const std::vector<GroupKey> & legacy_groups)
{
    j::parse_options parse_options; parse_options.numbers=j::number_precision::precise;
    const auto document=j::parse(text,{},parse_options); const auto & root=document.as_object();
    const int version=Read<int>(root,"version"); const bool legacy=version==1;
    if(version!=expected_version || (version!=1 && version!=2)) throw std::invalid_argument("Unsupported neutral stage format.");
    auto & data=ModelAnalysisData::Of(model); std::set<int> seen;
    for(const auto & v:root.at("atoms").as_array())
    {
        const auto & o=v.as_object(); const int id=Read<int>(o,"id"); auto * atom=model.FindAtomPtr(id);
        if(!atom || !seen.insert(id).second) throw std::invalid_argument("Invalid neutral atom identity.");
        auto & e=data.EnsureAtomLocalEntry(*atom);
        if(legacy)
        {
            for(const auto stage:{FittingStage::First,FittingStage::Second})
            {
                const auto incoming=Stage(o.at(stage==FittingStage::First ? "first":"second"));
                const auto & previous=e.StageEstimate(stage).point;
                if((incoming.source.method==EstimateMethod::JointComponents && previous) ||
                    (incoming.source.method!=EstimateMethod::JointComponents &&
                    (incoming.point.has_value()!=previous.has_value() ||
                    (incoming.point && incoming.point->ToVector()!=previous->ToVector()))))
                    throw std::invalid_argument("Legacy stage point mismatch.");
            }
            e.SetStageEstimate(FittingStage::First,Stage(o.at("first"))); e.SetStageEstimate(FittingStage::Second,Stage(o.at("second")));
        }
        else
        {
            RestoreNative(e,FittingStage::First,o.at("first_native"),Stage(o.at("first")));
            RestoreNative(e,FittingStage::Second,o.at("second_native"),Stage(o.at("second")));
        }
        if(!legacy && (!o.at("peeling").is_null() || e.StageEstimate(FittingStage::Second).source.method==EstimateMethod::JointComponents) &&
            (!o.at("peeled").is_null() || !o.at("neighbors").is_null()))
            throw std::invalid_argument("Canonical paired peeling contains a compact mirror.");
        const bool geometry=Read<bool>(o,"geometry"); e.SetSampleGeometryAvailable(geometry);
        const auto raw=Samples(o.at("raw"),geometry);
        if(legacy && Samples(raw,false)!=Samples(e.RawSamplingEntries(),false)) throw std::invalid_argument("Legacy raw samples mismatch.");
        if(legacy && Samples(Samples(o.at("peeled"),geometry),false)!=Samples(e.PeelingSamplingEntries(),false)) throw std::invalid_argument("Legacy peeling samples mismatch.");
        e.SetRawSamplingEntries(raw);
        if(!o.at("peeled").is_null() && e.StageEstimate(FittingStage::Second).source.method != EstimateMethod::JointComponents)
            e.SetPeelingSamplingEntries(Samples(o.at("peeled"),geometry));
        if(!legacy && !o.at("neighbors").is_null()) e.SetNeighborCountForPeeling(Read<int>(o,"neighbors"));
        if(!o.at("peeling").is_null()) { const auto & p=o.at("peeling").as_object(); PostFitPeelingResult r; if(legacy && Read<int>(p,"neighbors")!=e.NeighborCountForPeeling()) throw std::invalid_argument("Legacy peeling neighbor count mismatch."); r.source=Source(p.at("source")); r.mode=Read<std::string>(p,"mode"); r.neighbor_count=Read<std::size_t>(p,"neighbors"); for(const auto & sample:p.at("samples").as_array()) { const auto & s=sample.as_object(); PeelingSampleEstimate x; x.reason=Read<std::string>(s,"reason"); if(!s.at("response").is_null()) x.response=Number(s.at("response")); r.samples.push_back(x); } e.SetPostFitPeeling(std::move(r)); }
        if(legacy && e.PostFitPeeling() && Samples(e.PeelingSamplingEntries(),geometry)!=o.at("peeled"))
            throw std::invalid_argument("Legacy paired/compact peeling mismatch.");
        if(!o.at("evidence").is_null()) e.SetGroupEvidence(Evidence(o.at("evidence")));
        if(legacy && e.GroupMemberResult())
        {
            if(o.at("posterior").is_null()) throw std::invalid_argument("Legacy posterior missing.");
            const auto incoming=Member(o.at("posterior")); const auto & old=*e.GroupMemberResult();
            if(Pack(incoming.posterior)!=Pack(old.posterior) || incoming.is_outlier!=old.is_outlier || incoming.statistical_distance!=old.statistical_distance)
                throw std::invalid_argument("Legacy posterior mismatch.");
        }
        e.ClearGroupMemberResult(); if(!o.at("posterior").is_null()) e.SetGroupMemberResult(Member(o.at("posterior")));
    }
    if(seen.size()!=data.AtomLocalEntries().size()) throw std::invalid_argument("Incomplete neutral local records.");
    auto & groups=data.AtomGroupEntry(); groups.ClearMembers(); std::set<GroupKey> keys;
    for(const auto & v:root.at("groups").as_array())
    {
        const auto & o=v.as_object(); const auto key=Read<GroupKey>(o,"key");
        if(!keys.insert(key).second) throw std::invalid_argument("Duplicate neutral group.");
        std::set<int> members;
        for(const int id:Read<std::vector<int>>(o,"members")) {auto * atom=model.FindAtomPtr(id); if(!atom || !members.insert(id).second) throw std::invalid_argument("Invalid group member identity."); groups.AddMember(key,*atom);}
        if(!legacy && !o.at("alpha").is_null()) groups.SetAlphaG(key,Number(o.at("alpha")));
        if(!legacy && !o.at("native").is_null())
        {
            const auto & n=o.at("native").as_object(); GroupGaussianResult result;
            result.alpha_g=Number(n.at("alpha")); result.mean=Point(n.at("mean"),false);
            result.mdpde=Point(n.at("mdpde"),false); result.prior=Uncertain(n.at("prior"));
            groups.SetGaussianResult(key,result);
        }
        if(!o.at("summary").is_null()) { auto s=Summary(o.at("summary"),legacy);
            if(legacy && std::find(legacy_groups.begin(),legacy_groups.end(),key)!=legacy_groups.end())
            {
                const GaussianModel3D missing{0,0};
                if(Pack(groups.GetMean(key))!=Pack(s.descriptive_mean.value_or(missing)) ||
                    Pack(groups.GetMDPDE(key))!=Pack(s.inference ? s.inference->mdpde:missing) ||
                    Pack(groups.GetPriorWithUncertainty(key))!=Pack(s.inference ? s.inference->prior:GaussianModel3DWithUncertainty{missing,{}}) ||
                    (s.inference && groups.GetAlphaG(key)!=s.inference->alpha_g))
                    throw std::invalid_argument("Legacy group statistics mismatch.");
            }
            for(int id:s.member_ids) if(!members.contains(id)) throw std::invalid_argument("Posterior outside fitted group."); if (legacy && s.inference)
            {
                for (std::size_t index=0; index<s.member_ids.size(); ++index)
                {
                    const auto * entry=data.FindAtomLocalEntry(*model.FindAtomPtr(s.member_ids[index]));
                    if (!entry || !entry->GroupMemberResult() || Pack(*entry->GroupMemberResult()) != Pack(s.inference->member_results[index]))
                        throw std::invalid_argument("Legacy group/member posterior mismatch.");
                }
            }
            groups.SetParameterSummary(key,std::move(s)); }
    }
    for(const auto key:legacy_groups) if(!keys.contains(key)) throw std::invalid_argument("Legacy group missing from neutral document.");
    Validate(model);
}
void AdaptLegacy(ModelObject & model)
{
    auto & data=ModelAnalysisData::Of(model);
    for(auto & [id,e]:data.AtomLocalEntries())
    {
        e->SetSampleGeometryAvailable(false);
        for(const auto stage:{FittingStage::First,FittingStage::Second}) {auto s=e->StageEstimate(stage); s.source.atom_id=std::to_string(id); e->SetStageEstimate(stage,std::move(s));}
    }
    if(data.joint_result)
    {
        for(const auto & id:data.joint_result->atom_ids)
        {
            const auto * entry=data.FindAtomLocalEntry(*model.FindAtomPtr(std::stoi(id)));
            if(entry && entry->StageEstimate(FittingStage::Second).point)
                throw std::invalid_argument("Legacy snapshot conflicts with native Second point.");
        }
        MapSnapshot(model);
    }
    for(auto & [id,e]:data.AtomLocalEntries()) { (void)id; e->SetSampleGeometryAvailable(false); }
}
void MapSnapshot(ModelObject & model)
{
    auto & data=ModelAnalysisData::Of(model);
    auto snapshot = std::move(*data.joint_result);
    data_internal::ApplyJointStageEstimates(model,snapshot,"legacy-joint-snapshot");
    data.joint_result = std::move(snapshot);
    auto & groups=data.AtomGroupEntry(); groups.ClearMembers();
    for(const auto key:groups.CollectGroupKeys()) {GroupParameterSummary s; s.reason="not-recorded"; s.source_id="legacy-joint-snapshot"; groups.SetParameterSummary(key,s);}
    for(const auto & id:data.joint_result->atom_ids)
    {
        auto * atom=model.FindAtomPtr(std::stoi(id)); auto & e=data.EnsureAtomLocalEntry(*atom);
        e.ClearPeeling(); e.ClearGroupMemberResult();
        const auto & stage=e.StageEstimate(FittingStage::Second);
        GroupParameterEvidence evidence; evidence.atom_id=id; evidence.component_id=stage.source.component_id; evidence.source_id=stage.source.run_id; evidence.reason="not-recorded";
        e.SetGroupEvidence(evidence);
        if(stage.source.role==FittingRole::Target) {const auto key=data_internal::GetGroupKey(atom); groups.AddMember(key,*atom); GroupParameterSummary s; s.reason="not-recorded"; s.source_id=stage.source.run_id; groups.SetParameterSummary(key,s);}
    }
}

void Validate(const ModelObject & model)
{
    const auto & data=ModelAnalysisData::Of(model);
    for(const auto & [id,e]:data.AtomLocalEntries())
    {
        const auto & s=e->StageEstimate(FittingStage::Second);
        if(s.uncertainty.status==EvidenceStatus::Available && (!s.point || !s.uncertainty.covariance || !s.uncertainty.covariance->allFinite())) throw std::invalid_argument("Available uncertainty lacks finite point/covariance.");
        if(e->GroupEvidence() && e->GroupEvidence()->status==EvidenceStatus::Available && (!e->GroupEvidence()->estimate || !e->GroupEvidence()->covariance || !e->GroupEvidence()->estimate->allFinite() || !e->GroupEvidence()->covariance->allFinite())) throw std::invalid_argument("Available group evidence lacks finite parameters/covariance.");
        if(!s.source.atom_id.empty() && s.source.atom_id!=std::to_string(id)) throw std::invalid_argument("Stage atom identity mismatch.");
        if(e->PostFitPeeling()) {const auto & p=*e->PostFitPeeling(); if(p.samples.size()!=e->RawSamplingEntries().size() || p.source.run_id!=s.source.run_id || p.source.atom_id!=s.source.atom_id) throw std::invalid_argument("Peeling source/shape mismatch.");}
        if(e->GroupEvidence()) {const auto & g=*e->GroupEvidence(); if(g.atom_id!=s.source.atom_id || g.component_id!=s.source.component_id || g.source_id!=s.source.run_id) throw std::invalid_argument("Group evidence source mismatch.");}
        if(e->GroupMemberResult() && s.source.method==EstimateMethod::JointComponents && e->GroupMemberResult()->evidence_source_id!=s.source.run_id) throw std::invalid_argument("Posterior source mismatch.");
    }
    const auto & groups=data.AtomGroupEntry();
    for(const auto key:groups.CollectGroupKeys())
    {
        const auto & summary=groups.GetParameterSummary(key); if(!summary) continue;
        if(summary->eligible_count!=summary->member_ids.size()) throw std::invalid_argument("Group eligibility count mismatch.");
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
        }
    }
    if(!data.joint_result) return;
    const auto mapped=data_internal::BuildJointStageEstimates(*data.joint_result,"");
    for(const auto & [id,expected]:mapped)
    {
        const auto * e=data.FindAtomLocalEntry(*model.FindAtomPtr(id));
        if(!e || e->StageEstimate(FittingStage::Second).source.method!=EstimateMethod::JointComponents) throw std::invalid_argument("Joint snapshot lacks a matching Second record.");
        const auto & actual=e->StageEstimate(FittingStage::Second);
        if(actual.convergence!=expected.convergence || actual.source.component_id!=expected.source.component_id || actual.source.role!=expected.source.role || actual.point.has_value()!=expected.point.has_value() ||
            (actual.point && actual.point->ToVector()!=expected.point->ToVector())) throw std::invalid_argument("Joint snapshot/Second mismatch.");
    }
}
}
