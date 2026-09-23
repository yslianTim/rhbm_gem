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
V Pack(const GaussianModel3D & p) { return A{p.GetAmplitude(), p.GetWidth(), p.GetOffset()}; }
GaussianModel3D Point(const V & v, bool positive_width=true) { const auto & a = v.as_array(); if (a.size()!=3) throw std::invalid_argument("Invalid stage point."); GaussianModel3D p(Number(a[0]),Number(a[1]),Number(a[2])); GaussianModel3D::RequireFiniteModel(p); if(positive_width) GaussianModel3D::RequireFinitePositiveWidthModel(p); return p; }
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
    return O{{"point",Pack(p.GetModel())},{"sd",A{s.GetAmplitude(),s.GetWidth(),std::isfinite(s.GetOffset()) ? V(s.GetOffset()):V{}}}};
}
GaussianModel3DWithUncertainty Uncertain(const V & v)
{
    const auto & o=v.as_object(); const auto & a=o.at("sd").as_array();
    if(a.size()!=3) throw std::invalid_argument("Invalid uncertainty shape.");
    return {Point(o.at("point"),false),GaussianModel3DUncertainty(Number(a[0]),Number(a[1]),a[2].is_null() ? std::numeric_limits<double>::quiet_NaN():Number(a[2]))};
}
V Pack(const GroupGaussianMemberResult & m) { return O{{"posterior",Pack(m.posterior)},{"outlier",m.is_outlier},{"distance",m.statistical_distance},{"source",m.evidence_source_id},{"charge_inferred",m.charge_inferred},{"covariance",m.parameter_covariance ? Matrix(*m.parameter_covariance):V{}}}; }
GroupGaussianMemberResult Member(const V & v) { const auto & o=v.as_object(); GroupGaussianMemberResult m; m.posterior=Uncertain(o.at("posterior")); m.is_outlier=Read<bool>(o,"outlier"); m.statistical_distance=Number(o.at("distance")); m.evidence_source_id=Read<std::string>(o,"source"); m.charge_inferred=Read<bool>(o,"charge_inferred"); if(!o.at("covariance").is_null()) m.parameter_covariance=Matrix<2,2>(o.at("covariance")); return m; }
V Pack(const GroupParameterSummary & s)
{
    V inference;
    if(s.inference) { const auto & i=*s.inference; A members; for(const auto & m:i.member_results) members.push_back(Pack(m)); inference=O{{"alpha",i.alpha_g},{"mean",Pack(i.mean)},{"mdpde",Pack(i.mdpde)},{"prior",Pack(i.prior)},{"members",members}}; }
    return O{{"status",static_cast<int>(s.status)},{"reason",s.reason},{"source",s.source_id},{"correlation",s.correlation_policy},{"points",s.point_count},{"eligible",s.eligible_count},{"excluded",s.excluded_count},{"mean",s.descriptive_mean ? Pack(*s.descriptive_mean):V{}},{"inference",inference},{"ids",j::value_from(s.member_ids)}};
}
GroupParameterSummary Summary(const V & v)
{
    const auto & o=v.as_object(); GroupParameterSummary s; s.status=Enum<EvidenceStatus>(o.at("status"),3); s.reason=Read<std::string>(o,"reason"); s.source_id=Read<std::string>(o,"source"); s.correlation_policy=Read<std::string>(o,"correlation"); s.point_count=Read<std::size_t>(o,"points"); s.eligible_count=Read<std::size_t>(o,"eligible"); s.excluded_count=Read<std::size_t>(o,"excluded"); s.member_ids=Read<std::vector<int>>(o,"ids");
    if(!o.at("mean").is_null()) s.descriptive_mean=Point(o.at("mean"));
    if(!o.at("inference").is_null()) { const auto & i=o.at("inference").as_object(); GroupGaussianResult r; r.alpha_g=Number(i.at("alpha")); r.mean=Point(i.at("mean")); r.mdpde=Point(i.at("mdpde")); r.prior=Uncertain(i.at("prior")); for(const auto & m:i.at("members").as_array()) r.member_results.push_back(Member(m)); if(r.member_results.size()!=s.member_ids.size()) throw std::invalid_argument("Posterior member count mismatch."); s.inference=std::move(r); }
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
            {"geometry",e->SampleGeometryAvailable()},{"raw",Samples(e->RawSamplingEntries(),e->SampleGeometryAvailable())},{"peeled",Samples(e->PeelingSamplingEntries(),e->SampleGeometryAvailable())},
            {"peeling",peeling},{"evidence",e->GroupEvidence() ? Pack(*e->GroupEvidence()):V{}},{"posterior",e->GroupMemberResult() ? Pack(*e->GroupMemberResult()):V{}}});
    }
    const auto & g=data.AtomGroupEntry();
    for(const auto key:g.CollectGroupKeys()) { A ids; for(const auto * a:g.GetMembers(key)) ids.push_back(a->GetSerialID()); groups.push_back(O{{"key",static_cast<std::uint64_t>(key)},{"members",ids},{"summary",g.GetParameterSummary(key) ? Pack(*model.GetAnalysisView().GetGroupParameterSummary(key)):V{}}}); }
    return j::serialize(O{{"version",1},{"atoms",atoms},{"groups",groups}});
}
void Decode(ModelObject & model, const std::string & text)
{
    j::parse_options parse_options; parse_options.numbers=j::number_precision::precise;
    const auto document=j::parse(text,{},parse_options); const auto & root=document.as_object();
    if(Read<int>(root,"version")!=1) throw std::invalid_argument("Unsupported neutral stage format.");
    auto & data=ModelAnalysisData::Of(model); std::set<int> seen;
    for(const auto & v:root.at("atoms").as_array())
    {
        const auto & o=v.as_object(); const int id=Read<int>(o,"id"); auto * atom=model.FindAtomPtr(id);
        if(!atom || !seen.insert(id).second) throw std::invalid_argument("Invalid neutral atom identity.");
        auto & e=data.EnsureAtomLocalEntry(*atom);
        e.SetStageEstimate(FittingStage::First,Stage(o.at("first"))); e.SetStageEstimate(FittingStage::Second,Stage(o.at("second")));
        const bool geometry=Read<bool>(o,"geometry"); e.SetSampleGeometryAvailable(geometry);
        e.SetRawSamplingEntries(Samples(o.at("raw"),geometry));
        if(e.StageEstimate(FittingStage::Second).source.method != EstimateMethod::JointComponents)
            e.SetPeelingSamplingEntries(Samples(o.at("peeled"),geometry));
        if(!o.at("peeling").is_null()) { const auto & p=o.at("peeling").as_object(); PostFitPeelingResult r; r.source=Source(p.at("source")); r.mode=Read<std::string>(p,"mode"); r.neighbor_count=Read<std::size_t>(p,"neighbors"); for(const auto & sample:p.at("samples").as_array()) { const auto & s=sample.as_object(); PeelingSampleEstimate x; x.reason=Read<std::string>(s,"reason"); if(!s.at("response").is_null()) x.response=Number(s.at("response")); r.samples.push_back(x); } e.SetPostFitPeeling(std::move(r)); }
        if(!o.at("evidence").is_null()) e.SetGroupEvidence(Evidence(o.at("evidence")));
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
        if(!o.at("summary").is_null()) { auto s=Summary(o.at("summary")); for(int id:s.member_ids) if(!members.contains(id)) throw std::invalid_argument("Posterior outside fitted group."); if (s.inference)
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
    if(data.joint_result) MapSnapshot(model);
    for(auto & [id,e]:data.AtomLocalEntries()) { (void)id; e->SetSampleGeometryAvailable(false); }
}
void MapSnapshot(ModelObject & model)
{
    auto & data=ModelAnalysisData::Of(model);
    const auto snapshot = *data.joint_result;
    data_internal::ApplyJointStageEstimates(model,snapshot,"legacy-joint-snapshot");
    data.joint_result = snapshot;
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
        if(actual.source.component_id!=expected.source.component_id || actual.source.role!=expected.source.role || actual.point.has_value()!=expected.point.has_value() ||
            (actual.point && actual.point->ToVector()!=expected.point->ToVector())) throw std::invalid_argument("Joint snapshot/Second mismatch.");
    }
}
}
