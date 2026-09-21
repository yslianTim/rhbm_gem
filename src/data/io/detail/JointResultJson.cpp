#include "JointResultJson.hpp"
#include <boost/json.hpp>
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
j::value Number(std::optional<double> x)
{
    if(x && !std::isfinite(*x)) throw std::invalid_argument("Nonfinite joint result number.");
    return x ? j::value(*x) : j::value(nullptr);
}
std::optional<double> OptionalNumber(const j::value & v)
{return v.is_null() ? std::nullopt : std::optional<double>(j::value_to<double>(v));}
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
    return {Read<std::vector<double>>(o,"ac"),Read<std::vector<double>>(o,"b"),Read<std::vector<double>>(o,"log_b"),
        Read<std::vector<double>>(o,"width_gradient"),Read<double>(o,"objective")};
}
j::value OptionalState(const std::optional<JointState> & x) {return x ? j::value(State(*x)) : j::value(nullptr);}
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
        out.emplace_back(Object{{"name",x.name},{"scope",ScopeText(x.scope)},
            {"rank",x.rank},{"threshold",Number(x.threshold)},{"singular_values",j::value_from(x.singular_values)}});
    }
    return out;
}
std::vector<JointRankEvidence> ReadRanks(const j::value & v)
{
    std::vector<JointRankEvidence> out;
    for(const auto & item:v.as_array())
    {
        const auto & o=item.as_object();
        out.push_back({Read<std::string>(o,"name"),Scope(o),Read<std::size_t>(o,"rank"),Read<double>(o,"threshold"),
            Read<std::vector<double>>(o,"singular_values")});
    }
    return out;
}
Object Initialization(const JointInitialization & x)
{
    j::array atoms;
    for(const auto & a:x.atoms) atoms.emplace_back(Object{{"id",a.id},{"ols",Diagnostics(a.ols)},{"mdpde",Diagnostics(a.mdpde)},
        {"alpha",std::isfinite(a.alpha) ? j::value(a.alpha) : j::value(nullptr)},{"sample_count",a.sample_count},
        {"native_status",a.native_status ? j::value(*a.native_status) : j::value(nullptr)}});
    return {{"valid",x.valid},{"reason",x.reason},{"b",Diagnostics(x.b)},{"atoms",std::move(atoms)}};
}
JointInitialization ReadInitialization(const j::value & v)
{
    const auto & o=v.as_object(); JointInitialization x;
    x.valid=Read<bool>(o,"valid"); x.reason=Read<std::string>(o,"reason"); x.b=DiagnosticValues(o.at("b"));
    for(const auto & item:o.at("atoms").as_array())
    {
        const auto & a=item.as_object(); JointInitializationAtom atom;
        atom.id=Read<std::string>(a,"id");
        const auto ols=DiagnosticValues(a.at("ols")),mdpde=DiagnosticValues(a.at("mdpde"));
        if(ols.size()!=3 || mdpde.size()!=3) throw std::invalid_argument("Invalid joint initialization dimensions.");
        std::copy(ols.begin(),ols.end(),atom.ols.begin()); std::copy(mdpde.begin(),mdpde.end(),atom.mdpde.begin());
        atom.alpha=OptionalNumber(a.at("alpha")).value_or(std::numeric_limits<double>::quiet_NaN());
        atom.sample_count=Read<std::size_t>(a,"sample_count");
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
Object Units()
{
    return {{"contract","joint-kernel-map-units-v1"},{"map_value","fit-map-unit"},
        {"geometry","angstrom"},{"A","fit-map-unit*angstrom^3"},{"B","angstrom"},
        {"C","fit-map-unit*angstrom"}};
}
Object Metadata(const JointAnalysisMetadata & m)
{
    j::value normalization=nullptr,software=nullptr;
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
    return {{"model_path",m.model_path},{"map_path",m.map_path},
        {"grid_size",j::value_from(m.grid_size)},{"grid_spacing",j::value_from(m.grid_spacing)},
        {"origin",j::value_from(m.origin)},{"simulation",m.simulation},
        {"map_normalization",std::move(normalization)},{"model_sha256",OptionalText(m.model_sha256)},
        {"map_sha256",OptionalText(m.map_sha256)},{"software",std::move(software)},{"units",Units()}};
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
    for(double b:x.b) Require(b>0,"nonpositive width");
    Require(std::isfinite(x.objective),"nonfinite objective");
}
void Validate(const JointAnalysisResult & x)
{
    ValidateMetadata(x.metadata);
    const auto atoms=x.atom_ids.size(), rows=x.row_ids.size();
    Require(atoms>0 && std::set<std::string>(x.atom_ids.begin(),x.atom_ids.end()).size()==atoms,"atom identities");
    Require(std::set<std::string>(x.row_ids.begin(),x.row_ids.end()).size()==rows,"row identities");
    Require(x.available_row_mask.size()==rows,"row mask dimensions");
    Require(std::isfinite(x.observation_scale) && x.observation_scale>=1,"observation scale");
    Require(!x.objective || std::isfinite(*x.objective),"nonfinite objective");
    if(x.assembled_state) ValidateState(*x.assembled_state,atoms);
    if(x.initialization.valid)
    {
        Require(x.initialization.b.size()==atoms,"initial width dimensions");
        for(double b:x.initialization.b) Require(std::isfinite(b) && b>0,"invalid initial width");
    }
    Require(x.initialization.atoms.size()<=atoms,"initialization dimensions");
    std::set<std::size_t> used_atoms,used_rows;
    std::set<std::string> components;
    for(const auto & c:x.components)
    {
        Require(components.insert(c.id).second,"duplicate component");
        for(auto a:c.atoms) Require(a<atoms && used_atoms.insert(a).second,"component atom mapping");
        for(auto r:c.rows) Require(r<rows && used_rows.insert(r).second,"component row mapping");
        if(c.state) ValidateState(*c.state,c.atoms.size());
    }
}
}
std::string Encode(const JointAnalysisResult & x)
{
    Validate(x);
    const auto & m=x.metadata; const auto & t=x.costs;
    j::array components;
    for(const auto & c:x.components) components.emplace_back(Object{
        {"id",c.id},{"stop_reason",c.stop_reason},{"atoms",j::value_from(c.atoms)},{"rows",j::value_from(c.rows)},
        {"search_completed",c.search_completed},{"profile_evaluations",c.profile_evaluations},{"reference_evaluations",c.reference_evaluations},
        {"accepted_updates",c.accepted_updates},{"native_status",c.native_status},{"state",OptionalState(c.state)},
        {"evidence",Checks(c.evidence)},{"ranks",Ranks(c.ranks)},{"regular_certificate",StatusText(c.regular_certificate)},
        {"runtime_convergence",StatusText(c.runtime_convergence)}});
    Object out{{"schema_version",2},{"estimator","joint-components"},{"estimator_contract","guarded-joint-ls-v1"},{"objective_contract","parent-normalized-half-rss-v1"},
        {"support_contract","sphere-fma-v1"},{"metadata",Metadata(m)},
        {"atom_ids",j::value_from(x.atom_ids)},{"row_ids",j::value_from(x.row_ids)},{"initialization",Initialization(x.initialization)},
        {"costs",Object{{"initialization_seconds",t.initialization_seconds},{"search_seconds",t.search_seconds},
            {"search_reference_seconds",t.search_reference_seconds},{"assessment_seconds",t.assessment_seconds},{"assembly_seconds",t.assembly_seconds}}},
        {"components",std::move(components)},{"assembled_state",OptionalState(x.assembled_state)},{"objective",Number(x.objective)},
        {"available_row_mask",j::value_from(x.available_row_mask)},{"evidence",Checks(x.evidence)},{"ranks",Ranks(x.ranks)},
        {"search_completed",x.search_completed},{"observation_scale",x.observation_scale},
        {"regular_certificate",StatusText(x.regular_certificate)},{"runtime_convergence",StatusText(x.runtime_convergence)}};
    return j::serialize(out);
}
JointAnalysisResult Decode(std::string_view text)
{
    j::parse_options options; options.numbers=j::number_precision::precise;
    const auto parsed=j::parse(text,{},options); const auto & o=parsed.as_object();
    Require(Read<int>(o,"schema_version")==2,"unsupported result schema version (expected 2; regenerate older joint outcomes)");
    Require(Read<std::string>(o,"estimator")=="joint-components" &&
        Read<std::string>(o,"estimator_contract")=="guarded-joint-ls-v1" &&
        Read<std::string>(o,"objective_contract")=="parent-normalized-half-rss-v1" &&
        Read<std::string>(o,"support_contract")=="sphere-fma-v1","unsupported estimator contract");
    JointAnalysisResult x;
    x.metadata=ReadMetadata(o.at("metadata").as_object());
    x.atom_ids=Read<std::vector<std::string>>(o,"atom_ids"); x.row_ids=Read<std::vector<std::string>>(o,"row_ids");
    x.initialization=ReadInitialization(o.at("initialization"));
    const auto & t=o.at("costs").as_object();
    x.costs={Read<double>(t,"initialization_seconds"),Read<double>(t,"search_seconds"),Read<double>(t,"search_reference_seconds"),
        Read<double>(t,"assessment_seconds"),Read<double>(t,"assembly_seconds")};
    for(const auto & item:o.at("components").as_array())
    {
        const auto & c=item.as_object(); JointAnalysisComponent component;
        component.id=Read<std::string>(c,"id"); component.stop_reason=Read<std::string>(c,"stop_reason");
        component.atoms=Read<std::vector<std::size_t>>(c,"atoms"); component.rows=Read<std::vector<std::size_t>>(c,"rows");
        component.search_completed=Read<bool>(c,"search_completed"); component.profile_evaluations=Read<int>(c,"profile_evaluations");
        component.reference_evaluations=Read<int>(c,"reference_evaluations"); component.accepted_updates=Read<int>(c,"accepted_updates");
        component.native_status=Read<int>(c,"native_status");
        if(!c.at("state").is_null()) component.state=ReadState(c.at("state"));
        component.evidence=ReadChecks(c.at("evidence")); component.ranks=ReadRanks(c.at("ranks"));
        component.regular_certificate=Status(c.at("regular_certificate")); component.runtime_convergence=Status(c.at("runtime_convergence"));
        x.components.push_back(std::move(component));
    }
    if(!o.at("assembled_state").is_null()) x.assembled_state=ReadState(o.at("assembled_state"));
    x.objective=OptionalNumber(o.at("objective")); x.available_row_mask=Read<std::vector<bool>>(o,"available_row_mask");
    x.evidence=ReadChecks(o.at("evidence")); x.ranks=ReadRanks(o.at("ranks"));
    x.search_completed=Read<bool>(o,"search_completed"); x.observation_scale=Read<double>(o,"observation_scale");
    x.regular_certificate=Status(o.at("regular_certificate")); x.runtime_convergence=Status(o.at("runtime_convergence"));
    Validate(x); return x;
}
}
