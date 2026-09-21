#include "utils/domain/FileFingerprint.hpp"
#include "support/JointOfflineTools.hpp"
#include "support/JointComponentChecks.hpp"
#include "support/JointOfflineAudit.hpp"
#include "support/JointLocalEvidence.hpp"
#include "support/JointPrecisionAudit.hpp"
#include <filesystem>
#include <fstream>

namespace second_stage_test::matched::joint_abc {
namespace {
namespace j=boost::json;
namespace fs=std::filesystem;
using Vector=Eigen::VectorXd;
j::value Read(const fs::path & path)
{
    std::ifstream in(path); if(!in) throw std::runtime_error("Missing component input: "+path.string());
    j::parse_options options; options.numbers=j::number_precision::precise;
    return j::parse(std::string(std::istreambuf_iterator<char>(in),{}),{},options);
}
void Write(const fs::path & path,const j::value & value)
{std::ofstream out(path); out.exceptions(std::ios::failbit|std::ios::badbit); out<<j::serialize(value)<<'\n';}
Vector Parse(const j::value & v)
{Vector out(static_cast<Eigen::Index>(v.as_array().size())); for(Eigen::Index k=0;k<out.size();++k) out(k)=j::value_to<double>(v.at(static_cast<std::size_t>(k))); return out;}
j::array Values(const Vector & v) {j::array out; for(double x:v) out.push_back(x); return out;}
void AuditFit(const Domain & domain,const Vector & y,const j::object & search_fit,const EvaluationContext & context,
    const fs::path & output)
{
    fs::create_directories(output); auto fit=search_fit;
    // PrepareLocalAudit already assessed the actual state and local directions.
    fit["audit_assessment_refreshed"]=false;
    fit["audit_assessment_seconds"]=0.;
    Write(output/"fit.json",fit);
    certification::Audit(domain,y,fit,output/"audit.json",&context);
    auto audit=Read(output/"audit.json").as_object();
    const bool verified=fit.if_contains("derivative_verified") && fit.at("derivative_verified").as_bool();
    const bool requested=context.audit.expanded_if_unverified && (!verified || context.audit.precision);
    audit["scope_evidence"]=j::object{{"legacy_derivative_verified",verified},{"expanded_requested",requested},
        {"expanded_executed",audit.contains("ladders") && !audit.at("ladders").as_array().empty()},
        {"trigger",!requested ? "legacy-evidence-or-unregistered-expansion" : context.audit.precision ? "registered-high-precision-target" : "unverified-legacy-derivative"},
        {"precision_requested",context.audit.precision},{"boundary_requested",context.audit.boundary}};
    if(context.audit.precision && context.observations) audit["precision_normalization"]=certification::PrecisionNormalization(*context.observations);
    Write(output/"audit.json",audit);
    Write(output/"context.json",ContextEvidence(context));
    j::array direction_diagnostics;
    if(fit.contains("primary") && fit.at("primary").at("valid").as_bool())
    {
        const auto & state=fit.at("primary");
        const auto e=AtState(domain,y,Parse(state.at("eta")),Parse(state.at("beta")),context);
        const auto d=Differentiate(e,context.scale,&context);
        for(Eigen::Index k=0;k<context.audit.directions.cols();++k)
        {
            const double norm=context.audit.directions.col(k).norm();
            direction_diagnostics.push_back(j::object{{"direction",k},{"norm",norm},{"zero_direction",norm==0},
                {"jacobian_direction_norm",d.valid ? j::value((d.jacobian*context.audit.directions.col(k)).norm()) : j::value(nullptr)}});
        }
    }
    Write(output/"direction-diagnostics.json",direction_diagnostics);
    if(context.audit.boundary && fit.contains("primary") && fit.at("primary").at("valid").as_bool())
        Write(output/"boundary.json",certification::BoundaryAudit(domain,y,Parse(fit.at("primary").at("eta")),Parse(fit.at("primary").at("beta")),&context));
}
}
namespace {
void LocalAuditFit(const ComponentView & view,const Vector & parent_y,const j::object & fit,
    const EvaluationContext & parent,const fs::path & output)
{
    const auto y=Select(parent_y,view.rows);
    auto local=certification::PrepareLocalAudit(view.domain,y,fit,ComponentContext(parent,view,true));
    certification::ResetPrecisionCache();
    AuditFit(view.domain,y,local.fit,local.context,output);
    Write(output/"scope.json",local.scope);
}
}
void WriteLocalAudit(const ComponentView & view,const Vector & y,const j::object & fit,
    const EvaluationContext & context,const fs::path & output)
{LocalAuditFit(view,y,fit,context,output);}
void ComponentLocalBundleRerun(const std::string & bundle_path,const std::string & output_path)
{
    Eigen::setNbThreads(1); const auto bundle=Read(bundle_path); const fs::path output(output_path);
    if(fs::exists(output)) throw std::runtime_error("Use a fresh isolated component directory.");
    const auto & saved=bundle.at("search_context"),registration=saved.at("audit");
    const auto y=Parse(bundle.at("parent_observations"));
    AuditPlan plan;
    plan.trial_details=registration.at("trial_details").as_bool();
    plan.expanded_if_unverified=registration.at("expanded_if_unverified").as_bool();
    plan.precision=registration.at("precision").as_bool(); plan.boundary=registration.at("boundary").as_bool();
    plan.block_precision=registration.at("block_precision_reference").as_bool();
    plan.cache_precision=registration.at("precision_cache").as_string()=="exact-input, current-audit-case-only";
    for(const auto & atom:registration.at("boundary_atoms").as_array()) plan.boundary_atoms.push_back(j::value_to<Eigen::Index>(atom));
    const auto count=static_cast<Eigen::Index>(saved.at("atom_ids").as_array().size());
    const auto & directions=registration.at("directions").as_array();
    if(!directions.empty())
    {
        plan.directions.resize(count,static_cast<Eigen::Index>(directions.size()));
        for(std::size_t k=0;k<directions.size();++k)
        {
            const auto direction=Parse(directions[k]);
            if(direction.size()!=count) throw std::invalid_argument("Invalid parent audit direction.");
            plan.directions.col(static_cast<Eigen::Index>(k))=direction;
        }
    }
    auto context=runtime::CreateContext(y,count,j::value_to<std::string>(saved.at("parent_snapshot_sha256")),plan);
    std::vector<std::string> atom_ids,row_ids;
    for(const auto & id:saved.at("atom_ids").as_array()) atom_ids.push_back(j::value_to<std::string>(id));
    for(const auto & id:saved.at("row_ids").as_array()) row_ids.push_back(j::value_to<std::string>(id));
    context.atom_ids=std::move(atom_ids); context.row_ids=std::move(row_ids);
    if(ContextEvidence(context)!=saved) throw std::invalid_argument("Isolated parent context differs from its observations or numerical policy.");
    const auto & input=bundle.at("component_input"); ComponentView view;
    view.id=j::value_to<std::string>(input.at("id"));
    auto mappings=std::make_shared<runtime::PartitionMappings>(); view.mappings=mappings;
    mappings->atom_to_local.assign(static_cast<std::size_t>(count),-1); mappings->row_to_local.assign(static_cast<std::size_t>(y.size()),-1);
    mappings->atom_component.assign(static_cast<std::size_t>(count),-1); mappings->row_component.assign(static_cast<std::size_t>(y.size()),-1);
    for(const auto & atom:input.at("atoms").as_array())
    {
        const auto index=j::value_to<Eigen::Index>(atom);
        mappings->atom_to_local.at(static_cast<std::size_t>(index))=static_cast<Eigen::Index>(view.atoms.size()); mappings->atom_component.at(static_cast<std::size_t>(index))=0; view.atoms.push_back(index);
    }
    for(const auto & row:input.at("rows").as_array())
    {
        const auto index=j::value_to<Eigen::Index>(row);
        if(mappings->row_to_local.at(static_cast<std::size_t>(index))!=-1) throw std::invalid_argument("Duplicate component row.");
        mappings->row_to_local.at(static_cast<std::size_t>(index))=static_cast<Eigen::Index>(view.rows.size()); mappings->row_component.at(static_cast<std::size_t>(index))=0; view.rows.push_back(index);
    }
    std::vector<std::vector<Support>> support(view.atoms.size());
    for(const auto & entry:input.at("memberships").as_array())
        support.at(j::value_to<std::size_t>(entry.at(1))).push_back({j::value_to<Eigen::Index>(entry.at(0)),j::value_to<double>(entry.at(2))});
    view.domain=Domain(static_cast<Eigen::Index>(view.rows.size()),std::move(support));
    std::vector<std::string> ids; for(auto atom:view.atoms) ids.push_back(context.atom_ids.at(static_cast<std::size_t>(atom)));
    const auto partition=BuildPartition(view.domain,ids);
    if(partition.components.size()!=1 || partition.components[0].id!=view.id) throw std::invalid_argument("Bundle must contain exactly the requested structural component.");
    const auto local_b=Parse(bundle.at("component_initial_b"));
    if(local_b.size()!=static_cast<Eigen::Index>(view.atoms.size()) || !local_b.allFinite() || (local_b.array()<=0).any())
        throw std::invalid_argument("Invalid isolated component initialization.");
    Vector initial=Vector::Zero(count);
    for(std::size_t k=0;k<view.atoms.size();++k) initial(view.atoms[k])=local_b(static_cast<Eigen::Index>(k));
    auto fit=FitComponent(view,y,initial,context);
    fit["dataset"]=bundle.at("dataset"); fit["case"]=bundle.at("case"); fit["initial_b"]=Values(local_b);
    fit["observation_snapshot_sha256"]=context.snapshot_hash;
    fs::create_directories(output); Write(output/"fit.json",fit);
    LocalAuditFit(view,y,fit,context,output/"audit");
    Write(output/"completion.json",j::object{{"complete",true},{"component_id",view.id},
        {"snapshot_sha256",context.snapshot_hash},{"input_bundle_sha256",rhbm_gem::FileSha256(bundle_path)}});
}
}
