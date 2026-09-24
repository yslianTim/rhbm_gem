#include "SimulationBuildInfo.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/utils/math/EigenHelper.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <set>
#include <string_view>
#include <stdexcept>

namespace rhbm_gem::core {
namespace {
namespace n=joint_component;
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
std::vector<double> Values(const n::Vector & v) {return {v.data(),v.data()+v.size()};}
std::vector<const AtomObject *> EligibleAtoms(const ModelObject & model)
{
    std::vector<const AtomObject *> atoms;
    for(const auto & atom:model.GetAtomList()) if(atom->GetElement()!=Element::HYDROGEN) atoms.push_back(atom.get());
    return atoms;
}
template<class Visit> void VisitSupport(const MapObject & map,const AtomObject & atom,Visit visit)
{
    const auto dims=map.GetGridSize(); const auto spacing=map.GetGridSpacing(),origin=map.GetOrigin();
    const auto position=atom.GetPosition(); std::array<int,3> lo{},hi{};
    for(std::size_t k=0;k<3;++k)
    {
        if(!std::isfinite(position[k])) throw std::invalid_argument("Nonfinite atom position.");
        lo[k]=static_cast<int>(std::clamp(std::floor((position[k]-2.5-origin[k])/spacing[k]),0.,static_cast<double>(dims[k])));
        hi[k]=static_cast<int>(std::clamp(std::floor((position[k]+2.5-origin[k])/spacing[k]),-1.,static_cast<double>(dims[k]-1)));
    }
    for(int z=lo[2];z<=hi[2];++z) for(int y=lo[1];y<=hi[1];++y) for(int x=lo[0];x<=hi[0];++x)
    {
        const auto index=static_cast<std::size_t>(x)+static_cast<std::size_t>(dims[0])*(static_cast<std::size_t>(y)+static_cast<std::size_t>(dims[1])*static_cast<std::size_t>(z));
        const double square=simulation::SupportSquare(simulation::GridPosition({x,y,z},spacing,origin),position);
        if(square<=6.25) visit(index,square);
    }
}
JointState State(const n::Endpoint & e,double scale)
{return {Values(e.beta),Values(e.eta.array().exp()),Values(e.eta),Values(e.gradient),e.certificate.objective/(scale*scale)};}
}
namespace joint_component {
std::vector<JointCheck> AssessmentEvidence(const Assessment & a,JointEvidenceScope scope)
{
    using Status=JointCheckStatus;
    auto check=[scope](std::string name,bool available,bool pass,std::optional<double> value={},std::optional<double> threshold={}) {
        return JointCheck{std::move(name),available ? (pass ? Status::Passed : Status::Failed) : Status::Unavailable,scope,value,threshold,{}};
    };
    const bool valid=a.primary.valid && a.reference.valid;
    std::vector<JointCheck> out;
    out.push_back(check("inner",a.design.has_value(),a.inner,
        std::isfinite(a.coefficient_difference) ? std::optional<double>(a.coefficient_difference) : std::nullopt,1e-10));
    out.push_back(check("kkt",a.primary.certificate.available,a.primary.certificate.kkt_passed,
        a.primary.certificate.available ? std::optional<double>(a.primary.certificate.projected_kkt) : std::nullopt,1e-10));
    out.push_back(check("width-stationarity",valid && a.widths.has_value(),a.gradient,
        valid && a.primary.gradient.size() ? std::optional<double>(std::max(a.primary.gradient.lpNorm<Eigen::Infinity>(),
            a.reference.gradient.lpNorm<Eigen::Infinity>())) : std::nullopt,1e-12));
    out.push_back(check("local-correction",a.jacobian.has_value(),a.local,
        a.correction.size() ? std::optional<double>(a.correction.lpNorm<Eigen::Infinity>()) : std::nullopt,1e-10));
    out.push_back(check("numerical-identifiability",a.widths.has_value(),a.identified));
    for(const char * name:{"two-step-derivative","richardson","precision-50-100","boundary-audit","regular-certificate"})
        out.push_back({name,Status::NotRun,scope,{},{},"offline-audit-not-run"});
    for(auto & evidence:out) if(evidence.status==Status::Unavailable)
        evidence.reason=a.failure.empty() ? "missing-trusted-state" : a.failure;
    return out;
}
JointCheckStatus MergeConvergenceStatus(JointCheckStatus a,JointCheckStatus b)
{
    for(const auto status:{JointCheckStatus::Failed,JointCheckStatus::Unavailable,JointCheckStatus::NotRun})
        if(a==status || b==status) return status;
    return JointCheckStatus::Passed;
}
JointCheckStatus ConvergenceStatus(const std::vector<JointCheck> & evidence,JointEvidenceScope scope,bool assembled)
{
    auto status=JointCheckStatus::Passed;
    for(const char * name:{"inner","kkt","width-stationarity","local-correction","numerical-identifiability","assembled-profile"})
    {
        if(!assembled && std::string_view(name)=="assembled-profile") continue;
        const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & check){return check.name==name && check.scope==scope;});
        status=MergeConvergenceStatus(status,found==evidence.end() ? JointCheckStatus::Unavailable : found->status);
    }
    return status;
}
}
JointCheckStatus JointComponentResult::RuntimeConvergence() const
{
    return state ? joint_component::ConvergenceStatus(evidence,JointEvidenceScope::ComponentLocal) : JointCheckStatus::Unavailable;
}
JointCheckStatus JointFitResult::RuntimeConvergence() const
{
    if(!assembled_state || !prediction || !objective || components.empty() ||
        std::any_of(components.begin(),components.end(),[](const auto & c){return !c.state;}) ||
        !std::all_of(available_row_mask.begin(),available_row_mask.end(),[](bool available){return available;}))
        return JointCheckStatus::Unavailable;
    auto status=joint_component::ConvergenceStatus(evidence,JointEvidenceScope::AssembledGlobal,true);
    for(const auto & component:components) status=joint_component::MergeConvergenceStatus(status,component.RuntimeConvergence());
    return status;
}
std::vector<JointRankEvidence> joint_component::AssessmentRanks(const n::Assessment & a,JointEvidenceScope scope)
{
    std::vector<JointRankEvidence> out;
    for(const auto & [name,spectrum]:std::vector<std::pair<std::string,const std::optional<n::Spectrum> *>>{
        {"design",&a.design},{"projected-width",&a.widths},{"normalized-width",&a.normalized_widths},{"jacobian",&a.jacobian}})
        if(*spectrum) out.push_back({name,scope,static_cast<std::size_t>((*spectrum)->rank),(*spectrum)->threshold,Values((*spectrum)->singular_values)});
    return out;
}
JointProblem::JointProblem(JointProblemInput input)
{
    if(input.atom_ids.empty() || input.support.size()!=input.atom_ids.size() || input.row_ids.size()!=input.observations.size() ||
        std::set<std::string>(input.row_ids.begin(),input.row_ids.end()).size()!=input.row_ids.size())
        throw std::invalid_argument("Invalid joint problem identities or dimensions.");
    if(input.selection_domain && !input.selection_domain->IsValid(input.atom_ids.size()))
        throw std::invalid_argument("Invalid joint selection domain.");
    for(const auto & atom:input.support) for(const auto & s:atom)
        if(s.row>=input.observations.size()) throw std::invalid_argument("Invalid joint contributor row.");
    auto data=std::make_shared<n::ProblemData>(std::make_shared<const JointProblemInput>(std::move(input)));
    if(!data->y.allFinite()) throw std::invalid_argument("Nonfinite joint observations.");
    data->partition=n::Partition(data->domain,data->input->atom_ids);
    data->context=n::CreateContext(data->input);
    data->layout=n::BuildParameterLayout(*data->input);
    m_data=std::move(data);
}
const JointProblemInput & JointProblem::Input() const {return *m_data->input;}
double JointProblem::ObservationScale() const {return m_data->context.scale;}
const JointParameterLayout & JointProblem::ParameterLayout() const {return m_data->layout;}
JointProblem BuildJointProblem(const MapObject & map,const ModelObject & model)
{
    const auto atoms=EligibleAtoms(model); JointProblemInput input;
    std::set<const AtomObject *> selected;
    for(const auto * atom:model.GetSelectedAtoms()) if(atom->GetElement()!=Element::HYDROGEN) selected.insert(atom);
    if(selected.empty()) throw std::invalid_argument("Joint fitting requires at least one selected non-hydrogen atom.");
    const auto dims=map.GetGridSize(); const auto spacing=map.GetGridSpacing(),origin=map.GetOrigin();
    for(std::size_t k=0;k<3;++k) if(dims[k]<=0 || !std::isfinite(spacing[k]) || spacing[k]<=0 || !std::isfinite(origin[k]))
        throw std::invalid_argument("Invalid joint map geometry.");
    std::vector<std::size_t> rows;
    for(const auto * atom:atoms) if(selected.contains(atom))
        VisitSupport(map,*atom,[&](auto index,double) {rows.push_back(index);});
    std::sort(rows.begin(),rows.end()); rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    for(auto index:rows) {input.row_ids.push_back(std::to_string(index)); input.observations.push_back(map.GetMapValue(index));}
    input.selection_domain.emplace();
    for(const auto * atom:atoms)
    {
        std::vector<JointSupport> support;
        VisitSupport(map,*atom,[&](auto index,double square) {
            const auto row=std::lower_bound(rows.begin(),rows.end(),index);
            if(row!=rows.end() && *row==index) support.push_back({static_cast<std::size_t>(row-rows.begin()),square});
        });
        if(selected.contains(atom) || !support.empty())
        {
            if(selected.contains(atom)) input.selection_domain->target_indices.push_back(input.atom_ids.size());
            input.atom_ids.push_back(std::to_string(atom->GetSerialID()));
            input.support.push_back(std::move(support));
        }
    }
    return JointProblem(std::move(input));
}
JointFitResult FitJointComponents(const JointProblem & problem,const std::vector<double> & initial_b)
{return n::FitWithSearchPolicy(problem,initial_b,{});}
JointFitResult n::FitWithSearchPolicy(const JointProblem & problem,const std::vector<double> & initial_b,const n::SearchPolicy & search_policy)
{
    eigen_helper::ScopedEigenThreadCount eigen_thread_guard{1};
    if(!problem.ParameterLayout().groups.empty()) return n::FitObservableComponents(problem,initial_b,search_policy);
    const auto & data=JointProblemAccess::Get(problem); JointFitResult out; out.problem=problem;
    out.observation_scale=data.context.scale; out.initialization.b=initial_b;
    out.initialization.valid=initial_b.size()==data.domain.atoms.size() && std::all_of(initial_b.begin(),initial_b.end(),[](double b){return std::isfinite(b) && b>0;});
    out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
    if(initial_b.size()!=data.domain.atoms.size())
    {
        out.available_row_mask.assign(data.input->observations.size(),false);
        for(auto row:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(row)]=true;
        return out;
    }
    auto search_context=data.context; search_context.search=search_policy;
    const n::Vector b=Eigen::Map<const n::Vector>(initial_b.data(),static_cast<Eigen::Index>(initial_b.size()));
    std::vector<n::ComponentResult> results;
    for(const auto & view:data.partition.components)
    {
        const auto component_start=Clock::now();
        n::ComponentResult result;
        const bool valid=std::all_of(view.atoms.begin(),view.atoms.end(),[&](auto a) {
            return std::isfinite(initial_b[static_cast<std::size_t>(a)]) && initial_b[static_cast<std::size_t>(a)]>0;
        });
        if(valid) result=n::SolveComponent(view,data.y,b,search_context);
        else {result.search.stopped=true; result.search.stop_reason="invalid-initial-widths";}

        JointComponentResult component; component.id=view.id; component.stop_reason=result.search.stop_reason;
        for(auto a:view.atoms) component.atoms.push_back(static_cast<std::size_t>(a));
        for(auto row:view.rows) component.rows.push_back(static_cast<std::size_t>(row));
        component.search_completed=result.search_success; component.native_status=result.search.lm_status;
        component.profile_evaluations=result.search.evaluations; component.reference_evaluations=result.search.references;
        component.accepted_updates=result.search.accepted;
        if(result.trusted_state)
        {
            component.state=State(*result.trusted_state,data.context.scale);
            component.evidence=n::AssessmentEvidence(*result.trusted_assessment,JointEvidenceScope::ComponentLocal);
            component.ranks=n::AssessmentRanks(*result.trusted_assessment,JointEvidenceScope::ComponentLocal);
        }
        else
        {
            n::Assessment missing;
            missing.failure=result.search.stop_reason;
            component.evidence=n::AssessmentEvidence(missing,JointEvidenceScope::ComponentLocal);
        }
        out.costs.search_seconds+=result.search.seconds;
        out.costs.search_reference_seconds+=result.search.reference_seconds;
        out.costs.assessment_seconds+=Seconds(component_start)-result.search.seconds;
        out.components.push_back(std::move(component)); results.push_back(std::move(result));
    }
    const auto assembly_start=Clock::now();
    std::optional<n::EvaluationContext> reuse_context;
    std::optional<n::AssessmentReuse> reuse;
    if(results.size()==1 && results[0].trusted_assessment && data.partition.constant_rows.empty())
    {
        const auto & view=data.partition.components[0];
        // A single full component has identity row/atom mappings in this immutable snapshot.
        reuse_context=n::ChildContext(data.context,view,true);
        reuse.emplace(n::AssessmentReuse{view.domain,data.y,*reuse_context,*results[0].trusted_assessment});
    }
    const auto assembly=n::AssembleComponents(data.domain,data.y,data.partition,data.context,results,reuse ? &*reuse : nullptr);
    out.search_completed=assembly.completed; out.available_row_mask=assembly.row_mask;
    out.evidence=n::AssessmentEvidence(assembly.assessment,JointEvidenceScope::AssembledGlobal);
    out.ranks=n::AssessmentRanks(assembly.assessment,JointEvidenceScope::AssembledGlobal);
    if(assembly.available)
    {
        const auto & state=assembly.raw;
        out.assembled_state=State(state,data.context.scale);
        if(state.valid) {out.prediction=Values(assembly.prediction); out.objective=assembly.objective/(data.context.scale*data.context.scale);}
        out.evidence.push_back({"assembled-profile",assembly.profile_agrees ? JointCheckStatus::Passed : JointCheckStatus::Failed,
            JointEvidenceScope::AssembledGlobal,assembly.profile_difference,1e-10,{}});
    }
    out.costs.assembly_seconds=Seconds(assembly_start);
    return out;
}
JointAnalysisResult CaptureJointAnalysisResult(const JointFitResult & fit, JointAnalysisMetadata metadata)
{
    if (!fit.problem) throw std::invalid_argument("Cannot capture joint result without problem identities.");
    JointAnalysisResult out;
    out.layout=fit.layout;
    if(out.layout) out.parameterization_contract="singleton-halo-profile-v1";
    metadata.software=JointSoftwareProvenance{RHBM_GEM_SIMULATION_VERSION,
        RHBM_GEM_SIMULATION_SOURCE_SHA256,RHBM_GEM_SIMULATION_CONFIG_SHA256,
        RHBM_GEM_SIMULATION_BUILD_SHA256};
    out.metadata=std::move(metadata);
    out.atom_ids=fit.problem->Input().atom_ids; out.row_ids=fit.problem->Input().row_ids;
    out.selection_domain=fit.problem->Input().selection_domain;
    out.initialization=fit.initialization; out.costs=fit.costs;
    for (const auto & component:fit.components)
    {
        JointAnalysisComponent saved;
        static_cast<JointComponentData &>(saved)=component;
        saved.runtime_convergence=component.RuntimeConvergence();
        out.components.push_back(std::move(saved));
    }
    out.assembled_state=fit.assembled_state; out.objective=fit.objective;
    out.available_row_mask=fit.available_row_mask; out.evidence=fit.evidence; out.ranks=fit.ranks;
    out.search_completed=fit.search_completed; out.observation_scale=fit.observation_scale;
    out.runtime_convergence=fit.RuntimeConvergence(); out.regular_certificate=fit.regular_certificate;
    return out;
}
JointInitialization joint_component::InitializeContributors(MapObject & map,ModelObject & model,const JointProblem & problem)
{
    ModelObject initializer(model);
    const auto workset = detail::MakeJointFittingWorkset(initializer, problem);
    FitOptions options; options.thread_size = 1; options.quiet_mode = true; options.exclude_hydrogen = true;
    const auto initialization = detail::RunFirstStage(initializer, workset, options, detail::FirstStageMode::SampleContributorsIsolated, &map);
    for (std::size_t i = 0; i < workset.contributors.size(); ++i)
    {
        if (!workset.target_mask[i] || initialization.atoms[i].reason != "valid-width") continue;
        const auto * source = workset.contributors[i];
        const auto * target = model.FindAtomPtr(source->GetSerialID());
        const auto view = AtomLocalPotentialView::For(*source);
        model.EditAnalysis().SetAtomLocalRawSamplingEntries(*target, view.GetRawSamplingEntries(false));
        model.EditAnalysis().SetAtomLocalGaussianResult(FittingStage::First, *target,
            view.GetGaussianResult(FittingStage::First));
    }
    return initialization;
}
JointFitResult EstimateJointComponents(MapObject & map,ModelObject & model)
{
    const auto construction_start=Clock::now();
    const auto problem=BuildJointProblem(map,model);
    const double construction_seconds=Seconds(construction_start);
    const auto initialization_start=Clock::now();
    auto initialization=n::InitializeContributors(map,model,problem);
    const double initialization_seconds=Seconds(initialization_start);
    auto out=FitJointComponents(problem,initialization.b); out.costs.initialization_seconds=initialization_seconds;
    out.costs.construction_seconds=construction_seconds;
    out.initialization=std::move(initialization);
    return out;
}
}
