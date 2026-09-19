#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "core/detail/joint_component/Problem.hpp"
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
#include <limits>
#include <set>
#include <stdexcept>

namespace rhbm_gem::core {
namespace {
namespace n=joint_component;
using Clock=std::chrono::steady_clock;
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
std::vector<double> Values(const n::Vector & v) {return {v.data(),v.data()+v.size()};}
std::vector<const AtomObject *> Contributors(const ModelObject & model)
{
    std::set<const AtomObject *> selected;
    for(const auto * atom:model.GetSelectedAtoms()) if(atom->GetElement()!=Element::HYDROGEN) selected.insert(atom);
    std::vector<const AtomObject *> atoms;
    for(const auto & atom:model.GetAtomList()) if(atom->GetElement()!=Element::HYDROGEN)
    {
        if(!selected.contains(atom.get())) throw std::invalid_argument("Joint fitting requires all non-hydrogen contributors; partial selection is unsupported.");
        atoms.push_back(atom.get());
    }
    if(atoms.empty()) throw std::invalid_argument("Joint fitting requires at least one non-hydrogen atom.");
    return atoms;
}
JointState State(const n::Endpoint & e)
{return {Values(e.beta),Values(e.eta.array().exp()),Values(e.eta),Values(e.gradient),e.certificate.objective};}
std::vector<JointCheck> Evidence(const n::Assessment & a,JointEvidenceScope scope)
{
    using Status=JointCheckStatus;
    auto check=[scope](std::string name,bool available,bool pass,std::optional<double> value={},std::optional<double> threshold={}) {
        return JointCheck{std::move(name),available ? (pass ? Status::Passed : Status::Failed) : Status::Unavailable,scope,value,threshold,{}};
    };
    const bool valid=a.primary.valid && a.reference.valid;
    std::vector<JointCheck> out;
    out.push_back(check("inner",a.design.has_value(),a.design && a.design->rank==a.primary.beta.size() && a.coefficient_difference<=1e-10,
        std::isfinite(a.coefficient_difference) ? std::optional<double>(a.coefficient_difference) : std::nullopt,1e-10));
    out.push_back(check("kkt",a.primary.certificate.available,a.primary.certificate.kkt_passed,
        a.primary.certificate.available ? std::optional<double>(a.primary.certificate.projected_kkt) : std::nullopt,1e-10));
    out.push_back(check("width-stationarity",valid && a.widths.has_value(),a.gradient,
        valid && a.primary.gradient.size() ? std::optional<double>(std::max(a.primary.gradient.lpNorm<Eigen::Infinity>(),
            a.reference.gradient.lpNorm<Eigen::Infinity>())) : std::nullopt,1e-12));
    out.push_back(check("local-correction",a.jacobian.has_value(),a.local,
        a.correction.size() ? std::optional<double>(a.correction.lpNorm<Eigen::Infinity>()) : std::nullopt,1e-10));
    out.push_back(check("numerical-identifiability",a.widths.has_value(),a.identified));
    double derivative_error{}; bool same_face=true;
    for(const auto & d:a.derivatives)
    {derivative_error=std::isfinite(d.error) ? std::max(derivative_error,d.error) : std::numeric_limits<double>::infinity(); same_face &= d.same_face;}
    out.push_back(check("two-step-derivative",!a.derivatives.empty(),a.derivative_verified,
        !a.derivatives.empty() && std::isfinite(derivative_error) ? std::optional<double>(derivative_error) : std::nullopt,1e-6));
    if(!a.derivatives.empty() && !same_face) out.back().reason="active-face-not-preserved";
    for(const char * name:{"richardson","precision-50-100","boundary-audit","regular-certificate"})
        out.push_back({name,Status::NotRun,scope,{},{},"offline-audit-not-run"});
    for(auto & evidence:out) if(evidence.status==Status::Unavailable)
        evidence.reason=a.failure.empty() ? "missing-trusted-state" : a.failure;
    return out;
}
std::vector<JointRankEvidence> Ranks(const n::Assessment & a,JointEvidenceScope scope)
{
    std::vector<JointRankEvidence> out;
    for(const auto & [name,spectrum]:std::vector<std::pair<std::string,const std::optional<n::Spectrum> *>>{
        {"design",&a.design},{"projected-width",&a.widths},{"normalized-width",&a.normalized_widths},{"jacobian",&a.jacobian}})
        if(*spectrum) out.push_back({name,scope,static_cast<std::size_t>((*spectrum)->rank),(*spectrum)->threshold,Values((*spectrum)->singular_values)});
    return out;
}
}
JointProblem::JointProblem(JointProblemInput input)
{
    if(input.atom_ids.empty() || input.support.size()!=input.atom_ids.size() || input.row_ids.size()!=input.observations.size() ||
        std::set<std::string>(input.row_ids.begin(),input.row_ids.end()).size()!=input.row_ids.size())
        throw std::invalid_argument("Invalid joint problem identities or dimensions.");
    auto data=std::make_shared<n::ProblemData>(); data->input=std::move(input);
    data->y=Eigen::Map<const n::Vector>(data->input.observations.data(),static_cast<Eigen::Index>(data->input.observations.size()));
    if(!data->y.allFinite()) throw std::invalid_argument("Nonfinite joint observations.");
    data->domain.rows=data->y.size(); data->domain.atoms.resize(data->input.support.size());
    for(std::size_t a=0;a<data->input.support.size();++a) for(const auto & s:data->input.support[a])
    {
        if(s.row>=data->input.observations.size()) throw std::invalid_argument("Invalid joint contributor row.");
        data->domain.atoms[a].push_back({static_cast<Eigen::Index>(s.row),s.squared_distance});
    }
    data->partition=n::Partition(data->domain,data->input.atom_ids);
    data->context=n::CreateContext(data->y,static_cast<Eigen::Index>(data->input.atom_ids.size()));
    data->context.atom_ids=data->input.atom_ids; data->context.row_ids=data->input.row_ids;
    m_data=std::move(data);
}
const JointProblemInput & JointProblem::Input() const {return m_data->input;}
double JointProblem::ObservationScale() const {return m_data->context.scale;}
JointProblem BuildJointProblem(const MapObject & map,const ModelObject & model)
{
    const auto atoms=Contributors(model); JointProblemInput input;
    const auto dims=map.GetGridSize(); const auto spacing=map.GetGridSpacing(),origin=map.GetOrigin();
    for(std::size_t k=0;k<3;++k) if(dims[k]<=0 || !std::isfinite(spacing[k]) || spacing[k]<=0 || !std::isfinite(origin[k]))
        throw std::invalid_argument("Invalid joint map geometry.");
    std::vector<std::vector<std::pair<std::size_t,double>>> support(atoms.size());
    std::vector<Eigen::Index> rows(map.GetMapValueArraySize(),-1);
    for(std::size_t a=0;a<atoms.size();++a)
    {
        input.atom_ids.push_back(std::to_string(atoms[a]->GetSerialID()));
        const auto position=atoms[a]->GetPosition(); std::array<int,3> lo{},hi{};
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
            if(square<=6.25) {rows[index]=0; support[a].push_back({index,square});}
        }
    }
    for(std::size_t k=0;k<rows.size();++k) if(rows[k]==0)
    {rows[k]=static_cast<Eigen::Index>(input.observations.size()); input.row_ids.push_back(std::to_string(k)); input.observations.push_back(map.GetMapValue(k));}
    input.support.resize(atoms.size());
    for(std::size_t a=0;a<atoms.size();++a) for(const auto & [index,square]:support[a])
        input.support[a].push_back({static_cast<std::size_t>(rows[index]),square});
    return JointProblem(std::move(input));
}
JointFitResult FitJointComponents(const JointProblem & problem,const std::vector<double> & initial_b)
{
    eigen_helper::ScopedEigenThreadCount eigen_thread_guard{1};
    const auto & data=JointProblemAccess::Get(problem); JointFitResult out; out.problem=problem;
    out.observation_scale=data.context.scale; out.initialization.b=initial_b;
    out.initialization.valid=initial_b.size()==data.domain.atoms.size() && std::all_of(initial_b.begin(),initial_b.end(),[](double b){return std::isfinite(b) && b>0;});
    out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
    if(!out.initialization.valid)
    {
        out.available_row_mask.assign(data.input.observations.size(),false);
        for(auto row:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(row)]=true;
        return out;
    }
    const n::Vector b=Eigen::Map<const n::Vector>(initial_b.data(),static_cast<Eigen::Index>(initial_b.size()));
    std::vector<n::ComponentResult> results;
    for(const auto & view:data.partition.components)
    {
        const auto component_start=Clock::now();
        auto result=n::SolveComponent(view,data.y,b,data.context);
        JointComponentResult component; component.id=view.id; component.stop_reason=result.search.stop_reason;
        for(auto a:view.atoms) component.atoms.push_back(static_cast<std::size_t>(a));
        for(auto row:view.rows) component.rows.push_back(static_cast<std::size_t>(row));
        component.search_completed=result.search_success; component.native_status=result.search.lm_status;
        component.profile_evaluations=result.search.evaluations; component.reference_evaluations=result.search.references;
        component.accepted_updates=result.search.accepted;
        if(result.trusted_state)
        {
            component.state=State(*result.trusted_state);
            auto context=n::ChildContext(data.context,view,true); context.audit.directions.resize(0,0);
            const auto assessment=n::AssessProfile(view.domain,n::SelectValues(data.y,view.rows),result.trusted_state->eta,context,&result.trusted_state->beta);
            component.evidence=Evidence(assessment,JointEvidenceScope::ComponentLocal);
            component.ranks=Ranks(assessment,JointEvidenceScope::ComponentLocal);
        }
        else component.evidence=Evidence({},JointEvidenceScope::ComponentLocal);
        out.costs.search_seconds+=result.search.seconds;
        out.costs.search_reference_seconds+=result.search.reference_seconds;
        out.costs.assessment_seconds+=Seconds(component_start)-result.search.seconds;
        out.components.push_back(std::move(component)); results.push_back(std::move(result));
    }
    const auto assembly_start=Clock::now();
    const auto assembly=n::AssembleComponents(data.domain,data.y,data.partition,data.context,results);
    out.search_completed=assembly.completed; out.available_row_mask=assembly.row_mask;
    out.evidence=Evidence(assembly.assessment,JointEvidenceScope::AssembledGlobal);
    out.ranks=Ranks(assembly.assessment,JointEvidenceScope::AssembledGlobal);
    if(assembly.available)
    {
        const auto state=n::EvaluateState(data.domain,data.y,assembly.eta,assembly.beta,data.context);
        out.assembled_state=State(state);
        if(state.valid) {out.prediction=Values(assembly.prediction); out.objective=assembly.objective;}
        out.evidence.push_back({"assembled-profile",assembly.profile_agrees ? JointCheckStatus::Passed : JointCheckStatus::Failed,
            JointEvidenceScope::AssembledGlobal,assembly.profile_difference,1e-10,{}});
    }
    out.costs.assembly_seconds=Seconds(assembly_start);
    return out;
}
JointFitResult EstimateJointComponents(MapObject & map,ModelObject & model)
{
    const auto problem=BuildJointProblem(map,model);
    ModelObject initializer(model); const auto atoms=Contributors(initializer);
    JointInitialization initialization; const auto initialization_start=Clock::now();
    try
    {
        initializer.ApplyElementSelection(Element::HYDROGEN,true);
        initializer.EditAnalysis().InitializeFromSelection();
        RunPotentialSamplingWorkflow(map,initializer,SphereSamplingMethod::FibonacciDeterministic,1);
        initializer.EditAnalysis().InitializeLocalFittingSeedModels();
        FitOptions options; options.thread_size=1; options.quiet_mode=true; options.exclude_hydrogen=true;
        RunLocalAlphaTraining(initializer,options,FittingStage::First);
        RunFixedOffsetLocalFitting(initializer,options,FittingStage::First);
        for(const auto * atom:atoms)
        {
            const auto view=AtomLocalPotentialView::For(*atom); const auto & local=view.GetGaussianResult(FittingStage::First);
            JointInitializationAtom record; record.id=std::to_string(atom->GetSerialID());
            const auto ols=local.ols.GetModel().ToVector(),mdpde=local.mdpde.GetModel().ToVector();
            for(std::size_t k=0;k<3;++k) {record.ols[k]=ols(static_cast<Eigen::Index>(k)); record.mdpde[k]=mdpde(static_cast<Eigen::Index>(k));}
            record.alpha=local.alpha_r; record.sample_count=view.GetSamplingEntries(FittingStage::First).size();
            if(local.fit_result) record.native_status=static_cast<int>(local.fit_result->status);
            initialization.b.push_back(local.mdpde.GetModel().GetWidth()); initialization.atoms.push_back(std::move(record));
        }
    }
    catch(const std::exception & error)
    {
        JointFitResult out; out.problem=problem; initialization.reason=error.what(); out.initialization=std::move(initialization);
        out.costs.initialization_seconds=Seconds(initialization_start);
        out.observation_scale=problem.ObservationScale(); out.available_row_mask.assign(problem.Input().observations.size(),false); return out;
    }
    for(const auto * atom:atoms)
    {
        const auto view=AtomLocalPotentialView::For(*atom);
        const auto * destination=model.FindAtomPtr(atom->GetSerialID());
        model.EditAnalysis().SetAtomLocalRawSamplingEntries(*destination,view.GetRawSamplingEntries(false));
        model.EditAnalysis().SetAtomLocalGaussianResult(FittingStage::First,*destination,view.GetGaussianResult(FittingStage::First));
    }
    const double initialization_seconds=Seconds(initialization_start);
    auto out=FitJointComponents(problem,initialization.b); out.costs.initialization_seconds=initialization_seconds;
    out.initialization.atoms=std::move(initialization.atoms); return out;
}
}
