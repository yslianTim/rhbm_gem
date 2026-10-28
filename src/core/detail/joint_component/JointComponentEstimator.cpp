#include "SimulationBuildInfo.hpp"
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/WidthCalibration.hpp"
#include "core/detail/joint_component/ParameterStability.hpp"
#include "core/detail/FirstStageInitialization.hpp"
#include "core/detail/FittingWorkset.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
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
#include <numeric>
#include <set>
#include <string_view>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace rhbm_gem::core {
namespace {
namespace n=joint_component;
using Clock=std::chrono::steady_clock;
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
thread_local std::optional<std::size_t> problem_memory_budget_for_testing;
#endif
std::size_t JointProblemMemoryBudget()
{
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    if(problem_memory_budget_for_testing) return *problem_memory_budget_for_testing;
#endif
    return n::JointMemoryPlanningBudgetBytes;
}
double Seconds(Clock::time_point start) {return std::chrono::duration<double>(Clock::now()-start).count();}
std::vector<double> Values(const n::Vector & v) {return {v.data(),v.data()+v.size()};}
bool ValidJointState(const JointState & state,std::size_t atoms)
{
    if(atoms>std::numeric_limits<std::size_t>::max()/2 || state.ac.size()!=2*atoms ||
        state.b.size()!=atoms || state.log_b.size()!=atoms || state.width_gradient.size()!=atoms ||
        !std::isfinite(state.objective) || state.objective<0 ||
        !std::all_of(state.ac.begin(),state.ac.end(),[](double value){return std::isfinite(value);}) ||
        !std::all_of(state.b.begin(),state.b.end(),[](double value){return std::isfinite(value) && value>0;}) ||
        !std::all_of(state.log_b.begin(),state.log_b.end(),[](double value){return std::isfinite(value);}) ||
        !std::all_of(state.width_gradient.begin(),state.width_gradient.end(),[](double value){return std::isfinite(value);}))
        return false;
    for(std::size_t atom=0;atom<atoms;++atom)
    {
        const double amplitude=state.ac[2*atom],width=state.b[atom];
        if(amplitude<0) return false;
        const double saved_width=std::exp(state.log_b[atom]);
        if(!std::isfinite(saved_width) || !(saved_width>0) ||
            std::abs(saved_width-width)>8*std::numeric_limits<double>::epsilon()*width) return false;
    }
    return true;
}
void RemapAtomIndex(std::size_t & value,const std::vector<std::size_t> & mapping)
{
    if(value>=mapping.size()) throw std::invalid_argument("Invalid fitting-to-original atom mapping.");
    value=mapping[value];
}
void RemapRowIndex(std::size_t & value,const std::vector<std::size_t> & mapping)
{
    if(value>=mapping.size()) throw std::invalid_argument("Invalid fitting-to-original row mapping.");
    value=mapping[value];
}
void RemapTargetEvidence(std::optional<JointTargetEvidence> & evidence,
    const std::vector<std::size_t> & atom_mapping,std::size_t original_rows)
{
    if(!evidence) return;
    for(auto & point:evidence->point_atoms) RemapAtomIndex(point.atom,atom_mapping);
    evidence->original_rows=original_rows;
}
void RemapComponentForOriginalDomain(JointAnalysisComponent & component,
    const std::vector<std::size_t> & atom_mapping,const std::vector<std::size_t> & row_mapping)
{
    for(auto & atom:component.atoms) RemapAtomIndex(atom,atom_mapping);
    for(auto & row:component.rows) RemapRowIndex(row,row_mapping);
    if(component.layout)
    {
        for(auto & atom:component.layout->full_atoms) RemapAtomIndex(atom,atom_mapping);
        for(auto & row:component.layout->informative_rows) RemapRowIndex(row,row_mapping);
    }
    for(auto & block:component.local_conditioning)
        for(auto & atom:block.core_atoms) RemapAtomIndex(atom,atom_mapping);
    if(component.parameter_stability)
        for(auto & block:component.parameter_stability->blocks)
            for(auto & atom:block.core_atoms) RemapAtomIndex(atom,atom_mapping);
    RemapTargetEvidence(component.target_evidence,atom_mapping,component.rows.size());
}
std::vector<const AtomObject *> EligibleAtoms(const ModelObject & model)
{
    std::vector<const AtomObject *> atoms;
    for(const auto & atom:model.GetAtomList()) if(atom->GetElement()!=Element::HYDROGEN) atoms.push_back(atom.get());
    std::sort(atoms.begin(),atoms.end(),[](const auto * left,const auto * right) {
        return left->GetSerialID()<right->GetSerialID();
    });
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
struct TargetSpatialBounds
{
    std::array<double,3> minimum{},maximum{};
    bool available{};
};
TargetSpatialBounds BoundsForRows(const MapObject & map,const std::vector<std::size_t> & rows)
{
    TargetSpatialBounds bounds;
    if(rows.empty()) return bounds;
    bounds.minimum.fill(std::numeric_limits<double>::infinity());
    bounds.maximum.fill(-std::numeric_limits<double>::infinity());
    const auto spacing=map.GetGridSpacing(),origin=map.GetOrigin();
    for(const auto row:rows)
    {
        const auto point=simulation::GridPosition(map.GetGridIndex(row),spacing,origin);
        for(std::size_t axis=0;axis<3;++axis)
        {
            if(!std::isfinite(point[axis])) return {};
            bounds.minimum[axis]=std::min(bounds.minimum[axis],point[axis]);
            bounds.maximum[axis]=std::max(bounds.maximum[axis],point[axis]);
        }
    }
    bounds.available=true;
    return bounds;
}
bool MayOverlapTargetRows(const AtomObject & atom,const TargetSpatialBounds & bounds)
{
    if(!bounds.available) return true;
    const auto position=atom.GetPosition();
    for(const double coordinate:position) if(!std::isfinite(coordinate)) return true;
    std::array<double,3> nearest{};
    for(std::size_t axis=0;axis<3;++axis)
        nearest[axis]=std::clamp(position[axis],bounds.minimum[axis],bounds.maximum[axis]);
    return simulation::SupportSquare(nearest,position)<=6.25;
}
bool AddEstimatedBytes(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return n::CheckedMultiplyBytes(count,width,bytes) && n::CheckedAddBytes(total,bytes,total);
}
[[noreturn]] void RejectBuildPreflight(const n::MemoryShape & shape,const char * reason)
{
    auto estimate=n::PlanMemory(shape,n::MemoryStage::ProblemConstruction,JointProblemMemoryBudget());
    estimate.status=estimate.status==n::MemoryStatus::Overflow || std::string_view(reason)=="size-overflow" ?
        n::MemoryStatus::Overflow : n::MemoryStatus::Reject;
    estimate.reason=reason; n::LogMemoryEstimate(estimate);
    throw std::runtime_error(std::string("resource-limit:build-joint-problem:")+reason);
}
JointState State(const n::Endpoint & e,double scale)
{return {Values(e.beta),Values(e.eta.array().exp()),Values(e.eta),Values(e.gradient),e.certificate.objective/(scale*scale)};}
void ConfigureOptimizationDomain(n::ProblemData & data)
{
    const auto & input=*data.input;
    data.context=n::CreateContext(data.input);
    data.layout=n::BuildParameterLayout(input);
    data.partition=n::Partition(data.domain,data.context.atom_ids);
}

}
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
void SetJointProblemMemoryBudgetForTesting(std::optional<std::size_t> budget)
{
    problem_memory_budget_for_testing=budget;
}
#endif
namespace joint_component {
JointParameterLayout BuildParameterLayout(const JointProblemInput & input)
{
    JointParameterLayout out;
    out.full_atoms.resize(input.atom_ids.size());
    std::iota(out.full_atoms.begin(),out.full_atoms.end(),0);
    out.informative_rows.resize(input.observations.size());
    std::iota(out.informative_rows.begin(),out.informative_rows.end(),0);
    return out;
}

std::vector<JointCheck> AssessmentEvidence(const Assessment & a,JointEvidenceScope scope)
{
    using Status=JointCheckStatus;
    auto check=[scope](std::string name,bool available,bool pass,std::optional<double> value={},std::optional<double> threshold={}) {
        return JointCheck{std::move(name),available ? (pass ? Status::Passed : Status::Failed) : Status::Unavailable,scope,value,threshold,{}};
    };
    std::vector<JointCheck> out;
    const auto & replay=a.independent_replay;
    const bool resource_unavailable=a.failure.starts_with("resource-limit:");
    const bool evidence_available=replay.available || (!a.primary.valid && !resource_unavailable);
    out.push_back(check("ac-feasibility",evidence_available,replay.feasible,
        evidence_available ? std::optional<double>(replay.feasible ? 0.0 : 1.0) : std::nullopt,0.0));
    out.push_back(check("conditional-ac-optimality",evidence_available,
        replay.primary_kkt_passed && replay.feasible,
        std::isfinite(replay.primary_ac_kkt) ? std::optional<double>(replay.primary_ac_kkt) : std::nullopt,
        a.convergence_criteria.global_ac_kkt));
    out.push_back(check("independent-ac-optimality",evidence_available,
        replay.independent_kkt_passed && replay.feasible,
        std::isfinite(replay.independent_ac_kkt) ? std::optional<double>(replay.independent_ac_kkt) : std::nullopt,
        a.convergence_criteria.global_ac_kkt));
    out.push_back(check("independent-kkt-consistency",evidence_available,replay.kkt_agrees,
        std::isfinite(replay.kkt_difference) ? std::optional<double>(replay.kkt_difference) : std::nullopt,1e-13));
    out.push_back(check("independent-replay",evidence_available,replay.valid && replay.passed));
    out.back().reason=replay.reason.empty() ? a.failure : replay.reason;
    out.push_back(check("width-gradient-consistency",evidence_available,replay.gradient_agrees,
        std::isfinite(replay.gradient_difference) ? std::optional<double>(replay.gradient_difference) : std::nullopt));
    out.back().reason=replay.reason.empty() ? "per-coordinate-replay-enclosure" : replay.reason;
    const double width_stationarity_value=replay.width_bound_applied ?
        std::max(replay.primary_width_projected_kkt,replay.independent_width_projected_kkt) :
        replay.primary_width_gradient;
    const bool width_stationarity_passed=replay.width_bound_applied ?
        replay.width_bound_feasible && replay.primary_width_kkt_passed && replay.independent_width_kkt_passed &&
            std::isfinite(replay.primary_width_projected_kkt) &&
            std::isfinite(replay.independent_width_projected_kkt) &&
            replay.primary_width_projected_kkt<=a.convergence_criteria.global_width_gradient &&
            replay.independent_width_projected_kkt<=a.convergence_criteria.global_width_gradient :
        replay.primary_width_gradient<=a.convergence_criteria.global_width_gradient;
    out.push_back(check("width-stationarity",evidence_available && a.primary.valid && a.primary.gradient.size()>0,
        width_stationarity_passed,
        std::isfinite(width_stationarity_value) ? std::optional<double>(width_stationarity_value) : std::nullopt,
        a.convergence_criteria.global_width_gradient));
    if(replay.width_bound_applied) out.back().reason="independent-bound-aware-projected-width-kkt";
    if(replay.width_bound_applied)
    {
        out.push_back(check("width-bound-feasibility",evidence_available,
            replay.width_bound_feasible,
            replay.width_bound_feasible ? 0.0 : 1.0,0.0));
        out.back().reason="primary-and-independent-eta-upper-bound";
    }
    out.push_back(check("endpoint-verification",evidence_available,a.endpoint_verified));
    out.push_back(check("numerical-qualification",a.primary.valid,a.numerical_qualified));
    if(a.verification_memory_estimate.stage==n::MemoryStage::EndpointVerification)
        out.push_back({"endpoint-verification-memory",
            a.verification_memory_estimate.status==n::MemoryStatus::Pass ? Status::Passed : Status::Unavailable,
            scope,static_cast<double>(a.verification_memory_estimate.estimated_peak_bytes),
            static_cast<double>(a.verification_memory_estimate.budget_bytes),
            a.verification_memory_estimate.reason.empty() ? "" :
                std::string("resource-limit:endpoint-verification:")+a.verification_memory_estimate.reason});
    for(const char * name:{"two-step-derivative","richardson","precision-50-100","boundary-audit","regular-certificate"})
        out.push_back({name,Status::NotRun,scope,{},{},"offline-audit-not-run"});
    for(auto & evidence:out) if(evidence.status==Status::Unavailable && evidence.reason.empty())
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
    const bool reference_free=std::any_of(evidence.begin(),evidence.end(),[&](const auto & check) {
        return check.name=="endpoint-verification" && check.scope==scope;
    });
    const std::vector<const char *> required=reference_free ?
        std::vector<const char *>{"ac-feasibility","conditional-ac-optimality","independent-ac-optimality",
            "independent-kkt-consistency","independent-replay","width-gradient-consistency",
            "width-stationarity","endpoint-verification","numerical-qualification"} :
        std::vector<const char *>{"inner","kkt","width-stationarity","numerical-qualification"};
    for(const char * name:required)
    {
        const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & check){return check.name==name && check.scope==scope;});
        status=MergeConvergenceStatus(status,found==evidence.end() ? JointCheckStatus::Unavailable : found->status);
    }
    if(assembled)
    {
        const std::vector<const char *> assembly_required=reference_free ?
            std::vector<const char *>{"assembly-mapping","assembled-consistency"} :
            std::vector<const char *>{"assembled-profile"};

        for(const char * name:assembly_required)
        {
            const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & check){return check.name==name && check.scope==scope;});
            status=MergeConvergenceStatus(status,found==evidence.end() ? JointCheckStatus::Unavailable : found->status);
        }
    }
    if(std::any_of(evidence.begin(),evidence.end(),[&](const auto & check) {
        return check.name=="width-bound-feasibility" && check.scope==scope;
    }))
    {
        const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & check) {
            return check.name=="width-bound-feasibility" && check.scope==scope;
        });
        status=MergeConvergenceStatus(status,found->status);
    }
    return status;
}
}
JointCheckStatus JointComponentResult::RuntimeConvergence() const
{
    if(!state || !ValidJointState(*state,layout ? layout->full_atoms.size() : atoms.size()))
        return JointCheckStatus::Unavailable;
    auto status=joint_component::ConvergenceStatus(evidence,JointEvidenceScope::ComponentLocal);
    const auto stationary=std::find_if(evidence.begin(),evidence.end(),[](const auto & check) {
        return check.name=="complete-sweep-stationarity" && check.scope==JointEvidenceScope::ComponentLocal;
    });
    status=joint_component::MergeConvergenceStatus(status,stationary==evidence.end() ?
        JointCheckStatus::Unavailable : stationary->status);
    const auto active_set=std::find_if(evidence.begin(),evidence.end(),[](const auto & check) {
        return check.name=="active-set-lm-subproblems" && check.scope==JointEvidenceScope::ComponentLocal;
    });
    if(active_set!=evidence.end()) status=joint_component::MergeConvergenceStatus(status,active_set->status);
    return status;
}
JointCheckStatus JointFitResult::RuntimeConvergence() const
{
    if(!assembled_state || !prediction || !objective || components.empty() ||
        !ValidJointState(*assembled_state,assembled_state->b.size()) || !std::isfinite(*objective) || *objective<0 ||
        std::abs(assembled_state->objective-*objective)>8*std::numeric_limits<double>::epsilon()*
            std::max(1.0,std::abs(*objective)) ||
        prediction->size()!=available_row_mask.size() ||
        !std::all_of(prediction->begin(),prediction->end(),[](double value){return std::isfinite(value);}) ||
        !std::all_of(available_row_mask.begin(),available_row_mask.end(),[](bool available){return available;}))
        return JointCheckStatus::Unavailable;
    auto status=joint_component::ConvergenceStatus(evidence,JointEvidenceScope::AssembledGlobal,true);
    for(const auto & component:components) status=joint_component::MergeConvergenceStatus(status,component.RuntimeConvergence());
    return status;
}
JointProblem::JointProblem(JointProblemInput input)
{
    if(input.atom_ids.empty() || input.support.size()!=input.atom_ids.size() || input.row_ids.size()!=input.observations.size())
        throw std::invalid_argument("Invalid joint problem identities or dimensions.");
    if(input.selection_domain && !input.selection_domain->IsValid(input.atom_ids.size()))
        throw std::invalid_argument("Invalid joint selection domain.");
    if(!std::all_of(input.observations.begin(),input.observations.end(),[](double value){return std::isfinite(value);}))
        throw std::invalid_argument("Nonfinite joint observations.");
    for(const auto & atom:input.support) for(const auto & s:atom)
        if(s.row>=input.observations.size() || !std::isfinite(s.squared_distance) || s.squared_distance<0 || s.squared_distance>6.25)
            throw std::invalid_argument("Invalid joint contributor row.");
    n::MemoryShape shape;
    shape.component_atoms=input.atom_ids.size();
    shape.component_observations=input.observations.size();
    shape.problem_scope=true;
    const auto support_entries=n::CountSupportEntries(input);
    const auto resident=n::EstimateSnapshotBytes(input);
    const auto construction_scratch=n::EstimateConstructionScratchBytes(input);
    std::size_t free_columns{};
    const bool dimensions_ok=n::CheckedMultiplyBytes(input.atom_ids.size(),2,free_columns);
    n::MemoryEstimate construction_memory;
    if(!support_entries || !resident || !construction_scratch || !dimensions_ok)
    {
        construction_memory.shape=shape;
        construction_memory.stage=n::MemoryStage::ProblemConstruction;
        construction_memory.status=n::MemoryStatus::Overflow;
        construction_memory.reason="size-overflow";
        construction_memory.budget_bytes=n::JointMemoryPlanningBudgetBytes;
    }
    else
    {
        shape.component_support_entries=*support_entries; shape.resident_bytes=*resident;
        shape.identity_validation_bytes=*construction_scratch;
        shape.component_free_columns=free_columns;
        construction_memory=n::PlanMemory(shape,n::MemoryStage::ProblemConstruction);
    }
    if(construction_memory.status!=n::MemoryStatus::Pass ||
        construction_memory.estimated_peak_bytes>n::JointMemoryPlanningBudgetBytes/2)
        n::LogMemoryEstimate(construction_memory);
    if(construction_memory.status==n::MemoryStatus::Pass &&
        std::set<std::string>(input.row_ids.begin(),input.row_ids.end()).size()!=input.row_ids.size())
        throw std::invalid_argument("Invalid joint problem identities or dimensions.");
    auto data=std::make_shared<n::ProblemData>(std::make_shared<const JointProblemInput>(std::move(input)));
    data->construction_memory=construction_memory;
    if(data->construction_memory.status==n::MemoryStatus::Pass)
        ConfigureOptimizationDomain(*data);
    m_data=std::move(data);
}
const JointProblemInput & JointProblem::Input() const {return *m_data->input;}
const JointProblemInput & JointProblem::OriginalInput() const {return *m_data->original_input;}
const std::vector<std::size_t> & JointProblem::FittingToOriginalContributors() const
{return m_data->fitting_to_original_contributors;}
const std::vector<std::size_t> & JointProblem::FittingToOriginalRows() const
{return m_data->fitting_to_original_rows;}
JointProblem::JointProblem(std::shared_ptr<const n::ProblemData> data):m_data(std::move(data)) {}
JointProblem JointProblem::WithReconstructionDomain(
    JointProblemInput original, std::vector<std::size_t> fitting_to_original_contributors,
    std::vector<std::size_t> fitting_to_original_rows) const
{
    const auto & fitting=*m_data->input;
    if (original.atom_ids.empty() || original.support.size()!=original.atom_ids.size() ||
        original.row_ids.size()!=original.observations.size() ||
        fitting_to_original_contributors.size()!=fitting.atom_ids.size() ||
        fitting_to_original_rows.size()!=fitting.row_ids.size())
        throw std::invalid_argument("Invalid Joint reconstruction domain dimensions.");
    for (std::size_t atom=0;atom<fitting.atom_ids.size();++atom)
    {
        const auto original_atom=fitting_to_original_contributors[atom];
        if (original_atom>=original.atom_ids.size() || fitting.atom_ids[atom]!=original.atom_ids[original_atom])
            throw std::invalid_argument("Invalid fitting-to-original contributor mapping.");
    }
    for (std::size_t row=0;row<fitting.row_ids.size();++row)
    {
        const auto original_row=fitting_to_original_rows[row];
        if (original_row>=original.row_ids.size() || fitting.row_ids[row]!=original.row_ids[original_row])
            throw std::invalid_argument("Invalid fitting-to-original row mapping.");
    }
    if (!std::all_of(original.observations.begin(),original.observations.end(),
            [](double value){return std::isfinite(value);}))
        throw std::invalid_argument("Nonfinite Joint reconstruction observations.");
    auto data=std::make_shared<n::ProblemData>(*m_data);
    data->original_input=std::make_shared<const JointProblemInput>(std::move(original));
    data->fitting_to_original_contributors=std::move(fitting_to_original_contributors);
    data->fitting_to_original_rows=std::move(fitting_to_original_rows);
    const Eigen::Map<const Eigen::VectorXd> original_values(
        data->original_input->observations.data(),static_cast<Eigen::Index>(data->original_input->observations.size()));
    const double original_norm=original_values.norm();
    data->context.scale=std::max(1.0,original_norm);
    data->context.linear.release_response_norm=original_norm;
    return JointProblem(std::move(data));
}
double JointProblem::ObservationScale() const {return m_data->context.scale;}
const JointParameterLayout & JointProblem::ParameterLayout() const {return m_data->layout;}
namespace {
struct SelectedObservationClosure
{
    JointProblemInput original;
    std::vector<std::size_t> target_indices;
    std::size_t candidate_halos{},retained_halos{},zero_support_halos{};
};

SelectedObservationClosure BuildSelectedObservationClosure(const MapObject & map,const ModelObject & model)
{
    const auto atoms=EligibleAtoms(model); SelectedObservationClosure closure;
    auto & input=closure.original;
    std::set<const AtomObject *> selected;
    for(const auto * atom:model.GetSelectedAtoms()) if(atom->GetElement()!=Element::HYDROGEN) selected.insert(atom);
    if(selected.empty()) throw std::invalid_argument("Joint fitting requires at least one selected non-hydrogen atom.");
    const auto dims=map.GetGridSize(); const auto spacing=map.GetGridSpacing(),origin=map.GetOrigin();
    for(std::size_t k=0;k<3;++k) if(dims[k]<=0 || !std::isfinite(spacing[k]) || spacing[k]<=0 || !std::isfinite(origin[k]))
        throw std::invalid_argument("Invalid joint map geometry.");
    std::size_t selected_support_entries{};
    for(const auto * atom:atoms) if(selected.contains(atom))
        VisitSupport(map,*atom,[&](auto,double) {
            if(!n::CheckedAddBytes(selected_support_entries,1,selected_support_entries))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
        });
    std::size_t selected_row_buffer{};
    if(!n::CheckedMultiplyBytes(selected_support_entries,sizeof(std::size_t),selected_row_buffer))
        throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
    std::size_t builder_resident=selected_row_buffer,builder_tracking{};
    if(!AddEstimatedBytes(builder_tracking,atoms.capacity(),sizeof(const AtomObject *)) ||
        !AddEstimatedBytes(builder_tracking,selected.size(),4*sizeof(void *)+32) ||
        !AddEstimatedBytes(builder_tracking,selected_support_entries,sizeof(std::size_t)) ||
        !AddEstimatedBytes(builder_tracking,atoms.size(),sizeof(std::size_t)) ||
        !n::CheckedAddBytes(builder_resident,builder_tracking,builder_resident))
        throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
    if(builder_resident>JointProblemMemoryBudget())
    {
        n::MemoryShape shape; shape.component_atoms=atoms.size(); shape.component_observations=selected_support_entries;
        shape.component_support_entries=selected_support_entries; shape.resident_bytes=builder_resident;
        shape.problem_scope=true;
        RejectBuildPreflight(shape,"raw-target-row-buffer-exceeds-budget");
    }
    std::vector<std::size_t> rows; rows.reserve(selected_support_entries);
    for(const auto * atom:atoms) if(selected.contains(atom))
        VisitSupport(map,*atom,[&](auto index,double) {rows.push_back(index);});
    std::sort(rows.begin(),rows.end()); rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    const auto target_bounds=BoundsForRows(map,rows);

    std::size_t contributors=selected.size(),support_entries{},target_support_entries{},maximum_atom_support{};
    std::size_t candidate_halos{},retained_halos{},zero_support_halos{};
    std::vector<std::size_t> atom_support_counts;
    atom_support_counts.reserve(atoms.size());
    for(const auto * atom:atoms)
    {
        const bool selected_target=selected.contains(atom);
        if(!selected_target && !n::CheckedAddBytes(candidate_halos,1,candidate_halos))
            throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
        std::size_t memberships{};
        if(!selected_target && !MayOverlapTargetRows(*atom,target_bounds))
        {
            if(!n::CheckedAddBytes(zero_support_halos,1,zero_support_halos))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
        }
        else
        {
            VisitSupport(map,*atom,[&](auto index,double) {
                if(std::binary_search(rows.begin(),rows.end(),index) &&
                    !n::CheckedAddBytes(memberships,1,memberships))
                    throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
            });
            if(!selected_target && memberships==0 &&
                !n::CheckedAddBytes(zero_support_halos,1,zero_support_halos))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
        }
        atom_support_counts.push_back(memberships);
        maximum_atom_support=std::max(maximum_atom_support,memberships);
        if(selected_target)
        {
            if(!n::CheckedAddBytes(support_entries,memberships,support_entries) ||
                !n::CheckedAddBytes(target_support_entries,memberships,target_support_entries))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
            continue;
        }
        if(memberships>0 && (!n::CheckedAddBytes(retained_halos,1,retained_halos) ||
            !n::CheckedAddBytes(contributors,1,contributors) ||
            !n::CheckedAddBytes(support_entries,memberships,support_entries)))
            throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
    }
    n::MemoryShape shape; shape.problem_scope=true; shape.component_atoms=selected.size();
    shape.component_observations=rows.size(); shape.component_support_entries=support_entries;
    std::size_t free_columns{};
    if(!n::CheckedMultiplyBytes(selected.size(),2,free_columns))
        RejectBuildPreflight(shape,"size-overflow");
    shape.component_free_columns=free_columns;
    std::size_t planned_snapshot{};
    std::size_t identity_characters{},identity_scratch{},identity_map_entries{};
    const auto identity_node=sizeof(std::string)+4*sizeof(void *)+32;
    if(!AddEstimatedBytes(identity_characters,contributors,32) || !AddEstimatedBytes(identity_characters,rows.size(),32) ||
        !AddEstimatedBytes(identity_scratch,rows.size(),identity_node+32) ||
        !AddEstimatedBytes(identity_scratch,contributors,2*sizeof(std::string)+identity_node+64) ||
        !AddEstimatedBytes(identity_scratch,maximum_atom_support,sizeof(std::size_t)+4*sizeof(void *)+32) ||
        !n::CheckedAddBytes(contributors,selected.size(),identity_map_entries) ||
        !AddEstimatedBytes(planned_snapshot,identity_characters,2*sizeof(char)))
        RejectBuildPreflight(shape,"size-overflow");
    std::size_t all_support_entries{},all_identity_entries{};
    if(!n::CheckedAddBytes(support_entries,target_support_entries,all_support_entries) ||
        !n::CheckedAddBytes(contributors,selected.size(),all_identity_entries))
        RejectBuildPreflight(shape,"size-overflow");
    shape.identity_validation_bytes=identity_scratch;
    if(!AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(double)) ||
        !AddEstimatedBytes(planned_snapshot,all_identity_entries,sizeof(std::vector<JointSupport>)) ||
        !AddEstimatedBytes(planned_snapshot,all_identity_entries,sizeof(std::string)) ||
        !AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(std::string)) ||
        !AddEstimatedBytes(planned_snapshot,all_support_entries,sizeof(JointSupport)) ||
        !AddEstimatedBytes(planned_snapshot,selected.size(),sizeof(std::size_t)) ||
        !AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(std::size_t)) ||
        !AddEstimatedBytes(planned_snapshot,rows.size(),sizeof(std::uint8_t)) ||
        !AddEstimatedBytes(planned_snapshot,selected.size(),sizeof(std::vector<n::LocalSupportIndex>)) ||
        !AddEstimatedBytes(planned_snapshot,selected_support_entries,sizeof(n::LocalSupportIndex)) ||
        !AddEstimatedBytes(planned_snapshot,identity_map_entries,sizeof(std::size_t)) ||
        !n::CheckedAddBytes(planned_snapshot,builder_resident,planned_snapshot))
        RejectBuildPreflight(shape,"size-overflow");
    shape.resident_bytes=planned_snapshot;
    const auto build_estimate=n::PlanMemory(shape,n::MemoryStage::ProblemConstruction,
        JointProblemMemoryBudget());
    if(build_estimate.status!=n::MemoryStatus::Pass)
        RejectBuildPreflight(shape,build_estimate.reason.c_str());

    for(auto index:rows) {input.row_ids.push_back(std::to_string(index)); input.observations.push_back(map.GetMapValue(index));}
    for(std::size_t atom_index=0;atom_index<atoms.size();++atom_index)
    {
        const auto * atom=atoms[atom_index];
        bool selected_target=selected.contains(atom);
        const auto memberships=atom_support_counts[atom_index];
        if(!selected_target && memberships==0) continue;
        std::vector<JointSupport> support;
        support.reserve(memberships);
        VisitSupport(map,*atom,[&](auto index,double square) {
            const auto row=std::lower_bound(rows.begin(),rows.end(),index);
            if(row!=rows.end() && *row==index)
                support.push_back({static_cast<std::size_t>(row-rows.begin()),square});
        });
        if(selected_target || !support.empty())
        {
            if(selected_target) closure.target_indices.push_back(input.atom_ids.size());
            input.atom_ids.push_back(std::to_string(atom->GetSerialID()));
            input.support.push_back(std::move(support));
        }
    }
    closure.candidate_halos=candidate_halos;
    closure.retained_halos=retained_halos;
    closure.zero_support_halos=zero_support_halos;
    return closure;
}
}

JointProblem BuildPureHaloSupportExclusionProblem(const MapObject & map,const ModelObject & model)
{
    auto closure=BuildSelectedObservationClosure(map,model);
    const auto & source=closure.original;
    const auto & targets=closure.target_indices;
    std::vector<bool> retained(source.observations.size(),true);
    for (std::size_t atom=0;atom<source.support.size();++atom)
    {
        if (std::binary_search(targets.begin(),targets.end(),atom)) continue;
        for (const auto & support:source.support[atom])
        {
            if (support.row>=retained.size()) throw std::invalid_argument("Invalid Halo support incidence.");
            retained[support.row]=false;
        }
    }
    const auto retained_count=static_cast<std::size_t>(std::count(retained.begin(),retained.end(),true));
    if (!retained_count) throw std::invalid_argument("Pure Halo exclusion leaves an empty fitting domain.");
    JointProblemInput fitting;
    fitting.selection_domain.emplace();
    fitting.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    fitting.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        source.observations.size(),retained_count,source.observations.size()-retained_count,0,0,0};
    std::vector<std::size_t> fitting_rows;
    fitting_rows.reserve(retained_count);
    std::vector<std::size_t> original_to_fitting(source.observations.size(),std::numeric_limits<std::size_t>::max());
    for (std::size_t row=0;row<retained.size();++row) if (retained[row])
    {
        original_to_fitting[row]=fitting_rows.size();
        fitting_rows.push_back(row);
        fitting.row_ids.push_back(source.row_ids[row]);
        fitting.observations.push_back(source.observations[row]);
    }
    fitting.atom_ids.reserve(targets.size());
    fitting.support.reserve(targets.size());
    fitting.selection_domain->target_indices.resize(targets.size());
    std::iota(fitting.selection_domain->target_indices.begin(),fitting.selection_domain->target_indices.end(),0);
    for (const auto original_atom:targets)
    {
        fitting.atom_ids.push_back(source.atom_ids[original_atom]);
        std::vector<JointSupport> support;
        support.reserve(source.support[original_atom].size());
        for (const auto entry:source.support[original_atom]) if (retained[entry.row])
            support.push_back({original_to_fitting[entry.row],entry.squared_distance});
        if (support.empty()) throw std::invalid_argument("Pure Halo exclusion removes all observations for a Target.");
        fitting.support.push_back(std::move(support));
    }
    auto original=std::move(closure.original);
    original.selection_domain.emplace();
    original.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    original.selection_domain->target_indices=targets;
    original.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        original.observations.size(),retained_count,original.observations.size()-retained_count,
        closure.candidate_halos,closure.retained_halos,closure.zero_support_halos};
    std::vector<std::size_t> fitting_to_original_contributors=targets;
    return JointProblem(std::move(fitting)).WithReconstructionDomain(
        std::move(original),std::move(fitting_to_original_contributors),std::move(fitting_rows));
}
JointFitResult FitJointComponentsInitialized(const JointProblem & problem,const JointInitialization & initialization)
{
    return detail::FitJointComponentsInitializedWithProgress(problem,initialization,{});
}
JointFitResult detail::FitJointComponentsInitializedWithProgress(
    const JointProblem & problem,const JointInitialization & initialization,
    const joint_component::JointProgressObserver & observer)
{
    const auto & input=problem.Input();
    const auto invalid=[&](std::string reason) {
        JointFitResult out; out.problem=problem; out.initialization=initialization;
        out.initialization.valid=false; out.initialization.reason=std::move(reason);
        out.observation_scale=problem.ObservationScale();
        out.available_row_mask.assign(input.observations.size(),false);
        out.evidence.push_back({"joint-initialization",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},out.initialization.reason});
        return out;
    };
    if(!input.selection_domain ||
        input.selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy)
        return invalid("joint-initialization-unavailable:missing-pure-exclusion-policy");
    if(initialization.atoms.size()!=input.atom_ids.size() || initialization.b.size()!=input.atom_ids.size())
        return invalid("joint-initialization-unavailable:target-mapping-mismatch");
    if(!initialization.valid) return invalid(initialization.reason);
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
        if(initialization.atoms[atom].id!=input.atom_ids[atom] || !initialization.atoms[atom].target)
            return invalid("joint-initialization-unavailable:target-mapping-mismatch");

    const auto calibration=joint_component::CalibrateGlobalWidth(initialization,"joint-target-widths");
    const auto bounded=joint_component::MakeBoundedInitialization(initialization,calibration);
    if(!bounded.available)
    {
        JointFitResult out; out.problem=problem; out.initialization=initialization;
        out.observation_scale=problem.ObservationScale();
        out.available_row_mask.assign(input.observations.size(),false);
        out.evidence.push_back({"global-width-calibration",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,static_cast<double>(calibration.qualified_donor_count),
            static_cast<double>(calibration.target_count),"calibration-unavailable:"+bounded.reason});
        return out;
    }
    joint_component::FixedNeighborSearchPolicy solver_policy;
    solver_policy.eta_upper_bound=calibration.eta_max4;
    auto out=n::FitFixedNeighborComponents(problem,bounded.widths,solver_policy,observer);
    out.initialization=initialization;
    out.evidence.push_back({"global-width-calibration",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,calibration.p95,bounded.bmax,
        "global-first-stage-p95-k4; qualified-donors="+
            std::to_string(calibration.qualified_donor_count)+"/"+std::to_string(calibration.target_count)});
    out.evidence.push_back({"qualified-first-stage-donor-coverage",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,calibration.qualified_donor_coverage,0.95,
        "target-only-qualified-first-stage-widths"});
    out.evidence.push_back({"bounded-initialization-projection",JointCheckStatus::Passed,
        JointEvidenceScope::AssembledGlobal,static_cast<double>(bounded.projection_count),{},
        "projected-target-widths="+std::to_string(bounded.projection_count)+
            "; maximum-projection-angstroms="+std::to_string(bounded.maximum_projection_magnitude)});
    return out;
}
JointFitResult n::FitFixedNeighborComponentsImpl(const JointProblem & problem,const std::vector<double> & initial_b,
    const n::FixedNeighborSearchPolicy & production_policy,const n::JointProgressObserver & observer)
{
    eigen_helper::ScopedEigenThreadCount eigen_thread_guard{1};
    constexpr auto fixed_neighbor_core_atoms=n::ProductionFixedNeighborCoreAtoms;
    constexpr auto fixed_neighbor_maximum_sweeps=n::ProductionFixedNeighborMaximumSweeps;
    constexpr auto convergence_criteria=n::PracticalJointConvergenceCriteria;
    if(observer)
    {
        n::JointProgressEvent event;
        event.phase=n::JointProgressPhase::SolverConfigured;
        event.solver_configuration=n::JointSolverConfiguration{
            fixed_neighbor_core_atoms};
        observer(event);
    }
    auto provenance=[&] {
        JointSolverProvenance out;
        out.search_method="fixed-neighbor";
        out.contract_version=std::string(JointSolverProvenanceContractVersion);
        out.sparse_backend="EIGEN";
        out.fixed_neighbor_core_atoms=fixed_neighbor_core_atoms;
        out.fixed_neighbor_policy_version=std::string(FixedNeighborPolicyContractVersion);
        out.fixed_neighbor_maximum_sweeps=fixed_neighbor_maximum_sweeps;
        out.fixed_neighbor_order="forward";
        out.fixed_neighbor_local_search=std::string(FixedNeighborBoundedLocalSearchName);
        return out;
    };
    const auto & data=JointProblemAccess::Get(problem);
    if(data.construction_memory.status!=n::MemoryStatus::Pass)
    {
        JointFitResult out; out.problem=problem; out.solver_provenance=provenance();
        out.observation_scale=data.context.scale; out.initialization.b=initial_b;
        out.initialization.valid=initial_b.size()==data.domain.atoms.size() &&
            std::all_of(initial_b.begin(),initial_b.end(),[](double b){return std::isfinite(b) && b>0;});
        out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
        out.available_row_mask.assign(data.input->observations.size(),false);
        out.evidence.push_back({"memory-preflight",JointCheckStatus::Unavailable,JointEvidenceScope::AssembledGlobal,
            data.construction_memory.estimated_peak_bytes ?
                std::optional<double>(static_cast<double>(data.construction_memory.estimated_peak_bytes)) : std::nullopt,
            static_cast<double>(data.construction_memory.budget_bytes),
            "resource-limit:problem-construction:"+data.construction_memory.reason});
        return out;
    }
    JointFitResult out; out.problem=problem;
    out.solver_provenance=provenance();
    out.observation_scale=data.context.scale; out.initialization.b=initial_b;
    out.initialization.valid=initial_b.size()==data.domain.atoms.size() &&
        std::all_of(initial_b.begin(),initial_b.end(),[](double b){return std::isfinite(b) && b>0;});
    out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
    if(initial_b.size()!=data.domain.atoms.size())
    {
        out.available_row_mask.assign(data.input->observations.size(),false);
        for(auto row:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(row)]=true;
        return out;
    }
    if(!production_policy.eta_upper_bound || !std::isfinite(*production_policy.eta_upper_bound))
    {
        out.available_row_mask.assign(data.input->observations.size(),false);
        for(auto row:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(row)]=true;
        out.evidence.push_back({"global-width-bound",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},"finite-width-upper-bound-required"});
        return out;
    }
    const auto component_count=data.partition.components.size();
    auto bounded_policy=production_policy;
    auto assembly_context=data.context;
    assembly_context.eta_upper_bound=bounded_policy.eta_upper_bound;
    std::size_t remaining_perturbation_blocks=n::MaximumPerturbedBlocksPerFit;
    std::size_t remaining_general_perturbation_blocks=
        n::MaximumPerturbedBlocksPerFit>1 ? n::MaximumPerturbedBlocksPerFit-1 : 0;
    bool reserved_bound_active_probe=n::MaximumPerturbedBlocksPerFit>0;
    const auto request_perturbation_block=production_policy.claim_parameter_stability_block;
    bounded_policy.claim_parameter_stability_block=[&remaining_perturbation_blocks,
        &remaining_general_perturbation_blocks,&reserved_bound_active_probe,
        request_perturbation_block](int priority) {
        if(!remaining_perturbation_blocks) return false;
        if(priority==0 && reserved_bound_active_probe)
        {
            if(request_perturbation_block && !request_perturbation_block(priority)) return false;
            reserved_bound_active_probe=false;
            --remaining_perturbation_blocks;
            return true;
        }
        if(!remaining_general_perturbation_blocks ||
            (request_perturbation_block && !request_perturbation_block(priority))) return false;
        --remaining_general_perturbation_blocks;
        --remaining_perturbation_blocks;
        return true;
    };
    auto component_batch=n::OrchestrateFixedNeighborComponents(problem,initial_b,bounded_policy,observer);
    out.components=std::move(component_batch.components);
    out.costs.search_seconds+=component_batch.costs.search_seconds;
    out.costs.assessment_seconds+=component_batch.costs.assessment_seconds;
    const auto assembly_start=Clock::now();
    n::NotifyJointAssemblyProgress(observer,n::JointProgressPhase::AssemblyStarted,component_count);
    const auto assembly=n::AssembleComponents(data.domain,data.y,data.partition,assembly_context,out.components,
        convergence_criteria);
    out.search_completed=assembly.completed; out.available_row_mask=assembly.row_mask;
    out.evidence=n::AssessmentEvidence(assembly.assessment,JointEvidenceScope::AssembledGlobal);
    out.evidence.push_back({"assembly-mapping",assembly.mapping_valid ? JointCheckStatus::Passed :
        JointCheckStatus::Unavailable,JointEvidenceScope::AssembledGlobal,{},{},
        assembly.mapping_valid ? "exact-component-identity-mapping" :
            assembly.mapping_reason.empty() ? "missing-component-state" : assembly.mapping_reason});
    if(assembly.available)
    {
        const auto & state=assembly.raw;
        const auto consistency_status=assembly.consistency_passed ? JointCheckStatus::Passed :
            !state.valid && assembly.verification_memory_estimate.status!=n::MemoryStatus::Pass ?
                JointCheckStatus::Unavailable : JointCheckStatus::Failed;
        out.evidence.push_back({"assembled-consistency",consistency_status,
            JointEvidenceScope::AssembledGlobal,
            std::isfinite(assembly.replay.prediction_difference) ?
                std::optional<double>(assembly.replay.prediction_difference) : std::nullopt,2e-12,
            assembly.replay.reason});
        if(assembly.verification_memory_estimate.status!=n::MemoryStatus::Pass)
            out.evidence.push_back({"assembly-verification-memory",JointCheckStatus::Unavailable,
                JointEvidenceScope::AssembledGlobal,
                static_cast<double>(assembly.verification_memory_estimate.estimated_peak_bytes),
                static_cast<double>(assembly.verification_memory_estimate.budget_bytes),
                std::string("resource-limit:endpoint-verification:")+assembly.verification_memory_estimate.reason});
        if(assembly.assessment.numerical_qualified && assembly.consistency_passed)
        {
            out.assembled_state=State(state,data.context.scale);
            out.prediction=Values(assembly.prediction);
            out.objective=assembly.objective/(data.context.scale*data.context.scale);
        }
    }
    out.costs.assembly_seconds=Seconds(assembly_start);
    n::NotifyJointAssemblyProgress(observer,n::JointProgressPhase::AssemblyCompleted,component_count,
        out.costs.assembly_seconds,assembly.available && assembly.raw.valid);
    return out;
}
namespace {
template<class Fit,class Value>
decltype(auto) CaptureTransfer(Value & value)
{
    if constexpr(std::is_rvalue_reference_v<Fit &&>) return std::move(value);
    else return (value);
}
template<class Fit>
JointAnalysisResult CaptureJointAnalysisResultImpl(Fit && fit,JointAnalysisMetadata metadata)
{
    if (!fit.problem) throw std::invalid_argument("Cannot capture joint result without problem identities.");
    const auto & fitting_input=fit.problem->Input();
    const auto & original_input=fit.problem->OriginalInput();
    const bool pure_exclusion_policy=fitting_input.selection_domain &&
        fitting_input.selection_domain->contributor_policy==PureHaloSupportExclusionContributorPolicy;
    const auto captured_runtime=fit.RuntimeConvergence();
    const auto captured_target_runtime=fit.TargetRuntimeConvergence();
    if constexpr(std::is_rvalue_reference_v<Fit &&>) fit.prediction.reset();
    JointAnalysisResult out;
    metadata.software=JointSoftwareProvenance{RHBM_GEM_SIMULATION_VERSION,
        RHBM_GEM_SIMULATION_SOURCE_SHA256,RHBM_GEM_SIMULATION_CONFIG_SHA256,
        RHBM_GEM_SIMULATION_BUILD_SHA256};
    out.metadata=std::move(metadata);
    if(!fit.solver_provenance.search_method.empty())
        out.metadata.solver=CaptureTransfer<Fit>(fit.solver_provenance);
    out.atom_ids=pure_exclusion_policy ? original_input.atom_ids : fitting_input.atom_ids;
    out.row_ids=pure_exclusion_policy ? original_input.row_ids : fitting_input.row_ids;
    out.selection_domain=pure_exclusion_policy ? original_input.selection_domain : fitting_input.selection_domain;
    if(pure_exclusion_policy)
    {
        out.parameterization_contract=std::string(JointPureHaloSupportExclusionParameterizationContract);
        out.fitting_row_ids=fitting_input.row_ids;
        out.fitting_to_original_contributors=fit.problem->FittingToOriginalContributors();
        out.fitting_to_original_rows=fit.problem->FittingToOriginalRows();
    }
    out.initialization=CaptureTransfer<Fit>(fit.initialization);
    out.costs=CaptureTransfer<Fit>(fit.costs);
    for(auto & component:fit.components)
    {
        const auto runtime=component.RuntimeConvergence();
        const auto target_runtime=component.TargetRuntimeConvergence();
        JointAnalysisComponent saved;
        static_cast<JointComponentData &>(saved)=CaptureTransfer<Fit>(component);
        if(pure_exclusion_policy)
            RemapComponentForOriginalDomain(saved,fit.problem->FittingToOriginalContributors(),
                fit.problem->FittingToOriginalRows());
        saved.runtime_convergence=runtime;
        saved.target_runtime_convergence=target_runtime;
        out.components.push_back(std::move(saved));
    }
    out.assembled_state=CaptureTransfer<Fit>(fit.assembled_state); out.objective=fit.objective;
    if(pure_exclusion_policy)
    {
        out.available_row_mask.assign(original_input.row_ids.size(),false);
        if(fit.available_row_mask.size()!=fit.problem->FittingToOriginalRows().size())
            throw std::invalid_argument("Invalid Pure Exclusion fitting row mask.");
        for(std::size_t row=0;row<fit.available_row_mask.size();++row)
        {
            const auto original_row=fit.problem->FittingToOriginalRows()[row];
            if(original_row>=out.available_row_mask.size())
                throw std::invalid_argument("Invalid Pure Exclusion fitting row mapping.");
            out.available_row_mask[original_row]=fit.available_row_mask[row];
        }
    }
    else out.available_row_mask=CaptureTransfer<Fit>(fit.available_row_mask);
    out.evidence=CaptureTransfer<Fit>(fit.evidence);
    out.search_completed=fit.search_completed; out.observation_scale=fit.observation_scale;
    out.runtime_convergence=captured_runtime; out.regular_certificate=fit.regular_certificate;
    out.target_evidence=CaptureTransfer<Fit>(fit.target_evidence);
    if(pure_exclusion_policy)
        RemapTargetEvidence(out.target_evidence,fit.problem->FittingToOriginalContributors(),original_input.row_ids.size());
    out.target_runtime_convergence=captured_target_runtime;
    return out;
}
}
JointAnalysisResult CaptureJointAnalysisResult(const JointFitResult & fit,JointAnalysisMetadata metadata)
{return CaptureJointAnalysisResultImpl(fit,std::move(metadata));}
JointAnalysisResult CaptureJointAnalysisResult(JointFitResult && fit,JointAnalysisMetadata metadata)
{return CaptureJointAnalysisResultImpl(std::move(fit),std::move(metadata));}
JointInitialization joint_component::InitializeContributors(MapObject & map,ModelObject & model,const JointProblem & problem)
{
    ModelObject initializer(model);
    const auto workset = detail::MakeJointFittingWorkset(initializer, problem);
    FitOptions options; options.thread_size = 1; options.quiet_mode = true;
    const auto initialization = detail::PrepareJointWidthInitialization(map, initializer, workset, options);
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
    const auto problem=BuildPureHaloSupportExclusionProblem(map,model);
    const double construction_seconds=Seconds(construction_start);
    const auto initialization_start=Clock::now();
    auto initialization=n::InitializeContributors(map,model,problem);
    const double initialization_seconds=Seconds(initialization_start);
    auto out=FitJointComponentsInitialized(problem,initialization); out.costs.initialization_seconds=initialization_seconds;
    out.costs.construction_seconds=construction_seconds;
    out.initialization=std::move(initialization);
    return out;
}
}
