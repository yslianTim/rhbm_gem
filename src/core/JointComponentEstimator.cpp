#include "SimulationBuildInfo.hpp"
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/HaloEligibility.hpp"
#include "core/detail/joint_component/TargetEvidence.hpp"
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"
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

namespace rhbm_gem::core {
namespace {
namespace n=joint_component;
using Clock=std::chrono::steady_clock;
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
bool AddEstimatedBytes(std::size_t & total,std::size_t count,std::size_t width)
{
    std::size_t bytes{}; return n::CheckedMultiplyBytes(count,width,bytes) && n::CheckedAddBytes(total,bytes,total);
}
[[noreturn]] void RejectBuildPreflight(const n::MemoryShape & shape,const char * reason)
{
    auto estimate=n::PlanMemory(shape,n::MemoryStage::ProblemConstruction);
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
    if(!input.selection_domain || input.selection_domain->contributor_policy!=FixedHaloContributorPolicy)
    {
        data.target_to_contributor.resize(input.atom_ids.size());
        std::iota(data.target_to_contributor.begin(),data.target_to_contributor.end(),0);
        data.contributor_to_target.resize(input.atom_ids.size());
        std::iota(data.contributor_to_target.begin(),data.contributor_to_target.end(),0);
        data.layout=n::BuildParameterLayout(input);
        data.partition=n::Partition(data.domain,input.atom_ids);
        return;
    }

    const auto & targets=input.selection_domain->target_indices;
    n::Indices selected; selected.reserve(targets.size());
    data.target_to_contributor=targets;
    data.contributor_to_target.assign(input.atom_ids.size(),std::numeric_limits<std::size_t>::max());
    for(std::size_t local=0;local<targets.size();++local)
    {
        selected.push_back(static_cast<Eigen::Index>(targets[local]));
        data.contributor_to_target[targets[local]]=local;
    }
    data.domain=data.domain.Select(selected,static_cast<Eigen::Index>(input.observations.size()),
        n::CompactRowMapping::Identity(static_cast<Eigen::Index>(input.observations.size())));
    data.context.atom_ids=data.context.atom_ids.Select(selected);
    const auto rows=static_cast<Eigen::Index>(input.observations.size());
    const auto atoms=static_cast<Eigen::Index>(targets.size());
    data.context.rank={rows,2*atoms,atoms};
    data.context.linear.rank_relative=data.context.rank.Relative(2*atoms);
    data.layout.full_atoms.resize(targets.size());
    std::iota(data.layout.full_atoms.begin(),data.layout.full_atoms.end(),0);
    data.layout.informative_rows.resize(input.observations.size());
    std::iota(data.layout.informative_rows.begin(),data.layout.informative_rows.end(),0);
    data.partition=n::Partition(data.domain,data.context.atom_ids);
}
struct FullReplay
{
    bool passed{},residual_passed{},halo_passed{},objective_passed{},mapping_passed{},
        observation_identities_passed{},fixed_parameters_passed{},target_feasible{},
        target_ac_passed{},target_width_passed{};
    std::vector<double> prediction;
    double residual_difference{},halo_difference{},objective_difference{},objective{},tolerance{},
        objective_tolerance{},target_ac_kkt{},target_width_gradient{};
    std::string reason;
};
FullReplay ReplayFixedHaloFullModel(const n::ProblemData & data,const n::AssemblyResult & assembly,
    const n::JointConvergenceCriteria & convergence_criteria)
{
    FullReplay out;
    out.fixed_parameters_passed=true;
    const auto & input=*data.input;
    if(!assembly.available || !assembly.raw.valid || assembly.prediction.size()!=data.y.size() ||
        data.target_to_contributor.size()!=data.domain.atoms.size() ||
        data.contributor_to_target.size()!=input.atom_ids.size())
    {out.reason="missing-target-assembly"; return out;}
    out.mapping_passed=true;
    for(std::size_t target=0;target<data.target_to_contributor.size();++target)
    {
        const auto contributor=data.target_to_contributor[target];
        if(contributor>=input.atom_ids.size() || data.contributor_to_target[contributor]!=target ||
            data.context.atom_ids.size()!=data.target_to_contributor.size() ||
            data.context.atom_ids[target]!=input.atom_ids[contributor])
        {out.mapping_passed=false; break;}
    }
    out.observation_identities_passed=data.context.row_ids.size()==input.row_ids.size();
    for(std::size_t row=0;out.observation_identities_passed && row<input.row_ids.size();++row)
        out.observation_identities_passed=data.context.row_ids[row]==input.row_ids[row];
    out.prediction.assign(input.observations.size(),0.0);
    std::vector<double> target_prediction(input.observations.size()),halo_prediction(input.observations.size());
    std::size_t fixed_cursor{};
    const auto unassigned=std::numeric_limits<std::size_t>::max();
    for(std::size_t contributor=0;contributor<input.atom_ids.size();++contributor)
    {
        const auto target=data.contributor_to_target[contributor];
        double amplitude{},width{},charge{};
        if(target!=unassigned)
        {
            if(target>=static_cast<std::size_t>(assembly.raw.eta.size()) ||
                2*target+1>=static_cast<std::size_t>(assembly.raw.beta.size()))
            {out.reason="target-state-mapping-mismatch"; return out;}
            amplitude=assembly.raw.beta(static_cast<Eigen::Index>(2*target));
            charge=assembly.raw.beta(static_cast<Eigen::Index>(2*target+1));
            width=std::exp(assembly.raw.eta(static_cast<Eigen::Index>(target)));
        }
        else
        {
            if(fixed_cursor>=data.fixed_halos.size() ||
                data.fixed_halos[fixed_cursor].contributor_index!=contributor)
            {out.reason="fixed-halo-state-mapping-mismatch"; return out;}
            const auto & fixed=data.fixed_halos[fixed_cursor++];
            out.fixed_parameters_passed &= fixed.atom_id==input.atom_ids[contributor] &&
                std::isfinite(fixed.amplitude) && fixed.amplitude>=0 &&
                std::isfinite(fixed.width) && fixed.width>0 && fixed.charge==0.0 &&
                IsQualifiedFirstStageFixedHalo(fixed) &&
                fixed.source_contract==FixedHaloParameterSourceContractVersion;
            amplitude=fixed.amplitude; width=fixed.width; charge=0;
        }
        if(!std::isfinite(amplitude) || amplitude<0 || !std::isfinite(width) || !(width>0) ||
            !std::isfinite(charge))
        {out.reason="invalid-full-model-parameters"; return out;}
        for(const auto & support:input.support[contributor])
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            const double value=amplitude*basis.gaussian+charge*basis.charge;
            const auto row=support.row;
            const double updated=out.prediction[row]+value;
            if(!std::isfinite(value) || !std::isfinite(updated))
            {out.reason="nonfinite-full-model-prediction"; return out;}
            out.prediction[row]=updated;
            auto & part=target==unassigned ? halo_prediction[row] : target_prediction[row];
            part+=target==unassigned ? amplitude*basis.gaussian : value;
            if(!std::isfinite(part)) {out.reason="nonfinite-full-model-prediction"; return out;}
        }
    }
    if(fixed_cursor!=data.fixed_halos.size()) {out.reason="unused-fixed-halo-state"; return out;}
    if(!data.fixed_halos.empty() && !out.fixed_parameters_passed)
    {out.reason="invalid-fixed-halo-provenance"; return out;}
    if(!data.fixed_halos.empty()) out.fixed_parameters_passed=true;
    double full_rss{},effective_rss{},full_compensation{},effective_compensation{};
    double maximum_scale=1;
    for(std::size_t row=0;row<input.observations.size();++row)
    {
        const double raw_residual=out.prediction[row]-input.observations[row];
        const double effective_residual=target_prediction[row]-data.y(static_cast<Eigen::Index>(row));
        const double expected_full=target_prediction[row]+halo_prediction[row];
        out.residual_difference=std::max(out.residual_difference,std::abs(raw_residual-effective_residual));
        out.halo_difference=std::max(out.halo_difference,
            std::abs(halo_prediction[row]-data.fixed_prediction(static_cast<Eigen::Index>(row))));
        out.halo_difference=std::max(out.halo_difference,
            std::abs(expected_full-assembly.prediction(static_cast<Eigen::Index>(row))-
                data.fixed_prediction(static_cast<Eigen::Index>(row))));
        maximum_scale=std::max({maximum_scale,std::abs(raw_residual),std::abs(effective_residual),
            std::abs(out.prediction[row]),std::abs(input.observations[row])});
        const double full_term=raw_residual*raw_residual-full_compensation;
        const double full_next=full_rss+full_term;
        full_compensation=(full_next-full_rss)-full_term; full_rss=full_next;
        const double effective_term=effective_residual*effective_residual-effective_compensation;
        const double effective_next=effective_rss+effective_term;
        effective_compensation=(effective_next-effective_rss)-effective_term; effective_rss=effective_next;
    }
    out.tolerance=2e-12*maximum_scale;
    const double assembled_objective=.5*effective_rss;
    const double replay_objective=.5*full_rss;
    out.objective_tolerance=2e-12*std::max(1.0,assembled_objective);
    out.objective_difference=std::max(std::abs(replay_objective-assembled_objective),
        std::abs(assembled_objective-assembly.objective));
    out.objective=replay_objective/(data.context.scale*data.context.scale);
    out.residual_passed=out.residual_difference<=out.tolerance;
    out.halo_passed=out.halo_difference<=out.tolerance;
    out.objective_passed=std::isfinite(replay_objective) && out.objective_difference<=out.objective_tolerance;

    const auto target_count=data.target_to_contributor.size();
    std::vector<long double> ac_gradient(2*target_count),column_norm_squared(2*target_count),
        gaussian_width_gradient(target_count),charge_width_gradient(target_count);
    out.target_feasible=true;
    for(std::size_t target=0;target<target_count;++target)
    {
        const auto contributor=data.target_to_contributor[target];
        const double amplitude=assembly.raw.beta(static_cast<Eigen::Index>(2*target));
        const double charge=assembly.raw.beta(static_cast<Eigen::Index>(2*target+1));
        const double width=std::exp(assembly.raw.eta(static_cast<Eigen::Index>(target)));
        out.target_feasible &= std::isfinite(amplitude) && amplitude>=0 && std::isfinite(charge) &&
            std::isfinite(width) && width>0;
        for(const auto & support:input.support[contributor])
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,width,2.5);
            const long double residual=out.prediction[support.row]-input.observations[support.row];
            ac_gradient[2*target]+=static_cast<long double>(basis.gaussian)*residual;
            ac_gradient[2*target+1]+=static_cast<long double>(basis.charge)*residual;
            column_norm_squared[2*target]+=static_cast<long double>(basis.gaussian)*basis.gaussian;
            column_norm_squared[2*target+1]+=static_cast<long double>(basis.charge)*basis.charge;
            gaussian_width_gradient[target]+=static_cast<long double>(basis.gaussian_log_width)*residual;
            charge_width_gradient[target]+=static_cast<long double>(basis.charge_log_width)*residual;
        }
    }
    const long double scale=data.context.scale;
    out.target_ac_kkt=0;
    out.target_width_gradient=0;
    for(std::size_t target=0;target<target_count;++target)
    {
        const double amplitude=assembly.raw.beta(static_cast<Eigen::Index>(2*target));
        const double charge=assembly.raw.beta(static_cast<Eigen::Index>(2*target+1));
        for(std::size_t coefficient=0;coefficient<2;++coefficient)
        {
            const auto index=2*target+coefficient;
            const double norm=std::sqrt(static_cast<double>(column_norm_squared[index]));
            if(!std::isfinite(norm) || !(norm>0))
            {out.target_feasible=false; out.target_ac_kkt=std::numeric_limits<double>::infinity(); continue;}
            const double beta=coefficient==0 ? amplitude : charge;
            const double scaled_beta=norm*beta/static_cast<double>(scale);
            const double scaled_gradient=static_cast<double>(ac_gradient[index]/norm/scale);
            double projected=scaled_beta-scaled_gradient;
            if(coefficient==0) projected=std::max(0.0,projected);
            out.target_ac_kkt=std::max(out.target_ac_kkt,std::abs(scaled_beta-projected));
        }
        const double width_gradient=static_cast<double>((amplitude*gaussian_width_gradient[target]+
            charge*charge_width_gradient[target])/(scale*scale));
        out.target_width_gradient=std::max(out.target_width_gradient,std::abs(width_gradient));
    }
    out.target_ac_passed=out.target_feasible && std::isfinite(out.target_ac_kkt) &&
        out.target_ac_kkt<=convergence_criteria.global_ac_kkt;
    out.target_width_passed=std::isfinite(out.target_width_gradient) &&
        out.target_width_gradient<=convergence_criteria.global_width_gradient;
    out.passed=out.mapping_passed && out.observation_identities_passed && out.fixed_parameters_passed &&
        out.residual_passed && out.halo_passed && out.objective_passed && out.target_ac_passed && out.target_width_passed;
    if(!out.passed) out.reason=!out.mapping_passed ? "target-contributor-mapping-mismatch" :
        !out.observation_identities_passed ? "original-observation-identity-mismatch" :
        !out.fixed_parameters_passed ? "invalid-fixed-halo-provenance" :
        !out.residual_passed ? "full-and-effective-residuals-disagree" :
        !out.halo_passed ? "fixed-halo-replay-disagrees" :
        !out.objective_passed ? "full-and-effective-objectives-disagree" :
        !out.target_feasible ? "target-parameters-infeasible" :
        !out.target_ac_passed ? "target-ac-kkt-failed" : "target-width-gradient-failed";
    return out;
}
}
namespace joint_component {
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
    out.push_back(check("width-stationarity",evidence_available && a.primary.valid && a.primary.gradient.size()>0,
        replay.primary_width_gradient<=a.convergence_criteria.global_width_gradient,
        std::isfinite(replay.primary_width_gradient) ? std::optional<double>(replay.primary_width_gradient) : std::nullopt,
        a.convergence_criteria.global_width_gradient));
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
        std::vector<const char *> assembly_required=reference_free ?
            std::vector<const char *>{"assembly-mapping","assembled-consistency"} :
            std::vector<const char *>{"assembled-profile"};
        const bool fixed_halo_replay=std::any_of(evidence.begin(),evidence.end(),[&](const auto & check) {
            return check.name=="fixed-halo-replay" && check.scope==scope;
        });
        if(fixed_halo_replay)
            for(const char * name:{"full-observation-replay","fixed-halo-replay","full-objective-replay",
                "target-ac-feasibility","target-ac-optimality","target-width-stationarity",
                "target-contributor-mapping","original-observation-identities","fixed-halo-provenance"})
                assembly_required.push_back(name);
        for(const char * name:assembly_required)
        {
            const auto found=std::find_if(evidence.begin(),evidence.end(),[&](const auto & check){return check.name==name && check.scope==scope;});
            status=MergeConvergenceStatus(status,found==evidence.end() ? JointCheckStatus::Unavailable : found->status);
        }
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
    return joint_component::MergeConvergenceStatus(status,stationary==evidence.end() ?
        JointCheckStatus::Unavailable : stationary->status);
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
    if(input.selection_domain && input.selection_domain->contributor_policy==FullABCHaloContributorPolicy)
    {
        const auto & provenance=*input.selection_domain->halo_eligibility;
        for(const auto & excluded:provenance.excluded_halos)
            if(std::find(input.atom_ids.begin(),input.atom_ids.end(),excluded.atom_id)!=input.atom_ids.end())
                throw std::invalid_argument("Excluded halo is present in the retained contributor list.");
        const auto & targets=input.selection_domain->target_indices;
        for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
        {
            if(std::binary_search(targets.begin(),targets.end(),atom)) continue;
            auto support=input.support[atom];
            if(!n::EvaluateFullABCHaloEligibility(support,provenance.minimum_observations).eligible)
                throw std::invalid_argument("Ineligible halo in a FullABC-pruned joint problem.");
        }
    }
    if(input.selection_domain && input.selection_domain->contributor_policy==FixedHaloContributorPolicy)
    {
        const auto & targets=input.selection_domain->target_indices;
        for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
            if(!std::binary_search(targets.begin(),targets.end(),atom) && input.support[atom].empty())
                throw std::invalid_argument("Zero-support halo is present in a fixed-halo joint problem.");
    }
    if(!std::all_of(input.observations.begin(),input.observations.end(),[](double value){return std::isfinite(value);}))
        throw std::invalid_argument("Nonfinite joint observations.");
    for(const auto & atom:input.support) for(const auto & s:atom)
        if(s.row>=input.observations.size() || !std::isfinite(s.squared_distance) || s.squared_distance<0 || s.squared_distance>6.25)
            throw std::invalid_argument("Invalid joint contributor row.");
    const bool fixed_halo_policy=input.selection_domain &&
        input.selection_domain->contributor_policy==FixedHaloContributorPolicy;
    const std::size_t optimization_atoms=fixed_halo_policy ?
        input.selection_domain->target_indices.size() : input.atom_ids.size();
    n::MemoryShape shape;
    shape.component_atoms=optimization_atoms;
    shape.component_observations=input.observations.size();
    shape.problem_scope=true;
    const auto support_entries=n::CountSupportEntries(input);
    const auto resident=n::EstimateSnapshotBytes(input);
    const auto construction_scratch=n::EstimateConstructionScratchBytes(input);
    std::size_t free_columns{};
    const bool dimensions_ok=n::CheckedMultiplyBytes(optimization_atoms,2,free_columns);
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
        if(fixed_halo_policy)
        {
            std::size_t future_fixed_bytes{},optimization_view_bytes{},target_support_entries{};
            for(const auto target:input.selection_domain->target_indices)
                if(!n::CheckedAddBytes(target_support_entries,input.support[target].size(),target_support_entries))
                {future_fixed_bytes=std::numeric_limits<std::size_t>::max(); break;}
            const auto halo_count=input.atom_ids.size()-optimization_atoms;
            if(future_fixed_bytes!=std::numeric_limits<std::size_t>::max() &&
                (!AddEstimatedBytes(future_fixed_bytes,halo_count,sizeof(JointFixedHaloParameter)) ||
                    !AddEstimatedBytes(future_fixed_bytes,input.observations.size(),2*sizeof(double)) ||
                    !AddEstimatedBytes(future_fixed_bytes,input.observations.size(),2*sizeof(Eigen::Index)) ||
                    !AddEstimatedBytes(future_fixed_bytes,optimization_atoms,sizeof(std::vector<n::LocalSupportIndex>)) ||
                    !AddEstimatedBytes(future_fixed_bytes,target_support_entries,sizeof(n::LocalSupportIndex)) ||
                    !AddEstimatedBytes(future_fixed_bytes,optimization_atoms+input.atom_ids.size(),sizeof(std::size_t))))
                future_fixed_bytes=std::numeric_limits<std::size_t>::max();
            if(future_fixed_bytes==std::numeric_limits<std::size_t>::max() ||
                !n::CheckedAddBytes(shape.resident_bytes,future_fixed_bytes,optimization_view_bytes))
            {
                construction_memory.shape=shape;
                construction_memory.stage=n::MemoryStage::ProblemConstruction;
                construction_memory.status=n::MemoryStatus::Overflow;
                construction_memory.reason="size-overflow";
                construction_memory.budget_bytes=n::JointMemoryPlanningBudgetBytes;
            }
            else shape.resident_bytes=optimization_view_bytes;
        }
        shape.identity_validation_bytes=*construction_scratch;
        shape.component_free_columns=free_columns;
        if(construction_memory.status!=n::MemoryStatus::Overflow)
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
const std::vector<JointFixedHaloParameter> & JointProblem::FixedHalos() const {return m_data->fixed_halos;}
JointProblem JointProblem::WithFixedHalos(std::vector<JointFixedHaloParameter> fixed_halos) const
{
    const auto & input=*m_data->input;
    if(!input.selection_domain || input.selection_domain->contributor_policy!=FixedHaloContributorPolicy)
        throw std::invalid_argument("Fixed halos require the fixed-halo contributor policy.");
    const auto & targets=input.selection_domain->target_indices;
    if(fixed_halos.size()!=input.atom_ids.size()-targets.size())
        throw std::invalid_argument("Fixed-halo initialization is incomplete.");
    std::size_t previous{};
    for(std::size_t k=0;k<fixed_halos.size();++k)
    {
        const auto & fixed=fixed_halos[k];
        if(fixed.contributor_index>=input.atom_ids.size() ||
            std::binary_search(targets.begin(),targets.end(),fixed.contributor_index) ||
            input.atom_ids[fixed.contributor_index]!=fixed.atom_id || input.support[fixed.contributor_index].empty() ||
            (k && fixed.contributor_index<=previous) || !std::isfinite(fixed.amplitude) || fixed.amplitude<0 ||
            !std::isfinite(fixed.width) || !(fixed.width>0) || fixed.charge!=0.0 ||
            !IsQualifiedFirstStageFixedHalo(fixed) ||
            fixed.source_contract!=FixedHaloParameterSourceContractVersion)
            throw std::invalid_argument("Invalid fixed-halo parameter snapshot.");
        previous=fixed.contributor_index;
    }
    Eigen::VectorXd fixed_prediction=Eigen::VectorXd::Zero(static_cast<Eigen::Index>(input.observations.size()));
    for(const auto & fixed:fixed_halos)
        for(const auto & support:input.support[fixed.contributor_index])
        {
            const auto basis=n::EvaluateKernel(support.squared_distance,fixed.width,2.5);
            const double contribution=fixed.amplitude*basis.gaussian;
            const auto row=static_cast<Eigen::Index>(support.row);
            const double updated=fixed_prediction(row)+contribution;
            if(!std::isfinite(contribution) || !std::isfinite(updated))
                throw std::invalid_argument("Nonfinite fixed-halo prediction.");
            fixed_prediction(row)=updated;
        }
    auto effective=std::make_shared<std::vector<double>>(input.observations.size());
    for(std::size_t row=0;row<effective->size();++row)
    {
        (*effective)[row]=input.observations[row]-fixed_prediction(static_cast<Eigen::Index>(row));
        if(!std::isfinite((*effective)[row]))
            throw std::invalid_argument("Nonfinite effective observation.");
    }
    std::shared_ptr<const std::vector<double>> effective_snapshot=effective;
    auto data=std::make_shared<n::ProblemData>(m_data->input,std::move(fixed_halos),effective_snapshot);
    data->target_to_contributor=m_data->target_to_contributor;
    data->contributor_to_target=m_data->contributor_to_target;
    data->domain=m_data->domain;
    data->context=m_data->context;
    data->context.observations=n::Observe(effective_snapshot);
    data->layout=m_data->layout;
    data->partition=m_data->partition;
    data->construction_memory=m_data->construction_memory;
    data->fixed_prediction=std::move(fixed_prediction);
    return JointProblem(std::move(data));
}
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
JointProblem BuildJointProblem(const MapObject & map,const ModelObject & model)
{
    const auto atoms=EligibleAtoms(model); JointProblemInput input;
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
        !n::CheckedAddBytes(builder_resident,builder_tracking,builder_resident))
        throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
    if(builder_resident>n::JointMemoryPlanningBudgetBytes)
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

    std::size_t maximum_candidate_support{};
    for(const auto * atom:atoms)
    {
        std::size_t memberships{};
        VisitSupport(map,*atom,[&](auto index,double) {
            if(std::binary_search(rows.begin(),rows.end(),index) &&
                !n::CheckedAddBytes(memberships,1,memberships))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
        });
        maximum_candidate_support=std::max(maximum_candidate_support,memberships);
    }

    std::size_t eligibility_scratch_bytes{},eligibility_builder_resident{};
    if(!n::CheckedMultiplyBytes(maximum_candidate_support,sizeof(JointSupport),eligibility_scratch_bytes) ||
        !n::CheckedAddBytes(builder_resident,eligibility_scratch_bytes,eligibility_builder_resident))
        throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
    if(eligibility_builder_resident>n::JointMemoryPlanningBudgetBytes)
    {
        n::MemoryShape scan_shape; scan_shape.component_atoms=atoms.size();
        scan_shape.component_observations=rows.size(); scan_shape.component_support_entries=maximum_candidate_support;
        scan_shape.resident_bytes=eligibility_builder_resident; scan_shape.problem_scope=true;
        RejectBuildPreflight(scan_shape,"halo-eligibility-scratch-exceeds-budget");
    }

    std::size_t contributors=selected.size(),support_entries{},maximum_atom_support{};
    std::size_t candidate_halos{},retained_halos{},zero_support_halos{};
    std::vector<JointSupport> eligibility_support;
    eligibility_support.reserve(maximum_candidate_support);
    for(const auto * atom:atoms)
    {
        if(selected.contains(atom))
        {
            std::size_t memberships{};
            VisitSupport(map,*atom,[&](auto index,double) {
                if(std::binary_search(rows.begin(),rows.end(),index) &&
                    !n::CheckedAddBytes(memberships,1,memberships))
                    throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
            });
            maximum_atom_support=std::max(maximum_atom_support,memberships);
            if(!n::CheckedAddBytes(support_entries,memberships,support_entries))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
            continue;
        }
        if(!n::CheckedAddBytes(candidate_halos,1,candidate_halos))
            throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
        eligibility_support.clear();
        VisitSupport(map,*atom,[&](auto index,double square) {
            const auto row=std::lower_bound(rows.begin(),rows.end(),index);
            if(row!=rows.end() && *row==index)
                eligibility_support.push_back({static_cast<std::size_t>(row-rows.begin()),square});
        });
        const auto memberships=eligibility_support.size();
        if(memberships==0)
        {
            if(!n::CheckedAddBytes(zero_support_halos,1,zero_support_halos))
                throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
            continue;
        }
        maximum_atom_support=std::max(maximum_atom_support,memberships);
        if(!n::CheckedAddBytes(retained_halos,1,retained_halos) ||
            !n::CheckedAddBytes(contributors,1,contributors) ||
            !n::CheckedAddBytes(support_entries,memberships,support_entries))
            throw std::runtime_error("resource-limit:build-joint-problem:size-overflow");
    }
    eligibility_support.clear();
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
        !n::CheckedAddBytes(identity_scratch,eligibility_scratch_bytes,identity_scratch) ||
        !n::CheckedAddBytes(contributors,selected.size(),identity_map_entries) ||
        !AddEstimatedBytes(planned_snapshot,identity_characters,2*sizeof(char)))
        RejectBuildPreflight(shape,"size-overflow");
    shape.identity_validation_bytes=identity_scratch;
    if(!AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(double)) ||
        !AddEstimatedBytes(planned_snapshot,contributors,2*sizeof(std::vector<JointSupport>)) ||
        !AddEstimatedBytes(planned_snapshot,contributors,2*sizeof(std::string)) ||
        !AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(std::string)) ||
        !AddEstimatedBytes(planned_snapshot,support_entries,2*sizeof(JointSupport)) ||
        !AddEstimatedBytes(planned_snapshot,selected.size(),2*sizeof(std::size_t)) ||
        !AddEstimatedBytes(planned_snapshot,retained_halos,sizeof(JointFixedHaloParameter)) ||
        !AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(double)) ||
        !AddEstimatedBytes(planned_snapshot,rows.size(),2*sizeof(Eigen::Index)) ||
        !AddEstimatedBytes(planned_snapshot,selected.size(),sizeof(std::vector<n::LocalSupportIndex>)) ||
        !AddEstimatedBytes(planned_snapshot,selected_support_entries,sizeof(n::LocalSupportIndex)) ||
        !AddEstimatedBytes(planned_snapshot,identity_map_entries,sizeof(std::size_t)) ||
        !n::CheckedAddBytes(planned_snapshot,builder_resident,planned_snapshot))
        RejectBuildPreflight(shape,"size-overflow");
    shape.resident_bytes=planned_snapshot;
    const auto build_estimate=n::PlanMemory(shape,n::MemoryStage::ProblemConstruction);
    if(build_estimate.status!=n::MemoryStatus::Pass)
        RejectBuildPreflight(shape,build_estimate.reason.c_str());

    for(auto index:rows) {input.row_ids.push_back(std::to_string(index)); input.observations.push_back(map.GetMapValue(index));}
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(FixedHaloContributorPolicy);
    input.selection_domain->fixed_halo_selection=JointFixedHaloSelectionProvenance{
        1,candidate_halos,retained_halos,zero_support_halos};
    for(const auto * atom:atoms)
    {
        bool selected_target=selected.contains(atom);
        std::vector<JointSupport> support;
        if(selected_target)
        {
            VisitSupport(map,*atom,[&](auto index,double square) {
                const auto row=std::lower_bound(rows.begin(),rows.end(),index);
                if(row!=rows.end() && *row==index) support.push_back({static_cast<std::size_t>(row-rows.begin()),square});
            });
        }
        else
        {
            eligibility_support.clear();
            VisitSupport(map,*atom,[&](auto index,double square) {
                const auto row=std::lower_bound(rows.begin(),rows.end(),index);
                if(row!=rows.end() && *row==index)
                    eligibility_support.push_back({static_cast<std::size_t>(row-rows.begin()),square});
            });
            support=eligibility_support;
        }
        if(selected_target || !support.empty())
        {
            if(selected_target) input.selection_domain->target_indices.push_back(input.atom_ids.size());
            input.atom_ids.push_back(std::to_string(atom->GetSerialID()));
            input.support.push_back(std::move(support));
        }
    }
    return JointProblem(std::move(input));
}
JointProblem BuildPureHaloSupportExclusionProblem(const MapObject & map,const ModelObject & model)
{
    const auto fixed_problem=BuildJointProblem(map,model);
    const auto & source=fixed_problem.Input();
    if (!source.selection_domain || source.selection_domain->contributor_policy!=FixedHaloContributorPolicy ||
        !source.selection_domain->fixed_halo_selection)
        throw std::invalid_argument("Pure Halo exclusion requires the selected voxel closure.");
    const auto & targets=source.selection_domain->target_indices;
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
    auto original=source;
    original.selection_domain->contributor_policy=std::string(PureHaloSupportExclusionContributorPolicy);
    original.selection_domain->fixed_halo_selection.reset();
    const auto fixed=&*source.selection_domain->fixed_halo_selection;
    original.selection_domain->pure_halo_support_exclusion=JointPureHaloSupportExclusionProvenance{
        source.observations.size(),retained_count,source.observations.size()-retained_count,
        fixed->candidate_halos,fixed->retained_halos,fixed->zero_support_halos};
    std::vector<std::size_t> fitting_to_original_contributors=targets;
    return JointProblem(std::move(fitting)).WithReconstructionDomain(
        std::move(original),std::move(fitting_to_original_contributors),std::move(fitting_rows));
}
JointFitResult FitJointComponents(const JointProblem & problem,const std::vector<double> & initial_b)
{
    const auto & input=problem.Input();
    if(input.selection_domain && input.selection_domain->contributor_policy==FixedHaloContributorPolicy &&
        problem.FixedHalos().size()!=input.atom_ids.size()-input.selection_domain->target_indices.size())
    {
        JointFitResult out; out.problem=problem; out.initialization.b=initial_b;
        out.initialization.valid=false; out.initialization.reason="fixed-halo-initialization-unavailable";
        out.observation_scale=problem.ObservationScale();
        out.available_row_mask.assign(input.observations.size(),false);
        out.evidence.push_back({"fixed-halo-initialization",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},out.initialization.reason});
        return out;
    }
    return n::FitFixedNeighborComponents(problem,initial_b,{});
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
        out.evidence.push_back({"fixed-halo-initialization",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{}, {},out.initialization.reason});
        return out;
    };
    if(!input.selection_domain || input.selection_domain->contributor_policy!=FixedHaloContributorPolicy)
        return invalid("fixed-halo-initialization-unavailable:missing-fixed-halo-policy");
    if(initialization.atoms.size()!=input.atom_ids.size() || initialization.b.size()!=input.atom_ids.size())
        return invalid("fixed-halo-initialization-unavailable:contributor-mapping-mismatch");
    if(!initialization.valid)
        return invalid(initialization.reason.starts_with("fixed-halo-initialization-unavailable:") ?
            initialization.reason : "fixed-halo-initialization-unavailable:"+initialization.reason);

    const auto & targets=input.selection_domain->target_indices;
    std::vector<JointFixedHaloParameter> fixed_halos;
    fixed_halos.reserve(input.atom_ids.size()-targets.size());
    for(std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        const auto & record=initialization.atoms[atom];
        if(record.id!=input.atom_ids[atom])
            return invalid("fixed-halo-initialization-unavailable:contributor-identity-mismatch");
        const bool target=std::binary_search(targets.begin(),targets.end(),atom);
        if(record.target!=target)
            return invalid("fixed-halo-initialization-unavailable:contributor-role-mismatch");
        if(target) continue;
        std::string failure;
        auto fixed=detail::MakeFixedHaloParameter(record,atom,input.atom_ids[atom],failure);
        if(!fixed) return invalid("fixed-halo-initialization-unavailable:"+failure);
        fixed_halos.push_back(std::move(*fixed));
    }
    JointProblem initialized_problem=problem;
    try {initialized_problem=problem.WithFixedHalos(std::move(fixed_halos));}
    catch(const std::exception & error)
    {return invalid(std::string("fixed-halo-initialization-unavailable:")+error.what());}
    std::vector<double> target_widths; target_widths.reserve(targets.size());
    for(const auto contributor:targets) target_widths.push_back(initialization.b[contributor]);
    auto out=n::FitFixedNeighborComponents(initialized_problem,target_widths,{},observer);
    out.initialization=initialization;
    return out;
}
JointFitResult n::FitFixedNeighborComponentsImpl(const JointProblem & problem,const std::vector<double> & initial_b,
    const n::FixedNeighborSearchPolicy & production_policy,const n::JointProgressObserver & observer)
{
    eigen_helper::ScopedEigenThreadCount eigen_thread_guard{1};
    const auto configuration=n::ResolveJointSolverConfiguration(production_policy);
    n::NotifyJointSolverConfigured(observer,configuration);
    auto provenance=[&] {
        JointSolverProvenance out;
        out.search_method="fixed-neighbor";
        out.contract_version=std::string(JointSolverProvenanceContractVersion);
        out.sparse_backend="EIGEN";
        out.fixed_neighbor_core_atoms=production_policy.core_atoms;
        out.fixed_neighbor_policy_version=std::string(FixedNeighborPolicyContractVersion);
        out.fixed_neighbor_maximum_sweeps=production_policy.maximum_sweeps;
        out.fixed_neighbor_order=std::string(n::FixedNeighborBlockOrderName(production_policy.order));
        out.fixed_neighbor_local_search=std::string(FixedNeighborLocalSearchName);
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
            data.construction_memory.estimated_peak_bytes ? std::optional<double>(static_cast<double>(data.construction_memory.estimated_peak_bytes)) : std::nullopt,
            static_cast<double>(data.construction_memory.budget_bytes),
            "resource-limit:problem-construction:"+data.construction_memory.reason});
        return out;
    }
    const bool fixed_halo_policy=data.input->selection_domain &&
        data.input->selection_domain->contributor_policy==FixedHaloContributorPolicy;
    if(fixed_halo_policy && data.fixed_halos.size()!=
        data.input->atom_ids.size()-data.input->selection_domain->target_indices.size())
    {
        JointFitResult out; out.problem=problem; out.solver_provenance=provenance();
        out.observation_scale=data.context.scale; out.initialization.b=initial_b;
        out.initialization.valid=false;
        out.initialization.reason="fixed-halo-initialization-unavailable:missing-validated-first-stage-snapshot";
        out.available_row_mask.assign(data.input->observations.size(),false);
        out.evidence.push_back({"fixed-halo-initialization",JointCheckStatus::Unavailable,
            JointEvidenceScope::AssembledGlobal,{},{},out.initialization.reason});
        return out;
    }
    JointFitResult out; out.problem=problem;
    out.solver_provenance=provenance();
    out.observation_scale=data.context.scale; out.initialization.b=initial_b;
    out.initialization.valid=initial_b.size()==data.domain.atoms.size() && std::all_of(initial_b.begin(),initial_b.end(),[](double b){return std::isfinite(b) && b>0;});
    out.initialization.reason=out.initialization.valid ? "valid-widths" : "invalid-widths";
    if(initial_b.size()!=data.domain.atoms.size())
    {
        out.available_row_mask.assign(data.input->observations.size(),false);
        for(auto row:data.partition.constant_rows) out.available_row_mask[static_cast<std::size_t>(row)]=true;
        return out;
    }
    const auto component_count=data.partition.components.size();
    auto bounded_policy=production_policy;
    std::size_t remaining_perturbation_blocks=n::MaximumPerturbedBlocksPerFit;
    const auto request_perturbation_block=production_policy.claim_parameter_stability_block;
    bounded_policy.claim_parameter_stability_block=[&remaining_perturbation_blocks,request_perturbation_block] {
        if(!remaining_perturbation_blocks ||
            (request_perturbation_block && !request_perturbation_block())) return false;
        --remaining_perturbation_blocks; return true;
    };
    auto component_batch=n::OrchestrateFixedNeighborComponents(problem,initial_b,bounded_policy,observer);
    out.components=std::move(component_batch.components);
    out.costs.search_seconds+=component_batch.costs.search_seconds;
    out.costs.assessment_seconds+=component_batch.costs.assessment_seconds;
    const auto assembly_start=Clock::now();
    n::NotifyJointAssemblyProgress(observer,n::JointProgressPhase::AssemblyStarted,component_count);
    std::optional<n::EvaluationContext> reuse_context;
    std::optional<n::AssessmentReuse> reuse;
    if(component_batch.single_component_assessment)
    {
        const auto & view=data.partition.components[0];
        // A single full component has identity row/atom mappings in this immutable snapshot.
        reuse_context=n::ChildContext(data.context,view,true);
        reuse.emplace(n::AssessmentReuse{view.domain,data.y,*reuse_context,*component_batch.single_component_assessment});
    }
    const auto assembly=n::AssembleComponents(data.domain,data.y,data.partition,data.context,out.components,
        reuse ? &*reuse : nullptr,bounded_policy.convergence_criteria);
    out.search_completed=assembly.completed; out.available_row_mask=assembly.row_mask;
    out.evidence=n::AssessmentEvidence(assembly.assessment,JointEvidenceScope::AssembledGlobal);
    out.evidence.push_back({"assembly-mapping",assembly.mapping_valid ? JointCheckStatus::Passed :
        JointCheckStatus::Unavailable,JointEvidenceScope::AssembledGlobal,{},{},
        assembly.mapping_valid ? "exact-component-identity-mapping" :
            assembly.mapping_reason.empty() ? "missing-component-state" : assembly.mapping_reason});
    std::optional<FullReplay> full_replay;
    if(assembly.available)
    {
        const auto & state=assembly.raw;
        if(fixed_halo_policy)
        {
            full_replay=ReplayFixedHaloFullModel(data,assembly,bounded_policy.convergence_criteria);
            out.evidence.push_back({"full-observation-replay",full_replay->residual_passed ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->residual_difference,full_replay->tolerance,full_replay->residual_passed ?
                    "full-and-effective-residuals-agree" : full_replay->reason});
            out.evidence.push_back({"fixed-halo-replay",full_replay->halo_passed ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->halo_difference,full_replay->tolerance,full_replay->halo_passed ?
                    "immutable-first-stage-snapshot" : full_replay->reason});
            out.evidence.push_back({"full-objective-replay",full_replay->objective_passed ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->objective_difference,full_replay->objective_tolerance,full_replay->objective_passed ?
                    "full-and-effective-objectives-agree" : full_replay->reason});
            out.evidence.push_back({"target-ac-feasibility",full_replay->target_feasible ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->target_feasible ? 0.0 : 1.0,0.0,full_replay->target_feasible ?
                    "target-parameters-feasible" : full_replay->reason});
            out.evidence.push_back({"target-ac-optimality",full_replay->target_ac_passed ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->target_ac_kkt,bounded_policy.convergence_criteria.global_ac_kkt,full_replay->target_ac_passed ?
                    "target-only-projected-kkt" : full_replay->reason});
            out.evidence.push_back({"target-width-stationarity",full_replay->target_width_passed ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->target_width_gradient,bounded_policy.convergence_criteria.global_width_gradient,
                full_replay->target_width_passed ?
                    "target-only-log-width-gradient" : full_replay->reason});
            out.evidence.push_back({"target-contributor-mapping",full_replay->mapping_passed ? JointCheckStatus::Passed :
                JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->mapping_passed ? 0.0 : 1.0,0.0,full_replay->mapping_passed ?
                    "target-contributor-identity" : full_replay->reason});
            out.evidence.push_back({"original-observation-identities",full_replay->observation_identities_passed ?
                JointCheckStatus::Passed : JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->observation_identities_passed ? 0.0 : 1.0,0.0,full_replay->observation_identities_passed ?
                    "original-row-identity" : full_replay->reason});
            out.evidence.push_back({"fixed-halo-provenance",full_replay->fixed_parameters_passed ?
                JointCheckStatus::Passed : JointCheckStatus::Failed,JointEvidenceScope::AssembledGlobal,
                full_replay->fixed_parameters_passed ? 0.0 : 1.0,0.0,full_replay->fixed_parameters_passed ?
                    "qualified-first-stage-fixed-parameters" : full_replay->reason});
        }
        const bool replay_passed=assembly.consistency_passed && (!full_replay || full_replay->passed);
        const auto consistency_status=replay_passed ? JointCheckStatus::Passed :
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
        if(assembly.assessment.numerical_qualified && replay_passed)
        {
            out.assembled_state=State(state,data.context.scale);
            out.prediction=full_replay ? std::move(full_replay->prediction) : Values(assembly.prediction);
            out.objective=full_replay ? full_replay->objective :
                assembly.objective/(data.context.scale*data.context.scale);
            if(full_replay) out.assembled_state->objective=*out.objective;
        }
    }
    if(data.input->selection_domain && data.input->selection_domain->contributor_policy==FixedHaloContributorPolicy)
        for(auto & component:out.components)
        {
            const auto remap=[&](std::size_t & atom) {
                if(atom>=data.target_to_contributor.size()) return false;
                atom=data.target_to_contributor[atom]; return true;
            };
            for(auto & atom:component.atoms) remap(atom);
            const auto remap_blocks=[&](auto & blocks) {
                for(auto & block:blocks) for(auto & atom:block.core_atoms) remap(atom);
            };
            remap_blocks(component.local_conditioning);
            if(component.parameter_stability) remap_blocks(component.parameter_stability->blocks);
        }
    out.costs.assembly_seconds=Seconds(assembly_start);
    n::NotifyJointAssemblyProgress(observer,n::JointProgressPhase::AssemblyCompleted,component_count,
        out.costs.assembly_seconds,assembly.available && assembly.raw.valid);
    return out;
}
JointAnalysisResult CaptureJointAnalysisResult(const JointFitResult & fit, JointAnalysisMetadata metadata)
{
    if (!fit.problem) throw std::invalid_argument("Cannot capture joint result without problem identities.");
    JointAnalysisResult out;
    metadata.software=JointSoftwareProvenance{RHBM_GEM_SIMULATION_VERSION,
        RHBM_GEM_SIMULATION_SOURCE_SHA256,RHBM_GEM_SIMULATION_CONFIG_SHA256,
        RHBM_GEM_SIMULATION_BUILD_SHA256};
    out.metadata=std::move(metadata);
    out.metadata.solver=fit.solver_provenance.search_method.empty() ? std::nullopt :
        std::optional<JointSolverProvenance>{fit.solver_provenance};
    out.atom_ids=fit.problem->Input().atom_ids; out.row_ids=fit.problem->Input().row_ids;
    out.selection_domain=fit.problem->Input().selection_domain;
    out.fixed_halos=fit.problem->FixedHalos();
    if(out.selection_domain && out.selection_domain->contributor_policy==FixedHaloContributorPolicy)
        out.parameterization_contract=std::string(JointTargetOnlyFixedHaloParameterizationContract);
    out.initialization=fit.initialization;
    if(out.selection_domain && out.selection_domain->contributor_policy==FixedHaloContributorPolicy &&
        out.initialization.b.size()==out.atom_ids.size())
    {
        std::vector<double> target_widths; target_widths.reserve(out.selection_domain->target_indices.size());
        for(const auto target:out.selection_domain->target_indices) target_widths.push_back(out.initialization.b[target]);
        out.initialization.b=std::move(target_widths);
    }
    out.costs=fit.costs;
    for (const auto & component:fit.components)
    {
        JointAnalysisComponent saved;
        static_cast<JointComponentData &>(saved)=component;
        saved.runtime_convergence=component.RuntimeConvergence();
        saved.target_runtime_convergence=component.TargetRuntimeConvergence();
        out.components.push_back(std::move(saved));
    }
    out.assembled_state=fit.assembled_state; out.objective=fit.objective;
    out.available_row_mask=fit.available_row_mask; out.evidence=fit.evidence;
    out.search_completed=fit.search_completed; out.observation_scale=fit.observation_scale;
    out.runtime_convergence=fit.RuntimeConvergence(); out.regular_certificate=fit.regular_certificate;
    out.target_evidence=fit.target_evidence; out.target_runtime_convergence=fit.TargetRuntimeConvergence();
    return out;
}
JointInitialization joint_component::InitializeContributors(MapObject & map,ModelObject & model,const JointProblem & problem)
{
    ModelObject initializer(model);
    const auto workset = detail::MakeJointFittingWorkset(initializer, problem);
    FitOptions options; options.thread_size = 1; options.quiet_mode = true; options.exclude_hydrogen = true;
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
    const auto problem=BuildJointProblem(map,model);
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
