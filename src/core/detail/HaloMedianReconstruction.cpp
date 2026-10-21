#include "HaloMedianReconstruction.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include "joint_component/Problem.hpp"
#include "joint_component/WidthCalibration.hpp"
#include "PostFitPeeling.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace rhbm_gem::core::detail {
namespace {
struct TargetNode
{
    const HaloMedianTargetDonor * donor{};
    std::array<double,3> position{};
    const std::array<double,3> & GetPosition() const { return position; }
};

std::array<double,3> Point(const GaussianModel3D & point)
{
    return {point.GetAmplitude(),point.GetWidth(),point.GetOffset()};
}

double SquaredDistance(const std::array<double,3> & left,const std::array<double,3> & right)
{
    double square{};
    for (std::size_t axis=0;axis<3;++axis)
    {
        const auto delta=left[axis]-right[axis];
        square+=delta*delta;
    }
    return square;
}

std::array<double,3> Median(const std::vector<const HaloMedianTargetDonor *> & donors)
{
    if (donors.empty()) throw std::invalid_argument("Cannot compute Halo median from empty donors.");
    std::array<double,3> result{};
    for (std::size_t parameter=0;parameter<3;++parameter)
    {
        std::vector<double> values;
        values.reserve(donors.size());
        for (const auto * donor:donors)
        {
            const auto value=donor->point[parameter];
            if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite Halo median donor.");
            values.push_back(value);
        }
        result[parameter]=array_helper::ComputeMedian(values);
    }
    return result;
}

std::array<double,3> Dispersion(
    const std::vector<const HaloMedianTargetDonor *> & donors,const std::array<double,3> & median)
{
    std::array<double,3> result{};
    for (std::size_t parameter=0;parameter<3;++parameter)
        for (const auto * donor:donors)
            result[parameter]=std::max(result[parameter],std::abs(donor->point[parameter]-median[parameter]));
    return result;
}

std::map<std::size_t,GaussianModel3D> TargetPoints(const JointAnalysisResult & result)
{
    std::map<std::size_t,GaussianModel3D> points;
    for (const auto & component:result.components)
    {
        if (!component.state) continue;
        const auto & full=component.layout ? component.layout->full_atoms : component.atoms;
        if (component.state->ac.size()!=2*full.size() || component.state->b.size()!=full.size())
            throw std::invalid_argument("Invalid Joint target state for Halo median reconstruction.");
        for (std::size_t k=0;k<full.size();++k)
            if (!points.emplace(full[k],GaussianModel3D{component.state->ac[2*k],component.state->b[k],
                component.state->ac[2*k+1]}).second)
                throw std::invalid_argument("Duplicate Joint target state for Halo median reconstruction.");
    }
    return points;
}

bool QualifiedTarget(const JointAnalysisResult & result,std::size_t atom)
{
    if (result.runtime_convergence!=JointCheckStatus::Passed ||
        result.target_runtime_convergence!=JointCheckStatus::Passed || !result.target_evidence)
        return false;
    const auto found=std::find_if(result.target_evidence->point_atoms.begin(),result.target_evidence->point_atoms.end(),
        [atom](const auto & point){return point.atom==atom;});
    return found!=result.target_evidence->point_atoms.end() &&
        found->component_convergence==JointCheckStatus::Passed &&
        found->parameter_validity==JointCheckStatus::Passed &&
        found->local_conditioning=="Adequate";
}

bool SameSelectionDomain(const JointSelectionDomain & left,const JointSelectionDomain & right)
{
    const bool same_pure=left.pure_halo_support_exclusion.has_value()==right.pure_halo_support_exclusion.has_value() &&
        (!left.pure_halo_support_exclusion || *left.pure_halo_support_exclusion==*right.pure_halo_support_exclusion);
    const bool same_fixed=left.fixed_halo_selection.has_value()==right.fixed_halo_selection.has_value() &&
        (!left.fixed_halo_selection || *left.fixed_halo_selection==*right.fixed_halo_selection);
    return left.contract==right.contract && left.target_indices==right.target_indices &&
        left.observation_radius==right.observation_radius && left.support_radius==right.support_radius &&
        left.contributor_policy==right.contributor_policy &&
        left.halo_eligibility.has_value()==right.halo_eligibility.has_value() && same_fixed && same_pure;
}

bool ValidPoint(const GaussianModel3D & point)
{
    return std::isfinite(point.GetAmplitude()) && point.GetAmplitude()>=0 &&
        std::isfinite(point.GetWidth()) && point.GetWidth()>0 && std::isfinite(point.GetOffset());
}
}

struct GeometryLocalTargetDonorIndex::Impl
{
    const std::vector<HaloMedianTargetDonor> & donors;
    std::vector<TargetNode> nodes;
    std::unique_ptr<KDNode<TargetNode>> tree;

    explicit Impl(const std::vector<HaloMedianTargetDonor> & source):donors(source)
    {
        nodes.reserve(donors.size());
        for (const auto & donor:donors) nodes.push_back({&donor,donor.position});
        tree=KDTreeAlgorithm<TargetNode>::BuildKDTree(nodes);
    }
};

GeometryLocalTargetDonorIndex::GeometryLocalTargetDonorIndex(
    const std::vector<HaloMedianTargetDonor> & donors):m_impl(std::make_unique<Impl>(donors))
{
}

GeometryLocalTargetDonorIndex::~GeometryLocalTargetDonorIndex()=default;

std::vector<const HaloMedianTargetDonor *> GeometryLocalTargetDonorIndex::Query(
    const std::array<double,3> & position) const
{
    const auto & nodes=m_impl->nodes;
    if (nodes.empty()) return {};
    TargetNode query{nullptr,position};
    const auto count=std::min(GeometryLocalTargetMedianNeighborCount,nodes.size());
    const auto nearest=KDTreeAlgorithm<TargetNode>::KNearestNeighbors(m_impl->tree.get(),&query,count);
    if (nearest.empty()) return {};
    double radius{};
    for (const auto * donor:nearest)
        radius=std::max(radius,std::sqrt(SquaredDistance(position,donor->position)));
    auto candidates=KDTreeAlgorithm<TargetNode>::RangeSearch(m_impl->tree.get(),&query,radius+1.0e-12);
    std::sort(candidates.begin(),candidates.end(),[&](const auto * left,const auto * right) {
        const auto left_distance=SquaredDistance(position,left->position);
        const auto right_distance=SquaredDistance(position,right->position);
        if (left_distance!=right_distance) return left_distance<right_distance;
        return left->donor->atom_id<right->donor->atom_id;
    });
    if (candidates.size()>count) candidates.resize(count);
    std::vector<const HaloMedianTargetDonor *> selected;
    selected.reserve(candidates.size());
    for (const auto * candidate:candidates) selected.push_back(candidate->donor);
    return selected;
}

std::array<double,3> GeometryLocalTargetMedian(
    const std::vector<const HaloMedianTargetDonor *> & donors)
{
    return Median(donors);
}

std::vector<std::size_t> DistinctOriginalSupportRows(const std::vector<JointSupport> & support)
{
    std::vector<std::size_t> rows;
    rows.reserve(support.size());
    for(const auto & entry:support) rows.push_back(entry.row);
    std::sort(rows.begin(),rows.end());
    rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    return rows;
}

std::vector<HaloMedianTargetDonor> QualifiedJointTargetDonors(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result)
{
    const auto & input=problem.OriginalInput();
    if (!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    if (result.atom_ids!=input.atom_ids || !result.selection_domain ||
        !SameSelectionDomain(*result.selection_domain,*input.selection_domain))
        throw std::invalid_argument("halo-median-joint-run-identity-mismatch");
    const auto points=TargetPoints(result);
    std::vector<HaloMedianTargetDonor> out;
    for (const auto target:input.selection_domain->target_indices)
    {
        if (!QualifiedTarget(result,target)) continue;
        const auto found=points.find(target);
        if (found==points.end() || !ValidPoint(found->second)) continue;
        const auto atom_id=input.atom_ids[target];
        const auto * atom=model.FindAtomPtr(std::stoi(atom_id));
        if (!atom) throw std::invalid_argument("halo-median-unknown-target");
        const auto position=atom->GetPosition();
        if (!std::all_of(position.begin(),position.end(),[](double value){return std::isfinite(value);}))
            throw std::invalid_argument("halo-median-nonfinite-target-position");
        out.push_back({target,atom_id,position,Point(found->second)});
    }
    std::sort(out.begin(),out.end(),[](const auto & left,const auto & right) {
        return left.atom_id<right.atom_id;
    });
    return out;
}

HaloFullABCPreclassification BuildHaloFullABCPreclassification(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result)
{
    const auto & input=problem.OriginalInput();
    if(!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    const auto donors=QualifiedJointTargetDonors(model,problem,result);
    HaloFullABCPreclassification out;
    if(donors.empty()) out.unavailable_reason="no-qualified-joint-target-donors";
    else
    {
        std::vector<const HaloMedianTargetDonor *> donor_refs;
        donor_refs.reserve(donors.size());
        for(const auto & donor:donors) donor_refs.push_back(&donor);
        const auto median=GeometryLocalTargetMedian(donor_refs);
        if(!std::isfinite(median[0]) || median[0]<0 || !std::isfinite(median[1]) || !(median[1]>0) ||
            !std::isfinite(median[2]))
            out.unavailable_reason="invalid-qualified-joint-target-median";
        else
        {
            out.available=true;
            out.global_target_median=median;
        }
    }

    const auto contributor_universe=BuildJointPeelingContributorUniverse(model,problem,&result);
    out.halos.reserve(contributor_universe.contributors.size());
    for(const auto & contributor:contributor_universe.contributors)
    {
        if(!contributor.reconstruction_eligible || contributor.target) continue;
        HaloFullABCClassificationRecord record;
        record.atom_id=std::to_string(contributor.atom_id);
        if(contributor.contributor_index)
        {
            record.original_contributor_index=contributor.contributor_index;
            record.original_support_rows=DistinctOriginalSupportRows(input.support.at(*contributor.contributor_index));
            for(const auto row:record.original_support_rows)
            {
                if(row>=input.row_ids.size())
                    throw std::invalid_argument("Invalid Halo observation support row.");
            }
        }
        record.distinct_original_observations=record.original_support_rows.size();
        if(record.distinct_original_observations<3)
        {
            record.role=HaloFullABCRole::FixedGlobalTargetMedian;
            ++out.fixed_halo_count;
        }
        else
        {
            record.role=HaloFullABCRole::FreeFullABC;
            ++out.free_halo_count;
        }
        if(out.available) record.parameter_profile=*out.global_target_median;
        out.halos.push_back(std::move(record));
    }
    return out;
}

HaloFullABCReconstructionResult EstimateHaloFullABCJoint(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result,
    std::string_view joint_run_id,bool collect_solver_diagnostics)
{
    HaloFullABCReconstructionResult out;
    out.joint_run_id=std::string(joint_run_id);
    const auto & original=problem.OriginalInput();
    out.diagnostics.original_observation_rows=original.observations.size();
    if (result.parameterization_contract!=JointPureHaloSupportExclusionParameterizationContract ||
        !IsQualifiedJointAnalysisResult(result))
    {
        out.unavailable_reason="first-joint-endpoint-not-qualified-pure-exclusion";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    if (!original.selection_domain || result.atom_ids!=original.atom_ids || result.row_ids!=original.row_ids ||
        !result.selection_domain || !SameSelectionDomain(*result.selection_domain,*original.selection_domain))
    {
        out.unavailable_reason="first-joint-original-domain-identity-mismatch";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }

    const auto preclassification=BuildHaloFullABCPreclassification(model,problem,result);
    out.fixed_halo_count=preclassification.fixed_halo_count;
    out.free_halo_count=preclassification.free_halo_count;
    out.global_target_median=preclassification.global_target_median;
    out.halos.reserve(preclassification.halos.size());
    for (const auto & halo:preclassification.halos)
    {
        HaloFullABCReconstructionHalo reconstruction;
        reconstruction.classification=halo;
        if(halo.role==HaloFullABCRole::FixedGlobalTargetMedian)
        {
            reconstruction.reconstructed_parameters=halo.parameter_profile;
            reconstruction.parameter_source=HaloFullABCParameterSource::PreFixedGlobalMedian;
            reconstruction.qualified=halo.parameter_profile.has_value();
            reconstruction.qualification_reason=reconstruction.qualified ? "preclassified-before-solve" :
                "global-target-median-unavailable";
        }
        out.halos.push_back(std::move(reconstruction));
    }

    if (out.halos.empty())
    {
        out.available=true;
        out.diagnostics.reconstruction_qualified=true;
        out.global_target_median.reset();
        out.diagnostics.solver_status=JointCheckStatus::NotRun;
        return out;
    }
    if (!preclassification.available || !preclassification.global_target_median)
    {
        out.unavailable_reason=preclassification.unavailable_reason.empty() ?
            "global-target-median-unavailable" : preclassification.unavailable_reason;
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    if (out.free_halo_count==0)
    {
        out.available=true;
        out.diagnostics.reconstruction_qualified=true;
        out.diagnostics.solver_status=JointCheckStatus::NotRun;
        return out;
    }

    const auto calibration=joint_component::CalibrateGlobalWidth(result.initialization,
        "qualified-first-stage-target-widths");
    out.diagnostics.qualified_first_stage_width_donors=calibration.qualified_donor_count;
    out.diagnostics.first_stage_width_p95=calibration.p95;
    out.diagnostics.width_upper_bound=calibration.bmax4;
    if (calibration.status!="Available" || !calibration.eta_max4 || !calibration.bmax4 ||
        !calibration.p95)
    {
        out.unavailable_reason="first-stage-width-bound-unavailable:"+calibration.unavailable_reason;
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }

    const auto target_points=TargetPoints(result);
    const auto & targets=original.selection_domain->target_indices;
    for (const auto target:targets)
    {
        const auto point=target_points.find(target);
        if (point==target_points.end() || !ValidPoint(point->second))
        {
            out.unavailable_reason="first-joint-target-endpoint-missing-or-invalid";
            out.diagnostics.unavailable_reason=out.unavailable_reason;
            return out;
        }
    }

    JointProblemInput halo_input;
    halo_input.row_ids=original.row_ids;
    halo_input.observations=original.observations;
    std::vector<bool> free_sensitive_rows(original.observations.size(),false);
    std::vector<double> initial_widths;
    initial_widths.reserve(out.free_halo_count);
    const auto unique_support=[&](std::size_t contributor) {
        std::vector<JointSupport> support;
        std::map<std::size_t,double> by_row;
        for (const auto & entry:original.support.at(contributor))
        {
            const auto [found,inserted]=by_row.emplace(entry.row,entry.squared_distance);
            if (!inserted && found->second!=entry.squared_distance)
                throw std::invalid_argument("Conflicting duplicate Halo support row.");
        }
        support.reserve(by_row.size());
        for (const auto [row,square]:by_row) support.push_back({row,square});
        return support;
    };
    const auto subtract_profile=[&](std::size_t contributor,const std::array<double,3> & parameters) {
        if (contributor>=original.support.size())
            throw std::invalid_argument("Halo contributor outside OriginalInput.");
        for (const auto & entry:unique_support(contributor))
        {
            const auto basis=joint_component::EvaluateKernel(entry.squared_distance,parameters[1],2.5);
            halo_input.observations.at(entry.row)-=parameters[0]*basis.gaussian+parameters[2]*basis.charge;
        }
    };

    for (const auto target:targets)
        for (const auto & entry:unique_support(target))
        {
            const auto & point=target_points.at(target);
            const auto basis=joint_component::EvaluateKernel(entry.squared_distance,point.GetWidth(),2.5);
            halo_input.observations.at(entry.row)-=
                point.GetAmplitude()*basis.gaussian+point.GetOffset()*basis.charge;
        }

    for (auto & halo:out.halos)
    {
        if (halo.classification.role==HaloFullABCRole::FixedGlobalTargetMedian)
        {
            if (!halo.reconstructed_parameters)
            {
                out.unavailable_reason="fixed-halo-global-target-median-unavailable";
                out.diagnostics.unavailable_reason=out.unavailable_reason;
                return out;
            }
            if (halo.classification.original_contributor_index)
                subtract_profile(*halo.classification.original_contributor_index,*halo.reconstructed_parameters);
            continue;
        }
        const auto contributor=halo.classification.original_contributor_index;
        if (!contributor)
            throw std::invalid_argument("Free Halo has no OriginalInput contributor mapping.");
        halo_input.atom_ids.push_back(halo.classification.atom_id);
        halo_input.support.push_back(unique_support(*contributor));
        for (const auto & entry:halo_input.support.back()) free_sensitive_rows.at(entry.row)=true;
        initial_widths.push_back(std::min((*preclassification.global_target_median)[1],*calibration.bmax4));
    }

    if (!std::all_of(halo_input.observations.begin(),halo_input.observations.end(),
            [](double value){return std::isfinite(value);}))
    {
        out.unavailable_reason="nonfinite-effective-response";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    out.diagnostics.free_support_entries=0;
    for (const auto & support:halo_input.support)
        out.diagnostics.free_support_entries+=support.size();
    out.diagnostics.zero_sensitivity_rows=static_cast<std::size_t>(std::count(
        free_sensitive_rows.begin(),free_sensitive_rows.end(),false));
    out.diagnostics.solver_invoked=true;

    const auto free_atom_ids=halo_input.atom_ids;
    JointProblem halo_problem(std::move(halo_input));
    halo_problem=JointProblemAccess::WithReferenceNormalization(halo_problem,problem);
    std::map<std::string,std::vector<std::size_t>> component_atom_indices;
    for(const auto & component:JointProblemAccess::Get(halo_problem).partition.components)
        component_atom_indices.emplace(component.id,
            std::vector<std::size_t>(component.atoms.begin(),component.atoms.end()));
    joint_component::FixedNeighborSearchPolicy policy;
    policy.eta_upper_bound=calibration.eta_max4;
    policy.collect_diagnostics=collect_solver_diagnostics;
    if(collect_solver_diagnostics)
        policy.search_result_observer=[&out,&free_atom_ids,&component_atom_indices](std::string_view component_id,
            const joint_component::FixedNeighborSearchResult & search) {
            HaloFullABCComponentSearchDiagnostic component;
            component.component_id=std::string(component_id);
            component.reason=search.reason;
            component.search_converged=search.search_converged;
            component.sweep_count=search.sweep_count;
            if(std::isfinite(search.state.objective)) component.objective=search.state.objective;
            if(search.final_sweep)
            {
                component.global_ac_kkt=search.final_sweep->global_ac_kkt;
                component.global_width_gradient_inf_norm=search.final_sweep->global_width_gradient_inf_norm;
                component.eta_change_inf=search.final_sweep->eta_change_inf;
                component.objective_replay_error=search.final_sweep->objective_replay_error;
                component.cache_replay_error=search.final_sweep->cache_replay_error;
            }
            if(search.certification)
            {
                component.endpoint_verified=search.certification->verified();
                component.endpoint_qualified=search.certification->qualified();
                component.endpoint_failure=search.certification->assessment.failure;
                component.trust_failure=search.certification->trust.reason;
            }
            else component.endpoint_failure=search.reason;
            component.local_conditioning=search.local_conditioning;
            if(const auto atoms=component_atom_indices.find(std::string(component_id));
                atoms!=component_atom_indices.end() && atoms->second.size()==static_cast<std::size_t>(search.state.eta.size()) &&
                    search.state.beta.size()==2*search.state.eta.size())
                for(std::size_t k=0;k<atoms->second.size();++k)
                {
                    const auto parent=atoms->second[k];
                    if(parent>=free_atom_ids.size()) continue;
                    component.candidate_profiles.push_back({free_atom_ids[parent],{
                        search.state.beta(2*static_cast<Eigen::Index>(k)),
                        std::exp(search.state.eta(static_cast<Eigen::Index>(k))),
                        search.state.beta(2*static_cast<Eigen::Index>(k)+1)}});
                }
            for(const auto & block:search.blocks)
            {
                component.blocks.push_back({block.sweep,block.block,block.affected_rows,
                    block.status,block.reason,block.local_search_stop_reason,block.profile_evaluations,
                    block.accepted_updates,block.local_objective_before,block.local_objective_after,
                    block.objective_before,block.objective_after,block.global_replay_delta,
                    block.local_global_delta_error,block.local_final_ac_kkt});
                if(out.diagnostics.first_failing_block ||
                    (!block.failure_diagnostics && block.status!="failed" && block.status!="retry-pending")) continue;
                HaloFullABCBlockFailureDiagnostic captured;
                captured.component_index=block.failure_diagnostics ? block.failure_diagnostics->component_index : 0;
                captured.sweep_index=block.failure_diagnostics ? block.failure_diagnostics->sweep_index : block.sweep;
                captured.block_index=block.failure_diagnostics ? block.failure_diagnostics->block_index : block.block;
                captured.affected_rows=block.affected_rows;
                captured.component_id=block.failure_diagnostics && !block.failure_diagnostics->component_id.empty() ?
                    block.failure_diagnostics->component_id : std::string(component_id);
                captured.status=block.status; captured.reason=block.reason;
                captured.local_objective_before=block.local_objective_before;
                captured.local_objective_after=block.local_objective_after;
                captured.global_objective_before=block.objective_before;
                captured.global_objective_after=block.objective_after;
                captured.global_replay_delta=block.global_replay_delta;
                captured.local_global_delta_error=block.local_global_delta_error;
                captured.local_final_ac_kkt=block.local_final_ac_kkt;
                const auto sweep_index=captured.sweep_index ? captured.sweep_index-1 : 0;
                if(sweep_index<search.sweeps.size())
                {
                    captured.global_ac_kkt=search.sweeps[sweep_index].global_ac_kkt;
                    captured.global_width_gradient_inf_norm=search.sweeps[sweep_index].global_width_gradient_inf_norm;
                }
                if(const auto * failure=block.failure_diagnostics ? &*block.failure_diagnostics : nullptr)
                {
                    captured.core_atom_ids=failure->core_atom_ids;
                    captured.local_design_rows=static_cast<std::size_t>(std::max<Eigen::Index>(0,failure->local_design_rows));
                    captured.local_design_columns=static_cast<std::size_t>(std::max<Eigen::Index>(0,failure->local_design_columns));
                    if(failure->local_design_rank && *failure->local_design_rank>=0)
                        captured.local_design_rank=static_cast<std::size_t>(*failure->local_design_rank);
                    if(failure->first_invalid_coefficient && *failure->first_invalid_coefficient>=0)
                        captured.first_invalid_coefficient=static_cast<std::size_t>(*failure->first_invalid_coefficient);
                    captured.active_amplitude_indices.assign(failure->active_amplitude_indices.begin(),
                        failure->active_amplitude_indices.end());
                    captured.local_column_norms=failure->local_column_norms;
                    captured.local_design_singular_values=failure->local_design_singular_values;
                    captured.initial_widths=failure->initial_widths;
                    captured.width_derivative_column_norms=failure->width_derivative_column_norms;
                    captured.projected_width_singular_values=failure->projected_width_singular_values;
                    captured.projected_width_column_norms=failure->projected_width_column_norms;
                    captured.proposed_eta_step=failure->proposed_eta_step;
                    captured.local_design_rank_threshold=failure->local_design_rank_threshold;
                    captured.condition_indicator=failure->condition_indicator;
                    captured.local_response_norm=failure->local_response_norm;
                    captured.observation_scale=failure->parent_observation_scale;
                    captured.primal_feasibility=failure->primal_feasibility;
                    captured.scaled_projected_kkt=failure->scaled_projected_kkt;
                    captured.width_upper_bound=failure->width_upper_bound;
                    captured.width_step_inf=failure->width_step_inf;
                    captured.width_gradient_inf=failure->width_gradient_inf;
                    captured.lm_actual_decrease=failure->lm_actual_decrease;
                    captured.lm_predicted_decrease=failure->lm_predicted_decrease;
                    captured.lm_ratio=failure->lm_ratio;
                    captured.weak_width_step_fraction=failure->weak_width_step_fraction;
                    captured.projected_width_rank_threshold=failure->projected_width_rank_threshold;
                    captured.initial_profile_valid=failure->initial_profile_valid;
                    captured.initial_profile_trusted=failure->initial_profile_trusted;
                    captured.solver_valid=failure->solver_valid;
                    captured.width_bound_active=failure->width_bound_active;
                    captured.initial_profile_reason=failure->initial_profile_reason;
                    captured.local_search_stop_reason=failure->local_search_stop_reason;
                    captured.solver_reason=failure->solver_reason;
                    captured.lm_subproblem_status=failure->lm_subproblem_status;
                    captured.first_invalid_atom_id=failure->first_invalid_atom_id;
                    captured.first_invalid_parameter=failure->first_invalid_parameter;
                    captured.replay_prediction_discrepancy=failure->replay_prediction_discrepancy;
                    captured.replay_gradient_discrepancy=failure->replay_gradient_discrepancy;
                    captured.replay_kkt_discrepancy=failure->replay_kkt_discrepancy;
                }
                out.diagnostics.first_failing_block=std::move(captured);
                break;
            }
            out.diagnostics.component_searches.push_back(std::move(component));
        };
    auto fit=joint_component::FitFixedNeighborComponents(halo_problem,initial_widths,policy);
    out.diagnostics.solver_status=fit.RuntimeConvergence();
    out.diagnostics.solver_qualified=fit.search_completed &&
        out.diagnostics.solver_status==JointCheckStatus::Passed;
    out.diagnostics.objective=fit.objective;
    out.diagnostics.observation_scale=fit.observation_scale;
    out.diagnostics.solver_provenance=fit.solver_provenance;
    out.diagnostics.evidence=fit.evidence;
    const auto & memory=JointProblemAccess::Get(halo_problem).construction_memory;
    if (memory.estimated_peak_bytes)
        out.diagnostics.estimated_peak_memory_bytes=memory.estimated_peak_bytes;
    for (const auto & component:fit.components)
        if (component.local_conditioning_peak_memory_estimate_bytes &&
            (!out.diagnostics.estimated_peak_memory_bytes ||
                *component.local_conditioning_peak_memory_estimate_bytes>*out.diagnostics.estimated_peak_memory_bytes))
            out.diagnostics.estimated_peak_memory_bytes=component.local_conditioning_peak_memory_estimate_bytes;

    if (out.diagnostics.solver_qualified && fit.assembled_state &&
        fit.assembled_state->ac.size()==2*free_atom_ids.size() &&
        fit.assembled_state->b.size()==free_atom_ids.size() &&
        fit.assembled_state->log_b.size()==free_atom_ids.size())
    {
        for (std::size_t index=0;index<free_atom_ids.size();++index)
        {
            const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                return halo.classification.atom_id==free_atom_ids[index];
            });
            if (found==out.halos.end()) throw std::logic_error("FullABC Halo result mapping is incomplete.");
            const std::array<double,3> parameters{fit.assembled_state->ac[2*index],
                fit.assembled_state->b[index],fit.assembled_state->ac[2*index+1]};
            if (!std::isfinite(parameters[0]) || parameters[0]<0 || !std::isfinite(parameters[1]) ||
                !(parameters[1]>0) || fit.assembled_state->log_b[index]>*calibration.eta_max4 ||
                !std::isfinite(parameters[2]))
            {
                out.diagnostics.solver_qualified=false;
                out.diagnostics.solver_status=JointCheckStatus::Failed;
                out.diagnostics.unavailable_reason="invalid-fullabc-halo-endpoint";
                break;
            }
            found->reconstructed_parameters=parameters;
            found->parameter_source=HaloFullABCParameterSource::FullABC;
            found->qualified=true;
            found->qualification_reason="fullabc-endpoint-qualified";
        }
    }
    if (!out.diagnostics.solver_qualified && collect_solver_diagnostics)
        for (const auto & component:out.diagnostics.component_searches)
        {
            if (!component.endpoint_qualified) continue;
            for (const auto & candidate:component.candidate_profiles)
            {
                const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                    return halo.classification.atom_id==candidate.atom_id;
                });
                if (found==out.halos.end() || found->classification.role!=HaloFullABCRole::FreeFullABC)
                    continue;
                found->reconstructed_parameters=candidate.parameters;
                found->parameter_source=HaloFullABCParameterSource::FullABC;
                found->qualified=true;
                found->qualification_reason="independently-qualified-fullabc-component";
            }
        }
    if (out.diagnostics.solver_qualified && std::any_of(out.halos.begin(),out.halos.end(),[](const auto & halo) {
            return halo.classification.role==HaloFullABCRole::FreeFullABC && !halo.reconstructed_parameters;
        }))
    {
        out.diagnostics.solver_qualified=false;
        out.diagnostics.solver_status=JointCheckStatus::Unavailable;
        out.diagnostics.unavailable_reason="assembled-halo-endpoint-unavailable";
    }
    if (out.diagnostics.solver_qualified)
    {
        out.available=true;
        out.diagnostics.reconstruction_qualified=true;
        out.diagnostics.assembly_replay_qualified=true;
        return out;
    }
    out.unavailable_reason=out.diagnostics.unavailable_reason.empty() ?
        "halo-fullabc-solver-unqualified" : out.diagnostics.unavailable_reason;
    out.diagnostics.unavailable_reason=out.unavailable_reason;
    return out;
}

HaloFullABCReconstructionResult EstimateHaloFullABCJointWithRecovery(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result,
    std::string_view joint_run_id,HaloFullABCRecoveryMethod method)
{
    const auto recovery_method_name=[](HaloFullABCRecoveryMethod value) -> std::string {
        switch(value)
        {
        case HaloFullABCRecoveryMethod::FixedWidthACPerHalo: return "fixed-width-ac-per-halo";
        case HaloFullABCRecoveryMethod::FixedWidthACPerComponent: return "fixed-width-ac-per-component";
        case HaloFullABCRecoveryMethod::PostFailureGlobalMedian: return "post-failure-global-median";
        case HaloFullABCRecoveryMethod::IsolateFailedComponent: return "isolate-failed-component";
        case HaloFullABCRecoveryMethod::None: return "none";
        }
        return "unknown";
    };
    auto out=EstimateHaloFullABCJoint(model,problem,result,joint_run_id,true);
    out.diagnostics.recovery_method=recovery_method_name(method);
    if (out.available || method==HaloFullABCRecoveryMethod::None)
        return out;
    if (!out.diagnostics.solver_invoked)
    {
        out.unavailable_reason="recovery-fullabc-solver-not-invoked";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    if (!out.global_target_median)
    {
        out.unavailable_reason="recovery-global-target-median-unavailable";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    if (!out.diagnostics.width_upper_bound)
    {
        out.unavailable_reason="recovery-width-bound-unavailable";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }

    const auto & original=problem.OriginalInput();
    const auto preclassification=BuildHaloFullABCPreclassification(model,problem,result);
    if (!preclassification.available || !preclassification.global_target_median)
    {
        out.unavailable_reason="recovery-preclassification-unavailable:"+preclassification.unavailable_reason;
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    const auto target_points=TargetPoints(result);
    JointProblemInput halo_input;
    halo_input.row_ids=original.row_ids;
    halo_input.observations=original.observations;
    const auto unique_support=[&](std::size_t contributor) {
        std::map<std::size_t,double> by_row;
        for (const auto & entry:original.support.at(contributor))
        {
            const auto [found,inserted]=by_row.emplace(entry.row,entry.squared_distance);
            if (!inserted && found->second!=entry.squared_distance)
                throw std::invalid_argument("Conflicting duplicate Halo support row.");
        }
        std::vector<JointSupport> support;
        support.reserve(by_row.size());
        for (const auto [row,square]:by_row) support.push_back({row,square});
        return support;
    };
    const auto subtract_profile=[&](std::size_t contributor,const std::array<double,3> & parameters) {
        for (const auto & entry:unique_support(contributor))
        {
            const auto basis=joint_component::EvaluateKernel(entry.squared_distance,parameters[1],2.5);
            halo_input.observations.at(entry.row)-=parameters[0]*basis.gaussian+parameters[2]*basis.charge;
        }
    };
    for (const auto target:original.selection_domain->target_indices)
    {
        const auto & point=target_points.at(target);
        for (const auto & entry:unique_support(target))
        {
            const auto basis=joint_component::EvaluateKernel(entry.squared_distance,point.GetWidth(),2.5);
            halo_input.observations.at(entry.row)-=
                point.GetAmplitude()*basis.gaussian+point.GetOffset()*basis.charge;
        }
    }
    std::vector<std::size_t> free_contributors;
    for (const auto & halo:out.halos)
        if (halo.classification.role==HaloFullABCRole::FixedGlobalTargetMedian)
        {
            if (halo.classification.original_contributor_index)
                subtract_profile(*halo.classification.original_contributor_index,*out.global_target_median);
        }
        else
        {
            if (!halo.classification.original_contributor_index) return out;
            halo_input.atom_ids.push_back(halo.classification.atom_id);
            halo_input.support.push_back(unique_support(*halo.classification.original_contributor_index));
            free_contributors.push_back(*halo.classification.original_contributor_index);
        }
    if (!std::all_of(halo_input.observations.begin(),halo_input.observations.end(),
            [](double value){return std::isfinite(value);})) return out;

    const auto qualified_fullabc=[&](std::string_view atom_id) {
        const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
            return halo.classification.atom_id==atom_id &&
                halo.parameter_source==HaloFullABCParameterSource::FullABC && halo.qualified;
        });
        return found!=out.halos.end();
    };

    JointProblem halo_problem(std::move(halo_input));
    halo_problem=JointProblemAccess::WithReferenceNormalization(halo_problem,problem);
    const auto & parent_data=JointProblemAccess::Get(halo_problem);
    const double recovery_width=std::min((*out.global_target_median)[1],*out.diagnostics.width_upper_bound);
    const auto run_fixed_width_ac=[&](const std::vector<std::size_t> & parent_atoms,
        const std::vector<double> & parent_response) {
        HaloFullABCRecoveryDiagnostic diagnostic;
        diagnostic.attempted_halos=parent_atoms.size();
        diagnostic.fixed_width=recovery_width;
        JointProblemInput local_input;
        std::set<std::size_t> parent_rows;
        for (const auto atom:parent_atoms)
        {
            diagnostic.atom_ids.push_back(parent_data.input->atom_ids.at(atom));
            for (const auto & support:parent_data.input->support.at(atom)) parent_rows.insert(support.row);
        }
        if (parent_rows.empty())
        {
            diagnostic.status="unavailable";
            diagnostic.reason="component-has-no-support-rows";
            return std::pair{std::move(diagnostic),std::vector<std::array<double,3>>{}};
        }
        std::map<std::size_t,std::size_t> row_to_local;
        for (const auto row:parent_rows)
        {
            row_to_local.emplace(row,local_input.row_ids.size());
            local_input.row_ids.push_back(parent_data.input->row_ids.at(row));
            local_input.observations.push_back(parent_response.at(row));
        }
        for (const auto atom:parent_atoms)
        {
            local_input.atom_ids.push_back(parent_data.input->atom_ids.at(atom));
            std::vector<JointSupport> support;
            for (const auto & entry:parent_data.input->support.at(atom))
                support.push_back({row_to_local.at(entry.row),entry.squared_distance});
            local_input.support.push_back(std::move(support));
        }
        JointProblem local_problem(std::move(local_input));
        local_problem=JointProblemAccess::WithReferenceNormalization(local_problem,halo_problem);
        const auto & local_data=JointProblemAccess::Get(local_problem);
        const auto atoms=static_cast<Eigen::Index>(parent_atoms.size());
        const Eigen::VectorXd eta=Eigen::VectorXd::Constant(atoms,std::log(recovery_width));
        const auto profile=joint_component::EvaluateProfile(local_data.domain,local_data.y,eta,false,&local_data.context);
        if (profile.valid && profile.certificate.available)
        {
            const auto spectrum=joint_component::ComputeSpectrum(profile.x,local_data.context.rank,
                profile.x.cols(),true);
            diagnostic.ac_rank=static_cast<std::size_t>(std::max<Eigen::Index>(0,spectrum.rank));
            diagnostic.ac_columns=static_cast<std::size_t>(profile.x.cols());
            if (std::isfinite(spectrum.condition)) diagnostic.ac_condition=spectrum.condition;
            diagnostic.full_ac_rank=spectrum.available && !spectrum.threshold_sensitive &&
                spectrum.rank==profile.x.cols();
            const auto replay=joint_component::CheckIndependentReplay(local_data.domain,local_data.y,
                profile,local_data.context);
            diagnostic.independent_replay_passed=replay.passed;
            if (std::isfinite(replay.primary_ac_kkt)) diagnostic.primary_ac_kkt=replay.primary_ac_kkt;
            if (std::isfinite(replay.independent_ac_kkt)) diagnostic.independent_ac_kkt=replay.independent_ac_kkt;
            if (std::isfinite(replay.objective_difference)) diagnostic.objective_replay_error=replay.objective_difference;
            diagnostic.qualified=diagnostic.full_ac_rank && profile.certificate.kkt_passed && replay.passed;
            if (diagnostic.qualified)
            {
                diagnostic.qualified_halos=parent_atoms.size();
                diagnostic.status="qualified";
                diagnostic.reason="full-ac-rank-independent-replay-passed";
                std::vector<std::array<double,3>> parameters;
                parameters.reserve(parent_atoms.size());
                for (Eigen::Index atom=0;atom<atoms;++atom)
                    parameters.push_back({profile.beta(2*atom),recovery_width,profile.beta(2*atom+1)});
                return std::pair{std::move(diagnostic),std::move(parameters)};
            }
            diagnostic.reason=!diagnostic.full_ac_rank ? "ac-design-not-full-rank" :
                !profile.certificate.kkt_passed ? "primary-ac-kkt-failed" : replay.reason;
        }
        else diagnostic.reason=profile.reason.empty() ? "fixed-width-ac-profile-unavailable" : profile.reason;
        diagnostic.status="failed";
        return std::pair{std::move(diagnostic),std::vector<std::array<double,3>>{}};
    };

    out.diagnostics.recovery_method=recovery_method_name(method);
    const auto & components=parent_data.partition.components;
    for (const auto & component:components)
    {
        bool has_failed_halo{};
        for (const auto atom:component.atoms)
            has_failed_halo |= !qualified_fullabc(parent_data.input->atom_ids.at(static_cast<std::size_t>(atom)));
        if (!has_failed_halo) continue;
        HaloFullABCRecoveryDiagnostic diagnostic;
        diagnostic.component_id=component.id;
        for (const auto atom:component.atoms) diagnostic.atom_ids.push_back(
            parent_data.input->atom_ids.at(static_cast<std::size_t>(atom)));
        diagnostic.attempted_halos=component.atoms.size();
        diagnostic.method=out.diagnostics.recovery_method;
        if (method==HaloFullABCRecoveryMethod::FixedWidthACPerComponent)
        {
            std::vector<std::size_t> atoms;
            for (const auto atom:component.atoms) atoms.push_back(static_cast<std::size_t>(atom));
            auto solved=run_fixed_width_ac(atoms,halo_problem.OriginalInput().observations);
            diagnostic=std::move(solved.first);
            diagnostic.component_id=component.id;
            diagnostic.method=out.diagnostics.recovery_method;
            if (diagnostic.qualified)
                for (std::size_t k=0;k<atoms.size();++k)
                {
                    auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==diagnostic.atom_ids[k];
                    });
                    found->reconstructed_parameters=solved.second[k];
                    found->parameter_source=HaloFullABCParameterSource::ObservationACFixedWidthRecovery;
                    found->qualified=true;
                    found->qualification_reason="fixed-width-component-ac-independent-replay";
                }
            else
            {
                HaloFullABCRecoveryDiagnostic fallback;
                fallback.component_id=component.id;
                fallback.method=out.diagnostics.recovery_method+"-median-fallback";
                fallback.status="qualified";
                fallback.reason="component-ac-recovery-failed-then-fixed-at-global-median";
                fallback.attempted_halos=component.atoms.size();
                fallback.qualified_halos=component.atoms.size();
                fallback.qualified=true;
                for (const auto atom:component.atoms)
                {
                    const auto atom_id=parent_data.input->atom_ids.at(static_cast<std::size_t>(atom));
                    fallback.atom_ids.push_back(atom_id);
                    const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==atom_id;
                    });
                    found->reconstructed_parameters=*out.global_target_median;
                    found->parameter_source=HaloFullABCParameterSource::PostFailureGlobalMedianRecovery;
                    found->qualified=true;
                    found->qualification_reason="component-ac-failed-global-median-fallback";
                }
                out.diagnostics.recovery.push_back(std::move(fallback));
            }
        }
        else if (method==HaloFullABCRecoveryMethod::FixedWidthACPerHalo)
        {
            diagnostic.status="qualified";
            diagnostic.reason="all-per-halo-fixed-width-ac-profiles-qualified";
            diagnostic.qualified=true;
            diagnostic.full_ac_rank=true;
            diagnostic.independent_replay_passed=true;
            for (const auto atom:component.atoms)
            {
                const auto atom_index=static_cast<std::size_t>(atom);
                auto response=halo_problem.OriginalInput().observations;
                for (const auto other:component.atoms)
                {
                    if (other==atom) continue;
                    const auto global=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==parent_data.input->atom_ids.at(
                            static_cast<std::size_t>(other));
                    });
                    if (global==out.halos.end()) continue;
                    const auto contributor=global->classification.original_contributor_index;
                    if (!contributor) continue;
                    const auto & parameters=global->parameter_source==HaloFullABCParameterSource::FullABC &&
                        global->reconstructed_parameters ? *global->reconstructed_parameters :
                            *out.global_target_median;
                    for (const auto & entry:unique_support(*contributor))
                    {
                        const auto basis=joint_component::EvaluateKernel(entry.squared_distance,
                            parameters[1],2.5);
                        response.at(entry.row)-=parameters[0]*basis.gaussian+parameters[2]*basis.charge;
                    }
                }
                auto solved=run_fixed_width_ac({atom_index},response);
                diagnostic.qualified &= solved.first.qualified;
                diagnostic.qualified_halos+=solved.first.qualified_halos;
                if (!solved.first.qualified && diagnostic.status=="qualified")
                {diagnostic.status="partial"; diagnostic.reason=solved.first.reason;}
                if (solved.first.ac_rank)
                    diagnostic.ac_rank=diagnostic.ac_rank ?
                        std::min(*diagnostic.ac_rank,*solved.first.ac_rank) : solved.first.ac_rank;
                if (solved.first.ac_columns)
                    diagnostic.ac_columns=diagnostic.ac_columns ?
                        *diagnostic.ac_columns+*solved.first.ac_columns : solved.first.ac_columns;
                if (solved.first.ac_condition)
                    diagnostic.ac_condition=diagnostic.ac_condition ?
                        std::max(*diagnostic.ac_condition,*solved.first.ac_condition) : solved.first.ac_condition;
                if (solved.first.fixed_width) diagnostic.fixed_width=solved.first.fixed_width;
                if (solved.first.primary_ac_kkt)
                    diagnostic.primary_ac_kkt=diagnostic.primary_ac_kkt ?
                        std::max(*diagnostic.primary_ac_kkt,*solved.first.primary_ac_kkt) :
                        solved.first.primary_ac_kkt;
                if (solved.first.independent_ac_kkt)
                    diagnostic.independent_ac_kkt=diagnostic.independent_ac_kkt ?
                        std::max(*diagnostic.independent_ac_kkt,*solved.first.independent_ac_kkt) :
                        solved.first.independent_ac_kkt;
                if (solved.first.objective_replay_error)
                    diagnostic.objective_replay_error=diagnostic.objective_replay_error ?
                        std::max(*diagnostic.objective_replay_error,*solved.first.objective_replay_error) :
                        solved.first.objective_replay_error;
                diagnostic.full_ac_rank &= solved.first.full_ac_rank;
                diagnostic.independent_replay_passed &= solved.first.independent_replay_passed;
                if (solved.first.qualified)
                {
                    const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==parent_data.input->atom_ids.at(atom_index);
                    });
                    found->reconstructed_parameters=solved.second.front();
                    found->parameter_source=HaloFullABCParameterSource::ObservationACFixedWidthRecovery;
                    found->qualified=true;
                    found->qualification_reason="fixed-width-per-halo-ac-independent-replay";
                }
            }
            if (diagnostic.qualified_halos==diagnostic.attempted_halos) diagnostic.status="qualified";
        }
        else if (method==HaloFullABCRecoveryMethod::PostFailureGlobalMedian)
        {
            diagnostic.status="qualified";
            diagnostic.qualified=true;
            diagnostic.reason="failed-component-fixed-at-global-median";
            for (const auto atom:component.atoms)
            {
                const auto atom_id=parent_data.input->atom_ids.at(static_cast<std::size_t>(atom));
                auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                    return halo.classification.atom_id==atom_id;
                });
                if (!qualified_fullabc(atom_id))
                {
                    found->reconstructed_parameters=*out.global_target_median;
                    found->parameter_source=HaloFullABCParameterSource::PostFailureGlobalMedianRecovery;
                    found->qualified=true;
                    found->qualification_reason="failed-component-global-median";
                    ++diagnostic.qualified_halos;
                }
                else ++diagnostic.qualified_halos;
            }
        }
        else if (method==HaloFullABCRecoveryMethod::IsolateFailedComponent)
        {
            std::optional<std::size_t> isolated_atom;
            const auto set_isolated_from_id=[&](const std::string & atom_id) {
                for (const auto atom:component.atoms)
                    if (parent_data.input->atom_ids.at(static_cast<std::size_t>(atom))==atom_id)
                    {isolated_atom=static_cast<std::size_t>(atom); return;}
            };
            if (out.diagnostics.first_failing_block &&
                out.diagnostics.first_failing_block->component_id==component.id)
            {
                if (!out.diagnostics.first_failing_block->first_invalid_atom_id.empty())
                    set_isolated_from_id(out.diagnostics.first_failing_block->first_invalid_atom_id);
                if (!isolated_atom && out.diagnostics.first_failing_block->core_atom_ids.size()==1)
                    set_isolated_from_id(out.diagnostics.first_failing_block->core_atom_ids.front());
            }
            if (!isolated_atom && component.atoms.size()==1)
                isolated_atom=static_cast<std::size_t>(component.atoms.front());
            std::vector<std::size_t> remaining_atoms;
            for (const auto atom:component.atoms)
                if (!isolated_atom || static_cast<std::size_t>(atom)!=*isolated_atom)
                    remaining_atoms.push_back(static_cast<std::size_t>(atom));
            if (!isolated_atom)
            {
                diagnostic.reason="failed-halo-not-localized; failed-component-fixed-at-global-median";
                for (const auto atom:component.atoms)
                {
                    const auto atom_id=parent_data.input->atom_ids.at(static_cast<std::size_t>(atom));
                    const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                        return halo.classification.atom_id==atom_id;
                    });
                    found->reconstructed_parameters=*out.global_target_median;
                    found->parameter_source=HaloFullABCParameterSource::PostFailureGlobalMedianRecovery;
                    found->qualified=true;
                    found->qualification_reason="failed-component-global-median-not-localized";
                    ++diagnostic.qualified_halos;
                }
            }
            else
            {
                const auto atom_id=parent_data.input->atom_ids.at(*isolated_atom);
                const auto isolated=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                    return halo.classification.atom_id==atom_id;
                });
                isolated->reconstructed_parameters=*out.global_target_median;
                isolated->parameter_source=HaloFullABCParameterSource::PostFailureGlobalMedianRecovery;
                isolated->qualified=true;
                isolated->qualification_reason="failed-halo-isolated-at-global-median";
                ++diagnostic.qualified_halos;
                diagnostic.reason=remaining_atoms.empty() ? "failed-halo-isolated; no-other-halo-in-component" :
                    "failed-halo-isolated; remaining-component-refit";
                if (!remaining_atoms.empty())
                {
                    diagnostic.fullabc_refit_attempted=true;
                    JointProblemInput refit_input;
                    std::map<std::size_t,std::size_t> row_to_local;
                    for (const auto row:component.rows)
                    {
                        row_to_local.emplace(static_cast<std::size_t>(row),refit_input.row_ids.size());
                        refit_input.row_ids.push_back(parent_data.input->row_ids.at(static_cast<std::size_t>(row)));
                        refit_input.observations.push_back(parent_data.input->observations.at(
                            static_cast<std::size_t>(row)));
                    }
                    const auto isolated_contributor=isolated->classification.original_contributor_index;
                    if (isolated_contributor)
                        for (const auto & support:unique_support(*isolated_contributor))
                        {
                            const auto basis=joint_component::EvaluateKernel(support.squared_distance,
                                (*out.global_target_median)[1],2.5);
                            const auto local_row=row_to_local.find(support.row);
                            if (local_row!=row_to_local.end())
                                refit_input.observations.at(local_row->second)-=
                                    (*out.global_target_median)[0]*basis.gaussian+
                                    (*out.global_target_median)[2]*basis.charge;
                        }
                    std::vector<double> initial_widths;
                    for (const auto atom:remaining_atoms)
                    {
                        refit_input.atom_ids.push_back(parent_data.input->atom_ids.at(atom));
                        std::vector<JointSupport> support;
                        for (const auto & entry:parent_data.input->support.at(atom))
                            support.push_back({row_to_local.at(entry.row),entry.squared_distance});
                        refit_input.support.push_back(std::move(support));
                        initial_widths.push_back(recovery_width);
                    }
                    JointProblem refit_problem(std::move(refit_input));
                    refit_problem=JointProblemAccess::WithReferenceNormalization(refit_problem,halo_problem);
                    joint_component::FixedNeighborSearchPolicy refit_policy;
                    refit_policy.eta_upper_bound=std::log(*out.diagnostics.width_upper_bound);
                    const auto refit=joint_component::FitFixedNeighborComponents(
                        refit_problem,initial_widths,refit_policy);
                    const auto has_pass=[&](std::string_view name) {
                        return std::any_of(refit.evidence.begin(),refit.evidence.end(),[&](const auto & check) {
                            return check.name==name && check.scope==JointEvidenceScope::AssembledGlobal &&
                                check.status==JointCheckStatus::Passed;
                        });
                    };
                    diagnostic.fullabc_refit_qualified=refit.search_completed &&
                        refit.RuntimeConvergence()==JointCheckStatus::Passed && refit.assembled_state &&
                        refit.assembled_state->ac.size()==2*remaining_atoms.size() &&
                        refit.assembled_state->b.size()==remaining_atoms.size() &&
                        has_pass("independent-replay") && has_pass("width-stationarity") &&
                        has_pass("width-bound-feasibility") && has_pass("endpoint-verification") &&
                        has_pass("numerical-qualification") && has_pass("assembly-mapping") &&
                        has_pass("assembled-consistency");
                    for (std::size_t k=0;k<remaining_atoms.size();++k)
                    {
                        const auto remaining_id=parent_data.input->atom_ids.at(remaining_atoms[k]);
                        const auto found=std::find_if(out.halos.begin(),out.halos.end(),[&](const auto & halo) {
                            return halo.classification.atom_id==remaining_id;
                        });
                        if (diagnostic.fullabc_refit_qualified)
                        {
                            const auto & state=*refit.assembled_state;
                            const std::array<double,3> parameters{state.ac[2*k],state.b[k],state.ac[2*k+1]};
                            found->reconstructed_parameters=parameters;
                            found->parameter_source=HaloFullABCParameterSource::FullABC;
                            found->qualified=true;
                            found->qualification_reason="fullabc-refit-after-failed-halo-isolation";
                        }
                        else
                        {
                            found->reconstructed_parameters=*out.global_target_median;
                            found->parameter_source=HaloFullABCParameterSource::PostFailureGlobalMedianRecovery;
                            found->qualified=true;
                            found->qualification_reason="isolated-refit-failed-global-median-fallback";
                        }
                        ++diagnostic.qualified_halos;
                    }
                }
            }
            diagnostic.status="qualified";
            diagnostic.qualified=diagnostic.qualified_halos==diagnostic.attempted_halos;
        }
        out.diagnostics.recovery.push_back(std::move(diagnostic));
    }
    const auto & assembled_response=halo_problem.OriginalInput().observations;
    std::vector<double> primary_prediction(assembled_response.size(),0.0),
        independent_prediction(assembled_response.size(),0.0),prediction_compensation(assembled_response.size(),0.0);
    bool complete_profiles=true;
    const double sqrt_two=std::sqrt(2.0);
    for (const auto & halo:out.halos)
        if (halo.classification.role==HaloFullABCRole::FreeFullABC)
        {
            if (!halo.reconstructed_parameters)
            {
                complete_profiles=false;
                continue;
            }
            const auto & parameters=*halo.reconstructed_parameters;
            const auto contributor=halo.classification.original_contributor_index;
            if (!contributor)
            {
                complete_profiles=false;
                continue;
            }
            for (const auto & support:unique_support(*contributor))
            {
                const auto primary_basis=joint_component::EvaluateKernel(support.squared_distance,
                    parameters[1],2.5);
                primary_prediction.at(support.row)+=parameters[0]*primary_basis.gaussian+
                    parameters[2]*primary_basis.charge;
                const double radius=std::sqrt(support.squared_distance);
                const double gaussian=std::pow(2.0*std::acos(-1.0)*parameters[1]*parameters[1],-1.5)*
                    std::exp(-support.squared_distance/(2.0*parameters[1]*parameters[1]));
                const double charge=radius<1e-5 ? std::sqrt(2.0/std::acos(-1.0))/parameters[1] :
                    std::erf(radius/parameters[1]/sqrt_two)/radius;
                const auto row=support.row;
                const double value=parameters[0]*gaussian+parameters[2]*charge;
                const double sum=independent_prediction[row],next=sum+value;
                prediction_compensation[row]+=std::abs(sum)>=std::abs(value) ?
                    (sum-next)+value : (value-next)+sum;
                independent_prediction[row]=next;
            }
        }
    for (std::size_t row=0;row<independent_prediction.size();++row)
        independent_prediction[row]+=prediction_compensation[row];
    double maximum_prediction_error{};
    bool prediction_replay_passed=true;
    double primary_rss{},primary_rss_compensation{},independent_rss{},independent_rss_compensation{};
    const auto add_square=[](double residual,double & sum,double & compensation) {
        const double square=residual*residual;
        const double adjusted=square-compensation,next=sum+adjusted;
        compensation=(next-sum)-adjusted;
        sum=next;
    };
    for (std::size_t row=0;row<assembled_response.size();++row)
    {
        const double difference=std::abs(primary_prediction[row]-independent_prediction[row]);
        maximum_prediction_error=std::max(maximum_prediction_error,difference);
        prediction_replay_passed &= difference<=2e-12+2e-13*std::abs(primary_prediction[row]);
        add_square(primary_prediction[row]-assembled_response[row],primary_rss,primary_rss_compensation);
        add_square(independent_prediction[row]-assembled_response[row],independent_rss,independent_rss_compensation);
    }
    const double scale=parent_data.context.scale;
    const double primary_objective=.5*primary_rss/(scale*scale);
    const double independent_objective=.5*independent_rss/(scale*scale);
    const double objective_error=std::abs(primary_objective-independent_objective);
    out.diagnostics.assembled_prediction_replay_error=maximum_prediction_error;
    out.diagnostics.assembled_objective_replay_error=objective_error;
    out.diagnostics.original_domain_objective=independent_objective;
    out.diagnostics.original_domain_normalization_scale=scale;
    out.diagnostics.assembly_replay_qualified=complete_profiles && prediction_replay_passed &&
        objective_error<=1e-12+2e-12*std::abs(independent_objective);
    for (const auto & halo:out.halos)
        if (halo.classification.role==HaloFullABCRole::FreeFullABC)
        {
            if (!halo.reconstructed_parameters || !halo.qualified) out.diagnostics.assembly_replay_qualified=false;
            if (halo.parameter_source==HaloFullABCParameterSource::FullABC)
            {
                const auto component=std::find_if(out.diagnostics.component_searches.begin(),
                    out.diagnostics.component_searches.end(),[&](const auto & search) {
                        return search.endpoint_qualified && std::any_of(search.candidate_profiles.begin(),
                            search.candidate_profiles.end(),[&](const auto & candidate) {
                                return candidate.atom_id==halo.classification.atom_id;
                            });
                    });
                const bool original_component=component!=out.diagnostics.component_searches.end();
                const bool isolated_refit=std::any_of(out.diagnostics.recovery.begin(),
                    out.diagnostics.recovery.end(),[&](const auto & recovery) {
                        return recovery.fullabc_refit_qualified &&
                            std::find(recovery.atom_ids.begin(),recovery.atom_ids.end(),
                                halo.classification.atom_id)!=recovery.atom_ids.end();
                    });
                out.diagnostics.assembly_replay_qualified &= original_component || isolated_refit;
            }
            else if (halo.parameter_source==HaloFullABCParameterSource::ObservationACFixedWidthRecovery)
            {
                const auto recovery=std::find_if(out.diagnostics.recovery.begin(),out.diagnostics.recovery.end(),
                    [&](const auto & item) {
                        return std::find(item.atom_ids.begin(),item.atom_ids.end(),halo.classification.atom_id)!=
                            item.atom_ids.end() && item.method.starts_with("fixed-width-ac-");
                    });
                out.diagnostics.assembly_replay_qualified &= recovery!=out.diagnostics.recovery.end() &&
                    recovery->status=="qualified";
            }
            else if (halo.parameter_source!=HaloFullABCParameterSource::PostFailureGlobalMedianRecovery)
                out.diagnostics.assembly_replay_qualified=false;
        }
    out.diagnostics.reconstruction_qualified=out.diagnostics.assembly_replay_qualified;
    if (!out.diagnostics.reconstruction_qualified)
    {
        out.unavailable_reason="no-qualified-recovery-endpoint";
        out.diagnostics.unavailable_reason=out.unavailable_reason;
        return out;
    }
    out.available=true;
    out.unavailable_reason.clear();
    out.diagnostics.unavailable_reason.clear();
    return out;
}

HaloMedianBuildResult BuildHaloMedianReconstruction(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result,
    std::string_view joint_run_id)
{
    const auto & input=problem.OriginalInput();
    if (!input.selection_domain || input.selection_domain->target_indices.empty())
        throw std::invalid_argument("halo-median-no-selection");
    HaloMedianBuildResult out;
    out.target_donors=QualifiedJointTargetDonors(model,problem,result);
    if(out.target_donors.empty()) throw std::invalid_argument("halo-median-no-valid-targets");
    out.snapshot.joint_run_id=std::string(joint_run_id);
    out.snapshot.neighbor_count=GeometryLocalTargetMedianNeighborCount;
    out.snapshot.donor_policy=std::string(GeometryLocalTargetMedianDonorPolicy);
    std::vector<const HaloMedianTargetDonor *> all_donors;
    all_donors.reserve(out.target_donors.size());
    for(const auto & donor:out.target_donors) all_donors.push_back(&donor);
    out.snapshot.valid_target_count=all_donors.size();
    out.snapshot.global_median=Median(all_donors);
    out.audit.valid_target_count=all_donors.size();
    out.audit.global_median=out.snapshot.global_median;
    const GeometryLocalTargetDonorIndex donor_index(out.target_donors);

    const auto & targets=input.selection_domain->target_indices;
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if (std::binary_search(targets.begin(),targets.end(),atom)) continue;
        const auto * object=model.FindAtomPtr(std::stoi(input.atom_ids[atom]));
        if (!object) throw std::invalid_argument("halo-median-unknown-halo");
        const auto donors=donor_index.Query(object->GetPosition());
        if (donors.empty()) throw std::invalid_argument("halo-median-no-valid-targets");
        const auto median=GeometryLocalTargetMedian(donors);
        JointHaloMedianReconstruction halo;
        halo.halo_contributor_index=atom;
        halo.atom_id=input.atom_ids[atom];
        halo.amplitude=median[0];
        halo.width=median[1];
        halo.charge=median[2];
        halo.source=HaloMedianSource::GeometryLocalTargetMedian;
        halo.donor_count=donors.size();
        halo.joint_run_id=std::string(joint_run_id);
        halo.contract=std::string(HaloMedianReconstructionContractVersion);
        for (const auto * donor:donors) halo.donor_atom_ids.push_back(donor->atom_id);
        out.snapshot.halos.push_back(std::move(halo));

        std::set<std::size_t> rows;
        for (const auto & support:input.support[atom]) rows.insert(support.row);
        out.audit.halos.push_back({input.atom_ids[atom],donors.size(),rows.size(),median,
            Dispersion(donors,median)});
    }
    out.audit.geometry_local_halos=out.snapshot.halos.size();
    return out;
}
}
