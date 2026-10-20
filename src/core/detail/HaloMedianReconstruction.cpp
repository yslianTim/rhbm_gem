#include "HaloMedianReconstruction.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include "joint_component/Problem.hpp"
#include "joint_component/WidthCalibration.hpp"
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

    std::unordered_set<int> target_ids;
    for(const auto target:input.selection_domain->target_indices)
        target_ids.insert(std::stoi(input.atom_ids.at(target)));
    std::unordered_map<int,std::size_t> contributor_indices;
    contributor_indices.reserve(input.atom_ids.size());
    for(std::size_t index=0;index<input.atom_ids.size();++index)
        contributor_indices.emplace(std::stoi(input.atom_ids[index]),index);
    std::vector<const AtomObject *> candidates;
    for(const auto & atom:model.GetAtomList())
        if(atom->GetElement()!=Element::HYDROGEN && !target_ids.contains(atom->GetSerialID()))
            candidates.push_back(atom.get());
    std::sort(candidates.begin(),candidates.end(),[](const auto * left,const auto * right) {
        return left->GetSerialID()<right->GetSerialID();
    });
    out.halos.reserve(candidates.size());
    for(const auto * atom:candidates)
    {
        HaloFullABCClassificationRecord record;
        record.atom_id=std::to_string(atom->GetSerialID());
        if(const auto found=contributor_indices.find(atom->GetSerialID());found!=contributor_indices.end())
        {
            record.original_contributor_index=found->second;
            record.original_support_rows=DistinctOriginalSupportRows(input.support.at(found->second));
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
    std::string_view joint_run_id)
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
        out.halos.push_back({halo,halo.role==HaloFullABCRole::FixedGlobalTargetMedian ?
            halo.parameter_profile : std::nullopt});

    if (out.halos.empty())
    {
        out.available=true;
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
    joint_component::FixedNeighborSearchPolicy policy;
    policy.eta_upper_bound=calibration.eta_max4;
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
        return out;
    }
    out.unavailable_reason=out.diagnostics.unavailable_reason.empty() ?
        "halo-fullabc-solver-unqualified" : out.diagnostics.unavailable_reason;
    out.diagnostics.unavailable_reason=out.unavailable_reason;
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
