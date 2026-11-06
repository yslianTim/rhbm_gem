#include "PostFitPeeling.hpp"
#include "HaloMedianReconstruction.hpp"
#include "Numerics.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <cmath>
#include <set>

namespace rhbm_gem::core::detail {
namespace {
std::vector<bool> JointOutputMask(const JointProblemInput & input,
    std::span<const std::size_t> requested)
{
    std::vector<bool> mask(input.atom_ids.size(),false);
    for (auto atom:requested) mask.at(atom)=true;
    return mask;
}

bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role && a.fixed_from_first_stage==b.fixed_from_first_stage;
}

double SquaredDistance(const std::array<double, 3> & left, const std::array<double, 3> & right)
{
    double square{};
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const double delta = left[axis] - right[axis];
        square += delta * delta;
    }
    return square;
}

bool QualifiedTargetPoint(const JointAnalysisResult & result, std::size_t atom)
{
    if (result.runtime_convergence != JointCheckStatus::Passed ||
        result.target_runtime_convergence != JointCheckStatus::Passed || !result.target_evidence)
        return false;
    const auto found = std::find_if(result.target_evidence->point_atoms.begin(),
        result.target_evidence->point_atoms.end(),
        [atom](const auto & point) { return point.atom == atom; });
    return found != result.target_evidence->point_atoms.end() &&
        found->component_convergence == JointCheckStatus::Passed &&
        found->parameter_validity == JointCheckStatus::Passed;
}


}

std::vector<PostJointGridObservation> CollectPostJointGridObservations(
    const MapObject & map, const JointProblemInput & input, std::size_t target,
    const std::array<double,3> & target_position)
{
    auto support = input.support.at(target);
    std::sort(support.begin(), support.end(), [](const auto & a, const auto & b) { return a.row < b.row; });
    std::vector<PostJointGridObservation> samples;
    std::optional<std::size_t> previous;
    for (const auto & entry : support)
    {
        if (previous == entry.row) throw std::invalid_argument("Duplicate Original support row.");
        previous = entry.row;
        if (entry.row >= input.row_ids.size() || entry.row >= input.observations.size() ||
            !std::isfinite(entry.squared_distance) || entry.squared_distance < 0.0 || entry.squared_distance > 6.25)
            throw std::invalid_argument("Invalid Original support identity.");
        const auto grid = std::stoull(input.row_ids[entry.row]);
        if (std::to_string(grid) != input.row_ids[entry.row] || grid >= map.GetMapValueArraySize())
            throw std::invalid_argument("Invalid grid identity.");
        const auto position = map.GetGridPosition(grid);
        const double square = SquaredDistance(position, target_position);
        if (!std::isfinite(square) || std::abs(square - entry.squared_distance) > 1.0e-12 * std::max(1.0,square))
            throw std::invalid_argument("Original support geometry mismatch.");
        if (entry.squared_distance > PostJointCandidateRadius * PostJointCandidateRadius) continue;
        const double response = input.observations[entry.row];
        if (!std::isfinite(response) || !std::isfinite(map.GetMapValue(grid)))
            throw std::invalid_argument("missing-original-observation");
        if (response != map.GetMapValue(grid))
            throw std::invalid_argument("Original observation differs from current map.");
        samples.push_back({entry.row,grid,entry.squared_distance,response,position});
    }
    return samples;
}

JointPeelingContributorUniverse BuildJointPeelingContributorUniverse(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult & result)
{
    JointPeelingContributorUniverse universe;
    if (result.halo_median_reconstruction)
        universe.joint_run_id=result.halo_median_reconstruction->joint_run_id;
    const auto & input=problem.OriginalInput();
    // Every voxel in the 2 A candidate set lies in a Target's Original 2.5 A
    // support. The Original-domain builder already includes every non-H atom
    // contributing to those Original observations.
    universe.contributors.reserve(input.atom_ids.size());
    for (std::size_t index=0;index<input.atom_ids.size();++index)
    {
        const auto * object=model.FindAtomPtr(std::stoi(input.atom_ids[index]));
        JointPeelingContributor contributor;
        contributor.atom_id=object->GetSerialID();
        contributor.position=object->GetPosition();
        contributor.contributor_index=index;
        contributor.reconstruction_eligible=object->GetElement()!=Element::HYDROGEN;
        contributor.target=input.selection_domain && std::binary_search(input.selection_domain->target_indices.begin(),
            input.selection_domain->target_indices.end(),index);
        if (universe.joint_run_id.empty() && contributor.target)
        {
            const auto local=AtomLocalPotentialView::For(*object);
            if (local.IsAvailable()) universe.joint_run_id=local.GetStageEstimate(FittingStage::Second).source.run_id;
        }
        universe.contributors.push_back(std::move(contributor));
    }
    std::sort(universe.contributors.begin(),universe.contributors.end(),
        [](const auto & left,const auto & right){return left.atom_id<right.atom_id;});
    return universe;
}

LocalPotentialSampleList PostJointLocalSamples::BuildGaussianSamples() const
{
    LocalPotentialSampleList result;
    result.reserve(samples.size());
    for (const auto & sample : samples)
        if (sample.gaussian_response)
            result.push_back({*sample.gaussian_response,{sample.distance,sample.position,true}});
    return result;
}

JointPostprocessingContext::JointPostprocessingContext(
    const MapObject & geometry,const ModelObject & model,const JointProblem & problem,
    const JointAnalysisResult & snapshot,JointPostprocessingMetrics * metrics,
    const HaloMedianBuildResult * halo_medians,JointPeelingContributorUniverse contributors)
    : m_geometry(geometry), m_model(model), m_problem(problem), m_snapshot(snapshot), m_metrics(metrics),
      m_halo_medians(halo_medians),m_contributors(std::move(contributors))
{
    const auto & input = problem.OriginalInput();
    if (!input.selection_domain ||
        input.selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy)
        throw std::invalid_argument("Post-fit peeling requires a target-only Pure Exclusion problem.");
    if (snapshot.parameterization_contract!=JointPureHaloSupportExclusionParameterizationContract)
        throw std::invalid_argument("Post-fit peeling requires a target-only Pure Exclusion snapshot.");
    if (!IsQualifiedJointAnalysisResult(snapshot))
        throw std::invalid_argument("Post-fit peeling requires a qualified Joint endpoint.");
    if (snapshot.atom_ids != input.atom_ids)
        throw std::invalid_argument("Joint peeling snapshot identity mismatch.");
    if (snapshot.row_ids != input.row_ids)
        throw std::invalid_argument("Joint peeling snapshot identity mismatch.");
    if (!snapshot.selection_domain ||
        snapshot.selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy ||
        snapshot.available_row_mask.size()!=input.row_ids.size() ||
        snapshot.fitting_row_ids!=problem.Input().row_ids ||
        snapshot.fitting_to_original_contributors!=problem.FittingToOriginalContributors() ||
        snapshot.fitting_to_original_rows!=problem.FittingToOriginalRows())
        throw std::invalid_argument("Joint peeling fitting/reconstruction identity mismatch.");
    for (std::size_t row=0;row<snapshot.fitting_row_ids.size();++row)
    {
        const auto original_row=snapshot.fitting_to_original_rows[row];
        if (original_row>=input.row_ids.size() || snapshot.fitting_row_ids[row]!=input.row_ids[original_row])
            throw std::invalid_argument("Joint peeling fitting/reconstruction identity mismatch.");
    }
    std::size_t indexed_contributors{};
    int previous_contributor_id{};
    bool first_contributor=true;
    for (const auto & contributor:m_contributors.contributors)
    {
        if (contributor.atom_id<=0 || (!first_contributor && contributor.atom_id<=previous_contributor_id) ||
            (contributor.reconstruction_eligible && !std::all_of(contributor.position.begin(),
                contributor.position.end(),[](double value){return std::isfinite(value);})))
            throw std::invalid_argument("Invalid Joint Peeling contributor universe.");
        first_contributor=false;
        previous_contributor_id=contributor.atom_id;
        const auto index=contributor.contributor_index;
        if (index>=input.atom_ids.size() || std::stoi(input.atom_ids[index])!=contributor.atom_id)
            throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
        const bool target=std::binary_search(input.selection_domain->target_indices.begin(),
            input.selection_domain->target_indices.end(),index);
        if (contributor.target!=target || !contributor.reconstruction_eligible)
            throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
        ++indexed_contributors;
    }
    if (indexed_contributors!=input.atom_ids.size())
        throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
    std::set<std::size_t> grid_indices;
    for (std::size_t row=0;row<input.row_ids.size();++row)
    {
        const auto & row_id=input.row_ids[row];
        const auto grid = std::stoull(row_id);
        if (grid >= geometry.GetMapValueArraySize() || row_id != std::to_string(grid) ||
            !grid_indices.insert(grid).second)
            throw std::invalid_argument("Joint original-domain row mapping mismatch.");
        if (!std::isfinite(geometry.GetMapValue(grid)) || input.observations.at(row)!=geometry.GetMapValue(grid))
            throw std::invalid_argument("Original observations differ from current map.");
    }
    if (m_metrics)
    {
        m_metrics->estimated_context_lookup_bytes = m_contributors.contributors.size()*sizeof(JointPeelingContributor);
    }
}

std::vector<std::pair<int, PostFitPeelingResult>> JointPostprocessingContext::BuildPostFitPeelingBatch(
    std::span<const std::size_t> outputs) const
{
    const auto & input = m_problem.OriginalInput();
    const auto & result = m_snapshot;
    auto requested = JointOutputMask(input, outputs);

    // Internal representatives reconstruct the observation prediction even when
    // a nonunique halo has no public atom-wise point estimate.
    std::vector<std::optional<GaussianModel3D>> represented(input.atom_ids.size());
    if(!result.selection_domain ||
        result.selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy)
        throw std::invalid_argument("Invalid Pure Exclusion peeling contract.");
    const bool pure_has_halos=result.atom_ids.size()>result.selection_domain->target_indices.size();
    if (pure_has_halos && !result.halo_median_reconstruction)
        throw std::invalid_argument("Invalid Pure Exclusion peeling contract.");
    if (!pure_has_halos && result.halo_median_reconstruction)
        throw std::invalid_argument("Pure Exclusion Halo medians exist without original-domain Halos.");
    for(std::size_t atom=0;atom<requested.size();++atom)
        if(!std::binary_search(result.selection_domain->target_indices.begin(),
            result.selection_domain->target_indices.end(),atom)) requested[atom]=false;
    for (const auto & c : result.components) if (c.state)
    {
        const auto & full = c.layout ? c.layout->full_atoms : c.atoms;
        for (std::size_t k = 0; k < full.size(); ++k)
        {
            if(full[k]>=represented.size() || represented[full[k]])
                throw std::invalid_argument("Invalid target-only component mapping in peeling snapshot.");
            represented[full[k]] = GaussianModel3D{c.state->ac[2*k], c.state->tau[k], c.state->ac[2*k+1]};
        }
    }
    if(pure_has_halos)
    {
        const auto & snapshot=*result.halo_median_reconstruction;
        if (snapshot.contract!=HaloMedianReconstructionContractVersion || snapshot.joint_run_id.empty() ||
            snapshot.valid_target_count==0 || snapshot.neighbor_count!=GeometryLocalTargetMedianNeighborCount ||
            snapshot.donor_policy!=GeometryLocalTargetMedianDonorPolicy)
            throw std::invalid_argument("Invalid Pure Exclusion Halo median snapshot.");
        const auto valid_median=[](double amplitude,double width,double charge) {
            return std::isfinite(amplitude) && amplitude>=0 && std::isfinite(width) && width>0 &&
                std::isfinite(charge);
        };
        if (!valid_median(snapshot.global_median[0],snapshot.global_median[1],snapshot.global_median[2]))
            throw std::invalid_argument("Invalid Pure Exclusion global Halo median.");
        std::size_t previous{};
        bool first=true;
        for (const auto & halo:snapshot.halos)
        {
            if (halo.halo_contributor_index>=represented.size() || (!first && halo.halo_contributor_index<=previous) ||
                result.atom_ids[halo.halo_contributor_index]!=halo.atom_id ||
                std::binary_search(result.selection_domain->target_indices.begin(),
                    result.selection_domain->target_indices.end(),halo.halo_contributor_index) ||
                halo.contract!=HaloMedianReconstructionContractVersion || halo.joint_run_id!=snapshot.joint_run_id ||
                halo.source!=HaloMedianSource::GeometryLocalTargetMedian || halo.donor_count==0 ||
                halo.donor_count>snapshot.neighbor_count || halo.donor_atom_ids.size()!=halo.donor_count ||
                !valid_median(halo.amplitude,halo.width,halo.charge))
                throw std::invalid_argument("Invalid Pure Exclusion Halo median reconstruction.");
            std::set<std::string> donor_ids;
            for (const auto & donor_id:halo.donor_atom_ids)
                if (donor_id.empty() || !donor_ids.insert(donor_id).second)
                    throw std::invalid_argument("Invalid Pure Exclusion Halo donor identity.");
            previous=halo.halo_contributor_index;
            first=false;
            represented[halo.halo_contributor_index]=GaussianModel3D{halo.amplitude,halo.width,halo.charge};
        }

    }

    std::string expected_joint_run_id;
    if (result.target_evidence)
    {
        const auto target_index=result.selection_domain->target_indices.front();
        const auto * target=m_model.FindAtomPtr(std::stoi(result.atom_ids[target_index]));
        const auto stage=AtomLocalPotentialView::For(*target).GetStageEstimate(FittingStage::Second);
        if (stage.source.method!=EstimateMethod::JointComponents || stage.source.role!=FittingRole::Target ||
            stage.source.run_id.empty())
            throw std::invalid_argument("Invalid Pure Exclusion Joint run provenance.");
        expected_joint_run_id=stage.source.run_id;
        if (result.halo_median_reconstruction &&
            result.halo_median_reconstruction->joint_run_id!=expected_joint_run_id)
            throw std::invalid_argument("Invalid Pure Exclusion Halo median run provenance.");
        if (!m_contributors.joint_run_id.empty() && m_contributors.joint_run_id!=expected_joint_run_id)
            throw std::invalid_argument("Invalid Joint Peeling contributor universe run provenance.");
    }

    std::optional<HaloMedianBuildResult> rebuilt_halo_medians;
    const HaloMedianBuildResult * halo_medians=m_halo_medians;
    if (result.halo_median_reconstruction && (!halo_medians || !halo_medians->snapshot) && result.target_evidence)
    {
        try
        {
            rebuilt_halo_medians=BuildHaloMedianReconstruction(
                m_model,m_problem,result,expected_joint_run_id);
            halo_medians=&*rebuilt_halo_medians;
        }
        catch (const std::invalid_argument & error)
        {
            if (std::string_view(error.what())!="halo-median-no-valid-targets") throw;
            if (result.halo_median_reconstruction)
                throw std::invalid_argument("Saved Pure Exclusion Halo medians have no qualified donors.");
        }
    }
    if (result.halo_median_reconstruction && halo_medians && halo_medians->snapshot)
    {
        const auto & saved=*result.halo_median_reconstruction;
        const auto & rebuilt=*halo_medians->snapshot;
        if (saved.contract!=HaloMedianReconstructionContractVersion || saved.joint_run_id.empty() ||
            saved.joint_run_id!=rebuilt.joint_run_id || saved.valid_target_count!=rebuilt.valid_target_count ||
            saved.global_median!=rebuilt.global_median || saved.neighbor_count!=rebuilt.neighbor_count ||
            saved.donor_policy!=rebuilt.donor_policy || saved.halos!=rebuilt.halos)
            throw std::invalid_argument("Invalid Pure Exclusion Halo median run provenance.");
    }
    // Evaluate every recorded support contribution once. Sparse row incidence
    // permits direct neighbor sums without cancelling the target's own signal
    // against a shared total prediction (especially for very narrow widths).
    struct Contribution { std::size_t atom; std::optional<double> value; };
    std::vector<std::vector<Contribution>> rows(input.row_ids.size());
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        const bool target=std::binary_search(result.selection_domain->target_indices.begin(),
            result.selection_domain->target_indices.end(),atom);
        if (target && !QualifiedTargetPoint(result,atom)) represented[atom].reset();
        const auto & point=represented[atom];
        if (point && (!std::isfinite(point->GetAmplitude()) || point->GetAmplitude()<0 ||
            !std::isfinite(point->GetOffset()) || !std::isfinite(point->GetWidth()) || point->GetWidth()<=0))
            throw std::invalid_argument("Invalid reconstruction contributor parameters.");
        std::set<std::size_t> support_rows;
        for (const auto & support:input.support[atom])
        {
            if (!support_rows.insert(support.row).second)
                throw std::invalid_argument("Duplicate Original contributor support row.");
            std::optional<double> value;
            if (point)
            {
                const auto basis=joint_component::EvaluateKernel(support.squared_distance,point->GetWidth(),2.5);
                const double prediction=point->GetAmplitude()*basis.gaussian+point->GetOffset()*basis.charge;
                if (std::isfinite(prediction)) value=prediction;
            }
            rows.at(support.row).push_back({atom,value});
        }
    }
    std::vector<std::pair<int, PostFitPeelingResult>> results;
    for (std::size_t atom=0;atom<input.atom_ids.size();++atom)
    {
        if (!requested[atom]) continue;
        const auto id=std::stoi(input.atom_ids[atom]);
        const auto * object=m_model.FindAtomPtr(id);
        const auto view=AtomLocalPotentialView::For(*object);
        const auto & estimate=view.GetStageEstimate(FittingStage::Second);
        auto & output=results.emplace_back(id,PostFitPeelingResult{}).second;
        output.source=estimate.source;
        std::set<std::size_t> neighbors;
        for (const auto & observation:CollectPostJointGridObservations(m_geometry,input,atom,object->GetPosition()))
        {
            PeelingSampleEstimate sample;
            sample.observation=observation;
            if (!estimate.point) sample.reason=estimate.reason;
            else if (!represented[atom] || estimate.convergence!=JointCheckStatus::Passed)
                sample.reason="joint-endpoint-not-qualified";
            else if (estimate.point->ToVector()!=represented[atom]->ToVector())
                throw std::invalid_argument("Peeling Target endpoint differs from Joint snapshot.");
            long double prediction{};
            for (const auto & contributor:rows[observation.original_row_index])
            {
                if (contributor.atom==atom) continue;
                neighbors.insert(contributor.atom);
                if (!contributor.value) sample.reason="missing-contributor-state";
                else prediction+=*contributor.value;
            }
            if (sample.reason.empty())
            {
                const double response=static_cast<double>(static_cast<long double>(observation.raw_response)-prediction);
                if (std::isfinite(response)) sample.response=response;
                else sample.reason="nonfinite-peeled-response";
            }
            if (m_metrics)
            {
                ++m_metrics->grid_peeling_sample_count;
                if (sample.response) ++m_metrics->available_peeling_sample_count;
                else
                {
                    ++m_metrics->unavailable_peeling_sample_count;
                    if (sample.reason=="missing-contributor-state") ++m_metrics->missing_contributor_samples;
                    else ++m_metrics->other_peeling_failure_samples;
                }
            }
            output.samples.push_back(std::move(sample));
        }
        output.neighbor_count=neighbors.size();
    }
    if (m_metrics)
    {
        std::size_t entries{},samples{};
        for (const auto & row:rows) entries+=row.size();
        for (const auto & [id,output]:results) { (void)id; samples+=output.samples.size(); }
        m_metrics->estimated_peeling_result_bytes=samples*sizeof(PeelingSampleEstimate);
        m_metrics->estimated_peeling_workspace_bytes=rows.size()*sizeof(std::vector<Contribution>)+
            entries*sizeof(Contribution)+represented.size()*sizeof(std::optional<GaussianModel3D>);
        m_metrics->estimated_peak_postprocessing_payload_bytes=m_metrics->estimated_context_lookup_bytes+
            m_metrics->estimated_peeling_workspace_bytes+m_metrics->estimated_peeling_result_bytes;
    }
    std::sort(results.begin(), results.end(), [](const auto & left, const auto & right) {
        return left.first < right.first;
    });
    return results;
}

PostJointLocalSamples JointPostprocessingContext::BuildPostJointLocalSamples(
    int atom_id, const PostFitPeelingResult & peeling) const
{
    const auto & input = m_problem.OriginalInput();
    const auto atom_found=std::lower_bound(m_contributors.contributors.begin(),m_contributors.contributors.end(),atom_id,
        [](const JointPeelingContributor & contributor,int id){return contributor.atom_id<id;});
    if (atom_found==m_contributors.contributors.end() || atom_found->atom_id!=atom_id)
        throw std::invalid_argument("Post-Joint sample target is outside the Joint problem.");
    const auto atom_index=atom_found->contributor_index;
    if (!input.selection_domain || !std::binary_search(
            input.selection_domain->target_indices.begin(), input.selection_domain->target_indices.end(), atom_index))
        throw std::invalid_argument("Post-Joint samples require a selected Joint target.");

    const auto atom = m_model.FindAtomPtr(atom_id);
    const auto local = AtomLocalPotentialView::For(*atom);
    const auto & stage = local.GetStageEstimate(FittingStage::Second);
    if (stage.source.role != FittingRole::Target || !stage.point ||
        !SameSource(stage.source, peeling.source))
        throw std::invalid_argument("Post-Joint sample endpoint source mismatch.");
    const auto & point = *stage.point;
    PostJointLocalSamples output;
    output.source=stage.source;
    output.fixed_charge={point.GetOffset(),point.GetWidth()};
    const auto observations=CollectPostJointGridObservations(m_geometry,input,atom_index,atom->GetPosition());
    if (peeling.mode!="original-grid-v1" || peeling.samples.size()!=observations.size())
        throw std::invalid_argument("Post-Joint grid/peeling count mismatch.");
    for (std::size_t i=0;i<observations.size();++i)
    {
        const auto & observation=observations[i];
        const auto & peeled=peeling.samples[i];
        if (peeled.observation.original_row_index!=observation.original_row_index ||
            peeled.observation.global_grid_index!=observation.global_grid_index ||
            peeled.observation.position!=observation.position ||
            peeled.observation.squared_distance!=observation.squared_distance ||
            peeled.observation.raw_response!=observation.raw_response)
            throw std::invalid_argument("Post-Joint grid/peeling identity mismatch.");
        PostJointLocalSample sample;
        static_cast<PostJointGridObservation &>(sample)=observation;
        sample.distance=std::sqrt(observation.squared_distance);
        sample.peeled_response=peeled.response;
        if (!peeled.response) sample.reason=peeled.reason.empty() ? "peeling-unavailable" : peeled.reason;
        else if (!std::isfinite(*peeled.response)) sample.reason="nonfinite-peeled-response";
        else if (!std::isfinite(point.GetOffset()) || !std::isfinite(point.GetWidth()) || point.GetWidth()<=0)
            sample.reason="invalid-joint-charge-profile";
        else
        {
            const double charge=point.GetOffset()*joint_component::EvaluateKernel(
                observation.squared_distance,point.GetWidth(),2.5).charge;
            if (!std::isfinite(charge)) sample.reason="nonfinite-fixed-charge-response";
            else
            {
                sample.fixed_charge_response=charge;
                const double response=*peeled.response-charge;
                if (std::isfinite(response)) sample.gaussian_response=response;
                else sample.reason="nonfinite-charge-corrected-response";
            }
        }
        output.samples.push_back(std::move(sample));
    }
    if (m_metrics)
    {
        const auto payload=output.samples.size()*(sizeof(PostJointLocalSample)+sizeof(LocalPotentialSample));
        m_metrics->estimated_peak_postprocessing_payload_bytes=std::max(
            m_metrics->estimated_peak_postprocessing_payload_bytes,
            m_metrics->estimated_context_lookup_bytes+m_metrics->estimated_peeling_result_bytes+payload);
    }
    return output;
}

}
