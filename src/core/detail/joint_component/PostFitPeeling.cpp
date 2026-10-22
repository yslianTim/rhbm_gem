#include "PostFitPeeling.hpp"
#include "HaloMedianReconstruction.hpp"
#include "../JointPostprocessing.hpp"
#include "../MapInterpolation.hpp"
#include "Numerics.hpp"
#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

namespace rhbm_gem::core::detail {
namespace {
template<class HashMap>
std::size_t EstimateHashStorageBytes(const HashMap & values)
{
    return values.bucket_count() * sizeof(void *) + values.size() *
        (sizeof(typename HashMap::value_type) + 2 * sizeof(void *));
}

bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role && a.fixed_from_first_stage==b.fixed_from_first_stage;
}

struct StencilAudit
{
    std::size_t nonzero_nodes{};
    std::size_t map_boundary_nodes{};
    std::size_t original_support_boundary_nodes{};
    std::size_t fitting_mask_only_nodes{};
    bool map_boundary{};
    bool original_support_boundary{};
    bool mapping_mismatch{};
    bool fitting_mask_only{};
};

struct ReconstructionAtom
{
    int atom_id{};
    std::array<double, 3> position{};
    std::optional<std::size_t> contributor_index;
    std::optional<GaussianModel3D> point;

    const std::array<double, 3> & GetPosition() const { return position; }
};

struct ReconstructionNode
{
    double prediction{};
    std::size_t missing_contributors{};
};

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
        found->parameter_validity == JointCheckStatus::Passed &&
        found->local_conditioning == "Adequate";
}

StencilAudit AuditStencil(const TricubicStencil & stencil, const MapObject & geometry,
    const std::vector<std::string> & original_row_ids,
    const std::unordered_map<std::size_t, std::size_t> & original_rows,
    const std::unordered_map<std::size_t, std::size_t> & fitting_rows)
{
    const auto dimensions = geometry.GetGridSize();
    std::array<std::array<double, 4>, 3> weights{};
    for (std::size_t axis = 0; axis < 3; ++axis)
        for (std::size_t node = 0; node < 4; ++node)
            weights[axis][node] = CubicInterpolate(
                node == 0, node == 1, node == 2, node == 3, stencil.local[axis]);

    StencilAudit audit;
    for (int i = -1; i < 3; ++i) for (int j = -1; j < 3; ++j) for (int k = -1; k < 3; ++k)
    {
        const double weight = weights[0][static_cast<std::size_t>(i + 1)] *
            weights[1][static_cast<std::size_t>(j + 1)] *
            weights[2][static_cast<std::size_t>(k + 1)];
        if (weight == 0.0) continue;
        ++audit.nonzero_nodes;
        const auto raw_node = stencil.UnclampedNode(i, j, k);
        const bool in_bounds = raw_node[0] >= 0 && raw_node[1] >= 0 && raw_node[2] >= 0 &&
            raw_node[0] < dimensions[0] && raw_node[1] < dimensions[1] && raw_node[2] < dimensions[2];
        if (!in_bounds)
        {
            audit.map_boundary = true;
            ++audit.map_boundary_nodes;
            continue;
        }
        const auto grid = stencil.Index(raw_node);
        const auto found = original_rows.find(grid);
        if (found == original_rows.end())
        {
            audit.original_support_boundary = true;
            ++audit.original_support_boundary_nodes;
            continue;
        }
        if (found->second >= original_row_ids.size() || original_row_ids[found->second] != std::to_string(grid))
        {
            audit.mapping_mismatch = true;
            continue;
        }
        if (!fitting_rows.contains(grid))
        {
            audit.fitting_mask_only = true;
            ++audit.fitting_mask_only_nodes;
        }
    }
    return audit;
}
}

LocalPotentialSampleList PostJointLocalSamples::BuildGaussianSamples(
    const LocalPotentialSampleList & raw_samples) const
{
    if (samples.size() != raw_samples.size())
        throw std::invalid_argument("Post-Joint local/raw sample count mismatch.");
    LocalPotentialSampleList result;
    result.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        const auto & prepared = samples[i];
        const auto & raw = raw_samples[i];
        if (prepared.raw_sample_index != i || prepared.distance != raw.point.distance ||
            prepared.raw_response != raw.response)
            throw std::invalid_argument("Post-Joint local/raw sample identity mismatch.");
        if (!prepared.gaussian_response || !raw.point.is_selected) continue;
        auto sample = raw;
        sample.response = *prepared.gaussian_response;
        result.push_back(std::move(sample));
    }
    return result;
}

JointPostprocessingContext::JointPostprocessingContext(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    const JointAnalysisResult * snapshot, JointPostprocessingMetrics * metrics,
    const HaloMedianBuildResult * halo_medians)
    : JointPostprocessingContext(geometry,model,problem,snapshot,metrics,halo_medians,
        BuildJointPeelingContributorUniverse(model,problem,snapshot))
{
}

JointPostprocessingContext::JointPostprocessingContext(
    const MapObject & geometry,const ModelObject & model,const JointProblem & problem,
    const JointAnalysisResult * snapshot,JointPostprocessingMetrics * metrics,
    const HaloMedianBuildResult * halo_medians,JointPeelingContributorUniverse contributors)
    : m_geometry(geometry), m_model(model), m_problem(problem), m_snapshot(snapshot), m_metrics(metrics),
      m_halo_medians(halo_medians),m_contributors(std::move(contributors))
{
    const auto & input = problem.OriginalInput();
    const bool pure_exclusion=input.selection_domain &&
        input.selection_domain->contributor_policy==PureHaloSupportExclusionContributorPolicy;
    if (snapshot)
    {
        if (m_metrics) ++m_metrics->snapshot_identity_validation_count;
        if (snapshot->atom_ids != input.atom_ids)
            throw std::invalid_argument("Joint peeling snapshot identity mismatch.");
        if (m_metrics) ++m_metrics->full_row_identity_scan_count;
        if (snapshot->row_ids != input.row_ids)
            throw std::invalid_argument("Joint peeling snapshot identity mismatch.");
        if (pure_exclusion)
        {
            if (snapshot->parameterization_contract!=JointPureHaloSupportExclusionParameterizationContract ||
                !snapshot->selection_domain ||
                snapshot->selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy ||
                snapshot->available_row_mask.size()!=input.row_ids.size() ||
                snapshot->fitting_row_ids!=problem.Input().row_ids ||
                snapshot->fitting_to_original_contributors!=problem.FittingToOriginalContributors() ||
                snapshot->fitting_to_original_rows!=problem.FittingToOriginalRows())
                throw std::invalid_argument("Joint peeling fitting/reconstruction identity mismatch.");
            for (std::size_t row=0;row<snapshot->fitting_row_ids.size();++row)
            {
                const auto original_row=snapshot->fitting_to_original_rows[row];
                if (original_row>=input.row_ids.size() || snapshot->fitting_row_ids[row]!=input.row_ids[original_row])
                    throw std::invalid_argument("Joint peeling fitting/reconstruction identity mismatch.");
            }
        }
    }
    m_contributor_index_by_id.reserve(input.atom_ids.size());
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
        if (contributor.contributor_index)
        {
            const auto index=*contributor.contributor_index;
            if (index>=input.atom_ids.size() || std::stoi(input.atom_ids[index])!=contributor.atom_id)
                throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
            const bool target=input.selection_domain && std::binary_search(
                input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),index);
            if (contributor.target!=target || !contributor.reconstruction_eligible ||
                !m_contributor_index_by_id.emplace(contributor.atom_id,index).second)
                throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
        }
        else if (contributor.target)
            throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
    }
    if (m_contributor_index_by_id.size()!=input.atom_ids.size())
        throw std::invalid_argument("Joint Peeling contributor universe identity mismatch.");
    m_row_index_by_grid.reserve(input.row_ids.size());
    for (std::size_t row = 0; row < input.row_ids.size(); ++row)
    {
        const auto grid = std::stoull(input.row_ids[row]);
        if (grid >= geometry.GetMapValueArraySize()) throw std::invalid_argument("Joint row outside map geometry.");
        const auto [unused, inserted] = m_row_index_by_grid.emplace(grid, row);
        (void)unused;
        if (!inserted || input.row_ids[row] != std::to_string(grid))
            throw std::invalid_argument("Joint original-domain row mapping mismatch.");
    }
    if (m_metrics) ++m_metrics->observation_row_mapping_builds;

    const auto & fitting_input = problem.Input();
    m_fitting_row_index_by_grid.reserve(fitting_input.row_ids.size());
    for (std::size_t row = 0; row < fitting_input.row_ids.size(); ++row)
    {
        const auto grid = std::stoull(fitting_input.row_ids[row]);
        const auto [unused, inserted] = m_fitting_row_index_by_grid.emplace(grid, row);
        (void)unused;
        if (!inserted || fitting_input.row_ids[row] != std::to_string(grid))
            throw std::invalid_argument("Joint fitting-domain row mapping mismatch.");
    }

    if (m_metrics) ++m_metrics->contributor_mapping_builds;
    if (m_metrics)
        m_metrics->estimated_context_lookup_bytes = EstimateHashStorageBytes(m_row_index_by_grid) +
            EstimateHashStorageBytes(m_fitting_row_index_by_grid) +
            EstimateHashStorageBytes(m_contributor_index_by_id) +
            m_contributors.contributors.size()*sizeof(JointPeelingContributor);
}

std::vector<std::pair<int, PostFitPeelingResult>> JointPostprocessingContext::BuildPostFitPeelingBatch(
    std::optional<std::span<const std::size_t>> outputs) const
{
    const auto & input = m_problem.OriginalInput();
    const auto * result = m_snapshot;
    auto requested = JointOutputMask(input, outputs);
    std::vector<std::vector<std::size_t>> contributors(input.row_ids.size());
    std::vector<double> prediction(input.row_ids.size());
    std::vector<std::size_t> missing(input.row_ids.size());
    std::vector<std::unordered_map<std::size_t, double>> own(input.atom_ids.size());
    // Internal representatives reconstruct the observation prediction even when
    // a nonunique halo has no public atom-wise point estimate.
    std::vector<std::optional<GaussianModel3D>> represented(input.atom_ids.size());
    if (result && result->parameterization_contract!=JointPureHaloSupportExclusionParameterizationContract)
        throw std::invalid_argument("Post-fit peeling requires a target-only Pure Exclusion snapshot.");
    const bool pure_exclusion=result &&
        result->parameterization_contract==JointPureHaloSupportExclusionParameterizationContract;
    if (input.selection_domain && input.selection_domain->contributor_policy==PureHaloSupportExclusionContributorPolicy &&
        (!result || !pure_exclusion))
        throw std::invalid_argument("Pure Exclusion peeling requires a Pure Exclusion snapshot.");
    const bool pure_has_halos=pure_exclusion && result->selection_domain &&
        result->atom_ids.size()>result->selection_domain->target_indices.size();
    if(pure_exclusion && (!result->selection_domain ||
        result->selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy ||
        (pure_has_halos && !result->halo_median_reconstruction)))
        throw std::invalid_argument("Invalid Pure Exclusion peeling contract.");
    if (pure_exclusion && !pure_has_halos && result->halo_median_reconstruction)
        throw std::invalid_argument("Pure Exclusion Halo medians exist without original-domain Halos.");
    if(pure_exclusion)
        for(std::size_t atom=0;atom<requested.size();++atom)
            if(!std::binary_search(result->selection_domain->target_indices.begin(),
                result->selection_domain->target_indices.end(),atom)) requested[atom]=false;
    if (result) for (const auto & c : result->components) if (c.state)
    {
        const auto & full = c.layout ? c.layout->full_atoms : c.atoms;
        for (std::size_t k = 0; k < full.size(); ++k)
        {
            if(full[k]>=represented.size() || represented[full[k]])
                throw std::invalid_argument("Invalid target-only component mapping in peeling snapshot.");
            represented[full[k]] = GaussianModel3D{c.state->ac[2*k], c.state->b[k], c.state->ac[2*k+1]};
        }
    }
    if(pure_exclusion && pure_has_halos)
    {
        const auto & snapshot=*result->halo_median_reconstruction;
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
                result->atom_ids[halo.halo_contributor_index]!=halo.atom_id ||
                std::binary_search(result->selection_domain->target_indices.begin(),
                    result->selection_domain->target_indices.end(),halo.halo_contributor_index) ||
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
        if (result->assembled_state)
            for(std::size_t atom=0;atom<represented.size();++atom)
                if(!represented[atom])
                    throw std::invalid_argument("Fitted Pure Exclusion snapshot has an unmapped target or Halo.");
    }

    std::string expected_joint_run_id;
    if (pure_exclusion && result->target_evidence)
    {
        const auto target_index=result->selection_domain->target_indices.front();
        const auto * target=m_model.FindAtomPtr(std::stoi(result->atom_ids[target_index]));
        const auto stage=AtomLocalPotentialView::For(*target).GetStageEstimate(FittingStage::Second);
        if (stage.source.method!=EstimateMethod::JointComponents || stage.source.role!=FittingRole::Target ||
            stage.source.run_id.empty())
            throw std::invalid_argument("Invalid Pure Exclusion Joint run provenance.");
        expected_joint_run_id=stage.source.run_id;
        if (result->halo_median_reconstruction &&
            result->halo_median_reconstruction->joint_run_id!=expected_joint_run_id)
            throw std::invalid_argument("Invalid Pure Exclusion Halo median run provenance.");
        if (!m_contributors.joint_run_id.empty() && m_contributors.joint_run_id!=expected_joint_run_id)
            throw std::invalid_argument("Invalid Joint Peeling contributor universe run provenance.");
    }

    std::optional<HaloMedianBuildResult> rebuilt_halo_medians;
    const HaloMedianBuildResult * halo_medians=m_halo_medians;
    if (pure_exclusion && !halo_medians && result->target_evidence)
    {
        try
        {
            rebuilt_halo_medians=BuildHaloMedianReconstruction(
                m_model,m_problem,*result,expected_joint_run_id);
            halo_medians=&*rebuilt_halo_medians;
        }
        catch (const std::invalid_argument & error)
        {
            if (std::string_view(error.what())!="halo-median-no-valid-targets") throw;
            if (result->halo_median_reconstruction)
                throw std::invalid_argument("Saved Pure Exclusion Halo medians have no qualified donors.");
        }
    }
    if (pure_exclusion && result->halo_median_reconstruction && halo_medians && halo_medians->snapshot)
    {
        const auto & saved=*result->halo_median_reconstruction;
        const auto & rebuilt=*halo_medians->snapshot;
        if (saved.contract!=HaloMedianReconstructionContractVersion || saved.joint_run_id.empty() ||
            saved.joint_run_id!=rebuilt.joint_run_id || saved.valid_target_count!=rebuilt.valid_target_count ||
            saved.global_median!=rebuilt.global_median || saved.neighbor_count!=rebuilt.neighbor_count ||
            saved.donor_policy!=rebuilt.donor_policy || saved.halos!=rebuilt.halos)
            throw std::invalid_argument("Invalid Pure Exclusion Halo median run provenance.");
    }
    const std::string * halo_run_id=halo_medians && halo_medians->snapshot ?
        &halo_medians->snapshot->joint_run_id : result && result->halo_median_reconstruction ?
        &result->halo_median_reconstruction->joint_run_id : nullptr;
    if (halo_run_id && !expected_joint_run_id.empty() && *halo_run_id!=expected_joint_run_id)
        throw std::invalid_argument("Invalid Pure Exclusion reconstruction median run provenance.");

    const std::vector<HaloMedianTargetDonor> no_target_donors;
    const auto & target_donors=halo_medians ? halo_medians->target_donors : no_target_donors;
    const GeometryLocalTargetDonorIndex target_donor_index(target_donors);

    constexpr std::size_t reconstruction_node_cache_limit=4096;
    std::vector<ReconstructionAtom> reconstruction_atoms;
    std::unordered_map<int,std::size_t> reconstruction_atom_by_id;
    std::unique_ptr<KDNode<ReconstructionAtom>> reconstruction_tree;
    if (pure_exclusion)
    {
        reconstruction_atoms.reserve(m_contributors.contributors.size());
        reconstruction_atom_by_id.reserve(m_contributors.contributors.size());
        for (const auto & source:m_contributors.contributors)
        {
            if (!source.reconstruction_eligible) continue;
            ReconstructionAtom candidate;
            candidate.atom_id=source.atom_id;
            candidate.position=source.position;
            candidate.contributor_index=source.contributor_index;
            if (source.contributor_index)
            {
                const bool target=std::binary_search(result->selection_domain->target_indices.begin(),
                    result->selection_domain->target_indices.end(),*source.contributor_index);
                if (!target || QualifiedTargetPoint(*result,*source.contributor_index))
                    candidate.point=represented[*source.contributor_index];
            }
            if (!candidate.point && !candidate.contributor_index && !target_donors.empty())
            {
                const auto selected=target_donor_index.Query(candidate.position);
                if (!selected.empty())
                {
                    const auto median=GeometryLocalTargetMedian(selected);
                    candidate.point=GaussianModel3D{median[0],median[1],median[2]};
                }
            }
            if (candidate.point &&
                (!std::isfinite(candidate.point->GetAmplitude()) || candidate.point->GetAmplitude()<0 ||
                    !std::isfinite(candidate.point->GetWidth()) || candidate.point->GetWidth()<=0 ||
                    !std::isfinite(candidate.point->GetOffset())))
                throw std::invalid_argument("Invalid Pure Exclusion reconstruction contributor parameters.");
            reconstruction_atom_by_id.emplace(candidate.atom_id,reconstruction_atoms.size());
            reconstruction_atoms.push_back(std::move(candidate));
        }
        reconstruction_tree=KDTreeAlgorithm<ReconstructionAtom>::BuildKDTree(reconstruction_atoms);
        if (m_metrics)
        {
            m_metrics->reconstruction_lookup_atom_count=reconstruction_atoms.size();
            std::size_t donor_index_workspace_bytes{};
            for (const auto & donor:target_donors)
                donor_index_workspace_bytes+=sizeof(HaloMedianTargetDonor)+donor.atom_id.capacity()+
                    sizeof(std::array<double,3>)+sizeof(HaloMedianTargetDonor *)+8*sizeof(void *);
            m_metrics->estimated_reconstruction_workspace_bytes=
                reconstruction_atoms.size()*(sizeof(ReconstructionAtom)+sizeof(ReconstructionAtom *)+
                    sizeof(KDNode<ReconstructionAtom>))+donor_index_workspace_bytes+
                reconstruction_node_cache_limit*(sizeof(std::pair<const std::size_t,ReconstructionNode>)+
                    2*sizeof(void *))+
                EstimateHashStorageBytes(reconstruction_atom_by_id)+
                64*(sizeof(std::pair<const std::size_t,double>)+2*sizeof(void *));
        }
    }
    std::unordered_map<std::size_t,ReconstructionNode> reconstruction_node_cache;
    if (pure_exclusion) reconstruction_node_cache.reserve(reconstruction_node_cache_limit);
    const auto reconstruct_node=[&](std::size_t grid) {
        if (const auto found=reconstruction_node_cache.find(grid);found!=reconstruction_node_cache.end())
        {
            if (m_metrics) ++m_metrics->reconstruction_node_cache_hit_count;
            return found->second;
        }
        ReconstructionNode node;
        if (reconstruction_tree)
        {
            ReconstructionAtom query;
            query.position=m_geometry.GetGridPosition(grid);
            const auto neighbors=KDTreeAlgorithm<ReconstructionAtom>::RangeSearch(
                reconstruction_tree.get(),&query,2.5);
            for (const auto * candidate:neighbors)
            {
                const double square=SquaredDistance(query.position,candidate->position);
                if (square>6.25) continue;
                if (!candidate->point)
                {
                    ++node.missing_contributors;
                    continue;
                }
                const auto basis=joint_component::EvaluateKernel(square,candidate->point->GetWidth(),2.5);
                node.prediction+=candidate->point->GetAmplitude()*basis.gaussian+
                    candidate->point->GetOffset()*basis.charge;
            }
        }
        if (m_metrics) ++m_metrics->reconstruction_only_node_evaluation_count;
        if (reconstruction_node_cache.size()<reconstruction_node_cache_limit)
        {
            reconstruction_node_cache.emplace(grid,node);
            if (m_metrics) m_metrics->reconstruction_node_cache_peak_entries=
                std::max(m_metrics->reconstruction_node_cache_peak_entries,reconstruction_node_cache.size());
        }
        return node;
    };

    std::vector<bool> has_point;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        const auto view = AtomLocalPotentialView::For(*m_model.FindAtomPtr(std::stoi(input.atom_ids[atom])));
        const auto & point = result ? represented[atom] : view.GetStageEstimate(FittingStage::Second).point;
        has_point.push_back(point.has_value());
        for (const auto & support : input.support[atom])
        {
            contributors[support.row].push_back(atom);
            double value = 0;
            if (point)
            {
                const auto basis = joint_component::EvaluateKernel(support.squared_distance, point->GetWidth(), 2.5);
                value = point->GetAmplitude() * basis.gaussian + point->GetOffset() * basis.charge;
                prediction[support.row] += value;
            }
            else ++missing[support.row];
            if (requested[atom]) own[atom].emplace(support.row, value);
        }
    }
    std::vector<std::pair<int, PostFitPeelingResult>> results;
    results.reserve(static_cast<std::size_t>(std::count(requested.begin(), requested.end(), true)));
    const auto record_sample = [&](const PeelingSampleEstimate & sample, const StencilAudit & audit) {
        if (!m_metrics) return;
        ++m_metrics->raw_peeling_sample_count;
        m_metrics->nonzero_stencil_node_count += audit.nonzero_nodes;
        m_metrics->map_geometry_boundary_nodes += audit.map_boundary_nodes;
        m_metrics->original_support_boundary_nodes += audit.original_support_boundary_nodes;
        m_metrics->fitting_mask_only_intersection_nodes += audit.fitting_mask_only_nodes;
        if (audit.map_boundary) ++m_metrics->map_geometry_boundary_samples;
        if (audit.mapping_mismatch) ++m_metrics->mapping_mismatch_samples;
        if (audit.original_support_boundary) ++m_metrics->original_support_boundary_samples;
        if (sample.reason == "missing-contributor-state") ++m_metrics->missing_contributor_samples;
        else if (!sample.response) ++m_metrics->other_peeling_failure_samples;
        if (audit.fitting_mask_only) ++m_metrics->fitting_mask_only_intersection_samples;
        if (sample.response) ++m_metrics->available_peeling_sample_count;
        else ++m_metrics->unavailable_peeling_sample_count;
    };
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        if (!requested[atom]) continue;
        const auto id = std::stoi(input.atom_ids[atom]);
        const auto view = AtomLocalPotentialView::For(*m_model.FindAtomPtr(id));
        auto & output = results.emplace_back(id, PostFitPeelingResult{}).second;
        output.source = view.GetStageEstimate(FittingStage::Second).source;
        std::set<int> neighbors;
        for (const auto & raw : view.GetRawSamplingEntriesView())
        {
            PeelingSampleEstimate sample;
            StencilAudit audit;
            if (result)
            {
                const auto & estimate = view.GetStageEstimate(FittingStage::Second);
                if (!estimate.point)
                { sample.reason = estimate.reason; record_sample(sample, audit); output.samples.push_back(sample); continue; }
            }
            if (!view.HasSampleGeometry())
            { sample.reason = "sample-geometry-unavailable"; record_sample(sample, audit); output.samples.push_back(sample); continue; }
            const auto stencil = MakeTricubicStencil(m_geometry, raw.point.position);
            audit = AuditStencil(stencil, m_geometry, input.row_ids, m_row_index_by_grid,
                m_fitting_row_index_by_grid);
            std::unordered_map<std::size_t,double> node_predictions;
            node_predictions.reserve(64);
            for (const auto & [grid, weight] : TricubicWeights(stencil))
            {
                if (weight == 0) continue;
                const auto found = m_row_index_by_grid.find(grid);
                if (found == m_row_index_by_grid.end())
                {
                    if (!pure_exclusion) { sample.reason = "outside-joint-domain"; break; }
                    const auto reconstructed=reconstruct_node(grid);
                    auto target_candidate=reconstruction_atom_by_id.find(id);
                    const auto grid_position=m_geometry.GetGridPosition(grid);
                    const bool self_missing=target_candidate!=reconstruction_atom_by_id.end() &&
                        !reconstruction_atoms[target_candidate->second].point &&
                        SquaredDistance(reconstruction_atoms[target_candidate->second].position,grid_position)<=6.25;
                    if (reconstructed.missing_contributors>static_cast<std::size_t>(self_missing))
                    { sample.reason = "missing-contributor-state"; break; }
                    double self_prediction{};
                    if (target_candidate!=reconstruction_atom_by_id.end())
                    {
                        const auto & candidate=reconstruction_atoms[target_candidate->second];
                        if (candidate.point)
                        {
                            const auto basis=joint_component::EvaluateKernel(
                                SquaredDistance(candidate.position,grid_position),candidate.point->GetWidth(),2.5);
                            self_prediction=candidate.point->GetAmplitude()*basis.gaussian+
                                candidate.point->GetOffset()*basis.charge;
                        }
                    }
                    node_predictions.emplace(grid,reconstructed.prediction-self_prediction);
                    ReconstructionAtom query;
                    query.position=grid_position;
                    for (const auto * neighbor:KDTreeAlgorithm<ReconstructionAtom>::RangeSearch(
                            reconstruction_tree.get(),&query,2.5))
                        if (neighbor->atom_id!=id)
                            neighbors.insert(neighbor->atom_id);
                    continue;
                }
                const auto row = found->second;
                for (const auto neighbor : contributors[row])
                    if (neighbor != atom) neighbors.insert(std::stoi(input.atom_ids[neighbor]));
                const auto self_missing = !has_point[atom] && own[atom].contains(row) ? 1u : 0u;
                if (missing[row] > self_missing) { sample.reason = "missing-contributor-state"; break; }
                const auto self = own[atom].find(row);
                node_predictions.emplace(grid,prediction[row]-(self==own[atom].end() ? 0.0 : self->second));
            }
            if (sample.reason.empty())
            {
                const double neighbor_prediction = InterpolateTricubic(stencil, [&](const auto & node) {
                    const auto grid=stencil.Index(node);
                    const auto found=node_predictions.find(grid);
                    return found==node_predictions.end() ? 0.0 : found->second;
                });
                const double response = raw.response - neighbor_prediction;
                if (std::isfinite(response)) sample.response = response;
                else sample.reason = "nonfinite-response";
            }
            record_sample(sample, audit);
            output.samples.push_back(std::move(sample));
        }
        output.neighbor_count = neighbors.size();
    }
    if (m_metrics)
    {
        std::size_t support_entries{};
        std::size_t requested_support_entries{};
        std::size_t output_sample_count{};
        for (std::size_t atom = 0; atom < input.support.size(); ++atom)
        {
            support_entries += input.support[atom].size();
            if (requested[atom]) requested_support_entries += input.support[atom].size();
        }
        for (const auto & [unused, output] : results)
        {
            (void)unused;
            output_sample_count += output.samples.size();
        }
        m_metrics->estimated_peeling_result_bytes = output_sample_count * sizeof(PeelingSampleEstimate);
        const auto row_workspace = input.row_ids.size() *
            (sizeof(std::vector<std::size_t>) + sizeof(double) + sizeof(std::size_t));
        const auto contributor_workspace = input.atom_ids.size() *
            (sizeof(std::unordered_map<std::size_t, double>) + sizeof(std::optional<GaussianModel3D>));
        const auto support_workspace = support_entries * sizeof(std::size_t);
        const auto own_charge_workspace = requested_support_entries *
            (sizeof(std::pair<const std::size_t, double>) + 2 * sizeof(void *));
        m_metrics->estimated_peak_postprocessing_payload_bytes = row_workspace +
            m_metrics->estimated_context_lookup_bytes + contributor_workspace + support_workspace + own_charge_workspace +
            m_metrics->estimated_peeling_result_bytes + m_metrics->estimated_reconstruction_workspace_bytes;
    }
    std::sort(results.begin(), results.end(), [](const auto & left, const auto & right) {
        return left.first < right.first;
    });
    return results;
}

std::map<int, PostFitPeelingResult> JointPostprocessingContext::BuildPostFitPeelingSamples(
    std::optional<std::span<const std::size_t>> outputs) const
{
    std::map<int, PostFitPeelingResult> results;
    for (auto & [id, peeling] : BuildPostFitPeelingBatch(outputs))
        results.emplace(id, std::move(peeling));
    return results;
}

PostJointLocalSamples JointPostprocessingContext::BuildPostJointLocalSamples(
    int atom_id, const PostFitPeelingResult & peeling) const
{
    if (!m_snapshot) throw std::invalid_argument("Post-Joint samples require a trusted Joint snapshot.");
    const auto & input = m_problem.OriginalInput();
    const auto atom_found = m_contributor_index_by_id.find(atom_id);
    if (atom_found == m_contributor_index_by_id.end())
        throw std::invalid_argument("Post-Joint sample target is outside the Joint problem.");
    const auto atom_index = atom_found->second;
    if (!input.selection_domain || !std::binary_search(
            input.selection_domain->target_indices.begin(), input.selection_domain->target_indices.end(), atom_index))
        throw std::invalid_argument("Post-Joint samples require a selected Joint target.");
    if (m_metrics) ++m_metrics->target_preparation_count;

    const auto atom = m_model.FindAtomPtr(atom_id);
    const auto local = AtomLocalPotentialView::For(*atom);
    const auto & stage = local.GetStageEstimate(FittingStage::Second);
    if (stage.source.role != FittingRole::Target || !stage.point ||
        !SameSource(stage.source, peeling.source))
        throw std::invalid_argument("Post-Joint sample endpoint source mismatch.");
    const auto & point = *stage.point;
    const auto & raw = local.GetRawSamplingEntriesView();
    if (raw.size() != peeling.samples.size())
        throw std::invalid_argument("Post-Joint sample/peeling count mismatch.");
    if (!local.HasSampleGeometry()) throw std::invalid_argument("Post-Joint sample geometry is unavailable.");

    std::unordered_map<std::size_t, double> own_charge;
    own_charge.reserve(input.support[atom_index].size());
    for (const auto & support : input.support[atom_index])
    {
        const auto basis = joint_component::EvaluateKernel(support.squared_distance, point.GetWidth(), 2.5);
        own_charge.emplace(support.row, point.GetOffset() * basis.charge);
    }

    PostJointLocalSamples output;
    output.source = stage.source;
    output.fixed_charge = {point.GetOffset(), point.GetWidth()};
    output.samples.reserve(raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i)
    {
        const auto & raw_sample = raw[i];
        const auto & peeled = peeling.samples[i];
        PostJointLocalSample sample;
        sample.raw_sample_index = i;
        sample.distance = raw_sample.point.distance;
        sample.raw_response = raw_sample.response;
        sample.peeled_response = peeled.response;
        if (!peeled.response)
        {
            sample.reason = peeled.reason.empty() ? "peeling-unavailable" : peeled.reason;
            output.samples.push_back(std::move(sample));
            continue;
        }
        if (!std::all_of(raw_sample.point.position.begin(), raw_sample.point.position.end(),
                [](double coordinate) { return std::isfinite(coordinate); }))
        {
            sample.reason = "sample-geometry-unavailable";
            output.samples.push_back(std::move(sample));
            continue;
        }
        const auto stencil = MakeTricubicStencil(m_geometry, raw_sample.point.position);
        const double charge_response = InterpolateTricubic(stencil, [&](const auto & node) {
            const auto grid=stencil.Index(node);
            const auto found = m_row_index_by_grid.find(grid);
            if (found != m_row_index_by_grid.end())
            {
                const auto charge = own_charge.find(found->second);
                return charge == own_charge.end() ? 0.0 : charge->second;
            }
            const auto basis=joint_component::EvaluateKernel(
                SquaredDistance(atom->GetPosition(),m_geometry.GetGridPosition(grid)),point.GetWidth(),2.5);
            return point.GetOffset()*basis.charge;
        });
        if (std::isfinite(charge_response))
        {
            sample.fixed_charge_response = charge_response;
            const double gaussian_response = *peeled.response - charge_response;
            if (std::isfinite(gaussian_response)) sample.gaussian_response = gaussian_response;
            else sample.reason = "nonfinite-charge-corrected-response";
        }
        else sample.reason = "nonfinite-fixed-charge-response";
        output.samples.push_back(std::move(sample));
    }
    if (m_metrics)
    {
        m_metrics->total_prepared_sample_count += output.samples.size();
        for (const auto & sample : output.samples)
            if (sample.gaussian_response) ++m_metrics->valid_local_sample_count;
            else ++m_metrics->invalid_local_sample_count;
        const auto local_payload = output.samples.size() * sizeof(PostJointLocalSample) +
            raw.size() * (sizeof(LocalPotentialSample) + 5 * sizeof(double)) +
            input.support[atom_index].size() *
                (sizeof(std::pair<const std::size_t, double>) + 2 * sizeof(void *));
        m_metrics->estimated_peak_postprocessing_payload_bytes = std::max(
            m_metrics->estimated_peak_postprocessing_payload_bytes,
            m_metrics->estimated_peeling_result_bytes + local_payload +
                m_metrics->estimated_context_lookup_bytes +
                m_metrics->estimated_reconstruction_workspace_bytes);
    }
    return output;
}

std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    std::optional<std::span<const std::size_t>> outputs, const JointAnalysisResult * result)
{
    return JointPostprocessingContext(geometry, model, problem, result)
        .BuildPostFitPeelingSamples(outputs);
}

PostJointLocalSamples BuildPostJointLocalSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    int atom_id, const PostFitPeelingResult & peeling, const JointAnalysisResult & result)
{
    return JointPostprocessingContext(geometry, model, problem, &result)
        .BuildPostJointLocalSamples(atom_id, peeling);
}
}
