#include "PostFitPeeling.hpp"
#include "JointPostprocessing.hpp"
#include "MapInterpolation.hpp"
#include "joint_component/Numerics.hpp"
#include "data/detail/AtomClassifier.hpp"
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
    const JointAnalysisResult * snapshot, JointPostprocessingMetrics * metrics)
    : m_geometry(geometry), m_model(model), m_problem(problem), m_snapshot(snapshot), m_metrics(metrics)
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

    m_contributor_index_by_id.reserve(input.atom_ids.size());
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
        m_contributor_index_by_id.emplace(std::stoi(input.atom_ids[atom]), atom);
    if (m_metrics) ++m_metrics->contributor_mapping_builds;
    if (m_metrics)
        m_metrics->estimated_context_lookup_bytes = EstimateHashStorageBytes(m_row_index_by_grid) +
            EstimateHashStorageBytes(m_fitting_row_index_by_grid) +
            EstimateHashStorageBytes(m_contributor_index_by_id);
}

std::map<int, PostFitPeelingResult> JointPostprocessingContext::BuildPostFitPeelingSamples(
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
    const bool fixed_target_only=result &&
        result->parameterization_contract==JointTargetOnlyFixedHaloParameterizationContract;
    const bool pure_exclusion=result &&
        result->parameterization_contract==JointPureHaloSupportExclusionParameterizationContract;
    if (input.selection_domain && input.selection_domain->contributor_policy==PureHaloSupportExclusionContributorPolicy &&
        (!result || !pure_exclusion))
        throw std::invalid_argument("Pure Exclusion peeling requires a Pure Exclusion snapshot.");
    if(fixed_target_only && (!result->selection_domain ||
        result->selection_domain->contributor_policy!=FixedHaloContributorPolicy || result->layout))
        throw std::invalid_argument("Invalid fixed-halo peeling contract.");
    const bool pure_has_halos=pure_exclusion && result->selection_domain &&
        result->atom_ids.size()>result->selection_domain->target_indices.size();
    if(pure_exclusion && (!result->selection_domain ||
        result->selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy ||
        result->layout || !result->fixed_halos.empty() || (pure_has_halos && !result->halo_median_reconstruction)))
        throw std::invalid_argument("Invalid Pure Exclusion peeling contract.");
    if(fixed_target_only || pure_exclusion)
        for(std::size_t atom=0;atom<requested.size();++atom)
            if(!std::binary_search(result->selection_domain->target_indices.begin(),
                result->selection_domain->target_indices.end(),atom)) requested[atom]=false;
    if (result) for (const auto & c : result->components) if (c.state)
    {
        const auto & full = c.layout ? c.layout->full_atoms : c.atoms;
        for (std::size_t k = 0; k < full.size(); ++k)
        {
            if(full[k]>=represented.size() || (fixed_target_only &&
                !std::binary_search(result->selection_domain->target_indices.begin(),
                    result->selection_domain->target_indices.end(),full[k])) || represented[full[k]])
                throw std::invalid_argument("Invalid target-only component mapping in peeling snapshot.");
            represented[full[k]] = GaussianModel3D{c.state->ac[2*k], c.state->b[k], c.state->ac[2*k+1]};
        }
    }
    if(fixed_target_only)
    {
        std::size_t previous{};
        for(std::size_t k=0;k<result->fixed_halos.size();++k)
        {
            const auto & fixed=result->fixed_halos[k];
            if(fixed.contributor_index>=represented.size() || (k && fixed.contributor_index<=previous) ||
                result->atom_ids[fixed.contributor_index]!=fixed.atom_id || represented[fixed.contributor_index] ||
                std::binary_search(result->selection_domain->target_indices.begin(),
                    result->selection_domain->target_indices.end(),fixed.contributor_index) ||
                !std::isfinite(fixed.amplitude) || fixed.amplitude<0 ||
                !std::isfinite(fixed.width) || fixed.width<=0 || fixed.charge!=0.0 ||
                !IsQualifiedFirstStageFixedHalo(fixed) ||
                fixed.source_contract!=FixedHaloParameterSourceContractVersion)
                throw std::invalid_argument("Invalid fixed FirstStage halo in peeling snapshot.");
            previous=fixed.contributor_index;
            represented[fixed.contributor_index]=GaussianModel3D{fixed.amplitude,fixed.width,0.0};
        }
        if(result->assembled_state)
            for(std::size_t atom=0;atom<represented.size();++atom)
                if(!represented[atom])
                    throw std::invalid_argument("Fitted peeling snapshot has an unmapped target or fixed halo.");
    }
    if(pure_exclusion && pure_has_halos)
    {
        const auto & snapshot=*result->halo_median_reconstruction;
        if (snapshot.contract!=HaloMedianReconstructionContractVersion || snapshot.joint_run_id.empty() ||
            snapshot.valid_target_count==0)
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
                halo.donor_count==0 ||
                (halo.source==HaloMedianSource::GlobalTargetMedian && halo.donor_count!=snapshot.valid_target_count) ||
                !valid_median(halo.amplitude,halo.width,halo.charge))
                throw std::invalid_argument("Invalid Pure Exclusion Halo median reconstruction.");
            const auto * atom=m_model.FindAtomPtr(std::stoi(halo.atom_id));
            if (!atom || data_internal::GetGroupKey(atom)!=halo.group_key)
                throw std::invalid_argument("Invalid Pure Exclusion Halo GroupKey provenance.");
            previous=halo.halo_contributor_index;
            first=false;
            represented[halo.halo_contributor_index]=GaussianModel3D{halo.amplitude,halo.width,halo.charge};
        }
        if (result->assembled_state)
            for(std::size_t atom=0;atom<represented.size();++atom)
                if(!represented[atom])
                    throw std::invalid_argument("Fitted Pure Exclusion snapshot has an unmapped target or Halo.");
    }
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
    if (result) for (const auto & component : result->components) if (component.state && component.layout)
        for (std::size_t k = 0; k < component.layout->groups.size(); ++k)
        {
            const auto & group = component.layout->groups[k];
            prediction.at(group.row) += component.state->nuisance_amplitudes.at(k);
            missing.at(group.row) -= group.atoms.size();
        }
    std::map<int, PostFitPeelingResult> results;
    const auto record_sample = [&](const PeelingSampleEstimate & sample, const StencilAudit & audit) {
        if (!m_metrics) return;
        ++m_metrics->raw_peeling_sample_count;
        m_metrics->nonzero_stencil_node_count += audit.nonzero_nodes;
        m_metrics->map_geometry_boundary_nodes += audit.map_boundary_nodes;
        m_metrics->original_support_boundary_nodes += audit.original_support_boundary_nodes;
        m_metrics->fitting_mask_only_intersection_nodes += audit.fitting_mask_only_nodes;
        if (audit.map_boundary) ++m_metrics->map_geometry_boundary_samples;
        else if (audit.mapping_mismatch) ++m_metrics->mapping_mismatch_samples;
        else if (audit.original_support_boundary) ++m_metrics->original_support_boundary_samples;
        else if (sample.reason == "missing-contributor-state") ++m_metrics->missing_contributor_samples;
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
        if (fixed_target_only && view.GetStageEstimate(FittingStage::Second).source.component_id.empty()) continue;
        auto & output = results[id];
        output.source = view.GetStageEstimate(FittingStage::Second).source;
        std::set<std::size_t> neighbors;
        for (const auto & raw : view.GetRawSamplingEntriesView())
        {
            PeelingSampleEstimate sample;
            StencilAudit audit;
            if (result && result->layout && !std::binary_search(
                    result->layout->full_atoms.begin(), result->layout->full_atoms.end(), atom))
            { sample.reason = "observable-contribution-only"; record_sample(sample, audit); output.samples.push_back(sample); continue; }
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
            for (const auto & [grid, weight] : TricubicWeights(stencil))
            {
                if (weight == 0) continue;
                const auto found = m_row_index_by_grid.find(grid);
                if (found == m_row_index_by_grid.end()) { sample.reason = "outside-joint-domain"; break; }
                const auto row = found->second;
                for (const auto neighbor : contributors[row]) if (neighbor != atom) neighbors.insert(neighbor);
                const auto self_missing = !has_point[atom] && own[atom].contains(row) ? 1u : 0u;
                if (missing[row] > self_missing) { sample.reason = "missing-contributor-state"; break; }
            }
            if (sample.reason.empty())
            {
                const double neighbor_prediction = InterpolateTricubic(stencil, [&](const auto & node) {
                    const auto found = m_row_index_by_grid.find(stencil.Index(node));
                    if (found == m_row_index_by_grid.end()) return 0.0; // Zero-weight node, checked above.
                    const auto row = found->second;
                    const auto self = own[atom].find(row);
                    return prediction[row] - (self == own[atom].end() ? 0 : self->second);
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
            contributor_workspace + support_workspace + own_charge_workspace +
            m_metrics->estimated_peeling_result_bytes;
    }
    return results;
}

PostJointLocalSamples JointPostprocessingContext::BuildPostJointLocalSamples(
    int atom_id, const PostFitPeelingResult & peeling) const
{
    if (!m_snapshot) throw std::invalid_argument("Post-Joint samples require a trusted Joint snapshot.");
    const auto & input = m_problem.Input();
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
        for (const auto & [grid, weight] : TricubicWeights(stencil))
        {
            if (weight == 0.0) continue;
            if (!m_row_index_by_grid.contains(grid)) { sample.reason = "outside-joint-domain"; break; }
        }
        if (sample.reason.empty())
        {
            const double charge_response = InterpolateTricubic(stencil, [&](const auto & node) {
                const auto found = m_row_index_by_grid.find(stencil.Index(node));
                if (found == m_row_index_by_grid.end()) return 0.0;
                const auto charge = own_charge.find(found->second);
                return charge == own_charge.end() ? 0.0 : charge->second;
            });
            if (std::isfinite(charge_response))
            {
                sample.fixed_charge_response = charge_response;
                const double gaussian_response = *peeled.response - charge_response;
                if (std::isfinite(gaussian_response)) sample.gaussian_response = gaussian_response;
                else sample.reason = "nonfinite-charge-corrected-response";
            }
            else sample.reason = "nonfinite-fixed-charge-response";
        }
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
            m_metrics->estimated_peeling_result_bytes + local_payload);
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
