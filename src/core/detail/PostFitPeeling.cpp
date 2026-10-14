#include "PostFitPeeling.hpp"
#include "JointPostprocessing.hpp"
#include "MapInterpolation.hpp"
#include "joint_component/Numerics.hpp"
#include "data/detail/JointStageAdapter.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

namespace rhbm_gem::core::detail {
namespace {
bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role;
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

std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    std::optional<std::span<const std::size_t>> outputs,const JointAnalysisResult * result)
{
    const auto & input = problem.Input();
    if(result && (result->atom_ids!=input.atom_ids || result->row_ids!=input.row_ids))
        throw std::invalid_argument("Joint peeling snapshot identity mismatch.");
    const auto estimates=result ? data_internal::BuildJointStageEstimates(*result,{}) : std::map<int,LocalStageEstimate>{};
    const auto requested = JointOutputMask(input, outputs);
    std::unordered_map<std::size_t, std::size_t> rows;
    for (std::size_t row = 0; row < input.row_ids.size(); ++row)
    {
        const auto grid = std::stoull(input.row_ids[row]);
        if (grid >= geometry.GetMapValueArraySize()) throw std::invalid_argument("Joint row outside map geometry.");
        rows.emplace(grid, row);
    }
    std::vector<std::vector<std::size_t>> contributors(input.row_ids.size());
    std::vector<double> prediction(input.row_ids.size());
    std::vector<std::size_t> missing(input.row_ids.size());
    std::vector<std::unordered_map<std::size_t, double>> own(input.atom_ids.size());
    // Internal representatives reconstruct the observation prediction even when
    // a nonunique halo has no public atom-wise point estimate.
    std::vector<std::optional<GaussianModel3D>> represented(input.atom_ids.size());
    if(result) for(const auto & c:result->components) if(c.state)
    {
        const auto & full=c.layout ? c.layout->full_atoms : c.atoms;
        for(std::size_t k=0;k<full.size();++k)
            represented[full[k]]=GaussianModel3D{c.state->ac[2*k],c.state->b[k],c.state->ac[2*k+1]};
    }
    std::vector<bool> has_point;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        const auto view = AtomLocalPotentialView::For(*model.FindAtomPtr(std::stoi(input.atom_ids[atom])));
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
    if(result) for(const auto & component:result->components) if(component.state && component.layout)
        for(std::size_t k=0;k<component.layout->groups.size();++k)
        {
            const auto & group=component.layout->groups[k];
            prediction.at(group.row)+=component.state->nuisance_amplitudes.at(k);
            missing.at(group.row)-=group.atoms.size();
        }
    std::map<int, PostFitPeelingResult> results;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        if (!requested[atom]) continue;
        const auto id = std::stoi(input.atom_ids[atom]);
        const auto view = AtomLocalPotentialView::For(*model.FindAtomPtr(id));
        auto & output = results[id];
        output.source = view.GetStageEstimate(FittingStage::Second).source;
        std::set<std::size_t> neighbors;
        for (const auto & raw : view.GetRawSamplingEntriesView())
        {
            PeelingSampleEstimate sample;
            if(result && result->layout && !std::binary_search(result->layout->full_atoms.begin(),result->layout->full_atoms.end(),atom))
            {sample.reason="observable-contribution-only"; output.samples.push_back(sample); continue;}
            if(result && !estimates.at(id).point)
            {sample.reason=estimates.at(id).reason; output.samples.push_back(sample); continue;}
            if (!view.HasSampleGeometry()) { sample.reason = "sample-geometry-unavailable"; output.samples.push_back(sample); continue; }
            const auto stencil = MakeTricubicStencil(geometry, raw.point.position);
            for (const auto & [grid, weight] : TricubicWeights(stencil))
            {
                if (weight == 0) continue;
                const auto found = rows.find(grid);
                if (found == rows.end()) { sample.reason = "outside-joint-domain"; break; }
                const auto row = found->second;
                for (const auto neighbor : contributors[row]) if (neighbor != atom) neighbors.insert(neighbor);
                const auto self_missing = !has_point[atom] && own[atom].contains(row) ? 1u : 0u;
                if (missing[row] > self_missing) { sample.reason = "missing-contributor-state"; break; }
            }
            if (sample.reason.empty())
            {
                const double neighbor_prediction = InterpolateTricubic(stencil, [&](const auto & node) {
                    const auto found = rows.find(stencil.Index(node));
                    if (found == rows.end()) return 0.0; // Zero-weight node, checked above.
                    const auto row = found->second;
                    const auto self = own[atom].find(row);
                    return prediction[row] - (self == own[atom].end() ? 0 : self->second);
                });
                const double response = raw.response - neighbor_prediction;
                if (std::isfinite(response)) sample.response = response;
                else sample.reason = "nonfinite-response";
            }
            output.samples.push_back(std::move(sample));
        }
        output.neighbor_count = neighbors.size();
    }
    return results;
}

PostJointLocalSamples BuildPostJointLocalSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    int atom_id, const PostFitPeelingResult & peeling, const JointAnalysisResult & result)
{
    const auto & input = problem.Input();
    if (result.atom_ids != input.atom_ids || result.row_ids != input.row_ids)
        throw std::invalid_argument("Post-Joint sample snapshot identity mismatch.");
    const auto atom_found = std::find(input.atom_ids.begin(), input.atom_ids.end(), std::to_string(atom_id));
    if (atom_found == input.atom_ids.end())
        throw std::invalid_argument("Post-Joint sample target is outside the Joint problem.");
    const auto atom_index = static_cast<std::size_t>(std::distance(input.atom_ids.begin(), atom_found));
    if (!input.selection_domain || !std::binary_search(
            input.selection_domain->target_indices.begin(), input.selection_domain->target_indices.end(), atom_index))
        throw std::invalid_argument("Post-Joint samples require a selected Joint target.");

    const auto estimates = data_internal::BuildJointStageEstimates(result, peeling.source.run_id);
    const auto stage_found = estimates.find(atom_id);
    if (stage_found == estimates.end() || stage_found->second.source.role != FittingRole::Target ||
        !stage_found->second.point || !SameSource(stage_found->second.source, peeling.source))
        throw std::invalid_argument("Post-Joint sample endpoint source mismatch.");
    const auto & point = *stage_found->second.point;
    const auto & raw = AtomLocalPotentialView::For(*model.FindAtomPtr(atom_id)).GetRawSamplingEntriesView();
    if (raw.size() != peeling.samples.size())
        throw std::invalid_argument("Post-Joint sample/peeling count mismatch.");
    if (!AtomLocalPotentialView::For(*model.FindAtomPtr(atom_id)).HasSampleGeometry())
        throw std::invalid_argument("Post-Joint sample geometry is unavailable.");

    std::unordered_map<std::size_t, std::size_t> rows;
    rows.reserve(input.row_ids.size());
    for (std::size_t row = 0; row < input.row_ids.size(); ++row)
    {
        const auto grid = std::stoull(input.row_ids[row]);
        if (grid >= geometry.GetMapValueArraySize()) throw std::invalid_argument("Joint row outside map geometry.");
        rows.emplace(grid, row);
    }
    std::unordered_map<std::size_t, double> own_charge;
    own_charge.reserve(input.support[atom_index].size());
    for (const auto & support : input.support[atom_index])
    {
        const auto basis = joint_component::EvaluateKernel(support.squared_distance, point.GetWidth(), 2.5);
        own_charge.emplace(support.row, point.GetOffset() * basis.charge);
    }

    PostJointLocalSamples output;
    output.source = stage_found->second.source;
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
        const auto stencil = MakeTricubicStencil(geometry, raw_sample.point.position);
        for (const auto & [grid, weight] : TricubicWeights(stencil))
        {
            if (weight == 0.0) continue;
            if (!rows.contains(grid)) { sample.reason = "outside-joint-domain"; break; }
        }
        if (sample.reason.empty())
        {
            const double charge_response = InterpolateTricubic(stencil, [&](const auto & node) {
                const auto found = rows.find(stencil.Index(node));
                if (found == rows.end()) return 0.0;
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
    return output;
}
}
