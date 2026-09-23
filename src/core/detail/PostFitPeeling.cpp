#include "PostFitPeeling.hpp"
#include "MapInterpolation.hpp"
#include "joint_component/Numerics.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <cmath>
#include <set>
#include <unordered_map>

namespace rhbm_gem::core::detail {
std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem)
{
    const auto & input = problem.Input();
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
    std::vector<bool> has_point;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        const auto view = AtomLocalPotentialView::For(*model.FindAtomPtr(std::stoi(input.atom_ids[atom])));
        const auto & point = view.GetStageEstimate(FittingStage::Second).point;
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
            own[atom].emplace(support.row, value);
        }
    }
    std::map<int, PostFitPeelingResult> results;
    for (std::size_t atom = 0; atom < input.atom_ids.size(); ++atom)
    {
        const auto id = std::stoi(input.atom_ids[atom]);
        const auto view = AtomLocalPotentialView::For(*model.FindAtomPtr(id));
        auto & output = results[id];
        output.source = view.GetStageEstimate(FittingStage::Second).source;
        std::set<std::size_t> neighbors;
        for (const auto & raw : view.GetRawSamplingEntries(false))
        {
            PeelingSampleEstimate sample;
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
}
