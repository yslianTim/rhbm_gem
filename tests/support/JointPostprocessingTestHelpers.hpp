#pragma once

#include "core/detail/joint_component/PostFitPeeling.hpp"

#include <stdexcept>
#include <utility>

namespace second_stage_test::postprocessing {
using rhbm_gem::core::detail::JointPostprocessingContext;
using rhbm_gem::core::detail::JointPostprocessingMetrics;
using rhbm_gem::core::JointProblem;
using rhbm_gem::MapObject;
using rhbm_gem::ModelObject;
using rhbm_gem::core::detail::PostJointLocalSamples;
using rhbm_gem::JointAnalysisResult;
using rhbm_gem::PostFitPeelingResult;

// Small abstract grid evidence for editor/eligibility tests that do not own a map.
inline void AssignGridObservationsForTesting(PostFitPeelingResult & peeling,
    const LocalPotentialSampleList & raw, const std::array<double,3> & center = {})
{
    if (peeling.samples.size()!=raw.size()) throw std::invalid_argument("Test sample count mismatch.");
    for (std::size_t i=0;i<raw.size();++i)
    {
        auto position=center;
        position[0]+=raw[i].point.distance;
        peeling.samples[i].observation={i,i,raw[i].point.distance*raw[i].point.distance,raw[i].response,position};
    }
}

inline JointPostprocessingContext MakeJointPostprocessingContext(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    const JointAnalysisResult & snapshot, JointPostprocessingMetrics * metrics = nullptr,
    const rhbm_gem::core::detail::HaloMedianBuildResult * halo_medians = nullptr)
{
    return JointPostprocessingContext(geometry,model,problem,snapshot,metrics,halo_medians,
        rhbm_gem::core::detail::BuildJointPeelingContributorUniverse(model,problem,snapshot));
}

inline std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const JointPostprocessingContext & context,
    std::span<const std::size_t> outputs)
{
    std::map<int, PostFitPeelingResult> result;
    for (auto & [id, peeling] : context.BuildPostFitPeelingBatch(outputs))
        result.emplace(id, std::move(peeling));
    return result;
}

inline std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    std::optional<std::span<const std::size_t>> outputs = std::nullopt,
    const JointAnalysisResult * snapshot = nullptr)
{
    if (!snapshot)
        throw std::invalid_argument("Pure Exclusion peeling requires a Pure Exclusion snapshot.");
    auto context=MakeJointPostprocessingContext(geometry,model,problem,*snapshot);
    const auto & selection=problem.OriginalInput().selection_domain;
    if (!outputs && !selection)
        throw std::invalid_argument("Pure Exclusion peeling requires a target selection.");
    const std::span<const std::size_t> selected=outputs ? *outputs :
        std::span<const std::size_t>(selection->target_indices);
    return BuildPostFitPeelingSamples(context,selected);
}

inline PostJointLocalSamples BuildPostJointLocalSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    int atom_id, const PostFitPeelingResult & peeling, const JointAnalysisResult & snapshot)
{
    return MakeJointPostprocessingContext(geometry,model,problem,snapshot)
        .BuildPostJointLocalSamples(atom_id,peeling);
}
}
