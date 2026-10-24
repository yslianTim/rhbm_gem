#pragma once

#include "core/detail/joint_component/PostFitPeeling.hpp"

#include <stdexcept>
#include <utility>

namespace second_stage_test::postprocessing {
using rhbm_gem::core::detail::JointPostprocessingContext;
using rhbm_gem::core::JointProblem;
using rhbm_gem::MapObject;
using rhbm_gem::ModelObject;
using rhbm_gem::core::detail::PostJointLocalSamples;
using rhbm_gem::JointAnalysisResult;
using rhbm_gem::PostFitPeelingResult;

inline std::map<int, PostFitPeelingResult> BuildPostFitPeelingSamples(
    const JointPostprocessingContext & context,
    std::optional<std::span<const std::size_t>> outputs = std::nullopt)
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
    return BuildPostFitPeelingSamples(
        JointPostprocessingContext(geometry, model, problem, *snapshot), outputs);
}

inline PostJointLocalSamples BuildPostJointLocalSamples(
    const MapObject & geometry, const ModelObject & model, const JointProblem & problem,
    int atom_id, const PostFitPeelingResult & peeling, const JointAnalysisResult & snapshot)
{
    return JointPostprocessingContext(geometry, model, problem, snapshot)
        .BuildPostJointLocalSamples(atom_id, peeling);
}
}
