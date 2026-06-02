#pragma once

#include <array>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

class GridSampler;

namespace rhbm_gem {

class AtomObject;
class MapObject;
class ModelObject;

namespace core {

LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const GridSampler & sampler,
    const std::array<float, 3> & position,
    const std::array<float, 3> & direction);

LocalPotentialSampleList SampleAtomMapValues(
    const MapObject & map_object,
    const AtomObject & atom,
    SphereSamplingMethod sampling_method);

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    SphereSamplingMethod sampling_method,
    int thread_count);

} // namespace core

} // namespace rhbm_gem
