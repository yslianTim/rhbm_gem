#include <rhbm_gem/core/MapSampler.hpp>
#include "core/detail/MapSampler.hpp"
#include "core/detail/MapInterpolation.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/math/GridSampler.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

namespace rhbm_gem::core::detail {
LocalPotentialSampleList BuildLocalPotentialSampleList(
    const MapObject & map_object,
    const SamplingPointList & sample_point_list)
{
    LocalPotentialSampleList sampling_data_list;
    sampling_data_list.reserve(sample_point_list.size());
    for (const auto & sampling_point : sample_point_list)
    {
        const auto stencil{MakeTricubicStencil(map_object, sampling_point.position)};
        const auto map_value{InterpolateTricubic(stencil,
            [&](const auto & node) { return map_object.GetMapValue(node[0], node[1], node[2]); })};
        sampling_data_list.emplace_back(LocalPotentialSample{
            map_value,
            sampling_point
        });
    }
    return sampling_data_list;
}
} // namespace rhbm_gem::core::detail

namespace rhbm_gem::core {

LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const GridSampler & sampler,
    const std::array<double, 3> & position,
    const std::array<double, 3> & direction)
{
    const auto sample_point_list{ sampler.GenerateSamplingPoints(position, direction) };
    return detail::BuildLocalPotentialSampleList(map_object, sample_point_list);
}

LocalPotentialSampleList SampleAtomMapValues(
    const MapObject & map_object,
    const AtomObject & atom,
    SphereSamplingMethod sampling_method)
{
    const auto local_position{ atom.GetPosition() };
    auto sample_point_list{
        sphere_sampler::GenerateSamplingPointList(local_position, sampling_method)
    };
    const auto neighbor_atom_list{ atom.FindNeighborAtoms() };
    std::vector<std::array<double, 3>> reject_position_list;
    reject_position_list.reserve(neighbor_atom_list.size());
    for (const auto * neighbor_atom : neighbor_atom_list)
    {
        reject_position_list.emplace_back(neighbor_atom->GetPosition());
    }
    sample_filter::FilterSamplingPointList(sample_point_list, local_position, reject_position_list);
    return detail::BuildLocalPotentialSampleList(map_object, sample_point_list);
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    SphereSamplingMethod sampling_method,
    int thread_count)
{
    ScopeTimer timer("MapSampler::RunPotentialSamplingWorkflow");
    size_t atom_count{ 0 };
    std::vector<LocalPotentialSampleList> raw_sampling_entries_list(atom_list.size());
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_count)
#endif
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        raw_sampling_entries_list[i] =
            SampleAtomMapValues(map_object, *atom_list[i], sampling_method);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_list.size());
        }
    }

    auto analysis{ model_object.EditAnalysis() };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        analysis.SetAtomLocalRawSamplingEntries(
            *atom_list[i], std::move(raw_sampling_entries_list[i]));
    }
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    SphereSamplingMethod sampling_method,
    int thread_count)
{
    RunPotentialSamplingWorkflow(
        map_object,
        model_object,
        model_object.GetSelectedAtoms(),
        sampling_method,
        thread_count);
}

} // namespace rhbm_gem::core
