#include <rhbm_gem/core/MapSampler.hpp>
#include "core/detail/MapInterpolation.hpp"
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include "support/ForwardModelExperiment.hpp"
#endif

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

namespace rhbm_gem::core {
namespace {

double InterpolateMapValue(const MapObject & data_object, const std::array<double, 3> & position)
{
    return detail::InterpolateTricubic(detail::MakeTricubicStencil(data_object, position),
        [&](const auto & node) { return data_object.GetMapValue(node[0], node[1], node[2]); });
}

LocalPotentialSampleList BuildLocalPotentialSampleList(
    const MapObject & map_object,
    const SamplingPointList & sample_point_list)
{
    LocalPotentialSampleList sampling_data_list;
    sampling_data_list.reserve(sample_point_list.size());
    for (const auto & sampling_point : sample_point_list)
    {
        auto map_value{
            InterpolateMapValue(map_object, sampling_point.position)
        };
        sampling_data_list.emplace_back(LocalPotentialSample{
            map_value,
            sampling_point
        });
    }
    return sampling_data_list;
}

} // namespace

LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const GridSampler & sampler,
    const std::array<double, 3> & position,
    const std::array<double, 3> & direction)
{
    const auto sample_point_list{ sampler.GenerateSamplingPoints(position, direction) };
    return BuildLocalPotentialSampleList(map_object, sample_point_list);
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
    return BuildLocalPotentialSampleList(map_object, sample_point_list);
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    SphereSamplingMethod sampling_method,
    int thread_count)
{
    ScopeTimer timer("MapSampler::RunPotentialSamplingWorkflow");
    const auto & atom_list{ model_object.GetSelectedAtoms() };
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

} // namespace rhbm_gem::core

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
LocalPotentialSampleList second_stage_test::SampleExperimentPoints(
    const rhbm_gem::MapObject & map, const SamplingPointList & points)
{
    return rhbm_gem::core::BuildLocalPotentialSampleList(map, points);
}
#endif
