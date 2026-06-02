#include <rhbm_gem/core/MapSampler.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include <rhbm_gem/data/object/AtomLocalPotentialEditor.hpp>
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

float InterpolateMapValue(const MapObject & data_object, const std::array<float, 3> & position)
{
    auto index{ data_object.GetIndexFromPosition(position) };
    auto origin{ data_object.GetOrigin() };
    auto grid_spacing{ data_object.GetGridSpacing() };

    std::array<float, 3> local;
    for (size_t i = 0; i < 3; i++)
    {
        local.at(i) =
            (position.at(i) - origin.at(i) - static_cast<float>(index.at(i)) * grid_spacing.at(i))
            / grid_spacing.at(i);
    }

    auto cubic_interpolate = [](float p0, float p1, float p2, float p3, float t)
    {
        float a0{ p1 };
        float a1{ 0.5f * (p2 - p0) };
        float a2{ 0.5f * (-p3 + 4.0f * p2 - 5.0f * p1 + 2.0f * p0) };
        float a3{ 0.5f * (p3 - 3.0f * p2 + 3.0f * p1 - p0) };
        return a3 * t * t * t + a2 * t * t + a1 * t + a0;
    };

    std::array<std::array<std::array<float, 4>, 4>, 4> values;
    const auto grid_size{ data_object.GetGridSize() };
    for (int i = -1; i < 3; i++)
    {
        for (int j = -1; j < 3; j++)
        {
            for (int k = -1; k < 3; k++)
            {
                size_t i_next{ static_cast<size_t>(i + 1) };
                size_t j_next{ static_cast<size_t>(j + 1) };
                size_t k_next{ static_cast<size_t>(k + 1) };
                int xi{ std::clamp(index.at(0) + i, 0, grid_size.at(0) - 1) };
                int yj{ std::clamp(index.at(1) + j, 0, grid_size.at(1) - 1) };
                int zk{ std::clamp(index.at(2) + k, 0, grid_size.at(2) - 1) };
                values[i_next][j_next][k_next] = data_object.GetMapValue(xi, yj, zk);
            }
        }
    }

    std::array<std::array<float, 4>, 4> temp_y;
    for (size_t j = 0; j < 4; j++)
    {
        for (size_t k = 0; k < 4; k++)
        {
            temp_y[j][k] = cubic_interpolate(
                values[0][j][k], values[1][j][k], values[2][j][k], values[3][j][k], local.at(0));
        }
    }

    std::array<float, 4> temp_z;
    for (size_t k = 0; k < 4; k++)
    {
        temp_z[k] = cubic_interpolate(
            temp_y[0][k], temp_y[1][k], temp_y[2][k], temp_y[3][k], local.at(1));
    }

    return cubic_interpolate(temp_z[0], temp_z[1], temp_z[2], temp_z[3], local.at(2));
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

std::vector<AtomLocalPotentialEditor> BuildSelectedAtomLocalEditors(ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }
    return local_editor_list;
}

} // namespace

LocalPotentialSampleList SampleMapValues(
    const MapObject & map_object,
    const GridSampler & sampler,
    const std::array<float, 3> & position,
    const std::array<float, 3> & direction)
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
    std::vector<std::array<float, 3>> reject_position_list;
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
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_count)
#endif
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        auto sampling_entries{
            SampleAtomMapValues(map_object, *atom_list[i], sampling_method)
        };
        local_editor_list[i].SetSamplingEntries(std::move(sampling_entries));

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_list.size());
        }
    }
}

} // namespace rhbm_gem::core
