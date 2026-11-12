#include "core/sampling/MapSampler.hpp"
#include "core/sampling/MapInterpolation.hpp"

#include <array>
#include <cstddef>
#include <exception>
#include <optional>
#include <string>
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

LocalPotentialSampleList BuildAtomInterpolatedSamplesFromPoints(
    const MapObject & map_object,
    const AtomObject & atom,
    SamplingPointList sampling_points)
{
    const auto local_position{ atom.GetPosition() };
    const auto neighbor_atom_list{ atom.FindNeighborAtoms() };
    std::vector<std::array<double, 3>> reject_position_list;
    reject_position_list.reserve(neighbor_atom_list.size());
    for (const auto * neighbor_atom : neighbor_atom_list)
    {
        if (neighbor_atom->GetElement() == Element::HYDROGEN) continue;
        reject_position_list.emplace_back(neighbor_atom->GetPosition());
    }
    sample_filter::FilterSamplingPointList(sampling_points, local_position, reject_position_list);
    return detail::BuildInterpolatedSamplingEntries(map_object, sampling_points);
}

} // namespace

LocalPotentialSampleList BuildInterpolatedMapSamples(
    const MapObject & map_object,
    const GridSampler & sampler,
    const std::array<double, 3> & position,
    const std::array<double, 3> & direction)
{
    const auto sampling_points{ sampler.GenerateSamplingPoints(position, direction) };
    return detail::BuildInterpolatedSamplingEntries(map_object, sampling_points);
}

LocalPotentialSampleList BuildAtomInterpolatedMapSamples(
    const MapObject & map_object,
    const AtomObject & atom)
{
    return BuildAtomInterpolatedSamplesFromPoints(
        map_object,
        atom,
        sphere_sampler::GenerateSamplingPointList(atom.GetPosition()));
}

LocalPotentialSampleList BuildAtomInterpolatedMapSamples(
    const MapObject & map_object,
    const AtomObject & atom,
    SphereSamplingMethod sampling_method)
{
    return BuildAtomInterpolatedSamplesFromPoints(
        map_object,
        atom,
        sphere_sampler::GenerateSamplingPointList(atom.GetPosition(), sampling_method));
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    int thread_count,
    std::string_view progress_label,
    bool quiet_mode)
{
    RunPotentialSamplingWorkflow(map_object, model_object, atom_list, thread_count,
        progress_label, quiet_mode, SphereSamplingMethod::FibonacciDeterministic);
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    int thread_count,
    std::string_view progress_label,
    bool quiet_mode,
    SphereSamplingMethod sampling_method)
{
    std::optional<ScopeTimer> timer;
    if (!quiet_mode) timer.emplace("MapSampler::RunPotentialSamplingWorkflow");
    const std::string progress_message{progress_label};
    size_t atom_count{ 0 };
    std::vector<LocalPotentialSampleList> sampling_entries_by_atom(atom_list.size());
    std::vector<std::exception_ptr> failures(atom_list.size());
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_count)
#endif
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        try { sampling_entries_by_atom[i] = BuildAtomInterpolatedMapSamples(map_object, *atom_list[i], sampling_method); }
        catch (...) { failures[i] = std::current_exception(); }

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            if (!quiet_mode) Logger::ProgressPercent(atom_count, atom_list.size(), 50, progress_message);
        }
    }

    for (const auto & failure : failures) if (failure) std::rethrow_exception(failure);
    auto analysis{ model_object.EditAnalysis() };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        analysis.SetAtomLocalSamplingEntries(
            *atom_list[i], std::move(sampling_entries_by_atom[i]));
    }
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    int thread_count)
{
    RunPotentialSamplingWorkflow(
        map_object,
        model_object,
        model_object.GetSelectedAtoms(),
        thread_count);
}

} // namespace rhbm_gem::core
