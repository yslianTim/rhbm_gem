#include <rhbm_gem/core/MapSampler.hpp>
#include "core/detail/MapInterpolation.hpp"
#include "core/detail/JointSeedSampling.hpp"
#include "core/detail/gaussian_fit/FittingRanges.hpp"
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include "support/MapSamplerTestSupport.hpp"
#endif

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
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

namespace rhbm_gem::core::detail {
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
PotentialSamplingObserver & PotentialSamplingObserverForTesting()
{
    static thread_local PotentialSamplingObserver observer;
    return observer;
}
#endif

} // namespace rhbm_gem::core::detail

namespace rhbm_gem::core {
namespace {

LocalPotentialSampleList SampleAtomMapValuesFromPoints(
    const MapObject & map_object,
    const AtomObject & atom,
    SamplingPointList sample_point_list,
    detail::MapInterpolationMethod interpolation = detail::MapInterpolationMethod::Tricubic)
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
    sample_filter::FilterSamplingPointList(sample_point_list, local_position, reject_position_list);
    return detail::BuildLocalPotentialSampleList(map_object, sample_point_list, interpolation);
}

template<class SampleAtom>
void RunPotentialSamplingWorkflowImpl(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    int thread_count,
    std::string_view progress_label,
    bool quiet_mode,
    SampleAtom sample_atom)
{
    std::optional<ScopeTimer> timer;
    if (!quiet_mode) timer.emplace("MapSampler::RunPotentialSamplingWorkflow");
    const std::string progress_message{progress_label};
    size_t atom_count{ 0 };
    std::vector<LocalPotentialSampleList> raw_sampling_entries_list(atom_list.size());
    std::vector<std::exception_ptr> failures(atom_list.size());
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_count)
#endif
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        try { raw_sampling_entries_list[i] = sample_atom(i); }
        catch (...) { failures[i] = std::current_exception(); }
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
        if (detail::PotentialSamplingObserverForTesting())
            detail::PotentialSamplingObserverForTesting()(atom_list[i]->GetSerialID());
#endif

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
        analysis.SetAtomLocalRawSamplingEntries(
            *atom_list[i], std::move(raw_sampling_entries_list[i]));
    }
}

} // namespace

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
    const AtomObject & atom)
{
    return SampleAtomMapValuesFromPoints(
        map_object,
        atom,
        sphere_sampler::GenerateSamplingPointList(atom.GetPosition()));
}

LocalPotentialSampleList SampleAtomMapValues(
    const MapObject & map_object,
    const AtomObject & atom,
    SphereSamplingMethod sampling_method)
{
    return SampleAtomMapValuesFromPoints(
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
    RunPotentialSamplingWorkflowImpl(
        model_object, atom_list, thread_count, progress_label, quiet_mode,
        [&](std::size_t i) { return SampleAtomMapValues(map_object, *atom_list[i]); });
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

std::string_view detail::JointSeedSamplingName(JointSeedSamplingMethod method)
{
    switch (method)
    {
        case JointSeedSamplingMethod::FibonacciShells: return "FibonacciShells";
        case JointSeedSamplingMethod::QScoreSpiralShells: return "QScoreSpiralShells";
    }
    throw std::invalid_argument("Unsupported Joint seed sampling method.");
}

std::string_view detail::MapInterpolationName(MapInterpolationMethod method)
{
    switch (method)
    {
        case MapInterpolationMethod::Tricubic: return "Tricubic";
        case MapInterpolationMethod::Trilinear: return "Trilinear";
    }
    throw std::invalid_argument("Unsupported map interpolation method.");
}

void detail::RunJointSeedSamplingWorkflow(
    MapObject & map, ModelObject & model, const std::vector<AtomObject *> & atoms,
    int thread_count, std::string_view progress_label, bool quiet_mode,
    const JointSeedSamplingPolicy & policy, JointSeedSamplingReport * report)
{
    // Validate before entering OpenMP and snapshot the policy for this run.
    JointSeedSamplingName(policy.sampling);
    MapInterpolationName(policy.interpolation);
    const auto run_policy = policy;
    if (report) { report->policy = run_policy; report->atoms.assign(atoms.size(), {}); }
    RunPotentialSamplingWorkflowImpl(model, atoms, std::clamp(thread_count, 1, 4),
        progress_label, quiet_mode, [&](std::size_t index) {
            const auto & atom = *atoms[index];
            LocalPotentialSampleList samples;
            if (run_policy.sampling == JointSeedSamplingMethod::FibonacciShells)
                samples = SampleAtomMapValuesFromPoints(map, atom,
                    sphere_sampler::GenerateSamplingPointList(atom.GetPosition()), run_policy.interpolation);
            else
            {
                auto seed = BuildQScoreSeedSamples(atom, model);
                samples = BuildLocalPotentialSampleList(map, seed.points, run_policy.interpolation);
                if (report) report->atoms[index].shells = std::move(seed.shells);
            }
            if (report)
            {
                auto & record = report->atoms[index];
                record.atom_id = atom.GetSerialID();
                record.raw_count = samples.size();
                for (const auto & sample : samples)
                {
                    if (!sample.point.is_selected) continue;
                    ++record.selected_count;
                    const bool positive = std::isfinite(sample.response) && sample.response > 0.;
                    if (positive) ++record.positive_count;
                    if (IsSignalDistance(sample.point.distance))
                    {
                        ++record.fitting_count;
                        if (positive) ++record.positive_fitting_count;
                    }
                }
            }
            return samples;
        });
}

} // namespace rhbm_gem::core
