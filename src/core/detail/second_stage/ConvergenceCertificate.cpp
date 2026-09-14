#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr double kConvergencePercentile{ 0.99 };

static void ValidateBlockActivitySize(
    std::size_t atom_count,
    const SuspiciousBlockActivity & block_activity)
{
    if (block_activity.shape_fixed_atom_mask.size() != atom_count ||
        block_activity.offset_fixed_atom_mask.size() != atom_count ||
        block_activity.hard_failure_atom_mask.size() != atom_count)
    {
        throw std::invalid_argument("Active-coordinate convergence inputs are inconsistent.");
    }
}

} // namespace

bool ConvergenceCertificate::StrictOperatorPassed() const
{
    return solver_qualified && operator_complete &&
        IsTransformedPercentileConverged(operator_nominal_p99);
}

bool ConvergenceCertificate::ProductionConverged() const
{
    return StrictOperatorPassed() &&
        IsTransformedPercentileConverged(accepted_active_p99) &&
        !objective_domain_changed && !quarantine_transition &&
        !suspicious_block_fallback && !rejected_cluster;
}

ActiveCoordinatePopulation BuildActiveCoordinatePopulation(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousBlockActivity & block_activity)
{
    ValidateBlockActivitySize(block_activity.shape_fixed_atom_mask.size(), block_activity);
    ActiveCoordinatePopulation result;
    for (const auto atom_index : atom_index_list)
    {
        if (block_activity.HasActiveShape(atom_index))
        {
            result.active_shape_atom_index_list.emplace_back(atom_index);
        }
        if (block_activity.HasActiveOffset(atom_index))
        {
            result.active_offset_atom_index_list.emplace_back(atom_index);
        }
    }
    return result;
}

TransformedChangeSummary SummarizeActiveDofChanges(
    const std::vector<TransformedChange> & change_list,
    const ActiveCoordinatePopulation & population)
{
    std::array<std::vector<double>, GaussianModel3D::TransformedCoordinateSize()> parameter_change_lists;
    for (const auto parameter_index : std::array<std::size_t, 2>{
        GaussianModel3D::LogPeakHeightCoordinateIndex(),
        GaussianModel3D::LogWidthCoordinateIndex() })
    {
        auto & values{ parameter_change_lists.at(parameter_index) };
        values.reserve(population.active_shape_atom_index_list.size());
        for (const auto atom_index : population.active_shape_atom_index_list)
        {
            if (atom_index >= change_list.size())
            {
                throw std::invalid_argument("Active-coordinate shape change input is inconsistent.");
            }
            values.emplace_back(change_list.at(atom_index).at(parameter_index));
        }
    }

    auto & offset_change_list{
        parameter_change_lists.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex())
    };
    offset_change_list.reserve(population.active_offset_atom_index_list.size());
    for (const auto atom_index : population.active_offset_atom_index_list)
    {
        if (atom_index >= change_list.size())
        {
            throw std::invalid_argument("Active-coordinate offset change input is inconsistent.");
        }
        const auto value{
            change_list.at(atom_index).at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex())
        };
        offset_change_list.emplace_back(
            std::isfinite(value) ? std::abs(value) : std::numeric_limits<double>::infinity());
    }

    TransformedChangeSummary result;
    for (std::size_t parameter_index = 0; parameter_index < parameter_change_lists.size(); parameter_index++)
    {
        const auto & values{ parameter_change_lists.at(parameter_index) };
        result.population_size_list.at(parameter_index) = values.size();
        result.percentile_list.at(parameter_index) = array_helper::ComputePercentile(values, kConvergencePercentile);
        result.maximum_list.at(parameter_index) = values.empty() ? 0.0 : *std::ranges::max_element(values);
    }
    return result;
}

TransformedChangeSummary SummarizeActiveDofChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const ActiveCoordinatePopulation & population)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument("Active-coordinate transformed state sizes are inconsistent.");
    }
    std::vector<TransformedChange> change_list;
    change_list.reserve(current_state.size());
    for (std::size_t atom_index = 0; atom_index < current_state.size(); atom_index++)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, atom_index),
            GetFitModel(previous_state, atom_index)));
    }
    return SummarizeActiveDofChanges(change_list, population);
}

FixedPointOperatorSummary SummarizeFixedPointOperator(
    const FixedPointOperatorEvidence & evidence,
    const FitState & previous_state,
    const std::vector<std::size_t> & atom_index_list)
{
    FixedPointOperatorSummary result;
    const ActiveCoordinatePopulation operator_nominal_population{
        atom_index_list,
        atom_index_list
    };

    auto & change_list{ result.change_list };
    change_list.reserve(previous_state.size());
    for (std::size_t atom_index = 0; atom_index < previous_state.size(); atom_index++)
    {
        auto change{ CalculateTransformedChange(
            GetFitModel(evidence.state, atom_index),
            GetFitModel(previous_state, atom_index)) };
        if (evidence.shape_available_atom_mask.at(atom_index) == 0)
        {
            change.at(GaussianModel3D::LogPeakHeightCoordinateIndex()) = std::numeric_limits<double>::infinity();
            change.at(GaussianModel3D::LogWidthCoordinateIndex()) = std::numeric_limits<double>::infinity();
        }
        if (evidence.offset_available_atom_mask.at(atom_index) == 0)
        {
            change.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) = std::numeric_limits<double>::infinity();
        }
        change_list.emplace_back(std::move(change));
    }
    result.nominal_residual = SummarizeActiveDofChanges(change_list, operator_nominal_population);
    result.operator_complete = std::ranges::all_of(
        atom_index_list,
        [&](const auto atom_index)
        {
            return evidence.shape_available_atom_mask.at(atom_index) != 0 &&
                evidence.offset_available_atom_mask.at(atom_index) != 0;
        });
    return result;
}

bool AreActiveCoordinatesSolverQualified(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key)
{
    const auto atom_count{ block_activity.shape_fixed_atom_mask.size() };
    ValidateBlockActivitySize(atom_count, block_activity);
    if (local_refit_status_by_atom.size() != atom_count)
    {
        throw std::invalid_argument("Convergence audit qualification inputs are inconsistent.");
    }

    for (const auto atom_index : atom_index_list)
    {
        if (!block_activity.HasActiveShape(atom_index)) continue;
        const auto & status{ local_refit_status_by_atom[atom_index] };
        if (!status.has_value() || !IsLocalRefitStatusSolverQualified(*status))
        {
            return false;
        }
    }

    for (const auto atom_index : atom_index_list)
    {
        if (!block_activity.HasActiveOffset(atom_index)) continue;
        const auto owner{ std::ranges::find_if(cluster_key_list, [&](const auto & key)
        {
            return std::ranges::binary_search(key, atom_index);
        }) };
        if (owner == cluster_key_list.end()) return false;
        const auto health_iter{ health_by_key.find(*owner) };
        if (health_iter == health_by_key.end() ||
            health_iter->second.joint_offset_status != JointOffsetSolveStatus::Converged)
        {
            return false;
        }
    }

    return true;
}

} // namespace rhbm_gem::core::detail
