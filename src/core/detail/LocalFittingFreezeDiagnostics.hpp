#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

#include <rhbm_gem/utils/algorithm/ParameterChange.hpp>

namespace rhbm_gem::core::detail {

enum class LocalFittingFreezeOutcome
{
    Ineligible,
    AboveThreshold,
    Stabilizing,
    NewlyFrozen
};

struct LocalFittingFreezeEvidenceDiagnostic
{
    LocalFittingFreezeOutcome outcome{ LocalFittingFreezeOutcome::Ineligible };
    std::size_t dominant_parameter_index{ 0 };
    double maximum_evidence{ 0.0 };
};

struct LocalFittingSelfPeerBlockerDiagnostic
{
    bool self{ false };
    bool peer{ false };
};

inline LocalFittingSelfPeerBlockerDiagnostic ClassifyLocalFittingSelfPeerBlocker(
    std::size_t atom_index,
    const std::vector<std::size_t> & blocker_atom_index_list)
{
    const auto self{
        std::find(
            blocker_atom_index_list.begin(),
            blocker_atom_index_list.end(),
            atom_index) != blocker_atom_index_list.end()
    };
    const auto peer{
        std::any_of(
            blocker_atom_index_list.begin(),
            blocker_atom_index_list.end(),
            [&](std::size_t blocker_atom_index)
            {
                return blocker_atom_index != atom_index;
            })
    };
    return LocalFittingSelfPeerBlockerDiagnostic{ self, peer };
}

inline LocalFittingFreezeEvidenceDiagnostic ClassifyLocalFittingFreezeEvidence(
    bool stationarity_eligible,
    bool frozen_after_update,
    double freeze_threshold,
    const algorithm::ParameterChange & evidence)
{
    if (evidence.value_list.empty())
    {
        throw std::invalid_argument("Local fitting freeze evidence must not be empty.");
    }
    if (!std::isfinite(freeze_threshold) || freeze_threshold < 0.0)
    {
        throw std::invalid_argument("Local fitting freeze threshold must be finite and non-negative.");
    }
    const auto maximum_iter{
        std::max_element(evidence.value_list.begin(), evidence.value_list.end())
    };
    const auto dominant_parameter_index{
        static_cast<std::size_t>(std::distance(evidence.value_list.begin(), maximum_iter))
    };
    const auto maximum_evidence{ *maximum_iter };

    LocalFittingFreezeOutcome outcome{ LocalFittingFreezeOutcome::Ineligible };
    if (stationarity_eligible)
    {
        if (frozen_after_update)
        {
            outcome = LocalFittingFreezeOutcome::NewlyFrozen;
        }
        else if (maximum_evidence >= freeze_threshold)
        {
            outcome = LocalFittingFreezeOutcome::AboveThreshold;
        }
        else
        {
            outcome = LocalFittingFreezeOutcome::Stabilizing;
        }
    }
    return LocalFittingFreezeEvidenceDiagnostic{
        outcome,
        dominant_parameter_index,
        maximum_evidence
    };
}

} // namespace rhbm_gem::core::detail
