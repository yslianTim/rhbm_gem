#pragma once

#include "core/detail/CandidateSelection.hpp"

namespace rhbm_gem::core::detail {

void ReconcileSelectedBoundaries(
    const CandidateSelectionInputs & inputs,
    const std::map<ClusterKey, FitStatePatch> & rescue_patch_by_key,
    CandidateSelection & selection);

void ReauditFallbackSelection(const CandidateSelectionInputs & inputs, CandidateSelection & selection);

} // namespace rhbm_gem::core::detail
