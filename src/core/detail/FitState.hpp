#pragma once

#include <cstddef>
#include <vector>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem::core::detail {

using FitState = std::vector<LocalGaussianResult>;

inline const GaussianModel3D & GetFitModel(
    const FitState & state,
    std::size_t atom_index)
{
    return state.at(atom_index).mdpde.GetModel();
}

} // namespace rhbm_gem::core::detail
