#pragma once

#include <vector>

namespace rhbm_gem::algorithm {

template <typename ResultType, typename EstimationType>
struct IterationState
{
    std::vector<ResultType> result_list;
    std::vector<EstimationType> estimation_list;
};

} // namespace rhbm_gem::algorithm
