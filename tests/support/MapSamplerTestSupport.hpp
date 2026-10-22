#pragma once

#include <functional>

namespace rhbm_gem::core::detail {
using PotentialSamplingObserver = std::function<void(int)>;
PotentialSamplingObserver & PotentialSamplingObserverForTesting();
}
