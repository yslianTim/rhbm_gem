#pragma once

#include <cstddef>

namespace rhbm_gem::gui::internal {

// Exposes static wrapper pooling counters to tests without leaking into public API.
std::size_t DefaultServiceBuildCountForTesting();
std::size_t DefaultExecutorBuildCountForTesting();

} // namespace rhbm_gem::gui::internal
