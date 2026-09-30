#pragma once

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <array>
#include <vector>

namespace rhbm_gem::core::detail { struct FrozenBackground; }

namespace second_stage_test {

using BackgroundSnapshot = std::vector<std::array<double, 3>>;
using BackgroundSnapshots = std::vector<BackgroundSnapshot>;

void BeginBackgroundCapture();
BackgroundSnapshots EndBackgroundCapture();
void CaptureBackground(const rhbm_gem::core::detail::FrozenBackground &);

} // namespace second_stage_test

#define RHBM_TEST_CAPTURE_BACKGROUND(value) ::second_stage_test::CaptureBackground(value)
#else
#define RHBM_TEST_CAPTURE_BACKGROUND(value) ((void)0)
#endif
