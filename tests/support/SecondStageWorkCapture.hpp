#pragma once

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <array>
#include <cstddef>

namespace second_stage_test {

enum class Work { Solver, Operator, Objective, Snapshot };

using WorkCounts = std::array<std::size_t, 4>;

void BeginWorkCapture();
WorkCounts EndWorkCapture();
void CountWork(Work) noexcept;

} // namespace second_stage_test

#define RHBM_TEST_COUNT_WORK(kind) ::second_stage_test::CountWork(::second_stage_test::Work::kind)
#else
#define RHBM_TEST_COUNT_WORK(kind) ((void)0)
#endif
