#include "support/SecondStageWorkCapture.hpp"

#include "support/EndpointRefinementTestSupport.hpp"

#include <atomic>

namespace second_stage_test {
namespace {
std::atomic<bool> enabled{ false };
std::array<std::atomic<std::size_t>, 4> counts{};
}

void BeginWorkCapture()
{
    enabled = false;
    for (auto & count : counts) count = 0;
    enabled = true;
}

WorkCounts EndWorkCapture()
{
    enabled = false;
    WorkCounts result{};
    for (std::size_t i = 0; i < counts.size(); ++i) result[i] = counts[i].load();
    return result;
}

void CountWork(Work work) noexcept
{
    if (enabled.load(std::memory_order_relaxed) && !IsEndpointOperatorProbe())
    {
        counts[static_cast<std::size_t>(work)].fetch_add(1, std::memory_order_relaxed);
    }
}
}
