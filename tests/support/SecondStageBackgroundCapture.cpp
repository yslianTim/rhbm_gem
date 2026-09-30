#include "support/SecondStageBackgroundCapture.hpp"

#include "core/detail/second_stage/SecondStageState.hpp"

#include <utility>

namespace second_stage_test {
namespace {
bool enabled{};
BackgroundSnapshots backgrounds;
}

void BeginBackgroundCapture()
{
    enabled = false;
    backgrounds.clear();
    enabled = true;
}

BackgroundSnapshots EndBackgroundCapture()
{
    enabled = false;
    return std::move(backgrounds);
}

void CaptureBackground(const rhbm_gem::core::detail::FrozenBackground & background)
{
    if (!enabled) return;
    BackgroundSnapshot models;
    models.reserve(background.model_by_atom.size());
    for (const auto & model : background.model_by_atom)
    {
        models.push_back({ model.GetAmplitude(), model.GetWidth(), model.GetOffset() });
    }
    backgrounds.push_back(std::move(models));
}
}
