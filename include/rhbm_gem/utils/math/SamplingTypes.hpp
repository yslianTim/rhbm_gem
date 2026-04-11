#pragma once

#include <array>
#include <optional>
#include <vector>

struct SamplingPoint
{
    float distance{ 0.0f };
    std::array<float, 3> position{ 0.0f, 0.0f, 0.0f };
};

using SamplingPointList = std::vector<SamplingPoint>;

struct LocalPotentialSample
{
    float distance{ 0.0f };
    float response{ 0.0f };
    float weight{ 1.0f };
    std::optional<std::array<float, 3>> position{};

    bool HasPosition() const noexcept { return position.has_value(); }
};

using LocalPotentialSampleList = std::vector<LocalPotentialSample>;

struct SeriesPoint
{
    float x{ 0.0f };
    float y{ 0.0f };
    float weight{ 1.0f };
};

using SeriesPointList = std::vector<SeriesPoint>;
