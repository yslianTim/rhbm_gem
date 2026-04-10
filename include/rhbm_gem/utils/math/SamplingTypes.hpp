#pragma once

#include <array>
#include <tuple>
#include <vector>

using SamplingPoint = std::tuple<float, std::array<float, 3>>;
using SamplingPointList = std::vector<SamplingPoint>;
