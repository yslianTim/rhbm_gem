#pragma once

#include <string>

namespace rhbm_gem {

struct StyleSpec
{
    short color{ 1 };
    short solid_marker{ 1 };
    short open_marker{ 1 };
    short line_style{ 1 };
    std::string label{ "UNK" };
};

} // namespace rhbm_gem
