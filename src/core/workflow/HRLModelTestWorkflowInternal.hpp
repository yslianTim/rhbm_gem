#pragma once

#include <filesystem>
#include <functional>
#include <string_view>

#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>

namespace rhbm_gem::detail {

struct HRLModelTestWorkflowContext
{
    const HRLModelTestCommandOptions & options;
    std::function<std::filesystem::path(std::string_view)> build_output_path;
};

bool RunHRLModelTestWorkflow(const HRLModelTestWorkflowContext & context);

} // namespace rhbm_gem::detail
