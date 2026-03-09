#pragma once

#include <algorithm>
#include <optional>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/CommandBase.hpp>

namespace command_test {

inline const rhbm_gem::ValidationIssue * FindValidationIssue(
    const std::vector<rhbm_gem::ValidationIssue> & issues,
    std::string_view option_name,
    std::optional<rhbm_gem::ValidationPhase> phase = std::nullopt,
    std::optional<LogLevel> level = std::nullopt)
{
    const auto iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [&](const rhbm_gem::ValidationIssue & issue)
            {
                return issue.option_name == option_name
                    && (!phase.has_value() || issue.phase == phase.value())
                    && (!level.has_value() || issue.level == level.value());
            })
    };
    return iter == issues.end() ? nullptr : &(*iter);
}

template <typename CommandType>
const rhbm_gem::ValidationIssue * FindValidationIssue(
    const CommandType & command,
    std::string_view option_name,
    std::optional<rhbm_gem::ValidationPhase> phase = std::nullopt,
    std::optional<LogLevel> level = std::nullopt)
{
    return FindValidationIssue(command.GetValidationIssues(), option_name, phase, level);
}

} // namespace command_test
