#pragma once

#include <algorithm>
#include <optional>
#include <string_view>
#include <vector>

#include "command/detail/CommandBase.hpp"

namespace command_test {

inline const rhbm_gem::CommandDiagnostic * FindValidationIssue(
    const std::vector<rhbm_gem::CommandDiagnostic> & issues,
    std::string_view option_name)
{
    const auto iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [&](const rhbm_gem::CommandDiagnostic & issue)
            {
                return issue.option_name == option_name;
            })
    };
    return iter == issues.end() ? nullptr : &(*iter);
}

inline const rhbm_gem::ValidationIssueRecord * FindValidationIssue(
    const std::vector<rhbm_gem::ValidationIssueRecord> & issues,
    std::string_view option_name,
    std::optional<rhbm_gem::ValidationPhase> phase = std::nullopt,
    std::optional<LogLevel> level = std::nullopt)
{
    const auto iter{
        std::find_if(
            issues.begin(),
            issues.end(),
            [&](const rhbm_gem::ValidationIssueRecord & issue)
            {
                return issue.option_name == option_name
                    && (!phase.has_value() || issue.phase == phase.value())
                    && (!level.has_value() || issue.level == level.value());
            })
    };
    return iter == issues.end() ? nullptr : &(*iter);
}

} // namespace command_test
