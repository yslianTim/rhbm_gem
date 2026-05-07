#pragma once

#include <algorithm>
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

} // namespace command_test
