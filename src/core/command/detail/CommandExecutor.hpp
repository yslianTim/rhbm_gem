#pragma once

#include "CommandBase.hpp"

#include <vector>

namespace rhbm_gem::command_internal {

template <typename CommandType>
CommandResult ExecuteCommandInstance(const typename CommandType::RequestType & request)
{
    CommandType command{};
    command.ApplyRequest(request);

    CommandResult result;
    result.succeeded = command.Run();

    const auto & internal_issues{ command.GetValidationIssues() };
    result.issues.reserve(internal_issues.size());
    for (const auto & issue : internal_issues)
    {
        result.issues.push_back(CommandDiagnostic{ issue.option_name, issue.message });
    }
    return result;
}

} // namespace rhbm_gem::command_internal
