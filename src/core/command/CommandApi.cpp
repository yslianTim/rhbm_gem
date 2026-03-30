#include <rhbm_gem/core/command/CommandApi.hpp>

#include "internal/command/CommandRegistry.hpp"

namespace rhbm_gem {
namespace {

template <typename CommandType, typename RequestType>
ExecutionReport RunCommand(const RequestType & request)
{
    CommandType command{};
    command.ApplyRequest(request);

    ExecutionReport report;
    report.prepared = command.PrepareForExecution();
    report.validation_issues = command.GetValidationIssues();
    if (!report.prepared)
    {
        return report;
    }

    report.executed = command.Execute();
    const auto & post_execution_issues{ command.GetValidationIssues() };
    if (!post_execution_issues.empty())
    {
        report.validation_issues = post_execution_issues;
    }
    return report;
}

} // namespace

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    ExecutionReport Run##COMMAND_ID(const COMMAND_ID##Request & request)                       \
    {                                                                                          \
        return RunCommand<COMMAND_ID##Command>(request);                                       \
    }
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
