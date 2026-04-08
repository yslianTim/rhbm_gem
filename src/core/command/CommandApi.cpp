#include <rhbm_gem/core/command/CommandApi.hpp>

#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#endif

namespace rhbm_gem {
namespace {

std::vector<ValidationIssue> BuildPublicValidationIssues(
    const std::vector<ValidationIssueRecord> & internal_issues)
{
    std::vector<ValidationIssue> public_issues;
    public_issues.reserve(internal_issues.size());
    for (const auto & issue : internal_issues)
    {
        public_issues.push_back(ValidationIssue{
            issue.option_name,
            issue.message
        });
    }
    return public_issues;
}

template <typename CommandType, typename RequestType>
CommandResult RunCommand(const RequestType & request)
{
    CommandType command{};
    command.ApplyRequest(request);

    CommandResult result;
    result.succeeded = command.Run();
    result.issues = BuildPublicValidationIssues(command.GetValidationIssues());
    return result;
}

} // namespace

const std::vector<CommandInfo> & ListCommands()
{
    static const std::vector<CommandInfo> commands{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
        CommandInfo{ CLI_NAME, DESCRIPTION },
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND
    };
    return commands;
}

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    CommandResult Run##COMMAND_ID(const COMMAND_ID##Request & request)                         \
    {                                                                                          \
        return RunCommand<COMMAND_ID##Command>(request);                                       \
    }
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
