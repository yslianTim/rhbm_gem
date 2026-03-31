#include <rhbm_gem/core/command/CommandApi.hpp>

#include "internal/command/HRLModelTestCommand.hpp"
#include "internal/command/MapSimulationCommand.hpp"
#include "internal/command/PotentialAnalysisCommand.hpp"
#include "internal/command/PotentialDisplayCommand.hpp"
#include "internal/command/ResultDumpCommand.hpp"

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
#include "internal/command/MapVisualizationCommand.hpp"
#include "internal/command/PositionEstimationCommand.hpp"
#endif

namespace rhbm_gem {
namespace {

template <typename CommandType, typename RequestType>
CommandResult RunCommand(const RequestType & request)
{
    CommandType command{};
    command.ApplyRequest(request);

    const bool executed{ command.Run() };
    CommandResult result;
    if (executed)
    {
        result.outcome = CommandOutcome::Succeeded;
    }
    else if (command.WasPrepared())
    {
        result.outcome = CommandOutcome::ExecutionFailed;
    }
    else
    {
        result.outcome = CommandOutcome::ValidationFailed;
    }
    result.issues = command.GetValidationIssues();
    return result;
}

} // namespace

const std::vector<CommandInfo> & ListCommands()
{
    static const std::vector<CommandInfo> commands{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
        CommandInfo{ CLI_NAME, DESCRIPTION },
#include "CommandManifest.def"
#undef RHBM_GEM_COMMAND
    };
    return commands;
}

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    CommandResult Run##COMMAND_ID(const COMMAND_ID##Request & request)                         \
    {                                                                                          \
        return RunCommand<COMMAND_ID##Command>(request);                                       \
    }
#include "CommandManifest.def"
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
