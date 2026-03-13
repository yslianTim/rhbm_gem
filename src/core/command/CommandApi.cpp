#include <rhbm_gem/core/command/CommandApi.hpp>

#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

namespace rhbm_gem {
namespace {

template <typename CommandType, typename RequestType>
ExecutionReport RunCommand(CommonOptionProfile profile, const RequestType & request)
{
    CommandType command{ profile };
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

#define RHBM_GEM_COMMAND(COMMAND_ID, COMMAND_STEM, CLI_NAME, DESCRIPTION, PROFILE)             \
    ExecutionReport Run##COMMAND_STEM(const COMMAND_STEM##Request & request)                    \
    {                                                                                            \
        return RunCommand<COMMAND_STEM##Command>(CommonOptionProfile::PROFILE, request);        \
    }
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
