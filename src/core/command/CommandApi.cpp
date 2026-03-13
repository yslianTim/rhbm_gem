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

// BEGIN GENERATED: command-run-definitions
ExecutionReport RunPotentialAnalysis(const PotentialAnalysisRequest & request)
{
    return RunCommand<PotentialAnalysisCommand>(request);
}

ExecutionReport RunPotentialDisplay(const PotentialDisplayRequest & request)
{
    return RunCommand<PotentialDisplayCommand>(request);
}

ExecutionReport RunResultDump(const ResultDumpRequest & request)
{
    return RunCommand<ResultDumpCommand>(request);
}

ExecutionReport RunMapSimulation(const MapSimulationRequest & request)
{
    return RunCommand<MapSimulationCommand>(request);
}

ExecutionReport RunMapVisualization(const MapVisualizationRequest & request)
{
    return RunCommand<MapVisualizationCommand>(request);
}

ExecutionReport RunPositionEstimation(const PositionEstimationRequest & request)
{
    return RunCommand<PositionEstimationCommand>(request);
}

ExecutionReport RunHRLModelTest(const HRLModelTestRequest & request)
{
    return RunCommand<HRLModelTestCommand>(request);
}
// END GENERATED: command-run-definitions

} // namespace rhbm_gem
