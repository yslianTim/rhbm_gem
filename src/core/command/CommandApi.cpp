#include <rhbm_gem/core/command/CommandApi.hpp>

#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/command/MapSimulationCommand.hpp>
#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>
#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/command/PotentialDisplayCommand.hpp>
#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

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

} // namespace rhbm_gem
