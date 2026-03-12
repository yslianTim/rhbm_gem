#include <rhbm_gem/gui/GuiCommandExecutor.hpp>

#include <rhbm_gem/core/command/MapSimulationCommand.hpp>
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

namespace rhbm_gem::gui {

namespace {

template <typename CommandType>
void ApplyCommonOptions(
    CommandType & command,
    const CommonExecutionRequest & common)
{
    command.SetThreadSize(common.thread_size);
    command.SetVerboseLevel(common.verbose_level);
    command.SetFolderPath(common.folder_path);
    if constexpr (HasCommonOption(CommandType::kCommonOptions, CommonOption::Database))
    {
        command.SetDatabasePath(common.database_path);
    }
}

template <typename CommandType, typename ConfigureFn>
ExecutionResult ExecuteCommand(
    ConfigureFn && configure)
{
    CommandType command{};
    configure(command);

    ExecutionResult result;
    result.prepared = command.PrepareForExecution();
    result.validation_issues = command.GetValidationIssues();
    if (!result.prepared)
    {
        return result;
    }

    result.executed = command.Execute();
    const auto & execution_issues{ command.GetValidationIssues() };
    if (!execution_issues.empty())
    {
        result.validation_issues = execution_issues;
    }
    return result;
}

} // namespace

GuiCommandExecutor::GuiCommandExecutor() = default;

ExecutionResult GuiCommandExecutor::ExecuteMapSimulation(
    const MapSimulationRequest & request)
{
    return ExecuteCommand<MapSimulationCommand>([&](MapSimulationCommand & configured)
    {
        ApplyCommonOptions(configured, request.common);
        configured.SetModelFilePath(request.model_file_path);
        configured.SetMapFileName(request.map_file_name);
        configured.SetPotentialModelChoice(request.potential_model_choice);
        configured.SetPartialChargeChoice(request.partial_charge_choice);
        configured.SetCutoffDistance(request.cutoff_distance);
        configured.SetGridSpacing(request.grid_spacing);
        configured.SetBlurringWidthList(request.blurring_width_list);
    });
}

ExecutionResult GuiCommandExecutor::ExecutePotentialAnalysis(
    const PotentialAnalysisRequest & request)
{
    return ExecuteCommand<PotentialAnalysisCommand>([&](PotentialAnalysisCommand & configured)
    {
        ApplyCommonOptions(configured, request.common);
        configured.SetModelFilePath(request.model_file_path);
        configured.SetMapFilePath(request.map_file_path);
        configured.SetSimulationFlag(request.simulation_flag);
        configured.SetSimulatedMapResolution(request.simulated_map_resolution);
        configured.SetSavedKeyTag(request.saved_key_tag);
        configured.SetTrainingReportDir(request.training_report_dir);
        configured.SetTrainingAlphaFlag(request.training_alpha_flag);
        configured.SetAsymmetryFlag(request.asymmetry_flag);
        configured.SetSamplingSize(request.sampling_size);
        configured.SetSamplingRangeMinimum(request.sampling_range_min);
        configured.SetSamplingRangeMaximum(request.sampling_range_max);
        configured.SetSamplingHeight(request.sampling_height);
        configured.SetFitRangeMinimum(request.fit_range_min);
        configured.SetFitRangeMaximum(request.fit_range_max);
        configured.SetAlphaR(request.alpha_r);
        configured.SetAlphaG(request.alpha_g);
    });
}

ExecutionResult GuiCommandExecutor::ExecuteResultDump(
    const ResultDumpRequest & request)
{
    return ExecuteCommand<ResultDumpCommand>([&](ResultDumpCommand & configured)
    {
        ApplyCommonOptions(configured, request.common);
        configured.SetPrinterChoice(request.printer_choice);
        configured.SetModelKeyTagList(request.model_key_tag_list);
        configured.SetMapFilePath(request.map_file_path);
    });
}

} // namespace rhbm_gem::gui
