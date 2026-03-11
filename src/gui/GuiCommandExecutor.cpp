#include <rhbm_gem/gui/GuiCommandExecutor.hpp>

#include <rhbm_gem/core/command/MapSimulationCommand.hpp>
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/command/ResultDumpCommand.hpp>
#include <rhbm_gem/data/io/DataIoServices.hpp>
#include "internal/GuiCommandExecutorTestHooks.hpp"

#include <atomic>
#include <utility>

namespace rhbm_gem::gui {

namespace {

std::atomic<std::size_t> g_default_service_build_count{ 0U };
std::atomic<std::size_t> g_default_executor_build_count{ 0U };

const DataIoServices & SharedDataIoServices()
{
    static const DataIoServices data_io_services = []()
    {
        g_default_service_build_count.fetch_add(1U, std::memory_order_relaxed);
        return DataIoServices::BuildDefault();
    }();
    return data_io_services;
}

const GuiCommandExecutor & SharedThreadExecutor()
{
    thread_local const GuiCommandExecutor executor = []()
    {
        g_default_executor_build_count.fetch_add(1U, std::memory_order_relaxed);
        return GuiCommandExecutor{ SharedDataIoServices() };
    }();
    return executor;
}

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
    const DataIoServices & data_io_services,
    ConfigureFn && configure)
{
    CommandType command{ data_io_services };
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

GuiCommandExecutor::GuiCommandExecutor() :
    m_data_io_services{ DataIoServices::BuildDefault() }
{
}

GuiCommandExecutor::GuiCommandExecutor(DataIoServices data_io_services) :
    m_data_io_services{ std::move(data_io_services) }
{
}

ExecutionResult GuiCommandExecutor::RunMapSimulation(
    const MapSimulationRequest & request) const
{
    return ExecuteCommand<MapSimulationCommand>(m_data_io_services, [&](MapSimulationCommand & configured)
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

ExecutionResult GuiCommandExecutor::RunPotentialAnalysis(
    const PotentialAnalysisRequest & request) const
{
    return ExecuteCommand<PotentialAnalysisCommand>(m_data_io_services, [&](PotentialAnalysisCommand & configured)
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

ExecutionResult GuiCommandExecutor::RunResultDump(
    const ResultDumpRequest & request) const
{
    return ExecuteCommand<ResultDumpCommand>(m_data_io_services, [&](ResultDumpCommand & configured)
    {
        ApplyCommonOptions(configured, request.common);
        configured.SetPrinterChoice(request.printer_choice);
        configured.SetModelKeyTagList(request.model_key_tag_list);
        configured.SetMapFilePath(request.map_file_path);
    });
}

ExecutionResult GuiCommandExecutor::ExecuteMapSimulation(
    const MapSimulationRequest & request)
{
    return SharedThreadExecutor().RunMapSimulation(request);
}

ExecutionResult GuiCommandExecutor::ExecutePotentialAnalysis(
    const PotentialAnalysisRequest & request)
{
    return SharedThreadExecutor().RunPotentialAnalysis(request);
}

ExecutionResult GuiCommandExecutor::ExecuteResultDump(
    const ResultDumpRequest & request)
{
    return SharedThreadExecutor().RunResultDump(request);
}

namespace internal {

std::size_t DefaultServiceBuildCountForTesting()
{
    return g_default_service_build_count.load(std::memory_order_relaxed);
}

std::size_t DefaultExecutorBuildCountForTesting()
{
    return g_default_executor_build_count.load(std::memory_order_relaxed);
}

} // namespace internal

} // namespace rhbm_gem::gui
