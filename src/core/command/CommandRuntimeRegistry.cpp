#include "internal/CommandRuntimeRegistry.hpp"

#include "internal/CommandOptionBinding.hpp"
#include <rhbm_gem/core/command/CommandMetadata.hpp>

#include <memory>
#include <optional>
#include <string>

namespace {

using rhbm_gem::CommonCommandRequest;
using rhbm_gem::CommonOption;
using rhbm_gem::CommonOptionMask;
using rhbm_gem::HasCommonOption;

void BindCommonOptions(
    CLI::App * command,
    CommonCommandRequest & common,
    CommonOptionMask common_options)
{
    if (HasCommonOption(common_options, CommonOption::Threading))
    {
        rhbm_gem::command_cli::AddScalarOption<int>(
            command,
            "-j,--jobs",
            [&](int value) { common.thread_size = value; },
            "Number of threads",
            common.thread_size);
    }
    if (HasCommonOption(common_options, CommonOption::Verbose))
    {
        rhbm_gem::command_cli::AddScalarOption<int>(
            command,
            "-v,--verbose",
            [&](int value) { common.verbose_level = value; },
            "Verbose level",
            common.verbose_level);
    }
    if (HasCommonOption(common_options, CommonOption::Database))
    {
        rhbm_gem::command_cli::AddPathOption(
            command,
            "-d,--database",
            [&](const std::filesystem::path & value) { common.database_path = value; },
            "Database file path",
            common.database_path);
    }
    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        rhbm_gem::command_cli::AddPathOption(
            command,
            "-o,--folder",
            [&](const std::filesystem::path & value) { common.folder_path = value; },
            "folder path for output files",
            common.folder_path);
    }
}

} // namespace

namespace rhbm_gem {

CommandRunner BindPotentialAnalysisRuntime(CLI::App * command)
{
    auto request{ std::make_shared<PotentialAnalysisRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::DatabaseWorkflow));

    command_cli::AddPathOption(
        command,
        "-a,--model",
        [request](const std::filesystem::path & value) { request->model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "-m,--map",
        [request](const std::filesystem::path & value) { request->map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "--training-report-dir",
        [request](const std::filesystem::path & value) { request->training_report_dir = value; },
        "Optional output directory for alpha training reports");
    command_cli::AddScalarOption<bool>(
        command,
        "--simulation",
        [request](bool value) { request->simulation_flag = value; },
        "Simulation flag",
        request->simulation_flag);
    command_cli::AddScalarOption<double>(
        command,
        "-r,--sim-resolution",
        [request](double value) { request->simulated_map_resolution = value; },
        "Set simulated map's resolution (blurring width)",
        request->simulated_map_resolution);
    command_cli::AddStringOption(
        command,
        "-k,--save-key",
        [request](const std::string & value) { request->saved_key_tag = value; },
        "New key tag for saving ModelObject results into database",
        request->saved_key_tag);
    command_cli::AddScalarOption<bool>(
        command,
        "--training-alpha",
        [request](bool value) { request->training_alpha_flag = value; },
        "Turn On/Off alpha training flag",
        request->training_alpha_flag);
    command_cli::AddScalarOption<bool>(
        command,
        "--asymmetry",
        [request](bool value) { request->asymmetry_flag = value; },
        "Turn On/Off asymmetry flag",
        request->asymmetry_flag);
    command_cli::AddScalarOption<int>(
        command,
        "-s,--sampling",
        [request](int value) { request->sampling_size = value; },
        "Number of sampling points per atom",
        request->sampling_size);
    command_cli::AddScalarOption<double>(
        command,
        "--sampling-min",
        [request](double value) { request->sampling_range_min = value; },
        "Minimum sampling range",
        request->sampling_range_min);
    command_cli::AddScalarOption<double>(
        command,
        "--sampling-max",
        [request](double value) { request->sampling_range_max = value; },
        "Maximum sampling range",
        request->sampling_range_max);
    command_cli::AddScalarOption<double>(
        command,
        "--sampling-height",
        [request](double value) { request->sampling_height = value; },
        "Maximum sampling height",
        request->sampling_height);
    command_cli::AddScalarOption<double>(
        command,
        "--fit-min",
        [request](double value) { request->fit_range_min = value; },
        "Minimum fitting range",
        request->fit_range_min);
    command_cli::AddScalarOption<double>(
        command,
        "--fit-max",
        [request](double value) { request->fit_range_max = value; },
        "Maximum fitting range",
        request->fit_range_max);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha-r",
        [request](double value) { request->alpha_r = value; },
        "Alpha value for R",
        request->alpha_r);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha-g",
        [request](double value) { request->alpha_g = value; },
        "Alpha value for G",
        request->alpha_g);

    return [request]() { return RunPotentialAnalysis(*request); };
}

CommandRunner BindPotentialDisplayRuntime(CLI::App * command)
{
    auto request{ std::make_shared<PotentialDisplayRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::DatabaseWorkflow));

    command_cli::AddEnumOption<PainterType>(
        command,
        "-p,--painter",
        [request](PainterType value) { request->painter_choice = value; },
        "Painter choice",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        command,
        "-k,--model-keylist",
        [request](const std::string & value) { request->model_key_tag_list = value; },
        "List of model key tag to be display",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        command,
        "-r,--ref-model-keylist",
        [request](const std::string & value) { request->ref_model_key_tag_list = value; },
        "List of reference model key tag to be display",
        request->ref_model_key_tag_list);
    command_cli::AddStringOption(
        command,
        "--pick-chain",
        [request](const std::string & value) { request->pick_chain_id = value; },
        "Pick chain ID",
        request->pick_chain_id);
    command_cli::AddStringOption(
        command,
        "--pick-residue",
        [request](const std::string & value) { request->pick_residue = value; },
        "Pick residue type",
        request->pick_residue);
    command_cli::AddStringOption(
        command,
        "--pick-element",
        [request](const std::string & value) { request->pick_element = value; },
        "Pick element type",
        request->pick_element);
    command_cli::AddStringOption(
        command,
        "--veto-chain",
        [request](const std::string & value) { request->veto_chain_id = value; },
        "Veto chain ID",
        request->veto_chain_id);
    command_cli::AddStringOption(
        command,
        "--veto-residue",
        [request](const std::string & value) { request->veto_residue = value; },
        "Veto residue type",
        request->veto_residue);
    command_cli::AddStringOption(
        command,
        "--veto-element",
        [request](const std::string & value) { request->veto_element = value; },
        "Veto element type",
        request->veto_element);

    return [request]() { return RunPotentialDisplay(*request); };
}

CommandRunner BindResultDumpRuntime(CLI::App * command)
{
    auto request{ std::make_shared<ResultDumpRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::DatabaseWorkflow));

    command_cli::AddEnumOption<PrinterType>(
        command,
        "-p,--printer",
        [request](PrinterType value) { request->printer_choice = value; },
        "Printer choice",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        command,
        "-k,--model-keylist",
        [request](const std::string & value) { request->model_key_tag_list = value; },
        "List of model key tag to be display",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "-m,--map",
        [request](const std::filesystem::path & value) { request->map_file_path = value; },
        "Map file path",
        request->map_file_path);

    return [request]() { return RunResultDump(*request); };
}

CommandRunner BindMapSimulationRuntime(CLI::App * command)
{
    auto request{ std::make_shared<MapSimulationRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::FileWorkflow));

    // Keep CLI compatibility: omitted --blurring-width stays empty and fails validation.
    request->blurring_width_list.clear();

    command_cli::AddPathOption(
        command,
        "-a,--model",
        [request](const std::filesystem::path & value) { request->model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        command,
        "-n,--name",
        [request](const std::string & value) { request->map_file_name = value; },
        "File name for output map files",
        request->map_file_name);
    command_cli::AddEnumOption<PotentialModel>(
        command,
        "--potential-model",
        [request](PotentialModel value) { request->potential_model_choice = value; },
        "Atomic potential model option",
        PotentialModel::FIVE_GAUS_CHARGE);
    command_cli::AddEnumOption<PartialCharge>(
        command,
        "--charge",
        [request](PartialCharge value) { request->partial_charge_choice = value; },
        "Partial charge table option",
        PartialCharge::PARTIAL);
    command_cli::AddScalarOption<double>(
        command,
        "-c,--cut-off",
        [request](double value) { request->cutoff_distance = value; },
        "Cutoff distance",
        request->cutoff_distance);
    command_cli::AddScalarOption<double>(
        command,
        "-g,--grid-spacing",
        [request](double value) { request->grid_spacing = value; },
        "Grid spacing",
        request->grid_spacing);
    command_cli::AddStringOption(
        command,
        "--blurring-width",
        [request](const std::string & value) { request->blurring_width_list = value; },
        "Blurring width (list) setting",
        request->blurring_width_list);

    return [request]() { return RunMapSimulation(*request); };
}

CommandRunner BindMapVisualizationRuntime(CLI::App * command)
{
    auto request{ std::make_shared<MapVisualizationRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::FileWorkflow));

    command_cli::AddPathOption(
        command,
        "-a,--model",
        [request](const std::filesystem::path & value) { request->model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "-m,--map",
        [request](const std::filesystem::path & value) { request->map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    command_cli::AddScalarOption<int>(
        command,
        "-i,--atom-id",
        [request](int value) { request->atom_serial_id = value; },
        "Atom serial ID for visualization",
        request->atom_serial_id);
    command_cli::AddScalarOption<int>(
        command,
        "-s,--sampling",
        [request](int value) { request->sampling_size = value; },
        "Number of sampling points per atom",
        request->sampling_size);
    command_cli::AddScalarOption<double>(
        command,
        "--window-size",
        [request](double value) { request->window_size = value; },
        "Window size for sampling",
        request->window_size);

    return [request]() { return RunMapVisualization(*request); };
}

CommandRunner BindPositionEstimationRuntime(CLI::App * command)
{
    auto request{ std::make_shared<PositionEstimationRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::FileWorkflow));

    command_cli::AddPathOption(
        command,
        "-m,--map",
        [request](const std::filesystem::path & value) { request->map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    command_cli::AddScalarOption<int>(
        command,
        "--iter",
        [request](int value) { request->iteration_count = value; },
        "Iteration count for estimation",
        request->iteration_count);
    command_cli::AddScalarOption<int>(
        command,
        "--knn",
        [request](int value) { request->knn_size = value; },
        "KNN size for estimation",
        request->knn_size);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha",
        [request](double value) { request->alpha = value; },
        "Alpha value for robust regression",
        request->alpha);
    command_cli::AddScalarOption<double>(
        command,
        "--threshold",
        [request](double value) { request->threshold_ratio = value; },
        "Ratio of threshold of map values",
        request->threshold_ratio);
    command_cli::AddScalarOption<double>(
        command,
        "--dedup-tolerance",
        [request](double value) { request->dedup_tolerance = value; },
        "Tolerance for deduplicating points",
        request->dedup_tolerance);

    return [request]() { return RunPositionEstimation(*request); };
}

CommandRunner BindModelTestRuntime(CLI::App * command)
{
    auto request{ std::make_shared<HRLModelTestRequest>() };
    BindCommonOptions(
        command,
        request->common,
        CommonOptionMaskForProfile(CommonOptionProfile::FileWorkflow));

    command_cli::AddEnumOption<TesterType>(
        command,
        "-t,--tester",
        [request](TesterType value) { request->tester_choice = value; },
        "Tester option",
        TesterType::DATA_OUTLIER);
    command_cli::AddScalarOption<double>(
        command,
        "--fit-min",
        [request](double value) { request->fit_range_min = value; },
        "Minimum fitting range",
        request->fit_range_min);
    command_cli::AddScalarOption<double>(
        command,
        "--fit-max",
        [request](double value) { request->fit_range_max = value; },
        "Maximum fitting range",
        request->fit_range_max);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha-r",
        [request](double value) { request->alpha_r = value; },
        "Alpha value for R",
        request->alpha_r);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha-g",
        [request](double value) { request->alpha_g = value; },
        "Alpha value for G",
        request->alpha_g);

    return [request]() { return RunHRLModelTest(*request); };
}

} // namespace rhbm_gem
