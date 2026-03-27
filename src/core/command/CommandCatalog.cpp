#include "internal/command/CommandCatalog.hpp"

#include <CLI/CLI.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/CommandMetadata.hpp>
#include <rhbm_gem/core/command/OptionEnumClass.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>

namespace {

using rhbm_gem::CommonCommandRequest;
using rhbm_gem::CommonOption;
using rhbm_gem::CommonOptionProfile;
using rhbm_gem::HasCommonOption;
using rhbm_gem::PainterType;
using rhbm_gem::PartialCharge;
using rhbm_gem::PotentialModel;
using rhbm_gem::PrinterType;
using rhbm_gem::TesterType;

template <typename ValueType, typename Setter>
CLI::Option * AddScalarOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<ValueType> default_value = std::nullopt,
    bool required = false)
{
    auto * option{
        command->add_option_function<ValueType>(
            std::string(flags),
            [callback = std::forward<Setter>(setter)](const ValueType & value) mutable
            {
                std::invoke(callback, value);
            },
            std::string(help))
    };
    if (required)
    {
        option->required();
    }
    if (default_value.has_value())
    {
        option->default_val(*default_value);
    }
    return option;
}

template <typename Setter>
CLI::Option * AddStringOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<std::string> default_value = std::nullopt,
    bool required = false)
{
    return AddScalarOption<std::string>(
        command,
        flags,
        std::forward<Setter>(setter),
        help,
        std::move(default_value),
        required);
}

template <typename Setter>
CLI::Option * AddPathOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<std::filesystem::path> default_value = std::nullopt,
    bool required = false)
{
    auto * option{
        command->add_option_function<std::string>(
            std::string(flags),
            [callback = std::forward<Setter>(setter)](const std::string & value) mutable
            {
                std::invoke(callback, std::filesystem::path{ value });
            },
            std::string(help))
    };
    if (required)
    {
        option->required();
    }
    if (default_value.has_value())
    {
        option->default_val(default_value->string());
    }
    return option;
}

template <typename EnumType, typename Setter>
CLI::Option * AddEnumOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<EnumType> default_value = std::nullopt,
    bool required = false)
{
    auto * option{
        AddScalarOption<EnumType>(
            command,
            flags,
            std::forward<Setter>(setter),
            help,
            std::move(default_value),
            required)
    };
    option->transform(CLI::CheckedTransformer(
        rhbm_gem::BuildEnumCLIMap<EnumType>(), CLI::ignore_case));
    return option;
}

void BindCommonOptions(
    CLI::App * command,
    CommonCommandRequest & common,
    CommonOptionProfile profile)
{
    const auto common_options{ rhbm_gem::CommonOptionMaskForProfile(profile) };
    if (HasCommonOption(common_options, CommonOption::Threading))
    {
        AddScalarOption<int>(
            command,
            "-j,--jobs",
            [&](int value) { common.thread_size = value; },
            "Number of threads",
            common.thread_size);
    }
    if (HasCommonOption(common_options, CommonOption::Verbose))
    {
        AddScalarOption<int>(
            command,
            "-v,--verbose",
            [&](int value) { common.verbose_level = value; },
            "Verbose level",
            common.verbose_level);
    }
    if (HasCommonOption(common_options, CommonOption::Database))
    {
        AddPathOption(
            command,
            "-d,--database",
            [&](const std::filesystem::path & value) { common.database_path = value; },
            "Database file path",
            common.database_path);
    }
    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        AddPathOption(
            command,
            "-o,--folder",
            [&](const std::filesystem::path & value) { common.folder_path = value; },
            "folder path for output files",
            common.folder_path);
    }
}

template <typename RequestType, auto RunCommandFn, typename BindOptionsFn>
void RegisterCommand(
    CLI::App & app,
    const rhbm_gem::CommandDescriptor & descriptor,
    BindOptionsFn bind_options)
{
    CLI::App * command{
        app.add_subcommand(
            std::string(descriptor.name),
            std::string(descriptor.description))
    };
    auto request{ std::make_shared<RequestType>() };
    BindCommonOptions(command, request->common, descriptor.profile);
    std::invoke(bind_options, command, *request);
    command->callback([request]()
    {
        ScopeTimer timer("Command CLI callback");
        const auto report{ RunCommandFn(*request) };
        if (!report.executed)
        {
            throw CLI::RuntimeError(1);
        }
    });
}

void BindPotentialAnalysisRequestOptions(
    CLI::App * command,
    rhbm_gem::PotentialAnalysisRequest & request)
{
    AddPathOption(
        command,
        "-a,--model",
        [&](const std::filesystem::path & value) { request.model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    AddPathOption(
        command,
        "-m,--map",
        [&](const std::filesystem::path & value) { request.map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    AddPathOption(
        command,
        "--training-report-dir",
        [&](const std::filesystem::path & value) { request.training_report_dir = value; },
        "Optional output directory for alpha training reports");
    AddScalarOption<bool>(
        command,
        "--simulation",
        [&](bool value) { request.simulation_flag = value; },
        "Simulation flag",
        request.simulation_flag);
    AddScalarOption<double>(
        command,
        "-r,--sim-resolution",
        [&](double value) { request.simulated_map_resolution = value; },
        "Set simulated map's resolution (blurring width)",
        request.simulated_map_resolution);
    AddStringOption(
        command,
        "-k,--save-key",
        [&](const std::string & value) { request.saved_key_tag = value; },
        "New key tag for saving ModelObject results into database",
        request.saved_key_tag);
    AddScalarOption<bool>(
        command,
        "--training-alpha",
        [&](bool value) { request.training_alpha_flag = value; },
        "Turn On/Off alpha training flag",
        request.training_alpha_flag);
    AddScalarOption<bool>(
        command,
        "--asymmetry",
        [&](bool value) { request.asymmetry_flag = value; },
        "Turn On/Off asymmetry flag",
        request.asymmetry_flag);
    AddScalarOption<int>(
        command,
        "-s,--sampling",
        [&](int value) { request.sampling_size = value; },
        "Number of sampling points per atom",
        request.sampling_size);
    AddScalarOption<double>(
        command,
        "--sampling-min",
        [&](double value) { request.sampling_range_min = value; },
        "Minimum sampling range",
        request.sampling_range_min);
    AddScalarOption<double>(
        command,
        "--sampling-max",
        [&](double value) { request.sampling_range_max = value; },
        "Maximum sampling range",
        request.sampling_range_max);
    AddScalarOption<double>(
        command,
        "--sampling-height",
        [&](double value) { request.sampling_height = value; },
        "Maximum sampling height",
        request.sampling_height);
    AddScalarOption<double>(
        command,
        "--fit-min",
        [&](double value) { request.fit_range_min = value; },
        "Minimum fitting range",
        request.fit_range_min);
    AddScalarOption<double>(
        command,
        "--fit-max",
        [&](double value) { request.fit_range_max = value; },
        "Maximum fitting range",
        request.fit_range_max);
    AddScalarOption<double>(
        command,
        "--alpha-r",
        [&](double value) { request.alpha_r = value; },
        "Alpha value for R",
        request.alpha_r);
    AddScalarOption<double>(
        command,
        "--alpha-g",
        [&](double value) { request.alpha_g = value; },
        "Alpha value for G",
        request.alpha_g);
}

void BindPotentialDisplayRequestOptions(
    CLI::App * command,
    rhbm_gem::PotentialDisplayRequest & request)
{
    AddEnumOption<PainterType>(
        command,
        "-p,--painter",
        [&](PainterType value) { request.painter_choice = value; },
        "Painter choice",
        std::nullopt,
        true);
    AddStringOption(
        command,
        "-k,--model-keylist",
        [&](const std::string & value) { request.model_key_tag_list = value; },
        "List of model key tag to be display",
        std::nullopt,
        true);
    AddStringOption(
        command,
        "-r,--ref-model-keylist",
        [&](const std::string & value) { request.ref_model_key_tag_list = value; },
        "List of reference model key tag to be display",
        request.ref_model_key_tag_list);
    AddStringOption(
        command,
        "--pick-chain",
        [&](const std::string & value) { request.pick_chain_id = value; },
        "Pick chain ID",
        request.pick_chain_id);
    AddStringOption(
        command,
        "--pick-residue",
        [&](const std::string & value) { request.pick_residue = value; },
        "Pick residue type",
        request.pick_residue);
    AddStringOption(
        command,
        "--pick-element",
        [&](const std::string & value) { request.pick_element = value; },
        "Pick element type",
        request.pick_element);
    AddStringOption(
        command,
        "--veto-chain",
        [&](const std::string & value) { request.veto_chain_id = value; },
        "Veto chain ID",
        request.veto_chain_id);
    AddStringOption(
        command,
        "--veto-residue",
        [&](const std::string & value) { request.veto_residue = value; },
        "Veto residue type",
        request.veto_residue);
    AddStringOption(
        command,
        "--veto-element",
        [&](const std::string & value) { request.veto_element = value; },
        "Veto element type",
        request.veto_element);
}

void BindResultDumpRequestOptions(
    CLI::App * command,
    rhbm_gem::ResultDumpRequest & request)
{
    AddEnumOption<PrinterType>(
        command,
        "-p,--printer",
        [&](PrinterType value) { request.printer_choice = value; },
        "Printer choice",
        std::nullopt,
        true);
    AddStringOption(
        command,
        "-k,--model-keylist",
        [&](const std::string & value) { request.model_key_tag_list = value; },
        "List of model key tag to be display",
        std::nullopt,
        true);
    AddPathOption(
        command,
        "-m,--map",
        [&](const std::filesystem::path & value) { request.map_file_path = value; },
        "Map file path",
        request.map_file_path);
}

void BindMapSimulationRequestOptions(
    CLI::App * command,
    rhbm_gem::MapSimulationRequest & request)
{
    request.blurring_width_list.clear();

    AddPathOption(
        command,
        "-a,--model",
        [&](const std::filesystem::path & value) { request.model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    AddStringOption(
        command,
        "-n,--name",
        [&](const std::string & value) { request.map_file_name = value; },
        "File name for output map files",
        request.map_file_name);
    AddEnumOption<PotentialModel>(
        command,
        "--potential-model",
        [&](PotentialModel value) { request.potential_model_choice = value; },
        "Atomic potential model option",
        PotentialModel::FIVE_GAUS_CHARGE);
    AddEnumOption<PartialCharge>(
        command,
        "--charge",
        [&](PartialCharge value) { request.partial_charge_choice = value; },
        "Partial charge table option",
        PartialCharge::PARTIAL);
    AddScalarOption<double>(
        command,
        "-c,--cut-off",
        [&](double value) { request.cutoff_distance = value; },
        "Cutoff distance",
        request.cutoff_distance);
    AddScalarOption<double>(
        command,
        "-g,--grid-spacing",
        [&](double value) { request.grid_spacing = value; },
        "Grid spacing",
        request.grid_spacing);
    AddStringOption(
        command,
        "--blurring-width",
        [&](const std::string & value) { request.blurring_width_list = value; },
        "Blurring width (list) setting",
        request.blurring_width_list);
}

#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
void BindMapVisualizationRequestOptions(
    CLI::App * command,
    rhbm_gem::MapVisualizationRequest & request)
{
    AddPathOption(
        command,
        "-a,--model",
        [&](const std::filesystem::path & value) { request.model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    AddPathOption(
        command,
        "-m,--map",
        [&](const std::filesystem::path & value) { request.map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    AddScalarOption<int>(
        command,
        "-i,--atom-id",
        [&](int value) { request.atom_serial_id = value; },
        "Atom serial ID for visualization",
        request.atom_serial_id);
    AddScalarOption<int>(
        command,
        "-s,--sampling",
        [&](int value) { request.sampling_size = value; },
        "Number of sampling points per atom",
        request.sampling_size);
    AddScalarOption<double>(
        command,
        "--window-size",
        [&](double value) { request.window_size = value; },
        "Window size for sampling",
        request.window_size);
}

void BindPositionEstimationRequestOptions(
    CLI::App * command,
    rhbm_gem::PositionEstimationRequest & request)
{
    AddPathOption(
        command,
        "-m,--map",
        [&](const std::filesystem::path & value) { request.map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    AddScalarOption<int>(
        command,
        "--iter",
        [&](int value) { request.iteration_count = value; },
        "Iteration count for estimation",
        request.iteration_count);
    AddScalarOption<int>(
        command,
        "--knn",
        [&](int value) { request.knn_size = value; },
        "KNN size for estimation",
        request.knn_size);
    AddScalarOption<double>(
        command,
        "--alpha",
        [&](double value) { request.alpha = value; },
        "Alpha value for robust regression",
        request.alpha);
    AddScalarOption<double>(
        command,
        "--threshold",
        [&](double value) { request.threshold_ratio = value; },
        "Ratio of threshold of map values",
        request.threshold_ratio);
    AddScalarOption<double>(
        command,
        "--dedup-tolerance",
        [&](double value) { request.dedup_tolerance = value; },
        "Tolerance for deduplicating points",
        request.dedup_tolerance);
}
#endif

void BindHRLModelTestRequestOptions(
    CLI::App * command,
    rhbm_gem::HRLModelTestRequest & request)
{
    AddEnumOption<TesterType>(
        command,
        "-t,--tester",
        [&](TesterType value) { request.tester_choice = value; },
        "Tester option",
        TesterType::DATA_OUTLIER);
    AddScalarOption<double>(
        command,
        "--fit-min",
        [&](double value) { request.fit_range_min = value; },
        "Minimum fitting range",
        request.fit_range_min);
    AddScalarOption<double>(
        command,
        "--fit-max",
        [&](double value) { request.fit_range_max = value; },
        "Maximum fitting range",
        request.fit_range_max);
    AddScalarOption<double>(
        command,
        "--alpha-r",
        [&](double value) { request.alpha_r = value; },
        "Alpha value for R",
        request.alpha_r);
    AddScalarOption<double>(
        command,
        "--alpha-g",
        [&](double value) { request.alpha_g = value; },
        "Alpha value for G",
        request.alpha_g);
}

} // namespace

namespace rhbm_gem {

const std::vector<CommandDescriptor> & CommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
        CommandDescriptor{                                                                      \
            CommandId::COMMAND_ID,                                                              \
            CLI_NAME,                                                                           \
            DESCRIPTION,                                                                        \
            CommonOptionProfile::PROFILE                                                        \
        },
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
    };

    return catalog;
}

void ConfigureCommandCli(CLI::App & app)
{
    app.require_subcommand(1);
    RegisterCommandSubcommands(app);
}

void RegisterCommandSubcommands(CLI::App & app)
{
    for (const auto & descriptor : CommandCatalog())
    {
        switch (descriptor.id)
        {
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
        case CommandId::COMMAND_ID:                                                             \
            RegisterCommand<COMMAND_ID##Request, &Run##COMMAND_ID>(                             \
                app,                                                                            \
                descriptor,                                                                     \
                Bind##COMMAND_ID##RequestOptions);                                              \
            break;
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
        default:
            break;
        }
    }
}

} // namespace rhbm_gem
