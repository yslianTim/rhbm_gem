#include <rhbm_gem/core/command/ResultDumpCommand.hpp>
#include <rhbm_gem/core/command/CommandOptionBinding.hpp>
#include "workflow/ResultDumpWorkflowInternal.hpp"
#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace {
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kPrinterFlags{ "-p,--printer" };
constexpr std::string_view kPrinterOption{ "--printer" };
constexpr std::string_view kModelKeyListFlags{ "-k,--model-keylist" };
constexpr std::string_view kModelKeyListOption{ "--model-keylist" };
constexpr std::string_view kMapFlags{ "-m,--map" };
constexpr std::string_view kMapOption{ "--map" };
}

namespace rhbm_gem {

ResultDumpCommand::ResultDumpCommand(const DataIoServices & data_io_services) :
    CommandWithProfileOptions<
        ResultDumpCommandOptions,
        CommandId::ResultDump,
        CommonOptionProfile::DatabaseWorkflow>{ data_io_services },
    m_map_key_tag{ kMapKey }, m_map_object{ nullptr }
{
}

void ResultDumpCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddEnumOption<PrinterType>(
        cmd, kPrinterFlags,
        [&](PrinterType value) { SetPrinterChoice(value); },
        "Printer choice",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        cmd, kModelKeyListFlags,
        [&](const std::string & value) { SetModelKeyTagList(value); },
        "List of model key tag to be display",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        cmd, kMapFlags,
        [&](const std::filesystem::path & value) { SetMapFilePath(value); },
        "Map file path",
        m_options.map_file_path);
}

bool ResultDumpCommand::ExecuteImpl()
{
    detail::ResultDumpWorkflowContext context{
        m_data_manager,
        m_options,
        m_map_key_tag,
        m_selected_atom_list_map,
        m_model_object_list,
        m_map_object,
        [this]() { RequireDatabaseManager(); },
        [this](std::string_view stem, std::string_view extension)
        {
            return BuildOutputPath(stem, extension);
        },
    };
    return detail::ExecuteResultDumpWorkflow(context);
}

void ResultDumpCommand::ValidateOptions()
{
    ResetPrepareIssues(kMapOption);
    if (m_options.printer_choice == PrinterType::MAP_VALUE && m_options.map_file_path.empty())
    {
        AddValidationError(kMapOption,
            "A map file is required when '--printer map' is selected.");
    }
}

void ResultDumpCommand::ResetRuntimeState()
{
    m_selected_atom_list_map.clear();
    m_model_object_list.clear();
    m_map_object.reset();
}

void ResultDumpCommand::SetPrinterChoice(PrinterType value)
{
    SetValidatedEnumOption(
        m_options.printer_choice,
        value,
        kPrinterOption,
        PrinterType::GAUS_ESTIMATES,
        "Printer choice");
}

void ResultDumpCommand::SetMapFilePath(const std::filesystem::path & path)
{
    SetOptionalExistingPathOption(m_options.map_file_path, path, kMapOption, "Map file");
}

void ResultDumpCommand::SetModelKeyTagList(const std::string & value)
{
    MutateOptions([&]()
    {
        m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
        ResetParseIssues(kModelKeyListOption);
        if (m_options.model_key_tag_list.empty())
        {
            AddValidationError(
                kModelKeyListOption,
                "Model key list cannot be empty.",
                ValidationPhase::Parse);
        }
    });
}

} // namespace rhbm_gem
