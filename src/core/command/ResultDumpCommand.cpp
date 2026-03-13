#include "ResultDumpCommand.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

namespace rhbm_gem::detail {

bool ExecuteResultDumpWorkflow(
    DataObjectManager & data_manager,
    const ResultDumpCommandOptions & options,
    std::string & map_key_tag,
    std::unordered_map<std::string, std::vector<AtomObject *>> & selected_atom_list_map,
    std::vector<std::shared_ptr<ModelObject>> & model_object_list,
    std::shared_ptr<MapObject> & map_object);

} // namespace rhbm_gem::detail

namespace {
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kPrinterOption{ "--printer" };
constexpr std::string_view kModelKeyListOption{ "--model-keylist" };
constexpr std::string_view kMapOption{ "--map" };
}

namespace rhbm_gem {

ResultDumpCommand::ResultDumpCommand() :
    CommandWithProfileOptions<
        ResultDumpCommandOptions,
        CommandId::ResultDump,
        CommonOptionProfile::DatabaseWorkflow>{},
    m_map_key_tag{ kMapKey }, m_map_object{ nullptr }
{
}

void ResultDumpCommand::ApplyRequest(const ResultDumpRequest & request)
{
    ApplyCommonRequest(request.common);
    SetPrinterChoice(request.printer_choice);
    SetModelKeyTagList(request.model_key_tag_list);
    SetMapFilePath(request.map_file_path);
}

bool ResultDumpCommand::ExecuteImpl()
{
    return detail::ExecuteResultDumpWorkflow(
        m_data_manager,
        m_options,
        m_map_key_tag,
        m_selected_atom_list_map,
        m_model_object_list,
        m_map_object);
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
