#include "ResultDumpCommand.hpp"
#include "DataObjectManager.hpp"
#include "ResultDumpVisitor.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

REGISTER_COMMAND(ResultDumpCommand, "result_dump", "Run result dump");

ResultDumpCommand::ResultDumpCommand(void)
{

}

void ResultDumpCommand::RegisterCLIOptions(CLI::App * cmd)
{
    cmd->add_option("-p,--printer", m_options.printer_choice,
        "Printer choice")->required();
    cmd->add_option("-k,--model-keylist", m_options.model_key_tag_list,
        "List of model key tag to be display")->required();
    cmd->add_option("-m,--map", m_options.map_file_path,
        "Map file path")->default_val(m_options.map_file_path);
    RegisterCommandOptions(cmd);
}

bool ResultDumpCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "ResultDumpCommand::Execute() called.");
    SetModelKeyTagList(m_options.model_key_tag_list);
    Logger::Log(LogLevel::Info, "Total number of model object sets to be print: "
                + std::to_string(m_model_key_tag_list.size()));

    auto data_manager{ std::make_unique<DataObjectManager>(m_options.database_path) };
    try
    {
        if (m_options.map_file_path != "")
        {
            data_manager->ProcessFile(m_options.map_file_path, "map");
        }
        for (auto & key : m_model_key_tag_list)
        {
            data_manager->LoadDataObject(key);
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
        return false;
    }

    auto result_dump{ std::make_unique<ResultDumpVisitor>() };
    result_dump->SetMapObjectKeyTag("map");
    result_dump->SetModelObjectKeyTagList(m_model_key_tag_list);
    result_dump->SetFolderPath(m_options.folder_path);
    result_dump->SetPrinterChoice(m_options.printer_choice);
    data_manager->Accept(result_dump.get());
    return true;
}

void ResultDumpCommand::SetModelKeyTagList(const std::string & value)
{
    m_model_key_tag_list.clear();
    for (const auto & token : StringHelper::SplitStringLineFromDelimiter(value, ','))
    {
        m_model_key_tag_list.emplace_back(token);
    }
}
