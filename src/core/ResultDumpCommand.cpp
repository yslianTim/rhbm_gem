#include "ResultDumpCommand.hpp"
#include "DataObjectManager.hpp"
#include "ResultDumpVisitor.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

namespace {
CommandRegistrar<ResultDumpCommand> registrar_result_dump{
    "result_dump",
    "Run result dump"};
}

ResultDumpCommand::ResultDumpCommand(void)
{

}

void ResultDumpCommand::RegisterCLIOptions(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::RegisterCLIOptions() called.");
    cmd->add_option("-p,--printer", m_options.printer_choice,
        "Printer choice")->required();
    cmd->add_option("-k,--model-keylist", m_options.model_key_tag_list,
        "List of model key tag to be display")->required()->delimiter(',');
    cmd->add_option("-m,--map", m_options.map_file_path,
        "Map file path")->default_val(m_options.map_file_path);
    RegisterCommandOptions(cmd);
}

bool ResultDumpCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::Execute() called.");
    Logger::Log(LogLevel::Info, "Total number of model object sets to be print: "
                + std::to_string(m_options.model_key_tag_list.size()));

    auto data_manager{ std::make_unique<DataObjectManager>(m_options.database_path) };
    try
    {
        if (m_options.map_file_path != "")
        {
            data_manager->ProcessFile(m_options.map_file_path, "map");
        }
        for (auto & key : m_options.model_key_tag_list)
        {
            data_manager->LoadDataObject(key);
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, "ResultDumpCommand::Execute() : " + std::string(e.what()));
        return false;
    }

    auto result_dump{ std::make_unique<ResultDumpVisitor>(m_options) };
    std::vector<std::string> key_list{ m_options.model_key_tag_list };
    if (m_options.map_file_path != "")
    {
        key_list.emplace_back("map");
    }
    data_manager->Accept(result_dump.get(), key_list);
    result_dump->Finalize();
    return true;
}

void ResultDumpCommand::SetModelKeyTagList(const std::string & value)
{
    m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
}
