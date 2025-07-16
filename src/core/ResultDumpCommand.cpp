#include "ResultDumpCommand.hpp"
#include "DataObjectManager.hpp"
#include "ResultDumpVisitor.hpp"
#include "Logger.hpp"

#include <sstream>

ResultDumpCommand::ResultDumpCommand(void) :
    m_printer_choice{ 2 }, m_database_path{ "database.sqlite" },
    m_folder_path{ "" }, m_map_file_path{ "" }
{

}

ResultDumpCommand::ResultDumpCommand(
    const Options & options, const GlobalOptions & globals) :
    m_printer_choice{ options.printer_choice },
    m_database_path{ globals.database_path },
    m_folder_path{ globals.folder_path },
    m_map_file_path{ options.map_file_path }
{
    SetModelKeyTagList(options.model_key_tag_list);
}

void ResultDumpCommand::RegisterCLIOptions(CLI::App * cmd, Options & options)
{
    cmd->add_option("-p,--printer", options.printer_choice,
        "Printer choice")->required();
    cmd->add_option("-k,--model-keylist", options.model_key_tag_list,
        "List of model key tag to be display")->required();
    cmd->add_option("-m,--map", options.map_file_path,
        "Map file path")->default_val("");
}

void ResultDumpCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "ResultDumpCommand::Execute() called.");
    Logger::Log(LogLevel::Info, "Total number of model object sets to be print: "
                + std::to_string(m_model_key_tag_list.size()));

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    try
    {
        if (m_map_file_path != "")
        {
            data_manager->ProcessFile(m_map_file_path, "map");
        }
        for (auto & key : m_model_key_tag_list)
        {
            data_manager->LoadDataObject(key);
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
        return;
    }

    auto result_dump{ std::make_unique<ResultDumpVisitor>() };
    result_dump->SetMapObjectKeyTag("map");
    result_dump->SetModelObjectKeyTagList(m_model_key_tag_list);
    result_dump->SetFolderPath(m_folder_path);
    result_dump->SetPrinterChoice(m_printer_choice);
    data_manager->Accept(result_dump.get());
}

void ResultDumpCommand::SetModelKeyTagList(const std::string & value)
{
    m_model_key_tag_list.clear();
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_model_key_tag_list.emplace_back(segment);
    }
}
