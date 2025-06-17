#include "ResultDumpCommand.hpp"
#include "DataObjectManager.hpp"
#include "ResultDumpVisitor.hpp"

#include <iostream>
#include <sstream>

ResultDumpCommand::ResultDumpCommand(void)
{

}

ResultDumpCommand::~ResultDumpCommand()
{

}

void ResultDumpCommand::Execute(void)
{
    std::cout << "ResultDumpCommand::Execute() called." << std::endl;
    std::cout << "Total number of model object sets to be print: "
              << m_model_key_tag_list.size() << std::endl;

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
        std::cerr << e.what() << '\n';
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
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_model_key_tag_list.emplace_back(segment);
    }
}