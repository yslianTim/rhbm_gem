#pragma once

#include <memory>
#include <string>
#include <vector>
#include "CommandBase.hpp"

class ResultDumpCommand : public CommandBase
{
    int m_printer_choice;
    std::string m_database_path, m_folder_path, m_map_file_path;
    std::vector<std::string> m_model_key_tag_list;

public:
    ResultDumpCommand(void);
    ~ResultDumpCommand() = default;
    void Execute(void) override;

    void SetPrinterChoice(int value) { m_printer_choice = value; }
    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetFolderPath(const std::string & path) { m_folder_path = path; }
    void SetMapFilePath(const std::string & path) { m_map_file_path = path; }
    void SetModelKeyTagList(const std::string & value);

};
