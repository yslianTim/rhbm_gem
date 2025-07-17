#pragma once

#include <memory>
#include <string>
#include <vector>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "GlobalOptions.hpp"

class ResultDumpCommand : public CommandBase
{
public:
    struct Options
    {
        int printer_choice{ 2 };
        std::string model_key_tag_list;
        std::string map_file_path;
    };

private:
    Options m_options{};
    GlobalOptions m_globals{};
    std::vector<std::string> m_model_key_tag_list;

public:
    ResultDumpCommand(void);
    ~ResultDumpCommand() = default;
    void Execute(void) override;
    void RegisterCLIOptions(CLI::App * cmd) override;
    void SetGlobalOptions(const GlobalOptions & globals) override { m_globals = globals; }

    void SetPrinterChoice(int value) { m_options.printer_choice = value; }
    void SetDatabasePath(const std::string & path) { m_globals.database_path = path; }
    void SetFolderPath(const std::string & path) { m_globals.folder_path = path; }
    void SetMapFilePath(const std::string & path) { m_options.map_file_path = path; }
    void SetModelKeyTagList(const std::string & value);

};
