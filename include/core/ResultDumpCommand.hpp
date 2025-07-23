#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

class AtomObject;
class ModelObject;
class MapObject;

class ResultDumpCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        int printer_choice{ 2 };
        std::vector<std::string> model_key_tag_list{};
        std::filesystem::path map_file_path{ "" };
    };

private:
    Options m_options{};
    std::string m_map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<ModelObject *> m_model_object_list;
    MapObject * m_map_object;

public:
    ResultDumpCommand(void);
    ~ResultDumpCommand() = default;
    bool Execute(void) override;
    void RegisterCLIOptions(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }

    void SetPrinterChoice(int value) { m_options.printer_choice = value; }
    void SetDatabasePath(const std::filesystem::path & path) { m_options.database_path = path; }
    void SetFolderPath(const std::filesystem::path & path) { m_options.folder_path = path; }
    void SetMapFilePath(const std::filesystem::path & path) { m_options.map_file_path = path; }
    void SetModelKeyTagList(const std::string & value);

private:
    void BuildSelectedAtomList(void);
    void RunAtomPositionDumping(void);
    void RunMapValueDumping(void);
    void RunGausEstimatesDumping(void);
    void RunGroupGausEstimatesDumping(void);

};
