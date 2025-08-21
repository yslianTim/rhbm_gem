#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

class AtomObject;
class ModelObject;
class MapObject;

class ResultDumpCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
        std::vector<std::string> model_key_tag_list{};
        std::filesystem::path map_file_path{ "" };
    };

private:
    Options m_options{};
    std::string m_map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::shared_ptr<MapObject> m_map_object;

public:
    ResultDumpCommand(void);
    ~ResultDumpCommand();
    bool Execute(void) override;
    bool ValidateOptions(void) const override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetPrinterChoice(PrinterType value) { m_options.printer_choice = value; }
    void SetDatabasePath(const std::filesystem::path & path) { m_options.database_path = path; }
    void SetFolderPath(const std::filesystem::path & path) { m_options.SetFolderPath(path); }
    void SetMapFilePath(const std::filesystem::path & path) { m_options.map_file_path = path; }
    void SetModelKeyTagList(const std::string & value);

private:
    void BuildSelectedAtomList(void);
    void RunAtomPositionDumping(void);
    void RunMapValueDumping(void);
    void RunGausEstimatesDumping(void);
    void RunGroupGausEstimatesDumping(void);

};
