#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "OptionEnumClass.hpp"

namespace rhbm_gem {

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
    ResultDumpCommand();
    ~ResultDumpCommand();
    bool Execute() override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    const CommandOptions & GetOptions() const override { return m_options; }
    CommandOptions & GetOptions() override { return m_options; }

    void SetPrinterChoice(PrinterType value);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetModelKeyTagList(const std::string & value);

private:
    bool BuildDataObjectList();
    void RunResultDump();
    void RunAtomOutlierDumping();
    void RunAtomPositionDumping();
    void RunMapValueDumping();
    void RunGausEstimatesDumping();
    void RunGroupGausEstimatesDumping();

};

} // namespace rhbm_gem
