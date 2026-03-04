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

struct ResultDumpCommandOptions : public CommandOptions
{
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::vector<std::string> model_key_tag_list{};
    std::filesystem::path map_file_path{ "" };
};

class ResultDumpCommand
    : public CommandWithOptions<
          ResultDumpCommandOptions,
          CommandId::ResultDump,
          CommonOption::Threading
              | CommonOption::Verbose
              | CommonOption::Database
              | CommonOption::OutputFolder>
{
public:
    using Options = ResultDumpCommandOptions;

private:
    std::string m_map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::shared_ptr<MapObject> m_map_object;

public:
    ResultDumpCommand();
    ~ResultDumpCommand();
    void SetPrinterChoice(PrinterType value);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetModelKeyTagList(const std::string & value);

private:
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObjectList();
    void RunResultDump();
    void RunAtomOutlierDumping();
    void RunAtomPositionDumping();
    void RunMapValueDumping();
    void RunGausEstimatesDumping();
    void RunGroupGausEstimatesDumping();

};

} // namespace rhbm_gem
