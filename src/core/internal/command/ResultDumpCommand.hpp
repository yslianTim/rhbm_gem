#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

#include "CommandBase.hpp"
#include <rhbm_gem/core/command/CommandEnumClass.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class AtomObject;
class ModelObject;
class MapObject;
struct ResultDumpRequest;
void BindResultDumpRequestOptions(CLI::App * command, ResultDumpRequest & request);

struct ResultDumpCommandOptions
{
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::vector<std::string> model_key_tag_list{};
    std::filesystem::path map_file_path{ "" };
};

class ResultDumpCommand : public CommandBase
{
private:
    ResultDumpCommandOptions m_options{};
    std::string m_map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::shared_ptr<MapObject> m_map_object;

public:
    explicit ResultDumpCommand(CommonOptionProfile profile);
    ~ResultDumpCommand() override = default;
    void ApplyRequest(const ResultDumpRequest & request);

private:
    void SetPrinterChoice(PrinterType value);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetModelKeyTagList(const std::string & value);
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObjectList();
    void RunAtomOutlierDumping();
    void RunAtomPositionDumping();
    void RunMapValueDumping();
    void RunGausEstimatesDumping();
    void RunGroupGausEstimatesDumping();
    void RunDumpWorkflow();
};

} // namespace rhbm_gem
