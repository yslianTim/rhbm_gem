#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

#include "CommandBase.hpp"
#include <rhbm_gem/core/command/OptionEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;
class MapObject;
struct ResultDumpRequest;

struct ResultDumpCommandOptions : public CommandOptions
{
    PrinterType printer_choice{ PrinterType::GAUS_ESTIMATES };
    std::vector<std::string> model_key_tag_list{};
    std::filesystem::path map_file_path{ "" };
};

class ResultDumpCommand
    : public CommandWithProfileOptions<
          ResultDumpCommandOptions,
          CommandId::ResultDump,
          CommonOptionProfile::DatabaseWorkflow>
{
public:
    using Options = ResultDumpCommandOptions;

private:
    std::string m_map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::shared_ptr<MapObject> m_map_object;

public:
    explicit ResultDumpCommand();
    ~ResultDumpCommand() override = default;
    void ApplyRequest(const ResultDumpRequest & request);

private:
    void SetPrinterChoice(PrinterType value);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetModelKeyTagList(const std::string & value);
    void ValidateOptions() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
};

} // namespace rhbm_gem
