#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>

#include "CommandBase.hpp"

namespace rhbm_gem {

class AtomObject;
class ModelObject;
class MapObject;

class ResultDumpCommand : public CommandWithRequest<ResultDumpRequest>
{
private:
    std::string m_map_key_tag;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<std::shared_ptr<ModelObject>> m_model_object_list;
    std::shared_ptr<MapObject> m_map_object;

public:
    ResultDumpCommand();
    ~ResultDumpCommand() override = default;

private:
    void NormalizeRequest() override;
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
