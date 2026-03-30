#pragma once

#include <memory>
#include <string>
#include <filesystem>

#include "CommandBase.hpp"

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;

class MapVisualizationCommand : public CommandWithRequest<MapVisualizationRequest>
{
private:
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    explicit MapVisualizationCommand(CommonOptionProfile profile);
    ~MapVisualizationCommand() override = default;

private:
    void NormalizeRequest() override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunMapObjectPreprocessing();
    void RunModelObjectPreprocessing();
    bool RunAtomMapValueSampling();
    std::filesystem::path BuildOutputFilePath() const;

};

} // namespace rhbm_gem
