#pragma once

#include <memory>
#include <string>
#include <filesystem>

#include "CommandBase.hpp"

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;
struct MapVisualizationRequest;
void BindMapVisualizationRequestOptions(CLI::App * command, MapVisualizationRequest & request);

struct MapVisualizationCommandOptions
{
    int atom_serial_id{ 1 };
    int sampling_size{ 100 };
    double window_size{ 5.0 };
    std::filesystem::path model_file_path;
    std::filesystem::path map_file_path;
};

class MapVisualizationCommand : public CommandBase
{
private:
    MapVisualizationCommandOptions m_options{};
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    explicit MapVisualizationCommand(CommonOptionProfile profile);
    ~MapVisualizationCommand() override = default;
    void ApplyRequest(const MapVisualizationRequest & request);

private:
    void SetModelFilePath(const std::filesystem::path & path);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetAtomSerialID(int value);
    void SetSamplingSize(int value);
    void SetWindowSize(double value);
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunMapObjectPreprocessing();
    void RunModelObjectPreprocessing();
    bool RunAtomMapValueSampling();
    std::filesystem::path BuildOutputFilePath() const;

};

} // namespace rhbm_gem
