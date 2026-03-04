#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

namespace rhbm_gem {

class ModelObject;
class MapObject;
class AtomObject;

struct MapVisualizationCommandOptions : public CommandOptions
{
    int atom_serial_id{ 1 };
    int sampling_size{ 100 };
    double window_size{ 5.0 };
    std::filesystem::path model_file_path;
    std::filesystem::path map_file_path;
};

class MapVisualizationCommand
    : public CommandWithOptions<
          MapVisualizationCommandOptions,
          CommandId::MapVisualization,
          CommonOption::Threading
              | CommonOption::Verbose
              | CommonOption::OutputFolder>
{
public:
    using Options = MapVisualizationCommandOptions;

private:
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    MapVisualizationCommand();
    ~MapVisualizationCommand();
    void SetModelFilePath(const std::filesystem::path & path);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetAtomSerialID(int value);
    void SetSamplingSize(int value);
    void SetWindowSize(double value);

private:
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    void ResetRuntimeState() override;
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunMapObjectPreprocessing();
    void RunModelObjectPreprocessing();
    bool RunAtomMapValueSampling();
    std::filesystem::path BuildOutputFilePath() const;

};

} // namespace rhbm_gem
