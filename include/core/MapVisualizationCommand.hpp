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

class MapVisualizationCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        int atom_serial_id{ 1 };
        int sampling_size{ 100 };
        double window_size{ 5.0 };
        std::filesystem::path model_file_path;
        std::filesystem::path map_file_path;
    };

private:
    Options m_options;
    std::string m_model_key_tag, m_map_key_tag;
    std::shared_ptr<MapObject> m_map_object;
    std::shared_ptr<ModelObject> m_model_object;

public:
    MapVisualizationCommand();
    ~MapVisualizationCommand();
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    CommandId GetCommandId() const override { return CommandId::MapVisualization; }
    void ResetRuntimeState() override;
    const CommandOptions & GetOptions() const override { return m_options; }
    CommandOptions & GetOptions() override { return m_options; }

    void SetModelFilePath(const std::filesystem::path & path);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetAtomSerialID(int value);
    void SetSamplingSize(int value);
    void SetWindowSize(double value);

private:
    bool ExecuteImpl() override;
    bool BuildDataObject();
    void RunMapObjectPreprocessing();
    void RunModelObjectPreprocessing();
    bool RunAtomMapValueSampling();
    std::filesystem::path BuildOutputFilePath() const;

};

} // namespace rhbm_gem
