#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

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
    MapVisualizationCommand(void);
    ~MapVisualizationCommand();
    bool Execute(void) override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetModelFilePath(const std::filesystem::path & path);
    void SetMapFilePath(const std::filesystem::path & path);
    void SetAtomSerialID(int value) { m_options.atom_serial_id = value; }
    void SetSamplingSize(int value);
    void SetWindowSize(double value);

private:
    bool BuildDataObject(void);
    void RunMapObjectPreprocessing(void);
    void RunModelObjectPreprocessing(void);
    void RunAtomMapValueSampling(void);

};
