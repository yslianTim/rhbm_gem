#pragma once

#include <memory>
#include <string>
#include <filesystem>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

class MapObject;

class PositionEstimationCommand : public CommandBase
{
public:
    struct Options : public CommandOptions
    {
        int iteration_count{ 15 };
        int knn_size{ 20 };
        double alpha{ 2.0 };
        std::filesystem::path map_file_path;
        std::string saved_key_tag{"map"};
    };

private:
    Options m_options{};

public:
    PositionEstimationCommand(void);
    ~PositionEstimationCommand() = default;
    bool Execute(void) override;
    bool ValidateOptions(void) const override;
    void RegisterCLIOptionsExtend(CLI::App * cmd) override;
    const CommandOptions & GetOptions(void) const override { return m_options; }
    CommandOptions & GetOptions(void) override { return m_options; }

    void SetDatabasePath(const std::filesystem::path & path) { m_options.database_path = path; }
    void SetMapFilePath(const std::filesystem::path & path) { m_options.map_file_path = path; }
    void SetSavedKeyTag(const std::string & tag) { m_options.saved_key_tag = tag; }
    void SetThreadSize(int value) { m_options.thread_size = value; }

private:
    void RunMapValueConvergence(MapObject * map_object);

};
