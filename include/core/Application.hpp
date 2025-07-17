#pragma once

#include <memory>
#include <string>
#include <functional>
#include <unordered_map>

#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "GlobalOptions.hpp"

namespace CLI
{
    class App;
}

class CommandBase;

class Application
{
    CLI::App & m_cli_app;
    std::unordered_map<std::string, std::function<std::unique_ptr<CommandBase>()>> m_command_map;
    PotentialAnalysisCommand::Options m_potential_analysis_options{};
    PotentialDisplayCommand::Options  m_potential_display_options{};
    ResultDumpCommand::Options        m_result_dump_options{};
    MapSimulationCommand::Options     m_map_simulation_options{};
    GlobalOptions m_global_options{};

public:
    Application(CLI::App & app);
    ~Application() = default;
    void Run(void);

private:
    std::unique_ptr<CommandBase> CreateCommand(void);
    void RegisterAllCommands(void);
    template<typename Type>
    void RegisterCommand(const std::string & name,
                         const std::string & description,
                         typename Type::Options & options);
    void RegisterGlobalOptions(CLI::App * command);

};
