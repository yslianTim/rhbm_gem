#pragma once

#include <memory>
#include <string>

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
    CLI::App * m_potential_analysis_cmd;
    CLI::App * m_potential_display_cmd;
    CLI::App * m_result_dump_cmd;
    CLI::App * m_map_simulation_cmd;

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
    void RegisterCommands(void);
    void RegisterPotentialAnalysisCommand(void);
    void RegisterPotentialDisplayCommand(void);
    void RegisterResultDumpCommand(void);
    void RegisterMapSimulationCommand(void);
    void RegisterGlobalOptions(CLI::App * command);

};
