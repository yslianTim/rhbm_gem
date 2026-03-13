#pragma once

#include <functional>

#include <rhbm_gem/core/command/CommandApi.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

using CommandRunner = std::function<ExecutionReport()>;
using CommandRuntimeBinder = CommandRunner (*)(CLI::App *);

CommandRunner BindPotentialAnalysisRuntime(CLI::App * command);
CommandRunner BindPotentialDisplayRuntime(CLI::App * command);
CommandRunner BindResultDumpRuntime(CLI::App * command);
CommandRunner BindMapSimulationRuntime(CLI::App * command);
CommandRunner BindMapVisualizationRuntime(CLI::App * command);
CommandRunner BindPositionEstimationRuntime(CLI::App * command);
CommandRunner BindModelTestRuntime(CLI::App * command);

} // namespace rhbm_gem
