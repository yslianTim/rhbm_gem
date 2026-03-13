#pragma once

#include "CommandCatalog.hpp"

namespace rhbm_gem {

CommandRunner BindPotentialAnalysisRuntime(CLI::App * command);
CommandRunner BindPotentialDisplayRuntime(CLI::App * command);
CommandRunner BindResultDumpRuntime(CLI::App * command);
CommandRunner BindMapSimulationRuntime(CLI::App * command);
CommandRunner BindMapVisualizationRuntime(CLI::App * command);
CommandRunner BindPositionEstimationRuntime(CLI::App * command);
CommandRunner BindModelTestRuntime(CLI::App * command);

} // namespace rhbm_gem
