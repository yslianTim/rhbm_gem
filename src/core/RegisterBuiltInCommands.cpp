#include "RegisterBuiltInCommands.hpp"

#include "CommandRegistry.hpp"
#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

#include <mutex>
#include <utility>

namespace rhbm_gem {

namespace {

template <typename CommandType>
void RegisterBuiltInCommand(CommandRegistry & registry)
{
    registry.RegisterCommand(
        std::string(CommandType::CommandName()),
        std::string(CommandType::CommandDescription()),
        []() { return std::make_unique<CommandType>(); },
        CommandType::StaticCommandSurface());
}

} // namespace

void RegisterBuiltInCommands()
{
    static std::once_flag register_once;
    std::call_once(register_once, []()
    {
        auto & registry{ CommandRegistry::Instance() };
        RegisterBuiltInCommand<PotentialAnalysisCommand>(registry);
        RegisterBuiltInCommand<PotentialDisplayCommand>(registry);
        RegisterBuiltInCommand<ResultDumpCommand>(registry);
        RegisterBuiltInCommand<MapSimulationCommand>(registry);
        RegisterBuiltInCommand<MapVisualizationCommand>(registry);
        RegisterBuiltInCommand<PositionEstimationCommand>(registry);
        RegisterBuiltInCommand<HRLModelTestCommand>(registry);
    });
}

} // namespace rhbm_gem
