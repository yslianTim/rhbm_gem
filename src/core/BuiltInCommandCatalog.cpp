#include "BuiltInCommandCatalogInternal.hpp"

#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

#include <algorithm>
#include <stdexcept>

namespace rhbm_gem {

namespace {

CommandDescriptor MakeBuiltInDescriptor(
    CommandId id,
    std::string_view name,
    std::string_view description,
    CommandSurface surface,
    std::string_view python_binding_name,
    CommandFactory factory)
{
    if (python_binding_name.empty())
    {
        throw std::runtime_error(
            "Built-in command descriptors must provide a Python binding name.");
    }

    return CommandDescriptor{
        id,
        name,
        description,
        surface,
        python_binding_name,
        factory
    };
}

} // namespace

const std::vector<CommandDescriptor> & BuiltInCommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
        MakeBuiltInDescriptor(
            CommandId::PotentialAnalysis,
            "potential_analysis",
            "Run potential analysis",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::Database
                    | CommonOption::OutputFolder),
            "PotentialAnalysisCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<PotentialAnalysisCommand>();
            }),
        MakeBuiltInDescriptor(
            CommandId::PotentialDisplay,
            "potential_display",
            "Run potential display",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::Database
                    | CommonOption::OutputFolder),
            "PotentialDisplayCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<PotentialDisplayCommand>();
            }),
        MakeBuiltInDescriptor(
            CommandId::ResultDump,
            "result_dump",
            "Run result dump",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::Database
                    | CommonOption::OutputFolder),
            "ResultDumpCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<ResultDumpCommand>();
            }),
        MakeBuiltInDescriptor(
            CommandId::MapSimulation,
            "map_simulation",
            "Run map simulation command",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            "MapSimulationCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<MapSimulationCommand>();
            }),
        MakeBuiltInDescriptor(
            CommandId::MapVisualization,
            "map_visualization",
            "Run map visualization",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            "MapVisualizationCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<MapVisualizationCommand>();
            }),
        MakeBuiltInDescriptor(
            CommandId::PositionEstimation,
            "position_estimation",
            "Run atom position estimation",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            "PositionEstimationCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<PositionEstimationCommand>();
            }),
        MakeBuiltInDescriptor(
            CommandId::ModelTest,
            "model_test",
            "Run HRL model simulation test",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            "HRLModelTestCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<HRLModelTestCommand>();
            })
    };

    return catalog;
}

const CommandDescriptor & FindCommandDescriptor(CommandId command_id)
{
    const auto & catalog{ BuiltInCommandCatalog() };
    const auto iter{
        std::find_if(
            catalog.begin(),
            catalog.end(),
            [command_id](const CommandDescriptor & descriptor)
            {
                return descriptor.id == command_id;
            })
    };
    if (iter == catalog.end())
    {
        throw std::runtime_error("Unknown CommandId in built-in command catalog.");
    }
    return *iter;
}

} // namespace rhbm_gem
