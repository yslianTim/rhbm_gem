#include "BuiltInCommandCatalog.hpp"

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

const std::vector<CommandDescriptor> & BuiltInCommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
        CommandDescriptor{
            CommandId::PotentialAnalysis,
            "potential_analysis",
            "Run potential analysis",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::Database
                    | CommonOption::OutputFolder),
            DatabaseUsage::Required,
            BindingExposure::PythonPublic,
            "PotentialAnalysisCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<PotentialAnalysisCommand>();
            }
        },
        CommandDescriptor{
            CommandId::PotentialDisplay,
            "potential_display",
            "Run potential display",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::Database
                    | CommonOption::OutputFolder),
            DatabaseUsage::Required,
            BindingExposure::PythonPublic,
            "PotentialDisplayCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<PotentialDisplayCommand>();
            }
        },
        CommandDescriptor{
            CommandId::ResultDump,
            "result_dump",
            "Run result dump",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::Database
                    | CommonOption::OutputFolder),
            DatabaseUsage::Required,
            BindingExposure::PythonPublic,
            "ResultDumpCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<ResultDumpCommand>();
            }
        },
        CommandDescriptor{
            CommandId::MapSimulation,
            "map_simulation",
            "Run map simulation command",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            DatabaseUsage::NotUsed,
            BindingExposure::PythonPublic,
            "MapSimulationCommand",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<MapSimulationCommand>();
            }
        },
        CommandDescriptor{
            CommandId::MapVisualization,
            "map_visualization",
            "Run map visualization",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            DatabaseUsage::NotUsed,
            BindingExposure::CliOnly,
            "",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<MapVisualizationCommand>();
            }
        },
        CommandDescriptor{
            CommandId::PositionEstimation,
            "position_estimation",
            "Run atom position estimation",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            DatabaseUsage::NotUsed,
            BindingExposure::CliOnly,
            "",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<PositionEstimationCommand>();
            }
        },
        CommandDescriptor{
            CommandId::ModelTest,
            "model_test",
            "Run HRL model simulation test",
            MakeCommandSurface(
                CommonOption::Threading
                    | CommonOption::Verbose
                    | CommonOption::OutputFolder),
            DatabaseUsage::NotUsed,
            BindingExposure::CliOnly,
            "",
            []() -> std::unique_ptr<CommandBase>
            {
                return std::make_unique<HRLModelTestCommand>();
            }
        }
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
