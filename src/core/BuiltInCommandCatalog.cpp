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

template <typename CommandType>
CommandDescriptor MakeBuiltInDescriptor(
    std::string_view name,
    std::string_view description,
    std::string_view python_binding_name,
    CommandFactory factory = []() -> std::unique_ptr<CommandBase>
    {
        return std::make_unique<CommandType>();
    })
{
    if (python_binding_name.empty())
    {
        throw std::runtime_error(
            "Built-in command descriptors must provide a Python binding name.");
    }

    return CommandDescriptor{
        CommandType::kCommandId,
        name,
        description,
        CommandType::kCommonOptions,
        python_binding_name,
        factory
    };
}

} // namespace

const std::vector<CommandDescriptor> & BuiltInCommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
        MakeBuiltInDescriptor<PotentialAnalysisCommand>(
            "potential_analysis",
            "Run potential analysis",
            "PotentialAnalysisCommand"),
        MakeBuiltInDescriptor<PotentialDisplayCommand>(
            "potential_display",
            "Run potential display",
            "PotentialDisplayCommand"),
        MakeBuiltInDescriptor<ResultDumpCommand>(
            "result_dump",
            "Run result dump",
            "ResultDumpCommand"),
        MakeBuiltInDescriptor<MapSimulationCommand>(
            "map_simulation",
            "Run map simulation command",
            "MapSimulationCommand"),
        MakeBuiltInDescriptor<MapVisualizationCommand>(
            "map_visualization",
            "Run map visualization",
            "MapVisualizationCommand"),
        MakeBuiltInDescriptor<PositionEstimationCommand>(
            "position_estimation",
            "Run atom position estimation",
            "PositionEstimationCommand"),
        MakeBuiltInDescriptor<HRLModelTestCommand>(
            "model_test",
            "Run HRL model simulation test",
            "HRLModelTestCommand")
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
