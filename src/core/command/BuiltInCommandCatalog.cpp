#include "internal/BuiltInCommandCatalogInternal.hpp"

#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/command/MapSimulationCommand.hpp>
#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>
#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/command/PotentialDisplayCommand.hpp>
#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

#include <stdexcept>

namespace rhbm_gem {

namespace {

template <typename CommandType>
CommandDescriptor MakeBuiltInDescriptor(
    std::string_view name,
    std::string_view description,
    std::string_view python_binding_name)
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
        []() -> std::unique_ptr<CommandBase>
        {
            return std::make_unique<CommandType>();
        }
    };
}

} // namespace

const std::vector<CommandDescriptor> & BuiltInCommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
    #define RHBM_GEM_BUILTIN_COMMAND(COMMAND_TYPE, CLI_NAME, DESCRIPTION, PYTHON_BINDING_NAME) \
        MakeBuiltInDescriptor<COMMAND_TYPE>(CLI_NAME, DESCRIPTION, PYTHON_BINDING_NAME),
    #include "internal/BuiltInCommandList.def"
    #undef RHBM_GEM_BUILTIN_COMMAND
    };

    return catalog;
}

} // namespace rhbm_gem
