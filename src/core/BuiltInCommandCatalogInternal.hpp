#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include "CommandId.hpp"
#include "CommandMetadata.hpp"

namespace rhbm_gem {

class CommandBase;
using CommandFactory = std::unique_ptr<CommandBase>(*)();

struct CommandDescriptor
{
    CommandId id;
    std::string_view name;
    std::string_view description;
    CommandSurface surface;
    std::string_view python_binding_name;
    CommandFactory factory;
};

const std::vector<CommandDescriptor> & BuiltInCommandCatalog();
const CommandDescriptor & FindCommandDescriptor(CommandId command_id);

constexpr bool UsesDatabaseAtRuntime(const CommandSurface & surface)
{
    return HasCommonOption(surface.common_options, CommonOption::Database);
}

constexpr bool UsesOutputFolder(const CommandSurface & surface)
{
    return HasCommonOption(surface.common_options, CommonOption::OutputFolder);
}

} // namespace rhbm_gem
