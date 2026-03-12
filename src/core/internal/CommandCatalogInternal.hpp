#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandMetadata.hpp>

namespace rhbm_gem {

class CommandBase;
class DataIoServices;
using CommandFactory = std::unique_ptr<CommandBase>(*)(const DataIoServices &);

struct CommandDescriptor
{
    CommandId id;
    std::string_view name;
    std::string_view description;
    CommonOptionMask common_options;
    std::string_view python_binding_name;
    CommandFactory factory;
};

const std::vector<CommandDescriptor> & CommandCatalog();
std::string_view CommandPythonBindingName(CommandId id);

constexpr bool UsesDatabaseAtRuntime(CommonOptionMask common_options)
{
    return HasCommonOption(common_options, CommonOption::Database);
}

constexpr bool UsesOutputFolder(CommonOptionMask common_options)
{
    return HasCommonOption(common_options, CommonOption::OutputFolder);
}

} // namespace rhbm_gem
