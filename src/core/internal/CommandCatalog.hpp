#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandMetadata.hpp>

namespace rhbm_gem {

class CommandBase;
using CommandFactory = std::unique_ptr<CommandBase>(*)();

struct CommandDescriptor
{
    CommandId id;
    std::string_view name;
    std::string_view description;
    CommonOptionMask common_options;
    CommandFactory factory;
};

const std::vector<CommandDescriptor> & CommandCatalog();

} // namespace rhbm_gem
