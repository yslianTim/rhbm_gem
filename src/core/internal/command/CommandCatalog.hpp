#pragma once

#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandMetadata.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

struct CommandDescriptor
{
    CommandId id;
    std::string_view name;
    std::string_view description;
    CommonOptionProfile profile;
};

constexpr CommonOptionMask CommonOptionsForCommand(const CommandDescriptor & descriptor)
{
    return CommonOptionMaskForProfile(descriptor.profile);
}

const std::vector<CommandDescriptor> & CommandCatalog();
void RegisterCommandSubcommands(CLI::App & app);

} // namespace rhbm_gem
