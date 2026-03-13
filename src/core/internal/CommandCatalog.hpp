#pragma once

#include <functional>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandMetadata.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

using CommandRunner = std::function<ExecutionReport()>;
using CommandRuntimeBinder = CommandRunner (*)(CLI::App *);

struct CommandDescriptor
{
    CommandId id;
    std::string_view name;
    std::string_view description;
    CommonOptionMask common_options;
    CommandRuntimeBinder bind_runtime;
};

const std::vector<CommandDescriptor> & CommandCatalog();

} // namespace rhbm_gem
