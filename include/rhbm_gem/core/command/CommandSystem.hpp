#pragma once

#include <vector>

#include <rhbm_gem/core/command/CommandTypes.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

const std::vector<CommandInfo> & ListCommands();
void ConfigureCommandCLI(::CLI::App & app);

#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
CommandResult Run##COMMAND_ID(const COMMAND_ID##Request & request);
#include <rhbm_gem/core/command/CommandManifest.def>
#undef RHBM_GEM_COMMAND

} // namespace rhbm_gem
