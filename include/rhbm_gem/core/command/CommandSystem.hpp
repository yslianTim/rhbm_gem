#pragma once

#include <vector>

#include <rhbm_gem/core/command/CommandList.hpp>

namespace rhbm_gem {

const std::vector<CommandInfo> & ListCommands();
int RunCommandCLI(int argc, char * argv[]);

} // namespace rhbm_gem
