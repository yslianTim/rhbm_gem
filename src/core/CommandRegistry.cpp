#include "CommandRegistry.hpp"
#include "CommandBase.hpp"
#include "Logger.hpp"

#include <algorithm>

namespace rhbm_gem {

CommandRegistry & CommandRegistry::Instance()
{
    static CommandRegistry instance;
    return instance;
}

bool CommandRegistry::RegisterCommand(
    const std::string & name,
    const std::string & description,
    std::function<std::unique_ptr<CommandBase>()> factory)
{
    const auto duplicate_iter{
        std::find_if(
            m_commands.begin(),
            m_commands.end(),
            [&name](const CommandRegistry::CommandInfo & info)
            {
                return info.name == name;
            })
    };
    if (duplicate_iter != m_commands.end())
    {
        Logger::Log(LogLevel::Error, "Command [" + name + "] is already registered.");
        return false;
    }

    m_commands.push_back(CommandInfo{ name, description, std::move(factory) });
    return true;
}

} // namespace rhbm_gem
