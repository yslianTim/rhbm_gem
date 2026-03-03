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
    const CommandDescriptor & descriptor)
{
    const auto duplicate_iter{
        std::find_if(
            m_commands.begin(),
            m_commands.end(),
            [&descriptor](const CommandRegistry::CommandInfo & info)
            {
                return info.name == descriptor.name;
            })
    };
    if (duplicate_iter != m_commands.end())
    {
        Logger::Log(
            LogLevel::Error,
            "Command [" + std::string(descriptor.name) + "] is already registered.");
        return false;
    }

    m_commands.push_back(CommandInfo{
        descriptor.id,
        std::string(descriptor.name),
        std::string(descriptor.description),
        descriptor.factory,
        descriptor.surface,
        descriptor.database_usage,
        descriptor.binding_exposure,
        std::string(descriptor.python_binding_name)
    });
    return true;
}

} // namespace rhbm_gem
