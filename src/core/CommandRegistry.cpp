#include "CommandRegistry.hpp"
#include "CommandBase.hpp"
#include "Logger.hpp"

CommandRegistry & CommandRegistry::Instance(void)
{
    static CommandRegistry instance;
    return instance;
}

bool CommandRegistry::RegisterCommand(
    const std::string & name,
    const std::string & description,
    FactoryFunc factory)
{
    if (m_commands.find(name) != m_commands.end())
    {
        Logger::Log(LogLevel::Error, "Command [" + name + "] is already registered.");
        return false;
    }

    m_commands.emplace(name, CommandInfo{name, description, std::move(factory)});
    return true;
}
