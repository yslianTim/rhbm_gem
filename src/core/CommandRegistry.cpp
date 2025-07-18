#include "CommandRegistry.hpp"
#include "CommandBase.hpp"

CommandRegistry & CommandRegistry::Instance()
{
    static CommandRegistry instance;
    return instance;
}

bool CommandRegistry::RegisterCommand(const std::string & name,
                                      const std::string & description,
                                      FactoryFunc factory)
{
    m_commands.push_back({name, description, std::move(factory)});
    return true;
}
