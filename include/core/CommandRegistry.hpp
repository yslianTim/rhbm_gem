#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

class CommandBase;

class CommandRegistry
{
public:
    using FactoryFunc = std::function<std::unique_ptr<CommandBase>()>;
    struct CommandInfo
    {
        std::string name;
        std::string description;
        FactoryFunc factory;
    };

    static CommandRegistry & Instance(void);

    bool RegisterCommand(const std::string & name,
                         const std::string & description,
                         FactoryFunc factory);

    const std::unordered_map<std::string, CommandInfo> & GetCommands(void) const { return m_commands; }

private:
    std::unordered_map<std::string, CommandInfo> m_commands;
};

#define REGISTER_COMMAND(command_type, name, description) \
    namespace { \
        const bool registered_##command_type = \
            CommandRegistry::Instance().RegisterCommand(name, description, [](){ \
                return std::make_unique<command_type>(); \
            }); \
    }
