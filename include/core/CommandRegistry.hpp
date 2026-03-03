#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace rhbm_gem {

class CommandBase;

class CommandRegistry
{
public:
    struct CommandInfo
    {
        std::string name;
        std::string description;
        std::function<std::unique_ptr<CommandBase>()> factory;
    };

    static CommandRegistry & Instance();

    bool RegisterCommand(
        const std::string & name,
        const std::string & description,
        std::function<std::unique_ptr<CommandBase>()> factory
    );

    const std::vector<CommandInfo> &
    GetCommands() const { return m_commands; }

private:
    std::vector<CommandInfo> m_commands;
};

template <typename CommandType>
class CommandRegistrar
{
public:
    CommandRegistrar(const std::string & name, const std::string & description)
    {
        CommandRegistry::Instance().RegisterCommand(
            name,
            description,
            [](){ return std::make_unique<CommandType>(); }
        );
    }
};

} // namespace rhbm_gem
