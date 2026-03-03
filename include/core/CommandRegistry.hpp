#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "CommandMetadata.hpp"

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
        CommandSurface surface;
    };

    static CommandRegistry & Instance();

    bool RegisterCommand(
        const std::string & name,
        const std::string & description,
        std::function<std::unique_ptr<CommandBase>()> factory,
        CommandSurface surface
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
    CommandRegistrar(
        const std::string & name,
        const std::string & description,
        CommandSurface surface)
    {
        CommandRegistry::Instance().RegisterCommand(
            name,
            description,
            [](){ return std::make_unique<CommandType>(); },
            surface
        );
    }
};

} // namespace rhbm_gem
