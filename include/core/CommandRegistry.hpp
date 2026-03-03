#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "BuiltInCommandCatalog.hpp"
#include "CommandMetadata.hpp"

namespace rhbm_gem {

class CommandBase;

class CommandRegistry
{
public:
    struct CommandInfo
    {
        CommandId id;
        std::string name;
        std::string description;
        std::function<std::unique_ptr<CommandBase>()> factory;
        CommandSurface surface;
        DatabaseUsage database_usage;
        BindingExposure binding_exposure;
    };

    static CommandRegistry & Instance();

    bool RegisterCommand(const CommandDescriptor & descriptor);

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
        CommandRegistry::Instance().RegisterCommand(CommandDescriptor{
            CommandType{}.GetCommandId(),
            std::string_view{name},
            std::string_view{description},
            surface,
            DatabaseUsage::Required,
            BindingExposure::CliOnly,
            [](){ return std::make_unique<CommandType>(); }
        });
    }
};

} // namespace rhbm_gem
