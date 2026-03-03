#pragma once

#include <functional>
#include <memory>
#include <string>
#include <utility>
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
        std::string python_binding_name;
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
    explicit CommandRegistrar(CommandDescriptor descriptor)
    {
        CommandRegistry::Instance().RegisterCommand(std::move(descriptor));
    }
};

} // namespace rhbm_gem
