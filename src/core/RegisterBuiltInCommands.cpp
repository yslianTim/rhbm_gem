#include "RegisterBuiltInCommands.hpp"

#include "BuiltInCommandCatalog.hpp"
#include "CommandRegistry.hpp"

#include <mutex>

namespace rhbm_gem {

void RegisterBuiltInCommands()
{
    static std::once_flag register_once;
    std::call_once(register_once, []()
    {
        auto & registry{ CommandRegistry::Instance() };
        for (const auto & descriptor : BuiltInCommandCatalog())
        {
            registry.RegisterCommand(descriptor);
        }
    });
}

} // namespace rhbm_gem
