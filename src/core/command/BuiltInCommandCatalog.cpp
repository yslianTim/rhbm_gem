#include "internal/BuiltInCommandCatalogInternal.hpp"
#include "internal/BuiltInCommandCatalogIncludes.generated.inc"

#include <rhbm_gem/data/io/DataIoServices.hpp>

#include <stdexcept>

namespace rhbm_gem {

namespace {

template <typename CommandType>
CommandDescriptor MakeBuiltInDescriptor(
    std::string_view python_binding_name,
    std::string_view name,
    std::string_view description)
{
    if (python_binding_name.empty())
    {
        throw std::runtime_error(
            "Built-in command descriptors must provide a Python binding name.");
    }

    return CommandDescriptor{
        CommandType::kCommandId,
        name,
        description,
        CommandType::kCommonOptions,
        python_binding_name,
        [](const DataIoServices & data_io_services) -> std::unique_ptr<CommandBase>
        {
            return std::make_unique<CommandType>(data_io_services);
        }
    };
}

} // namespace

const std::vector<CommandDescriptor> & BuiltInCommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
    #include "internal/BuiltInCommandCatalogEntries.generated.inc"
    };

    return catalog;
}

std::string_view BuiltInPythonBindingName(CommandId id)
{
    for (const auto & descriptor : BuiltInCommandCatalog())
    {
        if (descriptor.id == id)
        {
            return descriptor.python_binding_name;
        }
    }

    throw std::runtime_error(
        "BuiltInPythonBindingName lookup failed: command id not found.");
}

} // namespace rhbm_gem
