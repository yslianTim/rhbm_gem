#include "internal/CommandCatalog.hpp"
#include "internal/CommandCatalogIncludes.generated.inc"

namespace rhbm_gem {

namespace {

template <typename CommandType>
CommandDescriptor MakeCommandDescriptor(
    std::string_view name,
    std::string_view description)
{
    return CommandDescriptor{
        CommandType::kCommandId,
        name,
        description,
        CommandType::kCommonOptions,
        []() -> std::unique_ptr<CommandBase>
        {
            return std::make_unique<CommandType>();
        }
    };
}

} // namespace

const std::vector<CommandDescriptor> & CommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
    #include "internal/CommandCatalogEntries.generated.inc"
    };

    return catalog;
}

} // namespace rhbm_gem
