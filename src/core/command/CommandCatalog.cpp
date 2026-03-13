#include "internal/CommandCatalog.hpp"
#include "internal/CommandRuntimeRegistry.hpp"

namespace rhbm_gem {

const std::vector<CommandDescriptor> & CommandCatalog()
{
    static const std::vector<CommandDescriptor> catalog{
    #include "internal/CommandCatalogEntries.generated.inc"
    };

    return catalog;
}

} // namespace rhbm_gem
