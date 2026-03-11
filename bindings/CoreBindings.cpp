#include "BindingHelpers.hpp"

namespace rhbm_gem {
#define RHBM_GEM_BUILTIN_COMMAND(COMMAND_TYPE, CLI_NAME, DESCRIPTION) \
    class COMMAND_TYPE;
#include "internal/BuiltInCommandList.def"
#undef RHBM_GEM_BUILTIN_COMMAND
} // namespace rhbm_gem

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommonTypes(module);
#define RHBM_GEM_BUILTIN_COMMAND(COMMAND_TYPE, CLI_NAME, DESCRIPTION) \
    rhbm_gem::bindings::BindCommand<rhbm_gem::COMMAND_TYPE>(module);
#include "internal/BuiltInCommandList.def"
#undef RHBM_GEM_BUILTIN_COMMAND
}
