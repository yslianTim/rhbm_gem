#include "BindingHelpers.hpp"

namespace rhbm_gem {
#define RHBM_GEM_COMMAND(COMMAND_ID, COMMAND_TYPE, CLI_NAME, DESCRIPTION, PROFILE, PYTHON_BINDING_NAME) \
    class COMMAND_TYPE;
#include "internal/CommandList.def"
#undef RHBM_GEM_COMMAND
} // namespace rhbm_gem

PYBIND11_MODULE(rhbm_gem_module, module)
{
    rhbm_gem::bindings::BindCommonTypes(module);
#define RHBM_GEM_COMMAND(COMMAND_ID, COMMAND_TYPE, CLI_NAME, DESCRIPTION, PROFILE, PYTHON_BINDING_NAME) \
    rhbm_gem::bindings::BindCommand<rhbm_gem::COMMAND_TYPE>(module);
#include "internal/CommandList.def"
#undef RHBM_GEM_COMMAND
}
