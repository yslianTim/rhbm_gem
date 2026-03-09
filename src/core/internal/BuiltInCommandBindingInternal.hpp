#pragma once

#include <string_view>

namespace rhbm_gem {

#define RHBM_GEM_BUILTIN_COMMAND(COMMAND_TYPE, CLI_NAME, DESCRIPTION, PYTHON_BINDING_NAME) \
    class COMMAND_TYPE;
#include "internal/BuiltInCommandList.def"
#undef RHBM_GEM_BUILTIN_COMMAND

template <typename CommandType>
struct BuiltInCommandBindingName;

#define RHBM_GEM_BUILTIN_COMMAND(COMMAND_TYPE, CLI_NAME, DESCRIPTION, PYTHON_BINDING_NAME) \
    template <>                                                                        \
    struct BuiltInCommandBindingName<COMMAND_TYPE>                                    \
    {                                                                                  \
        inline static constexpr std::string_view value{ PYTHON_BINDING_NAME };        \
    };
#include "internal/BuiltInCommandList.def"
#undef RHBM_GEM_BUILTIN_COMMAND

} // namespace rhbm_gem
