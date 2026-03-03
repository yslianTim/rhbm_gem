#pragma once

#include <cstdint>

namespace rhbm_gem {

enum class CommonOption : std::uint8_t
{
    Threading = 0,
    Verbose = 1,
    Database = 2,
    OutputFolder = 3
};

using CommonOptionMask = std::uint32_t;

constexpr CommonOptionMask ToMask(CommonOption option)
{
    return static_cast<CommonOptionMask>(1u << static_cast<std::uint8_t>(option));
}

constexpr CommonOptionMask operator|(CommonOption lhs, CommonOption rhs)
{
    return ToMask(lhs) | ToMask(rhs);
}

constexpr CommonOptionMask operator|(CommonOptionMask lhs, CommonOption rhs)
{
    return lhs | ToMask(rhs);
}

constexpr CommonOptionMask operator|(CommonOption lhs, CommonOptionMask rhs)
{
    return ToMask(lhs) | rhs;
}

constexpr bool HasCommonOption(CommonOptionMask mask, CommonOption option)
{
    return (mask & ToMask(option)) != 0u;
}

struct CommandSurface
{
    CommonOptionMask common_options{
        CommonOption::Threading
        | CommonOption::Verbose
        | CommonOption::Database
        | CommonOption::OutputFolder
    };
    CommonOptionMask deprecated_hidden_options{ 0u };
    bool requires_database_manager{ true };
    bool uses_output_folder{ true };
    bool exposed_to_python{ false };
};

constexpr CommandSurface MakeCommandSurface(
    CommonOptionMask common_options,
    CommonOptionMask deprecated_hidden_options,
    bool requires_database_manager,
    bool uses_output_folder,
    bool exposed_to_python)
{
    return CommandSurface{
        common_options,
        deprecated_hidden_options,
        requires_database_manager,
        uses_output_folder,
        exposed_to_python
    };
}

} // namespace rhbm_gem
