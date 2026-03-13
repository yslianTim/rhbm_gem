#pragma once

#include <cstdint>

namespace rhbm_gem {

enum class CommandId
{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    COMMAND_ID,
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
};

enum class CommonOption : std::uint8_t
{
    Threading = 0,
    Verbose = 1,
    Database = 2,
    OutputFolder = 3
};

using CommonOptionMask = std::uint32_t;

enum class CommonOptionProfile : std::uint8_t
{
    FileWorkflow = 0,
    DatabaseWorkflow = 1
};

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

constexpr CommonOptionMask CommonOptionMaskForProfile(CommonOptionProfile profile)
{
    switch (profile)
    {
    case CommonOptionProfile::FileWorkflow:
        return CommonOption::Threading | CommonOption::Verbose | CommonOption::OutputFolder;
    case CommonOptionProfile::DatabaseWorkflow:
        return CommonOption::Threading
            | CommonOption::Verbose
            | CommonOption::Database
            | CommonOption::OutputFolder;
    default:
        return CommonOption::Threading | CommonOption::Verbose | CommonOption::OutputFolder;
    }
}

} // namespace rhbm_gem
