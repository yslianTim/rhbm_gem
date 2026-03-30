#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <rhbm_gem/utils/domain/Logger.hpp>

namespace rhbm_gem {

// Default paths
inline std::filesystem::path GetDefaultDataRootPath()
{
    if (const char * configured_root{ std::getenv("RHBM_GEM_DATA_DIR") };
        configured_root != nullptr && configured_root[0] != '\0')
    {
        return std::filesystem::path(configured_root);
    }

    if (const char * home{ std::getenv("HOME") };
        home != nullptr && home[0] != '\0')
    {
        return std::filesystem::path(home) / ".rhbmgem" / "data";
    }

    return std::filesystem::path(".rhbmgem") / "data";
}

inline std::filesystem::path GetDefaultDatabasePath()
{
    return GetDefaultDataRootPath() / "database.sqlite";
}

// Command catalog
enum class CommandId
{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
    COMMAND_ID,
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
};

// Common options
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

struct CommandDescriptor
{
    CommandId id;
    std::string_view cli_name;
    std::string_view description;
    CommonOptionProfile common_option_profile;
};

constexpr auto GetCommandCatalog()
{
    return std::array{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION, PROFILE)                           \
        CommandDescriptor{                                                                     \
            CommandId::COMMAND_ID,                                                             \
            CLI_NAME,                                                                          \
            DESCRIPTION,                                                                       \
            CommonOptionProfile::PROFILE                                                       \
        },
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
    };
}

// Validation contract
enum class ValidationPhase : std::uint8_t
{
    Parse = 0,
    Prepare = 1
};

struct ValidationIssue
{
    std::string option_name;
    ValidationPhase phase;
    LogLevel level;
    std::string message;
    bool auto_corrected{ false };
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
