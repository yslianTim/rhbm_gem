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
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
    COMMAND_ID,
#include <rhbm_gem/core/command/CommandList.def>
#undef RHBM_GEM_COMMAND
};

struct CommandDescriptor
{
    CommandId id;
    std::string_view cli_name;
    std::string_view description;
};

constexpr auto GetCommandCatalog()
{
    return std::array{
#define RHBM_GEM_COMMAND(COMMAND_ID, CLI_NAME, DESCRIPTION)                                    \
        CommandDescriptor{                                                                     \
            CommandId::COMMAND_ID,                                                             \
            CLI_NAME,                                                                          \
            DESCRIPTION                                                                        \
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

} // namespace rhbm_gem
