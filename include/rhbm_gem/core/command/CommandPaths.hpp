#pragma once

#include <cstdlib>
#include <filesystem>

namespace rhbm_gem {

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

} // namespace rhbm_gem
