#pragma once

#include <filesystem>
#include "FilePathHelper.hpp"

struct CommandOptions
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::filesystem::path database_path{"database.sqlite"};
    std::filesystem::path folder_path{""};

    void SetFolderPath(const std::filesystem::path & path)
    {
        folder_path = std::filesystem::path(FilePathHelper::EnsureTrailingSlash(path));
    }
};
