#pragma once

#include <filesystem>

struct CommandOptions
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::filesystem::path database_path{"database.sqlite"};
    std::filesystem::path folder_path{""};
};
