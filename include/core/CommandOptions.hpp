#pragma once

#include <string>

struct CommandOptions
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::string database_path{"database.sqlite"};
    std::string folder_path{""};
};
