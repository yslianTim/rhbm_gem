#include "CommandBase.hpp"

#include <CLI/CLI.hpp>

void CommandBase::RegisterCommandOptions(CLI::App * command, CommandOptions & options)
{
    command->add_option("-d,--database", options.database_path,
        "Database file path")->default_val("database.sqlite");
    command->add_option("-o,--folder", options.folder_path,
        "folder path for output files")->default_val("");
    command->add_option("-j,--jobs", options.thread_size,
        "Number of threads")->default_val(1);
    command->add_option("-v,--verbose", options.verbose_level,
        "Verbose level")->default_val(3);
}