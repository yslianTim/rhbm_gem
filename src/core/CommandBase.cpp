#include "CommandBase.hpp"

#include <CLI/CLI.hpp>

DataObjectManager * CommandBase::GetDataManagerPtr(void)
{
    return &m_data_manager;
}

const DataObjectManager * CommandBase::GetDataManagerPtr(void) const
{
    return &m_data_manager;
}

void CommandBase::RegisterCLIOptionsBasic(CLI::App * command)
{
    auto & options{ GetOptions() };
    command->add_option("-d,--database", options.database_path,
        "Database file path")->default_val(options.database_path.string());
    command->add_option("-o,--folder", options.folder_path,
        "folder path for output files")->default_val(options.folder_path.string());
    command->add_option("-j,--jobs", options.thread_size,
        "Number of threads")->default_val(options.thread_size);
    command->add_option("-v,--verbose", options.verbose_level,
        "Verbose level")->default_val(options.verbose_level);
}

void CommandBase::RegisterCLIOptions(CLI::App * command)
{
    RegisterCLIOptionsExtend(command);
    RegisterCLIOptionsBasic(command);
}
