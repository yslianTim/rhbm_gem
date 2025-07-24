#include "CommandBase.hpp"
#include "DataObjectManager.hpp"

#include <CLI/CLI.hpp>

CommandBase::~CommandBase() = default;

void CommandBase::SetDataManager(std::unique_ptr<DataObjectManager> manager)
{
    m_data_manager = std::move(manager);
}

DataObjectManager * CommandBase::GetDataManagerPtr(void)
{
    if (!m_data_manager)
    {
        m_data_manager = std::make_unique<DataObjectManager>();
    }
    return m_data_manager.get();
}

const DataObjectManager * CommandBase::GetDataManagerPtr(void) const
{
    if (!m_data_manager)
    {
        m_data_manager = std::make_unique<DataObjectManager>();
    }
    return m_data_manager.get();
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
