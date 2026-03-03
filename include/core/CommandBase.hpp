#pragma once

#include <filesystem>

#include "DataObjectManager.hpp"

namespace CLI
{
    class App;
}

namespace rhbm_gem {

struct CommandOptions
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::filesystem::path database_path{"database.sqlite"};
    std::filesystem::path folder_path{""};
};

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    virtual bool Execute() = 0;
    virtual void RegisterCLIOptionsExtend(::CLI::App * command) = 0;
    virtual const CommandOptions & GetOptions() const = 0;
    virtual CommandOptions & GetOptions() = 0;

    void RegisterCLIOptions(::CLI::App * command);
    void SetThreadSize(int value);
    void SetVerboseLevel(int value);
    void SetDatabasePath(const std::filesystem::path & path);
    void SetFolderPath(const std::filesystem::path & path);
    bool IsValidateOptions() const { return m_valiate_options; }
    DataObjectManager * GetDataManagerPtr() { return &m_data_manager; }
    const DataObjectManager * GetDataManagerPtr() const { return &m_data_manager; }

protected:
    DataObjectManager m_data_manager;
    bool m_valiate_options{ true };

    CommandBase() = default;
    void RegisterCLIOptionsBasic(::CLI::App * command);
};

} // namespace rhbm_gem
