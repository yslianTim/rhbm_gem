#pragma once

#include <filesystem>

#include "DataObjectManager.hpp"

namespace CLI
{
    class App;
}

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
    virtual bool Execute(void) = 0;
    virtual bool ValidateOptions(void) const = 0;
    virtual void RegisterCLIOptionsExtend(CLI::App * command) = 0;
    virtual const CommandOptions & GetOptions(void) const = 0;
    virtual CommandOptions & GetOptions(void) = 0;

    void RegisterCLIOptions(CLI::App * command);
    void SetThreadSize(int value);
    void SetVerboseLevel(int value);
    void SetDatabasePath(const std::filesystem::path & path);
    void SetFolderPath(const std::filesystem::path & path);
    DataObjectManager * GetDataManagerPtr(void);
    const DataObjectManager * GetDataManagerPtr(void) const;

protected:
    DataObjectManager m_data_manager;

    CommandBase(void) = default;
    void RegisterCLIOptionsBasic(CLI::App * command);

};
