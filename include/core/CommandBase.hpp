#pragma once

#include <memory>

#include "CommandOptions.hpp"

namespace CLI
{
    class App;
}
class DataObjectManager;

class CommandBase
{
public:
    virtual ~CommandBase();
    virtual bool Execute(void) = 0;
    virtual bool ValidateOptions(void) const = 0;
    void RegisterCLIOptions(CLI::App * command);
    virtual void RegisterCLIOptionsExtend(CLI::App * command) = 0;
    virtual const CommandOptions & GetOptions(void) const = 0;
    virtual CommandOptions & GetOptions(void) = 0;

    void SetDataManager(std::unique_ptr<DataObjectManager> manager);
    DataObjectManager * GetDataManagerPtr(void);
    const DataObjectManager * GetDataManagerPtr(void) const;

protected:
    CommandBase(void);
    std::unique_ptr<DataObjectManager> m_data_manager;
    void RegisterCLIOptionsBasic(CLI::App * command);

};
