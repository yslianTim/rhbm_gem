#pragma once

#include <iostream>
#include <memory>
#include <string>
#include "CommandBase.hpp"

class PotentialDisplayCommand : public CommandBase
{
    std::string m_database_path, m_model_key_tag, m_folder_path;

public:
    PotentialDisplayCommand(void);
    ~PotentialDisplayCommand() = default;
    void Execute(void) override;

    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetModelKeyTag(const std::string & tag) { m_model_key_tag = tag; }
    void SetFolderPath(const std::string & path) { m_folder_path = path; }

};