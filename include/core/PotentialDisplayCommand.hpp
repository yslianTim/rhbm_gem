#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include "CommandBase.hpp"

class DataObjectManager;

class PotentialDisplayCommand : public CommandBase
{
    std::string m_database_path, m_folder_path;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;

public:
    PotentialDisplayCommand(void);
    ~PotentialDisplayCommand();
    void Execute(void) override;

    void SetDatabasePath(const std::string & path) { m_database_path = path; }
    void SetFolderPath(const std::string & path) { m_folder_path = path; }
    void SetModelKeyTagList(const std::string & value);
    void SetRefModelKeyTagList(const std::string & map_key, const std::string & value);
    void LoadModelObjects(DataObjectManager * data_manager);
    void LoadRefModelObjects(DataObjectManager * data_manager);

};