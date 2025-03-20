#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "DataObjectVisitorBase.hpp"

class PotentialDisplayVisitor : public DataObjectVisitorBase
{
    std::string m_folder_path;
    std::vector<std::string> m_model_key_tag_list;

public:
    PotentialDisplayVisitor(void);
    ~PotentialDisplayVisitor();

    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void RunModelPainter(DataObjectManager * data_manager);
    void SetModelObjectKeyTagList(const std::vector<std::string> & value) { m_model_key_tag_list = std::move(value); }
    void SetFolderPath(const std::string & value) { m_folder_path = value; }

};