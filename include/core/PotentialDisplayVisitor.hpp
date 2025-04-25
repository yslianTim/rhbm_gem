#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"

class DataObjectBase;

class PotentialDisplayVisitor : public DataObjectVisitorBase
{
    std::string m_folder_path;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    std::vector<DataObjectBase *> m_model_object_list;
    std::vector<DataObjectBase *> m_sim_no_charge_model_object_list;
    std::vector<DataObjectBase *> m_sim_with_charge_model_object_list;

public:
    PotentialDisplayVisitor(void);
    ~PotentialDisplayVisitor();

    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void RunAtomPainter(ModelObject * model_object);
    void RunModelPainter(DataObjectManager * data_manager);
    void RunComparisonPainter(DataObjectManager * data_manager);
    void BuildModelObjectList(DataObjectManager * data_manager, std::vector<DataObjectBase *> & model_object_list);
    void BuildReferenceModelObjectList(DataObjectManager * data_manager, const std::string & type_key, std::vector<DataObjectBase *> & model_object_list);
    void SetModelObjectKeyTagList(std::vector<std::string> value) { m_model_key_tag_list = std::move(value); }
    void SetRefModelObjectKeyTagListMap(std::unordered_map<std::string, std::vector<std::string>> value) { m_ref_model_key_tag_list_map = std::move(value); }
    void SetFolderPath(const std::string & value) { m_folder_path = value; }

};