#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "PotentialDisplayCommand.hpp"

class DataObjectBase;
class AtomSelector;
class PainterBase;
class ModelObject;

class PotentialDisplayVisitor : public DataObjectVisitorBase
{
    int m_painter_choice;
    std::string m_folder_path;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    AtomSelector * m_atom_selector;
    std::vector<ModelObject *> m_model_object_list;
    std::vector<ModelObject *> m_ordered_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ordered_ref_model_object_list_map;

public:
    PotentialDisplayVisitor(AtomSelector * atom_selector, const PotentialDisplayCommand::Options & options);
    ~PotentialDisplayVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;

    void Finalize(void);
    void SetRefModelObjectKeyTagListMap(const std::unordered_map<std::string, std::vector<std::string>> & value) { m_ref_model_key_tag_list_map = value; }

private:
    void BuildOrderedModelObjectList(void);
    void BuildOrderedRefModelObjectListMap(void);

};
