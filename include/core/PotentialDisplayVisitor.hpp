#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"

class DataObjectBase;
class AtomSelector;
class PainterBase;

class PotentialDisplayVisitor : public DataObjectVisitorBase
{
    int m_painter_choice;
    std::string m_folder_path;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<std::string>> m_ref_model_key_tag_list_map;
    AtomSelector * m_atom_selector;

public:
    PotentialDisplayVisitor(AtomSelector * atom_selector);
    ~PotentialDisplayVisitor();

    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void VisitDataObjectManager(DataObjectManager * data_manager) override;

    void SetModelObjectKeyTagList(const std::vector<std::string> & value) { m_model_key_tag_list = value; }
    void SetRefModelObjectKeyTagListMap(const std::unordered_map<std::string, std::vector<std::string>> & value) { m_ref_model_key_tag_list_map = value; }
    void SetFolderPath(const std::string & value) { m_folder_path = value; }
    void SetPainterChoice(int value) { m_painter_choice = value; }

private:
    void AddAtomObjectToPainter(DataObjectManager * data_manager, PainterBase * painter);
    void AddModelObjectToPainter(DataObjectManager * data_manager, PainterBase * painter);
    void AddReferenceModelObjectToPainter(DataObjectManager * data_manager, PainterBase * painter);

};
