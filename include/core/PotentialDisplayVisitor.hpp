#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "PotentialDisplayCommand.hpp"

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
    PotentialDisplayVisitor(AtomSelector * atom_selector, const PotentialDisplayCommand::Options & options);
    ~PotentialDisplayVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitDataObjectManager(DataObjectManager * data_manager) override;

    void SetRefModelObjectKeyTagListMap(const std::unordered_map<std::string, std::vector<std::string>> & value) { m_ref_model_key_tag_list_map = value; }

private:
    void AddAtomObjectToPainter(DataObjectManager * data_manager, PainterBase * painter);
    void AddModelObjectToPainter(DataObjectManager * data_manager, PainterBase * painter);
    void AddReferenceModelObjectToPainter(DataObjectManager * data_manager, PainterBase * painter);

};
