#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"

class DataObjectBase;

class ResultDumpVisitor : public DataObjectVisitorBase
{
    int m_printer_choice;
    std::string m_folder_path;
    std::string m_map_key_tag;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;

public:
    ResultDumpVisitor();
    ~ResultDumpVisitor();

    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetMapObjectKeyTag(const std::string & value) { m_map_key_tag = value; }
    void SetModelObjectKeyTagList(std::vector<std::string> value) { m_model_key_tag_list = std::move(value); }
    void SetFolderPath(const std::string & value);
    void SetPrinterChoice(int value) { m_printer_choice = value; }

private:
    void BuildSelectedAtomList(DataObjectManager * data_manager);
    void RunAtomPositionDumping(DataObjectManager * data_manager);
    void RunMapValueDumping(DataObjectManager * data_manager);
    void RunGausEstimatesDumping(DataObjectManager * data_manager);

};
