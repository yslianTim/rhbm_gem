#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "ResultDumpCommand.hpp"

class DataObjectBase;
class ModelObject;
class MapObject;

class ResultDumpVisitor : public DataObjectVisitorBase
{
    int m_printer_choice;
    std::string m_folder_path;
    std::string m_map_key_tag;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;
    std::vector<ModelObject *> m_model_object_list;
    MapObject * m_map_object;

public:
    explicit ResultDumpVisitor(const ResultDumpCommand::Options & options);
    ~ResultDumpVisitor();
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;

private:
    void BuildSelectedAtomList(void);
    void RunAtomPositionDumping(void);
    void RunMapValueDumping(void);
    void RunGausEstimatesDumping(void);
    void RunGroupGausEstimatesDumping(void);

};
