#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorAdapter.hpp"
#include "ResultDumpCommand.hpp"

class DataObjectBase;

class ResultDumpVisitor : public DataObjectVisitorAdapter
{
    int m_printer_choice;
    std::string m_folder_path;
    std::string m_map_key_tag;
    std::vector<std::string> m_model_key_tag_list;
    std::unordered_map<std::string, std::vector<AtomObject *>> m_selected_atom_list_map;

public:
    explicit ResultDumpVisitor(const ResultDumpCommand::Options & options);
    ~ResultDumpVisitor() = default;
    void VisitDataObjectManager(DataObjectManager * data_manager) override;

private:
    void BuildSelectedAtomList(DataObjectManager * data_manager);
    void RunAtomPositionDumping(DataObjectManager * data_manager);
    void RunMapValueDumping(DataObjectManager * data_manager);
    void RunGausEstimatesDumping(DataObjectManager * data_manager);
    void RunGroupGausEstimatesDumping(DataObjectManager * data_manager);

};
