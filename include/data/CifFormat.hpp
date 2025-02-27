#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "ModelFileFormatBase.hpp"

class CifFormat : public ModelFileFormatBase
{
    struct AtomSite
    {
        std::string group_PDB;
        int id;
        std::string type_symbol;
        std::string label_atom_id;
        std::string label_alt_id;
        std::string label_comp_id;
        std::string label_asym_id;
        std::string label_entity_id;
        std::string label_seq_id;
        std::string pdbx_PDB_ins_code;
        double position_x;
        double position_y;
        double position_z;
        double occupancy;
        double B_iso_or_equiv;
        std::string pdbx_formal_charge;
        std::string auth_seq_id;
        std::string auth_comp_id;
        std::string auth_asym_id;
        std::string auth_atom_id;
        int pdbx_PDB_model_num;
    };

    std::string m_map_id, m_model_id;
    std::vector<std::unique_ptr<AtomObject>> m_atom_object_list;

public:
    CifFormat(void);
    ~CifFormat();
    void LoadHeader(const std::string & filename) override;
    void PrintHeader(void) const override;
    void LoadDataArray(const std::string & filename) override;
    void BuildAtomObject(std::any atom_info, bool is_special_atom) override;
    std::vector<std::unique_ptr<AtomObject>> GetAtomObjectList(void) override;
    std::string GetPdbID(void) const override;
    std::string GetEmdID(void) const override;

private:
    void LoadPdbxData(const std::string & filename);
    void LoadAtomSiteData(const std::string & filename);

};