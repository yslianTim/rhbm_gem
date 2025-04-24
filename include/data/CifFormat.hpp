#pragma once

#include <string>
#include <vector>
#include <memory>
#include <functional>
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
        float position_x;
        float position_y;
        float position_z;
        float occupancy;
        float B_iso_or_equiv;
        std::string pdbx_formal_charge;
        std::string auth_seq_id;
        std::string auth_comp_id;
        std::string auth_asym_id;
        std::string auth_atom_id;
        int pdbx_PDB_model_num;
    };

    struct StructConf
    {
        std::string conf_type_id;
        std::string id;
        std::string pdbx_PDB_helix_id;
        std::string beg_label_comp_id;
        std::string beg_label_asym_id;
        std::string beg_label_seq_id;
        std::string pdbx_beg_PDB_ins_code;
        std::string end_label_comp_id;
        std::string end_label_asym_id;
        std::string end_label_seq_id;
        std::string pdbx_end_PDB_ins_code;
        std::string beg_auth_comp_id;
        std::string beg_auth_asym_id;
        std::string beg_auth_seq_id;
        std::string end_auth_comp_id;
        std::string end_auth_asym_id;
        std::string end_auth_seq_id;
        std::string pdbx_PDB_helix_class;
        std::string details;
        int         pdbx_PDB_helix_length;
    };

    std::string m_map_id, m_model_id;
    std::string m_resolution, m_resolution_method;
    std::vector<std::unique_ptr<AtomObject>> m_atom_object_list;
    std::vector<std::unique_ptr<StructConf>> m_struct_conf_list;

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
    double GetResolution(void) const override;
    std::string GetResolutionMethod(void) const override;

private:
    void LoadPdbxData(const std::string & filename);
    void LoadAtomSiteData(const std::string & filename);
    void LoadStructConfData(const std::string & filename);
    void ParseLoopBlock(std::ifstream & infile,
        const std::string & prefix,
        const std::function<void(const std::unordered_map<std::string, size_t> &,
                                 const std::vector<std::string> &)> & row_handler);

};