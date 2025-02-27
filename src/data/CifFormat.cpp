#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "StringHelper.hpp"

#include <fstream>
#include <iostream>
#include <cstring>
#include <sstream>
#include <stdexcept>

CifFormat::CifFormat(void)
{

}

CifFormat::~CifFormat()
{

}

void CifFormat::LoadHeader(const std::string & filename)
{
    LoadPdbxData(filename);
}

void CifFormat::PrintHeader(void) const
{
    
}

void CifFormat::LoadDataArray(const std::string & filename)
{
    LoadAtomSiteData(filename);
}

void CifFormat::LoadPdbxData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadPdbxDatabase failed!");
    }
    std::string line, header;
    auto found_model_id{ false };
    auto found_map_id{ false };
    while (std::getline(infile, line))
    {
        if (line.find("_pdbx_database_related.db_id") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header >> m_map_id;
            found_map_id = true;
        }
        if (line.find("_pdbx_database_status.entry_id") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header >> m_model_id;
            found_model_id = true;
        }
        if (found_model_id && found_map_id) break;
    }
}

void CifFormat::LoadAtomSiteData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadAtomSiteData failed!");
    }

    std::string line;
    auto in_atom_site_loop{ false };
    while (std::getline(infile, line))
    {
        if (line.find("loop_") != std::string::npos) in_atom_site_loop = false;
        if (line.find("_atom_site.") != std::string::npos)
        {
            in_atom_site_loop = true;
            continue;
        }
        auto is_special{ false };
        if (in_atom_site_loop && !line.empty() && line[0] != '#')
        {
            std::istringstream iss(line);
            auto atom{ std::make_shared<AtomSite>() };
            iss >> atom->group_PDB >> atom->id >> atom->type_symbol
                >> atom->label_atom_id >> atom->label_alt_id
                >> atom->label_comp_id >> atom->label_asym_id
                >> atom->label_entity_id >> atom->label_seq_id
                >> atom->pdbx_PDB_ins_code
                >> atom->position_x >> atom->position_y >> atom->position_z
                >> atom->occupancy >> atom->B_iso_or_equiv >> atom->pdbx_formal_charge
                >> atom->auth_seq_id >> atom->auth_comp_id
                >> atom->auth_asym_id >> atom->auth_atom_id >> atom->pdbx_PDB_model_num;
            if (atom->group_PDB == "HETATM") is_special = true;
            if (atom->pdbx_PDB_model_num > 1) break;  // end of reading file for multi-model PDB
            if (atom->label_alt_id != "." && atom->label_alt_id != "A") continue; // skip if there are other alternate position presented
            BuildAtomObject(atom, is_special);
        }
    }
}

void CifFormat::BuildAtomObject(std::any atom_info, bool is_special_atom)
{
    try
    {
        auto atom{ std::any_cast<std::shared_ptr<AtomSite>>(atom_info) };
        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetElement(atom->type_symbol);
        atom_object->SetRemoteness(StringHelper::ExtractCharAsString(atom->label_atom_id, atom->type_symbol.size()));
        atom_object->SetBranch(StringHelper::ExtractCharAsString(atom->label_atom_id, atom->type_symbol.size()+1));
        atom_object->SetResidue(atom->label_comp_id);
        atom_object->SetIndicator(atom->label_alt_id);
        auto residue_id{ (atom->label_seq_id == ".") ? -1 : std::stoi(atom->label_seq_id) };
        atom_object->SetResidueID(residue_id);
        atom_object->SetSerialID(atom->id);
        atom_object->SetChainID(atom->label_asym_id);
        atom_object->SetPosition(atom->position_x, atom->position_y, atom->position_z);
        atom_object->SetOccupancy(atom->occupancy);
        atom_object->SetTemperature(atom->B_iso_or_equiv);
        atom_object->SetSpecialAtomFlag(is_special_atom);
        //atom_object->SetCharge();
        m_atom_object_list.emplace_back(std::move(atom_object));
    }
    catch (const std::bad_any_cast &)
    {
        
    }
}

std::vector<std::unique_ptr<AtomObject>> CifFormat::GetAtomObjectList(void)
{
    return std::move(m_atom_object_list);
}

std::string CifFormat::GetPdbID(void) const
{
    return m_model_id;
}

std::string CifFormat::GetEmdID(void) const
{
    return m_map_id;
}