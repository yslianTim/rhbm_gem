#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "StringHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
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
    LoadElementTypeList(filename);
}

void CifFormat::PrintHeader(void) const
{
    
}

void CifFormat::LoadDataArray(const std::string & filename)
{
    LoadStructConfData(filename);
    LoadAtomSiteData(filename);
}

void CifFormat::LoadPdbxData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadPdbxData failed!");
    }
    std::string line, header;
    auto found_model_id{ false };
    auto found_map_id{ false };
    auto found_resolution{ false };
    auto found_resolution_method{ false };
    while (std::getline(infile, line))
    {
        if (line.find("_pdbx_database_related.db_id ") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header >> m_map_id;
            found_map_id = true;
        }
        if (line.find("_pdbx_database_status.entry_id ") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header >> m_model_id;
            found_model_id = true;
        }
        if (line.find("_em_3d_reconstruction.resolution ") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header >> m_resolution;
            found_resolution = true;
        }
        if (line.find("_em_3d_reconstruction.resolution_method ") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header;
            iss >> std::quoted(m_resolution_method, '\'');
            found_resolution_method = true;
        }
        if (found_model_id && found_map_id && found_resolution && found_resolution_method) break;
    }
}

void CifFormat::LoadElementTypeList(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadElementTypeList failed!");
    }

    ParseLoopBlock(infile, "_atom_type.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto element_type_string{ token_list[index_map.at("symbol")] };
            auto element{ AtomicInfoHelper::GetElementFromString(element_type_string) };
            m_element_type_list.emplace_back(element);
        }
    );
    std::cout <<"Number of element types in PDBx/mmCIF = "<< m_element_type_list.size() << std::endl;
}

void CifFormat::LoadAtomSiteData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadAtomSiteData failed!");
    }

    ParseLoopBlock(infile, "_atom_site.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto atom{ std::make_unique<AtomSite>() };
            atom->group_PDB          = token_list[index_map.at("group_PDB")];
            atom->id                 = std::stoi(token_list[index_map.at("id")]);
            atom->type_symbol        = token_list[index_map.at("type_symbol")];
            atom->label_atom_id      = token_list[index_map.at("label_atom_id")];
            atom->label_alt_id       = token_list[index_map.at("label_alt_id")];
            atom->label_comp_id      = token_list[index_map.at("label_comp_id")];
            atom->label_asym_id      = token_list[index_map.at("label_asym_id")];
            atom->label_entity_id    = token_list[index_map.at("label_entity_id")];
            atom->label_seq_id       = token_list[index_map.at("label_seq_id")];
            atom->pdbx_PDB_ins_code  = token_list[index_map.at("pdbx_PDB_ins_code")];
            atom->position_x         = std::stof(token_list[index_map.at("Cartn_x")]);
            atom->position_y         = std::stof(token_list[index_map.at("Cartn_y")]);
            atom->position_z         = std::stof(token_list[index_map.at("Cartn_z")]);
            atom->occupancy          = std::stof(token_list[index_map.at("occupancy")]);
            atom->B_iso_or_equiv     = std::stof(token_list[index_map.at("B_iso_or_equiv")]);
            atom->pdbx_formal_charge = token_list[index_map.at("pdbx_formal_charge")];
            atom->auth_seq_id        = token_list[index_map.at("auth_seq_id")];
            atom->auth_comp_id       = token_list[index_map.at("auth_comp_id")];
            atom->auth_asym_id       = token_list[index_map.at("auth_asym_id")];
            atom->auth_atom_id       = token_list[index_map.at("auth_atom_id")];
            atom->pdbx_PDB_model_num = std::stoi(token_list[index_map.at("pdbx_PDB_model_num")]);
            bool is_special{ (atom->group_PDB == "HETATM") ? true : false };
            if (atom->pdbx_PDB_model_num > 1) return; // end of reading file for multi-model PDB
            if (atom->label_alt_id != "." && atom->label_alt_id != "A") return; // skip if there are other alternate position presented
            BuildAtomObject(atom.get(), is_special);
        }
    );
}

void CifFormat::LoadStructConfData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadStructConfData failed!");
    }

    ParseLoopBlock(infile, "_struct_conf.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto entry{ std::make_unique<StructConf>() };
            entry->conf_type_id          = token_list[index_map.at("conf_type_id")];
            entry->id                    = token_list[index_map.at("id")];
            entry->beg_label_comp_id     = token_list[index_map.at("beg_label_comp_id")];
            entry->beg_label_asym_id     = token_list[index_map.at("beg_label_asym_id")];
            entry->beg_label_seq_id      = token_list[index_map.at("beg_label_seq_id")];
            entry->end_label_comp_id     = token_list[index_map.at("end_label_comp_id")];
            entry->end_label_asym_id     = token_list[index_map.at("end_label_asym_id")];
            entry->end_label_seq_id      = token_list[index_map.at("end_label_seq_id")];
            entry->pdbx_PDB_helix_length = std::stoi(token_list[index_map.at("pdbx_PDB_helix_length")]);
            m_struct_conf_list.emplace_back(std::move(entry));
        }
    );
}

void CifFormat::BuildAtomObject(std::any atom_info, bool is_special_atom)
{
    try
    {
        auto atom{ std::any_cast<AtomSite *>(atom_info) };
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
        SetStructureInfo(atom_object.get());
        m_atom_object_list.emplace_back(std::move(atom_object));
    }
    catch (const std::bad_any_cast &)
    {
        std::cout <<"Error: bad any cast in BuildAtomObject"<< std::endl;
    }
}

void CifFormat::SetStructureInfo(AtomObject * atom_object)
{
    auto chain_id{ atom_object->GetChainID() };
    auto residue_id{ atom_object->GetResidueID() };
    for (const auto & entry : m_struct_conf_list)
    {
        if (chain_id == entry->beg_label_asym_id || chain_id == entry->end_label_asym_id)
        {
            auto beg{ std::stoi(entry->beg_label_seq_id) };
            auto end{ std::stoi(entry->end_label_seq_id) };
            if (residue_id >= beg && residue_id <= end)
            {
                auto structure = AtomicInfoHelper::GetStructureFromString(entry->conf_type_id);
                atom_object->SetStructure(structure);
                return;
            }
        }
    }
    atom_object->SetStructure(Structure::FREE);
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

double CifFormat::GetResolution(void) const
{
    return std::stod(m_resolution);
}

std::string CifFormat::GetResolutionMethod(void) const
{
    return m_resolution_method;
}

void CifFormat::ParseLoopBlock(
    std::ifstream & infile,
    const std::string & block_prefix,
    const std::function<void(const std::unordered_map<std::string, size_t> &,
                             const std::vector<std::string> &)> & row_handler)
{
    std::string line;
    bool in_loop{ false };
    bool header_parsed{ false };
    std::vector<std::string> field_order;
    std::unordered_map<std::string, size_t> index_map;
    while (std::getline(infile, line))
    {
        if (in_loop == false)
        {
            if (line.find("loop_") != std::string::npos)
            {
                in_loop = true;
            }
            continue;
        }
        if (header_parsed == false)
        {
            if (line.rfind(block_prefix, 0) == 0)
            {
                std::istringstream iss(line);
                std::string full;
                iss >> full;
                field_order.emplace_back(full.substr(block_prefix.size()));
                continue;
            }
            if (field_order.empty() == false)
            {
                header_parsed = true;
                for (size_t i = 0; i < field_order.size(); i++)
                {
                    index_map[field_order[i]] = i;
                }
            }
        }
        if (header_parsed == true)
        {
            if (line.empty() || line[0] == '#') break;
            std::istringstream iss(line);
            std::vector<std::string> token_list;
            std::string token;
            while (iss >> token) token_list.emplace_back(token);
            row_handler(index_map, token_list);
        }
    }
}