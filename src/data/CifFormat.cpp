#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "StringHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"
#include "AtomicModelDataBlock.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cctype>

CifFormat::CifFormat(void) :
    m_data_block{ std::make_unique<AtomicModelDataBlock>() }
{

}

CifFormat::~CifFormat()
{

}

void CifFormat::LoadHeader(const std::string & filename)
{
    LoadDatabaseInfo(filename);
    LoadEntityInfo(filename);
    LoadPdbxData(filename);
    LoadElementTypeList(filename);
    LoadStructHelixInfo(filename);
    LoadStructSheetInfo(filename);
}

void CifFormat::PrintHeader(void) const
{
    std::cout <<"#Entities = "<< m_data_block->GetEntityTypeMap().size() << std::endl;
    for (auto & [entity_id, chain_id] : m_data_block->GetChainIDListMap())
    {
        std::cout << entity_id <<": ";
        for (auto & chain : chain_id)
        {
            std::cout << chain << ",";
        }
        std::cout << std::endl;
    }
    for (auto element : m_data_block->GetElementTypeList())
    {
        std::cout << AtomicInfoHelper::GetLabel(element) << ",";
    }
    std::cout << std::endl;
}

void CifFormat::LoadDataArray(const std::string & filename)
{
    LoadAtomSiteData(filename);
}

void CifFormat::LoadDatabaseInfo(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadDatabaseInfo failed!");
    }

    ParseLoopBlock(infile, "_database_2.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            std::unordered_map<std::string, std::string> data_map;
            auto key{ token_list[index_map.at("database_id")] };
            data_map[key] = token_list[index_map.at("database_code")];
            if (data_map.find("PDB") != data_map.end())
            {
                m_data_block->SetPdbID(data_map.at("PDB"));
            }
            if (data_map.find("EMDB") != data_map.end())
            {
                m_data_block->SetEmdID(data_map.at("EMDB"));
            }
        }
    );
}

void CifFormat::LoadEntityInfo(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadDatabaseInfo failed!");
    }

    ParseLoopBlock(infile, "_entity.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto entity_id{ token_list[index_map.at("id")] };
            auto entity_type{ token_list[index_map.at("type")] };
            auto molecules_size_string{ token_list[index_map.at("pdbx_number_of_molecules")] };
            int molecules_size{ -1 };
            try
            {
                molecules_size = std::stoi(molecules_size_string);
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() <<": "<< molecules_size_string << '\n';
            }
            
            m_data_block->AddEntityTypeInEntityMap(
                entity_id, AtomicInfoHelper::GetEntityFromString(entity_type));
            m_data_block->AddMoleculesSizeInEntityMap(entity_id, molecules_size);
        }
    );

    ParseLoopBlock(infile, "_struct_asym.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto entity_id{ token_list[index_map.at("entity_id")] };
            auto chain_id{ token_list[index_map.at("id")] };
            m_data_block->AddChainIDInEntityMap(entity_id, chain_id);
        }
    );
}

void CifFormat::LoadPdbxData(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadPdbxData failed!");
    }
    std::string line, header, resolution, resolution_method;
    auto found_resolution{ false };
    auto found_resolution_method{ false };
    while (std::getline(infile, line))
    {
        if (line.find("_em_3d_reconstruction.resolution ") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header >> resolution;
            m_data_block->SetResolution(resolution);
            found_resolution = true;
        }
        if (line.find("_em_3d_reconstruction.resolution_method ") != std::string::npos)
        {
            std::istringstream iss(line);
            iss >> header;
            iss >> std::quoted(resolution_method, '\'');
            m_data_block->SetResolutionMethod(resolution_method);
            found_resolution_method = true;
        }
        if (found_resolution && found_resolution_method) break;
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
            m_data_block->AddElementType(element);
        }
    );
}

void CifFormat::LoadStructHelixInfo(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadStructHelixInfo failed!");
    }

    ParseLoopBlock(infile, "_struct_conf.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto helix_id{ token_list[index_map.at("id")] };
            auto conf_type{ token_list[index_map.at("conf_type_id")] };
            auto chain_id_beg{ token_list[index_map.at("beg_label_asym_id")] };
            auto reisude_id_beg{ token_list[index_map.at("beg_label_seq_id")] };
            auto chain_id_end{ token_list[index_map.at("end_label_asym_id")] };
            auto reisude_id_end{ token_list[index_map.at("end_label_seq_id")] };
            m_data_block->AddHelixRange(
                helix_id,
                {chain_id_beg, reisude_id_beg, chain_id_end, reisude_id_end, conf_type});
        }
    );
}

void CifFormat::LoadStructSheetInfo(const std::string & filename)
{
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        std::cerr << "Cannot open the file: " << filename << std::endl;
        throw std::runtime_error("LoadStructSheetInfo failed!");
    }

    ParseLoopBlock(infile, "_struct_sheet.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto sheet_id{ token_list[index_map.at("id")] };
            auto strands_size{ std::stoi(token_list[index_map.at("number_strands")]) };
            m_data_block->AddSheetStrands(sheet_id, strands_size);
        }
    );

    ParseLoopBlock(infile, "_struct_sheet_range.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto sheet_id{ token_list[index_map.at("sheet_id")] };
            auto range_id{ token_list[index_map.at("id")] };
            auto composite_sheet_id{ sheet_id + range_id };
            auto chain_id_beg{ token_list[index_map.at("beg_label_asym_id")] };
            auto reisude_id_beg{ token_list[index_map.at("beg_label_seq_id")] };
            auto chain_id_end{ token_list[index_map.at("end_label_asym_id")] };
            auto reisude_id_end{ token_list[index_map.at("end_label_seq_id")] };
            m_data_block->AddSheetRange(
                composite_sheet_id,
                {chain_id_beg, reisude_id_beg, chain_id_end, reisude_id_end});
        }
    );
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
            static AtomObject * last_atom_object{ nullptr };
            auto group_type{ token_list[index_map.at("group_PDB")] };
            auto element_type{ token_list[index_map.at("type_symbol")] };
            if (element_type == "H") return; // Skip hydrogen atom
            auto is_special_atom{ (group_type == "HETATM") ? true : false };
            auto atom_object{ std::make_unique<AtomObject>() };
            atom_object->SetElement(element_type);
            atom_object->SetRemoteness(
                StringHelper::ExtractCharAsString(
                    token_list[index_map.at("label_atom_id")],
                    token_list[index_map.at("type_symbol")].size()
                ));
            atom_object->SetBranch(
                StringHelper::ExtractCharAsString(
                    token_list[index_map.at("label_atom_id")],
                    token_list[index_map.at("type_symbol")].size() + 1
                ));
            atom_object->SetResidue(token_list[index_map.at("label_comp_id")]);
            auto indicator{ token_list[index_map.at("label_alt_id")] };
            atom_object->SetIndicator(indicator);
            auto residue_id{ token_list[index_map.at("label_seq_id")] };
            atom_object->SetResidueID((residue_id == ".") ? -1 : std::stoi(residue_id));
            atom_object->SetSerialID(std::stoi(token_list[index_map.at("id")]));
            atom_object->SetChainID(token_list[index_map.at("label_asym_id")]);
            
            auto position_x{ std::stof(token_list[index_map.at("Cartn_x")]) };
            auto position_y{ std::stof(token_list[index_map.at("Cartn_y")]) };
            auto position_z{ std::stof(token_list[index_map.at("Cartn_z")]) };
            auto occupancy{ std::stof(token_list[index_map.at("occupancy")]) };
            auto temperature{ std::stof(token_list[index_map.at("B_iso_or_equiv")]) };
            atom_object->SetPosition(position_x, position_y, position_z);
            atom_object->SetOccupancy(occupancy);
            atom_object->SetTemperature(temperature);
            atom_object->SetSpecialAtomFlag(is_special_atom);
            m_data_block->SetStructureInfo(atom_object.get());

            auto model_number{ std::stoi(token_list[index_map.at("pdbx_PDB_model_num")]) };
            if (indicator == ".")
            {
                last_atom_object = nullptr;
                m_data_block->AddAtomObject(model_number, std::move(atom_object));
            }
            else if (indicator == "A")
            {
                last_atom_object = atom_object.get();
                m_data_block->AddAtomObject(model_number, std::move(atom_object));
            }
            else
            {
                if (last_atom_object == nullptr)
                {
                    std::cout <<"CifFormat::LoadAtomSiteData() atom_object is missing."<< std::endl;
                    return;
                }
                last_atom_object->AddAlternatePosition(indicator, {position_x, position_y, position_z});
                last_atom_object->AddAlternateOccupancy(indicator, occupancy);
                last_atom_object->AddAlternateTemperature(indicator, temperature);
            }
        }
    );
}

AtomicModelDataBlock * CifFormat::GetDataBlockPtr(void)
{
    return m_data_block.get();
}

void CifFormat::ParseLoopBlock(
    std::ifstream & infile, std::string_view data_block_prefix,
    const std::function<void(const std::unordered_map<std::string, size_t> &,
                             const std::vector<std::string> &)> & table_handler)
{
    std::string line;
    auto header_parsed{ false };
    std::vector<std::string> data_column_list;
    std::unordered_map<std::string, size_t> column_index_map;
    while (std::getline(infile, line))
    {
        if (header_parsed == false)
        {
            if (line.rfind(data_block_prefix, 0) == 0)
            {
                auto pos{ line.find_first_of(" \t") }; // extract the token up to the first whitespace
                std::string full_line{ (pos == std::string::npos) ? line : line.substr(0, pos) };
                data_column_list.emplace_back(full_line.substr(data_block_prefix.size()));
                continue;
            }
            
            if (data_column_list.empty() == true) continue;
            column_index_map.reserve(data_column_list.size());
            for (size_t i = 0; i < data_column_list.size(); i++)
            {
                column_index_map[data_column_list.at(i)] = i;
            }
            header_parsed = true;
            
        }
        if (header_parsed == true)
        {
            if (line.empty() || line[0] == '#') break;

            const auto expected_column_size{ column_index_map.size() };
            std::vector<std::string> token_list;
            token_list.reserve(expected_column_size);

            // initial tokens from this line
            auto initial{ StringHelper::SpliteStringLineAsTokens(line, expected_column_size) };
            token_list.insert(token_list.end(), initial.begin(), initial.end());

            // now read continuation lines until we have all fields
            bool in_multiline{ false };
            std::string multiline_content;
            std::string next_line;
            while (token_list.size() < expected_column_size && std::getline(infile, next_line))
            {
                if (next_line.empty()) continue;
                if (!in_multiline && next_line[0] == ';')
                {
                    // start of multiline literal for the next field
                    in_multiline = true;
                    multiline_content.clear();
                    // strip leading semicolon
                    auto content{ next_line.substr(1) };
                    if (!content.empty()) multiline_content = content;
                }
                else if (in_multiline)
                {
                    if (!next_line.empty() && next_line[0] == ';')
                    {
                        // end of multiline literal
                        token_list.emplace_back(std::move(multiline_content));
                        in_multiline = false;
                    }
                    else
                    {
                        // accumulate lines
                        multiline_content += "\n" + next_line;
                    }
                }
                else
                {
                    // normal continuation tokens
                    auto more{ StringHelper::SpliteStringLineAsTokens(next_line, expected_column_size) };
                    token_list.insert(token_list.end(), more.begin(), more.end());
                }
            }

            table_handler(column_index_map, token_list);
        }
    }
}