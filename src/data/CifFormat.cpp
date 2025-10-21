#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "ModelObject.hpp"
#include "StringHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "ChemicalDataHelper.hpp"
#include "ComponentHelper.hpp"
#include "AtomicModelDataBlock.hpp"
#include "LocalPotentialEntry.hpp"
#include "ChemicalComponentEntry.hpp"
#include "ComponentKeySystem.hpp"
#include "KDTreeAlgorithm.hpp"
#include "Logger.hpp"

#include <iomanip>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <ostream>
#include <algorithm>

CifFormat::CifFormat(void) :
    m_data_block{ std::make_unique<AtomicModelDataBlock>() }
{
    Logger::Log(LogLevel::Debug, "CifFormat::CifFormat() called");
}

CifFormat::~CifFormat()
{
    Logger::Log(LogLevel::Debug, "CifFormat::~CifFormat() called");
}

void CifFormat::LoadHeader(const std::string & filename)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadHeader() called");
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        Logger::Log(LogLevel::Error, "Cannot open the file: " + filename);
        throw std::runtime_error("LoadHeader failed!");
    }

    LoadChemicalComponentBlock(infile);
    LoadDatabaseBlock(infile);
    LoadEntityBlock(infile);
    LoadPdbxData(infile);
    LoadAtomTypeBlock(infile);
    LoadStructureConformationBlock(infile);
    LoadStructureSheetBlock(infile);

    if (m_find_component_bond_entry == false)
    {
        Logger::Log(LogLevel::Info,
            "No component bond entry found in the CIF file. Build default component bond entries.");
        BuildDefaultComponentBondEntry();
    }
}

void CifFormat::PrintHeader(void) const
{
    Logger::Log(LogLevel::Debug, "CifFormat::PrintHeader() called");
    std::ostringstream oss;
    oss << "CIF Header Information:\n";
    oss <<"#Entities = "<< m_data_block->GetEntityTypeMap().size() << "\n";
    for (auto & [entity_id, chain_id] : m_data_block->GetChainIDListMap())
    {
        oss <<"[" << entity_id <<"] : ";
        for (size_t i = 0; i < chain_id.size(); i++)
        {
            oss << chain_id.at(i);
            if (i < chain_id.size() - 1) oss << ",";
        }
        oss << "\n";
    }

    auto element_size{ m_data_block->GetElementTypeList().size() };
    oss <<"#Elementry types = "<< element_size << "\n";
    oss <<"Element type list : ";
    size_t count{ 0 };
    for (auto element : m_data_block->GetElementTypeList())
    {
        oss << ChemicalDataHelper::GetLabel(element);
        if (count < element_size - 1) oss << ",";
        count++;
    }
    oss << "\n";
    Logger::Log(LogLevel::Info, oss.str());
}

void CifFormat::LoadDataArray(const std::string & filename)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadDataArray() called");
    std::ifstream infile{ filename, std::ios::binary };
    if (!infile)
    {
        Logger::Log(LogLevel::Error, "Cannot open the file: " + filename);
        throw std::runtime_error("LoadDataArray failed!");
    }
    LoadAtomSiteBlock(infile);
    ConstructBondList();
    LoadStructureConnectionBlock(infile);
}

void CifFormat::LoadChemicalComponentBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadChemicalComponentBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_chem_comp.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto comp_id{ token_list[index_map.at("id")] };
            auto name{ token_list[index_map.at("name")] };
            auto type{ token_list[index_map.at("type")] };
            auto formula{ token_list[index_map.at("formula")] };
            auto formula_weight_str{ token_list[index_map.at("formula_weight")] };
            auto standard_flag_str{ token_list[index_map.at("mon_nstd_flag")] };
            auto formula_weight{ 0.0f };
            try
            {
                formula_weight = std::stof(formula_weight_str);
            }
            catch (const std::exception & e)
            {
                formula_weight = 0.0f;
            }
            auto entry{ std::make_unique<ChemicalComponentEntry>() };
            entry->SetComponentId(comp_id);
            entry->SetComponentName(name);
            entry->SetComponentType(type);
            entry->SetComponentFormula(formula);
            entry->SetComponentMolecularWeight(formula_weight);
            auto standard_flag{ false };
            if      (standard_flag_str == "n" || standard_flag_str == "no") standard_flag = false;
            else if (standard_flag_str == "y" || standard_flag_str == "yes") standard_flag = true;
            else standard_flag = false;
            entry->SetStandardMonomerFlag(standard_flag);

            m_data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
            auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
            m_data_block->AddChemicalComponentEntry(component_key, std::move(entry));

            m_find_chemical_component_entry = true;
        }
    );
    LoadChemicalComponentAtomBlock(infile);
    LoadChemicalComponentBondBlock(infile);
}

void CifFormat::LoadChemicalComponentAtomBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadChemicalComponentAtomBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_chem_comp_atom.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto comp_id{ token_list[index_map.at("comp_id")] };
            auto atom_id{ token_list[index_map.at("atom_id")] };
            auto element_symbol{ token_list[index_map.at("type_symbol")] };
            auto pdbx_aromatic_flag{ token_list[index_map.at("pdbx_aromatic_flag")] };
            auto pdbx_stereo_config{ token_list[index_map.at("pdbx_stereo_config")] };

            StringHelper::EraseCharFromString(atom_id, '\"');

            ComponentAtomEntry atom_entry;
            atom_entry.atom_id = atom_id;
            atom_entry.element_type = ChemicalDataHelper::GetElementFromString(element_symbol);
            atom_entry.aromatic_atom_flag = (pdbx_aromatic_flag == "Y") ? true : false;
            atom_entry.stereo_config = ChemicalDataHelper::GetStereoChemistryFromString(pdbx_stereo_config);

            m_data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
            auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
            auto atom_key{ m_data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id) };
            m_data_block->AddComponentAtomEntry(component_key, atom_key, atom_entry);

            m_find_component_atom_entry = true;
        }
    );
}

void CifFormat::LoadChemicalComponentBondBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadChemicalComponentBondBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_chem_comp_bond.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto comp_id{ token_list[index_map.at("comp_id")] };
            auto atom_id_1{ token_list[index_map.at("atom_id_1")] };
            auto atom_id_2{ token_list[index_map.at("atom_id_2")] };
            auto bond_order{ token_list[index_map.at("value_order")] };
            auto pdbx_aromatic_flag{ token_list[index_map.at("pdbx_aromatic_flag")] };
            auto pdbx_stereo_config{ token_list[index_map.at("pdbx_stereo_config")] };

            StringHelper::EraseCharFromString(atom_id_1, '\"');
            StringHelper::EraseCharFromString(atom_id_2, '\"');

            m_data_block->GetBondKeySystemPtr()->RegisterBond(atom_id_1, atom_id_2);
            auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
            auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(atom_id_1, atom_id_2) };
            auto bond_id{ m_data_block->GetBondKeySystemPtr()->GetBondId(bond_key) };

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = BondType::COVALENT; // The bonds in this block are all covalent bonds
            bond_entry.bond_order = ChemicalDataHelper::GetBondOrderFromString(bond_order);
            bond_entry.aromatic_atom_flag = (pdbx_aromatic_flag == "Y");
            bond_entry.stereo_config = ChemicalDataHelper::GetStereoChemistryFromString(pdbx_stereo_config);

            m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);

            m_find_component_bond_entry = true;
        }
    );
}

void CifFormat::LoadDatabaseBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadDatabaseBlock() called");
    infile.clear();
    infile.seekg(0);
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
            else
            {
                m_data_block->SetEmdID("X-RAY DIFF");
            }
        }
    );
}

void CifFormat::LoadEntityBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadEntityBlock() called");
    infile.clear();
    infile.seekg(0);
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
                Logger::Log(LogLevel::Error, "Invalid molecules size: " + molecules_size_string);
            }
            
            m_data_block->AddEntityTypeInEntityMap(
                entity_id, ChemicalDataHelper::GetEntityFromString(entity_type));
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

void CifFormat::LoadPdbxData(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadPdbxData() called");
    infile.clear();
    infile.seekg(0);
    std::string line, header, resolution, resolution_method;
    auto found_resolution{ false };
    auto found_resolution_method{ false };
    while (std::getline(infile, line))
    {
        StringHelper::StripCarriageReturn(line);
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
    if (!found_resolution || !found_resolution_method)
    {
        LoadXRayResolutionInfo(infile);
    }
}

void CifFormat::LoadXRayResolutionInfo(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadXRayResolutionInfo() called");
    infile.clear();
    infile.seekg(0);
    std::string line, header, resolution, resolution_method;
    auto found_resolution{ false };
    auto found_resolution_method{ false };
    while (std::getline(infile, line))
    {
        StringHelper::StripCarriageReturn(line);
        if (line.find("_refine.ls_d_res_high ") != std::string::npos && !found_resolution)
        {
            std::istringstream iss(line);
            iss >> header >> resolution;
            m_data_block->SetResolution(resolution);
            found_resolution = true;
        }

        if (line.find("_refine.pdbx_refine_id ") != std::string::npos && !found_resolution_method)
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

void CifFormat::LoadAtomTypeBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadAtomTypeBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_atom_type.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto element_type_string{ token_list[index_map.at("symbol")] };
            auto element{ ChemicalDataHelper::GetElementFromString(element_type_string) };
            m_data_block->AddElementType(element);
        }
    );
}

void CifFormat::LoadStructureConformationBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadStructureConformationBlock() called");
    infile.clear();
    infile.seekg(0);
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

void CifFormat::LoadStructureConnectionBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadStructureConnectionBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_struct_conn.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            auto conn_id{ token_list[index_map.at("id")] };
            auto conn_type_id{ token_list[index_map.at("conn_type_id")] };
            auto ptnr1_asym_id{ token_list[index_map.at("ptnr1_label_asym_id")] };
            auto ptnr1_comp_id{ token_list[index_map.at("ptnr1_label_comp_id")] };
            auto ptnr1_seq_id{ token_list[index_map.at("ptnr1_label_seq_id")] };
            auto ptnr1_atom_id{ token_list[index_map.at("ptnr1_label_atom_id")] };
            auto ptnr2_asym_id{ token_list[index_map.at("ptnr2_label_asym_id")] };
            auto ptnr2_comp_id{ token_list[index_map.at("ptnr2_label_comp_id")] };
            auto ptnr2_seq_id{ token_list[index_map.at("ptnr2_label_seq_id")] };
            auto ptnr2_atom_id{ token_list[index_map.at("ptnr2_label_atom_id")] };
            auto value_order_str{ token_list[index_map.at("pdbx_value_order")] };

            StringHelper::EraseCharFromString(ptnr1_atom_id, '\"');
            StringHelper::EraseCharFromString(ptnr2_atom_id, '\"');

            m_data_block->GetBondKeySystemPtr()->RegisterBond(ptnr1_atom_id, ptnr2_atom_id);
            auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(ptnr1_comp_id) };
            auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(ptnr1_atom_id, ptnr2_atom_id) };
            auto bond_id{ m_data_block->GetBondKeySystemPtr()->GetBondId(bond_key) };
            auto bond_type{ ChemicalDataHelper::GetBondTypeFromString(conn_type_id) };
            if (bond_type == BondType::HYDROGEN) return; // Skip hydrogen bond

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = bond_type;
            bond_entry.bond_order = (value_order_str == "?") ? BondOrder::UNK : ChemicalDataHelper::GetBondOrderFromString(value_order_str);
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;
            m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);

            auto model_number{ 1 };
            auto atom_1{ m_data_block->GetAtomObjectPtrInTuple(
                model_number, ptnr1_asym_id, ptnr1_comp_id, ptnr1_seq_id, ptnr1_atom_id) };
            auto atom_2{ m_data_block->GetAtomObjectPtrInTuple(
                model_number, ptnr2_asym_id, ptnr2_comp_id, ptnr2_seq_id, ptnr2_atom_id) };
            if (atom_1 == nullptr || atom_2 == nullptr)
            {
                Logger::Log(LogLevel::Warning, "Cannot find atom object for connection: " + conn_id);
                return;
            }

            auto bond_object{ std::make_unique<BondObject>(atom_1, atom_2) };
            bond_object->SetBondKey(bond_key);
            bond_object->SetBondType(bond_entry.bond_type);
            bond_object->SetBondOrder(bond_entry.bond_order);
            bond_object->SetSpecialBondFlag(atom_1->GetSpecialAtomFlag());
            m_data_block->AddBondObject(std::move(bond_object));
        }
    );
}

void CifFormat::LoadStructureSheetBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadStructureSheetBlock() called");
    infile.clear();
    infile.seekg(0);
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

void CifFormat::LoadAtomSiteBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadAtomSiteBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_atom_site.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            static AtomObject * last_atom_object{ nullptr };
            auto group_type{ token_list[index_map.at("group_PDB")] };
            auto element_type{ token_list[index_map.at("type_symbol")] };
            auto comp_id{ token_list[index_map.at("label_comp_id")] };
            auto atom_id{ token_list[index_map.at("label_atom_id")] };
            auto indicator{ token_list[index_map.at("label_alt_id")] };
            auto sequence_id{ token_list[index_map.at("label_seq_id")] };
            auto serial_id{ token_list[index_map.at("id")] };
            auto chain_id{ token_list[index_map.at("label_asym_id")] };
            auto position_x{ std::stof(token_list[index_map.at("Cartn_x")]) };
            auto position_y{ std::stof(token_list[index_map.at("Cartn_y")]) };
            auto position_z{ std::stof(token_list[index_map.at("Cartn_z")]) };
            auto occupancy{ std::stof(token_list[index_map.at("occupancy")]) };
            auto temperature{ std::stof(token_list[index_map.at("B_iso_or_equiv")]) };
            auto model_number{ token_list[index_map.at("pdbx_PDB_model_num")] };
            if (element_type == "H") return; // Skip hydrogen atom
            auto is_special_atom{ (group_type == "HETATM") ? true : false };
            StringHelper::EraseCharFromString(atom_id, '\"');

            if (m_find_chemical_component_entry == false)
            {
                BuildDefaultChemicalComponentEntry(comp_id);
            }

            if (m_find_component_atom_entry == false)
            {
                BuildDefaultComponentAtomEntry(comp_id, atom_id, element_type);
            }
            auto atom_object{ std::make_unique<AtomObject>() };
            atom_object->SetComponentID(comp_id);
            atom_object->SetComponentKey(m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id));
            atom_object->SetAtomID(atom_id);
            atom_object->SetAtomKey(m_data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id));
            atom_object->SetElement(element_type);
            atom_object->SetIndicator(indicator);
            atom_object->SetSequenceID((sequence_id == ".") ? -1 : std::stoi(sequence_id));
            atom_object->SetSerialID(std::stoi(serial_id));
            atom_object->SetChainID(chain_id);
            atom_object->SetPosition(position_x, position_y, position_z);
            atom_object->SetOccupancy(occupancy);
            atom_object->SetTemperature(temperature);
            atom_object->SetSpecialAtomFlag(is_special_atom);
            m_data_block->SetStructureInfo(atom_object.get());

            auto model_number_id{ std::stoi(model_number) };
            if (indicator == ".")
            {
                last_atom_object = nullptr;
                m_data_block->AddAtomObject(model_number_id, std::move(atom_object));
            }
            else if (indicator == "A")
            {
                last_atom_object = atom_object.get();
                m_data_block->AddAtomObject(model_number_id, std::move(atom_object));
            }
            else
            {
                if (last_atom_object == nullptr)
                {
                    Logger::Log(LogLevel::Error,
                        "CifFormat::LoadAtomSiteBlock() atom_object is missing.");
                    return;
                }
                last_atom_object->AddAlternatePosition(indicator, {position_x, position_y, position_z});
                last_atom_object->AddAlternateOccupancy(indicator, occupancy);
                last_atom_object->AddAlternateTemperature(indicator, temperature);
            }
        }
    );
}

void CifFormat::ConstructBondList(void)
{
    Logger::Log(LogLevel::Debug, "CifFormat::ConstructBondList() called");
    BuildPepetideBondEntry();
    BuildPhosphodiesterBondEntry();
    auto model_number{ 1 };
    auto & atom_object_list{ m_data_block->GetAtomObjectMap().at(model_number) };
    std::vector<AtomObject *> atom_ptr_list;
    atom_ptr_list.reserve(atom_object_list.size());
    for (auto & atom : atom_object_list)
    {
        atom_ptr_list.emplace_back(atom.get());
    }
    auto kd_tree_root{ KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0) };
    for (auto & atom : atom_object_list)
    {
        auto component_id_1{ atom->GetComponentID() };
        auto atom_id_1{ atom->GetAtomID() };
        auto sequence_id_1{ atom->GetSequenceID() };
        auto chain_id_1{ atom->GetChainID() };
        auto neighbor_atom_list{
            KDTreeAlgorithm<AtomObject>::RangeSearch(
                kd_tree_root.get(), atom.get(), m_bond_searching_radius)
        };

        for (auto neighbor_atom : neighbor_atom_list)
        {
            if (neighbor_atom == atom.get()) continue;
            auto component_id_2{ neighbor_atom->GetComponentID() };
            auto atom_id_2{ neighbor_atom->GetAtomID() };
            auto sequence_id_2{ neighbor_atom->GetSequenceID() };
            auto chain_id_2{ neighbor_atom->GetChainID() };
            auto bond_key_system{ m_data_block->GetBondKeySystemPtr() };
            if (bond_key_system->IsRegistedBond(atom_id_1, atom_id_2) == false) continue;
            auto component_key_1{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(component_id_1) };
            auto bond_key{ bond_key_system->GetBondKey(atom_id_1, atom_id_2) };
            if (m_data_block->HasComponentBondEntry(component_key_1, bond_key) == false) continue;

            bool is_in_same_component{ (component_id_1 == component_id_2) };
            bool is_in_same_chain{ (chain_id_1 == chain_id_2) };
            bool is_in_consecutive_sequence{ (sequence_id_1 + 1 == sequence_id_2) };

            // Peptide bond C-N between consecutive residues in the same chain
            bool is_peptide_bond{
                !is_in_same_component &&
                is_in_same_chain &&
                is_in_consecutive_sequence &&
                (bond_key == static_cast<BondKey>(Link::C_N))
            };

            // Phosphodiester bond P-O3' between consecutive residues in the same chain
            bool is_phosphodiester_bond{
                !is_in_same_component &&
                is_in_same_chain &&
                is_in_consecutive_sequence &&
                (bond_key == static_cast<BondKey>(Link::P_O3p))
            };

            if (is_in_same_component == false &&
                is_peptide_bond == false &&
                is_phosphodiester_bond == false) continue;

            auto bond_entry{
                m_data_block->GetComponentBondEntryPtr(component_key_1, bond_key)
            };
            
            auto bond_object{ std::make_unique<BondObject>(atom.get(), neighbor_atom) };
            bond_object->SetBondKey(bond_key);
            bond_object->SetBondType(bond_entry->bond_type);
            bond_object->SetBondOrder(bond_entry->bond_order);
            bond_object->SetSpecialBondFlag(false);
            m_data_block->AddBondObject(std::move(bond_object));
        }
    }
    Logger::Log(LogLevel::Info,
        "Construct " + std::to_string(m_data_block->GetBondObjectList().size()) + " bonds.");
}

AtomicModelDataBlock * CifFormat::GetDataBlockPtr(void)
{
    return m_data_block.get();
}

void CifFormat::SaveHeader(const ModelObject * model_object, std::ostream & stream)
{
    Logger::Log(LogLevel::Debug, "CifFormat::SaveHeader() called");
    if (model_object == nullptr) return;
    stream << "data_" << model_object->GetKeyTag() << '\n';
    stream << "#\n";
}

void CifFormat::SaveDataArray(
    const ModelObject * model_object, std::ostream & stream, int model_par)
{
    Logger::Log(LogLevel::Debug, "CifFormat::SaveDataArray() called");
    if (model_object == nullptr) return;

    WriteAtomSiteBlock(model_object, stream, model_par);
}

void CifFormat::WriteAtomSiteBlock(
    const ModelObject * model_object, std::ostream & stream, int model_par)
{
    Logger::Log(LogLevel::Debug, "CifFormat::WriteAtomSiteBlock() called");
    stream << "loop_\n";
    stream << "_atom_site.group_PDB\n";
    stream << "_atom_site.id\n";
    stream << "_atom_site.type_symbol\n";
    stream << "_atom_site.label_atom_id\n";
    stream << "_atom_site.label_alt_id\n";
    stream << "_atom_site.label_comp_id\n";
    stream << "_atom_site.label_asym_id\n";
    stream << "_atom_site.label_seq_id\n";
    stream << "_atom_site.Cartn_x\n";
    stream << "_atom_site.Cartn_y\n";
    stream << "_atom_site.Cartn_z\n";
    stream << "_atom_site.occupancy\n";
    stream << "_atom_site.B_iso_or_equiv\n";
    stream << "_atom_site.pdbx_PDB_model_num\n";

    const int model_number{ 1 };
    for (const auto & atom_ptr : model_object->GetAtomList())
    {
        const AtomObject * atom{ atom_ptr.get() };
        if (atom->GetLocalPotentialEntry() == nullptr) continue;
        auto model_entry{ atom->GetLocalPotentialEntry() };
        auto gaus_estimate{ model_entry->GetGausEstimateMDPDE(model_par) };
        auto position{ atom->GetPosition() };
        WriteAtomSiteBlockEntry(
            atom, position, atom->GetIndicator(), atom->GetOccupancy(),
            static_cast<float>(gaus_estimate), model_number, stream
        );
    }
    stream << "#\n";
}

void CifFormat::WriteAtomSiteBlockEntry(
    const AtomObject * atom,
    const std::array<float, 3> & position,
    const std::string & alt_id,
    float occupancy,
    float temperature,
    int model_number,
    std::ostream & stream)
{
    std::string group_PDB{ atom->GetSpecialAtomFlag() ? "HETATM" : "ATOM" };
    std::string type_symbol{ ChemicalDataHelper::GetLabel(atom->GetElement()) };
    std::string label_atom_id{ atom->GetAtomID() };
    std::string label_alt_id{ alt_id.empty() ? "." : alt_id };
    std::string label_comp_id{ atom->GetComponentID() };
    std::string label_asym_id{ atom->GetChainID() };
    std::string label_seq_id{
        (atom->GetSequenceID() == -1) ? "." : std::to_string(atom->GetSequenceID())
    };

    stream << group_PDB << ' '
           << atom->GetSerialID() << ' '
           << type_symbol << ' '
           << label_atom_id << ' '
           << label_alt_id << ' '
           << label_comp_id << ' '
           << label_asym_id << ' '
           << label_seq_id << ' '
           << std::fixed << std::setprecision(3)
           << position[0] << ' ' << position[1] << ' ' << position[2] << ' '
           << occupancy << ' ' << temperature << ' '
           << model_number << '\n';
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
        StringHelper::StripCarriageReturn(line);
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
            auto initial{ StringHelper::SplitStringLineAsTokens(line, expected_column_size) };
            token_list.insert(token_list.end(), initial.begin(), initial.end());

            // now read continuation lines until we have all fields
            bool in_multiline{ false };
            std::string multiline_content;
            std::string next_line;
            while (token_list.size() < expected_column_size && std::getline(infile, next_line))
            {
                StringHelper::StripCarriageReturn(next_line);
                if (next_line.empty() == true) continue;
                if (in_multiline == false && next_line[0] == ';')
                {
                    // start of multiline literal for the next field
                    in_multiline = true;
                    multiline_content.clear();
                    // strip leading semicolon
                    auto content{ next_line.substr(1) };
                    if (content.empty() == false) multiline_content = content;
                }
                else if (in_multiline == true)
                {
                    if (next_line.empty() == false && next_line[0] == ';')
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
                    auto more{ StringHelper::SplitStringLineAsTokens(next_line, expected_column_size) };
                    token_list.insert(token_list.end(), more.begin(), more.end());
                }
            }

            table_handler(column_index_map, token_list);
        }
    }
}

void CifFormat::BuildDefaultChemicalComponentEntry(const std::string & comp_id)
{
    Logger::Log(LogLevel::Debug, "CifFormat::BuildDefaultChemicalComponentEntry() called");
    auto entry{ std::make_unique<ChemicalComponentEntry>() };
    entry->SetComponentId(comp_id);
    entry->SetComponentName("");
    entry->SetComponentType("");
    entry->SetComponentFormula("");
    entry->SetComponentMolecularWeight(0.0f);
    auto standard_flag{ false };
    if (ChemicalDataHelper::GetResidueFromString(comp_id) != Residue::UNK) standard_flag = true;
    entry->SetStandardMonomerFlag(standard_flag);

    m_data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
    auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
    m_data_block->AddChemicalComponentEntry(component_key, std::move(entry));
}

void CifFormat::BuildDefaultComponentAtomEntry(
    const std::string & comp_id,
    const std::string & atom_id,
    const std::string & element_symbol)
{
    Logger::Log(LogLevel::Debug, "CifFormat::BuildDefaultComponentAtomEntry() called");
    auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };

    m_data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
    auto atom_key{ m_data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id) };

    ComponentAtomEntry atom_entry;
    atom_entry.atom_id = atom_id;
    atom_entry.element_type = ChemicalDataHelper::GetElementFromString(element_symbol);
    atom_entry.aromatic_atom_flag = false;
    atom_entry.stereo_config = StereoChemistry::NONE;
    
    m_data_block->AddComponentAtomEntry(component_key, atom_key, atom_entry);
}

void CifFormat::BuildDefaultComponentBondEntry(void)
{
    Logger::Log(LogLevel::Debug, "CifFormat::BuildDefaultComponentBondEntry() called");
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto comp_id{ ChemicalDataHelper::GetLabel(residue) };
        auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
        for (auto & link : ComponentHelper::GetLinkList(residue))
        {
            auto bond_id{ ChemicalDataHelper::GetLabel(link) };
            m_data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
            auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(bond_id) };

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = BondType::COVALENT;
            bond_entry.bond_order = BondOrder::UNK;
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;
            
            m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
        }
    }
}

void CifFormat::BuildPepetideBondEntry(void)
{
    Logger::Log(LogLevel::Debug, "CifFormat::BuildPepetideBondEntry() called");
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto comp_id{ ChemicalDataHelper::GetLabel(residue) };
        auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
        if (m_data_block->HasChemicalComponentEntry(component_key) == false) continue;
        auto bond_id{ ChemicalDataHelper::GetLabel(Link::C_N) };
        m_data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
        auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(bond_id) };

        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = BondType::COVALENT;
        bond_entry.bond_order = BondOrder::SINGLE;
        bond_entry.aromatic_atom_flag = false;
        bond_entry.stereo_config = StereoChemistry::NONE;
        
        m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
    }
}

void CifFormat::BuildPhosphodiesterBondEntry(void)
{
    Logger::Log(LogLevel::Debug, "CifFormat::BuildPhosphodiesterBondEntry() called");
    for (auto & residue : ChemicalDataHelper::GetStandardNucleotideList())
    {
        auto comp_id{ ChemicalDataHelper::GetLabel(residue) };
        auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id) };
        if (m_data_block->HasChemicalComponentEntry(component_key) == false) continue;
        auto bond_id{ ChemicalDataHelper::GetLabel(Link::P_O3p) };
        m_data_block->GetBondKeySystemPtr()->RegisterBond(bond_id);
        auto bond_key{ m_data_block->GetBondKeySystemPtr()->GetBondKey(bond_id) };

        ComponentBondEntry bond_entry;
        bond_entry.bond_id = bond_id;
        bond_entry.bond_type = BondType::COVALENT;
        bond_entry.bond_order = BondOrder::SINGLE;
        bond_entry.aromatic_atom_flag = false;
        bond_entry.stereo_config = StereoChemistry::NONE;
        
        m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);
    }
}
