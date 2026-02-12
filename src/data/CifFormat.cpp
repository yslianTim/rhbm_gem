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
#include <optional>
#include <initializer_list>
#include <map>

namespace
{
bool IsMmCifMissingValue(const std::string & value)
{
    return value.empty() || value == "." || value == "?";
}

std::vector<std::string> SplitMmCifTokens(const std::string & line)
{
    std::vector<std::string> token_list;
    for (size_t pos = 0; pos < line.size();)
    {
        while (pos < line.size() && std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
        if (pos >= line.size()) break;

        const char first_char{ line[pos] };
        if (first_char == '\'' || first_char == '"')
        {
            const char quote_char{ first_char };
            ++pos;
            const auto start{ pos };
            while (pos < line.size() && line[pos] != quote_char) ++pos;
            token_list.emplace_back(line.substr(start, pos - start));
            if (pos < line.size()) ++pos;
            continue;
        }

        const auto start{ pos };
        while (pos < line.size() && !std::isspace(static_cast<unsigned char>(line[pos]))) ++pos;
        token_list.emplace_back(line.substr(start, pos - start));
    }
    return token_list;
}

std::string BuildMmCifTokenPreview(const std::vector<std::string> & token_list, size_t max_items = 8)
{
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < token_list.size() && i < max_items; ++i)
    {
        if (i > 0) oss << ", ";
        oss << token_list[i];
    }
    if (token_list.size() > max_items) oss << ", ...";
    oss << "]";
    return oss.str();
}

std::optional<std::string> GetTokenOptional(
    const std::unordered_map<std::string, size_t> & index_map,
    const std::vector<std::string> & token_list,
    std::initializer_list<std::string_view> name_candidates)
{
    for (const auto & name : name_candidates)
    {
        auto iter{ index_map.find(std::string{name}) };
        if (iter == index_map.end()) continue;
        if (iter->second >= token_list.size()) continue;
        return token_list[iter->second];
    }
    return std::nullopt;
}

int ParseIntOrDefault(
    const std::string & value,
    int default_value,
    const std::string & field_name,
    const std::string & log_context)
{
    if (IsMmCifMissingValue(value)) return default_value;
    try
    {
        return std::stoi(value);
    }
    catch (const std::exception &)
    {
        Logger::Log(LogLevel::Warning,
            log_context + " Invalid integer in " + field_name + ": " + value
            + ", fallback = " + std::to_string(default_value));
        return default_value;
    }
}

std::optional<float> TryParseFloat(const std::string & value)
{
    if (IsMmCifMissingValue(value)) return std::nullopt;
    try
    {
        return std::stof(value);
    }
    catch (const std::exception &)
    {
        return std::nullopt;
    }
}

float ParseFloatOrDefault(
    const std::string & value,
    float default_value,
    const std::string & field_name,
    const std::string & log_context)
{
    auto parsed_value{ TryParseFloat(value) };
    if (parsed_value.has_value()) return *parsed_value;
    if (!IsMmCifMissingValue(value))
    {
        Logger::Log(LogLevel::Warning,
            log_context + " Invalid float in " + field_name + ": " + value
            + ", fallback = " + std::to_string(default_value));
    }
    return default_value;
}

std::string BuildAtomAltLocKey(
    int model_number,
    const std::string & chain_id,
    const std::string & comp_id,
    const std::string & sequence_id,
    const std::string & atom_id)
{
    return std::to_string(model_number) + "|" + chain_id + "|" + comp_id + "|" + sequence_id + "|" + atom_id;
}

std::optional<std::string> ParseMmCifDataItemValue(
    std::ifstream & infile, const std::string & line, std::string_view key)
{
    if (line.rfind(key, 0) != 0) return std::nullopt;
    auto token_list{ SplitMmCifTokens(line) };
    if (token_list.size() >= 2) return token_list[1];
    if (token_list.size() != 1) return std::string{};

    std::string next_line;
    while (std::getline(infile, next_line))
    {
        StringHelper::StripCarriageReturn(next_line);
        if (next_line.empty()) continue;
        if (next_line[0] == '#') return std::string{};
        if (next_line[0] != ';')
        {
            auto next_tokens{ SplitMmCifTokens(next_line) };
            return next_tokens.empty() ? std::string{} : next_tokens.front();
        }

        // Multi-line value starts with ';' and ends at the next line that starts with ';'.
        std::string value{ next_line.substr(1) };
        while (std::getline(infile, next_line))
        {
            StringHelper::StripCarriageReturn(next_line);
            if (next_line.empty() == false && next_line[0] == ';') break;
            if (!value.empty()) value += "\n";
            value += next_line;
        }
        return value;
    }
    return std::string{};
}
} // namespace

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
    std::unordered_map<std::string, std::string> data_map;
    ParseLoopBlock(infile, "_database_2.",
        [&data_map](const std::unordered_map<std::string, size_t> & index_map,
                    const std::vector<std::string> & token_list)
        {
            auto database_id{ GetTokenOptional(index_map, token_list, {"database_id"}) };
            auto database_code{ GetTokenOptional(index_map, token_list, {"database_code"}) };
            if (!database_id.has_value() || !database_code.has_value()) return;
            if (IsMmCifMissingValue(*database_id) || IsMmCifMissingValue(*database_code)) return;
            data_map[*database_id] = *database_code;
        }
    );
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

void CifFormat::LoadEntityBlock(std::ifstream & infile)
{
    Logger::Log(LogLevel::Debug, "CifFormat::LoadEntityBlock() called");
    infile.clear();
    infile.seekg(0);
    ParseLoopBlock(infile, "_entity.",
        [this](const std::unordered_map<std::string, size_t> & index_map,
               const std::vector<std::string> & token_list)
        {
            if (index_map.find("id") == index_map.end() ||
                index_map.find("type") == index_map.end())
            {
                return;
            }
            if (index_map.at("id") >= token_list.size() ||
                index_map.at("type") >= token_list.size())
            {
                return;
            }

            auto entity_id{ token_list[index_map.at("id")] };
            auto entity_type{ token_list[index_map.at("type")] };
            std::string molecules_size_string{ "-1" };
            if (index_map.find("pdbx_number_of_molecules") != index_map.end() &&
                index_map.at("pdbx_number_of_molecules") < token_list.size())
            {
                molecules_size_string = token_list[index_map.at("pdbx_number_of_molecules")];
            }
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
            if (index_map.find("entity_id") == index_map.end() ||
                index_map.find("id") == index_map.end())
            {
                return;
            }
            if (index_map.at("entity_id") >= token_list.size() ||
                index_map.at("id") >= token_list.size())
            {
                return;
            }
            auto entity_id{ token_list[index_map.at("entity_id")] };
            auto chain_id{ token_list[index_map.at("id")] };
            m_data_block->AddChainIDInEntityMap(entity_id, chain_id);
        }
    );

    const auto need_entity_fallback{ m_data_block->GetEntityTypeMap().empty() };
    const auto need_chain_fallback{ m_data_block->GetChainIDListMap().empty() };
    if (!need_entity_fallback && !need_chain_fallback) return;

    Logger::Log(LogLevel::Warning,
        "Cannot parse complete entity/chain metadata from loop blocks. "
        "Trying key-value fallback for _entity/_struct_asym.");

    infile.clear();
    infile.seekg(0);

    std::string entity_id;
    std::string entity_type;
    std::string molecules_size_raw;
    std::string struct_asym_id;
    std::string struct_asym_entity_id;
    std::string line;
    while (std::getline(infile, line))
    {
        StringHelper::StripCarriageReturn(line);

        if (need_entity_fallback)
        {
            if (auto value{ ParseMmCifDataItemValue(infile, line, "_entity.id") })
            {
                entity_id = *value;
                continue;
            }
            if (auto value{ ParseMmCifDataItemValue(infile, line, "_entity.type") })
            {
                entity_type = *value;
                continue;
            }
            if (auto value{ ParseMmCifDataItemValue(infile, line, "_entity.pdbx_number_of_molecules") })
            {
                molecules_size_raw = *value;
                continue;
            }
        }

        if (need_chain_fallback)
        {
            if (auto value{ ParseMmCifDataItemValue(infile, line, "_struct_asym.id") })
            {
                struct_asym_id = *value;
                continue;
            }
            if (auto value{ ParseMmCifDataItemValue(infile, line, "_struct_asym.entity_id") })
            {
                struct_asym_entity_id = *value;
                continue;
            }
        }
    }

    if (need_entity_fallback)
    {
        if (IsMmCifMissingValue(entity_id) || IsMmCifMissingValue(entity_type))
        {
            Logger::Log(LogLevel::Warning,
                "Key-value fallback cannot recover required _entity.id/_entity.type fields.");
        }
        else
        {
            int molecules_size{ -1 };
            if (!IsMmCifMissingValue(molecules_size_raw))
            {
                try
                {
                    molecules_size = std::stoi(molecules_size_raw);
                }
                catch (const std::exception &)
                {
                    Logger::Log(LogLevel::Warning,
                        "Invalid _entity.pdbx_number_of_molecules in key-value fallback: "
                        + molecules_size_raw);
                }
            }
            m_data_block->AddEntityTypeInEntityMap(
                entity_id, ChemicalDataHelper::GetEntityFromString(entity_type));
            m_data_block->AddMoleculesSizeInEntityMap(entity_id, molecules_size);
        }
    }

    if (need_chain_fallback)
    {
        if (IsMmCifMissingValue(struct_asym_id) || IsMmCifMissingValue(struct_asym_entity_id))
        {
            Logger::Log(LogLevel::Warning,
                "Key-value fallback cannot recover required _struct_asym.id/_struct_asym.entity_id fields.");
        }
        else
        {
            m_data_block->AddChainIDInEntityMap(struct_asym_entity_id, struct_asym_id);
        }
    }
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
            auto conn_id{ GetTokenOptional(index_map, token_list, {"id"}) };
            auto conn_type_id{ GetTokenOptional(index_map, token_list, {"conn_type_id"}) };
            auto ptnr1_asym_id{ GetTokenOptional(index_map, token_list, {"ptnr1_label_asym_id"}) };
            auto ptnr1_comp_id{ GetTokenOptional(index_map, token_list, {"ptnr1_label_comp_id"}) };
            auto ptnr1_seq_id{ GetTokenOptional(index_map, token_list, {"ptnr1_label_seq_id"}) };
            auto ptnr1_atom_id{ GetTokenOptional(index_map, token_list, {"ptnr1_label_atom_id"}) };
            auto ptnr2_asym_id{ GetTokenOptional(index_map, token_list, {"ptnr2_label_asym_id"}) };
            auto ptnr2_comp_id{ GetTokenOptional(index_map, token_list, {"ptnr2_label_comp_id"}) };
            auto ptnr2_seq_id{ GetTokenOptional(index_map, token_list, {"ptnr2_label_seq_id"}) };
            auto ptnr2_atom_id{ GetTokenOptional(index_map, token_list, {"ptnr2_label_atom_id"}) };
            auto value_order_str{
                GetTokenOptional(index_map, token_list, {"pdbx_value_order"})
            };

            if (!conn_type_id.has_value() ||
                !ptnr1_asym_id.has_value() || !ptnr1_comp_id.has_value() ||
                !ptnr1_seq_id.has_value() || !ptnr1_atom_id.has_value() ||
                !ptnr2_asym_id.has_value() || !ptnr2_comp_id.has_value() ||
                !ptnr2_seq_id.has_value() || !ptnr2_atom_id.has_value())
            {
                return;
            }
            auto conn_id_label{ conn_id.value_or("UNKNOWN_CONN") };
            auto ptnr1_atom_id_value{ *ptnr1_atom_id };
            auto ptnr2_atom_id_value{ *ptnr2_atom_id };
            auto ptnr1_seq_id_value{ IsMmCifMissingValue(*ptnr1_seq_id) ? "." : *ptnr1_seq_id };
            auto ptnr2_seq_id_value{ IsMmCifMissingValue(*ptnr2_seq_id) ? "." : *ptnr2_seq_id };

            StringHelper::EraseCharFromString(ptnr1_atom_id_value, '\"');
            StringHelper::EraseCharFromString(ptnr2_atom_id_value, '\"');

            m_data_block->GetBondKeySystemPtr()->RegisterBond(ptnr1_atom_id_value, ptnr2_atom_id_value);
            auto component_key{ m_data_block->GetComponentKeySystemPtr()->GetComponentKey(*ptnr1_comp_id) };
            auto bond_key{
                m_data_block->GetBondKeySystemPtr()->GetBondKey(ptnr1_atom_id_value, ptnr2_atom_id_value)
            };
            auto bond_id{ m_data_block->GetBondKeySystemPtr()->GetBondId(bond_key) };
            auto bond_type{ ChemicalDataHelper::GetBondTypeFromString(*conn_type_id) };
            if (bond_type == BondType::HYDROGEN) return; // Skip hydrogen bond

            ComponentBondEntry bond_entry;
            bond_entry.bond_id = bond_id;
            bond_entry.bond_type = bond_type;
            if (!value_order_str.has_value() || IsMmCifMissingValue(*value_order_str))
            {
                bond_entry.bond_order = BondOrder::UNK;
            }
            else
            {
                bond_entry.bond_order = ChemicalDataHelper::GetBondOrderFromString(*value_order_str);
            }
            bond_entry.aromatic_atom_flag = false;
            bond_entry.stereo_config = StereoChemistry::NONE;
            m_data_block->AddComponentBondEntry(component_key, bond_key, bond_entry);

            auto ptnr1_model_num_token{
                GetTokenOptional(index_map, token_list, {"pdbx_ptnr1_PDB_model_num", "ptnr1_PDB_model_num"})
            };
            auto ptnr2_model_num_token{
                GetTokenOptional(index_map, token_list, {"pdbx_ptnr2_PDB_model_num", "ptnr2_PDB_model_num"})
            };
            auto model_number_1{
                ParseIntOrDefault(ptnr1_model_num_token.value_or("1"), 1,
                    "_struct_conn.pdbx_ptnr1_PDB_model_num",
                    "LoadStructureConnectionBlock()[" + conn_id_label + "]")
            };
            auto model_number_2{
                ParseIntOrDefault(ptnr2_model_num_token.value_or(std::to_string(model_number_1)), model_number_1,
                    "_struct_conn.pdbx_ptnr2_PDB_model_num",
                    "LoadStructureConnectionBlock()[" + conn_id_label + "]")
            };

            auto atom_1{ m_data_block->GetAtomObjectPtrInTuple(
                model_number_1, *ptnr1_asym_id, *ptnr1_comp_id, ptnr1_seq_id_value, ptnr1_atom_id_value) };
            auto atom_2{ m_data_block->GetAtomObjectPtrInTuple(
                model_number_2, *ptnr2_asym_id, *ptnr2_comp_id, ptnr2_seq_id_value, ptnr2_atom_id_value) };
            if (atom_1 == nullptr)
            {
                atom_1 = m_data_block->GetAtomObjectPtrInAnyModel(
                    *ptnr1_asym_id, *ptnr1_comp_id, ptnr1_seq_id_value, ptnr1_atom_id_value, &model_number_1);
            }
            if (atom_2 == nullptr)
            {
                atom_2 = m_data_block->GetAtomObjectPtrInAnyModel(
                    *ptnr2_asym_id, *ptnr2_comp_id, ptnr2_seq_id_value, ptnr2_atom_id_value, &model_number_2);
            }
            if (atom_1 == nullptr || atom_2 == nullptr)
            {
                Logger::Log(LogLevel::Warning, "Cannot find atom object for connection: " + conn_id_label);
                return;
            }
            if (model_number_1 != model_number_2)
            {
                Logger::Log(LogLevel::Warning,
                    "Skip cross-model connection [" + conn_id_label + "] between model "
                    + std::to_string(model_number_1) + " and model "
                    + std::to_string(model_number_2) + ".");
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
    std::map<std::string, AtomObject *> altloc_primary_atom_map;
    size_t atom_site_row_count{ 0 };
    size_t atom_site_skip_count{ 0 };
    ParseLoopBlock(infile, "_atom_site.",
        [this, &altloc_primary_atom_map, &atom_site_row_count, &atom_site_skip_count](
            const std::unordered_map<std::string, size_t> & index_map,
            const std::vector<std::string> & token_list)
        {
            ++atom_site_row_count;
            const std::string context{
                "LoadAtomSiteBlock[row " + std::to_string(atom_site_row_count) + "]"
            };

            auto group_type{ GetTokenOptional(index_map, token_list, {"group_PDB"}).value_or("ATOM") };
            auto element_type{ GetTokenOptional(index_map, token_list, {"type_symbol"}).value_or("UNK") };
            auto comp_id{
                GetTokenOptional(index_map, token_list, {"label_comp_id", "auth_comp_id"}).value_or("UNK")
            };
            auto atom_id_opt{ GetTokenOptional(index_map, token_list, {"label_atom_id", "auth_atom_id"}) };
            auto indicator{
                GetTokenOptional(index_map, token_list, {"label_alt_id"}).value_or(".")
            };
            auto sequence_id{
                GetTokenOptional(index_map, token_list, {"label_seq_id", "auth_seq_id"}).value_or(".")
            };
            auto serial_id{ GetTokenOptional(index_map, token_list, {"id"}).value_or(
                std::to_string(atom_site_row_count)) };
            auto chain_id{
                GetTokenOptional(index_map, token_list, {"label_asym_id", "auth_asym_id"}).value_or("?")
            };
            auto position_x_str{ GetTokenOptional(index_map, token_list, {"Cartn_x"}) };
            auto position_y_str{ GetTokenOptional(index_map, token_list, {"Cartn_y"}) };
            auto position_z_str{ GetTokenOptional(index_map, token_list, {"Cartn_z"}) };
            auto occupancy_str{ GetTokenOptional(index_map, token_list, {"occupancy"}) };
            auto temperature_str{ GetTokenOptional(index_map, token_list, {"B_iso_or_equiv"}) };
            auto model_number_str{ GetTokenOptional(index_map, token_list, {"pdbx_PDB_model_num"}) };

            if (!atom_id_opt.has_value())
            {
                Logger::Log(LogLevel::Warning, context + " Missing atom id, skip this row.");
                ++atom_site_skip_count;
                return;
            }

            if (IsMmCifMissingValue(indicator)) indicator = ".";
            auto parsed_x{ TryParseFloat(position_x_str.value_or("?")) };
            auto parsed_y{ TryParseFloat(position_y_str.value_or("?")) };
            auto parsed_z{ TryParseFloat(position_z_str.value_or("?")) };
            if (!parsed_x.has_value() || !parsed_y.has_value() || !parsed_z.has_value())
            {
                Logger::Log(LogLevel::Warning,
                    context + " Invalid/missing Cartesian coordinates, skip this row.");
                ++atom_site_skip_count;
                return;
            }

            auto occupancy{
                ParseFloatOrDefault(occupancy_str.value_or("?"), 1.0f, "_atom_site.occupancy", context)
            };
            auto temperature{
                ParseFloatOrDefault(
                    temperature_str.value_or("?"), 0.0f, "_atom_site.B_iso_or_equiv", context)
            };
            auto model_number_id{
                ParseIntOrDefault(
                    model_number_str.value_or("1"), 1, "_atom_site.pdbx_PDB_model_num", context)
            };
            auto sequence_id_value{
                ParseIntOrDefault(sequence_id, -1, "_atom_site.(label/auth)_seq_id", context)
            };
            auto serial_id_value{ ParseIntOrDefault(serial_id, static_cast<int>(atom_site_row_count),
                "_atom_site.id", context) };

            auto atom_id{ *atom_id_opt };
            StringHelper::EraseCharFromString(atom_id, '\"');

            auto element_enum{ ChemicalDataHelper::GetElementFromString(element_type) };
            if (element_enum == Element::HYDROGEN) return; // Skip hydrogen atom

            auto is_special_atom{ (group_type == "HETATM") ? true : false };

            if (m_find_chemical_component_entry == false)
            {
                BuildDefaultChemicalComponentEntry(comp_id);
            }
            m_data_block->GetComponentKeySystemPtr()->RegisterComponent(comp_id);
            auto component_key{
                m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id)
            };
            if (m_data_block->HasChemicalComponentEntry(component_key) == false)
            {
                BuildDefaultChemicalComponentEntry(comp_id);
                component_key = m_data_block->GetComponentKeySystemPtr()->GetComponentKey(comp_id);
            }

            if (m_find_component_atom_entry == false)
            {
                BuildDefaultComponentAtomEntry(comp_id, atom_id, element_type);
            }
            m_data_block->GetAtomKeySystemPtr()->RegisterAtom(atom_id);
            auto atom_key{ m_data_block->GetAtomKeySystemPtr()->GetAtomKey(atom_id) };

            auto atom_altloc_key{
                BuildAtomAltLocKey(model_number_id, chain_id, comp_id, sequence_id, atom_id)
            };
            auto add_atom_as_primary{
                [&](std::unique_ptr<AtomObject> atom_object) -> AtomObject *
                {
                    auto raw_atom_ptr{ atom_object.get() };
                    m_data_block->AddAtomObject(model_number_id, std::move(atom_object));
                    altloc_primary_atom_map[atom_altloc_key] = raw_atom_ptr;
                    return raw_atom_ptr;
                }
            };

            if (indicator != "." && indicator != "A")
            {
                auto primary_iter{ altloc_primary_atom_map.find(atom_altloc_key) };
                if (primary_iter != altloc_primary_atom_map.end())
                {
                    primary_iter->second->AddAlternatePosition(indicator, {*parsed_x, *parsed_y, *parsed_z});
                    primary_iter->second->AddAlternateOccupancy(indicator, occupancy);
                    primary_iter->second->AddAlternateTemperature(indicator, temperature);
                    return;
                }
            }

            auto atom_object{ std::make_unique<AtomObject>() };
            atom_object->SetComponentID(comp_id);
            atom_object->SetComponentKey(component_key);
            atom_object->SetAtomID(atom_id);
            atom_object->SetAtomKey(atom_key);
            atom_object->SetElement(element_type);
            atom_object->SetIndicator(indicator);
            atom_object->SetSequenceID(sequence_id_value);
            atom_object->SetSerialID(serial_id_value);
            atom_object->SetChainID(chain_id);
            atom_object->SetPosition(*parsed_x, *parsed_y, *parsed_z);
            atom_object->SetOccupancy(occupancy);
            atom_object->SetTemperature(temperature);
            atom_object->SetSpecialAtomFlag(is_special_atom);
            m_data_block->SetStructureInfo(atom_object.get());

            if (indicator == "." || indicator == "A")
            {
                add_atom_as_primary(std::move(atom_object));
                return;
            }

            // First alternate indicator is not "A" (e.g. only "B"): treat as primary.
            add_atom_as_primary(std::move(atom_object));
        }
    );
    Logger::Log(LogLevel::Info,
        "LoadAtomSiteBlock parsed rows = " + std::to_string(atom_site_row_count)
        + ", skipped rows = " + std::to_string(atom_site_skip_count) + ".");
}

void CifFormat::ConstructBondList(void)
{
    Logger::Log(LogLevel::Debug, "CifFormat::ConstructBondList() called");
    BuildPepetideBondEntry();
    BuildPhosphodiesterBondEntry();
    const auto bond_count_before{ m_data_block->GetBondObjectList().size() };
    const auto & atom_object_map{ m_data_block->GetAtomObjectMap() };
    for (const auto & [model_number, atom_object_list] : atom_object_map)
    {
        if (atom_object_list.empty()) continue;
        std::vector<AtomObject *> atom_ptr_list;
        atom_ptr_list.reserve(atom_object_list.size());
        for (const auto & atom : atom_object_list)
        {
            atom_ptr_list.emplace_back(atom.get());
        }
        auto kd_tree_root{ KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0) };
        for (const auto & atom : atom_object_list)
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
                auto component_key_1{
                    m_data_block->GetComponentKeySystemPtr()->GetComponentKey(component_id_1)
                };
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
                if (bond_entry == nullptr) continue;
                
                auto bond_object{ std::make_unique<BondObject>(atom.get(), neighbor_atom) };
                bond_object->SetBondKey(bond_key);
                bond_object->SetBondType(bond_entry->bond_type);
                bond_object->SetBondOrder(bond_entry->bond_order);
                bond_object->SetSpecialBondFlag(false);
                m_data_block->AddBondObject(std::move(bond_object));
            }
        }
        Logger::Log(LogLevel::Debug,
            "ConstructBondList() processed model " + std::to_string(model_number)
            + " with " + std::to_string(atom_object_list.size()) + " atoms.");
    }
    Logger::Log(LogLevel::Info,
        "Construct " + std::to_string(m_data_block->GetBondObjectList().size() - bond_count_before)
        + " bonds (total " + std::to_string(m_data_block->GetBondObjectList().size()) + ").");
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
    size_t loop_row_number{ 0 };
    size_t file_line_number{ 0 };
    while (std::getline(infile, line))
    {
        ++file_line_number;
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
            if (line == "loop_" || line.rfind("data_", 0) == 0 || line[0] == '_') break;

            const auto expected_column_size{ column_index_map.size() };
            std::vector<std::string> token_list;
            token_list.reserve(expected_column_size);

            // initial tokens from this line
            auto initial{ SplitMmCifTokens(line) };
            token_list.insert(token_list.end(), initial.begin(), initial.end());

            // now read continuation lines until we have all fields
            bool in_multiline{ false };
            std::string multiline_content;
            std::string next_line;
            while (token_list.size() < expected_column_size && std::getline(infile, next_line))
            {
                ++file_line_number;
                StringHelper::StripCarriageReturn(next_line);
                if (next_line.empty() == true) continue;
                if (next_line == "loop_" || next_line.rfind("data_", 0) == 0 || next_line[0] == '_')
                {
                    break;
                }
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
                    auto more{ SplitMmCifTokens(next_line) };
                    token_list.insert(token_list.end(), more.begin(), more.end());
                }
            }

            ++loop_row_number;
            if (in_multiline)
            {
                Logger::Log(LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix)
                    + ") row " + std::to_string(loop_row_number)
                    + " has unterminated multiline token. Skip this row.");
                continue;
            }
            if (token_list.size() != expected_column_size)
            {
                Logger::Log(LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix)
                    + ") row " + std::to_string(loop_row_number)
                    + " has token count " + std::to_string(token_list.size())
                    + " but expects " + std::to_string(expected_column_size)
                    + ". tokens = " + BuildMmCifTokenPreview(token_list));
                continue;
            }

            try
            {
                table_handler(column_index_map, token_list);
            }
            catch (const std::exception & ex)
            {
                Logger::Log(LogLevel::Warning,
                    "ParseLoopBlock(" + std::string(data_block_prefix)
                    + ") row " + std::to_string(loop_row_number)
                    + " failed at file line " + std::to_string(file_line_number)
                    + ": " + std::string(ex.what())
                    + ". tokens = " + BuildMmCifTokenPreview(token_list));
                continue;
            }
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
