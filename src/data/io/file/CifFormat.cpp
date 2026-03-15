#include "internal/io/file/CifFormat.hpp"
#include "internal/io/file/MmCifLoopParser.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include "internal/object/AtomicModelDataBlock.hpp"
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <iomanip>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <ostream>
#include <algorithm>
#include <optional>
#include <initializer_list>
#include <map>
#include <fstream>
#include <unordered_set>
#include <cstdint>

namespace rhbm_gem {

namespace
{
bool IsMmCifMissingValue(const std::string & value)
{
    return value.empty() || value == "." || value == "?";
}

template <typename IndexMap>
std::optional<std::string> GetTokenOptional(
    const IndexMap & index_map,
    const std::vector<std::string> & token_list,
    std::initializer_list<std::string_view> name_candidates)
{
    for (const auto & name : name_candidates)
    {
        auto iter{ index_map.find(name) };
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

struct AtomAltLocKey
{
    int model_number;
    std::string chain_id;
    std::string comp_id;
    std::string sequence_id_token;
    std::string atom_id;

    bool operator==(const AtomAltLocKey & other) const
    {
        return model_number == other.model_number &&
               chain_id == other.chain_id &&
               comp_id == other.comp_id &&
               sequence_id_token == other.sequence_id_token &&
               atom_id == other.atom_id;
    }
};

struct AtomAltLocKeyHash
{
    size_t operator()(const AtomAltLocKey & key) const
    {
        size_t seed{ std::hash<int>{}(key.model_number) };
        seed ^= std::hash<std::string>{}(key.chain_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::string>{}(key.comp_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::string>{}(key.sequence_id_token) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<std::string>{}(key.atom_id) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

} // namespace

CifFormat::CifFormat() :
    m_data_block{ std::make_unique<AtomicModelDataBlock>() }
{
}

CifFormat::~CifFormat()
{
}

void CifFormat::ResetReadState()
{
    m_data_block = std::make_unique<AtomicModelDataBlock>();
    m_find_chemical_component_entry = false;
    m_find_component_atom_entry = false;
    m_find_component_bond_entry = false;
    m_loop_category_map.clear();
    m_data_item_map.clear();
}

std::optional<std::string> CifFormat::GetFirstDataItemValue(std::string_view key) const
{
    auto iter{ m_data_item_map.find(std::string{key}) };
    if (iter == m_data_item_map.end()) return std::nullopt;
    if (iter->second.empty()) return std::nullopt;
    return iter->second.front();
}

void CifFormat::ParseMmCifDocument(std::istream & stream, const std::string & source_name)
{
    if (!stream)
    {
        Logger::Log(LogLevel::Error, "Cannot open the file: " + source_name);
        throw std::runtime_error("ParseMmCifDocument failed!");
    }

    std::vector<std::string> line_list;
    std::string line;
    while (std::getline(stream, line))
    {
        StringHelper::StripCarriageReturn(line);
        line_list.emplace_back(std::move(line));
    }

    ResetReadState();
    auto parsed_document{ ParseMmCifDocumentLines(line_list, source_name) };

    m_data_item_map = std::move(parsed_document.data_item_map);
    m_loop_category_map = std::move(parsed_document.loop_category_map);
}

std::unique_ptr<ModelObject> CifFormat::ReadModel(
    std::istream & stream, const std::string & source_name)
{
    ParseMmCifDocument(stream, source_name);
    LoadChemicalComponentBlock();
    LoadDatabaseBlock();
    LoadEntityBlock();
    LoadPdbxData();
    LoadAtomTypeBlock();
    LoadStructureConformationBlock();
    LoadStructureSheetBlock();

    if (m_find_component_bond_entry == false)
    {
        Logger::Log(LogLevel::Info,
            "No component bond entry found in the CIF file. Build default component bond entries.");
        BuildDefaultComponentBondEntry();
    }
    LoadAtomSiteBlock();
    ConstructBondList();
    LoadStructureConnectionBlock();
    LogHeaderSummary();
    return m_data_block->TakeModelObject();
}

void CifFormat::LoadChemicalComponentBlock()
{
    ParseLoopBlock("_chem_comp.",
        [this](const ColumnIndexMap & index_map,
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
    LoadChemicalComponentAtomBlock();
    LoadChemicalComponentBondBlock();
}

void CifFormat::LoadChemicalComponentAtomBlock()
{
    ParseLoopBlock("_chem_comp_atom.",
        [this](const ColumnIndexMap & index_map,
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

void CifFormat::LoadChemicalComponentBondBlock()
{
    ParseLoopBlock("_chem_comp_bond.",
        [this](const ColumnIndexMap & index_map,
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

void CifFormat::LoadAtomSiteBlock()
{
    std::unordered_map<AtomAltLocKey, AtomObject *, AtomAltLocKeyHash> altloc_primary_atom_map;
    size_t atom_site_row_count{ 0 };
    size_t atom_site_skip_count{ 0 };
    ParseLoopBlock("_atom_site.",
        [this, &altloc_primary_atom_map, &atom_site_row_count, &atom_site_skip_count](
            const ColumnIndexMap & index_map,
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
            auto sequence_id_token{
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
            if (IsMmCifMissingValue(sequence_id_token)) sequence_id_token = ".";
            auto sequence_id_value{
                ParseIntOrDefault(
                    sequence_id_token, -1, "_atom_site.(label/auth)_seq_id", context)
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

            AtomAltLocKey atom_altloc_key{
                model_number_id, chain_id, comp_id, sequence_id_token, atom_id
            };
            auto add_atom_as_primary{
                [&](std::unique_ptr<AtomObject> atom_object, const AtomAltLocKey & key) -> AtomObject *
                {
                    auto raw_atom_ptr{ atom_object.get() };
                    m_data_block->AddAtomObject(
                        model_number_id, std::move(atom_object), sequence_id_token);
                    altloc_primary_atom_map[key] = raw_atom_ptr;
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
                add_atom_as_primary(std::move(atom_object), atom_altloc_key);
                return;
            }

            // First alternate indicator is not "A" (e.g. only "B"): treat as primary.
            add_atom_as_primary(std::move(atom_object), atom_altloc_key);
        }
    );
    Logger::Log(LogLevel::Info,
        "LoadAtomSiteBlock parsed rows = " + std::to_string(atom_site_row_count)
        + ", skipped rows = " + std::to_string(atom_site_skip_count) + ".");
}

} // namespace rhbm_gem
