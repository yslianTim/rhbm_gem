#include "internal/io/file/CifFormat.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomicModelDataBlock.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include <optional>
#include <string_view>
#include <unordered_map>
#include <vector>

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

std::optional<int> TryParseInt(const std::string & value)
{
    if (IsMmCifMissingValue(value)) return std::nullopt;
    try
    {
        return std::stoi(value);
    }
    catch (const std::exception &)
    {
        return std::nullopt;
    }
}

} // namespace

void CifFormat::LoadDatabaseBlock()
{
    std::unordered_map<std::string, std::string> data_map;
    ParseLoopBlock("_database_2.",
        [&data_map](const ColumnIndexMap & index_map,
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

void CifFormat::LoadEntityBlock()
{
    ParseLoopBlock("_entity.",
        [this](const ColumnIndexMap & index_map,
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

    ParseLoopBlock("_struct_asym.",
        [this](const ColumnIndexMap & index_map,
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

    std::string entity_id{ GetFirstDataItemValue("_entity.id").value_or("") };
    std::string entity_type{ GetFirstDataItemValue("_entity.type").value_or("") };
    std::string molecules_size_raw{
        GetFirstDataItemValue("_entity.pdbx_number_of_molecules").value_or("") };
    std::string struct_asym_id{ GetFirstDataItemValue("_struct_asym.id").value_or("") };
    std::string struct_asym_entity_id{ GetFirstDataItemValue("_struct_asym.entity_id").value_or("") };

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

void CifFormat::LoadPdbxData()
{
    auto found_resolution{ false };
    auto found_resolution_method{ false };

    if (auto resolution_value{ GetFirstDataItemValue("_em_3d_reconstruction.resolution") })
    {
        if (!IsMmCifMissingValue(*resolution_value))
        {
            m_data_block->SetResolution(*resolution_value);
            found_resolution = true;
        }
    }
    if (auto resolution_method_value{
            GetFirstDataItemValue("_em_3d_reconstruction.resolution_method") })
    {
        if (!IsMmCifMissingValue(*resolution_method_value))
        {
            m_data_block->SetResolutionMethod(*resolution_method_value);
            found_resolution_method = true;
        }
    }

    ParseLoopBlock("_em_3d_reconstruction.",
        [this, &found_resolution, &found_resolution_method](
            const ColumnIndexMap & index_map,
            const std::vector<std::string> & token_list)
        {
            if (!found_resolution)
            {
                auto resolution_value{ GetTokenOptional(index_map, token_list, {"resolution"}) };
                if (resolution_value.has_value() && !IsMmCifMissingValue(*resolution_value))
                {
                    m_data_block->SetResolution(*resolution_value);
                    found_resolution = true;
                }
            }
            if (!found_resolution_method)
            {
                auto resolution_method_value{
                    GetTokenOptional(index_map, token_list, {"resolution_method"}) };
                if (resolution_method_value.has_value() &&
                    !IsMmCifMissingValue(*resolution_method_value))
                {
                    m_data_block->SetResolutionMethod(*resolution_method_value);
                    found_resolution_method = true;
                }
            }
        });

    if (!found_resolution || !found_resolution_method)
    {
        LoadXRayResolutionInfo();
    }
}

void CifFormat::LoadXRayResolutionInfo()
{
    auto found_resolution{ false };
    auto found_resolution_method{ false };

    if (auto resolution_value{ GetFirstDataItemValue("_refine.ls_d_res_high") })
    {
        if (!IsMmCifMissingValue(*resolution_value))
        {
            m_data_block->SetResolution(*resolution_value);
            found_resolution = true;
        }
    }
    if (auto method_value{ GetFirstDataItemValue("_refine.pdbx_refine_id") })
    {
        if (!IsMmCifMissingValue(*method_value))
        {
            m_data_block->SetResolutionMethod(*method_value);
            found_resolution_method = true;
        }
    }

    ParseLoopBlock("_refine.",
        [this, &found_resolution, &found_resolution_method](
            const ColumnIndexMap & index_map,
            const std::vector<std::string> & token_list)
        {
            if (!found_resolution)
            {
                auto resolution_value{ GetTokenOptional(index_map, token_list, {"ls_d_res_high"}) };
                if (resolution_value.has_value() && !IsMmCifMissingValue(*resolution_value))
                {
                    m_data_block->SetResolution(*resolution_value);
                    found_resolution = true;
                }
            }
            if (!found_resolution_method)
            {
                auto method_value{ GetTokenOptional(index_map, token_list, {"pdbx_refine_id"}) };
                if (method_value.has_value() && !IsMmCifMissingValue(*method_value))
                {
                    m_data_block->SetResolutionMethod(*method_value);
                    found_resolution_method = true;
                }
            }
        });
}

void CifFormat::LoadAtomTypeBlock()
{
    ParseLoopBlock("_atom_type.",
        [this](const ColumnIndexMap & index_map,
               const std::vector<std::string> & token_list)
        {
            auto element_type_string{ token_list[index_map.at("symbol")] };
            auto element{ ChemicalDataHelper::GetElementFromString(element_type_string) };
            m_data_block->AddElementType(element);
        }
    );
}

void CifFormat::LoadStructureConformationBlock()
{
    ParseLoopBlock("_struct_conf.",
        [this](const ColumnIndexMap & index_map,
               const std::vector<std::string> & token_list)
        {
            auto helix_id_opt{ GetTokenOptional(index_map, token_list, {"id"}) };
            auto conf_type_opt{ GetTokenOptional(index_map, token_list, {"conf_type_id"}) };
            auto chain_id_beg_opt{ GetTokenOptional(index_map, token_list, {"beg_label_asym_id"}) };
            auto residue_id_beg_opt{ GetTokenOptional(index_map, token_list, {"beg_label_seq_id"}) };
            auto chain_id_end_opt{ GetTokenOptional(index_map, token_list, {"end_label_asym_id"}) };
            auto residue_id_end_opt{ GetTokenOptional(index_map, token_list, {"end_label_seq_id"}) };
            if (!helix_id_opt.has_value() || !conf_type_opt.has_value() ||
                !chain_id_beg_opt.has_value() || !residue_id_beg_opt.has_value() ||
                !chain_id_end_opt.has_value() || !residue_id_end_opt.has_value())
            {
                return;
            }

            auto residue_id_beg{ TryParseInt(*residue_id_beg_opt) };
            auto residue_id_end{ TryParseInt(*residue_id_end_opt) };
            if (!residue_id_beg.has_value() || !residue_id_end.has_value())
            {
                Logger::Log(LogLevel::Warning,
                    "LoadStructureConformationBlock(): invalid residue range in _struct_conf. "
                    "row is skipped.");
                return;
            }

            HelixRange range;
            range.chain_id_beg = *chain_id_beg_opt;
            range.seq_id_beg = *residue_id_beg;
            range.chain_id_end = *chain_id_end_opt;
            range.seq_id_end = *residue_id_end;
            range.conf_type = *conf_type_opt;

            m_data_block->AddHelixRange(*helix_id_opt, range);
        }
    );
}

void CifFormat::LoadStructureConnectionBlock()
{
    ParseLoopBlock("_struct_conn.",
        [this](const ColumnIndexMap & index_map,
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

void CifFormat::LoadStructureSheetBlock()
{
    ParseLoopBlock("_struct_sheet.",
        [this](const ColumnIndexMap & index_map,
               const std::vector<std::string> & token_list)
        {
            auto sheet_id_opt{ GetTokenOptional(index_map, token_list, {"id"}) };
            auto strands_size_opt{ GetTokenOptional(index_map, token_list, {"number_strands"}) };
            if (!sheet_id_opt.has_value() || !strands_size_opt.has_value()) return;
            auto strands_size{
                ParseIntOrDefault(*strands_size_opt, 0, "_struct_sheet.number_strands",
                    "LoadStructureSheetBlock()")
            };
            if (strands_size <= 0) return;
            auto sheet_id{ *sheet_id_opt };
            m_data_block->AddSheetStrands(sheet_id, strands_size);
        }
    );

    ParseLoopBlock("_struct_sheet_range.",
        [this](const ColumnIndexMap & index_map,
               const std::vector<std::string> & token_list)
        {
            auto sheet_id_opt{ GetTokenOptional(index_map, token_list, {"sheet_id"}) };
            auto range_id_opt{ GetTokenOptional(index_map, token_list, {"id"}) };
            auto chain_id_beg_opt{ GetTokenOptional(index_map, token_list, {"beg_label_asym_id"}) };
            auto residue_id_beg_opt{ GetTokenOptional(index_map, token_list, {"beg_label_seq_id"}) };
            auto chain_id_end_opt{ GetTokenOptional(index_map, token_list, {"end_label_asym_id"}) };
            auto residue_id_end_opt{ GetTokenOptional(index_map, token_list, {"end_label_seq_id"}) };
            if (!sheet_id_opt.has_value() || !range_id_opt.has_value() ||
                !chain_id_beg_opt.has_value() || !residue_id_beg_opt.has_value() ||
                !chain_id_end_opt.has_value() || !residue_id_end_opt.has_value())
            {
                return;
            }

            auto residue_id_beg{ TryParseInt(*residue_id_beg_opt) };
            auto residue_id_end{ TryParseInt(*residue_id_end_opt) };
            if (!residue_id_beg.has_value() || !residue_id_end.has_value())
            {
                Logger::Log(LogLevel::Warning,
                    "LoadStructureSheetBlock(): invalid residue range in _struct_sheet_range. "
                    "row is skipped.");
                return;
            }

            auto sheet_id{ *sheet_id_opt };
            auto range_id{ *range_id_opt };
            auto composite_sheet_id{ sheet_id + range_id };
            SheetRange range;
            range.chain_id_beg = *chain_id_beg_opt;
            range.seq_id_beg = *residue_id_beg;
            range.chain_id_end = *chain_id_end_opt;
            range.seq_id_end = *residue_id_end;
            m_data_block->AddSheetRange(composite_sheet_id, range);
        }
    );
}

} // namespace rhbm_gem
