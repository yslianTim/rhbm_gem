#include "ModelObjectDAOv2.hpp"

#include "AtomClassifier.hpp"
#include "AtomObject.hpp"
#include "BondClassifier.hpp"
#include "BondObject.hpp"
#include "ChemicalComponentEntry.hpp"
#include "ChemicalDataHelper.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"
#include "BondKeySystem.hpp"
#include "DataObjectBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "LocalPotentialEntry.hpp"
#include "ModelObject.hpp"
#include "SQLiteWrapper.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace
{
    constexpr std::string_view CREATE_MODEL_OBJECT_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_object (
            key_tag TEXT PRIMARY KEY,
            atom_size INTEGER,
            pdb_id TEXT,
            emd_id TEXT,
            map_resolution DOUBLE,
            resolution_method TEXT
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_CHAIN_MAP_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_chain_map (
            key_tag TEXT,
            entity_id TEXT,
            chain_ordinal INTEGER,
            chain_id TEXT,
            PRIMARY KEY (key_tag, entity_id, chain_ordinal)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_COMPONENT_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_component (
            key_tag TEXT,
            component_key INTEGER,
            id TEXT,
            name TEXT,
            type TEXT,
            formula TEXT,
            molecular_weight DOUBLE,
            is_standard_monomer INTEGER,
            PRIMARY KEY (key_tag, component_key)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_COMPONENT_ATOM_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_component_atom (
            key_tag TEXT,
            component_key INTEGER,
            atom_key INTEGER,
            atom_id TEXT,
            element_type INTEGER,
            aromatic_atom_flag INTEGER,
            stereo_config INTEGER,
            PRIMARY KEY (key_tag, component_key, atom_key)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_COMPONENT_BOND_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_component_bond (
            key_tag TEXT,
            component_key INTEGER,
            bond_key INTEGER,
            bond_id TEXT,
            bond_type INTEGER,
            bond_order INTEGER,
            aromatic_atom_flag INTEGER,
            stereo_config INTEGER,
            PRIMARY KEY (key_tag, component_key, bond_key)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_ATOM_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_atom (
            key_tag TEXT,
            serial_id INTEGER,
            sequence_id INTEGER,
            component_id TEXT,
            atom_id TEXT,
            chain_id TEXT,
            indicator TEXT,
            occupancy DOUBLE,
            temperature DOUBLE,
            element INTEGER,
            structure INTEGER,
            is_special_atom INTEGER,
            position_x DOUBLE,
            position_y DOUBLE,
            position_z DOUBLE,
            component_key INTEGER,
            atom_key INTEGER,
            PRIMARY KEY (key_tag, serial_id)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_BOND_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_bond (
            key_tag TEXT,
            atom_serial_id_1 INTEGER,
            atom_serial_id_2 INTEGER,
            bond_key INTEGER,
            bond_type INTEGER,
            bond_order INTEGER,
            is_special_bond INTEGER,
            PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_ATOM_LOCAL_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_atom_local_potential (
            key_tag TEXT,
            serial_id INTEGER,
            sampling_size INTEGER,
            distance_and_map_value_list BLOB,
            amplitude_estimate_ols DOUBLE,
            width_estimate_ols DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE,
            alpha_r DOUBLE,
            PRIMARY KEY (key_tag, serial_id)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_BOND_LOCAL_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_bond_local_potential (
            key_tag TEXT,
            atom_serial_id_1 INTEGER,
            atom_serial_id_2 INTEGER,
            sampling_size INTEGER,
            distance_and_map_value_list BLOB,
            amplitude_estimate_ols DOUBLE,
            width_estimate_ols DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE,
            alpha_r DOUBLE,
            PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_ATOM_POSTERIOR_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_atom_posterior (
            key_tag TEXT,
            class_key TEXT,
            serial_id INTEGER,
            amplitude_estimate_posterior DOUBLE,
            width_estimate_posterior DOUBLE,
            amplitude_variance_posterior DOUBLE,
            width_variance_posterior DOUBLE,
            outlier_tag INTEGER,
            statistical_distance DOUBLE,
            PRIMARY KEY (key_tag, class_key, serial_id)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_BOND_POSTERIOR_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_bond_posterior (
            key_tag TEXT,
            class_key TEXT,
            atom_serial_id_1 INTEGER,
            atom_serial_id_2 INTEGER,
            amplitude_estimate_posterior DOUBLE,
            width_estimate_posterior DOUBLE,
            amplitude_variance_posterior DOUBLE,
            width_variance_posterior DOUBLE,
            outlier_tag INTEGER,
            statistical_distance DOUBLE,
            PRIMARY KEY (key_tag, class_key, atom_serial_id_1, atom_serial_id_2)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_ATOM_GROUP_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_atom_group_potential (
            key_tag TEXT,
            class_key TEXT,
            group_key INTEGER,
            member_size INTEGER,
            amplitude_estimate_mean DOUBLE,
            width_estimate_mean DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE,
            amplitude_estimate_prior DOUBLE,
            width_estimate_prior DOUBLE,
            amplitude_variance_prior DOUBLE,
            width_variance_prior DOUBLE,
            alpha_g DOUBLE,
            PRIMARY KEY (key_tag, class_key, group_key)
        )
    )sql";

    constexpr std::string_view CREATE_MODEL_BOND_GROUP_TABLE_SQL = R"sql(
        CREATE TABLE IF NOT EXISTS model_bond_group_potential (
            key_tag TEXT,
            class_key TEXT,
            group_key INTEGER,
            member_size INTEGER,
            amplitude_estimate_mean DOUBLE,
            width_estimate_mean DOUBLE,
            amplitude_estimate_mdpde DOUBLE,
            width_estimate_mdpde DOUBLE,
            amplitude_estimate_prior DOUBLE,
            width_estimate_prior DOUBLE,
            amplitude_variance_prior DOUBLE,
            width_variance_prior DOUBLE,
            alpha_g DOUBLE,
            PRIMARY KEY (key_tag, class_key, group_key)
        )
    )sql";

    constexpr std::string_view UPSERT_MODEL_OBJECT_SQL = R"sql(
        INSERT INTO model_object (
            key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
        ) VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(key_tag) DO UPDATE SET
            atom_size = excluded.atom_size,
            pdb_id = excluded.pdb_id,
            emd_id = excluded.emd_id,
            map_resolution = excluded.map_resolution,
            resolution_method = excluded.resolution_method
    )sql";

    constexpr std::string_view INSERT_MODEL_CHAIN_MAP_SQL = R"sql(
        INSERT OR REPLACE INTO model_chain_map (
            key_tag, entity_id, chain_ordinal, chain_id
        ) VALUES (?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_COMPONENT_SQL = R"sql(
        INSERT OR REPLACE INTO model_component (
            key_tag, component_key, id, name, type, formula,
            molecular_weight, is_standard_monomer
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_COMPONENT_ATOM_SQL = R"sql(
        INSERT OR REPLACE INTO model_component_atom (
            key_tag, component_key, atom_key, atom_id,
            element_type, aromatic_atom_flag, stereo_config
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_COMPONENT_BOND_SQL = R"sql(
        INSERT OR REPLACE INTO model_component_bond (
            key_tag, component_key, bond_key, bond_id,
            bond_type, bond_order, aromatic_atom_flag, stereo_config
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_ATOM_SQL = R"sql(
        INSERT OR REPLACE INTO model_atom (
            key_tag, serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
            occupancy, temperature, element, structure, is_special_atom,
            position_x, position_y, position_z, component_key, atom_key
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_BOND_SQL = R"sql(
        INSERT OR REPLACE INTO model_bond (
            key_tag, atom_serial_id_1, atom_serial_id_2,
            bond_key, bond_type, bond_order, is_special_bond
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_ATOM_LOCAL_SQL = R"sql(
        INSERT OR REPLACE INTO model_atom_local_potential (
            key_tag, serial_id, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_BOND_LOCAL_SQL = R"sql(
        INSERT OR REPLACE INTO model_bond_local_potential (
            key_tag, atom_serial_id_1, atom_serial_id_2, sampling_size,
            distance_and_map_value_list, amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_ATOM_POSTERIOR_SQL = R"sql(
        INSERT OR REPLACE INTO model_atom_posterior (
            key_tag, class_key, serial_id, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_BOND_POSTERIOR_SQL = R"sql(
        INSERT OR REPLACE INTO model_bond_posterior (
            key_tag, class_key, atom_serial_id_1, atom_serial_id_2,
            amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_ATOM_GROUP_SQL = R"sql(
        INSERT OR REPLACE INTO model_atom_group_potential (
            key_tag, class_key, group_key, member_size,
            amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior, alpha_g
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view INSERT_MODEL_BOND_GROUP_SQL = R"sql(
        INSERT OR REPLACE INTO model_bond_group_potential (
            key_tag, class_key, group_key, member_size,
            amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior, alpha_g
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    )sql";

    constexpr std::string_view DELETE_ROWS_FOR_KEY_SQL_PREFIX = "DELETE FROM ";
    constexpr std::string_view DELETE_ROWS_FOR_KEY_SQL_SUFFIX = " WHERE key_tag = ?;";

    constexpr std::string_view SELECT_MODEL_OBJECT_SQL = R"sql(
        SELECT key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method
        FROM model_object WHERE key_tag = ? LIMIT 1;
    )sql";

    constexpr std::string_view SELECT_MODEL_CHAIN_MAP_SQL = R"sql(
        SELECT entity_id, chain_ordinal, chain_id
        FROM model_chain_map WHERE key_tag = ?
        ORDER BY entity_id, chain_ordinal;
    )sql";

    constexpr std::string_view SELECT_MODEL_COMPONENT_SQL = R"sql(
        SELECT component_key, id, name, type, formula, molecular_weight, is_standard_monomer
        FROM model_component WHERE key_tag = ?
        ORDER BY component_key;
    )sql";

    constexpr std::string_view SELECT_MODEL_COMPONENT_ATOM_SQL = R"sql(
        SELECT component_key, atom_key, atom_id, element_type, aromatic_atom_flag, stereo_config
        FROM model_component_atom WHERE key_tag = ?
        ORDER BY component_key, atom_key;
    )sql";

    constexpr std::string_view SELECT_MODEL_COMPONENT_BOND_SQL = R"sql(
        SELECT component_key, bond_key, bond_id, bond_type, bond_order, aromatic_atom_flag, stereo_config
        FROM model_component_bond WHERE key_tag = ?
        ORDER BY component_key, bond_key;
    )sql";

    constexpr std::string_view SELECT_MODEL_ATOM_SQL = R"sql(
        SELECT
            serial_id, sequence_id, component_id, atom_id, chain_id, indicator,
            occupancy, temperature, element, structure, is_special_atom,
            position_x, position_y, position_z, component_key, atom_key
        FROM model_atom WHERE key_tag = ?
        ORDER BY serial_id;
    )sql";

    constexpr std::string_view SELECT_MODEL_BOND_SQL = R"sql(
        SELECT atom_serial_id_1, atom_serial_id_2, bond_key, bond_type, bond_order, is_special_bond
        FROM model_bond WHERE key_tag = ?
        ORDER BY atom_serial_id_1, atom_serial_id_2;
    )sql";

    constexpr std::string_view SELECT_MODEL_ATOM_LOCAL_SQL = R"sql(
        SELECT
            serial_id, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
        FROM model_atom_local_potential WHERE key_tag = ?;
    )sql";

    constexpr std::string_view SELECT_MODEL_BOND_LOCAL_SQL = R"sql(
        SELECT
            atom_serial_id_1, atom_serial_id_2, sampling_size, distance_and_map_value_list,
            amplitude_estimate_ols, width_estimate_ols,
            amplitude_estimate_mdpde, width_estimate_mdpde, alpha_r
        FROM model_bond_local_potential WHERE key_tag = ?;
    )sql";

    constexpr std::string_view SELECT_MODEL_ATOM_POSTERIOR_SQL = R"sql(
        SELECT
            serial_id, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
        FROM model_atom_posterior WHERE key_tag = ? AND class_key = ?;
    )sql";

    constexpr std::string_view SELECT_MODEL_BOND_POSTERIOR_SQL = R"sql(
        SELECT
            atom_serial_id_1, atom_serial_id_2, amplitude_estimate_posterior, width_estimate_posterior,
            amplitude_variance_posterior, width_variance_posterior, outlier_tag, statistical_distance
        FROM model_bond_posterior WHERE key_tag = ? AND class_key = ?;
    )sql";

    constexpr std::string_view SELECT_MODEL_ATOM_GROUP_SQL = R"sql(
        SELECT
            group_key, member_size, amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior, alpha_g
        FROM model_atom_group_potential WHERE key_tag = ? AND class_key = ?;
    )sql";

    constexpr std::string_view SELECT_MODEL_BOND_GROUP_SQL = R"sql(
        SELECT
            group_key, member_size, amplitude_estimate_mean, width_estimate_mean,
            amplitude_estimate_mdpde, width_estimate_mdpde,
            amplitude_estimate_prior, width_estimate_prior,
            amplitude_variance_prior, width_variance_prior, alpha_g
        FROM model_bond_group_potential WHERE key_tag = ? AND class_key = ?;
    )sql";
}

namespace rhbm_gem {

ModelObjectDAOv2::ModelObjectDAOv2(SQLiteWrapper * db_manager) :
    m_database{ db_manager }
{
}

ModelObjectDAOv2::~ModelObjectDAOv2()
{
}

void ModelObjectDAOv2::EnsureSchema(SQLiteWrapper & database)
{
    database.Execute(std::string(CREATE_MODEL_OBJECT_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_CHAIN_MAP_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_COMPONENT_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_COMPONENT_ATOM_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_COMPONENT_BOND_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_ATOM_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_BOND_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_ATOM_LOCAL_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_BOND_LOCAL_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_ATOM_POSTERIOR_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_BOND_POSTERIOR_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_ATOM_GROUP_TABLE_SQL));
    database.Execute(std::string(CREATE_MODEL_BOND_GROUP_TABLE_SQL));
}

void ModelObjectDAOv2::Save(const DataObjectBase * obj, const std::string & key_tag)
{
    auto model_obj{ dynamic_cast<const ModelObject *>(obj) };
    if (!model_obj)
    {
        throw std::runtime_error(
            "ModelObjectDAOv2::Save() failed: object is not a ModelObject instance.");
    }

    EnsureSchema(*m_database);

    DeleteRowsForKey("model_chain_map", key_tag);
    DeleteRowsForKey("model_component", key_tag);
    DeleteRowsForKey("model_component_atom", key_tag);
    DeleteRowsForKey("model_component_bond", key_tag);
    DeleteRowsForKey("model_atom", key_tag);
    DeleteRowsForKey("model_bond", key_tag);
    DeleteRowsForKey("model_atom_local_potential", key_tag);
    DeleteRowsForKey("model_bond_local_potential", key_tag);
    DeleteRowsForKey("model_atom_posterior", key_tag);
    DeleteRowsForKey("model_bond_posterior", key_tag);
    DeleteRowsForKey("model_atom_group_potential", key_tag);
    DeleteRowsForKey("model_bond_group_potential", key_tag);

    SaveModelObjectRow(model_obj, key_tag);
    SaveChainMap(model_obj, key_tag);
    SaveChemicalComponentEntryList(model_obj, key_tag);
    SaveComponentAtomEntryList(model_obj, key_tag);
    SaveComponentBondEntryList(model_obj, key_tag);
    SaveAtomObjectList(model_obj, key_tag);
    SaveBondObjectList(model_obj, key_tag);
    SaveAtomLocalPotentialEntryList(model_obj, key_tag);
    SaveBondLocalPotentialEntryList(model_obj, key_tag);

    for (const auto & [class_key, group_entry] : model_obj->GetAtomGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveAtomLocalPotentialEntrySubList(model_obj, key_tag, class_key);
        SaveAtomGroupPotentialEntryList(group_entry.get(), key_tag, class_key);
    }
    for (const auto & [class_key, group_entry] : model_obj->GetBondGroupPotentialEntryMap())
    {
        if (group_entry == nullptr || group_entry->GetGroupKeySet().empty())
        {
            continue;
        }
        SaveBondLocalPotentialEntrySubList(model_obj, key_tag, class_key);
        SaveBondGroupPotentialEntryList(group_entry.get(), key_tag, class_key);
    }
}

std::unique_ptr<DataObjectBase> ModelObjectDAOv2::Load(const std::string & key_tag)
{
    auto model_object{ std::make_unique<ModelObject>() };
    LoadChemicalComponentEntryList(model_object.get(), key_tag);
    LoadComponentAtomEntryList(model_object.get(), key_tag);
    LoadComponentBondEntryList(model_object.get(), key_tag);
    model_object->SetAtomList(LoadAtomObjectList(key_tag));
    model_object->SetBondList(LoadBondObjectList(key_tag, model_object.get()));
    LoadModelObjectRow(model_object.get(), key_tag);
    LoadChainMap(model_object.get(), key_tag);

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddAtomGroupPotentialEntry(class_key, group_entry);
        LoadAtomGroupPotentialEntryList(model_object.get(), key_tag, class_key);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        auto class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_entry{ std::make_unique<GroupPotentialEntry>() };
        model_object->AddBondGroupPotentialEntry(class_key, group_entry);
        LoadBondGroupPotentialEntryList(model_object.get(), key_tag, class_key);
    }

    return model_object;
}

void ModelObjectDAOv2::DeleteRowsForKey(const std::string & table_name, const std::string & key_tag) const
{
    m_database->Prepare(
        std::string(DELETE_ROWS_FOR_KEY_SQL_PREFIX) + table_name + std::string(DELETE_ROWS_FOR_KEY_SQL_SUFFIX));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->StepOnce();
}

void ModelObjectDAOv2::SaveModelObjectRow(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(UPSERT_MODEL_OBJECT_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<int>(2, static_cast<int>(model_obj->GetNumberOfAtom()));
    m_database->Bind<std::string>(3, model_obj->GetPdbID());
    m_database->Bind<std::string>(4, model_obj->GetEmdID());
    m_database->Bind<double>(5, model_obj->GetResolution());
    m_database->Bind<std::string>(6, model_obj->GetResolutionMethod());
    m_database->StepOnce();
}

void ModelObjectDAOv2::SaveChainMap(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_CHAIN_MAP_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & [entity_id, chain_id_list] : model_obj->GetChainIDListMap())
    {
        for (size_t ordinal = 0; ordinal < chain_id_list.size(); ordinal++)
        {
            m_database->Bind<std::string>(1, key_tag);
            m_database->Bind<std::string>(2, entity_id);
            m_database->Bind<int>(3, static_cast<int>(ordinal));
            m_database->Bind<std::string>(4, chain_id_list.at(ordinal));
            m_database->StepOnce();
            m_database->Reset();
        }
    }
}

void ModelObjectDAOv2::SaveChemicalComponentEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_COMPONENT_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & [component_key, component_entry] : model_obj->GetChemicalComponentEntryMap())
    {
        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<ComponentKey>(2, component_key);
        m_database->Bind<std::string>(3, component_entry->GetComponentId());
        m_database->Bind<std::string>(4, component_entry->GetComponentName());
        m_database->Bind<std::string>(5, component_entry->GetComponentType());
        m_database->Bind<std::string>(6, component_entry->GetComponentFormula());
        m_database->Bind<double>(7, static_cast<double>(component_entry->GetComponentMolecularWeight()));
        m_database->Bind<int>(8, static_cast<int>(component_entry->IsStandardMonomer()));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveComponentAtomEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_COMPONENT_ATOM_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & [component_key, component_entry] : model_obj->GetChemicalComponentEntryMap())
    {
        for (const auto & [atom_key, atom_entry] : component_entry->GetComponentAtomEntryMap())
        {
            m_database->Bind<std::string>(1, key_tag);
            m_database->Bind<ComponentKey>(2, component_key);
            m_database->Bind<AtomKey>(3, atom_key);
            m_database->Bind<std::string>(4, atom_entry.atom_id);
            m_database->Bind<int>(5, static_cast<int>(atom_entry.element_type));
            m_database->Bind<int>(6, static_cast<int>(atom_entry.aromatic_atom_flag));
            m_database->Bind<int>(7, static_cast<int>(atom_entry.stereo_config));
            m_database->StepOnce();
            m_database->Reset();
        }
    }
}

void ModelObjectDAOv2::SaveComponentBondEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_COMPONENT_BOND_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & [component_key, component_entry] : model_obj->GetChemicalComponentEntryMap())
    {
        for (const auto & [bond_key, bond_entry] : component_entry->GetComponentBondEntryMap())
        {
            m_database->Bind<std::string>(1, key_tag);
            m_database->Bind<ComponentKey>(2, component_key);
            m_database->Bind<BondKey>(3, bond_key);
            m_database->Bind<std::string>(4, bond_entry.bond_id);
            m_database->Bind<int>(5, static_cast<int>(bond_entry.bond_type));
            m_database->Bind<int>(6, static_cast<int>(bond_entry.bond_order));
            m_database->Bind<int>(7, static_cast<int>(bond_entry.aromatic_atom_flag));
            m_database->Bind<int>(8, static_cast<int>(bond_entry.stereo_config));
            m_database->StepOnce();
            m_database->Reset();
        }
    }
}

void ModelObjectDAOv2::SaveAtomObjectList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_ATOM_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & atom_object : model_obj->GetAtomList())
    {
        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<int>(2, atom_object->GetSerialID());
        m_database->Bind<int>(3, atom_object->GetSequenceID());
        m_database->Bind<std::string>(4, atom_object->GetComponentID());
        m_database->Bind<std::string>(5, atom_object->GetAtomID());
        m_database->Bind<std::string>(6, atom_object->GetChainID());
        m_database->Bind<std::string>(7, atom_object->GetIndicator());
        m_database->Bind<double>(8, static_cast<double>(atom_object->GetOccupancy()));
        m_database->Bind<double>(9, static_cast<double>(atom_object->GetTemperature()));
        m_database->Bind<int>(10, static_cast<int>(atom_object->GetElement()));
        m_database->Bind<int>(11, static_cast<int>(atom_object->GetStructure()));
        m_database->Bind<int>(12, static_cast<int>(atom_object->GetSpecialAtomFlag()));
        m_database->Bind<double>(13, static_cast<double>(atom_object->GetPosition().at(0)));
        m_database->Bind<double>(14, static_cast<double>(atom_object->GetPosition().at(1)));
        m_database->Bind<double>(15, static_cast<double>(atom_object->GetPosition().at(2)));
        m_database->Bind<int>(16, static_cast<int>(atom_object->GetComponentKey()));
        m_database->Bind<int>(17, static_cast<int>(atom_object->GetAtomKey()));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveBondObjectList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_BOND_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & bond_object : model_obj->GetBondList())
    {
        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<int>(2, bond_object->GetAtomSerialID1());
        m_database->Bind<int>(3, bond_object->GetAtomSerialID2());
        m_database->Bind<int>(4, bond_object->GetBondKey());
        m_database->Bind<int>(5, static_cast<int>(bond_object->GetBondType()));
        m_database->Bind<int>(6, static_cast<int>(bond_object->GetBondOrder()));
        m_database->Bind<int>(7, static_cast<int>(bond_object->GetSpecialBondFlag()));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveAtomLocalPotentialEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_ATOM_LOCAL_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & atom_object : model_obj->GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<int>(2, atom_object->GetSerialID());
        m_database->Bind<int>(3, entry->GetDistanceAndMapValueListSize());
        m_database->Bind<std::vector<std::tuple<float, float>>>(4, entry->GetDistanceAndMapValueList());
        m_database->Bind<double>(5, entry->GetAmplitudeEstimateOLS());
        m_database->Bind<double>(6, entry->GetWidthEstimateOLS());
        m_database->Bind<double>(7, entry->GetAmplitudeEstimateMDPDE());
        m_database->Bind<double>(8, entry->GetWidthEstimateMDPDE());
        m_database->Bind<double>(9, entry->GetAlphaR());
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveBondLocalPotentialEntryList(const ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(INSERT_MODEL_BOND_LOCAL_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & bond_object : model_obj->GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<int>(2, bond_object->GetAtomSerialID1());
        m_database->Bind<int>(3, bond_object->GetAtomSerialID2());
        m_database->Bind<int>(4, entry->GetDistanceAndMapValueListSize());
        m_database->Bind<std::vector<std::tuple<float, float>>>(5, entry->GetDistanceAndMapValueList());
        m_database->Bind<double>(6, entry->GetAmplitudeEstimateOLS());
        m_database->Bind<double>(7, entry->GetWidthEstimateOLS());
        m_database->Bind<double>(8, entry->GetAmplitudeEstimateMDPDE());
        m_database->Bind<double>(9, entry->GetWidthEstimateMDPDE());
        m_database->Bind<double>(10, entry->GetAlphaR());
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveAtomLocalPotentialEntrySubList(
    const ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    m_database->Prepare(std::string(INSERT_MODEL_ATOM_POSTERIOR_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & atom_object : model_obj->GetAtomList())
    {
        auto entry{ atom_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<std::string>(2, class_key);
        m_database->Bind<int>(3, atom_object->GetSerialID());
        m_database->Bind<double>(4, entry->GetAmplitudeEstimatePosterior(class_key));
        m_database->Bind<double>(5, entry->GetWidthEstimatePosterior(class_key));
        m_database->Bind<double>(6, entry->GetAmplitudeVariancePosterior(class_key));
        m_database->Bind<double>(7, entry->GetWidthVariancePosterior(class_key));
        m_database->Bind<int>(8, static_cast<int>(entry->GetOutlierTag(class_key)));
        m_database->Bind<double>(9, entry->GetStatisticalDistance(class_key));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveBondLocalPotentialEntrySubList(
    const ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    m_database->Prepare(std::string(INSERT_MODEL_BOND_POSTERIOR_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & bond_object : model_obj->GetBondList())
    {
        auto entry{ bond_object->GetLocalPotentialEntry() };
        if (entry == nullptr) continue;

        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<std::string>(2, class_key);
        m_database->Bind<int>(3, bond_object->GetAtomSerialID1());
        m_database->Bind<int>(4, bond_object->GetAtomSerialID2());
        m_database->Bind<double>(5, entry->GetAmplitudeEstimatePosterior(class_key));
        m_database->Bind<double>(6, entry->GetWidthEstimatePosterior(class_key));
        m_database->Bind<double>(7, entry->GetAmplitudeVariancePosterior(class_key));
        m_database->Bind<double>(8, entry->GetWidthVariancePosterior(class_key));
        m_database->Bind<int>(9, static_cast<int>(entry->GetOutlierTag(class_key)));
        m_database->Bind<double>(10, entry->GetStatisticalDistance(class_key));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveAtomGroupPotentialEntryList(
    const GroupPotentialEntry * group_entry, const std::string & key_tag, const std::string & class_key)
{
    m_database->Prepare(std::string(INSERT_MODEL_ATOM_GROUP_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & group_key : group_entry->GetGroupKeySet())
    {
        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<std::string>(2, class_key);
        m_database->Bind<GroupKey>(3, group_key);
        m_database->Bind<int>(4, group_entry->GetAtomObjectPtrListSize(group_key));
        m_database->Bind<double>(5, std::get<0>(group_entry->GetGausEstimateMean(group_key)));
        m_database->Bind<double>(6, std::get<1>(group_entry->GetGausEstimateMean(group_key)));
        m_database->Bind<double>(7, std::get<0>(group_entry->GetGausEstimateMDPDE(group_key)));
        m_database->Bind<double>(8, std::get<1>(group_entry->GetGausEstimateMDPDE(group_key)));
        m_database->Bind<double>(9, std::get<0>(group_entry->GetGausEstimatePrior(group_key)));
        m_database->Bind<double>(10, std::get<1>(group_entry->GetGausEstimatePrior(group_key)));
        m_database->Bind<double>(11, std::get<0>(group_entry->GetGausVariancePrior(group_key)));
        m_database->Bind<double>(12, std::get<1>(group_entry->GetGausVariancePrior(group_key)));
        m_database->Bind<double>(13, group_entry->GetAlphaG(group_key));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::SaveBondGroupPotentialEntryList(
    const GroupPotentialEntry * group_entry, const std::string & key_tag, const std::string & class_key)
{
    m_database->Prepare(std::string(INSERT_MODEL_BOND_GROUP_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    for (const auto & group_key : group_entry->GetGroupKeySet())
    {
        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<std::string>(2, class_key);
        m_database->Bind<GroupKey>(3, group_key);
        m_database->Bind<int>(4, group_entry->GetBondObjectPtrListSize(group_key));
        m_database->Bind<double>(5, std::get<0>(group_entry->GetGausEstimateMean(group_key)));
        m_database->Bind<double>(6, std::get<1>(group_entry->GetGausEstimateMean(group_key)));
        m_database->Bind<double>(7, std::get<0>(group_entry->GetGausEstimateMDPDE(group_key)));
        m_database->Bind<double>(8, std::get<1>(group_entry->GetGausEstimateMDPDE(group_key)));
        m_database->Bind<double>(9, std::get<0>(group_entry->GetGausEstimatePrior(group_key)));
        m_database->Bind<double>(10, std::get<1>(group_entry->GetGausEstimatePrior(group_key)));
        m_database->Bind<double>(11, std::get<0>(group_entry->GetGausVariancePrior(group_key)));
        m_database->Bind<double>(12, std::get<1>(group_entry->GetGausVariancePrior(group_key)));
        m_database->Bind<double>(13, group_entry->GetAlphaG(group_key));
        m_database->StepOnce();
        m_database->Reset();
    }
}

void ModelObjectDAOv2::LoadModelObjectRow(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(SELECT_MODEL_OBJECT_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    const auto rc{ m_database->StepNext() };
    if (rc == SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    const auto atom_size{ m_database->GetColumn<int>(1) };
    model_obj->SetKeyTag(m_database->GetColumn<std::string>(0));
    model_obj->SetPdbID(m_database->GetColumn<std::string>(2));
    model_obj->SetEmdID(m_database->GetColumn<std::string>(3));
    model_obj->SetResolution(m_database->GetColumn<double>(4));
    model_obj->SetResolutionMethod(m_database->GetColumn<std::string>(5));
    if (atom_size != static_cast<int>(model_obj->GetNumberOfAtom()))
    {
        throw std::runtime_error(
            "The number of atoms in the model object does not match the database record.");
    }
}

void ModelObjectDAOv2::LoadChainMap(ModelObject * model_obj, const std::string & key_tag)
{
    std::unordered_map<std::string, std::vector<std::string>> chain_map;
    m_database->Prepare(std::string(SELECT_MODEL_CHAIN_MAP_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }
        auto & chain_id_list{ chain_map[m_database->GetColumn<std::string>(0)] };
        const auto ordinal{ m_database->GetColumn<int>(1) };
        if (static_cast<int>(chain_id_list.size()) <= ordinal)
        {
            chain_id_list.resize(static_cast<size_t>(ordinal + 1));
        }
        chain_id_list.at(static_cast<size_t>(ordinal)) = m_database->GetColumn<std::string>(2);
    }
    model_obj->SetChainIDListMap(chain_map);
}

void ModelObjectDAOv2::LoadChemicalComponentEntryList(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(SELECT_MODEL_COMPONENT_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto component_entry{ std::make_unique<ChemicalComponentEntry>() };
        const auto component_key{ m_database->GetColumn<ComponentKey>(0) };
        const auto component_id{ m_database->GetColumn<std::string>(1) };
        component_entry->SetComponentId(component_id);
        component_entry->SetComponentName(m_database->GetColumn<std::string>(2));
        component_entry->SetComponentType(m_database->GetColumn<std::string>(3));
        component_entry->SetComponentFormula(m_database->GetColumn<std::string>(4));
        component_entry->SetComponentMolecularWeight(static_cast<float>(m_database->GetColumn<double>(5)));
        component_entry->SetStandardMonomerFlag(static_cast<bool>(m_database->GetColumn<int>(6)));
        model_obj->AddChemicalComponentEntry(component_key, std::move(component_entry));
        model_obj->GetComponentKeySystemPtr()->RegisterComponent(component_id, component_key);
    }
}

void ModelObjectDAOv2::LoadComponentAtomEntryList(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(SELECT_MODEL_COMPONENT_ATOM_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto component_key{ m_database->GetColumn<ComponentKey>(0) };
        if (model_obj->GetChemicalComponentEntryMap().find(component_key)
            == model_obj->GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        auto atom_entry = ComponentAtomEntry{
            m_database->GetColumn<std::string>(2),
            static_cast<Element>(m_database->GetColumn<int>(3)),
            static_cast<bool>(m_database->GetColumn<int>(4)),
            static_cast<StereoChemistry>(m_database->GetColumn<int>(5))
        };
        const auto atom_key{ m_database->GetColumn<AtomKey>(1) };
        model_obj->GetChemicalComponentEntryMap().at(component_key)->AddComponentAtomEntry(atom_key, atom_entry);
        model_obj->GetAtomKeySystemPtr()->RegisterAtom(atom_entry.atom_id, atom_key);
    }
}

void ModelObjectDAOv2::LoadComponentBondEntryList(ModelObject * model_obj, const std::string & key_tag)
{
    m_database->Prepare(std::string(SELECT_MODEL_COMPONENT_BOND_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto component_key{ m_database->GetColumn<ComponentKey>(0) };
        if (model_obj->GetChemicalComponentEntryMap().find(component_key)
            == model_obj->GetChemicalComponentEntryMap().end())
        {
            continue;
        }
        ComponentBondEntry bond_entry;
        const auto bond_key{ m_database->GetColumn<BondKey>(1) };
        bond_entry.bond_id = m_database->GetColumn<std::string>(2);
        bond_entry.bond_type = static_cast<BondType>(m_database->GetColumn<int>(3));
        bond_entry.bond_order = static_cast<BondOrder>(m_database->GetColumn<int>(4));
        bond_entry.aromatic_atom_flag = static_cast<bool>(m_database->GetColumn<int>(5));
        bond_entry.stereo_config = static_cast<StereoChemistry>(m_database->GetColumn<int>(6));
        model_obj->GetChemicalComponentEntryMap().at(component_key)->AddComponentBondEntry(bond_key, bond_entry);
        model_obj->GetBondKeySystemPtr()->RegisterBond(bond_entry.bond_id, bond_key);
    }
}

std::vector<std::unique_ptr<AtomObject>> ModelObjectDAOv2::LoadAtomObjectList(const std::string & key_tag)
{
    auto local_entry_map{ LoadAtomLocalPotentialEntryMap(key_tag) };
    std::vector<std::unique_ptr<AtomObject>> atom_object_list;

    m_database->Prepare(std::string(SELECT_MODEL_ATOM_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto atom_object{ std::make_unique<AtomObject>() };
        atom_object->SetSerialID(m_database->GetColumn<int>(0));
        atom_object->SetSequenceID(m_database->GetColumn<int>(1));
        atom_object->SetComponentID(m_database->GetColumn<std::string>(2));
        atom_object->SetAtomID(m_database->GetColumn<std::string>(3));
        atom_object->SetChainID(m_database->GetColumn<std::string>(4));
        atom_object->SetIndicator(m_database->GetColumn<std::string>(5));
        atom_object->SetOccupancy(static_cast<float>(m_database->GetColumn<double>(6)));
        atom_object->SetTemperature(static_cast<float>(m_database->GetColumn<double>(7)));
        atom_object->SetElement(static_cast<Element>(m_database->GetColumn<int>(8)));
        atom_object->SetStructure(static_cast<Structure>(m_database->GetColumn<int>(9)));
        atom_object->SetSpecialAtomFlag(static_cast<bool>(m_database->GetColumn<int>(10)));
        atom_object->SetPosition(
            static_cast<float>(m_database->GetColumn<double>(11)),
            static_cast<float>(m_database->GetColumn<double>(12)),
            static_cast<float>(m_database->GetColumn<double>(13)));
        atom_object->SetComponentKey(m_database->GetColumn<ComponentKey>(14));
        atom_object->SetAtomKey(m_database->GetColumn<AtomKey>(15));

        const auto serial_id{ atom_object->GetSerialID() };
        auto iter{ local_entry_map.find(serial_id) };
        if (iter != local_entry_map.end())
        {
            atom_object->SetSelectedFlag(true);
            atom_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            atom_object->SetSelectedFlag(false);
        }
        atom_object_list.emplace_back(std::move(atom_object));
    }
    return atom_object_list;
}

std::vector<std::unique_ptr<BondObject>> ModelObjectDAOv2::LoadBondObjectList(
    const std::string & key_tag, const ModelObject * model_obj)
{
    auto local_entry_map{ LoadBondLocalPotentialEntryMap(key_tag) };
    std::vector<std::unique_ptr<BondObject>> bond_object_list;

    m_database->Prepare(std::string(SELECT_MODEL_BOND_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto atom_object_1{ model_obj->GetAtomPtr(m_database->GetColumn<int>(0)) };
        auto atom_object_2{ model_obj->GetAtomPtr(m_database->GetColumn<int>(1)) };
        auto bond_object{ std::make_unique<BondObject>(atom_object_1, atom_object_2) };
        bond_object->SetBondKey(m_database->GetColumn<BondKey>(2));
        bond_object->SetBondType(static_cast<BondType>(m_database->GetColumn<int>(3)));
        bond_object->SetBondOrder(static_cast<BondOrder>(m_database->GetColumn<int>(4)));
        bond_object->SetSpecialBondFlag(static_cast<bool>(m_database->GetColumn<int>(5)));

        const auto serial_id_pair{
            std::make_pair(atom_object_1->GetSerialID(), atom_object_2->GetSerialID()) };
        auto iter{ local_entry_map.find(serial_id_pair) };
        if (iter != local_entry_map.end())
        {
            bond_object->SetSelectedFlag(true);
            bond_object->AddLocalPotentialEntry(std::move(iter->second));
        }
        else
        {
            bond_object->SetSelectedFlag(false);
        }

        bond_object_list.emplace_back(std::move(bond_object));
    }
    return bond_object_list;
}

std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>>
ModelObjectDAOv2::LoadAtomLocalPotentialEntryMap(const std::string & key_tag)
{
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> entry_map;
    m_database->Prepare(std::string(SELECT_MODEL_ATOM_LOCAL_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto entry{ std::make_unique<LocalPotentialEntry>() };
        const auto serial_id{ m_database->GetColumn<int>(0) };
        entry->AddDistanceAndMapValueList(m_database->GetColumn<std::vector<std::tuple<float, float>>>(2));
        entry->AddGausEstimateOLS(m_database->GetColumn<double>(3), m_database->GetColumn<double>(4));
        entry->AddGausEstimateMDPDE(m_database->GetColumn<double>(5), m_database->GetColumn<double>(6));
        entry->SetAlphaR(m_database->GetColumn<double>(7));
        entry_map[serial_id] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        LoadAtomLocalPotentialEntrySubList(
            key_tag, ChemicalDataHelper::GetGroupAtomClassKey(i), entry_map);
    }
    return entry_map;
}

std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>>
ModelObjectDAOv2::LoadBondLocalPotentialEntryMap(const std::string & key_tag)
{
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> entry_map;
    m_database->Prepare(std::string(SELECT_MODEL_BOND_LOCAL_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        auto entry{ std::make_unique<LocalPotentialEntry>() };
        const auto key{ std::make_pair(m_database->GetColumn<int>(0), m_database->GetColumn<int>(1)) };
        entry->AddDistanceAndMapValueList(m_database->GetColumn<std::vector<std::tuple<float, float>>>(3));
        entry->AddGausEstimateOLS(m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        entry->AddGausEstimateMDPDE(m_database->GetColumn<double>(6), m_database->GetColumn<double>(7));
        entry->SetAlphaR(m_database->GetColumn<double>(8));
        entry_map[key] = std::move(entry);
    }

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        LoadBondLocalPotentialEntrySubList(
            key_tag, ChemicalDataHelper::GetGroupBondClassKey(i), entry_map);
    }
    return entry_map;
}

void ModelObjectDAOv2::LoadAtomLocalPotentialEntrySubList(
    const std::string & key_tag,
    const std::string & class_key,
    std::unordered_map<int, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    m_database->Prepare(std::string(SELECT_MODEL_ATOM_POSTERIOR_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto serial_id{ m_database->GetColumn<int>(0) };
        auto iter{ entry_map.find(serial_id) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(class_key, m_database->GetColumn<double>(1), m_database->GetColumn<double>(2));
        entry->AddGausVariancePosterior(class_key, m_database->GetColumn<double>(3), m_database->GetColumn<double>(4));
        entry->AddOutlierTag(class_key, static_cast<bool>(m_database->GetColumn<int>(5)));
        entry->AddStatisticalDistance(class_key, m_database->GetColumn<double>(6));
    }
}

void ModelObjectDAOv2::LoadBondLocalPotentialEntrySubList(
    const std::string & key_tag,
    const std::string & class_key,
    std::map<std::pair<int, int>, std::unique_ptr<LocalPotentialEntry>> & entry_map)
{
    m_database->Prepare(std::string(SELECT_MODEL_BOND_POSTERIOR_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto atom_pair{ std::make_pair(m_database->GetColumn<int>(0), m_database->GetColumn<int>(1)) };
        auto iter{ entry_map.find(atom_pair) };
        if (iter == entry_map.end()) continue;
        auto & entry{ iter->second };
        entry->AddGausEstimatePosterior(class_key, m_database->GetColumn<double>(2), m_database->GetColumn<double>(3));
        entry->AddGausVariancePosterior(class_key, m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        entry->AddOutlierTag(class_key, static_cast<bool>(m_database->GetColumn<int>(6)));
        entry->AddStatisticalDistance(class_key, m_database->GetColumn<double>(7));
    }
}

void ModelObjectDAOv2::LoadAtomGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    auto group_entry{ model_obj->GetAtomGroupPotentialEntry(class_key) };
    m_database->Prepare(std::string(SELECT_MODEL_ATOM_GROUP_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto group_key{ m_database->GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveAtomObjectPtrList(group_key, m_database->GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, m_database->GetColumn<double>(2), m_database->GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, m_database->GetColumn<double>(6), m_database->GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, m_database->GetColumn<double>(8), m_database->GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, m_database->GetColumn<double>(10));
    }

    for (auto & atom : model_obj->GetSelectedAtomList())
    {
        group_entry->AddAtomObjectPtr(AtomClassifier::GetGroupKeyInClass(atom, class_key), atom);
    }
}

void ModelObjectDAOv2::LoadBondGroupPotentialEntryList(
    ModelObject * model_obj, const std::string & key_tag, const std::string & class_key)
{
    auto group_entry{ model_obj->GetBondGroupPotentialEntry(class_key) };
    m_database->Prepare(std::string(SELECT_MODEL_BOND_GROUP_SQL));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, class_key);
    while (true)
    {
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            break;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
        }

        const auto group_key{ m_database->GetColumn<GroupKey>(0) };
        group_entry->InsertGroupKey(group_key);
        group_entry->ReserveBondObjectPtrList(group_key, m_database->GetColumn<int>(1));
        group_entry->AddGausEstimateMean(group_key, m_database->GetColumn<double>(2), m_database->GetColumn<double>(3));
        group_entry->AddGausEstimateMDPDE(group_key, m_database->GetColumn<double>(4), m_database->GetColumn<double>(5));
        group_entry->AddGausEstimatePrior(group_key, m_database->GetColumn<double>(6), m_database->GetColumn<double>(7));
        group_entry->AddGausVariancePrior(group_key, m_database->GetColumn<double>(8), m_database->GetColumn<double>(9));
        group_entry->AddAlphaG(group_key, m_database->GetColumn<double>(10));
    }

    for (auto & bond : model_obj->GetSelectedBondList())
    {
        group_entry->AddBondObjectPtr(BondClassifier::GetGroupKeyInClass(bond, class_key), bond);
    }
}

} // namespace rhbm_gem
