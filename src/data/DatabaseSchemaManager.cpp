#include "DatabaseSchemaManager.hpp"

#include "DataObjectBase.hpp"
#include "LegacyModelObjectDAO.hpp"
#include "Logger.hpp"
#include "MapObjectDAO.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"

#include "ChemicalDataHelper.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace
{
    constexpr std::string_view TABLE_EXISTS_SQL =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";
    constexpr std::string_view LIST_TABLES_SQL =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";
    constexpr std::string_view DATABASE_EMPTY_SQL =
        "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";
    constexpr std::string_view PRAGMA_USER_VERSION_SQL = "PRAGMA user_version;";
    constexpr std::string_view CREATE_OBJECT_METADATA_SQL =
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);";
    constexpr std::string_view UPSERT_OBJECT_METADATA_SQL = R"sql(
        INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?)
        ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;
    )sql";
    constexpr std::string_view LEGACY_MODEL_KEY_LIST_SQL =
        "SELECT key_tag FROM model_list;";
    constexpr std::string_view MODEL_METADATA_KEY_LIST_SQL =
        "SELECT key_tag FROM object_metadata WHERE object_type = 'model';";
    constexpr std::string_view MAP_METADATA_KEY_LIST_SQL =
        "SELECT key_tag FROM object_metadata WHERE object_type = 'map';";
    constexpr std::string_view MAP_KEY_LIST_SQL =
        "SELECT key_tag FROM map_list;";
    constexpr std::string_view MODEL_V2_KEY_LIST_SQL =
        "SELECT key_tag FROM model_object;";
    constexpr std::string_view DELETE_OBJECT_METADATA_SQL =
        "DELETE FROM object_metadata WHERE key_tag = ? AND object_type = ?;";

    constexpr std::array<std::string_view, 15> kNormalizedV2RequiredTables{
        "object_metadata",
        "model_object",
        "model_chain_map",
        "model_component",
        "model_component_atom",
        "model_component_bond",
        "model_atom",
        "model_bond",
        "model_atom_local_potential",
        "model_bond_local_potential",
        "model_atom_posterior",
        "model_bond_posterior",
        "model_atom_group_potential",
        "model_bond_group_potential",
        "map_list"
    };

    std::vector<std::string> QuerySingleStringColumn(
        rhbm_gem::SQLiteWrapper & database, const std::string & sql)
    {
        std::vector<std::string> values;
        database.Prepare(sql);
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        while (true)
        {
            const auto rc{ database.StepNext() };
            if (rc == rhbm_gem::SQLiteWrapper::StepDone())
            {
                break;
            }
            if (rc != rhbm_gem::SQLiteWrapper::StepRow())
            {
                throw std::runtime_error("Step failed: " + database.ErrorMessage());
            }
            values.push_back(database.GetColumn<std::string>(0));
        }
        return values;
    }

    int QueryUserVersion(rhbm_gem::SQLiteWrapper & database)
    {
        database.Prepare(PRAGMA_USER_VERSION_SQL.data());
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            return 0;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Failed to query schema version: " + database.ErrorMessage());
        }
        return database.GetColumn<int>(0);
    }

    void UpsertObjectMetadata(
        rhbm_gem::SQLiteWrapper & database,
        const std::vector<std::string> & key_list,
        const std::string & object_type)
    {
        if (key_list.empty())
        {
            return;
        }
        database.Prepare(std::string(UPSERT_OBJECT_METADATA_SQL));
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        for (const auto & key_tag : key_list)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, object_type);
            database.StepOnce();
            database.Reset();
        }
    }

    bool IsManagedTableName(std::string_view table_name)
    {
        return std::find(
                   kNormalizedV2RequiredTables.begin(),
                   kNormalizedV2RequiredTables.end(),
                   table_name) != kNormalizedV2RequiredTables.end();
    }

    std::string SanitizeLegacyModelKey(const std::string & key_tag)
    {
        std::string sanitized_key_tag;
        sanitized_key_tag.reserve(key_tag.size());
        for (char ch : key_tag)
        {
            if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_')
            {
                sanitized_key_tag.push_back(ch);
            }
            else
            {
                sanitized_key_tag.push_back('_');
            }
        }
        return sanitized_key_tag;
    }
}

namespace rhbm_gem {

DatabaseSchemaManager::DatabaseSchemaManager(SQLiteWrapper * database) :
    m_database{ database }
{
    if (m_database == nullptr)
    {
        throw std::runtime_error("DatabaseSchemaManager requires a valid database handle.");
    }
}

DatabaseSchemaVersion DatabaseSchemaManager::EnsureSchema()
{
    const auto raw_version{ QueryUserVersion(*m_database) };
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::LegacyV1))
    {
        MigrateLegacyV1ToNormalizedV2();
        return DatabaseSchemaVersion::NormalizedV2;
    }
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        ValidateNormalizedV2Schema();
        return DatabaseSchemaVersion::NormalizedV2;
    }
    if (raw_version != 0)
    {
        throw std::runtime_error("Unsupported database schema version: " + std::to_string(raw_version));
    }

    switch (InspectSchemaState())
    {
    case DatabaseSchemaState::Empty:
        EnsureNormalizedV2Tables();
        SetSchemaVersion(DatabaseSchemaVersion::NormalizedV2);
        ValidateNormalizedV2Schema();
        return DatabaseSchemaVersion::NormalizedV2;
    case DatabaseSchemaState::LegacyV1:
        MigrateLegacyV1ToNormalizedV2();
        return DatabaseSchemaVersion::NormalizedV2;
    case DatabaseSchemaState::ManagedButUnversioned:
        EnsureNormalizedV2Tables();
        RepairManagedMetadata();
        SetSchemaVersion(DatabaseSchemaVersion::NormalizedV2);
        ValidateNormalizedV2Schema();
        return DatabaseSchemaVersion::NormalizedV2;
    case DatabaseSchemaState::NormalizedV2:
        ValidateNormalizedV2Schema();
        return DatabaseSchemaVersion::NormalizedV2;
    case DatabaseSchemaState::MixedUnknown:
        throw std::runtime_error(
            "Database contains unmanaged or mixed schema tables and cannot be auto-upgraded safely.");
    }

    throw std::runtime_error("Unhandled database schema state.");
}

DatabaseSchemaVersion DatabaseSchemaManager::GetSchemaVersion() const
{
    return QueryUserVersion(*m_database) == static_cast<int>(DatabaseSchemaVersion::NormalizedV2)
        ? DatabaseSchemaVersion::NormalizedV2
        : DatabaseSchemaVersion::LegacyV1;
}

DatabaseSchemaState DatabaseSchemaManager::InspectSchemaState() const
{
    const auto raw_version{ QueryUserVersion(*m_database) };
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::LegacyV1))
    {
        return DatabaseSchemaState::LegacyV1;
    }
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        return DatabaseSchemaState::NormalizedV2;
    }
    if (IsDatabaseEmpty())
    {
        return DatabaseSchemaState::Empty;
    }
    if (HasLegacyModelSchema())
    {
        return DatabaseSchemaState::LegacyV1;
    }

    for (const auto & table_name : QuerySingleStringColumn(*m_database, std::string(LIST_TABLES_SQL)))
    {
        if (!IsManagedTableName(table_name))
        {
            return DatabaseSchemaState::MixedUnknown;
        }
    }
    return DatabaseSchemaState::ManagedButUnversioned;
}

bool DatabaseSchemaManager::IsDatabaseEmpty() const
{
    m_database->Prepare(DATABASE_EMPTY_SQL.data());
    SQLiteWrapper::StatementGuard guard(*m_database);
    const auto rc{ m_database->StepNext() };
    if (rc != SQLiteWrapper::StepRow())
    {
        if (rc == SQLiteWrapper::StepDone()) return true;
        throw std::runtime_error("Failed to inspect database tables: " + m_database->ErrorMessage());
    }
    return m_database->GetColumn<int>(0) == 0;
}

bool DatabaseSchemaManager::HasTable(const std::string & table_name) const
{
    m_database->Prepare(TABLE_EXISTS_SQL.data());
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, table_name);
    const auto rc{ m_database->StepNext() };
    if (rc == SQLiteWrapper::StepDone()) return false;
    if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Failed to inspect table existence: " + m_database->ErrorMessage());
    }
    return true;
}

bool DatabaseSchemaManager::HasLegacyModelSchema() const
{
    return HasTable("model_list");
}

void DatabaseSchemaManager::ValidateNormalizedV2Schema() const
{
    for (const auto required_table : kNormalizedV2RequiredTables)
    {
        if (!HasTable(std::string(required_table)))
        {
            throw std::runtime_error(
                "Normalized v2 schema is missing required table: " + std::string(required_table));
        }
    }
}

std::vector<std::string> DatabaseSchemaManager::GetLegacyModelKeyList() const
{
    if (!HasTable("model_list"))
    {
        return {};
    }

    const auto model_list_keys{
        QuerySingleStringColumn(*m_database, std::string(LEGACY_MODEL_KEY_LIST_SQL)) };
    if (!HasTable("object_metadata"))
    {
        return model_list_keys;
    }

    const auto metadata_keys{
        QuerySingleStringColumn(*m_database, std::string(MODEL_METADATA_KEY_LIST_SQL)) };
    std::unordered_set<std::string> metadata_key_set{ metadata_keys.begin(), metadata_keys.end() };
    std::unordered_set<std::string> model_list_key_set{ model_list_keys.begin(), model_list_keys.end() };

    for (const auto & key_tag : model_list_keys)
    {
        if (metadata_key_set.find(key_tag) == metadata_key_set.end())
        {
            Logger::Log(
                LogLevel::Warning,
                "Legacy database metadata is missing model key '" + key_tag +
                    "'. Migration will use model_list as the source of truth.");
        }
    }

    for (const auto & key_tag : metadata_keys)
    {
        if (model_list_key_set.find(key_tag) == model_list_key_set.end())
        {
            Logger::Log(
                LogLevel::Warning,
                "Legacy database metadata references model key '" + key_tag +
                    "' that is absent from model_list. The ghost metadata row will be ignored.");
        }
    }

    return model_list_keys;
}

std::vector<std::string> DatabaseSchemaManager::BuildLegacyOwnedTableNames(
    const std::vector<std::string> & model_key_list) const
{
    std::unordered_set<std::string> table_names;
    table_names.insert("model_list");

    for (const auto & key_tag : model_key_list)
    {
        const auto sanitized_key_tag{ SanitizeLegacyModelKey(key_tag) };
        table_names.insert("chemical_component_entry_in_" + sanitized_key_tag);
        table_names.insert("component_atom_entry_in_" + sanitized_key_tag);
        table_names.insert("component_bond_entry_in_" + sanitized_key_tag);
        table_names.insert("atom_list_in_" + sanitized_key_tag);
        table_names.insert("bond_list_in_" + sanitized_key_tag);
        table_names.insert("atom_local_potential_entry_in_" + sanitized_key_tag);
        table_names.insert("bond_local_potential_entry_in_" + sanitized_key_tag);

        for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
        {
            const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
            table_names.insert(
                class_key + "_atom_local_potential_entry_in_" + sanitized_key_tag);
            table_names.insert(
                class_key + "_atom_group_potential_entry_in_" + sanitized_key_tag);
        }
        for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
        {
            const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
            table_names.insert(
                class_key + "_bond_local_potential_entry_in_" + sanitized_key_tag);
            table_names.insert(
                class_key + "_bond_group_potential_entry_in_" + sanitized_key_tag);
        }
    }

    std::vector<std::string> owned_table_names{ table_names.begin(), table_names.end() };
    std::sort(owned_table_names.begin(), owned_table_names.end());
    return owned_table_names;
}

void DatabaseSchemaManager::EnsureNormalizedV2Tables() const
{
    m_database->Execute(std::string(CREATE_OBJECT_METADATA_SQL));
    ModelObjectDAOv2::EnsureSchema(*m_database);
    MapObjectDAO::EnsureSchema(*m_database);
}

void DatabaseSchemaManager::RepairManagedMetadata() const
{
    EnsureNormalizedV2Tables();

    std::vector<std::string> map_keys;
    if (HasTable("map_list"))
    {
        map_keys = QuerySingleStringColumn(*m_database, std::string(MAP_KEY_LIST_SQL));
        UpsertObjectMetadata(*m_database, map_keys, "map");
    }

    std::vector<std::string> model_keys;
    if (HasTable("model_object"))
    {
        model_keys = QuerySingleStringColumn(*m_database, std::string(MODEL_V2_KEY_LIST_SQL));
        UpsertObjectMetadata(*m_database, model_keys, "model");
    }

    if (HasTable("object_metadata"))
    {
        const std::unordered_set<std::string> map_key_set{ map_keys.begin(), map_keys.end() };
        const std::unordered_set<std::string> model_key_set{ model_keys.begin(), model_keys.end() };
        const auto map_metadata_keys{
            QuerySingleStringColumn(*m_database, std::string(MAP_METADATA_KEY_LIST_SQL)) };
        const auto model_metadata_keys{
            QuerySingleStringColumn(*m_database, std::string(MODEL_METADATA_KEY_LIST_SQL)) };
        m_database->Prepare(std::string(DELETE_OBJECT_METADATA_SQL));
        SQLiteWrapper::StatementGuard guard(*m_database);

        for (const auto & key_tag : map_metadata_keys)
        {
            if (map_key_set.find(key_tag) != map_key_set.end())
            {
                continue;
            }
            m_database->Bind<std::string>(1, key_tag);
            m_database->Bind<std::string>(2, "map");
            m_database->StepOnce();
            m_database->Reset();
        }

        for (const auto & key_tag : model_metadata_keys)
        {
            if (model_key_set.find(key_tag) != model_key_set.end())
            {
                continue;
            }
            m_database->Bind<std::string>(1, key_tag);
            m_database->Bind<std::string>(2, "model");
            m_database->StepOnce();
            m_database->Reset();
        }
    }
}

void DatabaseSchemaManager::MigrateLegacyV1ToNormalizedV2()
{
    Logger::Log(LogLevel::Info, "Migrating database schema from legacy v1 to normalized v2.");

    SQLiteWrapper::TransactionGuard transaction(*m_database);
    const auto legacy_model_keys{ GetLegacyModelKeyList() };
    EnsureNormalizedV2Tables();

    LegacyModelObjectDAO legacy_dao{ m_database };
    ModelObjectDAOv2 v2_dao{ m_database };

    for (const auto & key_tag : legacy_model_keys)
    {
        auto model_object{ legacy_dao.Load(key_tag) };
        v2_dao.Save(model_object.get(), key_tag);
        UpsertObjectMetadata(*m_database, { key_tag }, "model");
    }

    for (const auto & table_name : BuildLegacyOwnedTableNames(legacy_model_keys))
    {
        if (HasTable(table_name))
        {
            m_database->Execute("DROP TABLE IF EXISTS " + table_name + ";");
        }
    }

    RepairManagedMetadata();
    SetSchemaVersion(DatabaseSchemaVersion::NormalizedV2);
    ValidateNormalizedV2Schema();
}

void DatabaseSchemaManager::SetSchemaVersion(DatabaseSchemaVersion version) const
{
    m_database->Execute(
        "PRAGMA user_version = " + std::to_string(static_cast<int>(version)) + ";");
}

} // namespace rhbm_gem
