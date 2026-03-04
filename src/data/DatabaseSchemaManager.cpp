#include "DatabaseSchemaManager.hpp"

#include "Logger.hpp"
#include "MapObjectDAO.hpp"
#include "ModelObject.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectBase.hpp"

#include "ChemicalDataHelper.hpp"
#include "legacy/LegacyModelObjectReader.hpp"

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
    enum class UnversionedSchemaState
    {
        Empty,
        LegacyV1,
        ManagedButUnversioned,
        MixedUnknown
    };

    constexpr std::string_view kTableExistsSql =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";
    constexpr std::string_view kListTablesSql =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";
    constexpr std::string_view kDatabaseEmptySql =
        "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";
    constexpr std::string_view kPragmaUserVersionSql = "PRAGMA user_version;";
    constexpr std::string_view kCreateObjectMetadataSql =
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);";
    constexpr std::string_view kUpsertObjectMetadataSql = R"sql(
        INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?)
        ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;
    )sql";
    constexpr std::string_view kLegacyModelKeyListSql = "SELECT key_tag FROM model_list;";
    constexpr std::string_view kModelMetadataKeyListSql =
        "SELECT key_tag FROM object_metadata WHERE object_type = 'model';";
    constexpr std::string_view kMapMetadataKeyListSql =
        "SELECT key_tag FROM object_metadata WHERE object_type = 'map';";
    constexpr std::string_view kMapKeyListSql = "SELECT key_tag FROM map_list;";
    constexpr std::string_view kModelV2KeyListSql = "SELECT key_tag FROM model_object;";
    constexpr std::string_view kDeleteObjectMetadataSql =
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
        rhbm_gem::SQLiteWrapper & database, std::string_view sql)
    {
        std::vector<std::string> values;
        database.Prepare(std::string(sql));
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
        database.Prepare(std::string(kPragmaUserVersionSql));
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
        database.Prepare(std::string(kUpsertObjectMetadataSql));
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        for (const auto & key_tag : key_list)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, object_type);
            database.StepOnce();
            database.Reset();
        }
    }

    bool HasTable(rhbm_gem::SQLiteWrapper & database, const std::string & table_name)
    {
        database.Prepare(std::string(kTableExistsSql));
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        database.Bind<std::string>(1, table_name);
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            return false;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Failed to inspect table existence: " + database.ErrorMessage());
        }
        return true;
    }

    bool IsDatabaseEmpty(rhbm_gem::SQLiteWrapper & database)
    {
        database.Prepare(std::string(kDatabaseEmptySql));
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        const auto rc{ database.StepNext() };
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            if (rc == rhbm_gem::SQLiteWrapper::StepDone())
            {
                return true;
            }
            throw std::runtime_error("Failed to inspect database tables: " + database.ErrorMessage());
        }
        return database.GetColumn<int>(0) == 0;
    }

    bool HasLegacyModelSchema(rhbm_gem::SQLiteWrapper & database)
    {
        return HasTable(database, "model_list");
    }

    bool IsManagedTableName(std::string_view table_name)
    {
        return std::find(
                   kNormalizedV2RequiredTables.begin(),
                   kNormalizedV2RequiredTables.end(),
                   table_name) != kNormalizedV2RequiredTables.end();
    }

    UnversionedSchemaState InspectUnversionedSchemaState(rhbm_gem::SQLiteWrapper & database)
    {
        if (IsDatabaseEmpty(database))
        {
            return UnversionedSchemaState::Empty;
        }
        if (HasLegacyModelSchema(database))
        {
            return UnversionedSchemaState::LegacyV1;
        }

        for (const auto & table_name : QuerySingleStringColumn(database, kListTablesSql))
        {
            if (!IsManagedTableName(table_name))
            {
                return UnversionedSchemaState::MixedUnknown;
            }
        }
        return UnversionedSchemaState::ManagedButUnversioned;
    }

    void ValidateNormalizedV2Schema(rhbm_gem::SQLiteWrapper & database)
    {
        for (const auto required_table : kNormalizedV2RequiredTables)
        {
            if (!HasTable(database, std::string(required_table)))
            {
                throw std::runtime_error(
                    "Normalized v2 schema is missing required table: " +
                    std::string(required_table));
            }
        }
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

    std::vector<std::string> GetLegacyModelKeyList(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, "model_list"))
        {
            return {};
        }

        const auto model_list_keys{ QuerySingleStringColumn(database, kLegacyModelKeyListSql) };
        if (!HasTable(database, "object_metadata"))
        {
            return model_list_keys;
        }

        const auto metadata_keys{ QuerySingleStringColumn(database, kModelMetadataKeyListSql) };
        std::unordered_set<std::string> metadata_key_set{ metadata_keys.begin(), metadata_keys.end() };
        std::unordered_set<std::string> model_list_key_set{
            model_list_keys.begin(), model_list_keys.end() };

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

    std::vector<std::string> BuildLegacyOwnedTableNames(
        const std::vector<std::string> & model_key_list)
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
                table_names.insert(class_key + "_atom_local_potential_entry_in_" + sanitized_key_tag);
                table_names.insert(class_key + "_atom_group_potential_entry_in_" + sanitized_key_tag);
            }
            for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
            {
                const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
                table_names.insert(class_key + "_bond_local_potential_entry_in_" + sanitized_key_tag);
                table_names.insert(class_key + "_bond_group_potential_entry_in_" + sanitized_key_tag);
            }
        }

        std::vector<std::string> owned_table_names{ table_names.begin(), table_names.end() };
        std::sort(owned_table_names.begin(), owned_table_names.end());
        return owned_table_names;
    }

    void EnsureNormalizedV2Tables(rhbm_gem::SQLiteWrapper & database)
    {
        database.Execute(std::string(kCreateObjectMetadataSql));
        rhbm_gem::ModelObjectDAOv2::EnsureSchema(database);
        rhbm_gem::MapObjectDAO::EnsureSchema(database);
    }

    void RepairManagedMetadata(rhbm_gem::SQLiteWrapper & database)
    {
        EnsureNormalizedV2Tables(database);

        std::vector<std::string> map_keys;
        if (HasTable(database, "map_list"))
        {
            map_keys = QuerySingleStringColumn(database, kMapKeyListSql);
            UpsertObjectMetadata(database, map_keys, "map");
        }

        std::vector<std::string> model_keys;
        if (HasTable(database, "model_object"))
        {
            model_keys = QuerySingleStringColumn(database, kModelV2KeyListSql);
            UpsertObjectMetadata(database, model_keys, "model");
        }

        if (!HasTable(database, "object_metadata"))
        {
            return;
        }

        const std::unordered_set<std::string> map_key_set{ map_keys.begin(), map_keys.end() };
        const std::unordered_set<std::string> model_key_set{ model_keys.begin(), model_keys.end() };
        const auto map_metadata_keys{ QuerySingleStringColumn(database, kMapMetadataKeyListSql) };
        const auto model_metadata_keys{ QuerySingleStringColumn(database, kModelMetadataKeyListSql) };
        database.Prepare(std::string(kDeleteObjectMetadataSql));
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);

        for (const auto & key_tag : map_metadata_keys)
        {
            if (map_key_set.find(key_tag) != map_key_set.end())
            {
                continue;
            }
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, "map");
            database.StepOnce();
            database.Reset();
        }

        for (const auto & key_tag : model_metadata_keys)
        {
            if (model_key_set.find(key_tag) != model_key_set.end())
            {
                continue;
            }
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, "model");
            database.StepOnce();
            database.Reset();
        }
    }

    void SetSchemaVersion(
        rhbm_gem::SQLiteWrapper & database,
        rhbm_gem::DatabaseSchemaVersion version)
    {
        database.Execute(
            "PRAGMA user_version = " +
            std::to_string(static_cast<int>(version)) +
            ";");
    }

    void MigrateLegacyV1ToNormalizedV2(rhbm_gem::SQLiteWrapper & database)
    {
        Logger::Log(
            LogLevel::Info,
            "Migrating database schema from legacy v1 to normalized v2.");

        rhbm_gem::SQLiteWrapper::TransactionGuard transaction(database);
        const auto legacy_model_keys{ GetLegacyModelKeyList(database) };
        EnsureNormalizedV2Tables(database);

        rhbm_gem::LegacyModelObjectReader legacy_reader{ &database };
        rhbm_gem::ModelObjectDAOv2 v2_dao{ &database };
        for (const auto & key_tag : legacy_model_keys)
        {
            auto model_object{ legacy_reader.Load(key_tag) };
            v2_dao.Save(model_object.get(), key_tag);
            UpsertObjectMetadata(database, { key_tag }, "model");
        }

        for (const auto & table_name : BuildLegacyOwnedTableNames(legacy_model_keys))
        {
            if (HasTable(database, table_name))
            {
                database.Execute("DROP TABLE IF EXISTS " + table_name + ";");
            }
        }

        RepairManagedMetadata(database);
        SetSchemaVersion(database, rhbm_gem::DatabaseSchemaVersion::NormalizedV2);
        ValidateNormalizedV2Schema(database);
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
        MigrateLegacyV1ToNormalizedV2(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    }
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        ValidateNormalizedV2Schema(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    }
    if (raw_version != 0)
    {
        throw std::runtime_error(
            "Unsupported database schema version: " + std::to_string(raw_version));
    }

    switch (InspectUnversionedSchemaState(*m_database))
    {
    case UnversionedSchemaState::Empty:
        EnsureNormalizedV2Tables(*m_database);
        SetSchemaVersion(*m_database, DatabaseSchemaVersion::NormalizedV2);
        ValidateNormalizedV2Schema(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::LegacyV1:
        MigrateLegacyV1ToNormalizedV2(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::ManagedButUnversioned:
        EnsureNormalizedV2Tables(*m_database);
        RepairManagedMetadata(*m_database);
        SetSchemaVersion(*m_database, DatabaseSchemaVersion::NormalizedV2);
        ValidateNormalizedV2Schema(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::MixedUnknown:
        throw std::runtime_error(
            "Database contains unmanaged or mixed schema tables and cannot be auto-upgraded safely.");
    }

    throw std::runtime_error("Unhandled database schema state.");
}

} // namespace rhbm_gem
