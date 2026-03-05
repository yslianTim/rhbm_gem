#include "DatabaseSchemaManager.hpp"

#include "Logger.hpp"
#include "ModelObject.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"
#include "model_io/ModelSchemaSql.hpp"
#include "persistence/ManagedStoreRegistry.hpp"
#include "persistence/MapStoreSql.hpp"

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
    constexpr std::string_view kLegacyModelKeyListSql = "SELECT key_tag FROM model_list;";
    constexpr std::string_view kLegacyModelMetadataKeyListSql =
        "SELECT key_tag FROM object_metadata WHERE object_type = 'model' ORDER BY key_tag;";
    constexpr std::string_view kLegacyMapMetadataKeyListSql =
        "SELECT key_tag FROM object_metadata WHERE object_type = 'map' ORDER BY key_tag;";

    constexpr std::string_view kFinalCatalogTableName = "object_catalog";
    constexpr std::string_view kObjectCatalogShadowSuffix = "_final_v2";
    constexpr std::string_view kManagedObjectMetadataTableName = "object_metadata";

    constexpr std::array<std::string_view, 16> kKnownManagedTableNames{
        "object_catalog",
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
        rhbm_gem::SQLiteWrapper & database,
        const std::string & sql)
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

    int QuerySingleInt(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & sql)
    {
        database.Prepare(sql);
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        const auto rc{ database.StepNext() };
        if (rc == rhbm_gem::SQLiteWrapper::StepDone())
        {
            return 0;
        }
        if (rc != rhbm_gem::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Step failed: " + database.ErrorMessage());
        }
        return database.GetColumn<int>(0);
    }

    int QueryUserVersion(rhbm_gem::SQLiteWrapper & database)
    {
        return QuerySingleInt(database, std::string(kPragmaUserVersionSql));
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

    void DropTableIfExists(rhbm_gem::SQLiteWrapper & database, const std::string & table_name)
    {
        if (HasTable(database, table_name))
        {
            database.Execute("DROP TABLE " + table_name + ";");
        }
    }

    bool IsDatabaseEmpty(rhbm_gem::SQLiteWrapper & database)
    {
        return QuerySingleInt(database, std::string(kDatabaseEmptySql)) == 0;
    }

    bool HasLegacyModelSchema(rhbm_gem::SQLiteWrapper & database)
    {
        return HasTable(database, "model_list");
    }

    bool IsManagedTableName(std::string_view table_name)
    {
        return std::find(
                   kKnownManagedTableNames.begin(),
                   kKnownManagedTableNames.end(),
                   table_name) != kKnownManagedTableNames.end();
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

        for (const auto & table_name : QuerySingleStringColumn(database, std::string(kListTablesSql)))
        {
            if (!IsManagedTableName(table_name))
            {
                return UnversionedSchemaState::MixedUnknown;
            }
        }
        return UnversionedSchemaState::ManagedButUnversioned;
    }

    std::string BuildObjectCatalogTableName(const std::string & suffix = {})
    {
        return std::string(kFinalCatalogTableName) + suffix;
    }

    void CreateObjectCatalogTable(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & table_name = std::string(kFinalCatalogTableName))
    {
        database.Execute(
            "CREATE TABLE IF NOT EXISTS " + table_name +
            " (key_tag TEXT PRIMARY KEY, object_type TEXT NOT NULL, "
            "CHECK (object_type IN ('model', 'map')));");
    }

    void UpsertObjectCatalogRows(
        rhbm_gem::SQLiteWrapper & database,
        const std::vector<std::string> & key_list,
        const std::string & object_type,
        const std::string & table_name = std::string(kFinalCatalogTableName))
    {
        if (key_list.empty())
        {
            return;
        }
        database.Prepare(
            "INSERT INTO " + table_name + "(key_tag, object_type) VALUES (?, ?) "
            "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;");
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        for (const auto & key_tag : key_list)
        {
            database.Bind<std::string>(1, key_tag);
            database.Bind<std::string>(2, object_type);
            database.StepOnce();
            database.Reset();
        }
    }

    std::vector<std::string> QueryCatalogKeys(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & object_type,
        const std::string & table_name = std::string(kFinalCatalogTableName))
    {
        if (!HasTable(database, table_name))
        {
            return {};
        }
        database.Prepare(
            "SELECT key_tag FROM " + table_name + " WHERE object_type = ? ORDER BY key_tag;");
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        database.Bind<std::string>(1, object_type);
        std::vector<std::string> key_list;
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
            key_list.push_back(database.GetColumn<std::string>(0));
        }
        return key_list;
    }

    void ValidateFinalV2Schema(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, std::string(kFinalCatalogTableName)))
        {
            throw std::runtime_error("Normalized v2 schema is missing required table: object_catalog");
        }
        if (HasTable(database, std::string(kManagedObjectMetadataTableName)))
        {
            throw std::runtime_error(
                "Normalized v2 schema must not retain legacy object_metadata.");
        }

        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            descriptor.validate_schema_v2(database);
            const auto payload_keys{ descriptor.list_keys(database) };
            const auto catalog_keys{
                QueryCatalogKeys(database, std::string(descriptor.object_type))
            };
            const std::unordered_set<std::string> payload_key_set{
                payload_keys.begin(), payload_keys.end()
            };
            const std::unordered_set<std::string> catalog_key_set{
                catalog_keys.begin(), catalog_keys.end()
            };

            for (const auto & key_tag : payload_keys)
            {
                if (catalog_key_set.find(key_tag) == catalog_key_set.end())
                {
                    throw std::runtime_error(
                        "Normalized v2 catalog is missing " +
                        std::string(descriptor.object_type) +
                        " key: " + key_tag);
                }
            }

            for (const auto & key_tag : catalog_keys)
            {
                if (payload_key_set.find(key_tag) == payload_key_set.end())
                {
                    throw std::runtime_error(
                        "Normalized v2 catalog contains ghost " +
                        std::string(descriptor.object_type) +
                        " key without payload: " + key_tag);
                }
            }
        }
    }

    bool IsCatalogBasedV2Schema(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, std::string(kFinalCatalogTableName)))
        {
            return false;
        }
        if (HasTable(database, std::string(kManagedObjectMetadataTableName)))
        {
            return false;
        }

        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            descriptor.validate_schema_v2(database);
        }
        return true;
    }

    bool IsMetadataBasedLegacyV2Schema(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, std::string(kManagedObjectMetadataTableName)))
        {
            return false;
        }
        if (HasTable(database, std::string(kFinalCatalogTableName)))
        {
            return false;
        }
        for (const auto & table_name : QuerySingleStringColumn(database, std::string(kListTablesSql)))
        {
            if (!IsManagedTableName(table_name))
            {
                return false;
            }
        }
        return true;
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

    int CountRows(rhbm_gem::SQLiteWrapper & database, const std::string & table_name)
    {
        if (!HasTable(database, table_name))
        {
            return 0;
        }
        return QuerySingleInt(database, "SELECT COUNT(*) FROM " + table_name + ";");
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

        const auto model_list_keys{
            QuerySingleStringColumn(database, std::string(kLegacyModelKeyListSql))
        };
        if (!HasTable(database, std::string(kManagedObjectMetadataTableName)))
        {
            return model_list_keys;
        }

        const auto metadata_keys{
            QuerySingleStringColumn(database, std::string(kLegacyModelMetadataKeyListSql))
        };
        std::unordered_set<std::string> metadata_key_set{ metadata_keys.begin(), metadata_keys.end() };
        std::unordered_set<std::string> model_list_key_set{
            model_list_keys.begin(), model_list_keys.end()
        };

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

    void CreateFinalV2Tables(rhbm_gem::SQLiteWrapper & database)
    {
        CreateObjectCatalogTable(database);
        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            descriptor.ensure_schema_v2(database);
        }
    }

    void SeedCatalogFromLegacyMetadata(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & catalog_table_name)
    {
        if (!HasTable(database, std::string(kManagedObjectMetadataTableName)))
        {
            return;
        }

        const auto model_metadata_keys{
            QuerySingleStringColumn(database, std::string(kLegacyModelMetadataKeyListSql))
        };
        const auto map_metadata_keys{
            QuerySingleStringColumn(database, std::string(kLegacyMapMetadataKeyListSql))
        };
        const auto model_payload_keys{
            rhbm_gem::persistence::LookupManagedStoreDescriptor("model").list_keys(database)
        };
        const auto map_payload_keys{
            rhbm_gem::persistence::LookupManagedStoreDescriptor("map").list_keys(database)
        };
        const std::unordered_set<std::string> model_payload_key_set{
            model_payload_keys.begin(), model_payload_keys.end()
        };
        const std::unordered_set<std::string> map_payload_key_set{
            map_payload_keys.begin(), map_payload_keys.end()
        };

        std::vector<std::string> surviving_model_keys;
        for (const auto & key_tag : model_metadata_keys)
        {
            if (model_payload_key_set.find(key_tag) == model_payload_key_set.end())
            {
                Logger::Log(
                    LogLevel::Warning,
                    "Ignoring ghost model metadata row during normalized v2 upgrade: " + key_tag);
                continue;
            }
            surviving_model_keys.push_back(key_tag);
        }

        std::vector<std::string> surviving_map_keys;
        for (const auto & key_tag : map_metadata_keys)
        {
            if (map_payload_key_set.find(key_tag) == map_payload_key_set.end())
            {
                Logger::Log(
                    LogLevel::Warning,
                    "Ignoring ghost map metadata row during normalized v2 upgrade: " + key_tag);
                continue;
            }
            surviving_map_keys.push_back(key_tag);
        }

        UpsertObjectCatalogRows(database, surviving_model_keys, "model", catalog_table_name);
        UpsertObjectCatalogRows(database, surviving_map_keys, "map", catalog_table_name);
    }

    void SupplementCatalogFromPayload(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & catalog_table_name)
    {
        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            UpsertObjectCatalogRows(
                database,
                descriptor.list_keys(database),
                std::string(descriptor.object_type),
                catalog_table_name);
        }
    }

    std::vector<std::string> QueryRootKeysForObjectType(
        rhbm_gem::SQLiteWrapper & database,
        std::string_view object_type,
        const std::string & suffix)
    {
        if (object_type == "model")
        {
            return HasTable(database, "model_object" + suffix)
                ? QuerySingleStringColumn(
                    database,
                    "SELECT key_tag FROM model_object" + suffix + " ORDER BY key_tag;")
                : std::vector<std::string>{};
        }
        if (object_type == "map")
        {
            return HasTable(database, "map_list" + suffix)
                ? QuerySingleStringColumn(
                    database,
                    "SELECT key_tag FROM map_list" + suffix + " ORDER BY key_tag;")
                : std::vector<std::string>{};
        }
        throw std::runtime_error("Unsupported managed object type: " + std::string(object_type));
    }

    void ValidateShadowCatalogCoverage(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & suffix)
    {
        const auto shadow_catalog_table_name{ BuildObjectCatalogTableName(suffix) };
        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            const auto catalog_keys{
                QueryCatalogKeys(database, std::string(descriptor.object_type), shadow_catalog_table_name)
            };
            const auto payload_keys{
                QueryRootKeysForObjectType(database, descriptor.object_type, suffix)
            };
            const std::unordered_set<std::string> catalog_key_set{
                catalog_keys.begin(), catalog_keys.end()
            };
            const std::unordered_set<std::string> payload_key_set{
                payload_keys.begin(), payload_keys.end()
            };

            for (const auto & key_tag : catalog_keys)
            {
                if (payload_key_set.find(key_tag) == payload_key_set.end())
                {
                    throw std::runtime_error(
                        "Shadow object catalog contains ghost " +
                        std::string(descriptor.object_type) +
                        " key without payload: " + key_tag);
                }
            }

            for (const auto & key_tag : payload_keys)
            {
                if (catalog_key_set.find(key_tag) == catalog_key_set.end())
                {
                    throw std::runtime_error(
                        "Shadow object catalog is missing " +
                        std::string(descriptor.object_type) +
                        " key: " + key_tag);
                }
            }
        }
    }

    void ValidateShadowPayloadRowCounts(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & suffix)
    {
        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            if (!HasTable(database, std::string(table_name)))
            {
                continue;
            }
            const auto expected_count{ CountRows(database, std::string(table_name)) };
            const auto actual_count{ CountRows(database, std::string(table_name) + suffix) };
            if (expected_count != actual_count)
            {
                throw std::runtime_error(
                    "Shadow table row count mismatch for " +
                    std::string(table_name) +
                    ": expected " + std::to_string(expected_count) +
                    ", actual " + std::to_string(actual_count));
            }
        }

        for (const auto table_name : rhbm_gem::persistence::kMapCanonicalTableNames)
        {
            if (!HasTable(database, std::string(table_name)))
            {
                continue;
            }
            const auto expected_count{ CountRows(database, std::string(table_name)) };
            const auto actual_count{ CountRows(database, std::string(table_name) + suffix) };
            if (expected_count != actual_count)
            {
                throw std::runtime_error(
                    "Shadow table row count mismatch for " +
                    std::string(table_name) +
                    ": expected " + std::to_string(expected_count) +
                    ", actual " + std::to_string(actual_count));
            }
        }
    }

    void RebuildManagedSchemaToFinalV2(
        rhbm_gem::SQLiteWrapper & database,
        bool seed_from_legacy_metadata)
    {
        Logger::Log(
            LogLevel::Info,
            "Rebuilding managed database schema into final normalized v2 layout.");

        rhbm_gem::SQLiteWrapper::TransactionGuard transaction(database);
        const std::string suffix{ kObjectCatalogShadowSuffix };
        const auto shadow_catalog_table_name{ BuildObjectCatalogTableName(suffix) };

        CreateObjectCatalogTable(database, shadow_catalog_table_name);
        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            descriptor.create_shadow_tables_v2(database, suffix);
        }

        if (seed_from_legacy_metadata)
        {
            SeedCatalogFromLegacyMetadata(database, shadow_catalog_table_name);
        }
        SupplementCatalogFromPayload(database, shadow_catalog_table_name);

        int rebuilt_table_count{ 0 };
        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            descriptor.copy_into_shadow_tables_v2(database, suffix);
            ++rebuilt_table_count;
        }

        ValidateShadowCatalogCoverage(database, suffix);
        ValidateShadowPayloadRowCounts(database, suffix);

        DropTableIfExists(database, std::string(kManagedObjectMetadataTableName));
        for (const auto & descriptor : rhbm_gem::persistence::GetAllManagedStoreDescriptors())
        {
            descriptor.drop_old_and_rename_shadow_tables_v2(database, suffix);
        }
        DropTableIfExists(database, std::string(kFinalCatalogTableName));
        database.Execute(
            "ALTER TABLE " + shadow_catalog_table_name +
            " RENAME TO " + std::string(kFinalCatalogTableName) + ";");

        SetSchemaVersion(database, rhbm_gem::DatabaseSchemaVersion::NormalizedV2);
        ValidateFinalV2Schema(database);
        Logger::Log(
            LogLevel::Info,
            "Managed schema rebuild completed. Rebuilt store count = " +
                std::to_string(rebuilt_table_count) + ".");
    }

    void MigrateLegacyV1ToFinalV2(rhbm_gem::SQLiteWrapper & database)
    {
        Logger::Log(
            LogLevel::Info,
            "Migrating database schema from legacy v1 to final normalized v2.");

        rhbm_gem::SQLiteWrapper::TransactionGuard transaction(database);
        const auto legacy_model_keys{ GetLegacyModelKeyList(database) };

        CreateFinalV2Tables(database);

        rhbm_gem::LegacyModelObjectReader legacy_reader{ &database };
        rhbm_gem::ModelObjectDAOv2 v2_dao{ &database };
        for (const auto & key_tag : legacy_model_keys)
        {
            auto model_object{ legacy_reader.Load(key_tag) };
            UpsertObjectCatalogRows(database, { key_tag }, "model");
            v2_dao.Save(model_object.get(), key_tag);
        }

        const auto legacy_map_keys{
            rhbm_gem::persistence::LookupManagedStoreDescriptor("map").list_keys(database)
        };
        UpsertObjectCatalogRows(database, legacy_map_keys, "map");

        int dropped_legacy_table_count{ 0 };
        for (const auto & table_name : BuildLegacyOwnedTableNames(legacy_model_keys))
        {
            if (HasTable(database, table_name))
            {
                database.Execute("DROP TABLE " + table_name + ";");
                ++dropped_legacy_table_count;
            }
        }

        DropTableIfExists(database, std::string(kManagedObjectMetadataTableName));
        SetSchemaVersion(database, rhbm_gem::DatabaseSchemaVersion::NormalizedV2);
        ValidateFinalV2Schema(database);
        Logger::Log(
            LogLevel::Info,
            "Legacy v1 migration completed. Migrated model count = " +
                std::to_string(legacy_model_keys.size()) +
                ", migrated map count = " + std::to_string(legacy_map_keys.size()) +
                ", dropped legacy table count = " +
                std::to_string(dropped_legacy_table_count) + ".");
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
        MigrateLegacyV1ToFinalV2(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    }
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        if (IsCatalogBasedV2Schema(*m_database))
        {
            ValidateFinalV2Schema(*m_database);
            return DatabaseSchemaVersion::NormalizedV2;
        }
        if (IsMetadataBasedLegacyV2Schema(*m_database) ||
            HasTable(*m_database, std::string(kManagedObjectMetadataTableName)))
        {
            RebuildManagedSchemaToFinalV2(*m_database, true);
            return DatabaseSchemaVersion::NormalizedV2;
        }
        throw std::runtime_error(
            "Database reports schema version 2 but does not match a supported normalized v2 shape.");
    }
    if (raw_version != 0)
    {
        throw std::runtime_error(
            "Unsupported database schema version: " + std::to_string(raw_version));
    }

    switch (InspectUnversionedSchemaState(*m_database))
    {
    case UnversionedSchemaState::Empty:
        CreateFinalV2Tables(*m_database);
        SetSchemaVersion(*m_database, DatabaseSchemaVersion::NormalizedV2);
        ValidateFinalV2Schema(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::LegacyV1:
        MigrateLegacyV1ToFinalV2(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::ManagedButUnversioned:
        RebuildManagedSchemaToFinalV2(
            *m_database,
            HasTable(*m_database, std::string(kManagedObjectMetadataTableName)));
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::MixedUnknown:
        throw std::runtime_error(
            "Database contains unmanaged or mixed schema tables and cannot be auto-upgraded safely.");
    }

    throw std::runtime_error("Unhandled database schema state.");
}

} // namespace rhbm_gem
