#include "internal/migration/DatabaseSchemaManager.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "internal/io/sqlite/MapObjectDAO.hpp"
#include "internal/io/sqlite/ModelObjectDAOSqlite.hpp"
#include "internal/io/sqlite/SQLiteWrapper.hpp"
#include "io/sqlite/MapStoreSql.hpp"
#include "io/sqlite/ModelSchemaSql.hpp"

#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
#include "internal/migration/LegacyModelObjectReader.hpp"
#endif

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
        NonEmptyNonLegacy
    };

    struct TableColumnInfo
    {
        std::string name;
        int not_null;
        int primary_key_index;
    };

    struct ForeignKeyInfo
    {
        std::string referenced_table;
        std::string from_column;
        std::string to_column;
        std::string on_delete;
    };

    struct LegacyMapRow
    {
        std::string key_tag;
        std::array<int, 3> grid_size;
        std::array<double, 3> grid_spacing;
        std::array<double, 3> origin;
        std::vector<float> map_value_array;
    };

    constexpr std::string_view kTableExistsSql =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";
    constexpr std::string_view kDatabaseEmptySql =
        "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";
    constexpr std::string_view kPragmaUserVersionSql = "PRAGMA user_version;";
#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
    constexpr std::string_view kLegacyModelKeyListSql = "SELECT key_tag FROM model_list;";
#endif

    constexpr std::string_view kFinalCatalogTableName = "object_catalog";
    constexpr std::string_view kManagedObjectMetadataTableName = "object_metadata";

    std::string ToUpperCopy(std::string value)
    {
        std::transform(
            value.begin(),
            value.end(),
            value.begin(),
            [](const unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
        return value;
    }

#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
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
#endif

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

#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
    void DropTableIfExists(rhbm_gem::SQLiteWrapper & database, const std::string & table_name)
    {
        if (HasTable(database, table_name))
        {
            database.Execute("DROP TABLE " + table_name + ";");
        }
    }
#endif

    bool IsDatabaseEmpty(rhbm_gem::SQLiteWrapper & database)
    {
        return QuerySingleInt(database, std::string(kDatabaseEmptySql)) == 0;
    }

    bool HasLegacyModelSchema(rhbm_gem::SQLiteWrapper & database)
    {
        return HasTable(database, "model_list");
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
        return UnversionedSchemaState::NonEmptyNonLegacy;
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

#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
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
#endif

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

    std::vector<TableColumnInfo> QueryTableInfo(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & table_name)
    {
        std::vector<TableColumnInfo> column_info_list;
        database.Prepare("PRAGMA table_info(" + table_name + ");");
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
            column_info_list.push_back(
                TableColumnInfo{
                    database.GetColumn<std::string>(1),
                    database.GetColumn<int>(3),
                    database.GetColumn<int>(5)
                });
        }
        return column_info_list;
    }

    std::vector<ForeignKeyInfo> QueryForeignKeyList(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & table_name)
    {
        std::vector<ForeignKeyInfo> foreign_key_info_list;
        database.Prepare("PRAGMA foreign_key_list(" + table_name + ");");
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
            foreign_key_info_list.push_back(
                ForeignKeyInfo{
                    database.GetColumn<std::string>(2),
                    database.GetColumn<std::string>(3),
                    database.GetColumn<std::string>(4),
                    database.GetColumn<std::string>(6)
                });
        }
        return foreign_key_info_list;
    }

    std::vector<std::string> QueryKeyList(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & sql)
    {
        std::vector<std::string> key_list;
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
            key_list.push_back(database.GetColumn<std::string>(0));
        }
        return key_list;
    }

    void ValidatePrimaryKeyShape(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & table_name,
        const std::vector<std::string_view> & expected_primary_key_columns)
    {
        if (!HasTable(database, table_name))
        {
            throw std::runtime_error("Cannot validate missing table: " + table_name);
        }

        std::vector<std::pair<int, std::string>> actual_primary_key_columns;
        for (const auto & column_info : QueryTableInfo(database, table_name))
        {
            if (column_info.primary_key_index > 0)
            {
                actual_primary_key_columns.emplace_back(
                    column_info.primary_key_index,
                    column_info.name);
            }
        }
        std::sort(
            actual_primary_key_columns.begin(),
            actual_primary_key_columns.end(),
            [](const auto & left, const auto & right)
            {
                return left.first < right.first;
            });

        if (actual_primary_key_columns.size() != expected_primary_key_columns.size())
        {
            throw std::runtime_error(
                "Normalized v2 schema primary key shape mismatch for table: " + table_name);
        }

        for (size_t i = 0; i < expected_primary_key_columns.size(); i++)
        {
            if (actual_primary_key_columns[i].second != expected_primary_key_columns[i])
            {
                throw std::runtime_error(
                    "Normalized v2 schema primary key shape mismatch for table: " + table_name);
            }
        }
    }

    void ValidateForeignKey(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & table_name,
        std::string_view from_column,
        std::string_view referenced_table,
        std::string_view referenced_column,
        std::string_view on_delete)
    {
        const auto expected_on_delete{ ToUpperCopy(std::string(on_delete)) };
        for (const auto & foreign_key_info : QueryForeignKeyList(database, table_name))
        {
            if (foreign_key_info.from_column != from_column) continue;
            if (foreign_key_info.referenced_table != referenced_table) continue;
            if (foreign_key_info.to_column != referenced_column) continue;
            if (ToUpperCopy(foreign_key_info.on_delete) != expected_on_delete) continue;
            return;
        }
        throw std::runtime_error(
            "Normalized v2 schema is missing required foreign key on table: " + table_name);
    }

    template <std::size_t N>
    void ValidateRequiredTables(
        rhbm_gem::SQLiteWrapper & database,
        const std::array<std::string_view, N> & table_names,
        const std::string & object_type)
    {
        for (const auto table_name : table_names)
        {
            if (!HasTable(database, std::string(table_name)))
            {
                throw std::runtime_error(
                    "Normalized v2 schema is missing required " + object_type +
                    " table: " + std::string(table_name));
            }
        }
    }

    void EnsureModelSchema(rhbm_gem::SQLiteWrapper & database)
    {
        rhbm_gem::ModelObjectDAOSqlite::EnsureSchema(database);
    }

    void ValidateModelSchema(rhbm_gem::SQLiteWrapper & database)
    {
        ValidateRequiredTables(database, rhbm_gem::model_io::kModelCanonicalTableNames, "model");

        ValidatePrimaryKeyShape(database, "model_object", {"key_tag"});
        ValidatePrimaryKeyShape(database, "model_chain_map", {"key_tag", "entity_id", "chain_ordinal"});
        ValidatePrimaryKeyShape(database, "model_component", {"key_tag", "component_key"});
        ValidatePrimaryKeyShape(database, "model_component_atom", {"key_tag", "component_key", "atom_key"});
        ValidatePrimaryKeyShape(database, "model_component_bond", {"key_tag", "component_key", "bond_key"});
        ValidatePrimaryKeyShape(database, "model_atom", {"key_tag", "serial_id"});
        ValidatePrimaryKeyShape(database, "model_bond", {"key_tag", "atom_serial_id_1", "atom_serial_id_2"});
        ValidatePrimaryKeyShape(database, "model_atom_local_potential", {"key_tag", "serial_id"});
        ValidatePrimaryKeyShape(
            database,
            "model_bond_local_potential",
            {"key_tag", "atom_serial_id_1", "atom_serial_id_2"});
        ValidatePrimaryKeyShape(database, "model_atom_posterior", {"key_tag", "class_key", "serial_id"});
        ValidatePrimaryKeyShape(
            database,
            "model_bond_posterior",
            {"key_tag", "class_key", "atom_serial_id_1", "atom_serial_id_2"});
        ValidatePrimaryKeyShape(database, "model_atom_group_potential", {"key_tag", "class_key", "group_key"});
        ValidatePrimaryKeyShape(database, "model_bond_group_potential", {"key_tag", "class_key", "group_key"});

        ValidateForeignKey(database, "model_object", "key_tag", "object_catalog", "key_tag", "CASCADE");

        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            if (table_name == std::string_view("model_object"))
            {
                continue;
            }
            ValidateForeignKey(
                database,
                std::string(table_name),
                "key_tag",
                "model_object",
                "key_tag",
                "CASCADE");
        }
    }

    std::vector<std::string> ListModelKeys(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, "model_object"))
        {
            return {};
        }
        return QueryKeyList(database, "SELECT key_tag FROM model_object ORDER BY key_tag;");
    }

    void EnsureMapSchema(rhbm_gem::SQLiteWrapper & database)
    {
        rhbm_gem::MapObjectDAO::EnsureSchema(database);
    }

    void ValidateMapSchema(rhbm_gem::SQLiteWrapper & database)
    {
        ValidateRequiredTables(database, rhbm_gem::persistence::kMapCanonicalTableNames, "map");
        ValidatePrimaryKeyShape(database, "map_list", {"key_tag"});
        ValidateForeignKey(database, "map_list", "key_tag", "object_catalog", "key_tag", "CASCADE");
    }

    std::vector<std::string> ListMapKeys(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, "map_list"))
        {
            return {};
        }
        return QueryKeyList(database, "SELECT key_tag FROM map_list ORDER BY key_tag;");
    }

    void ValidateObjectCatalogShape(rhbm_gem::SQLiteWrapper & database)
    {
        const auto columns{ QueryTableInfo(database, std::string(kFinalCatalogTableName)) };

        const auto key_tag_iter{
            std::find_if(
                columns.begin(),
                columns.end(),
                [](const TableColumnInfo & column_info)
                {
                    return column_info.name == "key_tag";
                })
        };
        if (key_tag_iter == columns.end() || key_tag_iter->primary_key_index <= 0)
        {
            throw std::runtime_error(
                "Normalized v2 object_catalog must contain primary key column: key_tag");
        }

        const auto object_type_iter{
            std::find_if(
                columns.begin(),
                columns.end(),
                [](const TableColumnInfo & column_info)
                {
                    return column_info.name == "object_type";
                })
        };
        if (object_type_iter == columns.end())
        {
            throw std::runtime_error(
                "Normalized v2 object_catalog must contain column: object_type");
        }
        if (object_type_iter->not_null == 0)
        {
            throw std::runtime_error(
                "Normalized v2 object_catalog column object_type must be NOT NULL.");
        }

        const auto unexpected_object_type_count{
            QuerySingleInt(
                database,
                "SELECT COUNT(*) FROM object_catalog "
                "WHERE object_type IS NULL OR object_type NOT IN ('model', 'map');")
        };
        if (unexpected_object_type_count != 0)
        {
            throw std::runtime_error(
                "Normalized v2 object_catalog contains unsupported object_type values.");
        }
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

        ValidateObjectCatalogShape(database);

        const auto validate_catalog_consistency =
            [&database](std::string_view object_type, const std::vector<std::string> & payload_keys)
        {
            const auto catalog_keys{ QueryCatalogKeys(database, std::string(object_type)) };
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
                        std::string(object_type) +
                        " key: " + key_tag);
                }
            }

            for (const auto & key_tag : catalog_keys)
            {
                if (payload_key_set.find(key_tag) == payload_key_set.end())
                {
                    throw std::runtime_error(
                        "Normalized v2 catalog contains ghost " +
                        std::string(object_type) +
                        " key without payload: " + key_tag);
                }
            }
        };

        ValidateModelSchema(database);
        ValidateMapSchema(database);
        validate_catalog_consistency("model", ListModelKeys(database));
        validate_catalog_consistency("map", ListMapKeys(database));
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

#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
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
        return QuerySingleStringColumn(database, std::string(kLegacyModelKeyListSql));
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
#endif // RHBM_GEM_LEGACY_V1_SUPPORT

    void CreateFinalV2Tables(rhbm_gem::SQLiteWrapper & database)
    {
        CreateObjectCatalogTable(database);
        EnsureModelSchema(database);
        EnsureMapSchema(database);
    }

#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
    std::vector<LegacyMapRow> LoadLegacyMapRows(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, std::string(rhbm_gem::persistence::kMapTableName)))
        {
            return {};
        }

        std::vector<LegacyMapRow> rows;
        database.Prepare(
            "SELECT key_tag, "
            "grid_size_x, grid_size_y, grid_size_z, "
            "grid_spacing_x, grid_spacing_y, grid_spacing_z, "
            "origin_x, origin_y, origin_z, map_value_array "
            "FROM map_list ORDER BY key_tag;");
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

            LegacyMapRow row;
            row.key_tag = database.GetColumn<std::string>(0);
            row.grid_size = {
                database.GetColumn<int>(1),
                database.GetColumn<int>(2),
                database.GetColumn<int>(3)
            };
            row.grid_spacing = {
                database.GetColumn<double>(4),
                database.GetColumn<double>(5),
                database.GetColumn<double>(6)
            };
            row.origin = {
                database.GetColumn<double>(7),
                database.GetColumn<double>(8),
                database.GetColumn<double>(9)
            };
            row.map_value_array = database.GetColumn<std::vector<float>>(10);
            rows.emplace_back(std::move(row));
        }
        return rows;
    }

    std::vector<std::string> BuildLegacyMapKeyList(const std::vector<LegacyMapRow> & legacy_map_rows)
    {
        std::vector<std::string> key_list;
        key_list.reserve(legacy_map_rows.size());
        for (const auto & row : legacy_map_rows)
        {
            key_list.push_back(row.key_tag);
        }
        return key_list;
    }

    void ValidateNoLegacyRootKeyCollision(
        const std::vector<std::string> & legacy_model_keys,
        const std::vector<LegacyMapRow> & legacy_map_rows)
    {
        std::unordered_set<std::string> model_key_set{
            legacy_model_keys.begin(), legacy_model_keys.end()
        };
        for (const auto & map_row : legacy_map_rows)
        {
            if (model_key_set.find(map_row.key_tag) != model_key_set.end())
            {
                throw std::runtime_error(
                    "Legacy migration detected conflicting key_tag in both model_list and map_list: " +
                    map_row.key_tag);
            }
        }
    }

    void InsertLegacyMapRowsIntoFinalV2(
        rhbm_gem::SQLiteWrapper & database,
        const std::vector<LegacyMapRow> & legacy_map_rows)
    {
        if (legacy_map_rows.empty())
        {
            return;
        }

        database.Prepare(std::string(rhbm_gem::persistence::kInsertMapSql));
        rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
        for (const auto & row : legacy_map_rows)
        {
            database.Bind<std::string>(1, row.key_tag);
            database.Bind<int>(2, row.grid_size[0]);
            database.Bind<int>(3, row.grid_size[1]);
            database.Bind<int>(4, row.grid_size[2]);
            database.Bind<double>(5, row.grid_spacing[0]);
            database.Bind<double>(6, row.grid_spacing[1]);
            database.Bind<double>(7, row.grid_spacing[2]);
            database.Bind<double>(8, row.origin[0]);
            database.Bind<double>(9, row.origin[1]);
            database.Bind<double>(10, row.origin[2]);
            database.Bind<std::vector<float>>(11, row.map_value_array);
            database.StepOnce();
            database.Reset();
        }
    }

    void MigrateLegacyV1ToFinalV2(rhbm_gem::SQLiteWrapper & database)
    {
        Logger::Log(
            LogLevel::Info,
            "Migrating database schema from legacy v1 to final normalized v2.");

        rhbm_gem::SQLiteWrapper::TransactionGuard transaction(database);
        const auto legacy_model_keys{ GetLegacyModelKeyList(database) };
        const auto legacy_map_rows{ LoadLegacyMapRows(database) };
        ValidateNoLegacyRootKeyCollision(legacy_model_keys, legacy_map_rows);

        DropTableIfExists(database, std::string(rhbm_gem::persistence::kMapTableName));
        CreateFinalV2Tables(database);

        rhbm_gem::LegacyModelObjectReader legacy_reader{ &database };
        rhbm_gem::ModelObjectDAOSqlite v2_dao{ &database };
        for (const auto & key_tag : legacy_model_keys)
        {
            auto model_object{ legacy_reader.Load(key_tag) };
            UpsertObjectCatalogRows(database, { key_tag }, "model");
            v2_dao.Save(*model_object, key_tag);
        }

        const auto legacy_map_keys{ BuildLegacyMapKeyList(legacy_map_rows) };
        UpsertObjectCatalogRows(database, legacy_map_keys, "map");
        InsertLegacyMapRowsIntoFinalV2(database, legacy_map_rows);

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
                ", migrated map count = " + std::to_string(legacy_map_rows.size()) +
                ", dropped legacy table count = " +
                std::to_string(dropped_legacy_table_count) + ".");
    }
#endif // RHBM_GEM_LEGACY_V1_SUPPORT
} // namespace

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
#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
        MigrateLegacyV1ToFinalV2(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
#else
        throw std::runtime_error(
            "Detected legacy v1 database schema, but support is disabled. "
            "Rebuild with -DRHBM_GEM_LEGACY_V1_SUPPORT=ON to enable migration.");
#endif
    }
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        ValidateFinalV2Schema(*m_database);
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
        CreateFinalV2Tables(*m_database);
        SetSchemaVersion(*m_database, DatabaseSchemaVersion::NormalizedV2);
        ValidateFinalV2Schema(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
    case UnversionedSchemaState::LegacyV1:
#ifdef RHBM_GEM_LEGACY_V1_SUPPORT
        MigrateLegacyV1ToFinalV2(*m_database);
        return DatabaseSchemaVersion::NormalizedV2;
#else
        throw std::runtime_error(
            "Detected unversioned legacy v1 schema, but support is disabled. "
            "Rebuild with -DRHBM_GEM_LEGACY_V1_SUPPORT=ON to enable migration.");
#endif
    case UnversionedSchemaState::NonEmptyNonLegacy:
        throw std::runtime_error(
            "Unversioned non-empty database is unsupported unless it is a legacy v1 schema.");
    }

    throw std::runtime_error("Unhandled database schema state.");
}

} // namespace rhbm_gem
