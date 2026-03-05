#include "ManagedStoreRegistry.hpp"

#include "MapObjectDAO.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"
#include "../model_io/ModelSchemaSql.hpp"
#include "MapStoreSql.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{
    constexpr std::string_view kTableExistsSql =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";

    struct TableColumnInfo
    {
        std::string name;
        int primary_key_index;
    };

    struct ForeignKeyInfo
    {
        std::string referenced_table;
        std::string from_column;
        std::string to_column;
        std::string on_delete;
    };

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

    std::string ToUpperCopy(std::string value)
    {
        std::transform(
            value.begin(),
            value.end(),
            value.begin(),
            [](const unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
        return value;
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

    std::string BuildCreateModelShadowTableSql(
        std::string_view canonical_table_name,
        const std::string & suffix)
    {
        const auto table_name{ std::string(canonical_table_name) + suffix };
        const auto object_catalog_table_name{ "object_catalog" + suffix };
        const auto model_root_table_name{ "model_object" + suffix };

        if (canonical_table_name == "model_object")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT PRIMARY KEY REFERENCES " + object_catalog_table_name + "(key_tag) ON DELETE CASCADE, "
                "atom_size INTEGER, "
                "pdb_id TEXT, "
                "emd_id TEXT, "
                "map_resolution DOUBLE, "
                "resolution_method TEXT"
                ");";
        }
        if (canonical_table_name == "model_chain_map")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "entity_id TEXT, "
                "chain_ordinal INTEGER, "
                "chain_id TEXT, "
                "PRIMARY KEY (key_tag, entity_id, chain_ordinal), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_component")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "component_key INTEGER, "
                "id TEXT, "
                "name TEXT, "
                "type TEXT, "
                "formula TEXT, "
                "molecular_weight DOUBLE, "
                "is_standard_monomer INTEGER, "
                "PRIMARY KEY (key_tag, component_key), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_component_atom")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "component_key INTEGER, "
                "atom_key INTEGER, "
                "atom_id TEXT, "
                "element_type INTEGER, "
                "aromatic_atom_flag INTEGER, "
                "stereo_config INTEGER, "
                "PRIMARY KEY (key_tag, component_key, atom_key), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_component_bond")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "component_key INTEGER, "
                "bond_key INTEGER, "
                "bond_id TEXT, "
                "bond_type INTEGER, "
                "bond_order INTEGER, "
                "aromatic_atom_flag INTEGER, "
                "stereo_config INTEGER, "
                "PRIMARY KEY (key_tag, component_key, bond_key), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_atom")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "serial_id INTEGER, "
                "sequence_id INTEGER, "
                "component_id TEXT, "
                "atom_id TEXT, "
                "chain_id TEXT, "
                "indicator TEXT, "
                "occupancy DOUBLE, "
                "temperature DOUBLE, "
                "element INTEGER, "
                "structure INTEGER, "
                "is_special_atom INTEGER, "
                "position_x DOUBLE, "
                "position_y DOUBLE, "
                "position_z DOUBLE, "
                "component_key INTEGER, "
                "atom_key INTEGER, "
                "PRIMARY KEY (key_tag, serial_id), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_bond")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "atom_serial_id_1 INTEGER, "
                "atom_serial_id_2 INTEGER, "
                "bond_key INTEGER, "
                "bond_type INTEGER, "
                "bond_order INTEGER, "
                "is_special_bond INTEGER, "
                "PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_atom_local_potential")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "serial_id INTEGER, "
                "sampling_size INTEGER, "
                "distance_and_map_value_list BLOB, "
                "amplitude_estimate_ols DOUBLE, "
                "width_estimate_ols DOUBLE, "
                "amplitude_estimate_mdpde DOUBLE, "
                "width_estimate_mdpde DOUBLE, "
                "alpha_r DOUBLE, "
                "PRIMARY KEY (key_tag, serial_id), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_bond_local_potential")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "atom_serial_id_1 INTEGER, "
                "atom_serial_id_2 INTEGER, "
                "sampling_size INTEGER, "
                "distance_and_map_value_list BLOB, "
                "amplitude_estimate_ols DOUBLE, "
                "width_estimate_ols DOUBLE, "
                "amplitude_estimate_mdpde DOUBLE, "
                "width_estimate_mdpde DOUBLE, "
                "alpha_r DOUBLE, "
                "PRIMARY KEY (key_tag, atom_serial_id_1, atom_serial_id_2), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_atom_posterior")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "class_key TEXT, "
                "serial_id INTEGER, "
                "amplitude_estimate_posterior DOUBLE, "
                "width_estimate_posterior DOUBLE, "
                "amplitude_variance_posterior DOUBLE, "
                "width_variance_posterior DOUBLE, "
                "outlier_tag INTEGER, "
                "statistical_distance DOUBLE, "
                "PRIMARY KEY (key_tag, class_key, serial_id), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_bond_posterior")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "class_key TEXT, "
                "atom_serial_id_1 INTEGER, "
                "atom_serial_id_2 INTEGER, "
                "amplitude_estimate_posterior DOUBLE, "
                "width_estimate_posterior DOUBLE, "
                "amplitude_variance_posterior DOUBLE, "
                "width_variance_posterior DOUBLE, "
                "outlier_tag INTEGER, "
                "statistical_distance DOUBLE, "
                "PRIMARY KEY (key_tag, class_key, atom_serial_id_1, atom_serial_id_2), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_atom_group_potential")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "class_key TEXT, "
                "group_key INTEGER, "
                "member_size INTEGER, "
                "amplitude_estimate_mean DOUBLE, "
                "width_estimate_mean DOUBLE, "
                "amplitude_estimate_mdpde DOUBLE, "
                "width_estimate_mdpde DOUBLE, "
                "amplitude_estimate_prior DOUBLE, "
                "width_estimate_prior DOUBLE, "
                "amplitude_variance_prior DOUBLE, "
                "width_variance_prior DOUBLE, "
                "alpha_g DOUBLE, "
                "PRIMARY KEY (key_tag, class_key, group_key), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }
        if (canonical_table_name == "model_bond_group_potential")
        {
            return
                "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                "key_tag TEXT, "
                "class_key TEXT, "
                "group_key INTEGER, "
                "member_size INTEGER, "
                "amplitude_estimate_mean DOUBLE, "
                "width_estimate_mean DOUBLE, "
                "amplitude_estimate_mdpde DOUBLE, "
                "width_estimate_mdpde DOUBLE, "
                "amplitude_estimate_prior DOUBLE, "
                "width_estimate_prior DOUBLE, "
                "amplitude_variance_prior DOUBLE, "
                "width_variance_prior DOUBLE, "
                "alpha_g DOUBLE, "
                "PRIMARY KEY (key_tag, class_key, group_key), "
                "FOREIGN KEY(key_tag) REFERENCES " + model_root_table_name + "(key_tag) ON DELETE CASCADE"
                ");";
        }

        throw std::runtime_error(
            "Unsupported model table for shadow schema build: " + std::string(canonical_table_name));
    }

    std::string BuildCreateMapShadowTableSql(const std::string & suffix)
    {
        const auto table_name{ "map_list" + suffix };
        const auto object_catalog_table_name{ "object_catalog" + suffix };
        return
            "CREATE TABLE IF NOT EXISTS " + table_name + " ("
            "key_tag TEXT PRIMARY KEY REFERENCES " + object_catalog_table_name + "(key_tag) ON DELETE CASCADE, "
            "grid_size_x INTEGER, "
            "grid_size_y INTEGER, "
            "grid_size_z INTEGER, "
            "grid_spacing_x DOUBLE, "
            "grid_spacing_y DOUBLE, "
            "grid_spacing_z DOUBLE, "
            "origin_x DOUBLE, "
            "origin_y DOUBLE, "
            "origin_z DOUBLE, "
            "map_value_array BLOB"
            ");";
    }

    void EnsureModelSchema(rhbm_gem::SQLiteWrapper & database)
    {
        rhbm_gem::ModelObjectDAOv2::EnsureSchema(database);
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

    std::vector<std::string> ListModelKeysWithSuffix(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & suffix)
    {
        const auto table_name{ "model_object" + suffix };
        if (!HasTable(database, table_name))
        {
            return {};
        }
        return QueryKeyList(database, "SELECT key_tag FROM " + table_name + " ORDER BY key_tag;");
    }

    void CreateModelShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            database.Execute(BuildCreateModelShadowTableSql(table_name, suffix));
        }
    }

    void CopyModelIntoShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            if (!HasTable(database, std::string(table_name)))
            {
                continue;
            }
            database.Execute(
                "INSERT INTO " + std::string(table_name) + suffix +
                " SELECT * FROM " + std::string(table_name) + ";");
        }
    }

    void DropAndRenameModelShadowTables(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & suffix)
    {
        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            if (HasTable(database, std::string(table_name)))
            {
                database.Execute("DROP TABLE " + std::string(table_name) + ";");
            }
        }
        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            const auto shadow_table_name{ std::string(table_name) + suffix };
            if (!HasTable(database, shadow_table_name))
            {
                continue;
            }
            database.Execute(
                "ALTER TABLE " + shadow_table_name + " RENAME TO " + std::string(table_name) + ";");
        }
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

    std::vector<std::string> ListMapKeysWithSuffix(
        rhbm_gem::SQLiteWrapper & database,
        const std::string & suffix)
    {
        const auto table_name{ "map_list" + suffix };
        if (!HasTable(database, table_name))
        {
            return {};
        }
        return QueryKeyList(database, "SELECT key_tag FROM " + table_name + " ORDER BY key_tag;");
    }

    void CreateMapShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        database.Execute(BuildCreateMapShadowTableSql(suffix));
    }

    void CopyMapIntoShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        if (!HasTable(database, "map_list"))
        {
            return;
        }
        database.Execute("INSERT INTO map_list" + suffix + " SELECT * FROM map_list;");
    }

    void DropAndRenameMapShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        if (HasTable(database, "map_list"))
        {
            database.Execute("DROP TABLE map_list;");
        }
        if (HasTable(database, "map_list" + suffix))
        {
            database.Execute("ALTER TABLE map_list" + suffix + " RENAME TO map_list;");
        }
    }

    std::vector<std::string_view> BuildModelManagedTableNameList()
    {
        return {
            rhbm_gem::model_io::kModelCanonicalTableNames.begin(),
            rhbm_gem::model_io::kModelCanonicalTableNames.end()
        };
    }

    std::vector<std::string_view> BuildMapManagedTableNameList()
    {
        return {
            rhbm_gem::persistence::kMapCanonicalTableNames.begin(),
            rhbm_gem::persistence::kMapCanonicalTableNames.end()
        };
    }
}

namespace rhbm_gem::persistence {

const std::vector<ManagedStoreDescriptor> & GetAllManagedStoreDescriptors()
{
    static const std::vector<ManagedStoreDescriptor> k_descriptors{
        {
            "model",
            BuildModelManagedTableNameList(),
            EnsureModelSchema,
            ValidateModelSchema,
            ListModelKeys,
            ListModelKeysWithSuffix,
            CreateModelShadowTables,
            CopyModelIntoShadowTables,
            DropAndRenameModelShadowTables
        },
        {
            "map",
            BuildMapManagedTableNameList(),
            EnsureMapSchema,
            ValidateMapSchema,
            ListMapKeys,
            ListMapKeysWithSuffix,
            CreateMapShadowTables,
            CopyMapIntoShadowTables,
            DropAndRenameMapShadowTables
        }
    };
    return k_descriptors;
}

const ManagedStoreDescriptor & LookupManagedStoreDescriptor(std::string_view object_type)
{
    for (const auto & descriptor : GetAllManagedStoreDescriptors())
    {
        if (descriptor.object_type == object_type)
        {
            return descriptor;
        }
    }
    throw std::runtime_error("Unsupported managed store descriptor: " + std::string(object_type));
}

} // namespace rhbm_gem::persistence
