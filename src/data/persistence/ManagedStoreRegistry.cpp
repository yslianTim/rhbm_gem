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
} // namespace

namespace rhbm_gem::persistence {

const std::vector<ManagedStoreDescriptor> & GetAllManagedStoreDescriptors()
{
    static const std::vector<ManagedStoreDescriptor> k_descriptors{
        {
            "model",
            BuildModelManagedTableNameList(),
            EnsureModelSchema,
            ValidateModelSchema,
            ListModelKeys
        },
        {
            "map",
            BuildMapManagedTableNameList(),
            EnsureMapSchema,
            ValidateMapSchema,
            ListMapKeys
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
