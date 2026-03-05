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

namespace
{
    constexpr std::string_view kTableExistsSql =
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";

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

    bool IsIdentifierChar(const char ch)
    {
        return std::isalnum(static_cast<unsigned char>(ch)) != 0 || ch == '_';
    }

    std::string ReplaceIdentifier(
        const std::string & value,
        std::string_view identifier,
        const std::string & replacement)
    {
        if (identifier.empty())
        {
            return value;
        }

        std::string rewritten_sql;
        rewritten_sql.reserve(value.size() + replacement.size());

        std::size_t cursor{ 0 };
        while (cursor < value.size())
        {
            const auto position{ value.find(identifier, cursor) };
            if (position == std::string::npos)
            {
                rewritten_sql.append(value.substr(cursor));
                break;
            }

            rewritten_sql.append(value.substr(cursor, position - cursor));

            const bool has_left_boundary{
                position == 0 || !IsIdentifierChar(value[position - 1])
            };
            const auto end_position{ position + identifier.size() };
            const bool has_right_boundary{
                end_position >= value.size() || !IsIdentifierChar(value[end_position])
            };

            if (has_left_boundary && has_right_boundary)
            {
                rewritten_sql.append(replacement);
            }
            else
            {
                rewritten_sql.append(value.substr(position, identifier.size()));
            }

            cursor = end_position;
        }

        return rewritten_sql;
    }

    std::string BuildShadowModelSql(std::string_view sql, const std::string & suffix)
    {
        std::string rewritten_sql{ sql };
        std::vector<std::string> table_name_list;
        table_name_list.reserve(rhbm_gem::model_io::kModelCanonicalTableNames.size());
        for (const auto table_name : rhbm_gem::model_io::kModelCanonicalTableNames)
        {
            table_name_list.emplace_back(table_name);
        }
        std::sort(
            table_name_list.begin(),
            table_name_list.end(),
            [](const std::string & left, const std::string & right)
            {
                return left.size() > right.size();
            });

        for (const auto & table_name : table_name_list)
        {
            rewritten_sql = ReplaceIdentifier(
                std::move(rewritten_sql),
                table_name,
                table_name + suffix);
        }
        rewritten_sql = ReplaceIdentifier(
            std::move(rewritten_sql),
            "object_catalog",
            "object_catalog" + suffix);
        return rewritten_sql;
    }

    std::string BuildShadowMapSql(std::string_view sql, const std::string & suffix)
    {
        auto rewritten_sql{
            ReplaceIdentifier(std::string{ sql }, "map_list", "map_list" + suffix)
        };
        return ReplaceIdentifier(
            std::move(rewritten_sql),
            "object_catalog",
            "object_catalog" + suffix);
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
    }

    std::vector<std::string> ListModelKeys(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, "model_object"))
        {
            return {};
        }
        return QueryKeyList(database, "SELECT key_tag FROM model_object ORDER BY key_tag;");
    }

    void CreateModelShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        for (const auto create_sql : rhbm_gem::model_io::kCreateModelTableSqlList)
        {
            database.Execute(BuildShadowModelSql(create_sql, suffix));
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
    }

    std::vector<std::string> ListMapKeys(rhbm_gem::SQLiteWrapper & database)
    {
        if (!HasTable(database, "map_list"))
        {
            return {};
        }
        return QueryKeyList(database, "SELECT key_tag FROM map_list ORDER BY key_tag;");
    }

    void CreateMapShadowTables(rhbm_gem::SQLiteWrapper & database, const std::string & suffix)
    {
        database.Execute(BuildShadowMapSql(rhbm_gem::persistence::kCreateMapTableSql, suffix));
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
}

namespace rhbm_gem::persistence {

const std::vector<ManagedStoreDescriptor> & GetAllManagedStoreDescriptors()
{
    static const std::vector<ManagedStoreDescriptor> k_descriptors{
        {
            "model",
            EnsureModelSchema,
            ValidateModelSchema,
            ListModelKeys,
            CreateModelShadowTables,
            CopyModelIntoShadowTables,
            DropAndRenameModelShadowTables
        },
        {
            "map",
            EnsureMapSchema,
            ValidateMapSchema,
            ListMapKeys,
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
