#include <rhbm_gem/data/io/DataRepository.hpp>

#include "sqlite/ModelObjectStorage.hpp"
#include "sqlite/SQLiteWrapper.hpp"

#include <rhbm_gem/data/object/ModelObject.hpp>

#include <algorithm>
#include <array>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

using namespace std::literals;

constexpr int kCurrentSchemaVersion = 15;
constexpr std::string_view kUserSchemaObjectCountSql =
    "SELECT COUNT(*) FROM sqlite_master WHERE name NOT LIKE 'sqlite_%';";
constexpr std::string_view kTableNamesSql =
    "SELECT name FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%' "
    "ORDER BY name;";
constexpr std::string_view kPragmaUserVersionSql = "PRAGMA user_version;";

constexpr std::array<std::string_view, 10> kModelTableNames{
    "model_atom",
    "model_atom_group_potential",
    "model_atom_local_potential",
    "model_atom_posterior",
    "model_bond",
    "model_chain_map",
    "model_component",
    "model_component_atom",
    "model_component_bond",
    "model_object"
};

constexpr std::array<std::string_view, 9> kModelChildTableNames{
    "model_atom",
    "model_atom_group_potential",
    "model_atom_local_potential",
    "model_atom_posterior",
    "model_bond",
    "model_chain_map",
    "model_component",
    "model_component_atom",
    "model_component_bond"
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

int QuerySingleInt(rhbm_gem::SQLiteWrapper & database, const std::string & sql)
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

bool IsDatabaseEmpty(rhbm_gem::SQLiteWrapper & database)
{
    return QuerySingleInt(database, std::string(kUserSchemaObjectCountSql)) == 0;
}

std::vector<std::string> QueryTableNames(rhbm_gem::SQLiteWrapper & database)
{
    std::vector<std::string> table_names;
    database.Prepare(std::string(kTableNamesSql));
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
        table_names.push_back(database.GetColumn<std::string>(0));
    }
    return table_names;
}

std::vector<TableColumnInfo> QueryTableInfo(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view table_name)
{
    std::vector<TableColumnInfo> columns;
    database.Prepare("PRAGMA table_info(" + std::string(table_name) + ");");
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
        columns.push_back(TableColumnInfo{
            database.GetColumn<std::string>(1),
            database.GetColumn<int>(3),
            database.GetColumn<int>(5)
        });
    }
    return columns;
}

std::vector<ForeignKeyInfo> QueryForeignKeyList(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view table_name)
{
    std::vector<ForeignKeyInfo> foreign_keys;
    database.Prepare("PRAGMA foreign_key_list(" + std::string(table_name) + ");");
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
        foreign_keys.push_back(ForeignKeyInfo{
            database.GetColumn<std::string>(2),
            database.GetColumn<std::string>(3),
            database.GetColumn<std::string>(4),
            database.GetColumn<std::string>(6)
        });
    }
    return foreign_keys;
}

void ValidateColumns(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view table_name,
    std::initializer_list<std::string_view> expected_names)
{
    const auto columns{ QueryTableInfo(database, table_name) };
    if (columns.size() != expected_names.size())
    {
        throw std::runtime_error(
            "Schema v15 column mismatch for table: " + std::string(table_name));
    }

    auto column_iter{ columns.begin() };
    for (const auto expected_name : expected_names)
    {
        if (column_iter->name != expected_name)
        {
            throw std::runtime_error(
                "Schema v15 column mismatch for table: " + std::string(table_name));
        }
        ++column_iter;
    }
}

void ValidatePrimaryKey(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view table_name,
    std::initializer_list<std::string_view> expected_names)
{
    std::vector<std::pair<int, std::string>> primary_key_columns;
    for (const auto & column : QueryTableInfo(database, table_name))
    {
        if (column.primary_key_index > 0)
        {
            primary_key_columns.emplace_back(column.primary_key_index, column.name);
        }
    }
    std::sort(primary_key_columns.begin(), primary_key_columns.end());
    if (primary_key_columns.size() != expected_names.size())
    {
        throw std::runtime_error(
            "Schema v15 primary key mismatch for table: " + std::string(table_name));
    }

    auto primary_key_iter{ primary_key_columns.begin() };
    for (const auto expected_name : expected_names)
    {
        if (primary_key_iter->second != expected_name)
        {
            throw std::runtime_error(
                "Schema v15 primary key mismatch for table: " + std::string(table_name));
        }
        ++primary_key_iter;
    }
}

void ValidateSelectionColumn(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view table_name)
{
    for (const auto & column : QueryTableInfo(database, table_name))
    {
        if (column.name == "is_selected" && column.not_null != 0)
        {
            return;
        }
    }
    throw std::runtime_error(
        "Schema v15 requires a NOT NULL is_selected column on table: "
        + std::string(table_name));
}

void ValidateModelRootForeignKey(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view table_name)
{
    for (const auto & foreign_key : QueryForeignKeyList(database, table_name))
    {
        if (foreign_key.referenced_table == "model_object"
            && foreign_key.from_column == "key_tag"
            && foreign_key.to_column == "key_tag"
            && foreign_key.on_delete == "CASCADE")
        {
            return;
        }
    }
    throw std::runtime_error(
        "Schema v15 model root foreign key mismatch for table: "
        + std::string(table_name));
}

void ValidateCurrentSchema(rhbm_gem::SQLiteWrapper & database)
{
    if (QuerySingleInt(database, std::string(kUserSchemaObjectCountSql))
        != static_cast<int>(kModelTableNames.size()))
    {
        throw std::runtime_error("Schema v15 contains an unexpected schema object.");
    }

    const auto table_names{ QueryTableNames(database) };
    if (!std::equal(
            table_names.begin(),
            table_names.end(),
            kModelTableNames.begin(),
            kModelTableNames.end()))
    {
        throw std::runtime_error("Schema v15 contains an unexpected table set.");
    }

    ValidateColumns(database, "model_object", {
        "key_tag", "pdb_id", "emd_id", "map_resolution", "resolution_method",
        "standard_average_qscore", "reference_height", "reference_offset" });
    ValidateColumns(database, "model_chain_map", {
        "key_tag", "entity_id", "chain_ordinal", "chain_id" });
    ValidateColumns(database, "model_component", {
        "key_tag", "component_key", "id", "name", "type", "formula",
        "molecular_weight", "is_standard_monomer" });
    ValidateColumns(database, "model_component_atom", {
        "key_tag", "component_key", "atom_key", "atom_id", "element_type",
        "aromatic_atom_flag", "stereo_config" });
    ValidateColumns(database, "model_component_bond", {
        "key_tag", "component_key", "bond_key", "bond_id", "bond_type",
        "bond_order", "aromatic_atom_flag", "stereo_config" });
    ValidateColumns(database, "model_atom", {
        "key_tag", "serial_id", "sequence_id", "component_id", "atom_id",
        "chain_id", "indicator", "occupancy", "temperature", "element",
        "structure", "is_special_atom", "is_selected", "position_x",
        "position_y", "position_z", "component_key", "atom_key",
        "standard_qscore" });
    ValidateColumns(database, "model_bond", {
        "key_tag", "atom_serial_id_1", "atom_serial_id_2", "bond_key",
        "bond_type", "bond_order", "is_special_bond", "is_selected" });
    ValidateColumns(database, "model_atom_local_potential", {
        "key_tag", "serial_id", "raw_distance_and_map_value_list",
        "peeling_distance_and_map_value_list", "amplitude_estimate_ols_1st",
        "width_estimate_ols_1st", "intercept_estimate_ols_1st",
        "amplitude_estimate_mdpde_1st", "width_estimate_mdpde_1st",
        "intercept_estimate_mdpde_1st", "alpha_r_1st",
        "amplitude_estimate_ols_2nd", "width_estimate_ols_2nd",
        "intercept_estimate_ols_2nd", "amplitude_estimate_mdpde_2nd",
        "width_estimate_mdpde_2nd", "intercept_estimate_mdpde_2nd",
        "alpha_r_2nd", "amplitude_estimate_ols_3rd", "width_estimate_ols_3rd",
        "intercept_estimate_ols_3rd", "amplitude_estimate_mdpde_3rd",
        "width_estimate_mdpde_3rd", "intercept_estimate_mdpde_3rd",
        "alpha_r_3rd", "neighbor_count_for_peeling" });
    ValidateColumns(database, "model_atom_posterior", {
        "key_tag", "serial_id", "amplitude_estimate_posterior",
        "width_estimate_posterior", "intercept_estimate_posterior",
        "amplitude_variance_posterior", "width_variance_posterior",
        "intercept_variance_posterior", "outlier_tag", "statistical_distance" });
    ValidateColumns(database, "model_atom_group_potential", {
        "key_tag", "group_key", "amplitude_estimate_mean",
        "width_estimate_mean", "intercept_estimate_mean",
        "amplitude_estimate_mdpde", "width_estimate_mdpde",
        "intercept_estimate_mdpde", "amplitude_estimate_prior",
        "width_estimate_prior", "intercept_estimate_prior",
        "amplitude_variance_prior", "width_variance_prior",
        "intercept_variance_prior", "alpha_g" });

    ValidatePrimaryKey(database, "model_object", { "key_tag" });
    ValidatePrimaryKey(database, "model_chain_map", {
        "key_tag", "entity_id", "chain_ordinal" });
    ValidatePrimaryKey(database, "model_component", { "key_tag", "component_key" });
    ValidatePrimaryKey(database, "model_component_atom", {
        "key_tag", "component_key", "atom_key" });
    ValidatePrimaryKey(database, "model_component_bond", {
        "key_tag", "component_key", "bond_key" });
    ValidatePrimaryKey(database, "model_atom", { "key_tag", "serial_id" });
    ValidatePrimaryKey(database, "model_bond", {
        "key_tag", "atom_serial_id_1", "atom_serial_id_2" });
    ValidatePrimaryKey(database, "model_atom_local_potential", {
        "key_tag", "serial_id" });
    ValidatePrimaryKey(database, "model_atom_posterior", { "key_tag", "serial_id" });
    ValidatePrimaryKey(database, "model_atom_group_potential", {
        "key_tag", "group_key" });

    ValidateSelectionColumn(database, "model_atom");
    ValidateSelectionColumn(database, "model_bond");
    for (const auto table_name : kModelChildTableNames)
    {
        ValidateModelRootForeignKey(database, table_name);
    }
}

void EnsureCurrentSchema(rhbm_gem::SQLiteWrapper & database)
{
    const auto user_version{ QueryUserVersion(database) };
    if (user_version == 0 && IsDatabaseEmpty(database))
    {
        rhbm_gem::SQLiteWrapper::TransactionGuard transaction(database);
        rhbm_gem::model_storage::CreateTables(database);
        ValidateCurrentSchema(database);
        database.Execute("PRAGMA user_version = 15;");
        return;
    }
    if (user_version == kCurrentSchemaVersion)
    {
        ValidateCurrentSchema(database);
        return;
    }
    throw std::runtime_error(
        "Unsupported SQLite schema: expected an empty version-0 database or schema v15.");
}

} // namespace

namespace rhbm_gem {

DataRepository::DataRepository(const std::filesystem::path & database_path) :
    m_database{ nullptr }
{
    const auto resolved_database_path{
        database_path.empty() ? std::filesystem::path{ "database.sqlite" } : database_path
    };

    const auto parent_path{ resolved_database_path.parent_path() };
    if (!parent_path.empty())
    {
        std::error_code error_code;
        std::filesystem::create_directories(parent_path, error_code);
        if (error_code)
        {
            throw std::runtime_error(
                "Failed to create database parent directory '" + parent_path.string()
                + "': " + error_code.message());
        }
    }

    m_database = std::make_unique<SQLiteWrapper>(resolved_database_path);
    EnsureCurrentSchema(*m_database);
}

DataRepository::~DataRepository() = default;

std::unique_ptr<ModelObject> DataRepository::LoadModel(const std::string & key_tag) const
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    return model_storage::Load(*m_database, key_tag);
}

void DataRepository::SaveModel(
    const ModelObject & model_object,
    const std::string & key_tag) const
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    model_storage::Save(*m_database, model_object, key_tag);
}

} // namespace rhbm_gem
