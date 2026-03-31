#include "SQLitePersistence.hpp"

#include "../detail/TopLevelObjectRouting.hpp"
#include "ModelObjectStorage.hpp"
#include "SQLiteWrapper.hpp"
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_set>
#include <vector>

namespace {

using namespace std::literals;

constexpr int kCurrentSchemaVersion = 2;
constexpr std::string_view kCatalogTableName = "object_catalog";
constexpr std::string_view kMapTableName = "map_list";
constexpr std::string_view kUnsupportedMetadataTableName = "object_metadata";
constexpr std::string_view kTableExistsSql =
    "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;";
constexpr std::string_view kDatabaseEmptySql =
    "SELECT COUNT(*) FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%';";
constexpr std::string_view kPragmaUserVersionSql = "PRAGMA user_version;";

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

constexpr std::array<std::string_view, 13> kModelCanonicalTableNames{
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
    "model_bond_group_potential"};

constexpr std::array<std::string_view, 1> kMapCanonicalTableNames{
    kMapTableName};

constexpr std::string_view kCreateCatalogTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS object_catalog (
        key_tag TEXT PRIMARY KEY,
        object_type TEXT NOT NULL,
        CHECK (object_type IN ('model', 'map'))
    )
)sql";

constexpr std::string_view kCreateMapTableSql = R"sql(
    CREATE TABLE IF NOT EXISTS map_list (
        key_tag TEXT PRIMARY KEY REFERENCES object_catalog(key_tag) ON DELETE CASCADE,
        grid_size_x INTEGER,
        grid_size_y INTEGER,
        grid_size_z INTEGER,
        grid_spacing_x DOUBLE,
        grid_spacing_y DOUBLE,
        grid_spacing_z DOUBLE,
        origin_x DOUBLE,
        origin_y DOUBLE,
        origin_z DOUBLE,
        map_value_array BLOB
    )
)sql";

constexpr std::string_view kUpsertCatalogSql = R"sql(
    INSERT INTO object_catalog(key_tag, object_type) VALUES (?, ?)
    ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type
)sql";

constexpr std::string_view kLoadCatalogTypeSql =
    "SELECT object_type FROM object_catalog WHERE key_tag = ? LIMIT 1;";

constexpr std::string_view kInsertMapSql = R"sql(
    INSERT INTO map_list (
        key_tag,
        grid_size_x, grid_size_y, grid_size_z,
        grid_spacing_x, grid_spacing_y, grid_spacing_z,
        origin_x, origin_y, origin_z,
        map_value_array
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(key_tag) DO UPDATE SET
        key_tag = excluded.key_tag,
        grid_size_x = excluded.grid_size_x,
        grid_size_y = excluded.grid_size_y,
        grid_size_z = excluded.grid_size_z,
        grid_spacing_x = excluded.grid_spacing_x,
        grid_spacing_y = excluded.grid_spacing_y,
        grid_spacing_z = excluded.grid_spacing_z,
        origin_x = excluded.origin_x,
        origin_y = excluded.origin_y,
        origin_z = excluded.origin_z,
        map_value_array = excluded.map_value_array
)sql";

constexpr std::string_view kSelectMapSql = R"sql(
    SELECT
        key_tag,
        grid_size_x, grid_size_y, grid_size_z,
        grid_spacing_x, grid_spacing_y, grid_spacing_z,
        origin_x, origin_y, origin_z,
        map_value_array
    FROM map_list
    WHERE key_tag = ? LIMIT 1
)sql";

std::string ToUpperCopy(std::string value)
{
    std::transform(
        value.begin(),
        value.end(),
        value.begin(),
        [](const unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    return value;
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

int QueryUserVersion(rhbm_gem::SQLiteWrapper & database)
{
    return QuerySingleInt(database, std::string(kPragmaUserVersionSql));
}

bool IsDatabaseEmpty(rhbm_gem::SQLiteWrapper & database)
{
    return QuerySingleInt(database, std::string(kDatabaseEmptySql)) == 0;
}

void SetSchemaVersion(rhbm_gem::SQLiteWrapper & database)
{
    database.Execute("PRAGMA user_version = " + std::to_string(kCurrentSchemaVersion) + ";");
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
        column_info_list.push_back(TableColumnInfo{
            database.GetColumn<std::string>(1),
            database.GetColumn<int>(3),
            database.GetColumn<int>(5)});
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
        foreign_key_info_list.push_back(ForeignKeyInfo{
            database.GetColumn<std::string>(2),
            database.GetColumn<std::string>(3),
            database.GetColumn<std::string>(4),
            database.GetColumn<std::string>(6)});
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
            actual_primary_key_columns.emplace_back(column_info.primary_key_index, column_info.name);
        }
    }
    std::sort(actual_primary_key_columns.begin(), actual_primary_key_columns.end());

    if (actual_primary_key_columns.size() != expected_primary_key_columns.size())
    {
        throw std::runtime_error("Current schema primary key shape mismatch for table: " + table_name);
    }

    for (size_t i = 0; i < expected_primary_key_columns.size(); i++)
    {
        if (actual_primary_key_columns[i].second != expected_primary_key_columns[i])
        {
            throw std::runtime_error("Current schema primary key shape mismatch for table: " + table_name);
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
        if (foreign_key_info.from_column != from_column)
            continue;
        if (foreign_key_info.referenced_table != referenced_table)
            continue;
        if (foreign_key_info.to_column != referenced_column)
            continue;
        if (ToUpperCopy(foreign_key_info.on_delete) != expected_on_delete)
            continue;
        return;
    }
    throw std::runtime_error("Current schema is missing required foreign key on table: " + table_name);
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
                "Current schema is missing required " + object_type + " table: " + std::string(table_name));
        }
    }
}

void ValidateObjectCatalogShape(rhbm_gem::SQLiteWrapper & database)
{
    const auto columns{ QueryTableInfo(database, std::string(kCatalogTableName)) };

    const auto key_tag_iter{
        std::find_if(
            columns.begin(),
            columns.end(),
            [](const TableColumnInfo & column_info) {
                return column_info.name == "key_tag";
            }) };
    if (key_tag_iter == columns.end() || key_tag_iter->primary_key_index <= 0)
    {
        throw std::runtime_error("Current schema object_catalog must contain primary key column: key_tag");
    }

    const auto object_type_iter{
        std::find_if(
            columns.begin(),
            columns.end(),
            [](const TableColumnInfo & column_info) {
                return column_info.name == "object_type";
            }) };
    if (object_type_iter == columns.end())
    {
        throw std::runtime_error("Current schema object_catalog must contain column: object_type");
    }
    if (object_type_iter->not_null == 0)
    {
        throw std::runtime_error("Current schema object_catalog column object_type must be NOT NULL.");
    }

    const auto unexpected_object_type_count{
        QuerySingleInt(
            database,
            "SELECT COUNT(*) FROM object_catalog "
            "WHERE object_type IS NULL OR object_type NOT IN ('model', 'map');") };
    if (unexpected_object_type_count != 0)
    {
        throw std::runtime_error("Current schema object_catalog contains unsupported object_type values.");
    }
}

std::vector<std::string> QueryCatalogKeys(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & object_type)
{
    if (!HasTable(database, std::string(kCatalogTableName)))
    {
        return {};
    }
    database.Prepare("SELECT key_tag FROM object_catalog WHERE object_type = ? ORDER BY key_tag;");
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

std::vector<std::string> ListModelKeys(rhbm_gem::SQLiteWrapper & database)
{
    if (!HasTable(database, "model_object"))
    {
        return {};
    }
    return QueryKeyList(database, "SELECT key_tag FROM model_object ORDER BY key_tag;");
}

std::vector<std::string> ListMapKeys(rhbm_gem::SQLiteWrapper & database)
{
    if (!HasTable(database, std::string(kMapTableName)))
    {
        return {};
    }
    return QueryKeyList(database, "SELECT key_tag FROM map_list ORDER BY key_tag;");
}

void ValidateCatalogConsistency(
    rhbm_gem::SQLiteWrapper & database,
    std::string_view object_type,
    const std::vector<std::string> & payload_keys)
{
    const auto catalog_keys{ QueryCatalogKeys(database, std::string(object_type)) };
    const std::unordered_set<std::string> payload_key_set{ payload_keys.begin(), payload_keys.end() };
    const std::unordered_set<std::string> catalog_key_set{ catalog_keys.begin(), catalog_keys.end() };

    for (const auto & key_tag : payload_keys)
    {
        if (catalog_key_set.find(key_tag) == catalog_key_set.end())
        {
            throw std::runtime_error(
                "Current schema catalog is missing " + std::string(object_type) + " key: " + key_tag);
        }
    }

    for (const auto & key_tag : catalog_keys)
    {
        if (payload_key_set.find(key_tag) == payload_key_set.end())
        {
            throw std::runtime_error(
                "Current schema catalog contains ghost " + std::string(object_type)
                + " key without payload: " + key_tag);
        }
    }
}

void ValidateModelSchema(rhbm_gem::SQLiteWrapper & database)
{
    ValidateRequiredTables(database, kModelCanonicalTableNames, "model");

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
    for (const auto table_name : kModelCanonicalTableNames)
    {
        if (table_name == "model_object"sv)
        {
            continue;
        }
        ValidateForeignKey(database, std::string(table_name), "key_tag", "model_object", "key_tag", "CASCADE");
    }
}

void ValidateMapSchema(rhbm_gem::SQLiteWrapper & database)
{
    ValidateRequiredTables(database, kMapCanonicalTableNames, "map");
    ValidatePrimaryKeyShape(database, "map_list", {"key_tag"});
    ValidateForeignKey(database, "map_list", "key_tag", "object_catalog", "key_tag", "CASCADE");
}

void ValidateCurrentSchema(rhbm_gem::SQLiteWrapper & database)
{
    if (!HasTable(database, std::string(kCatalogTableName)))
    {
        throw std::runtime_error("Current schema is missing required table: object_catalog");
    }
    if (HasTable(database, std::string(kUnsupportedMetadataTableName)))
    {
        throw std::runtime_error("Unsupported schema/database state.");
    }

    ValidateObjectCatalogShape(database);
    ValidateModelSchema(database);
    ValidateMapSchema(database);
    ValidateCatalogConsistency(database, "model", ListModelKeys(database));
    ValidateCatalogConsistency(database, "map", ListMapKeys(database));
}

void CreateCurrentSchema(rhbm_gem::SQLiteWrapper & database)
{
    database.Execute(std::string(kCreateCatalogTableSql));
    rhbm_gem::ModelObjectStorage::CreateTables(database);
    database.Execute(std::string(kCreateMapTableSql));
    SetSchemaVersion(database);
    ValidateCurrentSchema(database);
}

void EnsureCurrentSchema(rhbm_gem::SQLiteWrapper & database)
{
    const auto raw_version{ QueryUserVersion(database) };
    if (raw_version == kCurrentSchemaVersion)
    {
        ValidateCurrentSchema(database);
        return;
    }
    if (raw_version == 0 && IsDatabaseEmpty(database))
    {
        CreateCurrentSchema(database);
        return;
    }
    throw std::runtime_error("Unsupported schema/database state.");
}

void SaveMapObject(
    rhbm_gem::SQLiteWrapper & database,
    const rhbm_gem::MapObject & map_object,
    const std::string & key_tag)
{
    database.Prepare(std::string(kInsertMapSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);

    const auto grid_size{ map_object.GetGridSize() };
    const auto grid_spacing{ map_object.GetGridSpacing() };
    const auto origin{ map_object.GetOrigin() };
    std::vector<float> values(map_object.GetMapValueArraySize());
    std::memcpy(values.data(), map_object.GetMapValueArray(), values.size() * sizeof(float));

    database.Bind<std::string>(1, key_tag);
    database.Bind<int>(2, grid_size[0]);
    database.Bind<int>(3, grid_size[1]);
    database.Bind<int>(4, grid_size[2]);
    database.Bind<double>(5, grid_spacing[0]);
    database.Bind<double>(6, grid_spacing[1]);
    database.Bind<double>(7, grid_spacing[2]);
    database.Bind<double>(8, origin[0]);
    database.Bind<double>(9, origin[1]);
    database.Bind<double>(10, origin[2]);
    database.Bind<std::vector<float>>(11, values);
    database.StepOnce();
}

std::unique_ptr<rhbm_gem::MapObject> LoadMapObject(
    rhbm_gem::SQLiteWrapper & database,
    const std::string & key_tag)
{
    database.Prepare(std::string(kSelectMapSql));
    rhbm_gem::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);

    const auto rc{ database.StepNext() };
    if (rc == rhbm_gem::SQLiteWrapper::StepDone())
    {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    }
    if (rc != rhbm_gem::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Step failed: " + database.ErrorMessage());
    }

    std::array<int, 3> grid_size{
        database.GetColumn<int>(1),
        database.GetColumn<int>(2),
        database.GetColumn<int>(3)};
    std::array<float, 3> grid_spacing{
        static_cast<float>(database.GetColumn<double>(4)),
        static_cast<float>(database.GetColumn<double>(5)),
        static_cast<float>(database.GetColumn<double>(6))};
    std::array<float, 3> origin{
        static_cast<float>(database.GetColumn<double>(7)),
        static_cast<float>(database.GetColumn<double>(8)),
        static_cast<float>(database.GetColumn<double>(9))};
    const auto values{ database.GetColumn<std::vector<float>>(10) };
    std::unique_ptr<float[]> data_array{ std::make_unique<float[]>(values.size()) };
    std::memcpy(data_array.get(), values.data(), values.size() * sizeof(float));

    auto map_object{
        std::make_unique<rhbm_gem::MapObject>(grid_size, grid_spacing, origin, std::move(data_array)) };
    map_object->SetKeyTag(database.GetColumn<std::string>(0));
    return map_object;
}

} // namespace

namespace rhbm_gem {

SQLitePersistence::SQLitePersistence(const std::filesystem::path & database_path) :
    m_database_path{ database_path }, m_database{ nullptr }, m_model_store{ nullptr }
{
    if (m_database_path.empty())
    {
        m_database_path = "database.sqlite";
    }

    const auto parent_path{ m_database_path.parent_path() };
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

    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
    m_model_store = std::make_unique<ModelObjectStorage>(m_database.get());
    EnsureCurrentSchema(*m_database);
}

SQLitePersistence::~SQLitePersistence() = default;

void SQLitePersistence::Save(const DataObjectBase & data_object, const std::string & key_tag)
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);

    const auto kind{
        io_detail::ResolveTopLevelObjectKind(data_object, "SQLitePersistence::Save()") };
    const auto type_name{ std::string(io_detail::GetCatalogTypeName(kind)) };
    m_database->Prepare(std::string(kUpsertCatalogSql));
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, type_name);
    m_database->StepOnce();

    switch (kind)
    {
    case io_detail::TopLevelObjectKind::Model:
        m_model_store->Save(static_cast<const ModelObject &>(data_object), key_tag);
        return;
    case io_detail::TopLevelObjectKind::Map:
        SaveMapObject(*m_database, static_cast<const MapObject &>(data_object), key_tag);
        return;
    }

    throw std::runtime_error("Unsupported data object type.");
}

std::unique_ptr<DataObjectBase> SQLitePersistence::Load(const std::string & key_tag)
{
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);

    m_database->Prepare(std::string(kLoadCatalogTypeSql));
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

    const auto kind{ io_detail::ParseCatalogTypeName(m_database->GetColumn<std::string>(0)) };
    switch (kind)
    {
    case io_detail::TopLevelObjectKind::Model:
        return m_model_store->Load(key_tag);
    case io_detail::TopLevelObjectKind::Map:
        return LoadMapObject(*m_database, key_tag);
    }

    throw std::runtime_error("Unsupported data object type.");
}

} // namespace rhbm_gem
