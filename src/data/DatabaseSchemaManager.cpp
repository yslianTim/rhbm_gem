#include "DatabaseSchemaManager.hpp"

#include "DataObjectBase.hpp"
#include "LegacyModelObjectDAO.hpp"
#include "MapObjectDAO.hpp"
#include "ModelObjectDAOv2.hpp"
#include "SQLiteWrapper.hpp"
#include "Logger.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

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
    constexpr std::string_view MAP_KEY_LIST_SQL =
        "SELECT key_tag FROM map_list;";
    constexpr std::string_view MODEL_V2_KEY_LIST_SQL =
        "SELECT key_tag FROM model_object;";

    bool IsLegacyModelTableName(const std::string & table_name)
    {
        if (table_name == "model_list") return true;
        constexpr std::string_view patterns[]{
            "chemical_component_entry_in_",
            "component_atom_entry_in_",
            "component_bond_entry_in_",
            "atom_list_in_",
            "bond_list_in_",
            "atom_local_potential_entry_in_",
            "bond_local_potential_entry_in_",
            "_atom_local_potential_entry_in_",
            "_bond_local_potential_entry_in_",
            "_atom_group_potential_entry_in_",
            "_bond_group_potential_entry_in_"
        };
        for (const auto pattern : patterns)
        {
            if (table_name.find(pattern) != std::string::npos)
            {
                return true;
            }
        }
        return false;
    }

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
    const auto raw_version{ [&]() -> int
    {
        m_database->Prepare(PRAGMA_USER_VERSION_SQL.data());
        SQLiteWrapper::StatementGuard guard(*m_database);
        const auto rc{ m_database->StepNext() };
        if (rc == SQLiteWrapper::StepDone())
        {
            return 0;
        }
        if (rc != SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Failed to query schema version: " + m_database->ErrorMessage());
        }
        return m_database->GetColumn<int>(0);
    }() };

    if (raw_version == 0)
    {
        if (IsDatabaseEmpty())
        {
            EnsureNormalizedV2Schema();
            SetSchemaVersion(DatabaseSchemaVersion::NormalizedV2);
            return DatabaseSchemaVersion::NormalizedV2;
        }
        if (HasLegacyModelSchema())
        {
            MigrateLegacyV1ToNormalizedV2();
            return DatabaseSchemaVersion::NormalizedV2;
        }
        EnsureNormalizedV2Schema();
        SetSchemaVersion(DatabaseSchemaVersion::NormalizedV2);
        return DatabaseSchemaVersion::NormalizedV2;
    }

    if (raw_version == static_cast<int>(DatabaseSchemaVersion::LegacyV1))
    {
        MigrateLegacyV1ToNormalizedV2();
        return DatabaseSchemaVersion::NormalizedV2;
    }

    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        EnsureNormalizedV2Schema();
        return DatabaseSchemaVersion::NormalizedV2;
    }

    throw std::runtime_error("Unsupported database schema version: " + std::to_string(raw_version));
}

DatabaseSchemaVersion DatabaseSchemaManager::GetSchemaVersion() const
{
    m_database->Prepare(PRAGMA_USER_VERSION_SQL.data());
    SQLiteWrapper::StatementGuard guard(*m_database);
    const auto rc{ m_database->StepNext() };
    if (rc == SQLiteWrapper::StepDone())
    {
        return DatabaseSchemaVersion::LegacyV1;
    }
    if (rc != SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Failed to query schema version: " + m_database->ErrorMessage());
    }

    const auto raw_version{ m_database->GetColumn<int>(0) };
    if (raw_version == static_cast<int>(DatabaseSchemaVersion::NormalizedV2))
    {
        return DatabaseSchemaVersion::NormalizedV2;
    }
    return DatabaseSchemaVersion::LegacyV1;
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

std::vector<std::string> DatabaseSchemaManager::GetLegacyModelKeyList() const
{
    if (HasTable("object_metadata"))
    {
        auto values{ QuerySingleStringColumn(*m_database, std::string(MODEL_METADATA_KEY_LIST_SQL)) };
        if (!values.empty())
        {
            return values;
        }
    }
    if (HasTable("model_list"))
    {
        return QuerySingleStringColumn(*m_database, std::string(LEGACY_MODEL_KEY_LIST_SQL));
    }
    return {};
}

std::vector<std::string> DatabaseSchemaManager::GetLegacyModelTableNameList() const
{
    std::vector<std::string> table_names;
    for (const auto & table_name : QuerySingleStringColumn(*m_database, std::string(LIST_TABLES_SQL)))
    {
        if (IsLegacyModelTableName(table_name))
        {
            table_names.push_back(table_name);
        }
    }
    return table_names;
}

void DatabaseSchemaManager::EnsureNormalizedV2Schema() const
{
    m_database->Execute(std::string(CREATE_OBJECT_METADATA_SQL));
    ModelObjectDAOv2::EnsureSchema(*m_database);
    MapObjectDAO::EnsureSchema(*m_database);

    if (HasTable("map_list"))
    {
        const auto key_list{ QuerySingleStringColumn(*m_database, std::string(MAP_KEY_LIST_SQL)) };
        m_database->Prepare(std::string(UPSERT_OBJECT_METADATA_SQL));
        SQLiteWrapper::StatementGuard guard(*m_database);
        for (const auto & key_tag : key_list)
        {
            m_database->Bind<std::string>(1, key_tag);
            m_database->Bind<std::string>(2, "map");
            m_database->StepOnce();
            m_database->Reset();
        }
    }

    if (HasTable("model_object"))
    {
        const auto key_list{ QuerySingleStringColumn(*m_database, std::string(MODEL_V2_KEY_LIST_SQL)) };
        m_database->Prepare(std::string(UPSERT_OBJECT_METADATA_SQL));
        SQLiteWrapper::StatementGuard guard(*m_database);
        for (const auto & key_tag : key_list)
        {
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
    EnsureNormalizedV2Schema();

    LegacyModelObjectDAO legacy_dao{ m_database };
    ModelObjectDAOv2 v2_dao{ m_database };

    for (const auto & key_tag : GetLegacyModelKeyList())
    {
        auto model_object{ legacy_dao.Load(key_tag) };
        v2_dao.Save(model_object.get(), key_tag);

        m_database->Prepare(std::string(UPSERT_OBJECT_METADATA_SQL));
        SQLiteWrapper::StatementGuard guard(*m_database);
        m_database->Bind<std::string>(1, key_tag);
        m_database->Bind<std::string>(2, "model");
        m_database->StepOnce();
    }

    for (const auto & table_name : GetLegacyModelTableNameList())
    {
        m_database->Execute("DROP TABLE IF EXISTS " + table_name + ";");
    }

    EnsureNormalizedV2Schema();
    SetSchemaVersion(DatabaseSchemaVersion::NormalizedV2);
}

void DatabaseSchemaManager::SetSchemaVersion(DatabaseSchemaVersion version) const
{
    m_database->Execute(
        "PRAGMA user_version = " + std::to_string(static_cast<int>(version)) + ";");
}

} // namespace rhbm_gem
