#pragma once

#include <string>
#include <vector>

namespace rhbm_gem {

class SQLiteWrapper;

enum class DatabaseSchemaVersion : int
{
    LegacyV1 = 1,
    NormalizedV2 = 2
};

class DatabaseSchemaManager
{
    SQLiteWrapper * m_database;

public:
    explicit DatabaseSchemaManager(SQLiteWrapper * database);
    DatabaseSchemaVersion EnsureSchema();
    DatabaseSchemaVersion GetSchemaVersion() const;

private:
    bool IsDatabaseEmpty() const;
    bool HasTable(const std::string & table_name) const;
    bool HasLegacyModelSchema() const;
    std::vector<std::string> GetLegacyModelKeyList() const;
    std::vector<std::string> GetLegacyModelTableNameList() const;
    void EnsureNormalizedV2Schema() const;
    void MigrateLegacyV1ToNormalizedV2();
    void SetSchemaVersion(DatabaseSchemaVersion version) const;
};

} // namespace rhbm_gem
