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

enum class DatabaseSchemaState
{
    Empty,
    LegacyV1,
    NormalizedV2,
    ManagedButUnversioned,
    MixedUnknown
};

class DatabaseSchemaManager
{
    SQLiteWrapper * m_database;

public:
    explicit DatabaseSchemaManager(SQLiteWrapper * database);
    DatabaseSchemaVersion EnsureSchema();
    DatabaseSchemaVersion GetSchemaVersion() const;

private:
    DatabaseSchemaState InspectSchemaState() const;
    bool IsDatabaseEmpty() const;
    bool HasTable(const std::string & table_name) const;
    bool HasLegacyModelSchema() const;
    void ValidateNormalizedV2Schema() const;
    std::vector<std::string> GetLegacyModelKeyList() const;
    std::vector<std::string> BuildLegacyOwnedTableNames(
        const std::vector<std::string> & model_key_list) const;
    void EnsureNormalizedV2Tables() const;
    void RepairManagedMetadata() const;
    void MigrateLegacyV1ToNormalizedV2();
    void SetSchemaVersion(DatabaseSchemaVersion version) const;
};

} // namespace rhbm_gem
