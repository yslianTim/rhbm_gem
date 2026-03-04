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
};

} // namespace rhbm_gem
