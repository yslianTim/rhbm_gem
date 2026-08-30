#pragma once

#include <filesystem>
#include <memory>
#include <mutex>
#include <string>

namespace rhbm_gem {

class ModelObject;
class SQLiteWrapper;

class DataRepository
{
    std::unique_ptr<SQLiteWrapper> m_database;
    mutable std::mutex m_db_mutex;

public:
    explicit DataRepository(const std::filesystem::path & database_path);
    ~DataRepository();
    DataRepository(const DataRepository &) = delete;
    DataRepository & operator=(const DataRepository &) = delete;

    std::unique_ptr<ModelObject> LoadModel(const std::string & key_tag) const;
    void SaveModel(const ModelObject & model_object, const std::string & key_tag) const;
};

} // namespace rhbm_gem
