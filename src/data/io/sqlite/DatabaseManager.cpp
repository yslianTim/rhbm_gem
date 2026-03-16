#include "internal/sqlite/DatabaseManager.hpp"
#include "internal/sqlite/DatabaseSchemaManager.hpp"
#include "internal/sqlite/MapObjectStorage.hpp"
#include "internal/sqlite/ModelObjectStorage.hpp"
#include "internal/sqlite/SQLiteWrapper.hpp"
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <system_error>

namespace rhbm_gem {

DatabaseManager::DatabaseManager(const std::filesystem::path& database_path) : m_database_path{database_path},
                                                                               m_database{nullptr},
                                                                               m_model_store{nullptr},
                                                                               m_map_store{nullptr},
                                                                               m_schema_version{DatabaseSchemaVersion::NormalizedV2} {
    if (m_database_path.empty())
        m_database_path = "database.sqlite";

    const auto parent_path{m_database_path.parent_path()};
    if (!parent_path.empty()) {
        std::error_code error_code;
        std::filesystem::create_directories(parent_path, error_code);
        if (error_code) {
            throw std::runtime_error(
                "Failed to create database parent directory '" + parent_path.string() + "': " + error_code.message());
        }
    }

    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
    m_model_store = std::make_unique<ModelObjectStorage>(m_database.get());
    m_map_store = std::make_unique<MapObjectStorage>(m_database.get());
    DatabaseSchemaManager schema_manager{m_database.get()};
    m_schema_version = schema_manager.EnsureSchema();
}

DatabaseManager::~DatabaseManager() {
}

void DatabaseManager::SaveDataObject(
    const DataObjectBase* data_object, const std::string& key_tag) {
    if (data_object == nullptr) {
        throw std::runtime_error("Null data object pointer provided.");
    }
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    auto type_name{GetCatalogTypeName(*data_object)};
    m_database->Prepare(
        "INSERT INTO object_catalog(key_tag, object_type) VALUES (?, ?) "
        "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;");
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    m_database->Bind<std::string>(2, type_name);
    m_database->StepOnce();

    if (type_name == "model") {
        m_model_store->Save(ExpectModelObject(*data_object, "DatabaseManager::SaveDataObject()"), key_tag);
        return;
    }
    if (type_name == "map") {
        m_map_store->Save(ExpectMapObject(*data_object, "DatabaseManager::SaveDataObject()"), key_tag);
        return;
    }
    throw std::runtime_error("Unsupported data object type: " + type_name);
}

std::unique_ptr<DataObjectBase> DatabaseManager::LoadDataObject(
    const std::string& key_tag) {
    std::lock_guard<std::mutex> lock(m_db_mutex);
    SQLiteWrapper::TransactionGuard transaction(*m_database);
    m_database->Prepare(
        "SELECT object_type FROM object_catalog WHERE key_tag = ? LIMIT 1;");
    SQLiteWrapper::StatementGuard guard(*m_database);
    m_database->Bind<std::string>(1, key_tag);
    auto rc{m_database->StepNext()};
    if (rc == SQLiteWrapper::StepDone()) {
        throw std::runtime_error("Cannot find the row with key_tag = " + key_tag);
    } else if (rc != SQLiteWrapper::StepRow()) {
        throw std::runtime_error("Step failed: " + m_database->ErrorMessage());
    }

    auto type_name{m_database->GetColumn<std::string>(0)};
    if (type_name == "model") {
        return m_model_store->Load(key_tag);
    }
    if (type_name == "map") {
        return m_map_store->Load(key_tag);
    }
    throw std::runtime_error("Unsupported data object type: " + type_name);
}

} // namespace rhbm_gem
