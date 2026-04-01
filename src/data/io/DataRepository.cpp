#include <rhbm_gem/data/io/DataRepository.hpp>

#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include "sqlite/SQLitePersistence.hpp"

#include <stdexcept>

namespace rhbm_gem {

DataRepository::DataRepository(const std::filesystem::path & database_path) :
    m_database{ std::make_unique<SQLitePersistence>(database_path) }
{
}

DataRepository::~DataRepository() = default;

std::unique_ptr<ModelObject> DataRepository::LoadModel(const std::string & key_tag) const
{
    return m_database->LoadModel(key_tag);
}

std::unique_ptr<MapObject> DataRepository::LoadMap(const std::string & key_tag) const
{
    return m_database->LoadMap(key_tag);
}

void DataRepository::SaveModel(const ModelObject & model_object, const std::string & key_tag) const
{
    m_database->SaveModel(model_object, key_tag);
}

void DataRepository::SaveMap(const MapObject & map_object, const std::string & key_tag) const
{
    m_database->SaveMap(map_object, key_tag);
}

} // namespace rhbm_gem
