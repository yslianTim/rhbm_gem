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
    auto base_object{ m_database->Load(key_tag) };
    auto typed_object{ dynamic_cast<ModelObject *>(base_object.get()) };
    if (!typed_object)
    {
        throw std::runtime_error("Invalid data type for " + key_tag);
    }
    (void)base_object.release();
    return std::unique_ptr<ModelObject>{ typed_object };
}

std::unique_ptr<MapObject> DataRepository::LoadMap(const std::string & key_tag) const
{
    auto base_object{ m_database->Load(key_tag) };
    auto typed_object{ dynamic_cast<MapObject *>(base_object.get()) };
    if (!typed_object)
    {
        throw std::runtime_error("Invalid data type for " + key_tag);
    }
    (void)base_object.release();
    return std::unique_ptr<MapObject>{ typed_object };
}

void DataRepository::SaveModel(const ModelObject & model_object, const std::string & key_tag) const
{
    m_database->Save(model_object, key_tag);
}

void DataRepository::SaveMap(const MapObject & map_object, const std::string & key_tag) const
{
    m_database->Save(map_object, key_tag);
}

} // namespace rhbm_gem
