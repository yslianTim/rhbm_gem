#include "DataObjectDAOFactoryRegistry.hpp"
#include "DataObjectDAOBase.hpp"
#include "Logger.hpp"
#include <stdexcept>

namespace rhbm_gem {

DataObjectDAOFactoryRegistry & DataObjectDAOFactoryRegistry::Instance()
{
    static DataObjectDAOFactoryRegistry instance;
    return instance;
}

bool DataObjectDAOFactoryRegistry::RegisterFactory(
    std::type_index type, const std::string & name,
    std::function<std::unique_ptr<DataObjectDAOBase>(SQLiteWrapper*)> factory)
{
    FactoryInfo info{ name, std::move(factory) };
    m_factory_map[type] = std::move(info);
    m_name_map.emplace(name, type);
    return true;
}

std::unique_ptr<DataObjectDAOBase> DataObjectDAOFactoryRegistry::CreateDAO(
    std::type_index type, SQLiteWrapper * db) const
{
    auto iter{ m_factory_map.find(type) };
    if (iter == m_factory_map.end())
    {
        throw std::runtime_error("Unsupported data object type.");
    }
    return (iter->second.factory)(db);
}

std::unique_ptr<DataObjectDAOBase> DataObjectDAOFactoryRegistry::CreateDAO(
    const std::string & name, SQLiteWrapper * db) const
{
    auto iter{ m_name_map.find(name) };
    if (iter == m_name_map.end())
    {
        throw std::runtime_error("Unsupported data object type: " + name);
    }
    return CreateDAO(iter->second, db);
}

std::string DataObjectDAOFactoryRegistry::GetTypeName(std::type_index type) const
{
    auto iter{ m_factory_map.find(type) };
    if (iter == m_factory_map.end())
    {
        throw std::runtime_error("Unsupported data object type.");
    }
    return iter->second.name;
}

std::type_index DataObjectDAOFactoryRegistry::GetTypeIndex(const std::string & name) const
{
    auto iter{ m_name_map.find(name) };
    if (iter == m_name_map.end())
    {
        throw std::runtime_error("Unsupported data object type: " + name);
    }
    return iter->second;
}

} // namespace rhbm_gem
