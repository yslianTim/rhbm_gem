#include "DataObjectDAOFactoryRegistry.hpp"
#include "DataObjectDAOBase.hpp"
#include "Logger.hpp"
#include <stdexcept>

DataObjectDAOFactoryRegistry & DataObjectDAOFactoryRegistry::Instance(void)
{
    static DataObjectDAOFactoryRegistry instance;
    return instance;
}

bool DataObjectDAOFactoryRegistry::RegisterFactory(
    std::type_index type,
    std::function<std::unique_ptr<DataObjectDAOBase>(SQLiteWrapper*)> factory)
{
    Logger::Log(LogLevel::Debug, "DataObjectDAOFactoryRegistry::RegisterFactory() called");
    m_factory_map[type] = std::move(factory);
    return true;
}

std::unique_ptr<DataObjectDAOBase> DataObjectDAOFactoryRegistry::CreateDAO(
    std::type_index type, SQLiteWrapper * db) const
{
    Logger::Log(LogLevel::Debug, "DataObjectDAOFactoryRegistry::CreateDAO() called");
    auto iter{ m_factory_map.find(type) };
    if (iter == m_factory_map.end())
    {
        throw std::runtime_error("Unsupported data object type.");
    }
    return (iter->second)(db);
}
