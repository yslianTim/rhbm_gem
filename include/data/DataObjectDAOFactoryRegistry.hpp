#pragma once

#include <functional>
#include <memory>
#include <typeindex>
#include <unordered_map>

class SQLiteWrapper;
class DataObjectDAOBase;

class DataObjectDAOFactoryRegistry
{
public:
    static DataObjectDAOFactoryRegistry & Instance(void);
    bool RegisterFactory(std::type_index type,
                         std::function<std::unique_ptr<DataObjectDAOBase>(SQLiteWrapper*)> factory);
    std::unique_ptr<DataObjectDAOBase> CreateDAO(std::type_index type, SQLiteWrapper * db) const;

private:
    std::unordered_map<std::type_index,
        std::function<std::unique_ptr<DataObjectDAOBase>(SQLiteWrapper*)>> m_factory_map;
};

template <typename DataObjectType, typename DAOType>
class DataObjectDAORegistrar
{
public:
    DataObjectDAORegistrar()
    {
        DataObjectDAOFactoryRegistry::Instance().RegisterFactory(
            typeid(DataObjectType),
            [](SQLiteWrapper* db){ return std::make_unique<DAOType>(db); });
    }
};
