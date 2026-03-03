#pragma once

#include <functional>
#include <memory>
#include <string>
#include <typeindex>
#include <unordered_map>

namespace rhbm_gem {

class SQLiteWrapper;
class DataObjectDAOBase;

class DataObjectDAOFactoryRegistry
{
    struct FactoryInfo
    {
        std::string name;
        std::function<std::unique_ptr<DataObjectDAOBase>(SQLiteWrapper*)> factory;
    };

    std::unordered_map<std::type_index, FactoryInfo> m_factory_map;
    std::unordered_map<std::string, std::type_index> m_name_map;

public:
    static DataObjectDAOFactoryRegistry & Instance();
    bool RegisterFactory(std::type_index type, const std::string & name,
                         std::function<std::unique_ptr<DataObjectDAOBase>(SQLiteWrapper*)> factory);
    std::unique_ptr<DataObjectDAOBase> CreateDAO(std::type_index type, SQLiteWrapper * db) const;
    std::unique_ptr<DataObjectDAOBase> CreateDAO(const std::string & name, SQLiteWrapper * db) const;
    std::string GetTypeName(std::type_index type) const;
    std::type_index GetTypeIndex(const std::string & name) const;

};

template <typename DataObjectType, typename DAOType>
class DataObjectDAORegistrar
{
public:
    explicit DataObjectDAORegistrar(const std::string & name)
    {
        DataObjectDAOFactoryRegistry::Instance().RegisterFactory(
            typeid(DataObjectType), name,
            [](SQLiteWrapper* db){ return std::make_unique<DAOType>(db); });
    }
};

} // namespace rhbm_gem
