#include "DatabaseManager.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectDAOBase.hpp"
#include "ModelObjectDAO.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"

DatabaseManager::DatabaseManager(const std::filesystem::path & database_path) :
    m_database_path{ database_path },
    m_database{ nullptr }
{
    if (m_database_path.empty()) m_database_path = "database.sqlite";
    m_database = std::make_unique<SQLiteWrapper>(m_database_path);
}

DatabaseManager::~DatabaseManager()
{

}

void DatabaseManager::SaveDataObject(const DataObjectBase * data_object, const std::string & key_tag)
{
    auto dao{ CreateDataObjectDAO(data_object) };
    dao->Save(data_object, key_tag);
}

std::unique_ptr<DataObjectBase> DatabaseManager::LoadDataObject(const std::string & key_tag)
{
    auto data_object{ std::make_unique<ModelObject>() };
    auto dao{ CreateDataObjectDAO(data_object.get()) };
    return dao->Load(key_tag);
}

std::unique_ptr<DataObjectDAOBase> DatabaseManager::CreateDataObjectDAO(const DataObjectBase * data_object)
{
    if (dynamic_cast<const ModelObject *>(data_object))
    {
        return std::make_unique<ModelObjectDAO>(m_database.get());
    }
    else
    {
        throw std::runtime_error("Unsupported data object type.");
    }
}
