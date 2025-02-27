#include "DatabaseManager.hpp"
#include "SQLiteWrapper.hpp"
#include "DataObjectDAOBase.hpp"
#include "ModelObjectDAO.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"

DatabaseManager::DatabaseManager(const std::string & database_path) :
    m_database{ std::make_unique<SQLiteWrapper>(database_path) }
{

}

DatabaseManager::~DatabaseManager()
{

}

void DatabaseManager::SaveDataObject(const DataObjectBase * obj)
{
    auto dao{ CreateDataObjectDAO(obj) };
    dao->Save(obj);
}

std::unique_ptr<DataObjectBase> DatabaseManager::LoadDataObject(const std::string & key_tag)
{
    auto obj{ std::make_unique<ModelObject>() };
    auto dao{ CreateDataObjectDAO(obj.get()) };
    return dao->Load(key_tag);
}

std::unique_ptr<DataObjectDAOBase> DatabaseManager::CreateDataObjectDAO(const DataObjectBase * obj)
{
    if (dynamic_cast<const ModelObject *>(obj))
    {
        return std::make_unique<ModelObjectDAO>(*this);
    }
    else
    {
        throw std::runtime_error("Unsupported data object type.");
    }
}