#include "DataObjectManager.hpp"
#include "FileProcessFactoryBase.hpp"
#include "FilePathHelper.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectVisitorBase.hpp"
#include "DatabaseManager.hpp"
#include "ModelObject.hpp"
#include "ScopeTimer.hpp"

DataObjectManager::DataObjectManager(void) :
    m_db_manager{ nullptr }
{

}

DataObjectManager::DataObjectManager(const std::string & dbname) :
    m_db_manager{ std::make_unique<DatabaseManager>(dbname) }
{

}

DataObjectManager::~DataObjectManager()
{

}

std::unique_ptr<FileProcessFactoryBase> DataObjectManager::CreateFactory(const std::string & file_extension)
{
    if (file_extension == ".pdb" || file_extension == ".cif")
    {
        return std::make_unique<ModelObjectFactory>();
    }
    else if (file_extension == ".mrc" || file_extension == ".map")
    {
        return std::make_unique<MapObjectFactory>();
    }
    else
    {
        throw std::runtime_error("Unsupported file format");
    }
}

void DataObjectManager::ProcessFile(const std::string & filename, const std::string & key_tag)
{
    ScopeTimer timer("DataObjectManager::ProcessFile");
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ CreateFactory(file_extension) };
    auto data_object{ factory->CreateDataObject(filename) };
    if (data_object == nullptr)
    {
        throw std::runtime_error("Failed to create data object");
    }
    data_object->SetKeyTag(key_tag);
    data_object->Display();

    if (m_data_object_map.find(key_tag) != m_data_object_map.end())
    {
        std::cout <<"The key tag: ["<< key_tag <<"] already presented"
                  <<", this data object will be replaced."<< std::endl;
    }
    m_data_object_map.insert_or_assign(key_tag, std::move(data_object));
}

void DataObjectManager::ProduceFile(const std::string & filename, const std::string & key_tag)
{
    ScopeTimer timer("DataObjectManager::ProduceFile");
    if (m_data_object_map.find(key_tag) == m_data_object_map.end())
    {
        std::cout <<"The data object with key tag: ["<< key_tag <<"] is not exists"
                  <<", no file will be produced."<< std::endl;
        return;
    }
    auto data_object{ m_data_object_map.at(key_tag).get() };
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ CreateFactory(file_extension) };
    factory->OutputDataObject(filename, data_object);
};

void DataObjectManager::AddDataObject(
    const std::string & key_tag, std::unique_ptr<DataObjectBase> data_object)
{
    if (m_data_object_map.find(key_tag) != m_data_object_map.end())
    {
        std::cout <<"[Warning] The key tag: ["<< key_tag <<"] already presented in the data object map"
                  <<", this data object will be replaced."<< std::endl;
    }
    m_data_object_map.insert_or_assign(key_tag, std::move(data_object));
}

void DataObjectManager::LoadDataObject(const std::string & key_tag)
{
    ScopeTimer timer("DataObjectManager::LoadDataObject");
    if (m_db_manager == nullptr)
    {
        throw std::runtime_error("Database manager is not initialized.");
    }
    if (m_data_object_map.find(key_tag) != m_data_object_map.end())
    {
        std::cout <<"The key tag: ["<< key_tag <<"] already presented"
                  <<", this data object will be replaced."<< std::endl;
    }
    auto data_object{ m_db_manager->LoadDataObject(key_tag) };
    m_data_object_map.insert_or_assign(key_tag, std::move(data_object));
}

void DataObjectManager::SaveDataObject(
    const std::string & key_tag, const std::string & renamed_key_tag) const
{
    ScopeTimer timer("DataObjectManager::SaveDataObject");
    if (m_db_manager == nullptr)
    {
        throw std::runtime_error("Database manager is not initialized.");
    }
    if (m_data_object_map.find(key_tag) == m_data_object_map.end())
    {
        std::cout <<"The key tag: ["<< key_tag <<"] isn't presented"
                  <<", skip saving data object."<< std::endl;
        return;
    }

    if (renamed_key_tag != "")
    {
        std::cout <<"The key tag of data object will be renamed as: "
                  <<"["<< renamed_key_tag <<"]"
                  <<" and saved into database: "<< m_db_manager->GetDatabasePath() << std::endl;
        m_data_object_map.at(key_tag)->SetKeyTag(renamed_key_tag);
    }
    else
    {
        std::cout <<"The data object ["<< key_tag <<"] will be saved into database: "
                  << m_db_manager->GetDatabasePath() << std::endl;
    }

    m_db_manager->SaveDataObject(m_data_object_map.at(key_tag).get());
}

void DataObjectManager::Accept(DataObjectVisitorBase * visitor)
{
    ScopeTimer timer("DataObjectManager::Accept");
    visitor->Analysis(this);
}

void DataObjectManager::PrintDataObjectInfo(const std::string & key_tag) const
{
    try
    {
        m_data_object_map.at(key_tag)->Display();
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

const std::unique_ptr<DataObjectBase> & DataObjectManager::GetDataObjectRef(const std::string & key_tag)
{
    if (m_data_object_map.find(key_tag) == m_data_object_map.end())
    {
        throw std::runtime_error("Cannot find the data object with key tag: " + key_tag);
    }
    return m_data_object_map.at(key_tag);
}
