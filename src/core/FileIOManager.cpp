#include "FileIOManager.hpp"
#include "FileProcessFactoryRegistry.hpp"
#include "FileProcessFactoryBase.hpp"
#include "FilePathHelper.hpp"
#include "DataObjectBase.hpp"
#include "Logger.hpp"

std::shared_ptr<DataObjectBase> FileIOManager::LoadDataObject(
    const std::filesystem::path & filename, const std::string & key_tag)
{
    Logger::Log(LogLevel::Debug, "FileIOManager::LoadDataObject() called");
    auto extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(extension) };
    auto data_object{ factory->CreateDataObject(filename) };
    if (!data_object)
    {
        throw std::runtime_error("Failed to create data object");
    }
    data_object->SetKeyTag(key_tag);
    data_object->Display();
    std::shared_ptr<DataObjectBase> shared_object{ std::move(data_object) };
    return shared_object;
}

void FileIOManager::WriteDataObject(const std::filesystem::path & filename,
                                    DataObjectBase * data_object)
{
    Logger::Log(LogLevel::Debug, "FileIOManager::WriteDataObject() called");
    auto extension{ FilePathHelper::GetExtension(filename) };
    auto factory{ FileProcessFactoryRegistry::Instance().CreateFactory(extension) };
    factory->OutputDataObject(filename, data_object);
}
