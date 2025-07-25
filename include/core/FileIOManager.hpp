#pragma once

#include <filesystem>
#include <memory>
#include <string>

class DataObjectBase;

class FileIOManager
{
public:
    FileIOManager(void) = delete;
    ~FileIOManager() = delete;

    static std::shared_ptr<DataObjectBase> LoadDataObject(const std::filesystem::path & filename,
                                                          const std::string & key_tag);
    static void WriteDataObject(const std::filesystem::path & filename,
                                DataObjectBase * data_object);
};
