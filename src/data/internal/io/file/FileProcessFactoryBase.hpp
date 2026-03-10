#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class DataObjectBase;
class FileFormatRegistry;

class FileProcessFactoryBase
{
public:
    virtual ~FileProcessFactoryBase() = default;
    virtual std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) = 0;
    virtual void OutputDataObject(
        const std::string & filename,
        const DataObjectBase & data_object) = 0;
};

class ModelObjectFactory : public FileProcessFactoryBase
{
    const FileFormatRegistry & m_file_format_registry;

public:
    explicit ModelObjectFactory(const FileFormatRegistry & file_format_registry);
    ~ModelObjectFactory() = default;
    std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) override;
    void OutputDataObject(
        const std::string & filename,
        const DataObjectBase & data_object) override;
};

class MapObjectFactory : public FileProcessFactoryBase
{
    const FileFormatRegistry & m_file_format_registry;

public:
    explicit MapObjectFactory(const FileFormatRegistry & file_format_registry);
    ~MapObjectFactory() = default;
    std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) override;
    void OutputDataObject(
        const std::string & filename,
        const DataObjectBase & data_object) override;
};

} // namespace rhbm_gem
