#pragma once

#include <memory>
#include <string>

class DataObjectBase;

class FileProcessFactoryBase
{
public:
    virtual ~FileProcessFactoryBase() = default;
    virtual std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) = 0;
    virtual void OutputDataObject(const std::string & filename, DataObjectBase * data_object) = 0;
};

class ModelObjectFactory : public FileProcessFactoryBase
{
public:
    ModelObjectFactory(void) = default;
    ~ModelObjectFactory() = default;
    std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) override;
    void OutputDataObject(const std::string & filename, DataObjectBase * data_object) override;
};

class MapObjectFactory : public FileProcessFactoryBase
{
public:
    MapObjectFactory(void) = default;
    ~MapObjectFactory() = default;
    std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) override;
    void OutputDataObject(const std::string & filename, DataObjectBase * data_object) override;
};