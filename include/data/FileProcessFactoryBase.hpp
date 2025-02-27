#pragma once

#include <iostream>
#include <memory>
#include <string>

class DataObjectBase;

class FileProcessFactoryBase
{
public:
    virtual ~FileProcessFactoryBase() = default;
    virtual std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) = 0;
};

class ModelObjectFactory : public FileProcessFactoryBase
{
public:
    ModelObjectFactory(void) = default;
    ~ModelObjectFactory() = default;
    std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) override;
};

class MapObjectFactory : public FileProcessFactoryBase
{
public:
    MapObjectFactory(void) = default;
    ~MapObjectFactory() = default;
    std::unique_ptr<DataObjectBase> CreateDataObject(const std::string & filename) override;
};