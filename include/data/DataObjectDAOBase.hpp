#pragma once

#include <memory>
#include <string>

class DataObjectBase;

class DataObjectDAOBase
{
public:
    virtual ~DataObjectDAOBase() = default;
    virtual void Save(const DataObjectBase * obj, const std::string & key_tag) = 0;
    virtual std::unique_ptr<DataObjectBase> Load(const std::string & key_tag) = 0;
};
