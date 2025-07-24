#pragma once

#include <memory>
#include <string>

class DataObjectVisitorBase;

class DataObjectBase
{
public:
    virtual ~DataObjectBase() = default;
    virtual std::unique_ptr<DataObjectBase> Clone(void) const = 0;
    virtual void Display(void) const = 0;
    virtual void Update(void) = 0;
    virtual void Accept(class DataObjectVisitorBase * visitor) = 0;
    virtual void SetKeyTag(const std::string & label) = 0;
    virtual std::string GetKeyTag(void) const = 0;
};
