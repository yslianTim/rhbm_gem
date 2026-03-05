#pragma once

#include <memory>
#include <string>

#include "ModelVisitMode.hpp"

namespace rhbm_gem {

class DataObjectVisitor;
class ConstDataObjectVisitor;

class DataObjectBase
{
public:
    virtual ~DataObjectBase() = default;
    virtual std::unique_ptr<DataObjectBase> Clone() const = 0;
    virtual void Display() const = 0;
    virtual void Update() = 0;
    virtual void Accept(DataObjectVisitor & visitor) = 0;
    virtual void Accept(ConstDataObjectVisitor & visitor) const = 0;
    virtual void Accept(DataObjectVisitor & visitor, ModelVisitMode model_mode) = 0;
    virtual void Accept(ConstDataObjectVisitor & visitor, ModelVisitMode model_mode) const = 0;
    virtual void SetKeyTag(const std::string & label) = 0;
    virtual std::string GetKeyTag() const = 0;
};

} // namespace rhbm_gem
