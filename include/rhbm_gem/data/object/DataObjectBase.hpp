#pragma once

#include <memory>
#include <string>

namespace rhbm_gem {

class DataObjectBase
{
public:
    virtual ~DataObjectBase() = default;
    virtual std::unique_ptr<DataObjectBase> Clone() const = 0;
    virtual void Display() const = 0;
    virtual void Update() = 0;
    virtual void SetKeyTag(const std::string & label) = 0;
    virtual std::string GetKeyTag() const = 0;
};

} // namespace rhbm_gem
