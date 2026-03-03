#pragma once

#include <string>

namespace rhbm_gem {

class DataObjectBase;
class PainterBase
{
public:
    virtual ~PainterBase() = default;
    virtual void SetFolder(const std::string & folder_path) = 0;
    virtual void AddDataObject(DataObjectBase * data_object) = 0;
    virtual void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) = 0;
    virtual void Painting(void) = 0;

};

} // namespace rhbm_gem
