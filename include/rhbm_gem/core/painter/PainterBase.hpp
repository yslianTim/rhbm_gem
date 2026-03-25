#pragma once

#include <string>

#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

namespace rhbm_gem {

class DataObjectBase;
class PainterBase
{
public:
    virtual ~PainterBase() = default;
    virtual void SetFolder(const std::string & folder_path)
    {
        m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
    }
    virtual void AddDataObject(DataObjectBase * data_object) = 0;
    virtual void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
    {
        (void)data_object;
        (void)label;
    }
    virtual void Painting() = 0;

protected:
    std::string m_folder_path{ "./" };
};

} // namespace rhbm_gem
