#pragma once

#include <string>

#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

namespace rhbm_gem {

class PainterBase
{
public:
    virtual ~PainterBase() = default;
    virtual void SetFolder(const std::string & folder_path)
    {
        m_folder_path = path_helper::EnsureTrailingSlash(folder_path);
    }
    virtual void Painting() = 0;

protected:
    std::string m_folder_path{ "./" };
};

} // namespace rhbm_gem
