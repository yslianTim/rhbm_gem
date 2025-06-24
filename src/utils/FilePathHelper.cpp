#include "FilePathHelper.hpp"

#include <algorithm>
#include <filesystem>

std::string FilePathHelper::GetExtension(const std::string & file_path)
{
    std::filesystem::path path{ file_path };
    return (path.has_extension()) ? path.extension().string() : std::string("");
}

std::string FilePathHelper::GetDirectory(const std::string & file_path)
{
    std::filesystem::path path{ file_path };
    auto parent{ path.parent_path() };
    if (!parent.empty())
    {
        std::string dir{ parent.string() };
        if (!dir.empty() && IsEndedWithSeparator(path) == false)
        {
            dir.push_back(std::filesystem::path::preferred_separator);
        }
        return dir;
    }
    return std::string("");
}

std::string FilePathHelper::GetFileName(const std::string & file_path)
{
    std::filesystem::path path{ file_path };
    return path.filename().string();
}

std::string FilePathHelper::EnsureTrailingSlash(const std::string & path)
{
    if (!path.empty() && IsEndedWithSeparator(path) == false)
    {
        return path + std::filesystem::path::preferred_separator;
    }
    return path;
}

bool FilePathHelper::IsEndedWithSeparator(const std::string & file_path)
{
    if (file_path.back() == '/' || file_path.back() == '\\') return true;
    if (file_path.back() == std::filesystem::path::preferred_separator) return true;
    return false;
}
