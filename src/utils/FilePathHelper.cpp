#include "FilePathHelper.hpp"
#include "Logger.hpp"

#include <algorithm>

std::string FilePathHelper::GetExtension(const std::filesystem::path & path)
{
    const std::string filename{ path.filename().string() };
    if (!filename.empty() && filename.front() == '.' && filename.find('.', 1) == std::string::npos)
    {
        return filename;
    }
    return (path.has_extension()) ? path.extension().string() : std::string("");
}

std::string FilePathHelper::GetDirectory(const std::filesystem::path & path)
{
    auto parent{ path.parent_path() };
    if (parent.empty()) return {};
    
    std::string dir{ parent.string() };
    if (!dir.empty() && IsEndedWithSeparator(dir) == false)
    {
        dir.push_back(std::filesystem::path::preferred_separator);
    }
    return dir;
}

std::string FilePathHelper::GetFileName(
    const std::filesystem::path & path, bool include_extension)
{
    std::string normalized{ path.string() };
    std::replace(normalized.begin(), normalized.end(), '\\', '/');
    if (include_extension == false)
    {
        auto pos{ normalized.find_last_of('.') };
        if (pos != std::string::npos)
        {
            normalized.erase(pos);
        }
    }
    return std::filesystem::path{ normalized }.filename().string();
}

std::string FilePathHelper::EnsureTrailingSlash(const std::filesystem::path & path)
{
    std::string result(path.string());
    if (!result.empty() && IsEndedWithSeparator(result) == false)
    {
        result.push_back(std::filesystem::path::preferred_separator);
    }
    return result;
}

bool FilePathHelper::EnsureFileExists(const std::filesystem::path & path,
                                      const std::string & log_prefix)
{
    if (std::filesystem::exists(path))
    {
        return true;
    }
    Logger::Log(LogLevel::Error, log_prefix + " does not exist: " + path.string());
    return false;
}

constexpr bool FilePathHelper::IsEndedWithSeparator(std::string_view file_path)
{
    if (file_path.empty()) return false;
    const char c{ file_path.back() };
    return c == '/' || c == '\\' || c == std::filesystem::path::preferred_separator;
}
