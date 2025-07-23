#include "FilePathHelper.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <filesystem>

std::string FilePathHelper::GetExtension(std::string_view file_path)
{
    const std::filesystem::path path{ file_path };
    return (path.has_extension()) ? path.extension().string() : std::string("");
}

std::string FilePathHelper::GetDirectory(std::string_view file_path)
{
    const std::filesystem::path path{ file_path };
    auto parent{ path.parent_path() };
    if (parent.empty()) return {};
    
    std::string dir{ parent.string() };
    if (!dir.empty() && IsEndedWithSeparator(dir) == false)
    {
        dir.push_back(std::filesystem::path::preferred_separator);
    }
    return dir;
}

std::string FilePathHelper::GetFileName(std::string_view file_path)
{
    return std::filesystem::path(file_path).filename().string();
}

std::string FilePathHelper::EnsureTrailingSlash(std::string_view path)
{
    std::string result(path);
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
