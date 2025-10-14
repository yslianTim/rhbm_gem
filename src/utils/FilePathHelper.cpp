#include "FilePathHelper.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <cctype>

std::string FilePathHelper::GetExtension(const std::filesystem::path & path)
{
    const std::string filename{ path.filename().string() };
    std::string extension;
    if (!filename.empty() && filename.front() == '.' && filename.find('.', 1) == std::string::npos)
    {
        extension = filename;
    }
    else if (path.has_extension())
    {
        extension = path.extension().string();
    }

    std::transform(
        extension.begin(), extension.end(), extension.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); }
    );
    return extension;
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
    std::filesystem::path normalized{ path };
    std::string native{ normalized.string() };
    if (native.find('\\') != std::string::npos)
    {
        std::replace(native.begin(), native.end(), '\\', '/');
        normalized = std::filesystem::path{ native };
    }

    const std::filesystem::path filename{ normalized.filename() };
    if (include_extension)
    {
        return filename.string();
    }

    if (filename.empty())
    {
        return filename.string();
    }

    const std::filesystem::path stem{ filename.stem() };
    if (stem.empty())
    {
        return filename.string();
    }

    return stem.string();
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
