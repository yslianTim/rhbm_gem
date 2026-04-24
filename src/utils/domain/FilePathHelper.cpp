#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <cctype>
#include <string_view>

namespace rhbm_gem {
namespace {

constexpr bool IsEndedWithSeparator(std::string_view file_path)
{
    if (file_path.empty()) return false;
    const char c{ file_path.back() };
    return c == '/' || c == '\\' || c == std::filesystem::path::preferred_separator;
}

} // namespace

namespace path_helper {

std::string GetExtension(const std::filesystem::path & path)
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

    string_helper::ToLowerCase(extension);
    return extension;
}

std::string GetDirectory(const std::filesystem::path & path)
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

std::string GetFileName(
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

std::string EnsureTrailingSlash(const std::filesystem::path & path)
{
    std::string result(path.string());
    if (!result.empty() && IsEndedWithSeparator(result) == false)
    {
        result.push_back(std::filesystem::path::preferred_separator);
    }
    return result;
}

bool EnsureFileExists(const std::filesystem::path & path,
                      const std::string & log_prefix)
{
    if (std::filesystem::exists(path))
    {
        return true;
    }
    Logger::Log(LogLevel::Error, log_prefix + " does not exist: " + path.string());
    return false;
}

} // namespace path_helper
} // namespace rhbm_gem
