#pragma once

#include <string>
#include <string_view>
#include <filesystem>

class FilePathHelper
{

public:
    static std::string GetExtension(const std::filesystem::path & path);
    static std::string GetDirectory(const std::filesystem::path & path);
    static std::string GetFileName(const std::filesystem::path & path, bool include_extension = true);
    static std::string EnsureTrailingSlash(const std::filesystem::path & path);
    static bool EnsureFileExists(const std::filesystem::path & path, const std::string & log_prefix);

private:
    static constexpr bool IsEndedWithSeparator(std::string_view path);
};
