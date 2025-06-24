#pragma once

#include <string>
#include <string_view>

class FilePathHelper
{

public:
    static std::string GetExtension(std::string_view file_path);
    static std::string GetDirectory(std::string_view file_path);
    static std::string GetFileName(std::string_view file_path);
    static std::string EnsureTrailingSlash(std::string_view path);

private:
    static constexpr bool IsEndedWithSeparator(std::string_view file_path);
};
