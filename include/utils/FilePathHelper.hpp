#pragma once

#include <string>

class FilePathHelper
{

public:
    static std::string GetExtension(const std::string & file_path);
    static std::string GetDirectory(const std::string & file_path);
    static std::string GetFileName(const std::string & file_path);
    static std::string EnsureTrailingSlash(const std::string & path);
};
