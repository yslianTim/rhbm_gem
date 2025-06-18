#include "FilePathHelper.hpp"

#include <algorithm>

std::string FilePathHelper::GetExtension(const std::string & file_path)
{
    std::size_t pos{ file_path.find_last_of('.') };
    if (pos != std::string::npos && pos + 1 < file_path.length())
    {
        return file_path.substr(pos);
    }
    return std::string("");
}

std::string FilePathHelper::GetDirectory(const std::string & file_path)
{
    std::size_t pos1{ file_path.find_last_of('/') };
    std::size_t pos2{ file_path.find_last_of('\\') };
    std::size_t pos{ std::max(pos1, pos2) };
    if (pos != std::string::npos)
    {
        return file_path.substr(0, pos + 1);
    }
    return std::string("");
}

std::string FilePathHelper::GetFileName(const std::string & file_path)
{
    std::size_t pos1{ file_path.find_last_of('/') };
    std::size_t pos2{ file_path.find_last_of('\\') };
    std::size_t pos{ std::max(pos1, pos2) };
    if (pos != std::string::npos && pos + 1 < file_path.length())
    {
        return file_path.substr(pos + 1);
    }
    return file_path;
}

std::string FilePathHelper::EnsureTrailingSlash(const std::string & path)
{
    if (!path.empty() && path.back() != '/' && path.back() != '\\')
    {
        if (path.find('\\') != std::string::npos)
        {
            return path + "\\";
        }
        else
        {
            return path + "/";
        }
    }
    return path;
}
