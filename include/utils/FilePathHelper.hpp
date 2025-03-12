#pragma once

#include <string>
#include <algorithm>

class FilePathHelper
{
public:
    
    static std::string GetExtension(const std::string & file_path)
    {
        std::size_t pos{ file_path.find_last_of('.') };
        if (pos != std::string::npos && pos + 1 < file_path.length())
        {
            return file_path.substr(pos);
        }
        return "";
    }
    
    static std::string GetDirectory(const std::string & file_path)
    {
        std::size_t pos1{ file_path.find_last_of('/') };
        std::size_t pos2{ file_path.find_last_of('\\') };
        std::size_t pos{ std::max(pos1, pos2) };
        if (pos != std::string::npos)
        {
            return file_path.substr(0, pos + 1);  // 包含最後的分隔符
        }
        return "";
    }
    
    static std::string GetFileName(const std::string & file_path)
    {
        std::size_t pos1{ file_path.find_last_of('/') };
        std::size_t pos2{ file_path.find_last_of('\\') };
        std::size_t pos{ std::max(pos1, pos2) };
        if (pos != std::string::npos && pos + 1 < file_path.length())
        {
            return file_path.substr(pos + 1);
        }
        return file_path;  // 如果沒有分隔符，整個字串就是檔案名稱
    }

    static std::string EnsureTrailingSlash(const std::string & path)
    {
        if (!path.empty() && path.back() != '/')
        {
            return path + "/";
        }
        return path;
    }
};