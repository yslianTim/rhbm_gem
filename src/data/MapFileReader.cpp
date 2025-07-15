#include "MapFileReader.hpp"
#include "FilePathHelper.hpp"
#include "MapFileFormatBase.hpp"
#include "MrcFormat.hpp"
#include "CCP4Format.hpp"
#include "Logger.hpp"

#include <stdexcept>
#include <algorithm>
#include <cctype>

MapFileReader::MapFileReader(const std::string & filename) :
    m_successfully_read_file{ false }, m_file_path{ filename }
{
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    std::transform(file_extension.begin(), file_extension.end(),
                   file_extension.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if      (file_extension == ".mrc")
    {
        m_file_format_helper = std::make_unique<MrcFormat>();
    }
    else if (file_extension == ".map")
    {
        m_file_format_helper = std::make_unique<CCP4Format>();
    }
    else
    {
        throw std::runtime_error("Unsupported file format");
    }
}

void MapFileReader::Read(void)
{
    std::ifstream file{ m_file_path, std::ios::binary };
    if (!file)
    {
        Logger::Log(LogLevel::Error, "Cannot open the file: " + m_file_path);
        m_successfully_read_file = false;
        return;
    }
    bool header_success{ ReadHeader(file) };
    bool data_success{ false };
    if (header_success == true)
    {
        data_success = ReadMapValueArray(file);
    }
    m_successfully_read_file = header_success && data_success;
}

bool MapFileReader::ReadHeader(std::ifstream & stream)
{
    try
    {
        m_file_format_helper->LoadHeader(stream);
        return true;
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, ex.what());
        return false;
    }
}

bool MapFileReader::ReadMapValueArray(std::ifstream & stream)
{
    try
    {
        m_file_format_helper->LoadDataArray(stream);
        return true;
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, ex.what());
        return false;
    }
}

std::unique_ptr<float[]> MapFileReader::GetMapValueArray(void)
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false) return nullptr;
    return m_file_format_helper->GetDataArray();
}

std::array<int, 3> MapFileReader::GetGridSizeArray(void) const
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false) return {};
    return m_file_format_helper->GetGridSize();
}

std::array<float, 3> MapFileReader::GetGridSpacingArray(void) const
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false) return {};
    return m_file_format_helper->GetGridSpacing();
}

std::array<float, 3> MapFileReader::GetOriginArray(void) const
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false) return {};
    return m_file_format_helper->GetOrigin();
}
