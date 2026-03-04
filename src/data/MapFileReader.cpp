#include "MapFileReader.hpp"
#include "FileFormatBackendFactory.hpp"
#include "FilePathHelper.hpp"
#include "FileFormatRegistry.hpp"
#include "MapFileFormatBase.hpp"
#include "Logger.hpp"

#include <stdexcept>
namespace rhbm_gem {

MapFileReader::MapFileReader(const std::string & filename) :
    m_successfully_read_file{ false }, m_file_path{ filename }
{
    const auto & descriptor{
        FileFormatRegistry::Instance().LookupForRead(FilePathHelper::GetExtension(filename)) };
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value())
    {
        throw std::runtime_error("Unsupported map file format");
    }
    m_file_format_helper = CreateMapFileFormatBackend(*descriptor.map_backend);
}

void MapFileReader::Read()
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
        Logger::Log(LogLevel::Error,
            "MapFileReader::ReadHeader : "
            "Failed to read header from file '" + m_file_path + "' -> " +
            ex.what()
        );
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
        Logger::Log(LogLevel::Error,
            "MapFileReader::ReadMapValueArray : "
            "Failed to read map value array from file '" + m_file_path + "' -> " +
            ex.what()
        );
        return false;
    }
}

std::unique_ptr<float[]> MapFileReader::GetMapValueArray()
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false)
    {
        Logger::Log(LogLevel::Error,
            "MapFileReader::GetMapValueArray : "
            "no data available or read failed."
        );
        return nullptr;
    }
    return m_file_format_helper->GetDataArray();
}

std::array<int, 3> MapFileReader::GetGridSizeArray() const
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false)
    {
        Logger::Log(LogLevel::Error,
            "MapFileReader::GetGridSizeArray : "
            "no grid size available or read failed."
        );
        return {};
    }
    return m_file_format_helper->GetGridSize();
}

std::array<float, 3> MapFileReader::GetGridSpacingArray() const
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false)
    {
        Logger::Log(LogLevel::Error,
            "MapFileReader::GetGridSpacingArray : "
            "no grid spacing available or read failed."
        );
        return {};
    }
    return m_file_format_helper->GetGridSpacing();
}

std::array<float, 3> MapFileReader::GetOriginArray() const
{
    if (m_file_format_helper == nullptr || m_successfully_read_file == false)
    {
        Logger::Log(LogLevel::Error,
            "MapFileReader::GetOriginArray : "
            "no origin available or read failed."
        );
        return {};
    }
    return m_file_format_helper->GetOrigin();
}

} // namespace rhbm_gem
