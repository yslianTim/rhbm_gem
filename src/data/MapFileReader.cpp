#include "MapFileReader.hpp"
#include "FileFormatBackendFactory.hpp"
#include "FilePathHelper.hpp"
#include "FileFormatRegistry.hpp"
#include "MapFileFormatBase.hpp"

#include <stdexcept>
namespace rhbm_gem {

MapFileReader::MapFileReader(const std::string & filename) :
    m_file_path{ filename },
    m_file_format_helper{ nullptr },
    m_has_loaded_data{ false }
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
    m_has_loaded_data = false;
    std::ifstream file{ m_file_path, std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + m_file_path);
    }
    if (m_file_format_helper == nullptr)
    {
        throw std::runtime_error("MapFileReader::Read(): file backend is not initialized.");
    }
    m_file_format_helper->Read(file, m_file_path);
    m_has_loaded_data = true;
}

std::unique_ptr<float[]> MapFileReader::GetMapValueArray()
{
    if (m_file_format_helper == nullptr || !m_has_loaded_data)
    {
        throw std::runtime_error("MapFileReader::GetMapValueArray(): no data available.");
    }
    return m_file_format_helper->GetDataArray();
}

std::array<int, 3> MapFileReader::GetGridSizeArray() const
{
    if (m_file_format_helper == nullptr || !m_has_loaded_data)
    {
        throw std::runtime_error("MapFileReader::GetGridSizeArray(): no data available.");
    }
    return m_file_format_helper->GetGridSize();
}

std::array<float, 3> MapFileReader::GetGridSpacingArray() const
{
    if (m_file_format_helper == nullptr || !m_has_loaded_data)
    {
        throw std::runtime_error("MapFileReader::GetGridSpacingArray(): no data available.");
    }
    return m_file_format_helper->GetGridSpacing();
}

std::array<float, 3> MapFileReader::GetOriginArray() const
{
    if (m_file_format_helper == nullptr || !m_has_loaded_data)
    {
        throw std::runtime_error("MapFileReader::GetOriginArray(): no data available.");
    }
    return m_file_format_helper->GetOrigin();
}

} // namespace rhbm_gem
