#include "MapFileWriter.hpp"
#include "FileFormatBackendFactory.hpp"
#include "FilePathHelper.hpp"
#include "FileFormatRegistry.hpp"
#include "MapFileFormatBase.hpp"
#include "MapObject.hpp"

#include <fstream>
#include <stdexcept>
namespace rhbm_gem {

MapFileWriter::MapFileWriter(const std::string & filename, const MapObject * map_object) :
    m_file_path{ filename }, m_map_object{ map_object }
{
    if (m_map_object == nullptr)
    {
        throw std::invalid_argument("MapFileWriter: map_object is null");
    }
    const auto & descriptor{
        FileFormatRegistry::Instance().LookupForWrite(FilePathHelper::GetExtension(filename)) };
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value())
    {
        throw std::runtime_error("Unsupported map file format");
    }
    m_file_format_helper = CreateMapFileFormatBackend(*descriptor.map_backend);
}

void MapFileWriter::Write()
{
    if (m_map_object == nullptr)
    {
        throw std::invalid_argument("MapFileWriter::Write(): map_object is null");
    }
    std::ofstream file{ m_file_path, std::ios::binary | std::ios::trunc };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + m_file_path);
    }
    if (m_file_format_helper == nullptr)
    {
        throw std::runtime_error("MapFileWriter::Write(): file backend is not initialized.");
    }
    m_file_format_helper->Write(*m_map_object, file);
}

} // namespace rhbm_gem
