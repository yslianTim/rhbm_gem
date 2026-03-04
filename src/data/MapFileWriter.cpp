#include "MapFileWriter.hpp"
#include "FileFormatBackendFactory.hpp"
#include "FilePathHelper.hpp"
#include "FileFormatRegistry.hpp"
#include "MapFileFormatBase.hpp"
#include "MapObject.hpp"
#include "Logger.hpp"

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
        Logger::Log(LogLevel::Error, "Cannot open the file: " + m_file_path);
        return;
    }
    WriteHeader(file);
    WriteMapValueArray(file);
}

void MapFileWriter::WriteHeader(std::ostream & stream)
{
    try
    {
        m_file_format_helper->SetHeader(
            m_map_object->GetGridSize(),
            m_map_object->GetGridSpacing(),
            m_map_object->GetOrigin()
        );
        m_file_format_helper->SaveHeader(stream);
        m_file_format_helper->PrintHeader();
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, ex.what());
    }
}

void MapFileWriter::WriteMapValueArray(std::ostream & stream)
{
    try
    {
        m_file_format_helper->SaveDataArray(
            m_map_object->GetMapValueArray(),
            m_map_object->GetMapValueArraySize(),
            stream);
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, ex.what());
    }
}

} // namespace rhbm_gem
