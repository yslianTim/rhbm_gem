#include "MapFileWriter.hpp"
#include "FilePathHelper.hpp"
#include "MapFileFormatBase.hpp"
#include "MrcFormat.hpp"
#include "CCP4Format.hpp"
#include "MapObject.hpp"
#include "Logger.hpp"

#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

MapFileWriter::MapFileWriter(const std::string & filename, const MapObject * map_object) :
    m_file_path{ filename }, m_map_object{ map_object }
{
    if (m_map_object == nullptr)
    {
        throw std::invalid_argument("MapFileWriter: map_object is null");
    }
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

void MapFileWriter::Write(void)
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
        m_file_format_helper->SetGridSize(m_map_object->GetGridSize());
        m_file_format_helper->SetGridSpacing(m_map_object->GetGridSpacing());
        m_file_format_helper->SetOrigin(m_map_object->GetOrigin());
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
