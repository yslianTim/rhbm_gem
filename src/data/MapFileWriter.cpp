#include "MapFileWriter.hpp"
#include "FilePathHelper.hpp"
#include "MapFileFormatBase.hpp"
#include "MrcFormat.hpp"
#include "CCP4Format.hpp"
#include "MapObject.hpp"

#include <iostream>

MapFileWriter::MapFileWriter(const std::string & filename, const MapObject * map_object) :
    m_file_path{ filename }, m_map_object{ map_object }
{
    auto file_extension{ FilePathHelper::GetExtension(filename) };
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

MapFileWriter::~MapFileWriter()
{

}

void MapFileWriter::Write(void)
{
    WriteHeader();
    WriteMapValueArray();
}

void MapFileWriter::WriteHeader(void)
{
    try
    {
        m_file_format_helper->SetGridSize(m_map_object->GetGridSize());
        m_file_format_helper->SetGridSpacing(m_map_object->GetGridSpacing());
        m_file_format_helper->SetOrigin(m_map_object->GetOrigin());
        m_file_format_helper->SaveHeader(m_file_path);
        m_file_format_helper->PrintHeader();
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

void MapFileWriter::WriteMapValueArray(void)
{
    try
    {
        m_file_format_helper->SetDataArray(
            m_map_object->GetMapValueArraySize(),
            m_map_object->GetMapValueArray());
        m_file_format_helper->SaveDataArray(m_file_path);
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}