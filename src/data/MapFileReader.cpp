#include "MapFileReader.hpp"
#include "FilePathHelper.hpp"
#include "MapFileFormatBase.hpp"
#include "MrcFormat.hpp"
#include "CCP4Format.hpp"

MapFileReader::MapFileReader(const std::string & filename) :
    m_file_path{ filename }
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

MapFileReader::~MapFileReader()
{

}

void MapFileReader::Read(void)
{
    ReadHeader();
    ReadMapValueArray();
}

void MapFileReader::ReadHeader(void)
{
    try
    {
        m_file_format_helper->LoadHeader(m_file_path);
        m_file_format_helper->PrintHeader();
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

void MapFileReader::ReadMapValueArray(void)
{
    try
    {
        m_file_format_helper->LoadDataArray(m_file_path);
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

std::unique_ptr<float[]> MapFileReader::GetMapValueArray(void)
{
    return m_file_format_helper->GetDataArray();
}

std::array<int, 3> MapFileReader::GetGridSizeArray(void)
{
    return m_file_format_helper->GetGridSize();
}

std::array<float, 3> MapFileReader::GetGridSpacingArray(void)
{
    return m_file_format_helper->GetGridSpacing();
}

std::array<float, 3> MapFileReader::GetOriginArray(void)
{
    return m_file_format_helper->GetOrigin();
}