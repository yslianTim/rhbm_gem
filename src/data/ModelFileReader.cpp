#include "ModelFileReader.hpp"
#include "FilePathHelper.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "AtomicModelDataBlock.hpp"
#include "Logger.hpp"

ModelFileReader::ModelFileReader(const std::string & filename) :
    m_successfully_read_file{ false }, m_file_path{ filename }
{
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    if      (file_extension == ".pdb")
    {
        m_file_object = std::make_unique<PdbFormat>();
    }
    else if (file_extension == ".cif")
    {
        m_file_object = std::make_unique<CifFormat>();
    }
    else
    {
        throw std::runtime_error("Unsupported file format");
    }
}

ModelFileReader::~ModelFileReader()
{

}

void ModelFileReader::Read(void)
{
    ReadHeader();
    ReadDataArray();
}

void ModelFileReader::ReadHeader(void)
{
    try
    {
        m_file_object->LoadHeader(m_file_path);
        m_file_object->PrintHeader();
        m_successfully_read_file = true;
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, ex.what());
        m_successfully_read_file = false;
    }
}

void ModelFileReader::ReadDataArray(void)
{
    try
    {
        m_file_object->LoadDataArray(m_file_path);
        m_successfully_read_file = true;
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, ex.what());
        m_successfully_read_file = false;
    }
}

AtomicModelDataBlock * ModelFileReader::GetDataBlockPtr(void)
{
    if (m_file_object == nullptr || m_successfully_read_file == false) return nullptr;
    return m_file_object->GetDataBlockPtr();
}
