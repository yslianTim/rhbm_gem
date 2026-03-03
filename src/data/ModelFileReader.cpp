#include "ModelFileReader.hpp"
#include "FilePathHelper.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "AtomicModelDataBlock.hpp"
#include "Logger.hpp"

namespace rhbm_gem {

ModelFileReader::ModelFileReader(const std::string & filename) :
    m_successfully_read_file{ false }, m_file_path{ filename }
{
    Logger::Log(LogLevel::Debug, "ModelFileReader::ModelFileReader() called");
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    if      (file_extension == ".pdb")
    {
        m_file_object = std::make_unique<PdbFormat>();
    }
    else if (file_extension == ".cif" || file_extension == ".mmcif" || file_extension == ".mcif")
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
    Logger::Log(LogLevel::Debug, "ModelFileReader::~ModelFileReader() called");
}

void ModelFileReader::Read(void)
{
    Logger::Log(LogLevel::Debug, "ModelFileReader::Read() called");
    ReadHeader();
    if (m_successfully_read_file == false)
    {
        Logger::Log(LogLevel::Error,
            "ModelFileReader::Read() header parsing failed; skip reading data array.");
        return;
    }
    ReadDataArray();
}

void ModelFileReader::ReadHeader(void)
{
    Logger::Log(LogLevel::Debug, "ModelFileReader::ReadHeader() called");
    try
    {
        m_file_object->LoadHeader(m_file_path);
        m_file_object->PrintHeader();
        m_successfully_read_file = true;
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, "ModelFileReader::ReadHeader() : " + std::string(ex.what()));
        m_successfully_read_file = false;
    }
}

void ModelFileReader::ReadDataArray(void)
{
    Logger::Log(LogLevel::Debug, "ModelFileReader::ReadDataArray() called");
    try
    {
        m_file_object->LoadDataArray(m_file_path);
        m_successfully_read_file = true;
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error, "ModelFileReader::ReadDataArray()" + std::string(ex.what()));
        m_successfully_read_file = false;
    }
}

AtomicModelDataBlock * ModelFileReader::GetDataBlockPtr(void)
{
    if (m_file_object == nullptr || m_successfully_read_file == false) return nullptr;
    return m_file_object->GetDataBlockPtr();
}

} // namespace rhbm_gem
