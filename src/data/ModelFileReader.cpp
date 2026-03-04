#include "ModelFileReader.hpp"
#include "FileFormatBackendFactory.hpp"
#include "FilePathHelper.hpp"
#include "FileFormatRegistry.hpp"
#include "ModelFileFormatBase.hpp"
#include "AtomObject.hpp"
#include "AtomicModelDataBlock.hpp"
#include "Logger.hpp"

namespace rhbm_gem {

ModelFileReader::ModelFileReader(const std::string & filename) :
    m_successfully_read_file{ false }, m_file_path{ filename }
{
    const auto & descriptor{
        FileFormatRegistry::Instance().LookupForRead(FilePathHelper::GetExtension(filename)) };
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value())
    {
        throw std::runtime_error("Unsupported model file format");
    }
    m_file_object = CreateModelFileFormatBackend(*descriptor.model_backend);
}

ModelFileReader::~ModelFileReader()
{
}

void ModelFileReader::Read()
{
    ReadHeader();
    if (m_successfully_read_file == false)
    {
        Logger::Log(LogLevel::Error,
            "ModelFileReader::Read() header parsing failed; skip reading data array.");
        return;
    }
    ReadDataArray();
}

void ModelFileReader::ReadHeader()
{
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

void ModelFileReader::ReadDataArray()
{
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

AtomicModelDataBlock * ModelFileReader::GetDataBlockPtr()
{
    if (m_file_object == nullptr || m_successfully_read_file == false) return nullptr;
    return m_file_object->GetDataBlockPtr();
}

} // namespace rhbm_gem
