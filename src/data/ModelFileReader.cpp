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
    m_file_path{ filename },
    m_file_object{ nullptr },
    m_has_loaded_data{ false }
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
    m_has_loaded_data = false;
    ReadHeader();
    ReadDataArray();
    m_has_loaded_data = true;
}

void ModelFileReader::ReadHeader()
{
    m_file_object->LoadHeader(m_file_path);
    m_file_object->PrintHeader();
}

void ModelFileReader::ReadDataArray()
{
    m_file_object->LoadDataArray(m_file_path);
}

AtomicModelDataBlock * ModelFileReader::GetDataBlockPtr()
{
    if (m_file_object == nullptr)
    {
        throw std::runtime_error("ModelFileReader::GetDataBlockPtr(): file backend is not initialized.");
    }
    if (!m_has_loaded_data)
    {
        throw std::runtime_error("ModelFileReader::GetDataBlockPtr(): file data has not been loaded.");
    }
    return m_file_object->GetDataBlockPtr();
}

} // namespace rhbm_gem
