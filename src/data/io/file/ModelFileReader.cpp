#include <rhbm_gem/data/io/ModelFileReader.hpp>
#include "internal/io/file/FileFormatBackendFactory.hpp"
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include "internal/io/file/FileFormatRegistry.hpp"
#include "internal/io/file/ModelFileFormatBase.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/AtomicModelDataBlock.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <fstream>

namespace rhbm_gem {

ModelFileReader::ModelFileReader(const std::string & filename) :
    ModelFileReader(filename, BuildDefaultFileFormatRegistry())
{
}

ModelFileReader::ModelFileReader(
    const std::string & filename,
    const FileFormatRegistry & file_format_registry) :
    m_file_path{ filename },
    m_file_object{ nullptr },
    m_has_loaded_data{ false }
{
    const auto & descriptor{
        file_format_registry.LookupForRead(FilePathHelper::GetExtension(filename)) };
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
    std::ifstream file{ m_file_path, std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + m_file_path);
    }
    m_file_object->Read(file, m_file_path);
    m_has_loaded_data = true;
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
