#include <rhbm_gem/data/ModelFileWriter.hpp>
#include "internal/FileFormatBackendFactory.hpp"
#include <rhbm_gem/utils/FilePathHelper.hpp>
#include "internal/FileFormatRegistry.hpp"
#include "internal/ModelFileFormatBase.hpp"
#include <rhbm_gem/data/ModelObject.hpp>

#include <stdexcept>
#include <fstream>

namespace rhbm_gem {

ModelFileWriter::ModelFileWriter(
    const std::string & filename, const ModelObject * model_object, int model_par) :
    m_file_path{ filename }, m_model_object{ model_object }, m_model_par{ model_par }
{
    const auto & descriptor{
        FileFormatRegistry::Instance().LookupForWrite(FilePathHelper::GetExtension(filename)) };
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value())
    {
        throw std::runtime_error("Unsupported model file format");
    }
    m_file_object = CreateModelFileFormatBackend(*descriptor.model_backend);
}

ModelFileWriter::~ModelFileWriter()
{
}

void ModelFileWriter::Write()
{
    if (m_model_object == nullptr)
    {
        throw std::invalid_argument("ModelFileWriter::Write(): model_object is null");
    }
    if (m_file_object == nullptr)
    {
        throw std::runtime_error("ModelFileWriter::Write(): file backend is not initialized");
    }
    std::ofstream outfile{ m_file_path, std::ios::binary };
    if (!outfile)
    {
        throw std::runtime_error("Cannot open the file: " + m_file_path);
    }
    m_file_object->Write(*m_model_object, outfile, m_model_par);
}

} // namespace rhbm_gem
