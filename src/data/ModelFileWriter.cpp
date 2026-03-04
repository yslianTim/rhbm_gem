#include "ModelFileWriter.hpp"
#include "FilePathHelper.hpp"
#include "FileFormatRegistry.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "ModelObject.hpp"
#include "Logger.hpp"

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

    if (descriptor.model_backend == ModelFormatBackend::Pdb)
    {
        m_file_object = std::make_unique<PdbFormat>();
    }
    else if (descriptor.model_backend == ModelFormatBackend::Cif)
    {
        m_file_object = std::make_unique<CifFormat>();
    }
    else
    {
        throw std::runtime_error("Unsupported file format");
    }
}

ModelFileWriter::~ModelFileWriter()
{
}

void ModelFileWriter::Write()
{
    if (m_model_object == nullptr || m_file_object == nullptr) return;
    std::ofstream outfile{ m_file_path, std::ios::binary };
    if (!outfile)
    {
        Logger::Log(LogLevel::Error, "Cannot open the file: " + m_file_path);
        return;
    }
    try
    {
        m_file_object->SaveHeader(m_model_object, outfile);
        m_file_object->SaveDataArray(m_model_object, outfile, m_model_par);
    }
    catch (const std::exception & ex)
    {
        Logger::Log(LogLevel::Error,
            "ModelFileWriter::Write : " + std::string(ex.what()));
    }
}

} // namespace rhbm_gem
