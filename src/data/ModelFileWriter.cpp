#include "ModelFileWriter.hpp"
#include "FilePathHelper.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "ModelObject.hpp"
#include "Logger.hpp"

#include <stdexcept>
#include <fstream>

ModelFileWriter::ModelFileWriter(
    const std::string & filename, const ModelObject * model_object, int model_par) :
    m_file_path{ filename }, m_model_object{ model_object }, m_model_par{ model_par }
{
    Logger::Log(LogLevel::Debug, "ModelFileWriter::ModelFileWriter() called");
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

ModelFileWriter::~ModelFileWriter()
{
    Logger::Log(LogLevel::Debug, "ModelFileWriter::~ModelFileWriter() called");
}

void ModelFileWriter::Write(void)
{
    Logger::Log(LogLevel::Debug, "ModelFileWriter::Write() called");
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
