#include "ModelFileWriter.hpp"
#include "FilePathHelper.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "ModelObject.hpp"

#include <iostream>

ModelFileWriter::ModelFileWriter(const std::string & filename, const ModelObject * model_object) :
    m_file_path{ filename }, m_model_object{ model_object }
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

ModelFileWriter::~ModelFileWriter()
{

}

void ModelFileWriter::Write(void)
{
    if (m_model_object == nullptr) return;
}