#include "ModelFileReader.hpp"
#include "FilePathHelper.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "AtomObject.hpp"
#include "AtomicModelDataBlock.hpp"

#include <iostream>

ModelFileReader::ModelFileReader(const std::string & filename) :
    m_file_path{ filename }
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
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

void ModelFileReader::ReadDataArray(void)
{
    try
    {
        m_file_object->LoadDataArray(m_file_path);
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

AtomicModelDataBlock * ModelFileReader::GetDataBlockPtr(void)
{
    return m_file_object->GetDataBlockPtr();
}