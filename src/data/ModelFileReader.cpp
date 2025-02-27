#include "ModelFileReader.hpp"
#include "FilePathHelper.hpp"
#include "ModelFileFormatBase.hpp"
#include "PdbFormat.hpp"
#include "CifFormat.hpp"
#include "AtomObject.hpp"

ModelFileReader::ModelFileReader(const std::string & filename) :
    m_file_path{ filename }
{
    auto file_extension{ FilePathHelper::GetExtension(filename) };
    if      (file_extension == ".pdb")
    {
        m_file_format_helper = std::make_unique<PdbFormat>();
    }
    else if (file_extension == ".cif")
    {
        m_file_format_helper = std::make_unique<CifFormat>();
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
        m_file_format_helper->LoadHeader(m_file_path);
        m_file_format_helper->PrintHeader();
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
        m_file_format_helper->LoadDataArray(m_file_path);
    }
    catch (const std::exception & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
}

std::vector<std::unique_ptr<AtomObject>> ModelFileReader::GetAtomObjectList(void)
{
    return m_file_format_helper->GetAtomObjectList();
}

std::string ModelFileReader::GetPdbID(void) const
{
    return m_file_format_helper->GetPdbID();
}

std::string ModelFileReader::GetEmdID(void) const
{
    return m_file_format_helper->GetEmdID();
}
