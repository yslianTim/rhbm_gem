#include "FileProcessFactoryBase.hpp"
#include "ModelFileReader.hpp"
#include "ModelFileWriter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"

#include <stdexcept>

std::unique_ptr<DataObjectBase> ModelObjectFactory::CreateDataObject(const std::string & filename)
{
    auto file_reader{ std::make_unique<ModelFileReader>(filename) };
    file_reader->Read();
    auto model_object{ std::make_unique<ModelObject>(file_reader->GetAtomObjectList()) };
    model_object->SetPdbID(file_reader->GetPdbID());
    model_object->SetEmdID(file_reader->GetEmdID());
    model_object->SetResolution(file_reader->GetResolution());
    model_object->SetResolutionMethod(file_reader->GetResolutionMethod());
    return model_object;
}

void ModelObjectFactory::OutputDataObject(const std::string & filename, DataObjectBase * data_object)
{
    auto model_object{ dynamic_cast<ModelObject *>(data_object) };
    auto file_writer{ std::make_unique<ModelFileWriter>(filename, model_object) };
    file_writer->Write();
}