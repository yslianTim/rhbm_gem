#include "FileProcessFactoryBase.hpp"
#include "ModelFileReader.hpp"
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
    return std::move(model_object);
}