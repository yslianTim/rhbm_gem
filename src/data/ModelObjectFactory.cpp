#include "FileProcessFactoryBase.hpp"
#include "ModelFileReader.hpp"
#include "ModelFileWriter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "AtomicModelDataBlock.hpp"
#include "Logger.hpp"

#include <stdexcept>

std::unique_ptr<DataObjectBase> ModelObjectFactory::CreateDataObject(const std::string & filename)
{
    Logger::Log(LogLevel::Debug, "ModelObjectFactory::CreateDataObject() called");
    auto file_reader{ std::make_unique<ModelFileReader>(filename) };
    file_reader->Read();
    if (file_reader->IsSuccessfullyRead() == false) return nullptr;
    auto data_block{ file_reader->GetDataBlockPtr() };
    auto model_object{ std::make_unique<ModelObject>(data_block->GetAtomObjectList()) };
    model_object->SetPdbID(data_block->GetPdbID());
    model_object->SetEmdID(data_block->GetEmdID());
    model_object->SetResolution(data_block->GetResolution());
    model_object->SetResolutionMethod(data_block->GetResolutionMethod());
    model_object->SetChainIDListMap(data_block->GetChainIDListMap());
    model_object->SetChemicalComponentEntryMap(data_block->GetChemicalComponentEntryMap());
    return model_object;
}

void ModelObjectFactory::OutputDataObject(const std::string & filename, DataObjectBase * data_object)
{
    Logger::Log(LogLevel::Debug, "ModelObjectFactory::OutputDataObject() called");
    auto model_object{ dynamic_cast<ModelObject *>(data_object) };
    if (model_object == nullptr)
    {
        throw std::runtime_error(
            "ModelObjectFactory::OutputDataObject(): invalid data_object type");
    }
    auto file_writer{ std::make_unique<ModelFileWriter>(filename, model_object) };
    file_writer->Write();
}
