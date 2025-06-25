#include "FileProcessFactoryBase.hpp"
#include "MapFileReader.hpp"
#include "MapFileWriter.hpp"
#include "MapObject.hpp"

#include <stdexcept>

std::unique_ptr<DataObjectBase> MapObjectFactory::CreateDataObject(const std::string & filename)
{
    auto file_reader{ std::make_unique<MapFileReader>(filename) };
    file_reader->Read();
    if (file_reader->IsSuccessfullyRead() == false) return nullptr;
    return std::make_unique<MapObject>(
        file_reader->GetGridSizeArray(),
        file_reader->GetGridSpacingArray(),
        file_reader->GetOriginArray(),
        file_reader->GetMapValueArray()
    );
}

void MapObjectFactory::OutputDataObject(const std::string & filename, DataObjectBase * data_object)
{
    auto map_object{ dynamic_cast<MapObject *>(data_object) };
    if (map_object == nullptr)
    {
        throw std::runtime_error(
            "MapObjectFactory::OutputDataObject(): invalid data_object type");
    }
    auto file_writer{ std::make_unique<MapFileWriter>(filename, map_object) };
    file_writer->Write();
}
