#include "FileProcessFactoryBase.hpp"
#include "MapFileReader.hpp"
#include "MapObject.hpp"

#include <stdexcept>

std::unique_ptr<DataObjectBase> MapObjectFactory::CreateDataObject(const std::string & filename)
{
    auto file_reader{ std::make_unique<MapFileReader>(filename) };
    file_reader->Read();
    return std::make_unique<MapObject>(
        file_reader->GetGridSizeArray(),
        file_reader->GetGridSpacingArray(),
        file_reader->GetOriginArray(),
        file_reader->GetMapValueArray()
    );
}