#include "FileProcessFactoryBase.hpp"
#include "MapFileReader.hpp"
#include "MapFileWriter.hpp"
#include "MapObject.hpp"
#include "DataObjectDispatch.hpp"
#include "Logger.hpp"

#include <stdexcept>

namespace rhbm_gem {

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

void MapObjectFactory::OutputDataObject(
    const std::string & filename,
    const DataObjectBase & data_object)
{
    const auto & map_object{
        ExpectMapObject(data_object, "MapObjectFactory::OutputDataObject()") };
    auto file_writer{ std::make_unique<MapFileWriter>(filename, &map_object) };
    file_writer->Write();
}

} // namespace rhbm_gem
