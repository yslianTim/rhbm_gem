#include "internal/io/file/FileProcessFactoryBase.hpp"
#include "internal/io/file/FileFormatRegistry.hpp"
#include <rhbm_gem/data/io/MapFileReader.hpp>
#include <rhbm_gem/data/io/MapFileWriter.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/dispatch/DataObjectDispatch.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <stdexcept>

namespace rhbm_gem {

MapObjectFactory::MapObjectFactory(const FileFormatRegistry & file_format_registry) :
    m_file_format_registry{ file_format_registry }
{
}

std::unique_ptr<DataObjectBase> MapObjectFactory::CreateDataObject(const std::string & filename)
{
    auto file_reader{ std::make_unique<MapFileReader>(filename, m_file_format_registry) };
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
    auto file_writer{
        std::make_unique<MapFileWriter>(filename, &map_object, m_file_format_registry) };
    file_writer->Write();
}

} // namespace rhbm_gem
