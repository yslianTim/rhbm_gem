#include <rhbm_gem/data/io/FileIO.hpp>

#include "internal/file/CCP4Format.hpp"
#include "internal/file/CifFormat.hpp"
#include "internal/file/FileFormatCatalog.hpp"
#include "internal/file/MrcFormat.hpp"
#include "internal/file/PdbFormat.hpp"
#include <rhbm_gem/data/object/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <fstream>
#include <stdexcept>

namespace rhbm_gem {

namespace {

const FileFormatCatalog & DefaultFileFormatCatalog()
{
    static const FileFormatCatalog catalog{ BuildDefaultFileFormatCatalog() };
    return catalog;
}

const FileFormatDescriptor & ResolveDescriptorForRead(const std::filesystem::path & filename)
{
    return DefaultFileFormatCatalog().LookupForRead(FilePathHelper::GetExtension(filename));
}

const FileFormatDescriptor & ResolveDescriptorForWrite(const std::filesystem::path & filename)
{
    return DefaultFileFormatCatalog().LookupForWrite(FilePathHelper::GetExtension(filename));
}

std::unique_ptr<ModelObject> ReadModelWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor)
{
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value())
    {
        throw std::runtime_error("Unsupported model file format.");
    }

    std::ifstream file{ filename, std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.model_backend)
    {
    case ModelFormatBackend::Pdb:
    {
        PdbFormat codec;
        return codec.ReadModel(file, filename.string());
    }
    case ModelFormatBackend::Cif:
    {
        CifFormat codec;
        return codec.ReadModel(file, filename.string());
    }
    }
    throw std::runtime_error("Unsupported model file format backend.");
}

void WriteModelWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor,
    const ModelObject & model_object,
    int model_parameter)
{
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value())
    {
        throw std::runtime_error("Unsupported model file format.");
    }

    std::ofstream outfile{ filename, std::ios::binary };
    if (!outfile)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.model_backend)
    {
    case ModelFormatBackend::Pdb:
    {
        PdbFormat codec;
        codec.WriteModel(model_object, outfile, model_parameter);
        return;
    }
    case ModelFormatBackend::Cif:
    {
        CifFormat codec;
        codec.WriteModel(model_object, outfile, model_parameter);
        return;
    }
    }
    throw std::runtime_error("Unsupported model file format backend.");
}

std::unique_ptr<MapObject> ReadMapWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor)
{
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value())
    {
        throw std::runtime_error("Unsupported map file format.");
    }

    std::ifstream file{ filename, std::ios::binary };
    if (!file)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.map_backend)
    {
    case MapFormatBackend::Mrc:
    {
        MrcFormat codec;
        return codec.ReadMap(file, filename.string());
    }
    case MapFormatBackend::Ccp4:
    {
        CCP4Format codec;
        return codec.ReadMap(file, filename.string());
    }
    }
    throw std::runtime_error("Unsupported map file format backend.");
}

void WriteMapWithDescriptor(
    const std::filesystem::path & filename,
    const FileFormatDescriptor & descriptor,
    const MapObject & map_object)
{
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value())
    {
        throw std::runtime_error("Unsupported map file format.");
    }

    std::ofstream outfile{ filename, std::ios::binary | std::ios::trunc };
    if (!outfile)
    {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.map_backend)
    {
    case MapFormatBackend::Mrc:
    {
        MrcFormat codec;
        codec.WriteMap(map_object, outfile);
        return;
    }
    case MapFormatBackend::Ccp4:
    {
        CCP4Format codec;
        codec.WriteMap(map_object, outfile);
        return;
    }
    }
    throw std::runtime_error("Unsupported map file format backend.");
}

} // namespace

std::unique_ptr<DataObjectBase> ReadDataObject(const std::filesystem::path & filename)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForRead(filename) };
        switch (descriptor.kind)
        {
        case DataObjectKind::Model:
            return ReadModelWithDescriptor(filename, descriptor);
        case DataObjectKind::Map:
            return ReadMapWithDescriptor(filename, descriptor);
        }
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read data object from '" + filename.string() + "': " + ex.what());
    }
    throw std::runtime_error("Failed to read data object from '" + filename.string() + "'.");
}

void WriteDataObject(const std::filesystem::path & filename, const DataObjectBase & data_object)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForWrite(filename) };
        switch (descriptor.kind)
        {
        case DataObjectKind::Model:
            WriteModelWithDescriptor(
                filename,
                descriptor,
                ExpectModelObject(data_object, "WriteDataObject()"),
                0);
            return;
        case DataObjectKind::Map:
            WriteMapWithDescriptor(
                filename,
                descriptor,
                ExpectMapObject(data_object, "WriteDataObject()"));
            return;
        }
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write data object to '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path & filename)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForRead(filename) };
        return ReadModelWithDescriptor(filename, descriptor);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read model file '" + filename.string() + "': " + ex.what());
    }
}

void WriteModel(
    const std::filesystem::path & filename,
    const ModelObject & model_object,
    int model_parameter)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForWrite(filename) };
        WriteModelWithDescriptor(filename, descriptor, model_object, model_parameter);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write model file '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path & filename)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForRead(filename) };
        return ReadMapWithDescriptor(filename, descriptor);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to read map file '" + filename.string() + "': " + ex.what());
    }
}

void WriteMap(const std::filesystem::path & filename, const MapObject & map_object)
{
    try
    {
        const auto & descriptor{ ResolveDescriptorForWrite(filename) };
        WriteMapWithDescriptor(filename, descriptor, map_object);
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to write map file '" + filename.string() + "': " + ex.what());
    }
}

} // namespace rhbm_gem
