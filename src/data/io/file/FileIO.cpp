#include <rhbm_gem/data/io/FileIO.hpp>

#include "internal/file/CCP4Format.hpp"
#include "internal/file/CifFormat.hpp"
#include "internal/file/MrcFormat.hpp"
#include "internal/file/PdbFormat.hpp"
#include <rhbm_gem/data/object/DataObjectDispatch.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem {

namespace {

enum class DataObjectKind {
    Model,
    Map
};

enum class ModelFormatBackend {
    Pdb,
    Cif
};

enum class MapFormatBackend {
    Mrc,
    Ccp4
};

struct FileFormatDescriptor {
    std::string extension;
    DataObjectKind kind;
    bool supports_read;
    bool supports_write;
    std::optional<ModelFormatBackend> model_backend;
    std::optional<MapFormatBackend> map_backend;
};

const std::vector<FileFormatDescriptor>& GetFileFormatDescriptors() {
    static const std::vector<FileFormatDescriptor> descriptors{
        {".pdb", DataObjectKind::Model, true, true, ModelFormatBackend::Pdb, std::nullopt},
        {".cif", DataObjectKind::Model, true, true, ModelFormatBackend::Cif, std::nullopt},
        {".mmcif", DataObjectKind::Model, true, false, ModelFormatBackend::Cif, std::nullopt},
        {".mcif", DataObjectKind::Model, true, false, ModelFormatBackend::Cif, std::nullopt},
        {".mrc", DataObjectKind::Map, true, true, std::nullopt, MapFormatBackend::Mrc},
        {".map", DataObjectKind::Map, true, true, std::nullopt, MapFormatBackend::Ccp4},
        {".ccp4", DataObjectKind::Map, true, true, std::nullopt, MapFormatBackend::Ccp4}};
    return descriptors;
}

const FileFormatDescriptor& ResolveFormatDescriptor(const std::string& extension) {
    auto normalized_extension{extension};
    StringHelper::ToLowerCase(normalized_extension);
    for (const auto& descriptor : GetFileFormatDescriptors()) {
        if (descriptor.extension == normalized_extension) {
            return descriptor;
        }
    }
    throw std::runtime_error("Unsupported file format: " + extension);
}

const FileFormatDescriptor& ResolveDescriptorForRead(const std::filesystem::path& filename) {
    const auto& descriptor{ResolveFormatDescriptor(FilePathHelper::GetExtension(filename))};
    if (!descriptor.supports_read) {
        throw std::runtime_error(
            "Unsupported file format for read: " + FilePathHelper::GetExtension(filename));
    }
    return descriptor;
}

const FileFormatDescriptor& ResolveDescriptorForWrite(const std::filesystem::path& filename) {
    const auto& descriptor{ResolveFormatDescriptor(FilePathHelper::GetExtension(filename))};
    if (!descriptor.supports_write) {
        throw std::runtime_error(
            "Unsupported file format for write: " + FilePathHelper::GetExtension(filename));
    }
    return descriptor;
}

std::unique_ptr<ModelObject> ReadModelWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor) {
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value()) {
        throw std::runtime_error("Unsupported model file format.");
    }

    std::ifstream file{filename, std::ios::binary};
    if (!file) {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.model_backend) {
    case ModelFormatBackend::Pdb: {
        PdbFormat codec;
        return codec.ReadModel(file, filename.string());
    }
    case ModelFormatBackend::Cif: {
        CifFormat codec;
        return codec.ReadModel(file, filename.string());
    }
    }
    throw std::runtime_error("Unsupported model file format backend.");
}

void WriteModelWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor,
    const ModelObject& model_object,
    int model_parameter) {
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value()) {
        throw std::runtime_error("Unsupported model file format.");
    }

    std::ofstream outfile{filename, std::ios::binary};
    if (!outfile) {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.model_backend) {
    case ModelFormatBackend::Pdb: {
        PdbFormat codec;
        codec.WriteModel(model_object, outfile, model_parameter);
        return;
    }
    case ModelFormatBackend::Cif: {
        CifFormat codec;
        codec.WriteModel(model_object, outfile, model_parameter);
        return;
    }
    }
    throw std::runtime_error("Unsupported model file format backend.");
}

std::unique_ptr<MapObject> ReadMapWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor) {
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value()) {
        throw std::runtime_error("Unsupported map file format.");
    }

    std::ifstream file{filename, std::ios::binary};
    if (!file) {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.map_backend) {
    case MapFormatBackend::Mrc: {
        MrcFormat codec;
        return codec.ReadMap(file, filename.string());
    }
    case MapFormatBackend::Ccp4: {
        CCP4Format codec;
        return codec.ReadMap(file, filename.string());
    }
    }
    throw std::runtime_error("Unsupported map file format backend.");
}

void WriteMapWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor,
    const MapObject& map_object) {
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value()) {
        throw std::runtime_error("Unsupported map file format.");
    }

    std::ofstream outfile{filename, std::ios::binary | std::ios::trunc};
    if (!outfile) {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    switch (*descriptor.map_backend) {
    case MapFormatBackend::Mrc: {
        MrcFormat codec;
        codec.WriteMap(map_object, outfile);
        return;
    }
    case MapFormatBackend::Ccp4: {
        CCP4Format codec;
        codec.WriteMap(map_object, outfile);
        return;
    }
    }
    throw std::runtime_error("Unsupported map file format backend.");
}

} // namespace

std::unique_ptr<DataObjectBase> ReadDataObject(const std::filesystem::path& filename) {
    try {
        const auto& descriptor{ResolveDescriptorForRead(filename)};
        switch (descriptor.kind) {
        case DataObjectKind::Model:
            return ReadModelWithDescriptor(filename, descriptor);
        case DataObjectKind::Map:
            return ReadMapWithDescriptor(filename, descriptor);
        }
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to read data object from '" + filename.string() + "': " + ex.what());
    }
    throw std::runtime_error("Failed to read data object from '" + filename.string() + "'.");
}

void WriteDataObject(const std::filesystem::path& filename, const DataObjectBase& data_object) {
    try {
        const auto& descriptor{ResolveDescriptorForWrite(filename)};
        switch (descriptor.kind) {
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
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to write data object to '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<ModelObject> ReadModel(const std::filesystem::path& filename) {
    try {
        const auto& descriptor{ResolveDescriptorForRead(filename)};
        return ReadModelWithDescriptor(filename, descriptor);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to read model file '" + filename.string() + "': " + ex.what());
    }
}

void WriteModel(
    const std::filesystem::path& filename,
    const ModelObject& model_object,
    int model_parameter) {
    try {
        const auto& descriptor{ResolveDescriptorForWrite(filename)};
        WriteModelWithDescriptor(filename, descriptor, model_object, model_parameter);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to write model file '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path& filename) {
    try {
        const auto& descriptor{ResolveDescriptorForRead(filename)};
        return ReadMapWithDescriptor(filename, descriptor);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to read map file '" + filename.string() + "': " + ex.what());
    }
}

void WriteMap(const std::filesystem::path& filename, const MapObject& map_object) {
    try {
        const auto& descriptor{ResolveDescriptorForWrite(filename)};
        WriteMapWithDescriptor(filename, descriptor, map_object);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to write map file '" + filename.string() + "': " + ex.what());
    }
}

} // namespace rhbm_gem
