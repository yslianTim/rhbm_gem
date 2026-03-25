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

const FileFormatDescriptor& ResolveDescriptor(
    const std::filesystem::path& filename,
    bool for_write) {
    const auto extension{FilePathHelper::GetExtension(filename)};
    const auto& descriptor{ResolveFormatDescriptor(extension)};
    const bool supported{for_write ? descriptor.supports_write : descriptor.supports_read};
    if (!supported) {
        throw std::runtime_error(
            "Unsupported file format for " + std::string(for_write ? "write" : "read") +
            ": " + extension);
    }
    return descriptor;
}

template <typename StreamType>
StreamType OpenBinaryFile(const std::filesystem::path& filename, std::ios::openmode mode) {
    StreamType stream{filename, mode};
    if (!stream) {
        throw std::runtime_error("Cannot open the file: " + filename.string());
    }
    return stream;
}

std::unique_ptr<ModelObject> ReadModelWithBackend(
    std::istream& file,
    const std::string& filename,
    ModelFormatBackend backend) {
    switch (backend) {
    case ModelFormatBackend::Pdb: {
        PdbFormat codec;
        return codec.ReadModel(file, filename);
    }
    case ModelFormatBackend::Cif: {
        CifFormat codec;
        return codec.ReadModel(file, filename);
    }
    }
    throw std::runtime_error("Unsupported model file format backend.");
}

void WriteModelWithBackend(
    std::ostream& outfile,
    const ModelObject& model_object,
    int model_parameter,
    ModelFormatBackend backend) {
    switch (backend) {
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

std::unique_ptr<MapObject> ReadMapWithBackend(
    std::istream& file,
    const std::string& filename,
    MapFormatBackend backend) {
    switch (backend) {
    case MapFormatBackend::Mrc: {
        MrcFormat codec;
        return codec.ReadMap(file, filename);
    }
    case MapFormatBackend::Ccp4: {
        CCP4Format codec;
        return codec.ReadMap(file, filename);
    }
    }
    throw std::runtime_error("Unsupported map file format backend.");
}

void WriteMapWithBackend(
    std::ostream& outfile,
    const MapObject& map_object,
    MapFormatBackend backend) {
    switch (backend) {
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

std::unique_ptr<ModelObject> ReadModelWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor) {
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value()) {
        throw std::runtime_error("Unsupported model file format.");
    }

    auto file{OpenBinaryFile<std::ifstream>(filename, std::ios::binary)};
    return ReadModelWithBackend(file, filename.string(), *descriptor.model_backend);
}

void WriteModelWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor,
    const ModelObject& model_object,
    int model_parameter) {
    if (descriptor.kind != DataObjectKind::Model || !descriptor.model_backend.has_value()) {
        throw std::runtime_error("Unsupported model file format.");
    }

    auto outfile{OpenBinaryFile<std::ofstream>(filename, std::ios::binary)};
    WriteModelWithBackend(outfile, model_object, model_parameter, *descriptor.model_backend);
}

std::unique_ptr<MapObject> ReadMapWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor) {
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value()) {
        throw std::runtime_error("Unsupported map file format.");
    }

    auto file{OpenBinaryFile<std::ifstream>(filename, std::ios::binary)};
    return ReadMapWithBackend(file, filename.string(), *descriptor.map_backend);
}

void WriteMapWithDescriptor(
    const std::filesystem::path& filename,
    const FileFormatDescriptor& descriptor,
    const MapObject& map_object) {
    if (descriptor.kind != DataObjectKind::Map || !descriptor.map_backend.has_value()) {
        throw std::runtime_error("Unsupported map file format.");
    }

    auto outfile{OpenBinaryFile<std::ofstream>(filename, std::ios::binary | std::ios::trunc)};
    WriteMapWithBackend(outfile, map_object, *descriptor.map_backend);
}

} // namespace

std::unique_ptr<DataObjectBase> ReadDataObject(const std::filesystem::path& filename) {
    try {
        const auto& descriptor{ResolveDescriptor(filename, false)};
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
        const auto& descriptor{ResolveDescriptor(filename, true)};
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
        const auto& descriptor{ResolveDescriptor(filename, false)};
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
        const auto& descriptor{ResolveDescriptor(filename, true)};
        WriteModelWithDescriptor(filename, descriptor, model_object, model_parameter);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to write model file '" + filename.string() + "': " + ex.what());
    }
}

std::unique_ptr<MapObject> ReadMap(const std::filesystem::path& filename) {
    try {
        const auto& descriptor{ResolveDescriptor(filename, false)};
        return ReadMapWithDescriptor(filename, descriptor);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to read map file '" + filename.string() + "': " + ex.what());
    }
}

void WriteMap(const std::filesystem::path& filename, const MapObject& map_object) {
    try {
        const auto& descriptor{ResolveDescriptor(filename, true)};
        WriteMapWithDescriptor(filename, descriptor, map_object);
    } catch (const std::exception& ex) {
        throw std::runtime_error(
            "Failed to write map file '" + filename.string() + "': " + ex.what());
    }
}

} // namespace rhbm_gem
