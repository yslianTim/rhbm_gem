#pragma once

#include <optional>
#include <string>
#include <vector>

namespace rhbm_gem {

enum class DataObjectKind
{
    Model,
    Map
};

enum class FileFormatOperation
{
    Read,
    Write
};

enum class ModelFormatBackend
{
    Pdb,
    Cif
};

enum class MapFormatBackend
{
    Mrc,
    Ccp4
};

struct FileFormatDescriptor
{
    std::string extension;
    DataObjectKind kind;
    bool supports_read;
    bool supports_write;
    std::optional<ModelFormatBackend> model_backend;
    std::optional<MapFormatBackend> map_backend;
};

class FileFormatRegistry
{
public:
    static const FileFormatRegistry & Instance();

    const FileFormatDescriptor & Lookup(const std::string & extension) const;
    const FileFormatDescriptor & LookupForRead(const std::string & extension) const;
    const FileFormatDescriptor & LookupForWrite(const std::string & extension) const;
    const std::vector<FileFormatDescriptor> & GetAllDescriptors() const;

private:
    FileFormatRegistry() = default;
};

} // namespace rhbm_gem
