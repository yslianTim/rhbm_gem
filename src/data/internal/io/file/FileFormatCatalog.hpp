#pragma once

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

enum class DataObjectKind
{
    Model,
    Map
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

class FileFormatCatalog
{
public:
    explicit FileFormatCatalog(std::vector<FileFormatDescriptor> descriptors);

    const FileFormatDescriptor & Lookup(const std::string & extension) const;
    const FileFormatDescriptor & LookupForRead(const std::string & extension) const;
    const FileFormatDescriptor & LookupForWrite(const std::string & extension) const;
    const std::vector<FileFormatDescriptor> & GetAllDescriptors() const;

private:
    std::vector<FileFormatDescriptor> m_descriptors;
    std::unordered_map<std::string, std::size_t> m_descriptor_index_map;
};

FileFormatCatalog BuildDefaultFileFormatCatalog();

} // namespace rhbm_gem
