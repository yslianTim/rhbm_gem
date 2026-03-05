#include "FileFormatRegistry.hpp"

#include "StringHelper.hpp"

#include <unordered_map>
#include <stdexcept>

namespace
{
    const std::vector<rhbm_gem::FileFormatDescriptor> k_descriptors{
        { ".pdb", rhbm_gem::DataObjectKind::Model, true, true, rhbm_gem::ModelFormatBackend::Pdb, std::nullopt },
        { ".cif", rhbm_gem::DataObjectKind::Model, true, true, rhbm_gem::ModelFormatBackend::Cif, std::nullopt },
        { ".mmcif", rhbm_gem::DataObjectKind::Model, true, false, rhbm_gem::ModelFormatBackend::Cif, std::nullopt },
        { ".mcif", rhbm_gem::DataObjectKind::Model, true, false, rhbm_gem::ModelFormatBackend::Cif, std::nullopt },
        { ".mrc", rhbm_gem::DataObjectKind::Map, true, true, std::nullopt, rhbm_gem::MapFormatBackend::Mrc },
        { ".map", rhbm_gem::DataObjectKind::Map, true, true, std::nullopt, rhbm_gem::MapFormatBackend::Ccp4 },
        { ".ccp4", rhbm_gem::DataObjectKind::Map, true, true, std::nullopt, rhbm_gem::MapFormatBackend::Ccp4 }
    };

    const std::unordered_map<std::string, std::size_t> k_descriptor_index_map{
        { ".pdb", 0 },
        { ".cif", 1 },
        { ".mmcif", 2 },
        { ".mcif", 3 },
        { ".mrc", 4 },
        { ".map", 5 },
        { ".ccp4", 6 }
    };
}

namespace rhbm_gem {

const FileFormatRegistry & FileFormatRegistry::Instance()
{
    static FileFormatRegistry instance;
    return instance;
}

const FileFormatDescriptor & FileFormatRegistry::Lookup(const std::string & extension) const
{
    auto normalized_extension{ extension };
    StringHelper::ToLowerCase(normalized_extension);
    const auto iter{ k_descriptor_index_map.find(normalized_extension) };
    if (iter != k_descriptor_index_map.end())
    {
        return k_descriptors.at(iter->second);
    }
    throw std::runtime_error("Unsupported file format: " + extension);
}

const FileFormatDescriptor & FileFormatRegistry::LookupForRead(const std::string & extension) const
{
    const auto & descriptor{ Lookup(extension) };
    if (!descriptor.supports_read)
    {
        throw std::runtime_error("Unsupported file format for read: " + extension);
    }
    return descriptor;
}

const FileFormatDescriptor & FileFormatRegistry::LookupForWrite(const std::string & extension) const
{
    const auto & descriptor{ Lookup(extension) };
    if (!descriptor.supports_write)
    {
        throw std::runtime_error("Unsupported file format for write: " + extension);
    }
    return descriptor;
}

const std::vector<FileFormatDescriptor> & FileFormatRegistry::GetAllDescriptors() const
{
    return k_descriptors;
}

} // namespace rhbm_gem
