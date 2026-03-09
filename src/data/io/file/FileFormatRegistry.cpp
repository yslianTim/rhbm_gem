#include "internal/io/file/FileFormatRegistry.hpp"

#include <rhbm_gem/utils/domain/StringHelper.hpp>

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

    const std::unordered_map<std::string, std::size_t> & GetDescriptorIndexMap()
    {
        static const auto k_descriptor_index_map{
            []()
            {
                std::unordered_map<std::string, std::size_t> index_map;
                index_map.reserve(k_descriptors.size());
                for (size_t i = 0; i < k_descriptors.size(); i++)
                {
                    const auto & descriptor{ k_descriptors.at(i) };
                    const auto [iter, inserted]{ index_map.emplace(descriptor.extension, i) };
                    if (!inserted)
                    {
                        throw std::runtime_error(
                            "Duplicate file format descriptor extension: " + descriptor.extension);
                    }
                    if (descriptor.kind == rhbm_gem::DataObjectKind::Model &&
                        !descriptor.model_backend.has_value())
                    {
                        throw std::runtime_error(
                            "Model descriptor must define model backend: " + descriptor.extension);
                    }
                    if (descriptor.kind == rhbm_gem::DataObjectKind::Map &&
                        !descriptor.map_backend.has_value())
                    {
                        throw std::runtime_error(
                            "Map descriptor must define map backend: " + descriptor.extension);
                    }
                }
                return index_map;
            }()
        };
        return k_descriptor_index_map;
    }
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
    const auto & descriptor_index_map{ GetDescriptorIndexMap() };
    const auto iter{ descriptor_index_map.find(normalized_extension) };
    if (iter != descriptor_index_map.end())
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
