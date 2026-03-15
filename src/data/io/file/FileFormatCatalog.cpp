#include "internal/io/file/FileFormatCatalog.hpp"

#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include <stdexcept>

namespace rhbm_gem {

namespace {

std::vector<FileFormatDescriptor> BuildDefaultDescriptors()
{
    return {
        { ".pdb", DataObjectKind::Model, true, true, ModelFormatBackend::Pdb, std::nullopt },
        { ".cif", DataObjectKind::Model, true, true, ModelFormatBackend::Cif, std::nullopt },
        { ".mmcif", DataObjectKind::Model, true, false, ModelFormatBackend::Cif, std::nullopt },
        { ".mcif", DataObjectKind::Model, true, false, ModelFormatBackend::Cif, std::nullopt },
        { ".mrc", DataObjectKind::Map, true, true, std::nullopt, MapFormatBackend::Mrc },
        { ".map", DataObjectKind::Map, true, true, std::nullopt, MapFormatBackend::Ccp4 },
        { ".ccp4", DataObjectKind::Map, true, true, std::nullopt, MapFormatBackend::Ccp4 }
    };
}

} // namespace

FileFormatCatalog::FileFormatCatalog(std::vector<FileFormatDescriptor> descriptors) :
    m_descriptors{ std::move(descriptors) },
    m_descriptor_index_map{}
{
    m_descriptor_index_map.reserve(m_descriptors.size());
    for (std::size_t i = 0; i < m_descriptors.size(); ++i)
    {
        const auto & descriptor{ m_descriptors.at(i) };
        const auto [iter, inserted]{ m_descriptor_index_map.emplace(descriptor.extension, i) };
        if (!inserted)
        {
            throw std::runtime_error(
                "Duplicate file format descriptor extension: " + descriptor.extension);
        }
        if (descriptor.kind == DataObjectKind::Model && !descriptor.model_backend.has_value())
        {
            throw std::runtime_error(
                "Model descriptor must define model backend: " + descriptor.extension);
        }
        if (descriptor.kind == DataObjectKind::Map && !descriptor.map_backend.has_value())
        {
            throw std::runtime_error(
                "Map descriptor must define map backend: " + descriptor.extension);
        }
    }
}

const FileFormatDescriptor & FileFormatCatalog::Lookup(const std::string & extension) const
{
    auto normalized_extension{ extension };
    StringHelper::ToLowerCase(normalized_extension);
    const auto iter{ m_descriptor_index_map.find(normalized_extension) };
    if (iter != m_descriptor_index_map.end())
    {
        return m_descriptors.at(iter->second);
    }
    throw std::runtime_error("Unsupported file format: " + extension);
}

const FileFormatDescriptor & FileFormatCatalog::LookupForRead(const std::string & extension) const
{
    const auto & descriptor{ Lookup(extension) };
    if (!descriptor.supports_read)
    {
        throw std::runtime_error("Unsupported file format for read: " + extension);
    }
    return descriptor;
}

const FileFormatDescriptor & FileFormatCatalog::LookupForWrite(const std::string & extension) const
{
    const auto & descriptor{ Lookup(extension) };
    if (!descriptor.supports_write)
    {
        throw std::runtime_error("Unsupported file format for write: " + extension);
    }
    return descriptor;
}

const std::vector<FileFormatDescriptor> & FileFormatCatalog::GetAllDescriptors() const
{
    return m_descriptors;
}

FileFormatCatalog BuildDefaultFileFormatCatalog()
{
    return FileFormatCatalog{ BuildDefaultDescriptors() };
}

} // namespace rhbm_gem
