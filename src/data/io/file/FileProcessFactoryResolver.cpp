#include "internal/io/file/FileProcessFactoryResolver.hpp"

#include "internal/io/file/FileFormatRegistry.hpp"
#include "internal/io/file/FileProcessFactoryBase.hpp"
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include <stdexcept>

namespace
{
    std::string NormalizeExtension(std::string extension)
    {
        StringHelper::ToLowerCase(extension);
        return extension;
    }
}

namespace rhbm_gem {

DefaultFileProcessFactoryResolver::DefaultFileProcessFactoryResolver(
    const FileFormatRegistry & file_format_registry) :
    m_file_format_registry{ file_format_registry }
{
}

std::unique_ptr<FileProcessFactoryBase> DefaultFileProcessFactoryResolver::CreateFactory(
    const std::string & extension) const
{
    const auto normalized_extension{ NormalizeExtension(extension) };
    const auto & descriptor{ m_file_format_registry.Lookup(normalized_extension) };
    switch (descriptor.kind)
    {
    case DataObjectKind::Model:
        return std::make_unique<ModelObjectFactory>(m_file_format_registry);
    case DataObjectKind::Map:
        return std::make_unique<MapObjectFactory>(m_file_format_registry);
    }

    throw std::runtime_error("Unsupported file format: " + extension);
}

OverrideableFileProcessFactoryResolver::OverrideableFileProcessFactoryResolver(
    std::shared_ptr<const FileProcessFactoryResolver> fallback_resolver) :
    m_fallback_resolver{ std::move(fallback_resolver) }
{
    if (!m_fallback_resolver)
    {
        throw std::runtime_error(
            "OverrideableFileProcessFactoryResolver requires a fallback resolver.");
    }
}

void OverrideableFileProcessFactoryResolver::RegisterFactory(
    const std::string & extension,
    std::function<std::unique_ptr<FileProcessFactoryBase>()> creator)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_factory_map[NormalizeExtension(extension)] = std::move(creator);
}

void OverrideableFileProcessFactoryResolver::UnregisterFactory(const std::string & extension)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_factory_map.erase(NormalizeExtension(extension));
}

std::unique_ptr<FileProcessFactoryBase> OverrideableFileProcessFactoryResolver::CreateFactory(
    const std::string & extension) const
{
    const auto normalized_extension{ NormalizeExtension(extension) };
    std::function<std::unique_ptr<FileProcessFactoryBase>()> creator;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        const auto iter{ m_factory_map.find(normalized_extension) };
        if (iter != m_factory_map.end())
        {
            creator = iter->second;
        }
    }
    if (creator)
    {
        return creator();
    }

    return m_fallback_resolver->CreateFactory(normalized_extension);
}

} // namespace rhbm_gem
