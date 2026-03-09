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

std::unique_ptr<FileProcessFactoryBase> DefaultFileProcessFactoryResolver::CreateFactory(
    const std::string & extension) const
{
    const auto normalized_extension{ NormalizeExtension(extension) };
    const auto & descriptor{ FileFormatRegistry::Instance().Lookup(normalized_extension) };
    switch (descriptor.kind)
    {
    case DataObjectKind::Model:
        return std::make_unique<ModelObjectFactory>();
    case DataObjectKind::Map:
        return std::make_unique<MapObjectFactory>();
    }

    throw std::runtime_error("Unsupported file format: " + extension);
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

    DefaultFileProcessFactoryResolver default_resolver;
    return default_resolver.CreateFactory(normalized_extension);
}

std::shared_ptr<const FileProcessFactoryResolver> CreateDefaultFileProcessFactoryResolver()
{
    static const auto resolver{
        std::make_shared<DefaultFileProcessFactoryResolver>()
    };
    return resolver;
}

} // namespace rhbm_gem
