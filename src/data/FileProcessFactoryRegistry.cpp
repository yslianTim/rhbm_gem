#include "FileProcessFactoryRegistry.hpp"
#include "FileProcessFactoryBase.hpp"
#include "FileFormatRegistry.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <stdexcept>

namespace rhbm_gem {

FileProcessFactoryRegistry & FileProcessFactoryRegistry::Instance()
{
    static FileProcessFactoryRegistry instance;
    return instance;
}

void FileProcessFactoryRegistry::RegisterDefaultFactories()
{
    m_factory_map.clear();
    for (const auto & descriptor : FileFormatRegistry::Instance().GetAllDescriptors())
    {
        switch (descriptor.kind)
        {
        case DataObjectKind::Model:
            RegisterFactory(descriptor.extension, []() { return std::make_unique<ModelObjectFactory>(); });
            break;
        case DataObjectKind::Map:
            RegisterFactory(descriptor.extension, []() { return std::make_unique<MapObjectFactory>(); });
            break;
        }
    }
}

void FileProcessFactoryRegistry::RegisterFactory(
    const std::string & extension, std::function<std::unique_ptr<FileProcessFactoryBase>()> creator)
{
    auto normalized_extension{ extension };
    StringHelper::ToLowerCase(normalized_extension);
    m_factory_map[normalized_extension] = std::move(creator);
}

std::unique_ptr<FileProcessFactoryBase> FileProcessFactoryRegistry::CreateFactory(
    const std::string & extension) const
{
    auto normalized_extension{ extension };
    StringHelper::ToLowerCase(normalized_extension);
    auto iter{ m_factory_map.find(normalized_extension) };
    if (iter != m_factory_map.end())
    {
        return (iter->second)();
    }

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

} // namespace rhbm_gem
