#include "FileProcessFactoryRegistry.hpp"
#include "FileProcessFactoryBase.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <stdexcept>

namespace rhbm_gem {

FileProcessFactoryRegistry & FileProcessFactoryRegistry::Instance(void)
{
    static FileProcessFactoryRegistry instance;
    return instance;
}

void FileProcessFactoryRegistry::RegisterDefaultFactories(void)
{
    RegisterFactory(".pdb",  []() { return std::make_unique<ModelObjectFactory>(); });
    RegisterFactory(".cif",  []() { return std::make_unique<ModelObjectFactory>(); });
    RegisterFactory(".mmcif",[]() { return std::make_unique<ModelObjectFactory>(); });
    RegisterFactory(".mcif", []() { return std::make_unique<ModelObjectFactory>(); });
    RegisterFactory(".mrc",  []() { return std::make_unique<MapObjectFactory>(); });
    RegisterFactory(".map",  []() { return std::make_unique<MapObjectFactory>(); });
    RegisterFactory(".ccp4", []() { return std::make_unique<MapObjectFactory>(); });
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
    if (iter == m_factory_map.end())
    {
        throw std::runtime_error("Unsupported file format: " + extension);
    }
    return (iter->second)();
}

} // namespace rhbm_gem
