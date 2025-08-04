#include "FileProcessFactoryRegistry.hpp"
#include "FileProcessFactoryBase.hpp"
#include "Logger.hpp"

#include <stdexcept>

FileProcessFactoryRegistry & FileProcessFactoryRegistry::Instance(void)
{
    Logger::Log(LogLevel::Debug, "FileProcessFactoryRegistry::Instance() called");
    static FileProcessFactoryRegistry instance;
    return instance;
}

void FileProcessFactoryRegistry::RegisterDefaultFactories(void)
{
    Logger::Log(LogLevel::Debug, "FileProcessFactoryRegistry::RegisterDefaultFactories() called");
    RegisterFactory(".pdb", []() { return std::make_unique<ModelObjectFactory>(); });
    RegisterFactory(".cif", []() { return std::make_unique<ModelObjectFactory>(); });
    RegisterFactory(".mrc", []() { return std::make_unique<MapObjectFactory>(); });
    RegisterFactory(".map", []() { return std::make_unique<MapObjectFactory>(); });
}

void FileProcessFactoryRegistry::RegisterFactory(
    const std::string & extension, std::function<std::unique_ptr<FileProcessFactoryBase>()> creator)
{
    Logger::Log(LogLevel::Debug,
        "FileProcessFactoryRegistry::RegisterFactory() called for extension: " + extension);
    m_factory_map[extension] = std::move(creator);
}

std::unique_ptr<FileProcessFactoryBase> FileProcessFactoryRegistry::CreateFactory(
    const std::string & extension) const
{
    Logger::Log(LogLevel::Debug,
        "FileProcessFactoryRegistry::CreateFactory() called for extension: " + extension);
    auto iter{ m_factory_map.find(extension) };
    if (iter == m_factory_map.end())
    {
        throw std::runtime_error("Unsupported file format: " + extension);
    }
    return (iter->second)();
}
