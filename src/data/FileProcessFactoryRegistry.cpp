#include "FileProcessFactoryRegistry.hpp"
#include "FileProcessFactoryBase.hpp"

namespace rhbm_gem {

FileProcessFactoryRegistry & FileProcessFactoryRegistry::Instance()
{
    static FileProcessFactoryRegistry instance;
    return instance;
}

void FileProcessFactoryRegistry::UnregisterFactory(const std::string & extension)
{
    m_resolver.UnregisterFactory(extension);
}

void FileProcessFactoryRegistry::RegisterFactory(
    const std::string & extension, std::function<std::unique_ptr<FileProcessFactoryBase>()> creator)
{
    m_resolver.RegisterFactory(extension, std::move(creator));
}

std::unique_ptr<FileProcessFactoryBase> FileProcessFactoryRegistry::CreateFactory(
    const std::string & extension) const
{
    return m_resolver.CreateFactory(extension);
}

} // namespace rhbm_gem
