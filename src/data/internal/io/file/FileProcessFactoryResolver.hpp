#pragma once

#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace rhbm_gem {

class FileProcessFactoryBase;
class FileFormatRegistry;

class FileProcessFactoryResolver
{
public:
    virtual ~FileProcessFactoryResolver() = default;
    virtual std::unique_ptr<FileProcessFactoryBase> CreateFactory(
        const std::string & extension) const = 0;
};

class DefaultFileProcessFactoryResolver : public FileProcessFactoryResolver
{
    const FileFormatRegistry & m_file_format_registry;

public:
    explicit DefaultFileProcessFactoryResolver(const FileFormatRegistry & file_format_registry);
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(
        const std::string & extension) const override;
};

class OverrideableFileProcessFactoryResolver : public FileProcessFactoryResolver
{
    std::shared_ptr<const FileProcessFactoryResolver> m_fallback_resolver;
    mutable std::mutex m_mutex;
    std::unordered_map<
        std::string,
        std::function<std::unique_ptr<FileProcessFactoryBase>()>> m_factory_map;

public:
    explicit OverrideableFileProcessFactoryResolver(
        std::shared_ptr<const FileProcessFactoryResolver> fallback_resolver);
    void RegisterFactory(
        const std::string & extension,
        std::function<std::unique_ptr<FileProcessFactoryBase>()> creator);
    void UnregisterFactory(const std::string & extension);
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(
        const std::string & extension) const override;
};

} // namespace rhbm_gem
