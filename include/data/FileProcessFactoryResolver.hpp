#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace rhbm_gem {

class FileProcessFactoryBase;

class FileProcessFactoryResolver
{
public:
    virtual ~FileProcessFactoryResolver() = default;
    virtual std::unique_ptr<FileProcessFactoryBase> CreateFactory(
        const std::string & extension) const = 0;
};

class DefaultFileProcessFactoryResolver : public FileProcessFactoryResolver
{
public:
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(
        const std::string & extension) const override;
};

class OverrideableFileProcessFactoryResolver : public FileProcessFactoryResolver
{
    std::unordered_map<
        std::string,
        std::function<std::unique_ptr<FileProcessFactoryBase>()>> m_factory_map;

public:
    void RegisterFactory(
        const std::string & extension,
        std::function<std::unique_ptr<FileProcessFactoryBase>()> creator);
    void UnregisterFactory(const std::string & extension);
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(
        const std::string & extension) const override;
};

std::shared_ptr<const FileProcessFactoryResolver> CreateDefaultFileProcessFactoryResolver();

} // namespace rhbm_gem
