#pragma once

#include <unordered_map>
#include <string>
#include <functional>
#include <memory>

namespace rhbm_gem {

class FileProcessFactoryBase;

class FileProcessFactoryRegistry
{
    std::unordered_map<
        std::string,
        std::function<std::unique_ptr<FileProcessFactoryBase>()>> m_factory_map;

public:
    static FileProcessFactoryRegistry & Instance();
    // Compatibility helper for tests or explicit bulk overrides. Normal built-in
    // dispatch should rely on FileFormatRegistry fallback in CreateFactory().
    void RegisterDefaultFactories();
    void RegisterFactory(
        const std::string & extension,
        std::function<std::unique_ptr<FileProcessFactoryBase>()> creator);
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(const std::string & extension) const;

private:
    FileProcessFactoryRegistry() = default;
};

} // namespace rhbm_gem
