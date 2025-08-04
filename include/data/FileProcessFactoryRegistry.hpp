#pragma once

#include <unordered_map>
#include <string>
#include <functional>
#include <memory>

class FileProcessFactoryBase;

class FileProcessFactoryRegistry
{
    std::unordered_map<
        std::string,
        std::function<std::unique_ptr<FileProcessFactoryBase>()>> m_factory_map;

public:
    static FileProcessFactoryRegistry & Instance(void);
    void RegisterDefaultFactories(void);
    void RegisterFactory(
        const std::string & extension,
        std::function<std::unique_ptr<FileProcessFactoryBase>()> creator);
    std::unique_ptr<FileProcessFactoryBase> CreateFactory(const std::string & extension) const;

private:
    FileProcessFactoryRegistry(void) = default;
    
};
