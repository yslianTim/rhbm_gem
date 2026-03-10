#pragma once

#include <memory>

namespace rhbm_gem {

class DataObjectDAOFactoryRegistry;
class FileFormatRegistry;
class FileProcessFactoryResolver;

class DataIoServices
{
public:
    DataIoServices(const DataIoServices &) = default;
    DataIoServices(DataIoServices &&) noexcept = default;
    DataIoServices & operator=(const DataIoServices &) = default;
    DataIoServices & operator=(DataIoServices &&) noexcept = default;
    ~DataIoServices();

    static DataIoServices BuildDefault();

private:
    const std::shared_ptr<const FileFormatRegistry> & FileFormatRegistryPtr() const;
    const std::shared_ptr<const FileProcessFactoryResolver> & FileProcessFactoryResolverPtr() const;
    const std::shared_ptr<const DataObjectDAOFactoryRegistry> & DaoFactoryRegistryPtr() const;

    struct Impl;
    std::shared_ptr<Impl> m_impl;

    explicit DataIoServices(std::shared_ptr<Impl> impl);

    friend class DataObjectManager;
};

} // namespace rhbm_gem
