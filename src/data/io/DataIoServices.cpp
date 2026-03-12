#include <rhbm_gem/data/io/DataIoServices.hpp>

#include "internal/io/file/FileFormatRegistry.hpp"
#include "internal/io/file/FileProcessFactoryResolver.hpp"
#include "internal/io/sqlite/DataObjectDAOFactoryRegistry.hpp"

namespace rhbm_gem {

struct DataIoServices::Impl
{
    std::shared_ptr<const FileFormatRegistry> file_format_registry;
    std::shared_ptr<const FileProcessFactoryResolver> file_process_factory_resolver;
    std::shared_ptr<const DataObjectDAOFactoryRegistry> dao_factory_registry;
};

DataIoServices::DataIoServices(std::shared_ptr<Impl> impl) :
    m_impl{ std::move(impl) }
{
}

DataIoServices::~DataIoServices() = default;

DataIoServices DataIoServices::BuildDefault()
{
    auto impl{ std::make_shared<Impl>() };

    auto file_format_registry{
        std::make_shared<FileFormatRegistry>(BuildDefaultFileFormatRegistry()) };
    auto default_resolver{
        std::make_shared<DefaultFileProcessFactoryResolver>(*file_format_registry) };
    auto dao_factory_registry{ std::make_shared<DataObjectDAOFactoryRegistry>() };
    RegisterDataObjectDaos(*dao_factory_registry);

    impl->file_format_registry = std::move(file_format_registry);
    impl->file_process_factory_resolver = std::move(default_resolver);
    impl->dao_factory_registry = std::move(dao_factory_registry);

    return DataIoServices{ std::move(impl) };
}

const std::shared_ptr<const FileFormatRegistry> & DataIoServices::FileFormatRegistryPtr() const
{
    return m_impl->file_format_registry;
}

const std::shared_ptr<const FileProcessFactoryResolver> &
DataIoServices::FileProcessFactoryResolverPtr() const
{
    return m_impl->file_process_factory_resolver;
}

const std::shared_ptr<const DataObjectDAOFactoryRegistry> &
DataIoServices::DaoFactoryRegistryPtr() const
{
    return m_impl->dao_factory_registry;
}

} // namespace rhbm_gem
