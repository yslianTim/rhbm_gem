#include <gtest/gtest.h>

#include <memory>
#include <stdexcept>

#include "command/detail/CommandObjectCache.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(CommandObjectCacheTest, StoresAndLoadsTypedTopLevelObjects)
{
    rg::command_detail::CommandObjectCache cache;
    auto model{ std::shared_ptr<rg::ModelObject>(data_test::MakeModelWithBond()) };
    auto map{ std::make_shared<rg::MapObject>(data_test::MakeMapObject()) };

    cache.PutModel("model", model);
    cache.PutMap("map", map);

    EXPECT_EQ(cache.GetKind("model"), rg::command_detail::CommandObjectCache::ObjectKind::Model);
    EXPECT_EQ(cache.GetKind("map"), rg::command_detail::CommandObjectCache::ObjectKind::Map);
    EXPECT_EQ(cache.GetModel("model"), model);
    EXPECT_EQ(cache.GetMap("map"), map);
}

TEST(CommandObjectCacheTest, TypedAccessThrowsOnMismatch)
{
    rg::command_detail::CommandObjectCache cache;
    auto model{ std::shared_ptr<rg::ModelObject>(data_test::MakeModelWithBond()) };
    cache.PutModel("model", model);

    EXPECT_THROW((void)cache.GetMap("model"), std::runtime_error);
    EXPECT_THROW((void)cache.GetKind("missing"), std::runtime_error);
}
