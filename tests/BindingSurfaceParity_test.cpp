#include <gtest/gtest.h>

#include "BuiltInCommandCatalog.hpp"

namespace rg = rhbm_gem;

TEST(BindingSurfaceParityTest, CatalogExposureMatchesPythonBindingMetadata)
{
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        EXPECT_EQ(
            rg::IsPythonPublic(descriptor.binding_exposure),
            !descriptor.python_binding_name.empty())
            << descriptor.name;
    }
}
