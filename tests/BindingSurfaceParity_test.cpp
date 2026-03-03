#include <gtest/gtest.h>

#include <unordered_set>

#include "BuiltInCommandCatalog.hpp"

namespace rg = rhbm_gem;

TEST(BindingSurfaceParityTest, CatalogExposureMatchesExpectedPublicCommandSet)
{
    const std::unordered_set<rg::CommandId> expected_python_public{
        rg::CommandId::PotentialAnalysis,
        rg::CommandId::PotentialDisplay,
        rg::CommandId::ResultDump,
        rg::CommandId::MapSimulation
    };

    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        const bool is_expected_public{
            expected_python_public.find(descriptor.id) != expected_python_public.end()
        };
        EXPECT_EQ(rg::IsPythonPublic(descriptor.binding_exposure), is_expected_public)
            << descriptor.name;
    }
}
