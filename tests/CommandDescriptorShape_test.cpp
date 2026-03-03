#include <gtest/gtest.h>

#include "BuiltInCommandCatalog.hpp"
#include "CommandBase.hpp"

namespace rg = rhbm_gem;

TEST(CommandDescriptorShapeTest, BuiltInFactoriesAreNonNullAndConstructCommands)
{
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        ASSERT_NE(descriptor.factory, nullptr) << descriptor.name;

        auto command{ descriptor.factory() };
        ASSERT_NE(command, nullptr) << descriptor.name;
        EXPECT_EQ(command->GetCommandId(), descriptor.id) << descriptor.name;
    }
}
