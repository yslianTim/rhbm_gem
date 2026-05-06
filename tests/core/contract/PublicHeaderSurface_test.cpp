#include <gtest/gtest.h>

#include <rhbm_gem/core/command/CommandTypes.hpp>

TEST(PublicHeaderSurfaceTest, CommandApiExposesStableDefaultDatabasePathHelper) {
    const auto default_path{ rhbm_gem::GetDefaultDatabasePath() };
    EXPECT_FALSE(default_path.empty());
    EXPECT_EQ(default_path.filename(), "database.sqlite");
}
