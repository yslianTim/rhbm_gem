#include <gtest/gtest.h>

#include <array>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "BuiltInCommandCatalogInternal.hpp"

namespace rg = rhbm_gem;

TEST(BuiltInCommandCatalogTest, BuiltInOrderIsStable)
{
    const auto & catalog{ rg::BuiltInCommandCatalog() };
    const std::vector<std::string> expected_names{
        "potential_analysis",
        "potential_display",
        "result_dump",
        "map_simulation",
        "map_visualization",
        "position_estimation",
        "model_test"
    };

    ASSERT_EQ(catalog.size(), expected_names.size());
    for (std::size_t index = 0; index < expected_names.size(); ++index)
    {
        EXPECT_EQ(catalog[index].name, expected_names[index]);
    }
}

TEST(BuiltInCommandCatalogTest, CommandIdsAndNamesAreUnique)
{
    const auto & catalog{ rg::BuiltInCommandCatalog() };
    std::unordered_set<int> unique_ids;
    std::unordered_set<std::string> unique_names;

    for (const auto & descriptor : catalog)
    {
        unique_ids.insert(static_cast<int>(descriptor.id));
        unique_names.emplace(descriptor.name);
    }

    EXPECT_EQ(unique_ids.size(), catalog.size());
    EXPECT_EQ(unique_names.size(), catalog.size());
}

TEST(BuiltInCommandCatalogTest, FactoriesAreNonNull)
{
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        EXPECT_NE(descriptor.factory, nullptr) << descriptor.name;
    }
}

TEST(BuiltInCommandCatalogTest, AllBuiltInCommandsProvidePythonBindingNames)
{
    static constexpr std::array<std::string_view, 7> kExpectedPythonCommandNames{
        "potential_analysis",
        "potential_display",
        "result_dump",
        "map_simulation",
        "map_visualization",
        "position_estimation",
        "model_test"
    };

    std::vector<std::string> actual_names;
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        EXPECT_FALSE(descriptor.python_binding_name.empty()) << descriptor.name;
        actual_names.emplace_back(descriptor.name);
    }

    ASSERT_EQ(actual_names.size(), kExpectedPythonCommandNames.size());
    for (std::size_t index = 0; index < kExpectedPythonCommandNames.size(); ++index)
    {
        EXPECT_EQ(actual_names[index], kExpectedPythonCommandNames[index]);
    }
}

TEST(BuiltInCommandCatalogTest, DatabaseEnabledCommandSetMatchesExpectedSurface)
{
    static constexpr std::array<std::string_view, 3> kExpectedDatabaseCommandNames{
        "potential_analysis",
        "potential_display",
        "result_dump"
    };

    std::vector<std::string> actual_names;
    for (const auto & descriptor : rg::BuiltInCommandCatalog())
    {
        if (!rg::UsesDatabaseAtRuntime(descriptor.surface)) continue;
        actual_names.emplace_back(descriptor.name);
    }

    ASSERT_EQ(actual_names.size(), kExpectedDatabaseCommandNames.size());
    for (std::size_t index = 0; index < kExpectedDatabaseCommandNames.size(); ++index)
    {
        EXPECT_EQ(actual_names[index], kExpectedDatabaseCommandNames[index]);
    }
}
