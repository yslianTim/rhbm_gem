#include <gtest/gtest.h>

#include <string>
#include <unordered_set>
#include <vector>

#include "BuiltInCommandCatalog.hpp"

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
