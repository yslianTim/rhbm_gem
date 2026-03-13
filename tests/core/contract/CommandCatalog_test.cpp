#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "internal/CommandCatalog.hpp"
#include <CLI/CLI.hpp>

namespace rg = rhbm_gem;

namespace {

struct ExpectedCommandMetadata
{
    rg::CommandId id;
    std::string_view name;
    rg::CommonOptionMask common_options;
    bool uses_database;
};

constexpr std::array<ExpectedCommandMetadata, 7> kExpectedCommandMetadata{{
    {
        rg::CommandId::PotentialAnalysis,
        "potential_analysis",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::DatabaseWorkflow),
        true
    },
    {
        rg::CommandId::PotentialDisplay,
        "potential_display",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::DatabaseWorkflow),
        true
    },
    {
        rg::CommandId::ResultDump,
        "result_dump",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::DatabaseWorkflow),
        true
    },
    {
        rg::CommandId::MapSimulation,
        "map_simulation",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::FileWorkflow),
        false
    },
    {
        rg::CommandId::MapVisualization,
        "map_visualization",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::FileWorkflow),
        false
    },
    {
        rg::CommandId::PositionEstimation,
        "position_estimation",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::FileWorkflow),
        false
    },
    {
        rg::CommandId::ModelTest,
        "model_test",
        rg::CommonOptionMaskForProfile(rg::CommonOptionProfile::FileWorkflow),
        false
    }
}};

const rg::CommandDescriptor * FindDescriptor(rg::CommandId command_id)
{
    const auto & catalog{ rg::CommandCatalog() };
    const auto iter{
        std::find_if(
            catalog.begin(),
            catalog.end(),
            [command_id](const rg::CommandDescriptor & descriptor)
            {
                return descriptor.id == command_id;
            })
    };
    return iter == catalog.end() ? nullptr : &(*iter);
}

void ExpectRuntimeBinderConstructs(rg::CommandId command_id)
{
    const auto * descriptor{ FindDescriptor(command_id) };
    ASSERT_NE(descriptor, nullptr) << static_cast<int>(command_id);
    ASSERT_NE(descriptor->bind_runtime, nullptr) << descriptor->name;

    CLI::App subcommand{ std::string(descriptor->name) };
    auto runner{ descriptor->bind_runtime(&subcommand) };
    EXPECT_TRUE(static_cast<bool>(runner));
}

} // namespace

TEST(CommandCatalogTest, CatalogMatchesExpectedMetadataAndOrder)
{
    const auto & catalog{ rg::CommandCatalog() };
    ASSERT_EQ(catalog.size(), kExpectedCommandMetadata.size());

    std::unordered_set<int> unique_ids;
    std::unordered_set<std::string> unique_names;
    for (std::size_t index = 0; index < kExpectedCommandMetadata.size(); ++index)
    {
        const auto & descriptor{ catalog[index] };
        const auto & expected{ kExpectedCommandMetadata[index] };

        EXPECT_EQ(descriptor.id, expected.id);
        EXPECT_EQ(std::string_view{ descriptor.name }, expected.name);
        EXPECT_EQ(descriptor.common_options, expected.common_options);
        EXPECT_EQ(
            rg::HasCommonOption(descriptor.common_options, rg::CommonOption::Database),
            expected.uses_database);
        EXPECT_TRUE(rg::HasCommonOption(descriptor.common_options, rg::CommonOption::Threading));
        EXPECT_TRUE(rg::HasCommonOption(descriptor.common_options, rg::CommonOption::Verbose));
        EXPECT_TRUE(rg::HasCommonOption(descriptor.common_options, rg::CommonOption::OutputFolder));

        unique_ids.insert(static_cast<int>(descriptor.id));
        unique_names.emplace(descriptor.name);
    }

    EXPECT_EQ(unique_ids.size(), catalog.size());
    EXPECT_EQ(unique_names.size(), catalog.size());
}

TEST(CommandCatalogTest, RuntimeBindersArePresentForAllCommands)
{
    ExpectRuntimeBinderConstructs(rg::CommandId::PotentialAnalysis);
    ExpectRuntimeBinderConstructs(rg::CommandId::PotentialDisplay);
    ExpectRuntimeBinderConstructs(rg::CommandId::ResultDump);
    ExpectRuntimeBinderConstructs(rg::CommandId::MapSimulation);
    ExpectRuntimeBinderConstructs(rg::CommandId::MapVisualization);
    ExpectRuntimeBinderConstructs(rg::CommandId::PositionEstimation);
    ExpectRuntimeBinderConstructs(rg::CommandId::ModelTest);
}
