#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

#include "internal/CommandCatalogInternal.hpp"
#include <rhbm_gem/core/command/CommandBase.hpp>
#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/command/MapSimulationCommand.hpp>
#include <rhbm_gem/core/command/MapVisualizationCommand.hpp>
#include <rhbm_gem/core/command/PositionEstimationCommand.hpp>
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/command/PotentialDisplayCommand.hpp>
#include <rhbm_gem/core/command/ResultDumpCommand.hpp>

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
        rg::PotentialAnalysisCommand::kCommonOptions,
        true
    },
    {
        rg::CommandId::PotentialDisplay,
        "potential_display",
        rg::PotentialDisplayCommand::kCommonOptions,
        true
    },
    {
        rg::CommandId::ResultDump,
        "result_dump",
        rg::ResultDumpCommand::kCommonOptions,
        true
    },
    {
        rg::CommandId::MapSimulation,
        "map_simulation",
        rg::MapSimulationCommand::kCommonOptions,
        false
    },
    {
        rg::CommandId::MapVisualization,
        "map_visualization",
        rg::MapVisualizationCommand::kCommonOptions,
        false
    },
    {
        rg::CommandId::PositionEstimation,
        "position_estimation",
        rg::PositionEstimationCommand::kCommonOptions,
        false
    },
    {
        rg::CommandId::ModelTest,
        "model_test",
        rg::HRLModelTestCommand::kCommonOptions,
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

template <typename CommandType>
void ExpectFactoryConstructs(rg::CommandId command_id)
{
    const auto * descriptor{ FindDescriptor(command_id) };
    ASSERT_NE(descriptor, nullptr) << static_cast<int>(command_id);
    ASSERT_NE(descriptor->factory, nullptr) << descriptor->name;
    EXPECT_EQ(descriptor->id, CommandType::kCommandId) << descriptor->name;
    EXPECT_EQ(descriptor->common_options, CommandType::kCommonOptions) << descriptor->name;
    auto command{ descriptor->factory() };
    ASSERT_NE(command, nullptr) << descriptor->name;
    EXPECT_NE(dynamic_cast<CommandType *>(command.get()), nullptr) << descriptor->name;
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

TEST(CommandCatalogTest, FactoriesConstructCommandsWithMatchingMetadata)
{
    ExpectFactoryConstructs<rg::PotentialAnalysisCommand>(rg::CommandId::PotentialAnalysis);
    ExpectFactoryConstructs<rg::PotentialDisplayCommand>(rg::CommandId::PotentialDisplay);
    ExpectFactoryConstructs<rg::ResultDumpCommand>(rg::CommandId::ResultDump);
    ExpectFactoryConstructs<rg::MapSimulationCommand>(rg::CommandId::MapSimulation);
    ExpectFactoryConstructs<rg::MapVisualizationCommand>(rg::CommandId::MapVisualization);
    ExpectFactoryConstructs<rg::PositionEstimationCommand>(rg::CommandId::PositionEstimation);
    ExpectFactoryConstructs<rg::HRLModelTestCommand>(rg::CommandId::ModelTest);
}
