#include <gtest/gtest.h>

#include <algorithm>

#include <rhbm_gem/core/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/MapSimulationCommand.hpp>
#include <rhbm_gem/core/MapVisualizationCommand.hpp>
#include <rhbm_gem/core/PositionEstimationCommand.hpp>
#include <rhbm_gem/core/PotentialAnalysisCommand.hpp>
#include <rhbm_gem/core/PotentialDisplayCommand.hpp>
#include <rhbm_gem/core/ResultDumpCommand.hpp>

#include "internal/BuiltInCommandCatalogInternal.hpp"
#include <rhbm_gem/core/CommandBase.hpp>

namespace rg = rhbm_gem;

namespace {

template <typename CommandType>
void ExpectDescriptorMatchesConcreteType(rg::CommandId command_id)
{
    const auto & catalog{ rg::BuiltInCommandCatalog() };
    const auto descriptor_iter{
        std::find_if(
            catalog.begin(),
            catalog.end(),
            [command_id](const rg::CommandDescriptor & descriptor)
            {
                return descriptor.id == command_id;
            })
    };
    ASSERT_NE(descriptor_iter, catalog.end()) << static_cast<int>(command_id);
    const auto & descriptor{ *descriptor_iter };
    ASSERT_NE(descriptor.factory, nullptr) << descriptor.name;
    EXPECT_EQ(descriptor.id, CommandType::kCommandId) << descriptor.name;
    EXPECT_EQ(descriptor.common_options, CommandType::kCommonOptions) << descriptor.name;

    auto command{ descriptor.factory() };
    ASSERT_NE(command, nullptr) << descriptor.name;
    EXPECT_NE(dynamic_cast<CommandType *>(command.get()), nullptr) << descriptor.name;
}

} // namespace

TEST(CommandDescriptorShapeTest, BuiltInFactoriesAreNonNullAndConstructCommandsWithMatchingMetadata)
{
    ExpectDescriptorMatchesConcreteType<rg::PotentialAnalysisCommand>(
        rg::CommandId::PotentialAnalysis);
    ExpectDescriptorMatchesConcreteType<rg::PotentialDisplayCommand>(
        rg::CommandId::PotentialDisplay);
    ExpectDescriptorMatchesConcreteType<rg::ResultDumpCommand>(
        rg::CommandId::ResultDump);
    ExpectDescriptorMatchesConcreteType<rg::MapSimulationCommand>(
        rg::CommandId::MapSimulation);
    ExpectDescriptorMatchesConcreteType<rg::MapVisualizationCommand>(
        rg::CommandId::MapVisualization);
    ExpectDescriptorMatchesConcreteType<rg::PositionEstimationCommand>(
        rg::CommandId::PositionEstimation);
    ExpectDescriptorMatchesConcreteType<rg::HRLModelTestCommand>(
        rg::CommandId::ModelTest);
}
