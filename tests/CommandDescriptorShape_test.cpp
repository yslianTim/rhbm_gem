#include <gtest/gtest.h>

#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

#include "BuiltInCommandCatalogInternal.hpp"
#include "CommandBase.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename CommandType>
void ExpectDescriptorMatchesConcreteType(rg::CommandId command_id)
{
    const auto & descriptor{ rg::FindCommandDescriptor(command_id) };
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
