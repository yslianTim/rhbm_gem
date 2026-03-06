#include <gtest/gtest.h>

#include <algorithm>
#include <string_view>

#include "BuiltInCommandBindingInternal.hpp"
#include "BuiltInCommandCatalogInternal.hpp"
#include "HRLModelTestCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "MapVisualizationCommand.hpp"
#include "PositionEstimationCommand.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename CommandType>
void ExpectBindingNameMatchesCatalog()
{
    const auto & catalog{ rg::BuiltInCommandCatalog() };
    const auto descriptor_iter{
        std::find_if(
            catalog.begin(),
            catalog.end(),
            [](const rg::CommandDescriptor & descriptor)
            {
                return descriptor.id == CommandType::kCommandId;
            })
    };
    ASSERT_NE(descriptor_iter, catalog.end());
    EXPECT_EQ(
        std::string_view{ descriptor_iter->python_binding_name },
        rg::BuiltInCommandBindingName<CommandType>::value);
}

} // namespace

TEST(BuiltInCommandBindingNameTest, BindingNamesMatchBuiltInCatalogMetadata)
{
    ExpectBindingNameMatchesCatalog<rg::PotentialAnalysisCommand>();
    ExpectBindingNameMatchesCatalog<rg::PotentialDisplayCommand>();
    ExpectBindingNameMatchesCatalog<rg::ResultDumpCommand>();
    ExpectBindingNameMatchesCatalog<rg::MapSimulationCommand>();
    ExpectBindingNameMatchesCatalog<rg::MapVisualizationCommand>();
    ExpectBindingNameMatchesCatalog<rg::PositionEstimationCommand>();
    ExpectBindingNameMatchesCatalog<rg::HRLModelTestCommand>();
}
